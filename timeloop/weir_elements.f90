module weir_elements

    use define_globals
    use define_keys
    use define_indexes
    use define_settings, only: setting
    use common_elements
    use adjust
    use utility, only: util_CLprint
    use utility_crash, only: util_crashpoint

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% Computes flow and head for weir elements
    !%----------------------------------------------------------------------------- 

    implicit none

    private

    public :: weir_toplevel
    public :: weir_set_setting

    contains
    !%==========================================================================
    !% PUBLIC
    !%==========================================================================
    !%
    subroutine weir_toplevel  (eIdx)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes diagnostic flow and head delta across a weir.
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: eIdx  !% must be a single element ID
        logical, pointer :: isSurcharged
        
        character(64) :: subroutine_name = 'weir_toplevel'
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%weir_elements) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------
        isSurcharged => elemYN(eIdx,eYN_isSurcharged)

        !% --- set the Setting for the fractional open
        ! call weir_set_setting (eIdx)  !% ss20220701 -- weir setting is already being set in control_update_setting subroutine

        !% get the flow direction and element head
        call  common_head_and_flowdirection_singular &
            (eIdx, esr_Weir_Zcrest, esr_Weir_NominalDownstreamHead, esi_Weir_FlowDirection)

        !% find effective head difference accross weir element
        call weir_effective_head_delta (eIdx)
        
        !% find flow on weir element
        if (isSurcharged) then
            call weir_surcharge_flow (eIdx)
        else
            call weir_flow (eIdx, esr_Weir_EffectiveHeadDelta, .true., .true.)
        endif
        
        !% update weir geometry from head
        call weir_geometry_update (eIdx)
        
        !% update velocity from flowrate and area
        call common_velocity_from_flowrate_singular (eIdx)
        
        if (setting%Debug%File%weir_elements)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine weir_toplevel    
!%
!%==========================================================================
!%==========================================================================   
!%
    subroutine weir_set_setting (eIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% adjusts weir values for 0 <= er_setting <= 1.0
        !% patterned after EPA-SWMM link.c/weir_setSetting
        !%------------------------------------------------------------------
        !% Declarations:
        integer, intent(in) :: eIdx

        real(8), pointer :: FullDepth, EffectiveFullDepth
        real(8), pointer :: CurrentSetting, TargetSetting
        !%------------------------------------------------------------------

        FullDepth          => elemSR(eIdx,esr_Weir_FullDepth)
        EffectiveFullDepth => elemSR(eIdx,esr_Weir_EffectiveFullDepth)
        CurrentSetting     => elemR(eIdx,er_Setting)
        TargetSetting      => elemR(eIdx,er_TargetSetting)

        !% --- instantaneous adjustment
        CurrentSetting = TargetSetting 

        !% --- error check
        !%     EPA-SWMM allows the weir setting to be between 0.0 and 1.0
        if (.not. ((CurrentSetting .ge. 0.0) .and. (CurrentSetting .le. 1.0))) then
            print *, 'CODE ERROR: orifice element has er_Setting that is not between 0.0 and 1.0'
            call util_crashpoint(668723)
        end if

        !% find effective weir opening
        EffectiveFullDepth = FullDepth * CurrentSetting

    end subroutine weir_set_setting
!% 
!%==========================================================================   
!% PRIVATE
!%==========================================================================       
!%  
     subroutine weir_effective_head_delta (eIdx)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes the effective head difference flowing over the top of a weir
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: eIdx !% single ID of element
        integer, pointer :: FlowDirection
        real(8), pointer :: EffectiveHeadDelta, Head, Zcrown, Zcrest
        real(8), pointer :: NominalDownstreamHead, CurrentSetting, EffectiveFullDepth
        logical, pointer :: CanSurcharge, IsSurcharged, hasFlapGate
        real(8) :: Zmidpt
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        !% input
        Head                  => elemR(eIdx,er_Head)
        FlowDirection         => elemSI(eIdx,esi_Weir_FlowDirection)
        hasFlapGate           => elemYN(eIdx,eYN_hasFlapGate)
        !% output
        EffectiveHeadDelta    => elemSR(eIdx,esr_Weir_EffectiveHeadDelta)
        EffectiveFullDepth    => elemSR(eIdx,esr_Weir_EffectiveFullDepth)
        Zcrown                => elemSR(eIdx,esr_Weir_Zcrown)
        Zcrest                => elemSR(eIdx,esr_Weir_Zcrest)
        NominalDownstreamHead => elemSR(eIdx,esr_Weir_NominalDownstreamHead)
        CanSurcharge          => elemYN(eIdx,eYN_canSurcharge)
        IsSurcharged          => elemYN(eIdx,eYN_isSurcharged)
        CurrentSetting        => elemR(eIdx,er_Setting)
        
        !% setting default surcharge condition as false
        IsSurcharged = .false.
        !%-----------------------------------------------------------------------------

        !% adjust weir crest height for partially open weir
        Zcrest = Zcrest + (oneR - CurrentSetting) * EffectiveFullDepth

        !% if the weir has a flapgate and the flow direction is reverse
        !% set EffectiveHeadDelta to zero, chich will result in zero flows
        if (hasFlapGate .and. (FlowDirection < zeroI)) then
            EffectiveHeadDelta = zeroR
        else
            if (Head <= Zcrest) then
                EffectiveHeadDelta = zeroR
            else
                EffectiveHeadDelta = Head - Zcrest
            end if
                
            if (Head > Zcrown) then
                !% use equivalent orifice head calculation if the weir can surcharge
                if (CanSurcharge) then
                    IsSurcharged = .true.
                    Zmidpt = (Zcrest + Zcrown) / twoR
                    if (NominalDownstreamHead < Zmidpt) then
                        EffectiveHeadDelta = Head - Zmidpt       
                    else
                        EffectiveHeadDelta = Head - NominalDownstreamHead    
                    endif  
                !% if the weir cannot surcharge, limit the head to height of weir opening
                else
                    EffectiveHeadDelta =  Zcrown - Zcrest
                end if      
            end if
        end if

    end subroutine weir_effective_head_delta
!%
!%========================================================================== 
!%==========================================================================    
!%  
     subroutine weir_surcharge_flow (eIdx)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% 
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: eIdx !% must be single element ID
        integer, pointer :: FlowDirection
        real(8), pointer :: Area, Flowrate, EffectiveFullDepth, Depth 
        real(8), pointer :: EffectiveHeadDelta, grav
        logical, pointer :: hasFlapGate
        real(8) :: CoeffOrifice, hLoss
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        Area               => elemR(eIdx,er_Area)
        Depth              => elemR(eIdx,er_Depth)
        Flowrate           => elemR(eIdx,er_Flowrate)
        hasFlapGate        => elemYN(eiDx,eYN_hasFlapGate) 
        FlowDirection      => elemSI(eIdx,esi_Weir_FlowDirection)
        EffectiveFullDepth => elemSR(eIdx,esr_Weir_EffectiveFullDepth)
        EffectiveHeadDelta => elemSR(eIdx,esr_Weir_EffectiveHeadDelta)
        grav               => setting%Constant%gravity
        !%-----------------------------------------------------------------------------
        ! get the flowrate for effective full depth without submergence correction
        call weir_flow(eIdx, esr_Weir_EffectiveFullDepth, .false., .false.)

        ! equivalent orifice flow coefficient for surcharge flow
        CoeffOrifice = Flowrate / sqrt(EffectiveFullDepth/twoR)
        
        !% old flowrate is overwritten by new surcharged flowrate
        Flowrate = FlowDirection * CoeffOrifice * sqrt(EffectiveHeadDelta)

        !% update the weir geometry to find the weir opening
        call weir_get_open_area (eIdx)

        !% apply aramco adjustments for flap gate head loss
        if (hasFlapGate) then
            call weir_get_flapgate_headLoss (eIdx, esr_Weir_EffectiveHeadDelta)
            !% recalculate the flowrate based on new adjusted head
            Flowrate = FlowDirection * CoeffOrifice * sqrt(EffectiveHeadDelta)
        end if

    end subroutine weir_surcharge_flow
!%
!%========================================================================== 
!%==========================================================================    
!%  
    subroutine weir_flow (eIdx, inCol, ApplySubmergenceCorrection, ApplyHeadlossCorrection)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes flow for standard weir types
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: eIdx, inCol
        logical, intent(in) :: ApplySubmergenceCorrection, ApplyHeadlossCorrection
        integer, pointer :: SpecificWeirType, EndContractions, FlowDirection
        real(8), pointer :: Flowrate, Head, EffectiveHeadDelta, CurrentSetting, fullDepth
        real(8), pointer :: RectangularBreadth, TrapezoidalBreadth
        real(8), pointer :: TriangularSideSlope, TrapezoidalLeftSlope, TrapezoidalRightSlope
        real(8), pointer :: CoeffTriangular, CoeffRectangular
        real(8), pointer :: WeirExponent, WeirExponentVNotch
        real(8), pointer :: WeirContractionFactor, VillemonteExponent, WeirCrestExponent
        real(8), pointer :: NominalDsHead, Zcrest
        logical, pointer :: hasFlapGate
        real(8) :: Zmidpt, CrestLength, SubCorrectionTriangular, SubCorrectionRectangular
        real(8) :: FlowRect, FlowTriang, ratio
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        SpecificWeirType => elemSI(eIdx,esi_Weir_SpecificType)
        EndContractions  => elemSI(eIdx,esi_Weir_EndContractions)
        FlowDirection    => elemSI(eIdx,esi_Weir_FlowDirection)
        
        Head                  => elemR(eIdx,er_Head)
        Flowrate              => elemR(eIdx,er_Flowrate)
        CurrentSetting        => elemR(eIdx,er_Setting)
        hasFlapGate           => elemYN(eIdx,eYN_hasFlapGate)
        EffectiveHeadDelta    => elemSR(eIdx,inCol)
        Zcrest                => elemSR(eIdx,esr_Weir_Zcrest)
        RectangularBreadth    => elemSR(eIdx,esr_Weir_RectangularBreadth)
        TrapezoidalBreadth    => elemSR(eIdx,esr_Weir_TrapezoidalBreadth)
        TriangularSideSlope   => elemSR(eIdx,esr_Weir_TriangularSideSlope)
        TrapezoidalLeftSlope  => elemSR(eIdx,esr_Weir_TrapezoidalLeftSlope)
        TrapezoidalRightSlope => elemSR(eIdx,esr_Weir_TrapezoidalRightSlope)
        CoeffTriangular       => elemSR(eIdx,esr_Weir_Triangular)
        CoeffRectangular      => elemSR(eIdx,esr_Weir_Rectangular)
        NominalDsHead         => elemSR(eIdx,esr_Weir_NominalDownstreamHead)
        FullDepth             => elemSR(eIdx,esr_Weir_FullDepth)
        !%-----------------------------------------------------------------------------
        !% initializing default local Villemonte submergence correction factors as 1
        !% These are changed below if needed
        SubCorrectionTriangular = oneR
        SubCorrectionRectangular = oneR
         
        select case (SpecificWeirType)
            case (transverse_weir)
                WeirExponent          => Setting%Weir%Transverse%WeirExponent
                WeirContractionFactor => Setting%Weir%Transverse%WeirContractionFactor
                VillemonteExponent    => Setting%Weir%Transverse%VillemonteCorrectionExponent

                !% effective crest length due to contraction for tranverse weir             
                CrestLength = max(zeroR, &
                        RectangularBreadth - WeirContractionFactor * real(EndContractions,8) * EffectiveHeadDelta)

                Flowrate = real(FlowDirection,8) * CrestLength * CoeffRectangular  * (EffectiveHeadDelta ** WeirExponent) 

                if (hasFlapGate .and. ApplyHeadlossCorrection) then
                    call weir_get_flapgate_headLoss (eIdx, inCol)
                    !% recalculate the flowrate based on new adjusted head
                    Flowrate = real(FlowDirection,8) * CrestLength * CoeffRectangular  * (EffectiveHeadDelta ** WeirExponent)
                end if

                !% correction factor for nominal downstream submergence
                if ((NominalDsHead > Zcrest) .and. (ApplySubmergenceCorrection)) then
                    ratio = (NominalDsHead - Zcrest) / (Head - Zcrest)        
                    SubCorrectionRectangular = ((oneR - (ratio ** WeirExponent)) ** VillemonteExponent)
                endif
                !% apply submergence correction
                Flowrate =  SubCorrectionRectangular * Flowrate   

            case (side_flow)
                WeirExponent          => Setting%Weir%SideFlow%WeirExponent
                WeirContractionFactor => Setting%Weir%SideFlow%WeirContractionFactor
                WeirCrestExponent     => Setting%Weir%SideFlow%SideFlowWeirCrestExponent
                VillemonteExponent    => Setting%Weir%SideFlow%VillemonteCorrectionExponent

                !% effective crest length due to contraction for sideflow weir   
                CrestLength = max(zeroR, &
                        RectangularBreadth - WeirContractionFactor * real(EndContractions,8) * EffectiveHeadDelta)
                        
                if (FlowDirection > zeroR) then
                
                    Flowrate = real(FlowDirection,8) * (CrestLength ** &
                        WeirCrestExponent) * CoeffRectangular * (EffectiveHeadDelta ** WeirExponent)

                    if (hasFlapGate .and. ApplyHeadlossCorrection) then
                        call weir_get_flapgate_headLoss (eIdx, inCol)
                        !% recalculate the flowrate based on new adjusted head
                        Flowrate = real(FlowDirection,8) * CrestLength * CoeffRectangular  * (EffectiveHeadDelta ** WeirExponent)
                    end if

                    !% correction factor for nominal downstream submergence
                    if ((NominalDsHead > Zcrest) .and. (ApplySubmergenceCorrection)) then
                        ratio = (NominalDsHead - Zcrest) / (Head - Zcrest) 
                        SubCorrectionRectangular = ((oneR - (ratio ** WeirExponent)) ** VillemonteExponent)
                    endif

                    !% apply submergence correction
                    Flowrate =  SubCorrectionRectangular * Flowrate 
                
                else
                    !% under reverse flow condition, sideflow weir behaves like a transverse weir
                    !% correction factor for nominal downstream submergence
                    !% note: flap-gate headloss computation is not needed because if a flap-gate is present
                    !% the flow will be zero anyway

                    WeirExponent => Setting%Weir%Transverse%WeirExponent

                    if ((NominalDsHead > Zcrest) .and. (ApplySubmergenceCorrection)) then
                        ratio = (NominalDsHead - Zcrest) / (Head - Zcrest)      
                        SubCorrectionRectangular = ((oneR - (ratio ** WeirExponent)) ** VillemonteExponent)
                    endif
                
                    Flowrate = real(FlowDirection,8) * SubCorrectionRectangular * CrestLength * &
                        CoeffRectangular  * (EffectiveHeadDelta ** WeirExponent)
                endif  

            case (trapezoidal_weir)
                WeirExponentVNotch    => Setting%Weir%VNotch%WeirExponent
                WeirExponent          => Setting%Weir%Trapezoidal%WeirExponent
                WeirContractionFactor => Setting%Weir%Trapezoidal%WeirContractionFactor
                WeirCrestExponent     => Setting%Weir%Trapezoidal%SideFlowWeirCrestExponent
                VillemonteExponent    => Setting%Weir%Trapezoidal%VillemonteCorrectionExponent

                !% effective crest length due for trapezoidal weir (changes if a control is present)
                CrestLength = TrapezoidalBreadth +  (oneR - CurrentSetting) * FullDepth &
                            * (TrapezoidalLeftSlope + TrapezoidalRightSlope)

                FlowRect    = real(FlowDirection,8) * (CoeffRectangular * CrestLength &
                            * (EffectiveHeadDelta ** WeirExponent))
                FlowTriang  = real(FlowDirection,8) * (CoeffTriangular * ((TrapezoidalLeftSlope &
                            + TrapezoidalRightSlope) / twoR) * (EffectiveHeadDelta ** WeirExponentVNotch))

                Flowrate    = FlowRect + FlowTriang

                if (hasFlapGate .and. ApplyHeadlossCorrection) then
                    call weir_get_flapgate_headLoss (eIdx, inCol)
                    !% recalculate the flowrate based on new adjusted head
                    FlowRect    = real(FlowDirection,8) * (CoeffRectangular * CrestLength &
                                * (EffectiveHeadDelta ** WeirExponent))
                    FlowTriang  = real(FlowDirection,8) * (CoeffTriangular * ((TrapezoidalLeftSlope &
                                + TrapezoidalRightSlope) / twoR) * (EffectiveHeadDelta ** WeirExponentVNotch))
                end if
                
                !% correction factor for nominal downstream submergence
                if ((NominalDsHead > Zcrest) .and. (ApplySubmergenceCorrection)) then
                    ratio = (NominalDsHead - Zcrest) / (Head - Zcrest)   
                    SubCorrectionRectangular = ((oneR - (ratio ** WeirExponent)) **  VillemonteExponent)
                    SubCorrectionTriangular  = ((oneR - (ratio ** WeirExponentVNotch)) ** VillemonteExponent)
                endif

                !% apply submergence correction
                FlowRect   = SubCorrectionRectangular * FlowRect
                FlowTriang = SubCorrectionTriangular  * FlowTriang
                Flowrate   = FlowRect + FlowTriang
                      
            case (vnotch_weir)
                WeirExponent          => Setting%Weir%VNotch%WeirExponent
                WeirContractionFactor => Setting%Weir%VNotch%WeirContractionFactor
                WeirCrestExponent     => Setting%Weir%VNotch%SideFlowWeirCrestExponent
                VillemonteExponent    => Setting%Weir%VNotch%VillemonteCorrectionExponent
                
                Flowrate = real(FlowDirection,8) * CoeffTriangular * &
                        TriangularSideSlope * (EffectiveHeadDelta ** WeirExponent) 

                if (hasFlapGate .and. ApplyHeadlossCorrection) then
                    call weir_get_flapgate_headLoss (eIdx, inCol)
                    !% recalculate the flowrate based on new adjusted head
                    Flowrate = real(FlowDirection,8) * CoeffTriangular * &
                        TriangularSideSlope * (EffectiveHeadDelta ** WeirExponent)
                end if

                !% correction factor for nominal downstream submergence
                if ((NominalDsHead > Zcrest) .and. (ApplySubmergenceCorrection)) then
                    ratio = (NominalDsHead - Zcrest) / (Head - Zcrest)
                    SubCorrectionTriangular = ((oneR - (ratio ** WeirExponent)) ** VillemonteExponent)
                endif

                !% apply submergence correction
                Flowrate = SubCorrectionTriangular * Flowrate
                
            case default
                print *, 'CODE ERROR: unknown weir type, ', specificWeirType,'  in network'
                print *, 'which has key ',trim(reverseKey(specificWeirType))
                stop 848555
        end Select

    end subroutine weir_flow
!%
!%========================================================================== 
!%==========================================================================    
!%  
    subroutine weir_geometry_update (eIdx)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% 
        !% HACK -- it is not clear as yet what geometries we actually need. There's
        !% an important difference between the geometry of the flow over the weir
        !% and the geometry surrounding the weir.
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: eIdx
        real(8), pointer :: FullDepth, Head, Length, Zbottom,  Zcrown
        real(8), pointer :: Depth, Area, Volume, Topwidth, HydRadius
        real(8), pointer :: Perimeter, HydDepth,  Zcrest, ell, CurrentSetting
        real(8), pointer :: RectangularBreadth, TrapezoidalBreadth
        real(8), pointer :: TriangularSideSlope, TrapezoidalLeftSlope, TrapezoidalRightSlope
        integer, pointer :: SpecificWeirType
        logical, pointer :: IsSurcharged
        real(8)          :: z, zY
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        !% pointers
        SpecificWeirType => elemSI(eIdx,esi_Weir_SpecificType)
        Area        => elemR(eIdx,er_Area)
        Depth       => elemR(eIdx,er_Depth)
        ell         => elemR(eIdx,er_ell)
        Head        => elemR(eIdx,er_Head)
        HydDepth    => elemR(eIdx,er_HydDepth)
        HydRadius   => elemR(eIdx,er_HydRadius)
        Length      => elemR(eIdx,er_Length)
        Perimeter   => elemR(eIdx,er_Perimeter)
        Topwidth    => elemR(eIdx,er_Topwidth)
        Volume      => elemR(eIdx,er_Volume)
        Zbottom     => elemR(eIdx,er_Zbottom)
        CurrentSetting          => elemR(eIdx,er_Setting)
        FullDepth               => elemSR(eIdx,esr_Weir_FullDepth)
        RectangularBreadth      => elemSR(eIdx,esr_Weir_RectangularBreadth)
        TrapezoidalBreadth      => elemSR(eIdx,esr_Weir_TrapezoidalBreadth)
        TriangularSideSlope     => elemSR(eIdx,esr_Weir_TriangularSideSlope)
        TrapezoidalLeftSlope    => elemSR(eIdx,esr_Weir_TrapezoidalLeftSlope)
        TrapezoidalRightSlope   => elemSR(eIdx,esr_Weir_TrapezoidalRightSlope)
        Zcrest                  => elemSR(eIdx,esr_Weir_Zcrest)
        Zcrown                  => elemSR(eIdx,esr_Weir_Zcrown)
        
        IsSurcharged => elemYN(eIdx,eYN_isSurcharged)
        !%-----------------------------------------------------------------------------     
        !% find depth on weir
        if (Head <= Zcrest) then
            Depth = zeroR
        elseif ((Head > Zcrest) .and. (Head < Zcrown)) then
            Depth =  Head - Zcrest
        else
            Depth = Zcrown - Zcrest
        endif

        !% find offset of weir cresr due to control setting
        z  = (oneR - CurrentSetting) * FullDepth
        zY = min(z+Depth,FullDepth) 
        
        !% set geometry variables for weir types
        select case (SpecificWeirType) 
            case (transverse_weir)
                Area      = RectangularBreadth * zY - RectangularBreadth * z
                Volume    = Area * Length  !% HACK this is not the correct volume in the element
                Topwidth  = RectangularBreadth
                HydDepth  = Depth !% HACK this is not the correct hydraulic depth in the element
                ell       = Head - Zbottom
                Perimeter = Topwidth + twoR * HydDepth
                HydRadius = Area / Perimeter
                
            case (side_flow)
                Area      = RectangularBreadth * zY - RectangularBreadth * z
                Volume    = Area * Length
                Topwidth  = RectangularBreadth
                HydDepth  = Depth
                ell       = Head - Zbottom
                Perimeter = Topwidth + twoR * HydDepth
                HydRadius = Area / Perimeter
            
            case (trapezoidal_weir)
                Area      = (TrapezoidalBreadth + onehalfR * (TrapezoidalLeftSlope + TrapezoidalRightSlope) * zY) * zY &
                        - (TrapezoidalBreadth + onehalfR * (TrapezoidalLeftSlope + TrapezoidalRightSlope) * z) * z 
                Volume    = Area * Length
                Topwidth  = TrapezoidalBreadth + Depth &
                            * (TrapezoidalLeftSlope + TrapezoidalRightSlope)
                HydDepth  = Area / Topwidth
                ell       = Head - Zbottom
                Perimeter = TrapezoidalBreadth + Depth &
                                * (sqrt(oneR + (TrapezoidalLeftSlope**twoR)) &
                                + sqrt(oneR + (TrapezoidalRightSlope**twoR)))
                HydRadius = Area / Perimeter
                
            case (vnotch_weir)
                Area      = TriangularSideSlope * zY ** twoR - TriangularSideSlope * z ** twoR
                Volume    = Area * Length
                Topwidth  = twoR * TriangularSideSlope * Depth
                HydDepth  = onehalfR * Depth
                Perimeter = twoR * Depth * sqrt(oneR + (TriangularSideSlope ** twoR))
                HydRadius = (TriangularSideSlope * Depth) &
                                / (twoR * sqrt(oneR + (TriangularSideSlope ** twoR)))
            case default
                print *, 'CODE ERROR: unknown weir type, ', SpecificWeirType,'  in network'
                print *, 'which has key ',trim(reverseKey(SpecificWeirType))
                stop 3358223
        end select

        !% apply geometry limiters
        call adjust_limit_by_zerovalues_singular (eIdx, er_Area,      setting%ZeroValue%Area,    .false.)
        call adjust_limit_by_zerovalues_singular (eIdx, er_Depth,     setting%ZeroValue%Depth,   .false.)
        call adjust_limit_by_zerovalues_singular (eIdx, er_HydDepth,  setting%ZeroValue%Depth,   .false.)
        call adjust_limit_by_zerovalues_singular (eIdx, er_HydRadius, setting%ZeroValue%Depth,   .false.)
        call adjust_limit_by_zerovalues_singular (eIdx, er_ell,       setting%ZeroValue%Depth,   .false.) 
        call adjust_limit_by_zerovalues_singular (eIdx, er_Topwidth,  setting%ZeroValue%Topwidth,.false.)
        call adjust_limit_by_zerovalues_singular (eIdx, er_Perimeter, setting%ZeroValue%Topwidth,.false.)
        call adjust_limit_by_zerovalues_singular (eIdx, er_Volume,    setting%ZeroValue%Volume,  .true.)

    end subroutine weir_geometry_update
!%
!%========================================================================== 
!%==========================================================================    
!%  
    subroutine weir_get_open_area (eIdx)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% specilized subroutine to get the flow area only
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: eIdx
        real(8), pointer :: Area, CurrentSetting, Depth, Head, FullDepth 
        real(8), pointer :: Zbottom,  Zcrown, Zcrest
        real(8), pointer :: RectangularBreadth, TrapezoidalBreadth
        real(8), pointer :: TriangularSideSlope, TrapezoidalLeftSlope, TrapezoidalRightSlope
        integer, pointer :: SpecificWeirType
        real(8)          :: z, zY
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        !% pointers
        SpecificWeirType => elemSI(eIdx,esi_Weir_SpecificType)
        Area           => elemR(eIdx,er_Area)
        Depth          => elemR(eIdx,er_Depth)
        Head           => elemR(eIdx,er_Head)
        CurrentSetting => elemR(eIdx,er_setting)
        Zbottom        => elemR(eIdx,er_Zbottom)
        FullDepth               => elemSR(eIdx,esr_Weir_FullDepth)
        RectangularBreadth      => elemSR(eIdx,esr_Weir_RectangularBreadth)
        TrapezoidalBreadth      => elemSR(eIdx,esr_Weir_TrapezoidalBreadth)
        TriangularSideSlope     => elemSR(eIdx,esr_Weir_TriangularSideSlope)
        TrapezoidalLeftSlope    => elemSR(eIdx,esr_Weir_TrapezoidalLeftSlope)
        TrapezoidalRightSlope   => elemSR(eIdx,esr_Weir_TrapezoidalRightSlope)
        Zcrest                  => elemSR(eIdx,esr_Weir_Zcrest)
        Zcrown                  => elemSR(eIdx,esr_Weir_Zcrown)
        !%-----------------------------------------------------------------------------     
        !% find depth on weir
        if (Head <= Zcrest) then
            Depth = zeroR
        elseif ((Head > Zcrest) .and. (Head < Zcrown)) then
            Depth =  Head - Zcrest
        else
            Depth = Zcrown - Zcrest
        endif
        
        !% find offset of weir cresr due to control setting
        z  = (oneR - CurrentSetting) * FullDepth
        zY = min(z+Depth,FullDepth) 
        
        !% set geometry variables for weir types
        select case (SpecificWeirType) 
            case (transverse_weir)
                Area      = RectangularBreadth * zY - RectangularBreadth * z
            case (side_flow)
                Area      = RectangularBreadth * zY - RectangularBreadth * z
            case (trapezoidal_weir)
                Area      = (TrapezoidalBreadth + onehalfR * (TrapezoidalLeftSlope + TrapezoidalRightSlope) * zY) * zY &
                        - (TrapezoidalBreadth + onehalfR * (TrapezoidalLeftSlope + TrapezoidalRightSlope) * z) * z
            case (vnotch_weir)
                Area      = TriangularSideSlope * zY ** twoR - TriangularSideSlope * z ** twoR
            case default
                print *, 'CODE ERROR: unknown weir type, ', SpecificWeirType,'  in network'
                print *, 'which has key ',trim(reverseKey(SpecificWeirType))
                stop 3358223
        end select

        !% apply geometry limiters
        call adjust_limit_by_zerovalues_singular (eIdx, er_Area, setting%ZeroValue%Area, .false.)

    end subroutine weir_get_open_area
!%
!%========================================================================== 
!%==========================================================================    
!%  
    subroutine weir_get_flapgate_headLoss (eIdx, inCol)
    !%-----------------------------------------------------------------------------
        !% Description:
        !%  computes the aramco headloss due to a flap gate
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: eIdx, inCol
        real(8), pointer :: Area, Flowrate, grav, Velocity, EffectiveHeadDelta, zeroArea
        real(8)          :: hLoss
        !%-----------------------------------------------------------------------------
        Area                  => elemR(eIdx,er_area)
        Flowrate              => elemR(eIdx,er_Flowrate)
        Velocity              => elemR(eIdx,er_Velocity)
        EffectiveHeadDelta    => elemSR(eIdx,inCol)
        zeroArea              => setting%ZeroValue%Area
        grav                  => setting%Constant%gravity
        !%-----------------------------------------------------------------------------

        call weir_get_open_area(eIdx) 

        if (Area > zeroArea) then
            Velocity = Flowrate / Area
            hLoss    = (fourR / grav) * Velocity * Velocity &
                     * exp(-1.15 * Velocity / sqrt(EffectiveHeadDelta))
            EffectiveHeadDelta = max(EffectiveHeadDelta - hLoss, zeroR)
        end if
        
    end subroutine weir_get_flapgate_headLoss
         
!%
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module weir_elements