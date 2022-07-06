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
            call weir_flow (eIdx, esr_Weir_EffectiveHeadDelta, .true.)
        endif
        
        !print *, 'in weir toplevel'
        !call util_CLprint()
        !stop 98374
        
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
        real(8), pointer :: EffectiveHeadDelta, Head, Zcrown, Zcrest
        real(8), pointer :: NominalDownstreamHead, CurrentSetting, EffectiveFullDepth
        logical, pointer :: CanSurcharge, IsSurcharged
        real(8) :: Zmidpt
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        !% input
        Head   => elemR(eIdx,er_Head)
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

        if (Head <= Zcrest) then
            EffectiveHeadDelta = zeroR
        else
            EffectiveHeadDelta = Head - Zcrest
        endif
            
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
        endif

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
        real(8), pointer :: Flowrate, EffectiveFullDepth, EffectiveHeadDelta
        real(8) :: CoeffOrifice
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        FlowDirection      => elemSI(eIdx,esi_Weir_FlowDirection)
        Flowrate           => elemR(eIdx,er_Flowrate) 
        EffectiveFullDepth => elemSR(eIdx,esr_Weir_EffectiveFullDepth)
        EffectiveHeadDelta => elemSR(eIdx,esr_Weir_EffectiveHeadDelta)
        !%-----------------------------------------------------------------------------
        ! get the flowrate for effective full depth without submergence correction
        call weir_flow(eIdx, esr_Weir_EffectiveFullDepth, .false.)

        ! equivalent orifice flow coefficient for surcharge flow
        CoeffOrifice = Flowrate / sqrt(EffectiveFullDepth/twoR)
        
        !% old flowrate is overwritten by new surcharged flowrate
        Flowrate = FlowDirection * CoeffOrifice * sqrt(EffectiveHeadDelta)

    end subroutine weir_surcharge_flow
!%
!%========================================================================== 
!%==========================================================================    
!%  
    subroutine weir_flow (eIdx, inCol, ApplySubmergenceCorrection)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Computes flow for standard weir types
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: eIdx, inCol
        logical, intent(in) :: ApplySubmergenceCorrection
        integer, pointer :: SpecificWeirType, EndContractions, FlowDirection
        real(8), pointer :: Flowrate, Head, EffectiveHeadDelta
        real(8), pointer :: RectangularBreadth, TrapezoidalBreadth
        real(8), pointer :: TriangularSideSlope, TrapezoidalLeftSlope, TrapezoidalRightSlope
        real(8), pointer :: CoeffTriangular, CoeffRectangular
        real(8), pointer :: WeirExponent, WeirExponentVNotch
        real(8), pointer :: WeirContractionFactor, VillemonteExponent, WeirCrestExponent
        real(8), pointer :: NominalDsHead, Zcrest
        real(8) :: Zmidpt, CrestLength, SubCorrectionTriangular, SubCorrectionRectangular
        real(8) :: ratio
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        SpecificWeirType => elemSI(eIdx,esi_Weir_SpecificType)
        EndContractions  => elemSI(eIdx,esi_Weir_EndContractions)
        FlowDirection    => elemSI(eIdx,esi_Weir_FlowDirection)

        Head                  => elemR(eIdx,er_Head)
        Flowrate              => elemR(eIdx,er_Flowrate)
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
        !%-----------------------------------------------------------------------------
        !% initializing default local Villemonte submergence correction factors as 1
        !% These are changed below if needed
        SubCorrectionTriangular = oneR
        SubCorrectionRectangular = oneR

        !print *, reverseKey(SpecificWeirType)
        !call util_CLprint()
        
    
        select case (SpecificWeirType)
        case (transverse_weir)
            WeirExponent          => Setting%Weir%Transverse%WeirExponent
            WeirContractionFactor => Setting%Weir%Transverse%WeirContractionFactor
            VillemonteExponent    => Setting%Weir%Transverse%VillemonteCorrectionExponent

            !% effective crest length due to contraction for tranverse weir             
            CrestLength = max(zeroR, &
                    RectangularBreadth - WeirContractionFactor * real(EndContractions,8) * EffectiveHeadDelta) 

            !% correction factor for nominal downstream submergence
            if ((NominalDsHead > Zcrest) .and. (ApplySubmergenceCorrection)) then
                ratio = (NominalDsHead - Zcrest) / (Head - Zcrest)        
                SubCorrectionRectangular = ((oneR - (ratio ** WeirExponent)) ** VillemonteExponent)
            endif

            Flowrate = real(FlowDirection,8) * SubCorrectionRectangular * CrestLength * &
                    CoeffRectangular  * (EffectiveHeadDelta ** WeirExponent)

            !print *, EffectiveHeadDelta        

        case (side_flow)
            WeirExponent          => Setting%Weir%SideFlow%WeirExponent
            WeirContractionFactor => Setting%Weir%SideFlow%WeirContractionFactor
            WeirCrestExponent     => Setting%Weir%SideFlow%SideFlowWeirCrestExponent
            VillemonteExponent    => Setting%Weir%SideFlow%VillemonteCorrectionExponent

            !% effective crest length due to contraction for sideflow weir   
            CrestLength = max(zeroR, &
                    RectangularBreadth - WeirContractionFactor * real(EndContractions,8) * EffectiveHeadDelta)
                    
            if (FlowDirection > zeroR) then
                
                !% correction factor for nominal downstream submergence
                if ((NominalDsHead > Zcrest) .and. (ApplySubmergenceCorrection)) then
                    ratio = (NominalDsHead - Zcrest) / (Head - Zcrest) 
                    SubCorrectionRectangular = ((oneR - (ratio ** WeirExponent)) ** VillemonteExponent)
                endif
            
                Flowrate = real(FlowDirection,8)  * SubCorrectionRectangular * (CrestLength ** &
                    WeirCrestExponent) * CoeffRectangular * (EffectiveHeadDelta ** WeirExponent)
            
            else
                !% under reverse flow condition, sideflow weir behaves like a transverse weir
                !% correction factor for nominal downstream submergence
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

            !% effective crest length due for trapezoidal weir
            !% HACK: the crest length changes if there is control present
            CrestLength = TrapezoidalBreadth
            
            !% correction factor for nominal downstream submergence
            if ((NominalDsHead > Zcrest) .and. (ApplySubmergenceCorrection)) then
                ratio = (NominalDsHead - Zcrest) / (Head - Zcrest)   
                SubCorrectionRectangular = ((oneR - (ratio ** WeirExponent)) **  VillemonteExponent)
                SubCorrectionTriangular  = ((oneR - (ratio ** WeirExponentVNotch)) ** VillemonteExponent)
            endif

            Flowrate = real(FlowDirection,8) * &
                    (SubCorrectionTriangular * CoeffTriangular * ((TrapezoidalLeftSlope + &
                    TrapezoidalRightSlope) / twoR) * (EffectiveHeadDelta ** WeirExponentVNotch) &
                    + &
                    SubCorrectionRectangular * CoeffRectangular * CrestLength * &
                    (EffectiveHeadDelta ** WeirExponent))
                            
        case (vnotch_weir)
            WeirExponent          => Setting%Weir%VNotch%WeirExponent
            WeirContractionFactor => Setting%Weir%VNotch%WeirContractionFactor
            WeirCrestExponent     => Setting%Weir%VNotch%SideFlowWeirCrestExponent
            VillemonteExponent    => Setting%Weir%VNotch%VillemonteCorrectionExponent

            !% correction factor for nominal downstream submergence
            if ((NominalDsHead > Zcrest) .and. (ApplySubmergenceCorrection)) then
                ratio = (NominalDsHead - Zcrest) / (Head - Zcrest)
                SubCorrectionTriangular = ((oneR - (ratio ** WeirExponent)) ** VillemonteExponent)
            endif
            
            Flowrate = real(FlowDirection,8) * SubCorrectionTriangular * CoeffTriangular * &
                    TriangularSideSlope * (EffectiveHeadDelta ** WeirExponent) 
            
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
        real(8), pointer :: Head, Length, Zbottom,  Zcrown
        real(8), pointer :: Depth, Area, Volume, Topwidth, HydRadius
        real(8), pointer :: Perimeter, HydDepth,  Zcrest, ell
        real(8), pointer :: RectangularBreadth, TrapezoidalBreadth
        real(8), pointer :: TriangularSideSlope, TrapezoidalLeftSlope, TrapezoidalRightSlope
        integer, pointer :: SpecificWeirType
        logical, pointer :: IsSurcharged
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
        
        !% set geometry variables for weir types
        select case (SpecificWeirType) 
        case (transverse_weir)
            Area      = RectangularBreadth * Depth
            Volume    = Area * Length  !% HACK this is not the correct volume in the element
            Topwidth  = RectangularBreadth
            HydDepth  = Depth !% HACK this is not the correct hydraulic depth in the element
            ell       = Head - Zbottom
            Perimeter = Topwidth + twoR * HydDepth
            HydRadius = Area / Perimeter
            
        case (side_flow)
            Area      = RectangularBreadth * Depth
            Volume    = Area * Length
            Topwidth  = RectangularBreadth
            HydDepth  = Depth
            ell       = Head - Zbottom
            Perimeter = Topwidth + twoR * HydDepth
            HydRadius = Area / Perimeter
        
        case (trapezoidal_weir)
            Area      =  (TrapezoidalBreadth + onehalfR * &
                            (TrapezoidalLeftSlope + TrapezoidalRightSlope) * Depth) * Depth 
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
            Area      =  TriangularSideSlope * Depth ** twoR
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
        !%-----------------------------------------------------------------------------
        !% Description:
        !% 
        !%-----------------------------------------------------------------------------

        !%-----------------------------------------------------------------------------
        !%  
!%
!%========================================================================== 
!%==========================================================================    
!%  
        !%-----------------------------------------------------------------------------
        !% Description:
        !% 
        !%-----------------------------------------------------------------------------

        !%-----------------------------------------------------------------------------
        !%  
!%
!%==========================================================================
!% END OF MODULE
!%+=========================================================================
end module weir_elements