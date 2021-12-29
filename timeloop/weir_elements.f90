module weir_elements

    use define_globals
    use define_keys
    use define_indexes
    use define_settings, only: setting
    use common_elements
    use adjust

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% Computes flow and head for weir elements
    !%----------------------------------------------------------------------------- 

    implicit none

    private

    public :: weir_toplevel

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
        if (icrash) return
        if (setting%Debug%File%weir_elements) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------
        isSurcharged => elemYN(eIdx,eYN_isSurcharged)
        !%  
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
        
        !% update weir geometry from head
        call weir_geometry_update (eIdx)
        
        !% update velocity from flowrate and area
        call common_velocity_from_flowrate_singular (eIdx)
        
        if (setting%Debug%File%weir_elements)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine weir_toplevel    
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
        real(8), pointer :: NominalDownstreamHead
        logical, pointer :: CanSurcharge, IsSurcharged
        real(8) :: Zmidpt
        !%-----------------------------------------------------------------------------
        if (icrash) return
        !% input
        Head   => elemR(eIdx,er_Head)
        !% output
        EffectiveHeadDelta    => elemSR(eIdx,esr_Weir_EffectiveHeadDelta)
        Zcrown                => elemSR(eIdx,esr_Weir_Zcrown)
        Zcrest                => elemSR(eIdx,esr_Weir_Zcrest)
        NominalDownstreamHead => elemSR(eIdx,esr_Weir_NominalDownstreamHead)
        CanSurcharge          => elemYN(eIdx,eYN_canSurcharge)
        IsSurcharged          => elemYN(eIdx,eYN_isSurcharged)
        
        !% setting default surcharge condition as false
        IsSurcharged = .false.
        !%-----------------------------------------------------------------------------

        if (Head <= Zcrest) then
            EffectiveHeadDelta = zeroR
        else
            EffectiveHeadDelta = Head - Zcrest
        endif
            
        if ((Head > Zcrown) .and. (CanSurcharge)) then
            IsSurcharged = .true.
            Zmidpt = (Zcrest + Zcrown) / twoR
            if (NominalDownstreamHead < Zmidpt) then
                EffectiveHeadDelta = Head - Zmidpt       
            else
                EffectiveHeadDelta = Head - NominalDownstreamHead    
            endif     
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
        if (icrash) return
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
        if (icrash) return
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
        real(8), pointer :: Perimeter, HydDepth,  Zcrest
        real(8), pointer :: RectangularBreadth, TrapezoidalBreadth
        real(8), pointer :: TriangularSideSlope, TrapezoidalLeftSlope, TrapezoidalRightSlope
        integer, pointer :: SpecificWeirType
        logical, pointer :: IsSurcharged
        !%-----------------------------------------------------------------------------
        if (icrash) return
        SpecificWeirType => elemSI(eIdx,esi_Weir_SpecificType)

        Head        => elemR(eIdx,er_Head)
        Length      => elemR(eIdx,er_Length)
        Zbottom     => elemR(eIdx,er_Zbottom)
        Depth       => elemR(eIdx,er_Depth)
        Area        => elemR(eIdx,er_Area)
        Volume      => elemR(eIdx,er_Volume)
        Topwidth    => elemR(eIdx,er_Topwidth)
        Perimeter   => elemR(eIdx,er_Perimeter)
        HydDepth    => elemR(eIdx,er_HydDepth)
        HydRadius   => elemR(eIdx,er_HydRadius)
        Zcrest                  => elemSR(eIdx,esr_Weir_Zcrest)
        Zcrown                  => elemSR(eIdx,esr_Weir_Zcrown)
        RectangularBreadth      => elemSR(eIdx,esr_Weir_RectangularBreadth)
        TrapezoidalBreadth      => elemSR(eIdx,esr_Weir_TrapezoidalBreadth)
        TriangularSideSlope     => elemSR(eIdx,esr_Weir_TriangularSideSlope)
        TrapezoidalLeftSlope    => elemSR(eIdx,esr_Weir_TrapezoidalLeftSlope)
        TrapezoidalRightSlope   => elemSR(eIdx,esr_Weir_TrapezoidalRightSlope)
        
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
                Perimeter = Topwidth + twoR * HydDepth
                HydRadius = Area / Perimeter
                
            case (side_flow)
                Area      = RectangularBreadth * Depth
                Volume    = Area * Length
                Topwidth  = RectangularBreadth
                HydDepth  = Depth
                Perimeter = Topwidth + twoR * HydDepth
                HydRadius = Area / Perimeter
            
            case (trapezoidal_weir)
                Area      =  (TrapezoidalBreadth + onehalfR * &
                             (TrapezoidalLeftSlope + TrapezoidalRightSlope) * Depth) * Depth 
                Volume    = Area * Length
                Topwidth  = TrapezoidalBreadth + Depth &
                            * (TrapezoidalLeftSlope + TrapezoidalRightSlope)
                HydDepth  = Area / Topwidth
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
                
        end select

        !% apply geometry limiters
        call adjust_limit_by_zerovalues_singular (eIdx, er_Area,      setting%ZeroValue%Area)
        call adjust_limit_by_zerovalues_singular (eIdx, er_Depth,     setting%ZeroValue%Depth)
        call adjust_limit_by_zerovalues_singular (eIdx, er_HydDepth,  setting%ZeroValue%Depth)
        call adjust_limit_by_zerovalues_singular (eIdx, er_HydRadius, setting%ZeroValue%Depth)
        call adjust_limit_by_zerovalues_singular (eIdx, er_Topwidth,  setting%ZeroValue%Topwidth)
        call adjust_limit_by_zerovalues_singular (eIdx, er_Perimeter, setting%ZeroValue%Topwidth)
        call adjust_limit_by_zerovalues_singular (eIdx, er_Volume,    setting%ZeroValue%Volume)

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