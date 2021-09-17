module orifice_elements

    use define_globals
    use define_keys
    use define_indexes
    use define_settings, only: setting
    use common_elements
    use adjust


    implicit none

    !%----------------------------------------------------------------------------- 
    !% Description:
    !% Computes diagnostic flow through orifice elements
    !%----------------------------------------------------------------------------- 

    private

    public :: orifice_toplevel

    real(8), pointer :: grav => setting%constant%gravity

    contains
    !%==========================================================================
    !% PUBLIC
    !%==========================================================================
    !%
    subroutine orifice_toplevel (eIdx)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% We need a subroutine to get the new full depth (eSr_EffectiveFullDepth) 
        !% crown (eSr_Zcrown) and crest (eSr_Zcrest) elevation from control setting.
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: eIdx  !% must be a single element ID
        !%-----------------------------------------------------------------------------
        !%  
        call common_head_and_flowdirection_singular &
            (eIdx, eSr_Orifice_Zcrest, eSr_Orifice_NominalDownstreamHead, eSi_Orifice_FlowDirection)
        
        !% find effective head on orifice element
         call orifice_effective_head_delta (eIdx)
        
        !% find flow on orifice element
        call orifice_flow (eIdx)
        
        !% update orifice geometry from head
        call orifice_geometry_update (eIdx)
        
         !% update velocity from flowrate and area
        call common_velocity_from_flowrate_singular (eIdx)

    end subroutine orifice_toplevel
    ! %
    !%==========================================================================
    !% PRIVATE
    !%========================================================================== 
    !%  
    subroutine orifice_effective_head_delta (eIdx)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% 
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: eIdx !% single ID of element
        real(8), pointer :: EffectiveHeadDelta, NominalDownstreamHead, Head
        real(8), pointer :: Zcrown, Zcrest
        integer, pointer :: SpecificOrificeType
        real(8) :: Zmidpt
        !%-----------------------------------------------------------------------------
        !% inputs
        SpecificOrificeType   => elemSI(eIdx,eSi_specific_orifice_type)
        Head                  => elemR(eIdx,er_Head)
        Zcrown                => elemR(eIdx,er_Zcrown)
        Zcrest                => elemSR(eIdx,eSr_Orifice_Zcrest)
        NominalDownstreamHead => elemSR(eIdx,eSr_Orifice_NominalDownstreamHead)
        !% output
        EffectiveHeadDelta         => elemSR(eIdx,eSr_Orifice_EffectiveHeadDelta)
        !%-----------------------------------------------------------------------------
        select case (SpecificOrificeType)
            case (BOTTOM_ORIFICE)
                if (Head <= Zcrest) then
                    EffectiveHeadDelta = zeroR
                elseif (NominalDownstreamHead > Zcrest) then
                    EffectiveHeadDelta = Head - NominalDownstreamHead
                else
                    EffectiveHeadDelta = Head - Zcrest
                end if
            case (SIDE_ORIFICE)
                if (Head <= Zcrest) then
                    EffectiveHeadDelta = zeroR
                elseif (Head < Zcrown) then
                    EffectiveHeadDelta = Head - Zcrest
                else
                    Zmidpt = (Zcrown + Zcrest)/2.0
                    if (NominalDownstreamHead < Zmidpt) then
                        EffectiveHeadDelta = Head - Zmidpt
                    else
                        EffectiveHeadDelta = Head - NominalDownstreamHead
                    end if
                end if
        end select
                    
    end subroutine orifice_effective_head_delta
    !%
    !%========================================================================== 
    !%========================================================================== 
    !%   
    subroutine orifice_flow (eIdx)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% 
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: eIdx
        integer, pointer :: SpecificOrificeType, FlowDirection, GeometryType
        real(8), pointer :: Flowrate, EffectiveHeadDelta, Zcrest, Head
        real(8), pointer :: RectangularBreadth, NominalDownstreamHead
        real(8), pointer :: DischargeCoeff, EffectiveFullDepth
        real(8), pointer :: WeirExponent, VillemonteExponent, SharpCrestedWeirCoeff
        real(8) :: CriticalDepth, AoverL, FractionCritDepth, Coef, FullArea
        real(8) :: ratio
        !%-----------------------------------------------------------------------------
        GeometryType          => elemI(eIdx,ei_geometryType)
        SpecificOrificeType   => elemSI(eIdx,eSi_specific_orifice_type)
        FlowDirection         => elemSI(eIdx,eSi_Orifice_FlowDirection)
        Flowrate              => elemR(eIdx,er_Flowrate)
        Head                  => elemR(eIdx,er_Head)
        EffectiveHeadDelta    => elemSR(eIdx,eSr_Orifice_EffectiveHeadDelta)
        Zcrest                => elemSR(eIdx,eSr_Orifice_Zcrest)
        RectangularBreadth    => elemSR(eIdx,eSr_Orifice_RectangularBreadth)
        DischargeCoeff        => elemSR(eIdx,eSr_Orifice_DischargeCoeff)
        EffectiveFullDepth    => elemSR(eIdx,eSr_Orifice_EffectiveFullDepth)
        NominalDownstreamHead => elemSR(eIdx,eSr_Orifice_NominalDownstreamHead)
        
        SharpCrestedWeirCoeff => Setting%Orifice%SharpCrestedWeirCoefficient
        WeirExponent          => Setting%Orifice%TransverseWeirExponent
        VillemonteExponent    => Setting%Orifice%VillemonteCorrectionExponent
        !%-----------------------------------------------------------------------------
        !% find full area for flow, and A/L for critical depth calculations
        select case (GeometryType)
            case (circular)
                FullArea = pi * (onehalfR * EffectiveFullDepth) ** twoR
                AoverL   = onefourthR * EffectiveFullDepth 
            case (rectangular)
                FullArea = EffectiveFullDepth * RectangularBreadth
                AoverL   = FullArea / (twoR * (EffectiveFullDepth + RectangularBreadth))
            case default
                print *, 'error, case default should not be reached.'
                stop 5983
        end select
        
        !% find critical depth to determine weir/orifice flow
        select case (SpecificOrificeType)
            case (BOTTOM_ORIFICE)
                !% find critical height above opening where orifice flow turns into 
                !% weir flow for Bottom orifice = (C_orifice/C_weir)*(Area/Length)
                !% where C_orifice = given orifice coeff, C_weir = weir_coeff/sqrt(2g),
                !% Area is the area of the opening, and Length = circumference
                !% of the opening. For a basic sharp crested weir, C_weir = 0.414.
                CriticalDepth = DischargeCoeff / SharpCrestedWeirCoeff * AoverL
                FractionCritDepth = min(EffectiveHeadDelta / CriticalDepth, oneR)
            case (SIDE_ORIFICE)
                CriticalDepth = EffectiveFullDepth
                FractionCritDepth = min(((Head - Zcrest) / CriticalDepth), oneR)
                !% another adjustment to critical depth is needed
                !% for weir coeff calculation for side orifice
                CriticalDepth = onehalfR * CriticalDepth 
            end select
        
        !% flow calculation conditions through an orifice
        if ((EffectiveHeadDelta == zeroR) .or. (FractionCritDepth <= zeroR)) then
            !% no flow case
            Flowrate = zeroR
        elseif (FractionCritDepth < oneR) then
            !% case where inlet depth is below critical depth thus,
            !% orifice behaves as a rectangular transverse weir
            Coef     = DischargeCoeff * sqrt(CriticalDepth)
            Flowrate = FlowDirection * Coef * (FractionCritDepth ** WeirExponent)
        else 
            !% standard orifice flow condition 
            Coef      = DischargeCoeff * FullArea * sqrt(twoR * grav)
            Flowrate  = FlowDirection * Coef * sqrt(EffectiveHeadDelta)       
        end if    
        
        !% applying Villemonte submergence correction for orifice having submerged weir flow
        if ((FractionCritDepth < oneR) .and. (NominalDownstreamHead > Zcrest)) then
            ratio = (NominalDownstreamHead - Zcrest) / (Head - Zcrest)
            Flowrate = Flowrate * ((oneR - (ratio ** WeirExponent)) ** VillemonteExponent) 
        end if
            
    end subroutine orifice_flow
    !%
    !%==========================================================================
    !%==========================================================================    
    !%  
    subroutine orifice_geometry_update (eIdx)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% HACK -- it is not clear as yet what geometries we actually need. There's
        !% an important difference between the geometry of the flow thru the orifice
        !% and the geometry surrounding the orifice.
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: eIdx
        real(8), pointer :: Head, Length, Zbottom,  Zcrown
        real(8), pointer :: Depth, Area, Volume, Topwidth
        real(8), pointer :: Perimeter, HydDepth, HydRadius,  Zcrest
        real(8), pointer :: RectangularBreadth, TrapezoidalBreadth
        real(8), pointer :: TriangularSideSlope, TrapezoidalLeftSlope, TrapezoidalRightSlope
        integer, pointer :: GeometryType
        !%-----------------------------------------------------------------------------
        GeometryType => elemI(eIdx,ei_geometryType)
        Head        => elemR(eIdx,er_Head)
        Length      => elemR(eIdx,er_Length)
        Zbottom     => elemR(eIdx,er_Zbottom)
        Zcrown      => elemR(eIdx,er_Zcrown)
        Depth       => elemR(eIdx,er_Depth)
        Area        => elemR(eIdx,er_Area)
        Volume      => elemR(eIdx,er_Volume)
        Topwidth    => elemR(eIdx,er_Topwidth)
        Perimeter   => elemR(eIdx,er_Perimeter)
        HydDepth    => elemR(eIdx,er_HydDepth)
        HydRadius   => elemR(eIdx,er_HydRadius)
        Zcrest                  => elemSR(eIdx,eSr_Orifice_Zcrest)
        RectangularBreadth      => elemSR(eIdx,eSr_Orifice_RectangularBreadth)
    
        !% find depth over bottom of orifice
        if (Head <= Zcrest) then
            Depth = zeroR
        elseif ((Head > Zcrest) .and. (Head < Zcrown)) then
            Depth =  Head - Zcrest
        else
            Depth = Zcrown - Zcrest
        end if
    
        !% set geometry
        select case (GeometryType)
            case (rectangular)
                Area      =  RectangularBreadth * Depth
                Volume    = Area * Length !% HACK this is not the correct volume in the element
                Topwidth  = RectangularBreadth
                HydDepth  = Depth !% HACK this is not the correct hydraulic depth in the element
                Perimeter = Topwidth + twoR * HydDepth
                HydRadius = Area / Perimeter
            case (circular)
                print *, 'error, the circular orifice is not yet implemented'
                stop 2087
            case default
                print *, 'error, the default case should not be reached'
                stop 9478
        end select
        
        !% apply geometry limiters
        call adjust_limit_by_zerovalues_singular (eIdx, er_Area,      setting%ZeroValue%Area)
        call adjust_limit_by_zerovalues_singular (eIdx, er_Depth,     setting%ZeroValue%Depth)
        call adjust_limit_by_zerovalues_singular (eIdx, er_HydDepth,  setting%ZeroValue%Depth)
        call adjust_limit_by_zerovalues_singular (eIdx, er_HydRadius, setting%ZeroValue%Depth)
        call adjust_limit_by_zerovalues_singular (eIdx, er_Topwidth,  setting%ZeroValue%Topwidth)
        call adjust_limit_by_zerovalues_singular (eIdx, er_Perimeter, setting%ZeroValue%Topwidth)
        call adjust_limit_by_zerovalues_singular (eIdx, er_Volume,    setting%ZeroValue%Volume)
    
    end subroutine  orifice_geometry_update
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
end module orifice_elements