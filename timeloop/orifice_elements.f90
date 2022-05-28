module orifice_elements

    use define_globals
    use define_keys
    use define_indexes
    use define_settings, only: setting
    use common_elements
    use adjust
    use define_xsect_tables
    use xsect_tables


    implicit none

    !%-----------------------------------------------------------------------------
    !% Description:
    !% Computes diagnostic flow through orifice elements
    !%-----------------------------------------------------------------------------

    private

    public :: orifice_toplevel

    contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine orifice_toplevel (eIdx)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% We need a subroutine to get the new full depth (esr_EffectiveFullDepth)
        !% crown (esr_Zcrown) and crest (esr_Zcrest) elevation from control setting.
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: eIdx  !% must be a single element ID

        character(64) :: subroutine_name = 'orifice_toplevel'
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%orifice_elements) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        call common_head_and_flowdirection_singular &
            (eIdx, esr_Orifice_Zcrest, esr_Orifice_NominalDownstreamHead, esi_Orifice_FlowDirection)

        !% find effective head on orifice element
         call orifice_effective_head_delta (eIdx)

        !% find flow on orifice element
        call orifice_flow (eIdx)

        !% update orifice geometry from head
        call orifice_geometry_update (eIdx)

         !% update velocity from flowrate and area
        call common_velocity_from_flowrate_singular (eIdx)

        if (setting%Debug%File%orifice_elements)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
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
        real(8), pointer    :: EffectiveHeadDelta, NominalDownstreamHead, Head
        real(8), pointer    :: Zcrown, Zcrest
        integer, pointer    :: SpecificOrificeType
        real(8)             :: Zmidpt

        character(64) :: subroutine_name = 'orifice_effective_head_delta'
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%orifice_elements) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% inputs
        SpecificOrificeType   => elemSI(eIdx,esi_Orifice_SpecificType)
        Head                  => elemR(eIdx,er_Head)
        Zcrown                => elemSR(eIdx,esr_Orifice_Zcrown)
        Zcrest                => elemSR(eIdx,esr_Orifice_Zcrest)
        NominalDownstreamHead => elemSR(eIdx,esr_Orifice_NominalDownstreamHead)
        !% output
        EffectiveHeadDelta    => elemSR(eIdx,esr_Orifice_EffectiveHeadDelta)
        !%-----------------------------------------------------------------------------
        select case (SpecificOrificeType)
        case (bottom_orifice)
            if (Head <= Zcrest) then
                EffectiveHeadDelta = zeroR
            elseif (NominalDownstreamHead > Zcrest) then
                EffectiveHeadDelta = Head - NominalDownstreamHead
            else
                EffectiveHeadDelta = Head - Zcrest
            end if
        case (side_orifice)
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
        case default
            print *, 'In ', subroutine_name
            print *, 'CODE ERROR: unknown orifice type, ', SpecificOrificeType,'  in network'
            print *, 'which has key ',trim(reverseKey(SpecificOrificeType))
            stop 862295
        end select

        if (setting%Debug%File%orifice_elements)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine orifice_effective_head_delta
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine orifice_flow (eIdx)
        !%-----------------------------------------------------------------------------
        !% Description: calculates the flow in an orifice elements
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: eIdx
        integer, pointer :: SpecificOrificeType, FlowDirection, GeometryType
        real(8), pointer :: Flowrate, EffectiveHeadDelta, Zcrest, Head, grav
        real(8), pointer :: RectangularBreadth, NominalDownstreamHead
        real(8), pointer :: DischargeCoeff, EffectiveFullDepth, FullArea
        real(8), pointer :: WeirExponent, VillemonteExponent, SharpCrestedWeirCoeff
        real(8) :: CriticalDepth, AoverL, FractionCritDepth, Coef
        real(8) :: ratio

        character(64) :: subroutine_name = 'orifice_flow'
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%orifice_elements) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        GeometryType          => elemI(eIdx,ei_geometryType)
        SpecificOrificeType   => elemSI(eIdx,esi_Orifice_SpecificType)
        FlowDirection         => elemSI(eIdx,esi_Orifice_FlowDirection)
        Flowrate              => elemR(eIdx,er_Flowrate)
        fullArea              => elemR(eIdx,er_FullArea)
        Head                  => elemR(eIdx,er_Head)
        EffectiveHeadDelta    => elemSR(eIdx,esr_Orifice_EffectiveHeadDelta)
        Zcrest                => elemSR(eIdx,esr_Orifice_Zcrest)
        RectangularBreadth    => elemSR(eIdx,esr_Orifice_RectangularBreadth)
        DischargeCoeff        => elemSR(eIdx,esr_Orifice_DischargeCoeff)
        EffectiveFullDepth    => elemSR(eIdx,esr_Orifice_EffectiveFullDepth)
        NominalDownstreamHead => elemSR(eIdx,esr_Orifice_NominalDownstreamHead)
        SharpCrestedWeirCoeff => Setting%Orifice%SharpCrestedWeirCoefficient
        WeirExponent          => Setting%Orifice%TransverseWeirExponent
        VillemonteExponent    => Setting%Orifice%VillemonteCorrectionExponent
        grav                  => setting%constant%gravity
        !%-----------------------------------------------------------------------------

        !% HACK: Hardcoded for Trajkovic cases
        ! if (eIdx == 146) then
        !     if ((setting%Time%Now .ge. 120.00) .and. (setting%Time%Now .lt. 150.00)) then
        !         EffectiveFullDepth = 0.00001
        !     else if (setting%Time%Now .ge. 150.00) then
        !         EffectiveFullDepth = 0.028
        !     end if
        ! end if


        !% find full area for flow, and A/L for critical depth calculations
        select case (GeometryType)
            case (circular)
                fullArea = pi * (onehalfR * EffectiveFullDepth) ** twoR
                AoverL   = onefourthR * EffectiveFullDepth
            case (rectangular_closed)
                FullArea = EffectiveFullDepth * RectangularBreadth
                AoverL   = FullArea / (twoR * (EffectiveFullDepth + RectangularBreadth))
            case default
                print *, 'element idx = ',eIdx
                print *, 'SpecificOrificeType = ',SpecificOrificeType, ' ',reverseKey(SpecificOrificeType)
                print *, 'CODE ERROR geometry type unknown for # ', GeometryType
                print *, 'which has key ',trim(reverseKey(GeometryType))
                stop 5983
        end select

        !% find critical depth to determine weir/orifice flow
        select case (SpecificOrificeType)
        case (bottom_orifice)
            !% find critical height above opening where orifice flow turns into
            !% weir flow for Bottom orifice = (C_orifice/C_weir)*(Area/Length)
            !% where C_orifice = given orifice coeff, C_weir = weir_coeff/sqrt(2g),
            !% Area is the area of the opening, and Length = circumference
            !% of the opening. For a basic sharp crested weir, C_weir = 0.414.
            CriticalDepth = DischargeCoeff / SharpCrestedWeirCoeff * AoverL
            FractionCritDepth = min(EffectiveHeadDelta / CriticalDepth, oneR)
        case (side_orifice)
            CriticalDepth = EffectiveFullDepth
            FractionCritDepth = min(((Head - Zcrest) / EffectiveFullDepth), oneR)
            !% another adjustment to critical depth is needed
            !% for weir coeff calculation for side orifice
            CriticalDepth = onehalfR * CriticalDepth
        case default
            print *, 'In ', subroutine_name
            print *, 'CODE ERROR: unknown orifice type, ', SpecificOrificeType,'  in network'
            print *, 'which has key ',trim(reverseKey(SpecificOrificeType))
            stop 8863411
        end select

        !% flow calculation conditions through an orifice
        if ((EffectiveHeadDelta == zeroR) .or. (FractionCritDepth <= zeroR)) then
            !% no flow case
            Flowrate = zeroR
        elseif (FractionCritDepth < oneR) then
            !% case where inlet depth is below critical depth thus,
            !% orifice behaves as a rectangular transverse weir
            Coef     = DischargeCoeff * FullArea * sqrt(twoR * grav * CriticalDepth)
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
        else
            !% no correction needed
        end if

        if (setting%Debug%File%orifice_elements) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
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
        real(8), pointer :: Perimeter, HydDepth, HydRadius,  Zcrest, Fullarea
        real(8), pointer :: RectangularBreadth, EffectiveFullDepth
        integer, pointer :: GeometryType
        real(8)          :: YoverYfull

        character(64) :: subroutine_name = 'orifice_geometry_update'
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%orifice_elements) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        GeometryType       => elemI(eIdx,ei_geometryType)
        Head               => elemR(eIdx,er_Head)
        Length             => elemR(eIdx,er_Length)
        Zbottom            => elemR(eIdx,er_Zbottom)
        Depth              => elemR(eIdx,er_Depth)
        Area               => elemR(eIdx,er_Area)
        Volume             => elemR(eIdx,er_Volume)
        Topwidth           => elemR(eIdx,er_Topwidth)
        Perimeter          => elemR(eIdx,er_Perimeter)
        HydDepth           => elemR(eIdx,er_HydDepth)
        HydRadius          => elemR(eIdx,er_HydRadius)
        FullArea           => elemR(eIdx,er_FullArea)
        Zcrown             => elemSR(eIdx,esr_Orifice_Zcrown)
        Zcrest             => elemSR(eIdx,esr_Orifice_Zcrest)
        RectangularBreadth => elemSR(eIdx,esr_Orifice_RectangularBreadth)
        EffectiveFullDepth => elemSR(eIdx,esr_Orifice_EffectiveFullDepth)

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
            case (rectangular_closed)
                Area      =  RectangularBreadth * Depth
                Volume    = Area * Length !% HACK this is not the correct volume in the element
                Topwidth  = RectangularBreadth
                HydDepth  = Depth !% HACK this is not the correct hydraulic depth in the element
                Perimeter = Topwidth + twoR * HydDepth
                HydRadius = Area / Perimeter
            case (circular)
                YoverYfull  = Depth / EffectiveFullDepth
                Area        = FullArea * &
                        xsect_table_lookup_singular (YoverYfull, ACirc)  !% 20220506brh removed NACirc
                Volume      = Area * Length
                Topwidth    = EffectiveFullDepth * &
                        xsect_table_lookup_singular (YoverYfull, TCirc) !% 20220506brh removed NTCirc
                HydDepth    = min(Area / Topwidth, EffectiveFullDepth)
                hydRadius   = onefourthR * EffectiveFullDepth * &
                        xsect_table_lookup_singular (YoverYfull, RCirc)  !% 20220506brh removed NRCirc
                Perimeter   = min(Area / hydRadius, &
                        FullArea / (onefourthR * EffectiveFullDepth))
            case default
                print *, 'CODE ERROR geometry type unknown for # ', GeometryType
                print *, 'which has key ',trim(reverseKey(GeometryType))
                stop 9478
        end select

        !% apply geometry limiters
        call adjust_limit_by_zerovalues_singular (eIdx, er_Area,      setting%ZeroValue%Area,     .false.)
        call adjust_limit_by_zerovalues_singular (eIdx, er_Depth,     setting%ZeroValue%Depth,    .false.)
        call adjust_limit_by_zerovalues_singular (eIdx, er_HydDepth,  setting%ZeroValue%Depth,    .false.)
        call adjust_limit_by_zerovalues_singular (eIdx, er_HydRadius, setting%ZeroValue%Depth,    .false.)
        call adjust_limit_by_zerovalues_singular (eIdx, er_Topwidth,  setting%ZeroValue%Topwidth, .false.)
        call adjust_limit_by_zerovalues_singular (eIdx, er_Perimeter, setting%ZeroValue%Topwidth, .false.)
        call adjust_limit_by_zerovalues_singular (eIdx, er_Volume,    setting%ZeroValue%Volume,   .true.)

        if (setting%Debug%File%orifice_elements) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine  orifice_geometry_update
!%==========================================================================
!% END OF MODULE
!%==========================================================================
end module orifice_elements