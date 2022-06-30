module orifice_elements

    use define_globals
    use define_keys
    use define_indexes
    use define_settings, only: setting
    use common_elements
    use adjust
    use define_xsect_tables
    use utility, only: util_sign_with_ones
    use xsect_tables
    use utility_crash, only: util_crashpoint


    implicit none

    !%-----------------------------------------------------------------------------
    !% Description:
    !% Computes diagnostic flow through orifice elements
    !%-----------------------------------------------------------------------------

    private

    public :: orifice_toplevel
    public :: orifice_set_setting

    contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine orifice_toplevel (eIdx)
        !%-----------------------------------------------------------------------------
        !% Description:
        !%  Calculate flow through an orifice
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: eIdx  !% must be a single element ID

        character(64) :: subroutine_name = 'orifice_toplevel'
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%orifice_elements) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        
        !% find the opening of the orifice due to control intervention
        print *, 'calling orifice_set_setting'
        call orifice_set_setting (eIdx)

        print *, 'calling common_head_and_flowdirection'
        call common_head_and_flowdirection_singular &
            (eIdx, esr_Orifice_Zcrest, esr_Orifice_NominalDownstreamHead, esi_Orifice_FlowDirection)

        print *, 'calling orifice_effective_head_delta'
        !% find effective head on orifice element
         call orifice_effective_head_delta (eIdx)

         print *, 'calling orifice_flow'
        !% find flow on orifice element
        call orifice_flow (eIdx)

        print *, 'calling orifice_geometry_update'
        !% update orifice geometry from head
        call orifice_geometry_update (eIdx)

        print *, 'calling common_velocity_from_flowrate_singular'
         !% update velocity from flowrate and area
        call common_velocity_from_flowrate_singular (eIdx)

        if (setting%Debug%File%orifice_elements)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine orifice_toplevel
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine orifice_set_setting (eIdx)
        !%-----------------------------------------------------------------------------
        !% Description:
        !%   evaluate the orifice setting based on control update
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: eIdx !% single ID of element

        real(8), pointer :: FullDepth, EffectiveFullDepth, dt
        real(8), pointer :: Orate, CurrentSetting, TargetSetting
        real(8) :: delta, step

        character(64) :: subroutine_name = 'orifice_set_settings'
        !%-----------------------------------------------------------------------------
        ! if (crashYN) return
        if (setting%Debug%File%orifice_elements) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        FullDepth          => elemSR(eIdx,esr_Orifice_FullDepth)
        EffectiveFullDepth => elemSR(eIdx,esr_Orifice_EffectiveFullDepth)
        Orate              => elemSR(eIdx,esr_Orifice_Orate)
        CurrentSetting     => elemR(eIdx,er_Setting)
        TargetSetting      => elemR(eIdx,er_TargetSetting)
        dt                 => setting%Time%Hydraulics%Dt
        
        !% case where adjustment is instantaneous
        if ((Orate == zeroR) .or. (dt == zeroR)) then
            CurrentSetting = TargetSetting
        else
            delta = TargetSetting - CurrentSetting
            step = dt / Orate
            if (step + oneOneThounsandthR >= abs(delta)) then
                CurrentSetting = TargetSetting
            else
                !CurrentSetting = CurrentSetting + util_sign_with_ones(delta) * step
                CurrentSetting = CurrentSetting + sign(step,delta)
            end if
        end if

        !% --- error check
        !%     EPA-SWMM allows the orifice setting to be between 0.0 and 1.0
        if (.not. ((CurrentSetting .ge. 0.0) .and. (CurrentSetting .le. 1.0))) then
            print *, 'CODE ERROR: orifice element has er_Setting that is not between 0.0 and 1.0'
            call util_crashpoint(623943)
        end if

        !% find effective orifice opening
        EffectiveFullDepth = FullDepth * CurrentSetting

        if (setting%Debug%File%orifice_elements)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine orifice_set_setting
!%
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
        real(8), pointer    :: EffectiveHeadDelta, NominalDsHead, Head
        real(8), pointer    :: EffectiveFullDepth, Zcrown, Zcrest
        integer, pointer    :: SpecificOrificeType
        real(8)             :: Zmidpt

        character(64) :: subroutine_name = 'orifice_effective_head_delta'
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%orifice_elements) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% inputs
        SpecificOrificeType   => elemSI(eIdx,esi_Orifice_SpecificType)
        EffectiveFullDepth    => elemSR(eIdx,esr_Orifice_EffectiveFullDepth)
        Head                  => elemR(eIdx,er_Head)
        NominalDsHead         => elemSR(eIdx,esr_Orifice_NominalDownstreamHead)
        Zcrest                => elemSR(eIdx,esr_Orifice_Zcrest) 
        Zcrown                => elemSR(eIdx,esr_Orifice_Zcrown)
        !% output
        EffectiveHeadDelta    => elemSR(eIdx,esr_Orifice_EffectiveHeadDelta)
        !%-----------------------------------------------------------------------------
        select case (SpecificOrificeType)
        case (bottom_orifice)
            if (Head <= Zcrest) then
                EffectiveHeadDelta = zeroR
            elseif (NominalDsHead > Zcrest) then
                EffectiveHeadDelta = Head - NominalDsHead
            else
                EffectiveHeadDelta = Head - Zcrest
            end if
        case (side_orifice)
            !% considering the effect of control intervention
            Zcrown = Zcrest + EffectiveFullDepth
            Zmidpt = (Zcrown + Zcrest)/2.0
            if (Head <= Zcrest) then
                EffectiveHeadDelta = zeroR
            elseif (Head < Zcrown) then
                EffectiveHeadDelta = Head - Zcrest
            else
                if (NominalDsHead < Zmidpt) then
                    EffectiveHeadDelta = Head - Zmidpt
                else
                    EffectiveHeadDelta = Head - NominalDsHead
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
        real(8), pointer :: RectangularBreadth, NominalDsHead, FullDepth, FullArea
        real(8), pointer :: DischargeCoeff, EffectiveFullDepth, EffectiveFullArea 
        real(8), pointer :: WeirExponent, VillemonteExponent, SharpCrestedWeirCoeff
        real(8) :: CriticalDepth, AoverL, FractionCritDepth, Coef
        real(8) :: ratio, YoverYfull

        character(64) :: subroutine_name = 'orifice_flow'
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%orifice_elements) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        FlowDirection         => elemSI(eIdx,esi_Orifice_FlowDirection)
        GeometryType          => elemSI(eIdx,esi_Orifice_GeometryType)
        SpecificOrificeType   => elemSI(eIdx,esi_Orifice_SpecificType)
        DischargeCoeff        => elemSR(eIdx,esr_Orifice_DischargeCoeff)
        EffectiveFullDepth    => elemSR(eIdx,esr_Orifice_EffectiveFullDepth)
        EffectiveFullArea     => elemSR(eIdx,esr_Orifice_EffectiveFullArea)
        EffectiveHeadDelta    => elemSR(eIdx,esr_Orifice_EffectiveHeadDelta)
        FullArea              => elemSR(eIdx,esr_Orifice_FullArea)
        FullDepth             => elemSR(eIdx,esr_Orifice_FullDepth)
        NominalDsHead         => elemSR(eIdx,esr_Orifice_NominalDownstreamHead)
        RectangularBreadth    => elemSR(eIdx,esr_Orifice_RectangularBreadth)
        Zcrest                => elemSR(eIdx,esr_Orifice_Zcrest)
        Flowrate              => elemR(eIdx,er_Flowrate)
        Head                  => elemR(eIdx,er_Head)
        SharpCrestedWeirCoeff => Setting%Orifice%SharpCrestedWeirCoefficient
        WeirExponent          => Setting%Orifice%TransverseWeirExponent
        VillemonteExponent    => Setting%Orifice%VillemonteCorrectionExponent
        grav                  => setting%constant%gravity
        !%-----------------------------------------------------------------------------
        !% find full area for flow, and A/L for critical depth calculations
        select case (GeometryType)
            case (circular)
                YoverYfull        = EffectiveFullDepth / FullDepth
                EffectiveFullArea = FullArea * xsect_table_lookup_singular (YoverYfull, ACirc) 
                AoverL            = onefourthR * EffectiveFullDepth
            case (rectangular_closed)
                EffectiveFullArea = EffectiveFullDepth * RectangularBreadth
                AoverL            = EffectiveFullArea / (twoR * (EffectiveFullDepth + RectangularBreadth))
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
            Coef     = DischargeCoeff * EffectiveFullArea * sqrt(twoR * grav * CriticalDepth)
            Flowrate = FlowDirection * Coef * (FractionCritDepth ** WeirExponent)
        else
            !% standard orifice flow condition
            Coef      = DischargeCoeff * EffectiveFullArea * sqrt(twoR * grav)
            Flowrate  = FlowDirection * Coef * sqrt(EffectiveHeadDelta)
        end if

        !% applying Villemonte submergence correction for orifice having submerged weir flow
        if ((FractionCritDepth < oneR) .and. (NominalDsHead > Zcrest)) then
            ratio    = (NominalDsHead - Zcrest) / (Head - Zcrest)
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
        real(8), pointer :: Perimeter, HydDepth, HydRadius,  Zcrest, EffectiveFullArea
        real(8), pointer :: RectangularBreadth, EffectiveFullDepth
        integer, pointer :: GeometryType
        real(8)          :: YoverYfull

        character(64) :: subroutine_name = 'orifice_geometry_update'
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%orifice_elements) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        print *, 'in ',trim(subroutine_name), ' with ',eIdx

        !% pointers
        GeometryType       => elemSI(eIdx,esi_Orifice_GeometryType)
        Area               => elemR(eIdx,er_Area)
        Depth              => elemR(eIdx,er_Depth)
        Head               => elemR(eIdx,er_Head)
        HydDepth           => elemR(eIdx,er_HydDepth)
        HydRadius          => elemR(eIdx,er_HydRadius)
        Length             => elemR(eIdx,er_Length)
        Perimeter          => elemR(eIdx,er_Perimeter)
        Topwidth           => elemR(eIdx,er_Topwidth)
        Volume             => elemR(eIdx,er_Volume)
        Zbottom            => elemR(eIdx,er_Zbottom)
        EffectiveFullArea  => elemSR(eIdx,esr_Orifice_EffectiveFullArea)
        EffectiveFullDepth => elemSR(eIdx,esr_Orifice_EffectiveFullDepth)
        RectangularBreadth => elemSR(eIdx,esr_Orifice_RectangularBreadth)
        Zcrown             => elemSR(eIdx,esr_Orifice_Zcrown)
        Zcrest             => elemSR(eIdx,esr_Orifice_Zcrest)

        !% find depth over bottom of orifice
        if (Head <= Zcrest) then
            Depth = zeroR
        elseif ((Head > Zcrest) .and. (Head < Zcrown)) then
            Depth =  Head - Zcrest
        else
            Depth = Zcrown - Zcrest
        end if

        print *, 'in ',trim(subroutine_name), ' with ',trim(reverseKey(GeometryType))

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
                print *, 'Depth              ',Depth
                print *, 'EffectiveFullDepth ',EffectiveFullDepth
                YoverYfull  = Depth / EffectiveFullDepth
                Area        = EffectiveFullArea * &
                        xsect_table_lookup_singular (YoverYfull, ACirc)  !% 20220506brh removed NACirc
                Volume      = Area * Length
                Topwidth    = EffectiveFullDepth * &
                        xsect_table_lookup_singular (YoverYfull, TCirc) !% 20220506brh removed NTCirc
                HydDepth    = min(Area / Topwidth, EffectiveFullDepth)
                hydRadius   = onefourthR * EffectiveFullDepth * &
                        xsect_table_lookup_singular (YoverYfull, RCirc)  !% 20220506brh removed NRCirc
                Perimeter   = min(Area / hydRadius, &
                        EffectiveFullArea / (onefourthR * EffectiveFullDepth))
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