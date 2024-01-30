module orifice_elements
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Computes diagnostic flow through orifice elements
    !%
    !% Methods:
    !% Follows methods of EPA-SWMM-C
    !%==========================================================================
    use define_globals
    use define_keys
    use define_indexes
    use define_settings, only: setting
    use common_elements
    use adjust
    use geometry_lowlevel, only: llgeo_elldepth_pure
    use define_xsect_tables
    use utility, only: util_sign_with_ones
    use xsect_tables
    use utility_crash, only: util_crashpoint

    implicit none

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
        !%------------------------------------------------------------------
        !% Description:
        !% Calculate flow through an orifice
        !%------------------------------------------------------------------
            integer, intent(in) :: eIdx  !% must be a single element ID
            character(64)       :: subroutine_name = 'orifice_toplevel'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%orifice_elements) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% --- NOTE the opening of the orifice due to control intervention
        !%     is already set in control_update_setting subroutine

        !% -- get the head and flow direction through orifice
        call common_head_and_flowdirection_singular &
            (eIdx, esr_Orifice_Zcrest, esr_Orifice_NominalDownstreamHead, esi_Orifice_FlowDirection)

        !% --- find effective head on orifice element
        call orifice_effective_head_delta (eIdx)

        !% --- find flow on orifice element
        call orifice_flow (eIdx)

        !% --- apply flap gate adjustment
        call orifice_flapgate_adjustment (eIdx)

        !% --- apply Villemonte submergence correction
        call orifice_submergence_correction (eIdx)

        !% --- update orifice geometry from head
        call orifice_geometry_update (eIdx)

        !% --- update velocity from flowrate and area
        call common_velocity_from_flowrate_singular (eIdx)

        !% --- compute downstream energy head
        call common_outflow_energyhead_singular &
         (eIdx, esr_Orifice_NominalDownstreamHead, esi_Orifice_FlowDirection)

        !%------------------------------------------------------------------
        !% Closing
        if (setting%Debug%File%orifice_elements)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine orifice_toplevel
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine orifice_set_setting (eIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% evaluate the orifice setting based on control update
        !% Note that the "Orate" the opening/closing rate is entered in the
        !% EPA-SWMM inp file in hours, but has been converted to seconds.
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: eIdx !% single ID of element
            real(8), pointer    :: FullDepth, EffectiveFullDepth, dt
            real(8), pointer    :: Orate, CurrentSetting, TargetSetting
            real(8) :: deltaRemaining, changeFraction

            character(64) :: subroutine_name = 'orifice_set_settings'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%orifice_elements) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% Aliases
            FullDepth          => elemSR(eIdx,esr_Orifice_FullDepth)
            EffectiveFullDepth => elemSR(eIdx,esr_Orifice_EffectiveFullDepth)
            Orate              => elemSR(eIdx,esr_Orifice_Orate)
            CurrentSetting     => elemR(eIdx,er_Setting)
            TargetSetting      => elemR(eIdx,er_TargetSetting)
            dt                 => setting%Time%Hydraulics%Dt
        !%-----------------------------------------------------------------------------
        !% --- case where adjustment is instantaneous
        if ((Orate == zeroR) .or. (dt == zeroR)) then
            CurrentSetting = TargetSetting
        else
            deltaRemaining = TargetSetting - CurrentSetting  !% fraction of orifice to open/close
            changeFraction = dt / Orate                      !% fraction we can complete in this step
            if (changefraction * (oneR + onefourthR) >= abs(deltaRemaining)) then
                !% --- if this change fraction is within 25% of opening/closing,
                !%     then we complete in the step rather than finish in the next step
                CurrentSetting = TargetSetting
            else
                CurrentSetting = CurrentSetting + sign(changeFraction,deltaRemaining)
            end if
        end if

        !% --- error check
        !%     EPA-SWMM allows the orifice setting to be between 0.0 and 1.0
        if (.not. ((CurrentSetting .ge. zeroR) .and. (CurrentSetting .le. oneR))) then
            print *, 'CODE ERROR orifice element has er_Setting that is not between 0.0 and 1.0'
            call util_crashpoint(623943)
        end if

        !% --- find effective orifice opening
        EffectiveFullDepth = FullDepth * CurrentSetting

        !%------------------------------------------------------------------
        !% Closing
        if (setting%Debug%File%orifice_elements)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine orifice_set_setting
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
     subroutine orifice_effective_head_delta (eIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% Compute the effective difference between head and crest
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: eIdx !% single ID of element
            real(8), pointer    :: CriticalDepth, CriticalHead, FullArea, FractionCritDepth
            real(8), pointer    :: EffectiveFullArea, EffectiveHeadDelta, FullDepth
            real(8), pointer    :: Head, NominalDsHead, RectangularBreadth
            real(8), pointer    :: EffectiveFullDepth, Zcrown, Zcrest
            real(8), pointer    :: DischargeCoeff, SharpCrestedWeirCoeff
            integer, pointer    :: SpecificOrificeType, FlowDirection, GeometryType
            logical, pointer    :: hasFlapGate
            real(8)             :: AoverL, YoverYfull, Zmidpt

            character(64) :: subroutine_name = 'orifice_effective_head_delta'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%orifice_elements) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------
        !% Aliases
            !% --- inputs
            SpecificOrificeType   => elemSI(eIdx,esi_Orifice_SpecificType)
            FlowDirection         => elemSI(eIdx,esi_Orifice_FlowDirection)
            GeometryType          => elemSI(eIdx,esi_Orifice_GeometryType)
            CriticalDepth         => elemSR(eIdx,esr_Orifice_CriticalDepth)
            DischargeCoeff        => elemSR(eIdx,esr_Orifice_DischargeCoeff)
            CriticalHead          => elemSR(eIdx,esr_Orifice_CriticalHead)
            EffectiveFullArea     => elemSR(eIdx,esr_Orifice_EffectiveFullArea)
            EffectiveFullDepth    => elemSR(eIdx,esr_Orifice_EffectiveFullDepth)
            FullArea              => elemSR(eIdx,esr_Orifice_FullArea)
            FullDepth             => elemSR(eIdx,esr_Orifice_FullDepth)
            FractionCritDepth     => elemSR(eIdx,esr_Orifice_FractionCriticalDepth)
            Head                  => elemR(eIdx,er_Head)
            NominalDsHead         => elemSR(eIdx,esr_Orifice_NominalDownstreamHead)
            RectangularBreadth    => elemSR(eIdx,esr_Orifice_RectangularBreadth)
            SharpCrestedWeirCoeff => Setting%Orifice%SharpCrestedWeirCoefficient
            Zcrest                => elemSR(eIdx,esr_Orifice_Zcrest) 
            Zcrown                => elemSR(eIdx,esr_Orifice_Zcrown)
            hasFlapGate           => elemYN(eIdx,eYN_hasFlapGate)
            !% --- output
            EffectiveHeadDelta    => elemSR(eIdx,esr_Orifice_EffectiveHeadDelta)
        !%------------------------------------------------------------------

        !% --- find the effective head delta
        if (hasFlapGate .and. (FlowDirection < zeroR)) then
            EffectiveHeadDelta = zeroR
        else
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
                    !% --- considering the effect of control intervention
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
                    print *, 'In ', trim(subroutine_name)
                    print *, 'CODE ERROR unknown orifice type, ', SpecificOrificeType,'  in network'
                    print *, 'which has key ',trim(reverseKey(SpecificOrificeType))
                    call util_crashpoint(7298734)
            end select
        end if

        !% --- find full area for flow, and A/L for critical depth calculations
        select case (GeometryType)
            case (circular)
                YoverYfull        = EffectiveFullDepth / FullDepth
                if (YoverYfull .le. zeroR) then
                    EffectiveFullArea =  zeroR
                    AoverL = zeroR
                else
                    EffectiveFullArea = FullArea * xsect_table_lookup_singular (YoverYfull, ACirc)
                    AoverL            = onefourthR * EffectiveFullDepth
                end if
            case (rectangular_closed)
                EffectiveFullArea = EffectiveFullDepth * RectangularBreadth
                AoverL            = EffectiveFullArea / (twoR * (EffectiveFullDepth + RectangularBreadth))
            case default
                print *, 'In ', trim(subroutine_name)
                print *, 'element idx = ',eIdx
                print *, 'SpecificOrificeType = ',SpecificOrificeType, ' ',reverseKey(SpecificOrificeType)
                print *, 'CODE ERROR geometry type unknown for # ', GeometryType
                print *, 'which has key ',trim(reverseKey(GeometryType))
                call util_crashpoint(7998734)
        end select

        !% --- find critical depth to determine weir/orifice flow
        select case (SpecificOrificeType)
            case (bottom_orifice)
                !% find critical height above opening where orifice flow turns into
                !% weir flow for Bottom orifice = (C_orifice/C_weir)*(Area/Length)
                !% where C_orifice = given orifice coeff, C_weir = weir_coeff/sqrt(2g),
                !% Area is the area of the opening, and Length = circumference
                !% of the opening. For a basic sharp crested weir, C_weir = 0.414.
                CriticalDepth = DischargeCoeff / SharpCrestedWeirCoeff * AoverL
                CriticalHead  = CriticalDepth
                FractionCritDepth = min(EffectiveHeadDelta / CriticalDepth, oneR)
            case (side_orifice)
                !% another adjustment to critical depth is needed
                !% for weir coeff calculation for side orifice
                CriticalDepth = EffectiveFullDepth
                CriticalHead  = CriticalDepth / twoR
                FractionCritDepth = min(((Head - Zcrest) / EffectiveFullDepth), oneR)
            case default
                print *, 'In ', trim(subroutine_name)
                print *, 'CODE ERROR unknown orifice type, ', SpecificOrificeType,'  in network'
                print *, 'which has key ',trim(reverseKey(SpecificOrificeType))
                call util_crashpoint(9298734)
        end select

        !%------------------------------------------------------------------
        !% Closing
        if (setting%Debug%File%orifice_elements)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine orifice_effective_head_delta
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine orifice_flow (eIdx)
        !%------------------------------------------------------------------
        !% Description: 
        !% calculates the flow in an orifice element
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: eIdx
            integer, pointer :: FlowDirection
            real(8), pointer :: Flowrate, EffectiveHeadDelta, Zcrest, grav
            real(8), pointer :: dQdH, DischargeCoeff, EffectiveFullArea
            real(8), pointer :: WeirExponent, CriticalHead, FractionCritDepth
            real(8) :: Coef

            character(64) :: subroutine_name = 'orifice_flow'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%orifice_elements) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------
        !% Aliases
            FlowDirection         => elemSI(eIdx,esi_Orifice_FlowDirection)
            CriticalHead          => elemSR(eIdx,esr_Orifice_CriticalHead)
            DischargeCoeff        => elemSR(eIdx,esr_Orifice_DischargeCoeff)
            EffectiveHeadDelta    => elemSR(eIdx,esr_Orifice_EffectiveHeadDelta)
            EffectiveFullArea     => elemSR(eIdx,esr_Orifice_EffectiveFullArea)
            FractionCritDepth     => elemSR(eIdx,esr_Orifice_FractionCriticalDepth)
            Zcrest                => elemSR(eIdx,esr_Orifice_Zcrest)
            dQdH                  => elemSR(eIdx,esr_Orifice_dQdHe)
            Flowrate              => elemR(eIdx,er_Flowrate)
            WeirExponent          => Setting%Orifice%TransverseWeirExponent
            grav                  => setting%constant%gravity
        !%------------------------------------------------------------------

        !% --- flow calculation conditions through an orifice
        if ((EffectiveHeadDelta == zeroR) .or. (FractionCritDepth <= zeroR)) then
            !% --- no flow case
            Flowrate = zeroR
            dQdH     = zeroR
        elseif (FractionCritDepth < oneR) then
            !% --- case where inlet depth is below critical depth thus,
            !%     orifice behaves as a rectangular transverse weir
            Coef     = DischargeCoeff * EffectiveFullArea * sqrt(twoR * grav * CriticalHead)
            Flowrate = FlowDirection * Coef * (FractionCritDepth ** WeirExponent)
            dQdH     = WeirExponent * Flowrate / (FractionCritDepth * CriticalHead)
        else
            !% --- standard orifice flow condition
            Coef      = DischargeCoeff * EffectiveFullArea * sqrt(twoR * grav)
            Flowrate  = FlowDirection * Coef * sqrt(EffectiveHeadDelta)
            dQdH      = onehalfR * (Flowrate / EffectiveHeadDelta)
        end if

        !%------------------------------------------------------------------
        !% Closing
        if (setting%Debug%File%orifice_elements) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine orifice_flow
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine orifice_geometry_update (eIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes geometry terms for orifice 
        !% requires either rectangular closed or circular geometry
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: eIdx
            real(8), pointer :: Head, Length, Zbottom,  Zcrown
            real(8), pointer :: Depth, Area, Volume, Topwidth, ellDepth
            real(8), pointer :: Perimeter, HydRadius,  Zcrest, EffectiveFullArea !, HydDepth
            real(8), pointer :: RectangularBreadth, EffectiveFullDepth
            integer, pointer :: GeometryType
            real(8)          :: YoverYfull

            integer, dimension(1) :: iA 
            real(8), dimension(1) :: outA

            character(64) :: subroutine_name = 'orifice_geometry_update'
        !%-----------------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%orifice_elements) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------
        !% Aliases
            GeometryType       => elemSI(eIdx,esi_Orifice_GeometryType)
            Area               => elemR(eIdx,er_Area)
            Depth              => elemR(eIdx,er_Depth)
            Head               => elemR(eIdx,er_Head)
            HydRadius          => elemR(eIdx,er_HydRadius)
            ellDepth           => elemR(eIdx,er_EllDepth)
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
        !%-----------------------------------------------------------------------------
        !% --- find depth over bottom of orifice
        if (Head <= Zcrest) then
            Depth = zeroR
        elseif ((Head > Zcrest) .and. (Head < Zcrown)) then
            Depth =  Head - Zcrest
        else
            Depth = Zcrown - Zcrest
        end if

        !% --- if the orifice is closed or depth is below crest, set all the geometry to zero
        if ((EffectiveFullDepth <= zeroR) .or. (Depth == zeroR)) then
            Area      = zeroR
            Volume    = zeroR
            Topwidth  = zeroR
            Perimeter = zeroR
            HydRadius = zeroR
            ellDepth  = zeroR
        else
            select case (GeometryType)

                case (rectangular_closed)
                    Area      = RectangularBreadth * Depth
                    Volume    = Area * Length !% HACK this is not the correct volume in the element
                    Topwidth  = RectangularBreadth
                    ellDepth  = Depth
                    Perimeter = Topwidth + twoR * Depth
                    HydRadius = Area / Perimeter

                case (circular)
                    YoverYfull  = Depth / EffectiveFullDepth
                    Area        = EffectiveFullArea * &
                            xsect_table_lookup_singular (YoverYfull, ACirc)
                    Volume      = Area * Length
                    Topwidth    = EffectiveFullDepth * &
                            xsect_table_lookup_singular (YoverYfull, TCirc)
                    hydRadius   = onefourthR * EffectiveFullDepth * &
                            xsect_table_lookup_singular (YoverYfull, RCirc)
                    if (hydRadius > zeroR) then
                        Perimeter   = min(Area / hydRadius, &
                            EffectiveFullArea / (onefourthR * EffectiveFullDepth))
                    else
                        Perimeter = setting%ZeroValue%Depth
                    end if
                    iA(1) = eIdx
                    outA = llgeo_elldepth_pure(iA)
                    ellDepth = outA(1)
    
                case default
                    print *, 'CODE ERROR geometry type unknown for # ', GeometryType
                    print *, 'which has key ',trim(reverseKey(GeometryType))
                    call util_crashpoint(2201777)
            end select
        end if

        !% apply geometry limiters
        call adjust_limit_by_zerovalues_singular (eIdx, er_Area,      setting%ZeroValue%Area,     .false.)
        call adjust_limit_by_zerovalues_singular (eIdx, er_Depth,     setting%ZeroValue%Depth,    .false.)
        call adjust_limit_by_zerovalues_singular (eIdx, er_EllDepth,  setting%ZeroValue%Depth,    .false.)
        call adjust_limit_by_zerovalues_singular (eIdx, er_HydRadius, setting%ZeroValue%Depth,    .false.)
        call adjust_limit_by_zerovalues_singular (eIdx, er_Topwidth,  setting%ZeroValue%Topwidth, .false.)
        call adjust_limit_by_zerovalues_singular (eIdx, er_Perimeter, setting%ZeroValue%Topwidth, .false.)
        call adjust_limit_by_zerovalues_singular (eIdx, er_Volume,    setting%ZeroValue%Volume,   .true.)

        !%------------------------------------------------------------------
        !% Closing
        if (setting%Debug%File%orifice_elements) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine  orifice_geometry_update
!%
!%==========================================================================
!%==========================================================================
!%
subroutine orifice_flapgate_adjustment (eIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% headloss adjustment for flap gates
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: eIdx
            real(8), pointer    :: Area, Flowrate, Velocity, CriticalDepth 
            real(8), pointer    :: EffectiveHeadDelta, FractionCritDepth, grav, zeroArea
            logical, pointer    :: hasFlapGate
            real(8)             :: hLoss

            character(64) :: subroutine_name = 'orifice_flapgate_adjustment'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%orifice_elements) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------
        !% Aliases
            Area               => elemR(eIdx,er_Area)
            Flowrate           => elemR(eIdx,er_Flowrate)
            Velocity           => elemR(eIdx,er_Velocity)
            hasFlapGate        => elemYN(eIdx,eYN_hasFlapGate)
            CriticalDepth      => elemSR(eIdx,esr_Orifice_CriticalDepth)
            FractionCritDepth  => elemSR(eIdx,esr_Orifice_FractionCriticalDepth)
            EffectiveHeadDelta => elemSR(eIdx,esr_Orifice_EffectiveHeadDelta)
            grav               => setting%constant%gravity
            zeroArea           => setting%ZeroValue%Area
        !%-----------------------------------------------------------------------------
        if (hasFlapGate) then
            !% --- find the flow area to calculate the velocity
            call orifice_flow_area (eIdx)

            if (Area > zeroArea) then
                Velocity = Flowrate / Area
                hLoss    = (fourR / grav) * Velocity * Velocity &
                        * exp(-1.15 * Velocity / sqrt(EffectiveHeadDelta))

                !% --- update the new headDelta
                if (FractionCritDepth < oneR) then
                    !% --- weir flow
                    FractionCritDepth = max(FractionCritDepth - hLoss / CriticalDepth, zeroR)
                else
                    EffectiveHeadDelta = max(EffectiveHeadDelta - hLoss, zeroR)
                end if

                call orifice_flow (eIdx)
            end if
        end if

        !%------------------------------------------------------------------
        !% Closing
        if (setting%Debug%File%orifice_elements) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine  orifice_flapgate_adjustment
!%
!%==========================================================================
!%==========================================================================
!%
subroutine orifice_flow_area (eIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% find orifice flow area
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: eIdx
            real(8), pointer    :: Area, Depth, Head, Zcrown
            real(8), pointer    :: EffectiveFullArea, Zcrest
            real(8), pointer    :: RectangularBreadth, EffectiveFullDepth
            integer, pointer    :: GeometryType
            real(8)             :: YoverYfull

            character(64) :: subroutine_name = 'orifice_flow_area'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%orifice_elements) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------
        !% Aliases
            GeometryType       => elemSI(eIdx,esi_Orifice_GeometryType)
            Area               => elemR(eIdx,er_Area)
            Depth              => elemR(eIdx,er_Depth)
            Head               => elemR(eIdx,er_Head)
            EffectiveFullArea  => elemSR(eIdx,esr_Orifice_EffectiveFullArea)
            EffectiveFullDepth => elemSR(eIdx,esr_Orifice_EffectiveFullDepth)
            RectangularBreadth => elemSR(eIdx,esr_Orifice_RectangularBreadth)
            Zcrown             => elemSR(eIdx,esr_Orifice_Zcrown)
            Zcrest             => elemSR(eIdx,esr_Orifice_Zcrest)
        !%-----------------------------------------------------------------------------
        !% --- find depth over bottom of orifice
        if (Head <= Zcrest) then
            Depth = zeroR
        elseif ((Head > Zcrest) .and. (Head < Zcrown)) then
            Depth =  Head - Zcrest
        else
            Depth = Zcrown - Zcrest
        end if

        !% --- if the orifice is closed or depth is below crest, set all the geometry to zero
        if ((EffectiveFullDepth <= zeroR) .or. (depth == zeroR)) then
            Area      = zeroR
        else
            select case (GeometryType)

                case (rectangular_closed)
                    Area      = RectangularBreadth * Depth
                case (circular)
                    YoverYfull  = Depth / EffectiveFullDepth
                    Area        = EffectiveFullArea * xsect_table_lookup_singular (YoverYfull, ACirc)
                case default
                    print *, 'CODE ERROR geometry type unknown for # ', GeometryType
                    print *, 'which has key ',trim(reverseKey(GeometryType))
                    call util_crashpoint(6222987)
            end select
        end if

        !% --- apply geometry limiters
        call adjust_limit_by_zerovalues_singular (eIdx, er_Area, setting%ZeroValue%Area, .false.)

        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%orifice_elements) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine  orifice_flow_area
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine orifice_submergence_correction (eIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% compute flow correction for submerged flow 
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: eIdx
            real(8), pointer    :: FractionCritDepth, NominalDsHead, Head, Zcrest
            real(8), pointer    :: Flowrate, WeirExponent, VillemonteExponent
            real(8)             :: ratio
            character(64) :: subroutine_name = 'orifice_submergence_correction'
        !%-----------------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%orifice_elements) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------
        !% Aliases
            Head               => elemR(eIdx,er_Head)
            Flowrate           => elemR(eIdx,er_Flowrate)
            FractionCritDepth  => elemSR(eIdx,esr_Orifice_FractionCriticalDepth)
            NominalDsHead      => elemSR(eIdx,esr_Orifice_NominalDownstreamHead)
            Zcrest             => elemSR(eIdx,esr_Orifice_Zcrest)
            WeirExponent       => Setting%Orifice%TransverseWeirExponent
            VillemonteExponent => Setting%Orifice%VillemonteCorrectionExponent
        !%-----------------------------------------------------------------------------
        !% --- applying Villemonte submergence correction for orifice having submerged flow
        if ((FractionCritDepth < oneR) .and. (NominalDsHead > Zcrest)) then
            ratio    = (NominalDsHead - Zcrest) / (Head - Zcrest)
            Flowrate = Flowrate * ((oneR - (ratio ** WeirExponent)) ** VillemonteExponent)
        else
            !% no correction needed
        end if

        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%orifice_elements) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
                
    end subroutine  orifice_submergence_correction
!%==========================================================================
!% END OF MODULE
!%==========================================================================
end module orifice_elements