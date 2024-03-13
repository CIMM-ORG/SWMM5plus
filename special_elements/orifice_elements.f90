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
    use geometry_lowlevel
    use irregular_channel, only: irregular_geometry_from_depth_singular
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
            integer, pointer    :: SpecificOrificeType
        !%------------------------------------------------------------------
        !% Aliases
            SpecificOrificeType   => elemSI(eIdx,esi_Orifice_SpecificType)
        !%------------------------------------------------------------------
        !% --- NOTE the opening of the orifice due to control intervention
        !%     is already set in control_update_setting subroutine

        !% -- get the head and flow direction through orifice
        call common_head_and_flowdirection_singular &
            (eIdx, esr_Orifice_Zcrest, esr_Orifice_NominalDownstreamHead, esi_Orifice_FlowDirection)

        !% --- find effective head difference across orifice element
        call orifice_effective_head_delta (eIdx)

        if ((SpecificOrificeType .eq. (equivalent_orifice_channel)) .or. &
            (SpecificOrificeType .eq. (equivalent_orifice_pipe))            ) then
            !% --- update geometry in elemR for channel/conduit based on
            !%     depth implied by upstream equivalent orifice head
            call orifice_equivalent_geometry_update (eIdx)
            !% --- compute the dynamic discharge coefficient
            call orifice_equivalent_dischargeCoef (eIdx)
            !% --- compute the equivalent orifice flowrate
            call orifice_equivalent_flow (eIdx)
        else
            !% --- find the effective full area of orifice
            call orifice_effective_full_area (eIdx)
            !% --- get the critical head and depth values
            call orifice_critical_head_and_depth (eIdx)
            !% --- find flow on orifice element
            call orifice_flow (eIdx)
        end if

        !% --- apply flap gate adjustment
        call orifice_flapgate_adjustment (eIdx)

        if ((SpecificOrificeType .ne. (equivalent_orifice_channel)) .and. &
            (SpecificOrificeType .ne. (equivalent_orifice_pipe))            ) then

            !% --- apply Villemonte submergence correction
            call orifice_submergence_correction (eIdx)

            !% --- update orifice geometry from head
            call orifice_geometry_update (eIdx)
        else
                !% --- no action for equivalent orifices
        end if

        !% --- update velocity from flowrate and area
        call common_velocity_from_flowrate_singular (eIdx)

        !% --- compute downstream energy head
        call common_outflow_energyhead_singular &
            (eIdx, esr_Orifice_NominalDownstreamHead, esi_Orifice_FlowDirection)

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
        !%-----------------------------------------------------------------------------
        !% Aliases
            !% --- inputs
            SpecificOrificeType   => elemSI(eIdx,esi_Orifice_SpecificType)
            FlowDirection         => elemSI(eIdx,esi_Orifice_FlowDirection)
            EffectiveFullDepth    => elemSR(eIdx,esr_Orifice_EffectiveFullDepth)
            Head                  => elemR(eIdx,er_Head)
            NominalDsHead         => elemSR(eIdx,esr_Orifice_NominalDownstreamHead)
            Zcrest                => elemSR(eIdx,esr_Orifice_Zcrest) 
            Zcrown                => elemSR(eIdx,esr_Orifice_Zcrown)
            hasFlapGate           => elemYN(eIdx,eYN_hasFlapGate)
            !% --- output
            EffectiveHeadDelta    => elemSR(eIdx,esr_Orifice_EffectiveHeadDelta)
        !%------------------------------------------------------------------

        !% --- find the effective head delta
        if (hasFlapGate .and. (FlowDirection < zeroI)) then
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
                    Zcrown = Zcrest + EffectiveFullDepth  !% HACK is effective full depth used before defined?  BRH20240306
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

                case (equivalent_orifice_channel,equivalent_orifice_pipe)
                    if (Head <= elemR(eIdx,er_Zbottom) + setting%ZeroValue%Depth) then
                        !% --- driving head is too small to generate flow
                        EffectiveHeadDelta = zeroR
                    elseif (NominalDsHead <= elemR(eIdx,er_Zbottom) + setting%ZeroValue%Depth) then
                        !% --- driving head is only height above bottom
                        EffectiveHeadDelta = Head - elemR(eIdx,er_Zbottom)
                    else
                        !% --- standard orifice head
                        EffectiveHeadDelta = Head - NominalDsHead
                    end if

                case default
                    print *, 'In ', trim(subroutine_name)
                    print *, 'CODE ERROR unknown orifice type, ', SpecificOrificeType,'  in network'
                    print *, 'which has key ',trim(reverseKey(SpecificOrificeType))
                    call util_crashpoint(7298734)
            end select
        end if

        !% MOVED TO SEPARATE SUBROUTINES BRH20240306
        ! !% --- find full area for flow, and A/L for critical depth calculations
        ! select case (GeometryType)
        !     case (circular)
        !         YoverYfull        = EffectiveFullDepth / FullDepth
        !         if (YoverYfull .le. zeroR) then
        !             EffectiveFullArea =  zeroR
        !             AoverL = zeroR
        !         else
        !             EffectiveFullArea = FullArea * xsect_table_lookup_singular (YoverYfull, ACirc)
        !             AoverL            = onefourthR * EffectiveFullDepth
        !         end if
        !     case (rectangular_closed)
        !         EffectiveFullArea = EffectiveFullDepth * RectangularBreadth
        !         AoverL            = EffectiveFullArea / (twoR * (EffectiveFullDepth + RectangularBreadth))
        !     case default
        !         if (SpecificOrificeType == equivalent_orifice) then 
        !             print *, 'need to fix this '
        !             stop 889874
        !         else
        !             print *, 'In ', trim(subroutine_name)
        !             print *, 'element idx = ',eIdx
        !             print *, 'SpecificOrificeType = ',SpecificOrificeType, ' ',reverseKey(SpecificOrificeType)
        !             print *, 'CODE ERROR geometry type unknown for # ', GeometryType
        !             print *, 'which has key ',trim(reverseKey(GeometryType))
        !             call util_crashpoint(7998734)
        !         end if
        ! end select

        ! !% --- find critical depth to determine weir/orifice flow
        ! select case (SpecificOrificeType)
        !     case (bottom_orifice)
        !         !% find critical height above opening where orifice flow turns into
        !         !% weir flow for Bottom orifice = (C_orifice/C_weir)*(Area/Length)
        !         !% where C_orifice = given orifice coeff, C_weir = weir_coeff/sqrt(2g),
        !         !% Area is the area of the opening, and Length = circumference
        !         !% of the opening. For a basic sharp crested weir, C_weir = 0.414.
        !         CriticalDepth = DischargeCoeff / SharpCrestedWeirCoeff * AoverL
        !         CriticalHead  = CriticalDepth
        !         FractionCritDepth = min(EffectiveHeadDelta / CriticalDepth, oneR)
        !     case (side_orifice)
        !         !% another adjustment to critical depth is needed
        !         !% for weir coeff calculation for side orifice
        !         CriticalDepth = EffectiveFullDepth
        !         CriticalHead  = CriticalDepth / twoR
        !         FractionCritDepth = min(((Head - Zcrest) / EffectiveFullDepth), oneR)
        !     case (equivalent_orifice)
        !         print *, 'need stuff here'
        !         stop 609874
        !     case default
        !         print *, 'In ', trim(subroutine_name)
        !         print *, 'CODE ERROR unknown orifice type, ', SpecificOrificeType,'  in network'
        !         print *, 'which has key ',trim(reverseKey(SpecificOrificeType))
        !         call util_crashpoint(9298734)
        ! end select

    end subroutine orifice_effective_head_delta
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine orifice_effective_full_area (eIdx)
        !%------------------------------------------------------------------
        !% Description: 
        !% calculates the effective full flow area for orifice
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: eIdx
            integer, pointer    :: SpecificOrificeType, GeometryType
            real(8), pointer    :: EffectiveFullArea, EffectiveFullDepth
            real(8), pointer    :: FullArea, FullDepth, RectangularBreadth
            real(8)             :: YoverYfull, AoverL
        !%------------------------------------------------------------------
        !% Aliases
            SpecificOrificeType   => elemSI(eIdx,esi_Orifice_SpecificType)
            GeometryType          => elemSI(eIdx,esi_Orifice_GeometryType)
            EffectiveFullArea     => elemSR(eIdx,esr_Orifice_EffectiveFullArea)
            EffectiveFullDepth    => elemSR(eIdx,esr_Orifice_EffectiveFullDepth)
            FullArea              => elemSR(eIdx,esr_Orifice_FullArea)
            FullDepth             => elemSR(eIdx,esr_Orifice_FullDepth)
            RectangularBreadth    => elemSR(eIdx,esr_Orifice_RectangularBreadth)
        !%------------------------------------------------------------------

        !% --- regular orifice geometry
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
                print *, 'element idx = ',eIdx
                print *, 'SpecificOrificeType = ',SpecificOrificeType, ' ',reverseKey(SpecificOrificeType)
                print *, 'CODE ERROR geometry type unknown for # ', GeometryType
                print *, 'which has key ',trim(reverseKey(GeometryType))
                call util_crashpoint(7998734)
        end select

    end subroutine orifice_effective_full_area
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine orifice_critical_head_and_depth (eIdx)
        !%------------------------------------------------------------------
        !% Description: 
        !% calculates the critical head, depth and fraction of critical depth
        !% for an orifice
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: eIdx
            integer, pointer    :: SpecificOrificeType
            real(8), pointer    :: CriticalDepth, DischargeCoeff, CriticalHead
            real(8), pointer    :: EffectiveFullDepth, FractionCritDepth
            real(8), pointer    :: EffectiveHeadDelta
            real(8), pointer    :: Head, SharpCrestedWeirCoeff, Zcrest
            real(8)             :: AoverL
        !%------------------------------------------------------------------
        !% Aliases
        !% --- inputs
        SpecificOrificeType   => elemSI(eIdx,esi_Orifice_SpecificType)
        CriticalDepth         => elemSR(eIdx,esr_Orifice_CriticalDepth)
        DischargeCoeff        => elemSR(eIdx,esr_Orifice_DischargeCoeff)
        CriticalHead          => elemSR(eIdx,esr_Orifice_CriticalHead)
        EffectiveFullDepth    => elemSR(eIdx,esr_Orifice_EffectiveFullDepth)
        EffectiveHeadDelta    => elemSR(eIdx,esr_Orifice_EffectiveHeadDelta)
        FractionCritDepth     => elemSR(eIdx,esr_Orifice_FractionCriticalDepth)
        Head                  => elemR(eIdx,er_Head)
        SharpCrestedWeirCoeff => Setting%Orifice%SharpCrestedWeirCoefficient
        Zcrest                => elemSR(eIdx,esr_Orifice_Zcrest) 

        !%------------------------------------------------------------------

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
            case (equivalent_orifice_channel, equivalent_orifice_pipe)
                !% --- these should not be used for equivalent orifice, nullvalues if accessed will cause blow up.
                CriticalDepth = nullValueR 
                CriticalHead  = nullValueR
                FractionCritDepth = nullvalueR
            case default
                print *, 'CODE ERROR unknown orifice type, ', SpecificOrificeType,'  in network'
                print *, 'which has key ',trim(reverseKey(SpecificOrificeType))
                call util_crashpoint(9298734)
        end select

    end subroutine orifice_critical_head_and_depth
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
            integer, pointer :: FlowDirection, SpecificOrificeType
            real(8), pointer :: Flowrate, EffectiveHeadDelta,  grav
            real(8), pointer :: dQdH, DischargeCoeff, EffectiveFullArea
            real(8), pointer :: WeirExponent, CriticalHead, FractionCritDepth
            real(8) :: Coef

            character(64) :: subroutine_name = 'orifice_flow'
        !%-----------------------------------------------------------------------------
        !% Aliases
            SpecificOrificeType   => elemSI(eIdx,esi_Orifice_SpecificType)
            FlowDirection         => elemSI(eIdx,esi_Orifice_FlowDirection)
            CriticalHead          => elemSR(eIdx,esr_Orifice_CriticalHead)
            DischargeCoeff        => elemSR(eIdx,esr_Orifice_DischargeCoeff)
            EffectiveHeadDelta    => elemSR(eIdx,esr_Orifice_EffectiveHeadDelta)
            EffectiveFullArea     => elemSR(eIdx,esr_Orifice_EffectiveFullArea)
            FractionCritDepth     => elemSR(eIdx,esr_Orifice_FractionCriticalDepth)
            dQdH                  => elemSR(eIdx,esr_Orifice_dQdHe)
            Flowrate              => elemR(eIdx,er_Flowrate)
            WeirExponent          => Setting%Orifice%TransverseWeirExponent
            grav                  => setting%constant%gravity
        !%------------------------------------------------------------------

        !% --- flow calculation conditions through an orifice
        if ((EffectiveHeadDelta <= zeroR) .or. (FractionCritDepth <= zeroR)) then
            !% --- no flow case
            Flowrate = zeroR
            dQdH     = zeroR
        elseif (FractionCritDepth < oneR) then
            !% --- case where inlet depth is below critical depth thus,
            !%     orifice behaves as a rectangular transverse weir
            Coef     = DischargeCoeff * EffectiveFullArea * sqrt(twoR * grav * CriticalHead)
            Flowrate = real(FlowDirection,real(8)) * Coef * (FractionCritDepth ** WeirExponent)
            dQdH     = WeirExponent * Flowrate / (FractionCritDepth * CriticalHead)
        else
            !% --- standard orifice flow condition
            Coef      = DischargeCoeff * EffectiveFullArea * sqrt(twoR * grav)
            Flowrate  = real(FlowDirection,real(8)) * Coef * sqrt(EffectiveHeadDelta)
            dQdH      = onehalfR * (Flowrate / EffectiveHeadDelta)
        end if

    end subroutine orifice_flow
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
            integer, pointer    :: SpecificOrificeType
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
            SpecificOrificeType=> elemSI(eIdx,esi_Orifice_SpecificType)
            Area               => elemR(eIdx,er_Area) !% --- orifice area, do not use AreaVelocity
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
            !% --- find the actual flow area to calculate the velocity (Not the effective)
            !%     only for true orifices.
            if ((SpecificOrificeType .eq. equivalent_orifice_channel) .or. &
                (SpecificOrificeType .eq. equivalent_orifice_pipe)            ) then
                !% --- no action required because Area is already stored for equivalent orifice
            else
                !% --- true orifice requires computing velocity through area
                call orifice_flow_area (eIdx)
            end if

            if (Area > zeroArea) then
                Velocity = Flowrate / Area
                hLoss    = (fourR / grav) * Velocity * Velocity &
                        * exp(-1.15 * Velocity / sqrt(EffectiveHeadDelta))

            !% --- update the new headDelta to account for flap gate loss
                if ((SpecificOrificeType .eq. equivalent_orifice_channel) .or. &
                    (SpecificOrificeType .eq. equivalent_orifice_pipe)            ) then

                        EffectiveHeadDelta = max(EffectiveHeadDelta - hloss, zeroR)

                else            
                    if (FractionCritDepth < oneR) then
                        !% --- weir flow
                        FractionCritDepth = max(FractionCritDepth - hLoss / CriticalDepth, zeroR)
                    else
                        EffectiveHeadDelta = max(EffectiveHeadDelta - hLoss, zeroR)
                    end if
                end if

                !% --- reset the orifice flow to account for headloss through flap
                if ((SpecificOrificeType .eq. equivalent_orifice_channel) .or. &
                    (SpecificOrificeType .eq. equivalent_orifice_pipe)            ) then
                    call orifice_equivalent_flow (eIdx)
                else
                    call orifice_flow (eIdx)
                end if

            end if
        end if

        !%------------------------------------------------------------------
        !% Closing

    end subroutine  orifice_flapgate_adjustment
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine orifice_flow_area (eIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% find orifice actual flow area (without effective area adjustments)
        !% used to estimate velocity for flapgate
        !% SHOULD NOT BE CALLED FOR AN EQUIVALENT ORIFICE ELEMENT
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
    end subroutine  orifice_submergence_correction
!%==========================================================================
!%
    subroutine orifice_geometry_update (eIdx)
        !%------------------------------------------------------------------
        !% Description:
        !% Computes geometry terms for orifice 
        !% requires either rectangular closed or circular geometry (not for equivalen orifice)
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
    subroutine orifice_equivalent_geometry_update (eIdx)
        !%------------------------------------------------------------------
        !% Description
        !% Updates the geometry for an equivalent orifice representation
        !% of a link
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: eIdx
            integer, pointer :: FlowDirection, thisGeoType, SpecificOrificeType
            integer, pointer :: fUp, fDn
            real(8), pointer :: Atable(:), Ttable(:), Rtable(:), Stable(:)
            real(8), pointer :: Head, ZbtmDn, ZbtmUp, FullArea, FullHydRadius
            real(8), dimension(1) :: Depth, Area, Perimeter, HydRadius
            integer, dimension(1) :: eIdxA
 
        !%------------------------------------------------------------------
        !% Aliases:
            FlowDirection         => elemSI(eIdx,esi_Orifice_FlowDirection)
            thisGeoType           => elemI (eIdx,ei_geometryType)
            SpecificOrificeType   => elemSI(eIdx,esi_Orifice_SpecificType)

            FullArea      => elemR (eIdx,er_FullArea)
            FullHydRadius => elemR (eIdx,er_FullHydRadius)
            Head          => elemR (eIdx,er_Head)  !% -- this will be the upstream head, depends on flow direction

            !% --- face aliases
            fUp           => elemI (eIdx,ei_Mface_uL)
            fDn           => elemI (eIdx,ei_Mface_dL)
            ZbtmDn        => faceR (fDn,fr_Zbottom)
            ZbtmUp        => faceR (fUp,fr_Zbottom)

        !%------------------------------------------------------------------
        !% -- store in size 1 arrays for pure functions
        !%    cannot simply use pointers due to shape matching rules with scalars.
        Depth         = elemR (eIdx,er_Depth)
        Area          = elemR (eIdx,er_Area)
        Perimeter     = elemR (eIdx,er_Perimeter)
        HydRadius     = elemR (eIdx,er_HydRadius)
        eIdxA         = eIdx

        if (FlowDirection < oneI ) then 
            !% --- reversed flow
            Depth = Head - ZbtmDn
        else
            !% --- nominal downstream flow
            Depth = Head - ZbtmUp
        end if    

        if (SpecificOrificeType .eq. equivalent_orifice_channel) then
            !% --- open channel types
            select case (thisGeoType)
                case (parabolic)
                    Area      = llgeo_parabolic_area_from_depth_pure      (eIdxA,Depth)
                    Perimeter = llgeo_parabolic_perimeter_from_depth_pure (eIdxA,Depth)
                case (power_function)
                    Area      = llgeo_powerfunction_area_from_depth_pure      (eIdxA,Depth)
                    Perimeter = llgeo_powerfunction_perimeter_from_depth_pure (eIdxA,Depth)
                case (rectangular)
                    Area      = llgeo_rectangular_area_from_depth_pure      (eIdxA,Depth)
                    Perimeter = llgeo_rectangular_perimeter_from_depth_pure (eIdxA,Depth)
                case (trapezoidal)
                    Area      = llgeo_trapezoidal_area_from_depth_pure      (eIdxA,Depth)
                    Perimeter = llgeo_trapezoidal_perimeter_from_depth_pure (eIdxA,Depth)
                case (triangular)
                    Area      = llgeo_triangular_area_from_depth_pure      (eIdxA,Depth)
                    Perimeter = llgeo_triangular_perimeter_from_depth_pure (eIdxA,Depth)
                case (irregular)
                    Area(1)      = irregular_geometry_from_depth_singular ( &
                        eIdx, tt_area, Depth(1), FullArea, setting%ZeroValue%Area)
                    HydRadius(1) = irregular_geometry_from_depth_singular ( &
                        eIdx, tt_hydradius, Depth(1), FullHydRadius, setting%ZeroValue%Depth)    
                case default
                    print *, 'CODE ERROR Unexpected case default'
                    call util_crashpoint(5709832)
            end select

        elseif (SpecificOrificeType .eq. equivalent_orifice_pipe) then 

            !% --- closed conduit types
            select case (thisGeoType)
                !% --- closed conduits with HydRadius (Rtable) lookup
                case (arch,basket_handle,circular,eggshaped,horiz_ellipse,horseshoe,vert_ellipse)

                    call llgeo_lookup_table_selection (eIdx, Atable, Rtable, Stable, Ttable)
                    Area(1) = llgeo_tabular_from_depth_singular &
                        (eIdx, Depth(1), FullArea, setting%ZeroValue%Depth, setting%ZeroValue%Area, Atable)
                    HydRadius(1) = llgeo_tabular_from_depth_singular &
                        (eIdx, Depth(1), FullHydRadius, setting%ZeroValue%Depth, setting%ZeroValue%Depth, Rtable)

                !% --- closed conduits with section factor (Stable) lookup
                case (catenary,gothic,semi_circular,semi_elliptical)

                    call llgeo_lookup_table_selection (eIdx, Atable, Rtable, Stable, Ttable)
                    Area(1)  = llgeo_tabular_from_depth_singular &
                        (eIdx, Depth(1), FullArea, setting%ZeroValue%Depth, setting%ZeroValue%Area, Atable)
                    HydRadius(1) = llgeo_tabular_hydradius_from_area_and_sectionfactor_singular &
                        (eIdx, Area(1), FullHydradius, setting%ZeroValue%Depth, Stable)

                case (custom)
                    print *, 'Custom conduit geometry not yet supported'
                    call util_crashpoint(70873)

                case (filled_circular)
                    Area(1)      = llgeo_filled_circular_area_from_depth_singular      (eIdx,Depth(1),setting%ZeroValue%Area)
                    Perimeter(1) = llgeo_filled_circular_perimeter_from_depth_singular (eIdx,Depth(1),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
                                    
                case (mod_basket)
                    Area(1)      = llgeo_mod_basket_area_from_depth_singular        (eIdx,Depth(1),setting%ZeroValue%Area)
                    Perimeter(1) = llgeo_mod_basket_perimeter_from_depth_singular   (eIdx,Depth(1),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
                    
                case (rectangular_closed)
                    Area(1)      = llgeo_rectangular_closed_area_from_depth_singular      (eIdx,Depth(1),setting%ZeroValue%Area)
                    Perimeter(1) = llgeo_rectangular_closed_perimeter_from_depth_singular (eIdx,Depth(1),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
                    
                case (rect_round)
                    Area (1)     = llgeo_rect_round_area_from_depth_singular       (eIdx,Depth(1),setting%ZeroValue%Area)
                    Perimeter(1) = llgeo_rect_round_perimeter_from_depth_singular  (eIdx,Depth(1),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
            
                case (rect_triang)
                    Area(1)      = llgeo_rectangular_triangular_area_from_depth_singular      (eIdx,Depth(1),setting%ZeroValue%Area)
                    Perimeter(1) = llgeo_rectangular_triangular_perimeter_from_depth_singular (eIdx,Depth(1),setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
            
                case default 
                    print *, 'CODE ERROR Unexpected case default'
                    call util_crashpoint(570983)
            end select

        else
            print *, 'CODE ERROR: unexpected else'
            call util_crashpoint(71098734)
        end if

        !% --- compute hydraulic radius for cases that only provide perimeter above
        select case (thisGeoType)
            case (parabolic, power_function, rectangular, trapezoidal, triangular, &
                  filled_circular, mod_basket, rectangular_closed, rect_round, rect_triang)
                Perimeter = max(Perimeter,setting%ZeroValue%Topwidth + setting%ZeroValue%Depth)
                HydRadius = llgeo_hydradius_from_area_and_perimeter_pure (eIdxA, Area, Perimeter)
            case default
                !% no action
        end select

        !% --- reset values
        elemR(eIdx,er_Depth)     = Depth(1)
        elemR(eIdx,er_Area)      = Area(1)
        elemR(eIdx,er_Perimeter) = Perimeter(1)
        elemR(eIdx,er_HydRadius) = HydRadius(1)



    end subroutine orifice_equivalent_geometry_update
!%    
!%==========================================================================    
!%========================================================================== 
!%
    subroutine orifice_equivalent_dischargeCoef (eIdx)
        !%------------------------------------------------------------------
        !% Description
        !% computes the discharge coefficient for channel/conduit treated
        !% as an equivalent orifice. Should only be called if eIdx is
        !% an equivalent orifice link
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: eIdx
            real(8), pointer    :: DischargeCoeff, HydRadius, ManningsN
            real(8), pointer    :: Length, grav
        !%------------------------------------------------------------------
        !% Alias 

            DischargeCoeff => elemSR(eIdx,esr_Orifice_DischargeCoeff)
            HydRadius      => elemR (eIdx,er_HydRadius)
            ManningsN      => elemR (eIdx,er_ManningsN)
            Length         => elemR (eIdx,er_Length)
            grav           => setting%Constant%gravity
        !%------------------------------------------------------------------
        
        DischargeCoeff = (HydRadius**twothirdR) &
                          / (ManningsN * sqrt( twoR * grav * Length))   

    end subroutine orifice_equivalent_dischargeCoef
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine orifice_equivalent_flow (eIdx)
        !%------------------------------------------------------------------
        !% Description
        !% computes flowrate for an equivalent orifice used for a conduit/channel
        !% Should only be called when SpecificOrificeType of eIdx is 
        !% equivalent_orifice_channel or ..._pipe
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: eIdx
            integer, pointer    :: FlowDirection
            real(8), pointer    :: DischargeCoeff, EffectiveHeadDelta, dQdH
            real(8), pointer    :: Flowrate, grav
        !%------------------------------------------------------------------
        !% Aliases
            FlowDirection         => elemSI(eIdx,esi_Orifice_FlowDirection)
            DischargeCoeff        => elemSR(eIdx,esr_Orifice_DischargeCoeff)
            EffectiveHeadDelta    => elemSR(eIdx,esr_Orifice_EffectiveHeadDelta)
            dQdH                  => elemSR(eIdx,esr_Orifice_dQdHe)
            Flowrate              => elemR(eIdx,er_Flowrate)
            grav                  => setting%constant%gravity
        !%------------------------------------------------------------------

        !% --- equivalent orifice does not have Critical Depth
        if (EffectiveHeadDelta <= zeroR) then 
            Flowrate = zeroR 
            dQdH     = zeroR
        else
            Flowrate = real(FlowDirection,real(8)) * DischargeCoeff * elemR(eIdx,er_Area) &
                       * sqrt(twoR * grav * EffectiveHeadDelta)
            
            dQdH = onehalfR * Flowrate / EffectiveHeadDelta
        end if

    end subroutine orifice_equivalent_flow
!%
!%==========================================================================
!% END OF MODULE
!%==========================================================================
end module orifice_elements