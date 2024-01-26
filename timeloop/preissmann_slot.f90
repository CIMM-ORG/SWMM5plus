module preissmann_slot
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Procedures for dynamic Preissmann Slot following paper by
    !% Sharior, Hodges, and Vasconcelos (2023)
    !%==========================================================================
    use define_globals
    use define_indexes
    use define_keys
    use define_settings, only: setting
    use utility_crash

    ! use utility_unit_testing, only: util_utest_CLprint

    implicit none

    private
    public :: slot_CC
    public :: slot_JM
    public :: slot_initialize
    !public :: slot_JM_head_PSadd
    !public :: slot_JM_head_PSremove
    public :: slot_CC_adjustments
    public :: slot_JM_adjustments
    public :: slot_JB_computation
    public :: slot_Vshaped_adjust

    integer :: printJM = 135
    integer :: printJB1 = 136
    integer :: printJB2 = 137
    integer :: printUp = 134
    integer :: printDn = 146
    integer :: stepCut = 57000


    contains
!%    
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine slot_initialize (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% initializes slot depth, depth, and PS logicals for surcharge
        !% or non-surcharged conditions
        !% 20230116 -- uses initial head for controlling value
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisP(:)
        !%------------------------------------------------------------------

        where (elemR(thisP,er_Head) < elemR(thisP,er_Zcrown))
            elemR(thisP,er_SlotDepth)         = zeroR   
            !% --- depth remains unchanged
            elemYN(thisP,eYN_isPSsurcharged)  = .false. 
            elemYN(thisP,eYN_isSurcharged)    = .false.    
        elsewhere
            elemR(thisP,er_SlotDepth)         = elemR(thisP,er_Head) - elemR(thisP,er_Zcrown)
            elemR(thisP,er_Depth)             = elemR(thisP,er_FullDepth)
            elemYN(thisP,eYN_isPSsurcharged)  = .true.
            elemYN(thisP,eYN_isSurcharged)    = .true.
        endwhere

    end subroutine slot_initialize
!%
!%==========================================================================


!%==========================================================================
!%
    subroutine slot_CC_adjustments (thisP)
        !%------------------------------------------------------------------
        !% Description:
        !% This subroutine adds back the slot geometry in all the closed elements
        !%------------------------------------------------------------------
        integer, intent(in) :: thisP(:)
        integer, pointer    :: SlotMethod
        real(8), pointer    :: SlotVolume(:), SlotDepth(:), dSlotDepth(:)
        real(8), pointer    :: volume(:), SlotArea(:) !, ell(:)
        real(8), pointer    :: SlotDepth_N0(:), SlotWidth(:)
        real(8), pointer    :: head(:),  fullDepth(:)
        real(8), pointer    :: Overflow(:), zbottom(:)
        logical, pointer    :: isSlot(:)

        character(64) :: subroutine_name = 'geo_CC_slot_adjustments'
        !%------------------------------------------------------------------
        !% Preliminaries
            !% --- exit if not using PS
            if (.not. setting%Solver%PreissmannSlot%useSlotTF) return

            if (setting%Debug%File%geometry) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% Aliases
            !area       => elemR(:,er_Area)
            !depth      => elemR(:,er_Depth)
            dSlotDepth => elemR(:,er_dSlotDepth)
            fullDepth  => elemR(:,er_FullDepth)
            !fullvolume => elemR(:,er_FullVolume)
            !fullArea   => elemR(:,er_FullArea)
            head       => elemR(:,er_Head)
            SlotWidth  => elemR(:,er_SlotWidth)
            SlotVolume => elemR(:,er_SlotVolume)
            SlotDepth  => elemR(:,er_SlotDepth)
            SlotDepth_N0 => elemR(:,er_SlotDepth_N0)
            SlotArea   => elemR(:,er_SlotArea)
            !SlotHydRad => elemR(:,er_SlotHydRadius)
            volume     => elemR(:,er_Volume)
            zbottom    => elemR(:,er_Zbottom)
            isSlot     => elemYN(:,eYN_isPSsurcharged)

            !% pointer to necessary settings struct
            SlotMethod => setting%Solver%PreissmannSlot%Method
        !%-----------------------------------------------------------------------------

            ! call util_utest_CLprint('       AAA slot  - - - - - - - - - - ')

        !% CC slot adjustment
        select case (SlotMethod)
            case (StaticSlot)
                where (isSlot(thisP)) 
                    volume(thisP)    = volume(thisP)  + SlotVolume(thisP)
                    head(thisP)      = zbottom(thisP) + fullDepth(thisP) + SlotDepth(thisP)
                end where 

                ! call util_utest_CLprint('       BBB slot  - - - - - - - - - - ')

            case (DynamicSlot,SplitDynamicSlot)
                where (isSlot(thisP)) 
                    volume(thisP)    = volume(thisP)  + SlotVolume(thisP)
                    SlotDepth(thisP) = max(SlotDepth_N0(thisP) + dSlotDepth(thisP), zeroR) 
                    head(thisP)      = max(zbottom(thisP) + fullDepth(thisP) + SlotDepth(thisP), zbottom(thisP))
                elsewhere
                    SlotDepth(thisP)  = zeroR
                end where 

                !% --- NOTE: for the dynamic slot the Slot width is updated in slot_CC,
                !%     but for the split slot it must be delayed until the slot depth
                !%     is set
                if (SlotMethod == SplitDynamicSlot) then 
                    where (isSlot(thisP)) 
                        SlotWidth(thisP) = SlotArea(thisP) / SlotDepth(thisP)
                    end where
                end if

                ! call util_utest_CLprint('       CCC slot  - - - - - - - - - - ')
            case default 
                print *, 'CODE ERROR unexpected case default'
                call util_crashpoint(6552981)
        end select

        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%geometry) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
                
    end subroutine slot_CC_adjustments
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine slot_JM_adjustments (thisCol, Npack)   
        !%------------------------------------------------------------------
        !% Description:
        !% Adjusts geometry in JM for Preissmann Slot
        !% Uses the Static Slot approach (only)
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisCol, Npack
            integer, pointer    :: thisP(:), SlotMethod
            real(8), pointer    :: volume(:), depth(:),  head(:), ellDepth(:)
            real(8), pointer    :: fullDepth(:), zbottom(:)
            real(8), pointer    :: SlotVolume(:), SlotDepth(:), SlotDepth_N0(:)
            real(8), pointer    :: dSlotDepth(:), SlotWidth(:), SlotArea(:)
            logical, pointer    :: isSlot(:)

            character(64) :: subroutine_name = 'slot_JM_adjustments'
        !%------------------------------------------------------------------
        !% Preliminaries:
            !% --- exit if not using PS
            if (.not. setting%Solver%PreissmannSlot%useSlotTF) return
            if (Npack < 1) return
            if (setting%Debug%File%geometry) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% Aliases
            thisP        => elemP(1:Npack,thisCol)
            volume       => elemR(:,er_Volume)
            depth        => elemR(:,er_Depth)
            dSlotDepth   => elemR(:,er_dSlotDepth)
            fullDepth    => elemR(:,er_FullDepth)
            head         => elemR(:,er_Head)
            ellDepth     => elemR(:,er_EllDepth)
            SlotArea     => elemR(:,er_SlotArea)
            SlotWidth    => elemR(:,er_SlotWidth)
            SlotVolume   => elemR(:,er_SlotVolume)
            SlotDepth    => elemR(:,er_SlotDepth)
            SlotDepth_N0 => elemR(:,er_SlotDepth_N0)
            isSlot       => elemYN(:,eYN_isPSsurcharged)
            zbottom      => elemR(:,er_Zbottom)
            !% pointer to necessary settings struct
            SlotMethod => setting%Solver%PreissmannSlot%Method
        !%------------------------------------------------------------------

        select case (SlotMethod)
            case (StaticSlot)
                where (isSlot(thisP)) 
                    !% brh 20230927 -- volume is not adjusted in Slot_JM
                    !volume(thisP)    = volume(thisP)  + SlotVolume(thisP)
                    head(thisP)      = zbottom(thisP) + fullDepth(thisP) + SlotDepth(thisP)
                end where 

            case (DynamicSlot,SplitDynamicSlot)
                where (isSlot(thisP)) 
                    SlotDepth(thisP) = max(SlotDepth_N0(thisP) + dSlotDepth(thisP), zeroR) 
                    head(thisP)      = max(zbottom(thisP) + fullDepth(thisP) + SlotDepth(thisP), zbottom(thisP))
                elsewhere
                    SlotDepth(thisP)  = zeroR
                end where 

                !% --- NOTE: for the dynamic slot the Slot width is updated in slot_CC,
                !%     but for the split slot it must be delayed until the slot depth
                !%     is set
                if (SlotMethod == SplitDynamicSlot) then 
                    where (isSlot(thisP)) 
                        SlotWidth(thisP) = SlotArea(thisP) / SlotDepth(thisP)
                    end where
                end if
            case default 
                print *, 'CODE ERROR unexpected case default'
                call util_crashpoint(655298)
        end select

        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%geometry) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine slot_JM_adjustments
!%
!%==========================================================================    
!%==========================================================================
!%
    subroutine slot_JB_computation (thisColP_JM)
        !%------------------------------------------------------------------
        !% Description:
        !% Slot computation for Junction Branches
        !% This both computes the slot values and also adjusts the volume
        !% and depth. The head for JB was already set to
        !% Note that overflow/ponding from JB is not allowed, so this does
        !% not consider the number of barrels.
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisColP_JM
            integer, pointer :: Npack, thisP(:), tM, BranchExists(:)
            real(8), pointer :: head(:), length(:), volume(:), zcrown(:)
            real(8), pointer :: fullArea(:), fullVolume(:), fPNumber(:), PNumber(:), PCelerity(:)
            real(8), pointer :: dSlotArea(:), dSlotVol(:), dSlotDepth(:)
            real(8), pointer :: SlotVolume(:), SlotDepth(:), SlotArea(:), SlotVolN0(:)
            real(8), pointer :: SlotWidth(:), SurchargeTime(:), PnumberInitial(:)
            real(8), pointer :: grav, TargetPCelerity, Dt, DecayRate
            logical, pointer :: isSlot(:) , isfSlot(:), isSurcharge(:), canSurcharge(:), isfBlocked(:)
            integer, pointer :: SlotMethod, fUp(:), fDn(:)
            integer :: tB, ii, kk

            character(64) :: subroutine_name = 'slot_JB_computation'
        !%------------------------------------------------------------------
        !% Preliminaries:
            !% --- exit if not using PS
            if (.not. setting%Solver%PreissmannSlot%useSlotTF) return
            Npack => npack_elemP(thisColP_JM)
            if (Npack < 1) return
            if (setting%Debug%File%geometry) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% Aliases
            thisP         => elemP(1:Npack,thisColP_JM)
            head          => elemR(:,er_Head)
            length        => elemR(:,er_Length)
            fullArea      => elemR(:,er_FullArea)
            fullVolume    => elemR(:,er_FullVolume)
            volume        => elemR(:,er_Volume)
            zcrown        => elemR(:,er_Zcrown)
            fUp           => elemI(:,ei_Mface_uL)
            fDn           => elemI(:,ei_Mface_dL)
            BranchExists  => elemSI(:,esi_JB_Exists)
            canSurcharge  => elemYN(:,eYN_canSurcharge)
            grav          => setting%Constant%gravity
        !% Slot Aliases
            dSlotArea      => elemR(:,er_dSlotArea)
            dSlotDepth     => elemR(:,er_dSlotDepth)
            dSlotVol       => elemR(:,er_dSlotVolume)
            PNumber        => elemR(:,er_Preissmann_Number)
            PCelerity      => elemR(:,er_Preissmann_Celerity)
            SlotVolN0      => elemR(:,er_SlotVolume_N0)
            SlotVolume     => elemR(:,er_SlotVolume)
            SlotDepth      => elemR(:,er_SlotDepth)
            SlotArea       => elemR(:,er_SlotArea)
            SlotWidth      => elemR(:,er_SlotWidth)
            SurchargeTime  => elemR(:,er_Surcharge_Time)
            PnumberInitial => elemR(:,er_Preissmann_Number_initial)
            fPNumber       => faceR(:,fr_Preissmann_Number)
            isSurcharge    => elemYN(:,eYN_isSurcharged)
            isSlot         => elemYN(:,eYN_isPSsurcharged)
            isfSlot        => faceYN(:,fYN_isPSsurcharged)
            isfBlocked     => faceYN(:,fYN_isAirflowBlocked)
            SlotMethod     => setting%Solver%PreissmannSlot%Method
            TargetPCelerity=> setting%Solver%PreissmannSlot%TargetCelerity
            DecayRate      => setting%Solver%PreissmannSlot%DecayRate
            Dt             => setting%Time%Hydraulics%Dt
            !Alpha           => setting%Solver%PreissmannSlot%Alpha
        !%------------------------------------------------------------------
        !% --- JB slot adjustment
        !%     cycle through the all the main junctions and each of its branches
        do ii=1,Npack
            tM => thisP(ii) !% junction main ID
            ! --- handle the upstream branches
            do kk=1,max_branch_per_node,2
                tB = tM + kk  !% JB branch ID
                
                if (BranchExists(tB)==1) then
                    !% --- initialize slot
                    isSlot(tB)     = .false.
                    isSurcharge(tB)= .false.
                    SlotDepth(tB)  = zeroR
                    SlotArea(tB)   = zeroR
                    
                    SlotVolume(tB) = zeroR
                    PCelerity(tB)  = zeroR

                    !% --- assuming a slot if the head is above the crown
                    !%     or the upstream CC is in a slot
                    if ((head(tB) .gt. zcrown(tB)) .and. (canSurcharge(tB))) then
                        !% --- classify face as slot ONLY if JM is slot
                        if (isSlot(tM)) then 
                            isfSlot(fUp(tB))    = .true.
                            isfBlocked(fUp(tB)) = .true.
                        else 
                            isfSlot(fUp(tB))    = .false.
                            isfBlocked(fUp(tB)) = .false.
                        end if
                        isSlot(tB)      = .true.
                        !isfSlot(fUp(tB)) = .true.
                        isSurcharge(tB) = .true.
                        !PNumber(tB)    = fPNumber(fUp(tB))
                        PNumber(tB)     = max(onehalfR * (fPNumber(fUp(tB)) + PNumber(tM)), oneR)
                        PCelerity(tB)   = min(TargetPCelerity / PNumber(tB), TargetPCelerity)
                        SlotDepth(tB)   = max(head(tB) - zcrown(tB), zeroR)   
                        SlotArea(tB)    = (SlotDepth(tB) * (PNumber(tB)**twoR) * grav * &
                                            fullArea(tB)) / (TargetPCelerity ** twoR)
                        SlotVolume(tB)  = SlotArea(tB) * length(tB)
                        
                        !% --- add the slot geometry back to previously solved geometry
                        volume(tB) = volume(tB)  + SlotVolume(tB)
                        !Overflow(tB) = zeroR  !% defined as zero (no oveflow from JB)
                    else 
                        !% --- excess head at JB in an open channel represents
                        !%     head from ponding of JM, so do not adjust.
                        !%     Volume and depth are retained at full levels
                        isfSlot(fUp(tB))    = .false.
                        isfBlocked(fUp(tB)) = .false.
                    end if  

                    !% --- increase surcharge time if edge and junction are both surcharge
                    if (isfSlot(fUp(tB)) .and. isSlot(tM)) then 
                        SurchargeTime(tB) = SurchargeTime(tB) + Dt / twoR
                    else
                        SurchargeTime(tB) = zeroR
                    end if

                    !% --- find the new preissmann number for all the closed elements
                    PNumber(tB) = (PnumberInitial(tB) - oneR) * exp((- SurchargeTime(tB) * tenR)/ DecayRate) + oneR

                end if
            end do
            !% --- handle the downstream branches
            do kk=2,max_branch_per_node,2
                tB = tM + kk
                if (BranchExists(tB)==1) then
                    !% --- initialize slot
                    isSlot(tB)     = .false.
                    isSurcharge(tB)= .false.
                    SlotDepth(tB)  = zeroR
                    SlotArea(tB)   = zeroR
                    SlotWidth(tB)  = zeroR
                    SlotVolume(tB) = zeroR
                    PCelerity(tB)  = zeroR

                    !%     or the downstream CC is in a slot
                    if ((head(tB) .gt. zcrown(tB)) .and. (canSurcharge(tB))) then
                        !% --- classify face as slot ONLY if JM is slot
                        if (isSlot(tM)) then 
                            isfSlot(fDn(tB))    = .true.
                            isfBlocked(fDn(tB)) = .true.
                        else 
                            isfSlot(fDn(tB))    = .false.
                            isfBlocked(fDn(tB)) = .false.
                        end if
                        isSlot(tB)     = .true.
                        isSurcharge(tB)= .true.
                        !ifSlot(fDn(tB)) = .true.
                        !PNumber(tB)    = fPNumber(fDn(tB))
                        PNumber(tB) = max(onehalfR * (fPNumber(fDn(tB)) + PNumber(tM)), oneR)
                        PCelerity(tB)  = min(TargetPCelerity / PNumber(tB), TargetPCelerity)

                        SlotDepth(tB)  = max(head(tB) - zcrown(tB), zeroR)    

                        SlotArea(tB)   = (SlotDepth(tB) * (PNumber(tB)**twoR) * grav * &
                                            fullArea(tB)) / (TargetPCelerity ** twoR)
                        SlotVolume(tB) = SlotArea(tB) * length(tB)

                        !% --- add the slot geometry back to previously solved geometry
                        volume(tB) = volume(tB)  + SlotVolume(tB)
                       ! Overflow(tB) = zeroR !% defined as zero (no overflow from JB)
                    else 
                        !% --- excess head at JB in an open channel represents
                        !%     head from ponding of JM, so do not adjust level.
                        !%     Volume and depth are retained at full levels
                        isfSlot(fDn(tB))    = .false.
                        isfBlocked(fDn(tB)) = .false.
                    end if

                    !% --- increase surcharge time if edge and junction are
                    if (isfSlot(fDn(tB)) .and. isSlot(tM)) then 
                        SurchargeTime(tB) = SurchargeTime(tB) + Dt / twoR
                    else
                        SurchargeTime(tB) = zeroR
                    end if

                    !% --- find the new preissmann number for all the closed elements
                    PNumber(tB) = (PnumberInitial(tB) - oneR) * exp((- SurchargeTime(tB) * tenR)/ DecayRate) + oneR

                end if
            end do
        end do
                  
        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%geometry) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine slot_JB_computation
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%    
    subroutine slot_CC (thisP, isSingularYN)
        !%------------------------------------------------------------------
        !% Description:
        !% Compute Preissmann slot for conduits in ETM methods
        !% 
        !% Note that conduits are not allowed to overflow/pond, so there
        !% is no need to consider the number of barrels.
        !%------------------------------------------------------------------
        !% Declarations
            !integer, intent(in) :: thisCol, Npack
            integer, intent(in) :: thisP(:)
            logical, intent(in) :: isSingularYN
            integer, pointer    :: SlotMethod, fUp(:), fDn(:)
            real(8), pointer    :: fullarea(:), PNumberOld(:) !, ellMax(:)
            real(8), pointer    :: fullVolume(:), length(:), PNumber(:), PCelerity(:) 
            real(8), pointer    :: SlotWidth(:), SlotVolume(:), SlotDepth(:), SlotArea(:)
            real(8), pointer    :: SlotDepth_N0(:)
            real(8), pointer    :: dSlotVol(:), dSlotArea(:), dSlotDepth(:), SlotVolN0(:), volume(:) 
            real(8), pointer    :: velocity(:), fPNumber(:), SurchargeTime(:), PnumberInitial(:)
            real(8), pointer    :: TargetPCelerity, grav, Dt, cfl, DecayRate !, Alpha
            logical, pointer    :: isSlot(:), isfSlot(:), isSurcharge(:), isfBlocked(:)
            character(64) :: subroutine_name = "slot_CC"
        !%------------------------------------------------------------------
        !% Preliminaries
            !% --- exit if not using PS
            if (.not. setting%Solver%PreissmannSlot%useSlotTF) return
            !% --- exit if no closed CC elements
            ! if (Npack < 1) return
        !%------------------------------------------------------------------
        !% Aliases
        !% HACK -- many of these aliases are unused. Need to clean up 20220913 brh
            !% --- elemR data
            dSlotArea  => elemR(:,er_dSlotArea)
            dSlotDepth => elemR(:,er_dSlotDepth)
            dSlotVol   => elemR(:,er_dSlotVolume)
            fullArea   => elemR(:,er_FullArea)
            fullVolume => elemR(:,er_FullVolume)
            length     => elemR(:,er_Length)
            SlotVolN0  => elemR(:,er_SlotVolume_N0)
            PNumber    => elemR(:,er_Preissmann_Number)
            PNumberOld  => elemR(:,er_Preissmann_Number_N0)
            PCelerity  => elemR(:,er_Preissmann_Celerity)
            SlotWidth  => elemR(:,er_SlotWidth)
            SlotVolume => elemR(:,er_SlotVolume)
            SlotDepth  => elemR(:,er_SlotDepth)
            SlotDepth_N0 => elemR(:,er_SlotDepth_N0)
            SlotArea   => elemR(:,er_SlotArea)
            volume     => elemR(:,er_Volume)
            velocity   => elemR(:,er_velocity)
            SurchargeTime  => elemR(:,er_Surcharge_Time)
            PnumberInitial => elemR(:,er_Preissmann_Number_initial)
            !% --- pointer to elemYN column
            isSurcharge=> elemYN(:,eYN_isSurcharged)
            isSlot     => elemYN(:,eYN_isPSsurcharged)
            isfSlot    => faceYN(:,fYN_isPSsurcharged)
            isfBlocked => faceYN(:,fYN_isAirflowBlocked)
            !% --- pointers to elemI columns
            fUp        => elemI(:,ei_Mface_uL)
            fDn        => elemI(:,ei_Mface_dL)
            !% --- pointer to faceR column
            fPNumber   => faceR(:,fr_Preissmann_Number)
            !% --- pointer to necessary settings struct
            SlotMethod          => setting%Solver%PreissmannSlot%Method
            TargetPCelerity     => setting%Solver%PreissmannSlot%TargetCelerity
            DecayRate           => setting%Solver%PreissmannSlot%DecayRate
            !Alpha               => setting%Solver%PreissmannSlot%Alpha
            cfl                 => setting%VariableDT%CFL_target
            grav                => setting%Constant%gravity
            Dt                  => setting%Time%Hydraulics%Dt
        !%------------------------------------------------------------------
        !% --- common initialization
        SlotVolume(thisP)   = zeroR
        !SlotArea(thisP)     = zeroR
        PCelerity(thisP)    = zeroR
        isSlot(thisP)       = .false.
        isSurcharge(thisP)  = .false.
        isfSlot(fUp(thisP)) = .false.
        isfSlot(fDn(thisP)) = .false.
        isfBlocked(fUp(thisP)) = .false.
        isfBlocked(fDn(thisP)) = .false.

        !% --- Select the type of slot method
        select case (SlotMethod)
            !% --- for a static slot, the preissmann number will always be one.
            case (StaticSlot)
                !% --- initialize static slot
                PNumber(thisP)    = oneR  !% static slot requirement
                SlotWidth(thisP)  = zeroR
                SlotDepth(thisP)  = zeroR
                !% --- find out the slot volume/ area/ and the faces that are surcharged
                where (volume(thisP) >= fullVolume(thisP))
                    !% --- find slot properties
                    SlotVolume(thisP) = max(volume(thisP) - fullVolume(thisP), zeroR)
                    SlotArea(thisP)   = SlotVolume(thisP) / length(thisP)
                    isfSlot(fUp(thisP)) = .true.
                    isfSlot(fDn(thisP)) = .true.
                    isSlot(thisP)       = .true.
                    isSurcharge(thisP)  = .true.
                    PCelerity(thisP)  = min(TargetPCelerity / PNumber(thisP), TargetPCelerity)
                    SlotWidth(thisP)  = (grav * fullarea(thisP)) / (PCelerity(thisP) ** twoR)
                    !%TEST VALUE: SlotWidth(thisP) = elemR(thisP,er_BreadthMax) * 0.1d0
                    SlotDepth(thisP)  = SlotArea(thisP) / SlotWidth(thisP)
                    !% Air passageway is blocked
                    isfBlocked(fUp(thisP)) = .true.
                    isfBlocked(fDn(thisP)) = .true.
                end where
            
            !% --- for dynamic slot, preissmann number is adjusted
            case (DynamicSlot,SplitDynamicSlot)

                !% --- initialize dynamic slot
                dSlotVol(thisP)   = zeroR
                dSlotArea(thisP)  = zeroR
                dSlotDepth(thisP) = zeroR

                ! if ((setting%Time%Step > 60486) .and. (.not. isSingularYN)) then 
                !     print *, 'VOLUME HERE ',volume(120)- fullVolume(120)
                ! end if

                !% ---find out the slot volume/ area/ and the faces that are surcharged
                where (volume(thisP) > fullVolume(thisP))
                    !% --- find slot properties
                    SlotVolume(thisP) = max(volume(thisP) - fullVolume(thisP), zeroR)
                    SlotArea(thisP)   = SlotVolume(thisP) / length(thisP)  
                    !% --- logicals
                    isfSlot(fUp(thisP)) = .true.
                    isfSlot(fDn(thisP)) = .true.
                    isSlot(thisP)       = .true.
                    isSurcharge(thisP)  = .true.
                    !% --- smooth out the preissmann number before celerity calculation
                    PNumber(thisP) = max(onehalfR * (fPNumber(fUP(thisP)) + fPNumber(fDn(thisP))), oneR)
                    !% --- find the preissmann celerity from the preissmann number
                    PCelerity(thisP) = min(TargetPCelerity / PNumber(thisP), TargetPCelerity)
                    !% --- find the change in slot volume
                    dSlotVol(thisP)   = SlotVolume(thisP) - SlotVolN0(thisP)
                    !% --- find the change in slot area
                    dSlotArea(thisP)  = dSlotVol(thisP) / length(thisP)
                    !% Air passageway is blocked
                    isfBlocked(fUp(thisP)) = .true.
                    isfBlocked(fDn(thisP)) = .true.
                endwhere

                if (SlotMethod == DynamicSlot) then
                    where (volume(thisP) > fullVolume(thisP))
                        !% --- find the change in slot depth
                        dSlotDepth(thisP) = (dSlotArea(thisP)  * (PCelerity(thisP) ** twoI)) / (grav * (fullArea(thisP)))
                        !% --- find the slot width
                        SlotWidth(thisP) =  dSlotArea(thisP) / dSlotDepth(thisP) 
                    endwhere 
                else !% SplitDynamicSlot
                    !% --- apply a different for dropping head
                    where (volume(thisP) > fullVolume(thisP))
                        where (dSlotVol > zeroR)
                            !% --- positive increase is the same as for standard dynamic slot
                            dSlotDepth(thisP) = (dSlotArea(thisP)  * (PCelerity(thisP) ** twoI)) / (grav * (fullArea(thisP)))   
                            !% --- find the slot width
                            SlotWidth(thisP) =  dSlotArea(thisP) / dSlotDepth(thisP)     
                        elsewhere
                            !% --- decrease is proportional to available depth (i.e., an average slot width)
                            dSlotDepth(thisP) = dSlotArea(thisP) * SlotDepth_N0(thisP) / SlotArea(thisP)
                            !% --- note width is not computed here!
                            !% --- positive increase is the same as for standard dynamic slot
                            !dSlotDepth(thisP) = (dSlotArea(thisP)  * (PCelerity(thisP) ** twoI)) / (grav * (fullArea(thisP)))   
                            !% --- find the slot width
                            SlotWidth(thisP) =  dSlotArea(thisP) / dSlotDepth(thisP)   
                        endwhere
                    endwhere
                endif

                !% --- reset isfSlot to find ventilated positions
                where (.not. isSlot(thisP))
                    isfSlot(fUp(thisP)) = .false.
                    isfSlot(fDn(thisP)) = .false.
                    !% Air passageway is not blocked
                    isfBlocked(fUp(thisP)) = .false.
                    isfBlocked(fDn(thisP)) = .false.
                end where

                !% --- not changing preissmann number at diagnostic adjacent elements
                where (faceYN(fUP(thisP),fYN_isDiag_adjacent_all))
                    isfSlot(fUp(thisP)) = .false.
                end where

                where (faceYN(fDn(thisP),fYN_isDiag_adjacent_all))
                    isfSlot(fDn(thisP)) = .false.
                end where

                where (isfSlot(fUp(thisP)) .and. isfSlot(fDn(thisP)))
                    !% --- calculate surcharge time
                    SurchargeTime(thisP) = SurchargeTime(thisP) + Dt / twoR 
                elsewhere
                    !% --- reset the surcharge timer
                    SurchargeTime(thisP) = zeroR
                end where

                !% --- find the new preissmann number for all the closed elements
                PNumber(thisP) = (PnumberInitial(thisP) - oneR) * exp((- SurchargeTime(thisP) * tenR)/ DecayRate) + oneR

            case default
                !% --- should not reach this stage
                print*, 'In ', subroutine_name
                print *, 'CODE ERROR Slot Method type unknown for # ', SlotMethod
                print *, 'which has key ',trim(reverseKey(SlotMethod))
                call util_crashpoint(668732)

        end select

    end subroutine slot_CC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine slot_Vshaped_adjust ()
        !%------------------------------------------------------------------
        !% Description:
        !% Adjust the slots where v-shaped head occurs to ensure the volume
        !% is preserved for the adjusted head.
        !% Is implemented over both CC and JM simultaneously
        !%------------------------------------------------------------------
        !% Declarations
            integer, pointer    :: SlotMethod, fUp(:), fDn(:)
            real(8), pointer    :: fullArea(:), fullDepth(:), PNumber(:), PCelerity(:)  
            real(8), pointer    :: SlotWidth(:), SlotDepth(:), SlotDepth_N0(:), SlotArea(:)
            real(8), pointer    :: dSlotArea(:), dSlotDepth(:), eHead(:) 
            real(8), pointer    :: PnumberInitial(:),  zCrown(:)
            real(8), pointer    :: TargetPCelerity, grav, Dt
            logical, pointer    :: isSurcharge(:), isfSlot(:), isJBup(:), isJBdn(:)
            character(64) :: subroutine_name = "slot_CC"
        !%------------------------------------------------------------------
        !% Preliminaries
        !% --- exit if not using PS
        if (.not. setting%Solver%PreissmannSlot%useSlotTF) return
        !% --- exit if no closed CC elements

        !%------------------------------------------------------------------
        !% Aliases
        !% HACK -- many of these aliases are unused. Need to clean up 20220913 brh
            !% --- elemR data
            eHead      => elemR(:,er_Head)
            dSlotArea  => elemR(:,er_dSlotArea)
            dSlotDepth => elemR(:,er_dSlotDepth)
            fullDepth  => elemR(:,er_FullDepth)
            fullArea   => elemR(:,er_FullArea)
            PNumber    => elemR(:,er_Preissmann_Number)
            PnumberInitial => elemR(:,er_Preissmann_Number_initial)
            PCelerity  => elemR(:,er_Preissmann_Celerity)
            SlotWidth  => elemR(:,er_SlotWidth)
            SlotDepth  => elemR(:,er_SlotDepth)
            SlotDepth_N0 => elemR(:,er_SlotDepth_N0)
            SlotArea   => elemR(:,er_SlotArea)
            zCrown     => elemR(:,er_Zcrown)
            !% --- pointer to elemYN column
            isSurcharge=> elemYN(:,eYN_isSurcharged)
            isfSlot    => faceYN(:,fYN_isPSsurcharged)
            isJBup     => faceYN(:,fYN_isUpstreamJBFace)
            isJBdn     => faceYN(:,fYN_isDownstreamJBFace)
            !% --- pointers to elemI columns
            fUp        => elemI(:,ei_Mface_uL)
            fDn        => elemI(:,ei_Mface_dL)
            !% --- pointer to necessary settings struct
            SlotMethod          => setting%Solver%PreissmannSlot%Method
            TargetPCelerity     => setting%Solver%PreissmannSlot%TargetCelerity
            grav                => setting%Constant%gravity
        !%------------------------------------------------------------------
        
        !% --- Select the type of slot method
        select case (SlotMethod)
            !% --- for a static slot, no furter adjustment is required
            case (StaticSlot)
                return
            !% --- for dynamic slot, slot adjustment is required
            case (DynamicSlot,SplitDynamicSlot)

                where (elemYN(:,eYN_isSurchargeHeadAdjusted))
                
                    !% only the slot if the upstream and downstream faces are surcharged
                    ! where ((isfSlot(fUP(thisP)) .and. isfSlot(fDn(thisP))) &
                    !  .and. (.not. isJBup(fDn(thisP))) .and. (.not. isJBdn(fUp(thisP))))

                    ! where (isfSlot(fUP(thisP)) .and. isfSlot(fDn(thisP)))
                        !% --- new slot depth -- required when V-shaped surcharge changes the head
                        SlotDepth(:) = max(eHead(:) - zCrown(:), zeroR)
                        !% --- new dslot depth corresponding to the new slot depth
                        dSlotDepth(:) = SlotDepth(:) - SlotDepth_N0(:)
                        !% --- maintain the same slot area (volume/length), which implies a new slot width
                        SlotWidth(:) = abs(dSlotArea(:) / dSlotDepth(:))
                        !% --- new preissmann celerity for this slot width
                        PCelerity(:) = min(sqrt((grav * fullArea(:)) / SlotWidth(:)), TargetPCelerity)
                        !% --- new preissmann number for this celerity
                        PNumber(:) =  min(TargetPCelerity / PCelerity(:), PnumberInitial(:))
                        !% ---- recalculate the preissmann celerity
                        PCelerity(:) = TargetPCelerity / PNumber(:)

                    ! !% reset the v-shaped head adjustment where the upstream and downstream ends are not surcharged
                    ! elsewhere 
                    !     eHead(thisP) = zCrown(thisP) + SlotDepth(thisP)
                    ! end where

                end where

            case default
                !% --- should not reach this stage
                print*, 'In ', subroutine_name
                print *, 'CODE ERROR Slot Method type unknown for # ', SlotMethod
                print *, 'which has key ',trim(reverseKey(SlotMethod))
                call util_crashpoint(668732)
        end select

    end subroutine slot_Vshaped_adjust
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine slot_JM (thisCol, Npack)
        !%------------------------------------------------------------------
        !% Description:
        !% Compute Preissmann slot for closed JM's in ETM methods
        !% Note this does NOT change the volume stored.
        !%
        !% NOTE: as this is JM only, which have exactly 1 barrel, there
        !% is no need to multiply overflow/ponding by number of barrels.
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisCol, Npack
            integer, pointer    :: thisP(:), SlotMethod
            real(8), pointer    :: dSlotArea(:), dSlotDepth(:), dSlotVol(:)
            real(8), pointer    :: fullarea(:), fullVolume(:), length(:)
            real(8), pointer    :: PNumber(:), PCelerity(:), pAreaSurcharge(:)
            real(8), pointer    :: volume(:) , SlotVolume(:), SlotDepth(:), SlotArea(:)
            real(8), pointer    :: SlotWidth(:), SurchargeTime(:)
            !real(8), pointer    :: maxSlotDepth(:)
            real(8), pointer    :: SlotVolN0(:),  SlotDepthN0(:), PnumberInitial(:)
            !real(8), pointer    :: VolumeExtra(:), VolumePonded(:), VolumeOverflow(:)
            real(8), pointer    :: TargetPCelerity, cfl, grav, Dt, DecayRate ! , Alpha
            logical, pointer    :: isSlot(:), isSurcharge(:), canSurcharge(:)
            integer ::  kk, mm, JMidx, bcount 
            real(8) :: PNadd

            character(64) :: subroutine_name = 'slot_JM'
        !%-----------------------------------------------------------------
        !% Preliminaries
            if (Npack < 1) return
        !%-----------------------------------------------------------------
        !% Aliases
        !%   HACK -- many of these aliases are unused. Need to clean up
            !% --- pointer packed element indexes
            thisP => elemP(1:Npack,thisCol)
            
            !% --- pointers to elemR columns
            dSlotArea  => elemR(:,er_dSlotArea)
            dSlotDepth => elemR(:,er_dSlotDepth)
            dSlotVol   => elemR(:,er_dSlotVolume)
            fullArea   => elemR(:,er_FullArea)
            fullVolume => elemR(:,er_FullVolume)
            length     => elemR(:,er_Length)
            PNumber    => elemR(:,er_Preissmann_Number)
            PCelerity  => elemR(:,er_Preissmann_Celerity)
            SlotVolN0  => elemR(:,er_SlotVolume_N0)
            SlotVolume => elemR(:,er_SlotVolume)
            SlotDepth  => elemR(:,er_SlotDepth)
            SlotDepthN0=> elemR(:,er_SlotDepth_N0)
            SlotWidth  => elemR(:,er_SlotWidth)
            SlotArea   => elemR(:,er_SlotArea)
            volume     => elemR(:,er_Volume)
            !maxSlotDepth   => elemSR(:,esr_JM_OverflowHeightAboveCrown)
            pAreaSurcharge => elemSR(:,esr_JM_Present_PlanArea)
            PnumberInitial => elemR(:,er_Preissmann_Number_initial)
            SurchargeTime  => elemR(:,er_Surcharge_Time)
            !VolumePonded   => elemR(:,er_VolumePondedTotal)
            !VolumeExtra    => elemR(:,er_Temp01)
            !VolumeOverflow => elemR(:,er_VolumeOverFlow)
            !% --- pointer to elemYN column
            isSlot         => elemYN(:,eYN_isPSsurcharged)
            isSurcharge    => elemYN(:,eYN_isSurcharged)
            canSurcharge   => elemYN(:,eYN_canSurcharge)
            !% --- pointer to necessary settings struct
            SlotMethod          => setting%Solver%PreissmannSlot%Method
            TargetPCelerity     => setting%Solver%PreissmannSlot%TargetCelerity
            !Alpha               => setting%Solver%PreissmannSlot%Alpha
            DecayRate           => setting%Solver%PreissmannSlot%DecayRate
            cfl                 => setting%VariableDT%CFL_target
            grav                => setting%Constant%gravity
            Dt                  => setting%Time%Hydraulics%Dt
        !%-----------------------------------------------------------------
        !% --- common initialization
        SlotVolume(thisP)     = zeroR
        pAreaSurcharge(thisP) = zeroR
        SlotDepth(thisP)      = zeroR
        PCelerity(thisP)      = zeroR
        SlotWidth(thisP)      = zeroR
        !VolumeExtra(thisP)    = zeroR
        isSlot(thisP)         = .false.
        isSurcharge(thisP)    = .false.

        !% --- initialize static slot
        select case (setting%Solver%PreissmannSlot%Method)
        case (StaticSlot)
            PNumber(thisP)        = oneR  !% required for static slot
            where ((volume(thisP) >= fullVolume(thisP)) .and. canSurcharge(thisP))
                isSlot(thisP)         = .true.
                isSurcharge(thisP)    = .true.
                SlotVolume(thisP)     = max(volume(thisP) - fullVolume(thisP), zeroR)
                PCelerity(thisP)      = min(TargetPCelerity / PNumber(thisP), TargetPCelerity)
                SlotWidth(thisP)      = (grav * fullarea(thisP)) / (PCelerity(thisP) ** twoR)
                pAreaSurcharge(thisP) = SlotWidth(thisP) * elemR(thisP,er_Length) 
                SlotDepth(thisP)      = SlotVolume(thisP) / pAreaSurcharge(thisP)

            end where
        case (DynamicSlot,SplitDynamicSlot)
            !% --- requires cycling through the junctions
            do mm=1,Npack
                JMidx = thisP(mm)
                !% --- initialize dynamic slot
                dSlotVol(JMidx)          = zeroR
                dSlotArea(JMidx)         = zeroR
                dSlotDepth(JMidx)        = zeroR

                if ((volume(JMidx) > fullVolume(JMidx)) .and. canSurcharge(JMidx)) then
                    SlotVolume(JMidx) = max(volume(JMidx) - fullVolume(JMidx), zeroR)


                    SlotArea(JMidx)   = SlotVolume(JMidx) / length(JMidx)  


                    !% --- logicals
                    isSlot(JMidx)       = .true.
                    isSurcharge(JMidx)  = .true.
                    !% --- smooth out the preissmann number before celerity calculation
                    !%     with adjacent branches
                    bcount = zeroI
                    PNadd  = zeroR
                    do kk=1,max_branch_per_node
                        if (elemSI(JMidx+kk,esi_JB_Exists) .ne. oneI) cycle 

                        PNadd = PNadd + PNumber(JMidx+kk)
                        bcount = bcount + oneI
                    end do

                    PNumber(JMidx) = max(PNadd/real(bcount,8), oneR)

                    ! if (printJM == JMidx) then 
                    !     print *, ' '
                    !     print *, 'PNumber here ',PNumber(JMidx), PNadd, bcount
                    !     !print *, ' '
                    ! end if

                    PCelerity(JMidx) = min(TargetPCelerity / PNumber(JMidx), TargetPCelerity)



                    !% --- find the change in slot volume
                    dSlotVol(JMidx)   = SlotVolume(JMidx) - SlotVolN0(JMidx)
                    !% --- find the change in slot area
                    dSlotArea(JMidx)  = dSlotVol(JMidx) / length(JMidx)


                    !% --- find the change in slot depth
                    if (SlotMethod == DynamicSlot) then
                        dSlotDepth(JMidx) = (dSlotArea(JMidx) * (PCelerity(JMidx) ** twoI)) &
                                         / (grav * (fullArea(JMidx)))
                    else !% SplitDynamicSlot
                        if (dSlotVol(JMidx) > zeroR) then 
                            dSlotDepth(JMidx) = (dSlotArea(JMidx) * (PCelerity(JMidx) ** twoI)) &
                                         / (grav * (fullArea(JMidx)))
                        else 
                            dSlotDepth(JMidx) = dSlotArea(JMidx) * SlotDepth(JMidx) / SlotArea(JMidx)
                            !dSlotDepth(JMidx) = (dSlotArea(JMidx) * (PCelerity(JMidx) ** twoI)) &
                            !             / (grav * (fullArea(JMidx)))
                        end if 
                    end if


                    !% -- update the plan area for surcharging at this level
                    pAreaSurcharge(JMidx) = dSlotVol(JMidx) / dSlotDepth(JMidx)

                    SlotWidth(JMidx) = pAreaSurcharge(JMidx) / length(JMidx)

       
                    !% --- surcharge time 
                    if (bcount > zeroI) then 
                        SurchargeTime(JMidx) = SurchargeTime(JMidx) + Dt / twoR   
                    else 
                        SurchargeTime(JMidx) = zeroR
                    end if               


                    PNumber(JMidx) = (PnumberInitial(JMidx) - oneR) &
                         * exp((- SurchargeTime(JMidx) * tenR)/ DecayRate) + oneR

                else 
                    !% --- not surcharged 
                    isSlot(JMidx)      = .false. 
                    isSurcharge(JMidx) = .false.
                end if
            end do

        case default 
            print *, 'CODE ERROR unexpected case default'
            call util_crashpoint(720984)
        end select

    end subroutine slot_JM
!%
!%==========================================================================
!% END OF MODULE
!%==========================================================================
end module preissmann_slot