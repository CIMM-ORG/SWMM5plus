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
    ! subroutine slot_JM_head_PSadd (thisColP_JM)
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Adjusts head on JM for surcharge and/or ponding effects
    !     !%------------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in) :: thisColP_JM
    !         integer, pointer :: Npack, thisP(:)
    !         !real(8), pointer :: PondedArea(:), PondedVolume(:), PondedHead(:)
    !         real(8), pointer :: Head(:), SlotDepth(:)
    !     !%------------------------------------------------------------------
    !     !% Preliminaries:
    !         !% --- exit if not using PS
    !         if (.not. setting%Solver%PreissmannSlot%useSlotTF) return
    !         !% --- exit if no JM elements
    !         Npack => npack_elemP(thisColP_JM)
    !         if (Npack < 1) return
    !     !%------------------------------------------------------------------
    !     !% Aliases
    !         thisP        => elemP(1:Npack,thisColP_JM)
    !         !PondedArea   => elemSR(:,esr_JunctionMain_PondedArea)
    !         !PondedVolume =>  elemR(:,er_VolumePonded)
    !         !PondedHead   =>  elemR(:,er_Temp01)
    !         Head         =>  elemR(:,er_Head)
    !         SlotDepth    =>  elemR(:,er_SlotDepth)

    !     !%------------------------------------------------------------------

    !     print *, 'OBSOLETE 20230628'
    !     stop 669873
    !     Head(thisP) = Head(thisP) + SlotDepth(thisP)

    ! end subroutine slot_JM_head_PSadd
!%   
!%==========================================================================
!%==========================================================================
!%
    ! subroutine slot_JM_head_PSremove (thisColP_JM)
    !     !%------------------------------------------------------------------
    !     !% Description:
    !     !% Adjusts head on JM for surcharge and/or ponding effects
    !     !%------------------------------------------------------------------
    !     !% Declarations:
    !         integer, intent(in) :: thisColP_JM
    !         integer, pointer :: Npack, thisP(:)
    !         !real(8), pointer :: PondedArea(:), PondedVolume(:), PondedHead(:)
    !         real(8), pointer :: Head(:), SlotDepth(:)
    !     !%------------------------------------------------------------------
    !     !% Preliminaries:
    !         !% --- exit if not using PS
    !         if (.not. setting%Solver%PreissmannSlot%useSlotTF) return
    !         Npack => npack_elemP(thisColP_JM)
    !         if (Npack < 1) return
    !     !%------------------------------------------------------------------
    !     !% Aliases:
    !         thisP        => elemP(1:Npack,thisColP_JM)
    !         !PondedArea   => elemSR(:,esr_JunctionMain_PondedArea)
    !         !PondedVolume => elemR (:,er_VolumePonded)
    !         !PondedHead   => elemR (:,er_Temp01)
    !         Head         => elemR (:,er_Head)
    !         SlotDepth    => elemR (:,er_SlotDepth)
    !     !%------------------------------------------------------------------

    !     print *, 'OBSOLETE 20230628'
    !     stop 77222987

    !     Head(thisP) = Head(thisP) - SlotDepth(thisP)
 

    ! end subroutine slot_JM_head_PSremove
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
        real(8), pointer    :: SlotWidth(:), SlotVolume(:), SlotDepth(:), dSlotDepth(:)
        real(8), pointer    :: volume(:), depth(:), area(:), SlotArea(:) !, ell(:)
        real(8), pointer    :: SlotDepth_N0(:)
        real(8), pointer    :: head(:), fullVolume(:), fullArea(:), fullDepth(:)
        real(8), pointer    :: Overflow(:), zbottom(:), SlotHydRad(:)!, ellMax(:)
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
            area       => elemR(:,er_Area)
            depth      => elemR(:,er_Depth)
            dSlotDepth => elemR(:,er_dSlotDepth)
            fullDepth  => elemR(:,er_FullDepth)
            fullvolume => elemR(:,er_FullVolume)
            fullArea   => elemR(:,er_FullArea)
            head       => elemR(:,er_Head)
            SlotWidth  => elemR(:,er_SlotWidth)
            SlotVolume => elemR(:,er_SlotVolume)
            SlotDepth  => elemR(:,er_SlotDepth)
            SlotDepth_N0 => elemR(:,er_SlotDepth_N0)
            SlotArea   => elemR(:,er_SlotArea)
            SlotHydRad => elemR(:,er_SlotHydRadius)
            volume     => elemR(:,er_Volume)
            zbottom    => elemR(:,er_Zbottom)
            isSlot     => elemYN(:,eYN_isPSsurcharged)

            !% pointer to necessary settings struct
            SlotMethod => setting%Solver%PreissmannSlot%Method
        !%-----------------------------------------------------------------------------

        !% CC slot adjustment
        select case (SlotMethod)
            case (StaticSlot)
                where (isSlot(thisP)) 
                    volume(thisP)    = volume(thisP)  + SlotVolume(thisP)
                    head(thisP)      = zbottom(thisP) + fullDepth(thisP) + SlotDepth(thisP)
                end where 

            case (DynamicSlot)
                where (isSlot(thisP)) 
                    volume(thisP)    = volume(thisP)  + SlotVolume(thisP)
                    SlotDepth(thisP) = max(SlotDepth_N0(thisP) + dSlotDepth(thisP), zeroR) 
                    head(thisP)      = max(zbottom(thisP) + fullDepth(thisP) + SlotDepth(thisP), zbottom(thisP))
                elsewhere
                    SlotDepth(thisP)  = zeroR
                end where 
            case default 
                print *, 'CODE ERROR: unexpected case default'
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
            real(8), pointer    :: dSlotDepth(:)
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
            thisP    => elemP(1:Npack,thisCol)
            volume   => elemR(:,er_Volume)
            depth    => elemR(:,er_Depth)
            dSlotDepth => elemR(:,er_dSlotDepth)
            fullDepth  => elemR(:,er_FullDepth)
            head     => elemR(:,er_Head)
            ellDepth => elemR(:,er_EllDepth)
            SlotVolume => elemR(:,er_SlotVolume)
            SlotDepth  => elemR(:,er_SlotDepth)
            SlotDepth_N0 => elemR(:,er_SlotDepth_N0)
            isSlot     => elemYN(:,eYN_isPSsurcharged)
            zbottom    => elemR(:,er_Zbottom)
            !% pointer to necessary settings struct
            SlotMethod => setting%Solver%PreissmannSlot%Method
        !%------------------------------------------------------------------

        select case (SlotMethod)
            case (StaticSlot)
                where (isSlot(thisP)) 
                    volume(thisP)    = volume(thisP)  + SlotVolume(thisP)
                    head(thisP)      = zbottom(thisP) + fullDepth(thisP) + SlotDepth(thisP)
                end where 

            case (DynamicSlot)
                where (isSlot(thisP)) 
                    SlotDepth(thisP) = max(SlotDepth_N0(thisP) + dSlotDepth(thisP), zeroR) 
                    head(thisP)      = max(zbottom(thisP) + fullDepth(thisP) + SlotDepth(thisP), zbottom(thisP))
                elsewhere
                    SlotDepth(thisP)  = zeroR
                end where 
            case default 
                print *, 'CODE ERROR: unexpected case default'
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
            real(8), pointer :: area(:), depth(:), head(:), length(:), volume(:), zcrown(:), zbottom(:)
            real(8), pointer :: fullDepth(:), fullArea(:), fPNumber(:), PNumber(:), PCelerity(:), ellDepth(:)
            real(8), pointer :: SlotWidth(:), SlotVolume(:), SlotDepth(:), SlotArea(:)!, ellMax(:)
            real(8), pointer :: overflow(:), grav, TargetPCelerity, Alpha
            logical, pointer :: isSlot(:) , fSlot(:), isSurcharge(:), canSurcharge(:)
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
            area          => elemR(:,er_Area)
            depth         => elemR(:,er_Depth)
            head          => elemR(:,er_Head)
            length        => elemR(:,er_Length)
            fullArea      => elemR(:,er_FullArea)
            fullDepth     => elemR(:,er_FullDepth)
            overflow      => elemR(:,er_VolumeOverFlow)
            volume        => elemR(:,er_Volume)
            zcrown        => elemR(:,er_Zcrown)
            zbottom       => elemR(:,er_Zbottom)
            ellDepth      => elemR(:,er_EllDepth)
            fUp           => elemI(:,ei_Mface_uL)
            fDn           => elemI(:,ei_Mface_dL)
            BranchExists  => elemSI(:,esi_JunctionBranch_Exists)
            canSurcharge  => elemYN(:,eYN_canSurcharge)
            grav          => setting%Constant%gravity
        !% Slot Aliases
            PNumber    => elemR(:,er_Preissmann_Number)
            PCelerity  => elemR(:,er_Preissmann_Celerity)
            SlotWidth  => elemR(:,er_SlotWidth)
            SlotVolume => elemR(:,er_SlotVolume)
            SlotDepth  => elemR(:,er_SlotDepth)
            SlotArea   => elemR(:,er_SlotArea)
            fPNumber   => faceR(:,fr_Preissmann_Number)
            isSurcharge=> elemYN(:,eYN_isSurcharged)
            isSlot     => elemYN(:,eYN_isPSsurcharged)
            fSlot      => faceYN(:,fYN_isPSsurcharged)
            SlotMethod      => setting%Solver%PreissmannSlot%Method
            TargetPCelerity => setting%Solver%PreissmannSlot%TargetCelerity
            Alpha           => setting%Solver%PreissmannSlot%Alpha
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
                    SlotWidth(tB)  = zeroR
                    SlotVolume(tB) = zeroR
                    PCelerity(tB)  = zeroR

                    !% --- HACK -- may want to define explicitly each JB
                    !%     as either closed/open  20220915 brh

                    !% --- assuming a slot if the head is above the crown
                    !%     or the upstream CC is in a slot
                    if ((head(tB) .gt. zcrown(tB)) .and. (canSurcharge(tB))) then
                        isSlot(tB)     = .true.
                        fSlot(fUp(tB)) = .true.
                        isSurcharge(tB)= .true.
                        PNumber(tB)    = fPNumber(fUp(tB))
                        PCelerity(tB)  = min(TargetPCelerity / PNumber(tB), TargetPCelerity)
                        SlotDepth(tB)  = max(head(tB) - zcrown(tB), zeroR)   
                        SlotArea(tB)   = (SlotDepth(tB) * (PNumber(tB)**twoR) * grav * &
                                            fullArea(tB)) / (TargetPCelerity ** twoR)
                        SlotVolume(tB) = SlotArea(tB) * length(tB)
                        
                        !% --- add the slot geometry back to previously solved geometry
                        volume(tB) = volume(tB)  + SlotVolume(tB)
                        Overflow(tB) = zeroR  !% defined as zero (no oveflow from JB)
                    else 
                        !% --- excess head at JB in an open channel represents
                        !%     head from ponding of JM, so do not adjust.
                        !%     Volume and depth are retained at full levels
                    end if  
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

                    !% --- assuming a slot if the head is above the crown
                    !%     or the downstream CC is in a slot
                    if ((head(tB) .gt. zcrown(tB)) .and. (canSurcharge(tB))) then
                        isSlot(tB)     = .true.
                        isSurcharge(tB)= .true.
                        fSlot(fDn(tB)) = .true.
                        PNumber(tB)    = fPNumber(fDn(tB))
                        PCelerity(tB)  = min(TargetPCelerity / PNumber(tB), TargetPCelerity)

                        SlotDepth(tB)  = max(head(tB) - zcrown(tB), zeroR)    

                        SlotArea(tB)   = (SlotDepth(tB) * (PNumber(tB)**twoR) * grav * &
                                            fullArea(tB)) / (TargetPCelerity ** twoR)
                        SlotVolume(tB) = SlotArea(tB) * length(tB)

                        !% --- add the slot geometry back to previously solved geometry
                        volume(tB) = volume(tB)  + SlotVolume(tB)
                        Overflow(tB) = zeroR !% defined as zero (no overflow from JB)
                    else 
                        !% --- excess head at JB in an open channel represents
                        !%     head from ponding of JM, so do not adjust level.
                        !%     Volume and depth are retained at full levels
                    end if
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
    subroutine slot_CC (thisP)
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
            integer, pointer    :: SlotMethod, fUp(:), fDn(:)
            real(8), pointer    :: fullarea(:), PNumberOld(:) !, ellMax(:)
            real(8), pointer    :: fullVolume(:), length(:), PNumber(:), PCelerity(:) 
            real(8), pointer    :: SlotWidth(:), SlotVolume(:), SlotDepth(:), SlotArea(:)
            real(8), pointer    :: dSlotVol(:), dSlotArea(:), dSlotDepth(:), SlotVolN0(:), volume(:) 
            real(8), pointer    :: velocity(:), fPNumber(:), SurchargeTime(:), PnumberInitial(:)
            real(8), pointer    :: TargetPCelerity, grav, Dt, Alpha, cfl, DecayRate
            logical, pointer    :: isSlot(:), isfSlot(:), isSurcharge(:)
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
            SlotArea   => elemR(:,er_SlotArea)
            volume     => elemR(:,er_Volume)
            velocity   => elemR(:,er_velocity)
            SurchargeTime  => elemR(:,er_Surcharge_Time)
            PnumberInitial => elemR(:,er_Preissmann_Number_initial)
            !% --- pointer to elemYN column
            isSurcharge=> elemYN(:,eYN_isSurcharged)
            isSlot     => elemYN(:,eYN_isPSsurcharged)
            isfSlot    => faceYN(:,fYN_isPSsurcharged)
            !% --- pointers to elemI columns
            fUp        => elemI(:,ei_Mface_uL)
            fDn        => elemI(:,ei_Mface_dL)
            !% --- pointer to faceR column
            fPNumber   => faceR(:,fr_Preissmann_Number)
            !% --- pointer to necessary settings struct
            SlotMethod          => setting%Solver%PreissmannSlot%Method
            TargetPCelerity     => setting%Solver%PreissmannSlot%TargetCelerity
            DecayRate           => setting%Solver%PreissmannSlot%DecayRate
            Alpha               => setting%Solver%PreissmannSlot%Alpha
            cfl                 => setting%VariableDT%CFL_target
            grav                => setting%Constant%gravity
            Dt                  => setting%Time%Hydraulics%Dt
        !%------------------------------------------------------------------
        !% --- common initialization
        SlotVolume(thisP)   = zeroR
        SlotArea(thisP)     = zeroR
        PCelerity(thisP)    = zeroR
        isSlot(thisP)       = .false.
        isSurcharge(thisP)  = .false.
        isfSlot(fUp(thisP)) = .false.
        isfSlot(fDn(thisP)) = .false.

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
                end where
            
            !% --- for dynamic slot, preissmann number is adjusted
            case (DynamicSlot)

                !% --- initialize dynamic slot
                dSlotVol          = zeroR
                dSlotArea         = zeroR
                dSlotDepth        = zeroR

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

                    !% --- find the change in slot depth
                    dSlotDepth(thisP) = (dSlotArea(thisP)  * (PCelerity(thisP) ** twoI)) / (grav * (fullArea(thisP)))
                end where

                !% --- reset isfSlot to find ventilated positions
                where (.not. isSlot(thisP))
                    isfSlot(fUp(thisP)) = .false.
                    isfSlot(fDn(thisP)) = .false.
                end where

                !% --- set any diagnostic face as venting point
                !%     thus not changing preissmann number here
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
            real(8), pointer    :: TargetPCelerity, cfl, grav, Alpha, Dt, DecayRate
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
            !maxSlotDepth   => elemSR(:,esr_JunctionMain_OverflowHeightAboveCrown)
            pAreaSurcharge => elemSR(:,esr_JunctionMain_Surcharge_Plan_Area)
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
            Alpha               => setting%Solver%PreissmannSlot%Alpha
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
                !% obsolete below 20230628 brh
                !% --- difference between SlotDepth and maximum is the "extra" volume
                !%     If negative, this is volume that could be filled by ponding
                !%     If positive, this is volume that is added to ponding or lost 
                !%     in overflow
                ! VolumeExtra(thisP) =  (SlotDepth(thisP) - maxSlotDepth(thisP)) &
                !                         * pAreaSurcharge(thisP)
            end where
        case (DynamicSlot)
            !% --- initialize dynamic slot
            dSlotVol          = zeroR
            dSlotArea         = zeroR
            dSlotDepth        = zeroR
            !% --- requires cycling through the junctions
            do mm=1,Npack
                JMidx = thisP(mm)

                if (volume(JMidx) > fullVolume(JMidx)) then
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
                        if (elemSI(JMidx+kk,esi_JunctionBranch_Exists) .ne. oneI) cycle 

                        PNadd = PNadd + PNumber(JMidx+kk)
                        bcount = bcount + oneI
                    end do

                    PNumber(JMidx) = max(PNadd/real(bcount,8), oneR)
                    PCelerity(JMidx) = min(TargetPCelerity / PNumber(JMidx), TargetPCelerity)

                    !% --- find the change in slot volume
                    dSlotVol(JMidx)   = SlotVolume(JMidx) - SlotVolN0(JMidx)
                    !% --- find the change in slot area
                    dSlotArea(JMidx)  = dSlotVol(JMidx) / length(JMidx)

                    !% --- find the change in slot depth
                    dSlotDepth(JMidx) = (dSlotArea(JMidx) * (PCelerity(JMidx) ** twoI)) &
                                         / (grav * (fullArea(JMidx)))


                    !% --- for case where change in slot depth is negative while the slot volume
                    !%     is still positive
                    if (dSlotDepth(JMidx) < SlotDepthN0(JMidx)) then 
                        dSlotDepth(JMidx) = (SlotVolume(JMidx) / elemSR(JMidx,esr_Storage_Plan_Area)) &
                                           - SlotDepthN0(JMidx)    
                    end if                 

                    !% -- update the plan area for surcharging at this level
                    pAreaSurcharge(JMidx) = dSlotVol(JMidx) / dSlotDepth(JMidx)
       
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
            print *, 'CODE ERROR: unexpected case default'
            call util_crashpoint(720984)
        end select

    end subroutine slot_JM
!%
!%==========================================================================
!% END OF MODULE
!%==========================================================================
end module preissmann_slot