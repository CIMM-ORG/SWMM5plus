module preissmann_slot

    use define_globals
    use define_indexes
    use define_keys
    use define_settings, only: setting
    use utility_crash
    use utility, only: util_CLprint


    implicit none

!%-----------------------------------------------------------------------------
!% Description:
!% Preissmann Slot computations for explicit time-marching of surcharged pipe
!%

    private
    public :: slot_toplevel 
    public :: slot_JM_head_PSadd
    public :: slot_JM_head_PSremove
    public :: slot_CC_adjustments
    public :: slot_JM_adjustments
    public :: slot_JB_computation

    contains
!%    
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine slot_toplevel (whichTM, colP_CC_PS, colP_JM_PS)
        !%------------------------------------------------------------------
        !% Description:
        !% Stores the slot volumes etc. and handles ponding/overflow
        !% associated with JM that are above maximum depth for Preissman Slot
        !% NOTE: need to ensure that "whichTM" is consistent with the packed
        !% columns of colP_CC_PS and colP_JM_PS
        !% input: thisSolver is the solver used for this set of elements
        !% whcih should be ETM or ALLtm.
        !% Pcol_CC_PS is the packed column of CC (closed only) elements
        !% where the Preissmann Slot is used
        !% Pcol_JM_PS is the packed column for JM elements where Preissmann
        !% slot is used.
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: whichTM, colP_CC_PS, colP_JM_PS
            integer, pointer    :: thisPackCol, Npack
        !%------------------------------------------------------------------
        !% Preliminaries
            !% --- exit if not using PS
            if (.not. setting%Solver%PreissmannSlot%useSlotTF) return
            !% --- compute Preissmann slot for conduits only if ETM solver is used
            if (.not. ((whichTM .eq. ETM) .or. (whichTM .eq. ALLtm))) return
        !%------------------------------------------------------------------
        !%    
        !% --- Handle Preissmann Slot for closed CC elements
        !%     with this time march type.
        thisPackCol => col_elemP(colP_CC_PS)
        Npack       => npack_elemP(thisPackCol)
        if (Npack > 0) then
            call slot_CC_ETM (thisPackCol, Npack)
        end if

        !% --- Handle Preissmann Slot for all JM elements
        !%    (although SurchargeExtraDepth=0 is effectively open)
        thisPackCol => col_elemP(colP_JM_PS)
        Npack => npack_elemP(thisPackCol)
        if (Npack > 0) then
            call slot_JM_ETM (thisPackCol, Npack)
        end if
  
    end subroutine slot_toplevel
!%
!%==========================================================================
!%==========================================================================
!%    
    subroutine slot_JM_head_PSadd (thisColP_JM)
        !%------------------------------------------------------------------
        !% Description:
        !% Adjusts head on JM for surcharge and/or ponding effects
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisColP_JM
            integer, pointer :: Npack, thisP(:)
            real(8), pointer :: PondedArea(:), PondedVolume(:), PondedHead(:)
            real(8), pointer :: Head(:), SlotDepth(:)

        !%------------------------------------------------------------------
        !% Preliminaries:
            !% --- exit if not using PS
            if (.not. setting%Solver%PreissmannSlot%useSlotTF) return
            !% --- exit if no JM elements
            Npack => npack_elemP(thisColP_JM)
            if (Npack < 1) return
        !%------------------------------------------------------------------
        !% Aliases
            thisP        => elemP(1:Npack,thisColP_JM)
            PondedArea   => elemSR(:,esr_JunctionMain_PondedArea)
            PondedVolume => elemSR(:,esr_JunctionMain_PondedVolume)
            PondedHead   =>  elemR(:,er_Temp01)
            Head         =>  elemR(:,er_Head)
            SlotDepth    =>  elemR(:,er_SlotDepth)

        !%------------------------------------------------------------------
        PondedHead = zeroR
        if (setting%SWMMinput%AllowPonding) then 
            PondedHead(thisP) = PondedVolume(thisP) / PondedArea(thisP)
            Head(thisP) = Head(thisP) + SlotDepth(thisP) + PondedHead(thisP)
        else
            Head(thisP) = Head(thisP) + SlotDepth(thisP)
        end if

    end subroutine slot_JM_head_PSadd
!%   
!%==========================================================================
!%==========================================================================
!%
    subroutine slot_JM_head_PSremove (thisColP_JM)
        !%------------------------------------------------------------------
        !% Description:
        !% Adjusts head on JM for surcharge and/or ponding effects
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisColP_JM
            integer, pointer :: Npack, thisP(:)
            real(8), pointer :: PondedArea(:), PondedVolume(:), PondedHead(:)
            real(8), pointer :: Head(:), SlotDepth(:)
        !%------------------------------------------------------------------
        !% Preliminaries:
            !% --- exit if not using PS
            if (.not. setting%Solver%PreissmannSlot%useSlotTF) return
            Npack => npack_elemP(thisColP_JM)
            if (Npack < 1) return
        !%------------------------------------------------------------------
        !% Aliases
            PondedArea   => elemSR(:,esr_JunctionMain_PondedArea)
            PondedVolume => elemSR(:,esr_JunctionMain_PondedVolume)
            PondedHead   =>  elemR(:,er_Temp01)
            Head         =>  elemR(:,er_Head)
            SlotDepth    =>  elemR(:,er_SlotDepth)
        !%------------------------------------------------------------------
        if (setting%SWMMinput%AllowPonding) then 
            PondedHead(thisP) = PondedVolume(thisP) / PondedArea(thisP)
            Head(thisP)       = Head(thisP) - SlotDepth(thisP) - PondedHead(thisP)
        else
            Head(thisP) = Head(thisP) - SlotDepth(thisP)
        end if

    end subroutine slot_JM_head_PSremove
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine slot_CC_adjustments (thisColP_closed_CC)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% This subroutine adds back the slot geometry in all the closed elements
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: thisColP_closed_CC
        integer, pointer    :: thisP(:), Npack, SlotMethod
        real(8), pointer    :: SlotWidth(:), SlotVolume(:), SlotDepth(:), dSlotDepth(:)
        real(8), pointer    :: volume(:), ell(:), depth(:), area(:), SlotArea(:)
        real(8), pointer    :: head(:), fullVolume(:), fullArea(:), fullDepth(:)
        real(8), pointer    :: Overflow(:), zbottom(:), ellMax(:), SlotHydRad(:)
        logical, pointer    :: isSlot(:)

        character(64) :: subroutine_name = 'geo_CC_slot_adjustments'
        !%-----------------------------------------------------------------------------
        !%------------------------------------------------------------------
        !% Preliminaries
            !% --- exit if not using PS
            if (.not. setting%Solver%PreissmannSlot%useSlotTF) return

            Npack => npack_elemP(thisColP_closed_CC)
            if (Npack < 1) return
            if (setting%Debug%File%geometry) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% Aliases
            thisP      => elemP(1:Npack,thisColP_closed_CC)
            area       => elemR(:,er_Area)
            depth      => elemR(:,er_Depth)
            dSlotDepth => elemR(:,er_dSlotDepth)
            ell        => elemR(:,er_ell)
            ellMax     => elemR(:,er_ell_max)
            fullDepth  => elemR(:,er_FullDepth)
            fullvolume => elemR(:,er_FullVolume)
            fullArea   => elemR(:,er_FullArea)
            head       => elemR(:,er_Head)
            !Overflow   => elemR(:,er_VolumeOverFlow) !% handled elsewhere
            SlotWidth  => elemR(:,er_SlotWidth)
            SlotVolume => elemR(:,er_SlotVolume)
            SlotDepth  => elemR(:,er_SlotDepth)
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
                    volume(thisP) = volume(thisP)  + SlotVolume(thisP)
                    ! area(thisP)   = area(thisP)    + SlotArea(thisP) !% KEEP AREA BASED ON CONDUIT
                    depth(thisP)  = depth(thisP)   + SlotDepth(thisP)
                    head(thisP)   = head(thisP)    + SlotDepth(thisP)
                    ell(thisP)    = ellMax(thisP)  !% IS THIS NEEDED, WHY?  20220915 brh
                    !Overflow(thisP) = zeroR !% Set elsewhsere
                end where 

            case (DynamicSlot)
                where (isSlot(thisP)) 
                    volume(thisP)    = volume(thisP)  + SlotVolume(thisP)
                    ! area(thisP)      = area(thisP)    + SlotArea(thisP)  !% KEEP AREA BASED ON CONDUIT
                    SlotDepth(thisP) = max(SlotDepth(thisP) + dSlotDepth(thisP), zeroR) 
                    head(thisP)      = max(head(thisP)  + SlotDepth(thisP), zbottom(thisP))
                    !% 20220915 brh WHY ISN'T the SLOTDEPTH ADDED TO DEPTH?
                    !Overflow(thisP)  = zeroR !% set elsewhere
                elsewhere
                    SlotDepth(thisP)  = zeroR
                    ! SlotVolume(thisP) = zeroR
                end where 
        end select


        if (setting%Debug%File%geometry) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine slot_CC_adjustments
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine slot_JM_adjustments (thisColP_JM)    
        !%------------------------------------------------------------------
        !% Description:
        !% Adjusts geometry in JM for Preissmann Slot
        !% Uses the Static Slot approach (only)
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisColP_JM
            integer, pointer    :: Npack, thisP(:)
            real(8), pointer    :: volume(:), depth(:), ell(:), ellMax(:)
            real(8), pointer    :: SlotVolume(:), SlotDepth(:)
        !%------------------------------------------------------------------
        !% Preliminaries:
            !% --- exit if not using PS
            if (.not. setting%Solver%PreissmannSlot%useSlotTF) return
            Npack => npack_elemP(thisColP_JM)
            if (Npack < 1) return
        !%------------------------------------------------------------------
        !% Aliases
            thisP => elemP(1:Npack,thisColP_JM)
            volume => elemR(:,er_Volume)
            depth  => elemR(:,er_Depth)
            ell    => elemR(:,er_ell)
            ellMax => elemR(:,er_ell_max) 
        !%------------------------------------------------------------------

        call slot_JM_head_PSadd (thisColP_JM)

        volume(thisp) = volume(thisP) + SlotVolume(thisP)
        depth(thisP)  = depth(thisP)  + SlotDepth(thisP)
        ell(thisP)    = ellMax(thisP)


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
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisColP_JM
            integer, pointer :: Npack, thisP(:), tM, BranchExists(:)
            real(8), pointer :: area(:), depth(:), head(:), length(:), volume(:), zcrown(:), zbottom(:)
            real(8), pointer :: fullDepth(:), fullArea(:), fPNumber(:), PNumber(:), PCelerity(:), ell(:)
            real(8), pointer :: SlotWidth(:), SlotVolume(:), SlotDepth(:), SlotArea(:), ellMax(:)
            real(8), pointer :: overflow(:), grav, TargetPCelerity, Alpha
            logical, pointer :: isSlot(:) , fSlot(:), isDnJB(:)
            integer, pointer :: SlotMethod, fUp(:), fDn(:)
            integer :: tB, ii, kk
        !%------------------------------------------------------------------
        !% Preliminaries:
            !% --- exit if not using PS
            if (.not. setting%Solver%PreissmannSlot%useSlotTF) return
            Npack => npack_elemP(thisColP_JM)
            if (Npack < 1) return
        !%------------------------------------------------------------------
        !% Aliases
            thisP        => elemP(1:Npack,thisColP_JM)
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
            ell           => elemR(:,er_ell)
            ellMax        => elemR(:,er_ell_max)
            fUp           => elemI(:,ei_Mface_uL)
            fDn           => elemI(:,ei_Mface_dL)
            BranchExists  => elemSI(:,esi_JunctionBranch_Exists)
            grav          => setting%Constant%gravity
        !% Slot Aliases
            PNumber    => elemR(:,er_Preissmann_Number)
            PCelerity  => elemR(:,er_Preissmann_Celerity)
            SlotWidth  => elemR(:,er_SlotWidth)
            SlotVolume => elemR(:,er_SlotVolume)
            SlotDepth  => elemR(:,er_SlotDepth)
            SlotArea   => elemR(:,er_SlotArea)
            fPNumber   => faceR(:,fr_Preissmann_Number)
            isSlot     => elemYN(:,eYN_isPSsurcharged)
            isDnJB     => elemYN(:,eYN_isDownstreamJB)
            fSlot      => faceYN(:,fYN_isPSsurcharged)
            SlotMethod      => setting%Solver%PreissmannSlot%Method
            TargetPCelerity => setting%Solver%PreissmannSlot%TargetCelerity
            Alpha           => setting%Solver%PreissmannSlot%Alpha
        !%------------------------------------------------------------------
        !% JB slot adjustment
        !% cycle through the all the main junctions and each of its branches
        do ii=1,Npack
            tM => thisP(ii) !% junction main ID
            ! handle the upstream branches
            do kk=1,max_branch_per_node,2
                tB = tM + kk  !% JB branch ID
                if (BranchExists(tB)==1) then
                    !% initialize slot
                    isSlot(tB)     = .false.
                    SlotDepth(tB)  = zeroR
                    SlotArea(tB)   = zeroR
                    SlotWidth(tB)  = zeroR
                    SlotVolume(tB) = zeroR
                    PCelerity(tB)  = zeroR

                    !% HACK -- may want to define explicitly each JB
                    !% as either closed/open  20220915 brh

                    !% assuming a slot if the head is above the crown
                    !% or the upstream CC is in a slot
                    if (head(tB) .gt. zcrown(tB)) then
                        isSlot(tB)     = .true.
                        fSlot(fUp(tB)) = .true.
                        PNumber(tB)    = fPNumber(fUp(tB))
                        PCelerity(tB)  = min(TargetPCelerity / PNumber(tB), TargetPCelerity)

                !% ERROR? the JB depth has already been limited to the full depth
                !% but the elem(:,er_Head) for JB should have the correct surcharged
                !% value.

                        SlotDepth(tB)  = max(depth(tB) - fulldepth(tB), zeroR)   

                        SlotArea(tB)   = (SlotDepth(tB) * (PNumber(tB)**twoR) * grav * &
                                            fullArea(tB)) / (TargetPCelerity ** twoR)
                        SlotVolume(tB) = SlotArea(tB) * length(tB)
                        
                        !% add the slot geometry back to previously solved geometry
                        volume(tB) = volume(tB)  + SlotVolume(tB)
                        !area(tB)   = area(tB)    + SlotArea(tB) !% 20220915 brh CONSISTENCY WITH CC adjustment
                        depth(tB)  = depth(tB)   + SlotDepth(tB)
                        Overflow(tB) = zeroR  !% defined as zero (no oveflow from JB)
                    end if  
                end if
            end do
            !% handle the downstream branches
            do kk=2,max_branch_per_node,2
                tB = tM + kk
                if (BranchExists(tB)==1) then
                    !% initialize slot
                    isSlot(tB)     = .false.
                    SlotDepth(tB)  = zeroR
                    SlotArea(tB)   = zeroR
                    SlotWidth(tB)  = zeroR
                    SlotVolume(tB) = zeroR
                    PCelerity(tB)  = zeroR

                    !% HACK -- may want to define explicitly each JB
                    !% as either closed/open  20220915 brh

                    !% assuming a slot if the head is above the crown
                    !% or the downstream CC is in a slot
                    if (head(tB) .gt. zcrown(tB)) then
                        isSlot(tB)     = .true.
                        fSlot(fDn(tB)) = .true.
                        PNumber(tB)    = fPNumber(fDn(tB))
                        PCelerity(tB)  = min(TargetPCelerity / PNumber(tB), TargetPCelerity)

                !% ERROR? the JB depth has already been limited to the full depth
                !% but the elem(:,er_Head) for JB should have the correct surcharged
                !% value.
                        SlotDepth(tB)  = max(depth(tB) - fulldepth(tB), zeroR)    

                        SlotArea(tB)   = (SlotDepth(tB) * (PNumber(tB)**twoR) * grav * &
                                            fullArea(tB)) / (TargetPCelerity ** twoR)
                        SlotVolume(tB) = SlotArea(tB) * length(tB)

                        !% add the slot geometry back to previously solved geometry
                        volume(tB) = volume(tB)  + SlotVolume(tB)
                        !area(tB)   = area(tB)    + SlotArea(tB) !% 20220915 brh CONSISTENCY WITH CC adjustment
                        depth(tB)  = depth(tB)   + SlotDepth(tB)
                        Overflow(tB) = zeroR !% defined as zero (no overflow from JB)
                    end if
                end if
            end do
        end do
                  
        !%------------------------------------------------------------------
        !% Closing
    end subroutine slot_JB_computation
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%    
    subroutine slot_CC_ETM (thisCol, Npack)
        !%------------------------------------------------------------------
        !% Description:
        !% Compute Preissmann slot for conduits in ETM methods
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisCol, Npack
            integer, pointer    :: thisP(:), SlotMethod, fUp(:), fDn(:)
            real(8), pointer    :: AreaN0(:), BreadthMax(:), ellMax(:), fullarea(:)
            real(8), pointer    :: fullVolume(:), length(:), PNumber(:), PCelerity(:) 
            real(8), pointer    :: SlotWidth(:), SlotVolume(:), SlotDepth(:), SlotArea(:), temp(:)
            real(8), pointer    :: dSlotVol(:), dSlotArea(:), dSlotDepth(:), oldSlotVOl(:), volume(:) 
            real(8), pointer    :: velocity(:), fPNumber(:)
            real(8), pointer    :: TargetPCelerity, grav, Dt, Alpha, cfl
            logical, pointer    :: isSlot(:), isfSlot(:)
            character(64) :: subroutine_name = "slot_CC_ETM"
        !%------------------------------------------------------------------
        !% Preliminaries
            !% --- exit if not using PS
            if (.not. setting%Solver%PreissmannSlot%useSlotTF) return
            !% --- exit if no closed CC elements
            if (Npack < 1) return
        !%------------------------------------------------------------------
        !% Aliases
        !% HACK -- many of these aliases are unused. Need to clean up 20220913 brh
            thisP      => elemP(1:Npack,thisCol)
            !% --- elemR data
            AreaN0     => elemR(:,er_Area_N0)
            BreadthMax => elemR(:,er_BreadthMax)
            dSlotArea  => elemR(:,er_dSlotArea)
            dSlotDepth => elemR(:,er_dSlotDepth)
            dSlotVol   => elemR(:,er_dSlotVolume)
            ellMax     => elemR(:,er_ell_max)
            fullArea   => elemR(:,er_FullArea)
            fullVolume => elemR(:,er_FullVolume)
            length     => elemR(:,er_Length)
            oldSlotVol => elemR(:,er_SlotVolumeOld)
            PNumber    => elemR(:,er_Preissmann_Number)
            PCelerity  => elemR(:,er_Preissmann_Celerity)
            SlotWidth  => elemR(:,er_SlotWidth)
            SlotVolume => elemR(:,er_SlotVolume)
            SlotDepth  => elemR(:,er_SlotDepth)
            SlotArea   => elemR(:,er_SlotArea)
            volume     => elemR(:,er_Volume)
            velocity   => elemR(:,er_velocity)
            temp       => elemR(:,er_Temp01)
            !% --- pointer to elemYN column
            isSlot     => elemYN(:,eYN_isPSsurcharged)
            isfSlot    => faceYN(:,fYN_isPSsurcharged)
            !% --- pointers to elemI columns
            fUp        => elemI(:,ei_Mface_uL)
            fDn        => elemI(:,ei_Mface_dL)
            !% --- pointer to faceR column
            fPNumber   => faceR(:,fr_Preissmann_Number)
            !% --- pointer to necessary settings struct
            SlotMethod          => setting%PreissmannSlot%Method
            TargetPCelerity     => setting%PreissmannSlot%TargetCelerity
            Alpha               => setting%PreissmannSlot%Alpha
            cfl                 => setting%VariableDT%CFL_target
            grav                => setting%Constant%gravity
            Dt                  => setting%Time%Hydraulics%Dt
        !%------------------------------------------------------------------
        !% --- common initialization
        SlotVolume(thisP)   = zeroR
        SlotArea(thisP)     = zeroR
        PCelerity(thisP)    = zeroR
        isSlot(thisP)       = .false.
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
                    isSlot(thisP)     = .true.
                    PCelerity(thisP)  = min(TargetPCelerity / PNumber(thisP), TargetPCelerity)
                    SlotWidth(thisP)  = (grav * fullarea(thisP)) / (PCelerity(thisP) ** twoR)
                    SlotDepth(thisP)  = SlotArea(thisP) / SlotWidth(thisP)
                end where
            
            !% --- for dynamic slot, preissmann number is adjusted
            case (DynamicSlot)
                !% --- initialize dynamic slot
                oldSlotVol(thisP) = SlotVolume(thisP)
                dSlotVol          = zeroR
                dSlotArea         = zeroR
                dSlotDepth        = zeroR
                temp              = zeroR
                !% ---find out the slot volume/ area/ and the faces that are surcharged
                where (volume(thisP) > fullVolume(thisP))
                    !% find slot properties
                    SlotVolume(thisP) = max(volume(thisP) - fullVolume(thisP), zeroR)
                    SlotArea(thisP)   = SlotVolume(thisP) / length(thisP)  
                    !% logicals
                    isfSlot(fUp(thisP)) = .true.
                    isfSlot(fDn(thisP)) = .true.
                    isSlot(thisP)       = .true.
                    !% smooth out preissmann number
                    ! PNumber = max(onehalfR * (fPNumber(fUP(thisP)) + fPNumber(fDn(thisP))), oneR)
                    !% find the preissmann celerity
                    PCelerity(thisP) = min(TargetPCelerity / PNumber(thisP), TargetPCelerity)
                    !% find the change in slot volume
                    dSlotVol(thisP)   = SlotVolume(thisP) - oldSlotVol(thisP)
                    !% find the change in slot area
                    dSlotArea(thisP)  = dSlotVol(thisP) / length(thisP)
                    !% find the change in slot depth
                    dSlotDepth(thisP) = (dSlotArea(thisP)  * (PCelerity(thisP) ** twoR)) / (grav * (fullArea(thisP)))
                end where

                !% --- reset isfSlot to find ventilated positions
                where (.not. isSlot(thisP))
                    isfSlot(fUp(thisP)) = .false.
                    isfSlot(fDn(thisP)) = .false.
                end where

                where (faceYN(fUP(thisP),fYN_isDiag_adjacent))
                    isfSlot(fUp(thisP)) = .false.
                end where

                where (faceYN(fDn(thisP),fYN_isDiag_adjacent))
                    isfSlot(fDn(thisP)) = .false.
                end where

                where (isfSlot(fUp(thisP)) .and. isfSlot(fDn(thisP)))
                    !% --- get a new decreased preissmann number for the next time step for elements having a slot
                    !    PNumber(thisP) = max((PNumber(thisP) ** twoR - PNumber(thisP) + oneR) / PNumber(thisP), oneR)
                    PNumber(thisP) = max(PNumber(thisP) - PNumber(thisP) * Dt, oneR)
                elsewhere
                    !% reset the preissmann number for every CC element not having a slot
                    PNumber(thisP) =  TargetPCelerity / (Alpha * sqrt(grav * ellMax(thisP)))
                end where

            case default
                !% --- should not reach this stage
                print*, 'In ', subroutine_name
                print *, 'CODE ERROR Slot Method type unknown for # ', SlotMethod
                print *, 'which has key ',trim(reverseKey(SlotMethod))
                call util_crashpoint(668732)

        end select

    end subroutine slot_CC_ETM
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine slot_JM_ETM (thisCol, Npack)
        !%------------------------------------------------------------------
        !% Description:
        !% Compute Preissmann slot for closed JM's in ETM methods
        !% Note this does NOT change the volume stored except in the case
        !% ponding or overflow
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisCol, Npack
            integer, pointer    :: thisP(:), SlotMethod
            real(8), pointer    :: fullarea(:), fullVolume(:), length(:)
            real(8), pointer    :: PNumber(:), PCelerity(:)
            real(8), pointer    :: volume(:) , SlotVolume(:), SlotDepth(:)
            real(8), pointer    :: SlotArea(:), SlotWidth(:), maxSlotDepth(:)
            real(8), pointer    :: VolumeExtra(:), VolumePonded(:), VolumeOverflow(:)
            real(8), pointer    :: TargetPCelerity, cfl, grav, Alpha, Dt
            logical, pointer    :: isSlot(:)
            integer :: ii

            character(64) :: subroutine_name = 'slot_JM_ETM'
        !%-----------------------------------------------------------------
        !% Preliminaries
            if (Npack < 1) return
        !%-----------------------------------------------------------------
        !% Aliases
        !% HACK -- many of these aliases are unused. Need to clean up 20220913 brh
            !% --- pointer packed element indexes
            thisP => elemP(1:Npack,thisCol)
            !% --- pointers to elemR columns
            fullArea   => elemR(:,er_FullArea)
            fullVolume => elemR(:,er_FullVolume)
            length     => elemR(:,er_Length)
            PNumber    => elemR(:,er_Preissmann_Number)
            PCelerity  => elemR(:,er_Preissmann_Celerity)
            SlotVolume => elemR(:,er_SlotVolume)
            SlotDepth  => elemR(:,er_SlotDepth)
            SlotArea   => elemR(:,er_SlotArea)
            SlotWidth  => elemR(:,er_SlotWidth)
            volume     => elemR(:,er_Volume)
            maxSlotDepth   => elemSR(:,esr_JunctionMain_SurchargeExtraDepth)
            VolumePonded   => elemSR(:,esr_JunctionMain_PondedVolume)
            VolumeExtra    => elemR(:,er_Temp01)
            VolumeOverflow => elemR(:,er_VolumeOverFlow)
            !% --- pointer to elemYN column
            isSlot     => elemYN(:,eYN_isPSsurcharged)
            !% --- pointer to necessary settings struct
            SlotMethod          => setting%Solver%PreissmannSlot%Method
            TargetPCelerity     => setting%Solver%PreissmannSlot%TargetCelerity
            Alpha               => setting%Solver%PreissmannSlot%Alpha
            cfl                 => setting%VariableDT%CFL_target
            grav                => setting%Constant%gravity
            Dt                  => setting%Time%Hydraulics%Dt
        !%-----------------------------------------------------------------
        ! !% Select the type of slot method
        ! select case (SlotMethod)
        !     !% for a static slot, the preissmann number will always be one.
        !     case (StaticSlot)

        !% --- HACK JM are always static slot. should test using JB average PNumber
        !%     for a dynamic slot

                !% initialize static slot
                PNumber(thisP)    = oneR  !% required for static slot
                SlotVolume(thisP) = zeroR
                SlotArea(thisP)   = zeroR
                SlotDepth(thisP)  = zeroR
                PCelerity(thisP)  = zeroR
                SlotWidth(thisP)  = zeroR
                VolumeExtra(thisP)= zeroR
                isSlot(thisP)     = .false.
                !% --- find the slot volume/ area/ and the faces that are surcharged

                !% --- initialize extra volume to zero so that non-surcharged JM will
                !%     not be affected. 
                !%     Note that a positive value is extra volume that would be in the 
                !%     slot but is lost to ponding/overflow. A negative value is the ability
                !%     take take in extra volume into the slot due to ponding (only). 
                !%     Remember the inflow from ponding into a non-surcharged JM is handled
                !%     in geo_ponding_inflow
                VolumeExtra(thisP) = zeroR

                !% --- compute baseline slot properties 
                !%     allows SlotDepth greater than maximum allowed
                where (volume(thisP) >= fullVolume(thisP))
                    isSlot(thisP)     = .true.
                    SlotVolume(thisP) = max(volume(thisP) - fullVolume(thisP), zeroR)
                    SlotArea(thisP)   = SlotVolume(thisP) / length(thisP)
                    PCelerity(thisP)  = min(TargetPCelerity / PNumber(thisP), TargetPCelerity)
                    SlotDepth(thisP)  = SlotArea(thisP)  * (PCelerity(thisP) ** twoR)  &
                                         / (grav * (fullArea(thisP)))
                    SlotWidth(thisP)  = SlotArea(thisP) / SlotDepth(thisP)
                    !% --- difference between SlotDepth and maximum is the "extra" volume
                    !%     If negative, this is volume that could be filled by ponding
                    !%     If positive, this is volume that is added to ponding or lost 
                    !%     in overflow
                    VolumeExtra(thisP) =  (SlotDepth(thisP) - maxSlotDepth(thisP)) &
                                         * SlotWidth(thisP) * length(thisP)
                end where

                !% --- if ponding is not allowed, there is no possibility of ponding inflow
                !%     so VolumeExtra must be positive or zero
                if (.not. setting%SWMMinput%AllowPonding) then
                    VolumeExtra(thisP) = max(VolumeExtra(thisP), zeroR)
                end if
    
                !% --- limit inflow volumes from ponding
                if (setting%SWMMinput%AllowPonding) then  
                    where (VolumeExtra(thisP) < zeroR)
                        !% --- inflow from ponding is the smaller volume of the extra available 
                        !%     in the slot and the amount in the ponded volume (negative values, 
                        !%     so use max). If there is no ponded volume, this returns 0 for extra
                        !%     volume
                        VolumeExtra(thisP) = max(VolumeExtra(thisP), -VolumePonded(thisp) )
                    endwhere
                end if

                where (VolumeExtra(thisP) .ne. zeroR) 
                    !% --- adjust slot for ponding inflow/outflow or overflow
                    !%     Note that non-surcharged will have VolumeExtra == zeroR
                    !%     so they are not affected. A negative volume extra is
                    !%     adding water to the slot.
                    SlotVolume(thisP) = SlotVolume(thisP) - VolumeExtra(thisP)
                    SlotArea(thisP)   = SlotVolume(thisP) / length(thisP)
                    PCelerity(thisP)  = min(TargetPCelerity / PNumber(thisP), TargetPCelerity)
                    SlotDepth(thisP)  = SlotArea(thisP)  * (PCelerity(thisP) ** twoR)  &
                                            / (grav * (fullArea(thisP)))
                    SlotWidth(thisP)  = SlotArea(thisP) / SlotDepth(thisP)
                    !% --- adjust the element volume for the extra volume leaving or entering
                    !%     (the latter from ponded only)
                    volume(thisP)     = volume(thisP) - VolumeExtra(thisP)
                end where

                if (setting%SWMMinput%AllowPonding) then 
                    VolumePonded(thisP)   = VolumePonded(thisP)   + VolumeExtra(thisP)
                else
                    !% --- note that VolumeExtra >= 0 is guaranteed if ponding isnot used
                    VolumeOverflow(thisP) = VolumeExtra(thisP)
                end if

            
            
        !     !% for dynamic slot, preissmann number is adjusted
        !     case (DynamicSlot)
                
        !         print *, 'In ', subroutine_name
        !         print *, 'Dynamic slot for junction mains is under development'
        !         call util_crashpoint(4098734)

        !     case default
        !         !% should not reach this stage
        !         print *, 'In ', subroutine_name
        !         print *, 'CODE ERROR Slot Method type unknown for # ', SlotMethod
        !         print *, 'which has key ',trim(reverseKey(SlotMethod))
        !         cal util_crashpoint(698723)

        ! end select

    end subroutine slot_JM_ETM
!%
!%==========================================================================
!%==========================================================================
!%
!%
!%==========================================================================
!%==========================================================================
!%
        !%------------------------------------------------------------------
        !% Description:
        !%
        !%------------------------------------------------------------------
        !% Declarations:
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------
        !% Aliases:
        !%------------------------------------------------------------------
    
    
        !%------------------------------------------------------------------
        !% Closing:

!%==========================================================================
!% END OF MODULE
!%==========================================================================
end module preissmann_slot