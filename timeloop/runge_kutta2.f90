module runge_kutta2

    use define_globals
    use define_indexes
    use define_keys
    use define_settings, only: setting
    use update
    use face
    use lowlevel_rk2
    use pack_mask_arrays, only: pack_small_and_zero_depth_elements
    use adjust
    use diagnostic_elements
    use utility, only: util_CLprint, util_syncwrite
    use utility_crash

    implicit none

    !%-----------------------------------------------------------------------------
    !% Description:
    !% Runge-Kutta time-marching method for both ETM and AC approaches with RK2
    !%-----------------------------------------------------------------------------

    private

    public :: rk2_toplevel_ETM
    public :: rk2_toplevel_AC
    public :: rk2_toplevel_ETMAC

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
    subroutine rk2_toplevel_ETM ()
        !%------------------------------------------------------------------
        !% Description:
        !% single RK2 step for explicit time advance of SVE
        !%------------------------------------------------------------------
        !% Declarations:
            integer :: istep, ii
            integer :: whichTM
            character(64) :: subroutine_name = 'rk2_toplevel_ETM'
        !%------------------------------------------------------------------
        !% Preliminaries
            !if (crashYN) return
            if (setting%Debug%File%runge_kutta2) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
            !% --- set the time-marching type for this routine
            whichTM = ETM
            !% --- reset the overflow counter for this time level
            elemR(:,er_VolumeOverFlow) = zeroR        
        !%-----------------------------------------------------------------
        !% Aliases
        !%-----------------------------------------------------------------

        ! print *, ' '
        ! call util_CLprint ('======= AAA  start of RK2 ==============================')

        !% --- compute the dynamic mannings N (DISABLED AS OF 20220817 brh)
        if (setting%Solver%ManningsN%useDynamicManningsN) then
            call rk2_dynamic_ManningsN (whichTM)
        else
            !% --- arguably not needed, but here to prevent bugs
            !%     if dynamic Mannings N is used in place of standard Mannings N
            elemR(:,er_ManningsN_Dynamic) = elemR(:,er_ManningsN)
        end if

        !% --- compute Force Main Manning's N for non-submerged FM elements
        if (setting%Solver%ForceMain%UseForceMainTF) call rk2_ForceMain_ManningsN ()

        !% --- RK2 solution step -- single time advance step for CC and JM
        istep=1
        call rk2_step_ETM (istep)

            ! call util_CLprint ('BBB after volume/momentum step 1---------------------------')
   
        !% --- RK2 solution step -- update all non-diagnostic aux variables
        !%     Note, these updates CANNOT depend on face values
        call update_auxiliary_variables (whichTM)

            ! call util_CLprint ('CCC  after update aux step 1-----------------------')

        !% --- set the flagged zero and small depth cells (allow depth to change)
        !%     This does not reset the zero/small depth packing
        call adjust_zero_and_small_depth_elem (whichTM, .false.)
        call util_crashstop(340927)

            ! call util_CLprint ('DDD  after adjust zero/small elem-----------------')
     
        !% --- RK2 solution step  -- all face interpolation
        sync all
        call face_interpolation(fp_all,whichTM)

            ! call util_CLprint ('EEE  after face interpolation step 1---------------')

        !% --- set the zero and small depth fluxes
        call adjust_zero_and_small_depth_face (whichTM, .true.)
        call util_crashstop(440223)

            ! call util_CLprint ('FFF  after zero/small face step 1-----------------')

        !% --- RK2 solution step  -- update diagnostic elements and faces
        call diagnostic_toplevel (.true.)
        call util_crashstop(402873)

            ! call util_CLprint ('GGG  after diagnostic step 1')

        !% --- RK2 solution step  -- make ad hoc adjustments
        call adjust_Vfilter (whichTM) ! brh20220211 this is useful in lateral flow induced oscillations
        call util_crashstop(13987)

            ! call util_CLprint ('HHH  after Vfilter step 1 -------------------------')
        
        !% -- the conservative fluxes from N to N_1 are the values just before the second RK2 step
        call rk2_store_conservative_fluxes (whichTM)

            ! call util_CLprint ('III  after consQ store step 1 ----------------------')

        !% --- reset the overflow counter (we only save conservation in the 2nd step)
        elemR(:,er_VolumeOverFlow) = zeroR

        !% --------------------------------------------------------------------------
        !% --- RK2 solution step -- RK2 second step for ETM 

        istep=2
        call rk2_step_ETM (istep)
        
            ! call util_CLprint ('JJJ  after volume rk2 step 2 -----------------------')

        !% --- RK2 solution step -- update non-diagnostic auxiliary variables
        !%     Note, these updates CANNOT depend on face values
        call update_auxiliary_variables(whichTM)  

            ! call util_CLprint ('KKK  after update aux step 2 --------------------------')

        !% --- set the flagged zero and small depth cells (allow depth to change)
        call adjust_zero_and_small_depth_elem (whichTM, .false.)
        call util_crashstop(12973)

            ! call util_CLprint ('LLL  after zero/small elem step 2 -------------------')

        !% --- RK2 solution step -- update all faces
        sync all
        call face_interpolation(fp_all,whichTM)

            ! call util_CLprint ('MMM  after face interp step 2 --------------------------')

        !% --- set the zero and small depth fluxes
        call adjust_zero_and_small_depth_face (whichTM, .false.)

            ! call util_CLprint ('NNN  after zero/small face step 2 ---------------------')

        !% --- RK2 solution step -- update diagnostic elements and faces
        call diagnostic_toplevel (.false.)
        call util_crashstop(662398)

            ! call util_CLprint ('OOO  after diagnostic step 2')
        
        !% --- RK2 solution step -- make ad hoc adjustments (V filter)
        call adjust_Vfilter (whichTM)
        call util_crashstop(449872)

            ! call util_CLprint ('PPP  after Vfilter step 2-----------------------------')

        !% --- ensures that the Vfilter hasn't affected the zero/small depth cells        
        call adjust_zero_and_small_depth_elem (whichTM, .true.)
        call util_crashstop(64987)

            ! call util_CLprint ('QQQ  after zero/small elem step 2 (2nd time)')

        !% --- accumulate the volume overflow
        elemR(:,er_VolumeOverFlowTotal) = elemR(:,er_VolumeOverFlowTotal) + elemR(:,er_VolumeOverFlow)

        
            !  call util_CLprint ('ZZZ  after accumulate overflow step 2')
            !  print *, '==================================================='
        

        !%-----------------------------------------------------------------
        !% closing
            if (setting%Debug%File%runge_kutta2)  &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine rk2_toplevel_ETM
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rk2_toplevel_AC ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !%
        !%-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'rk2_toplevel_AC'
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%runge_kutta2) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%runge_kutta2)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"


        print *, "need rk2_toplevel_AC to be written"
        stop 57839

    end subroutine rk2_toplevel_AC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rk2_toplevel_ETMAC ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !%
        !%-----------------------------------------------------------------------------
        integer :: istep, faceMaskCol, thisCol
        integer, pointer :: Npack
        character(64) :: subroutine_name = 'rk2_toplevel_ETMAC'
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%runge_kutta2) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        
        !% step 1 -- RK2 step 1 for ETM

        istep=1
        if (N_etm > 0) then
            call rk2_step_ETM (istep)
        end if

        !% step 2 -- RK2 step 1 for AC
        if (N_ac > 0) then
            !% step 2(a)
            call rk2_step_AC (istep)
            if (N_etm > 0) then
                !% step 2(b,c) create time n+1(1*) consistency for AC
                call rk2_extrapolate_to_fullstep_ETM()
            end if
        end if

        !% step 3 -- all aux variables for non-diagnostic
        call update_auxiliary_variables(ALLtm)

        !% step 4 -- all face interpolation
        sync all
        call face_interpolation(fp_all,ALLtm)
        

        !% step 5 -- update diagnostic elements and faces
        call diagnostic_toplevel (.true.)
        call util_crashstop(66234)

        !% step X -- make ad hoc adjustments
        call adjust_Vfilter (ALLtm)
        call util_crashstop(53434)

        !% step 6 -- RK2 step 2 for AC
        istep=2
        if (N_ac > 0) then
            !% step 6(a)
            call rk2_step_AC (istep)
            if (N_etm > 0) then
                !% step 6(b)
                call rk2_restore_to_midstep_ETM()
                !% step 6(c,d)
                call rk2_interpolate_to_halfstep_AC()
            end if
            !% step 6(e)
            call update_auxiliary_variables (AC)

            !% step 6(f,g) -- update faces for AC elements
            !% HACK -- need to sync some of the processors before this, but cannot sync all because some have N_ac <1
            call face_interpolation (fp_AC,AC)
            

        end if

        !% step 7 -- update diagnostic elements and faces
        call diagnostic_toplevel(.false.)
        call util_crashstop(12293)


        !% step 8 -- RK2 step 2 for ETM
        if (N_etm > 0) then
            !% step 8(a)
            call rk2_step_ETM (istep)

            if (N_ac > 0) then
                !% step 8(b)
                call rk2_restore_to_fullstep_AC ()
            end if

            !% step 8(c)
            call update_auxiliary_variables(ETM)

            !% step 8(d,e) -- update all faces
            !% HACK -- need to sync some of the processors before this, but cannot sync all because some have N_ac <1
            call face_interpolation (fp_all,ALLtm)
        end if

        !% step 9 -- update diagnostic elements and faces
        call diagnostic_toplevel (.false.)
        call util_crashstop(23422)

        !% step X -- make ad hoc adjustments
        call adjust_Vfilter (ALLtm)
        call util_crashstop(99287)

        if (setting%Debug%File%runge_kutta2)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine rk2_toplevel_ETMAC
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    subroutine rk2_step_ETM (istep)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single rk2 step for ETM
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: istep
        integer :: tmType
        !%-----------------------------------------------------------------------------
        !%
        !% perform the continuity step of the rk2 for ETM
        call rk2_continuity_step_ETM(istep)

            ! print *, this_image(),'    aaaa  after rk2 continuity step etm',this_image()
            !  call util_CLprint ('after rk2 continuity step etm')

        !% only adjust extremely small element volumes that have been introduced
        call adjust_limit_by_zerovalues (er_Volume, setting%ZeroValue%Volume/twentyR, col_elemP(ep_CCJM_H_ETM), .true.)

            ! print *, this_image(),'    bbbb  after rk2 call to adjust limit by zero',this_image()
            !  call util_CLprint ('after rk2 call to adjust limit by zero')

        !% perform the momentum step of the rk2 for ETM
        call rk2_momentum_step_ETM(istep)

            ! print *, this_image(),'    cccc  after rk2 call to rk2_momentum_step_ETM',this_image()
            !    call util_CLprint (' after rk2 call to rk2_momentum_step_ETM')

    end subroutine rk2_step_ETM
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rk2_step_AC (istep)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Single step of RK2 for AC method
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: istep
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        !% AC continuity
        call rk2_continuity_step_AC(istep)
        !% AC momentum
        call rk2_momentum_step_AC(istep)

    end subroutine rk2_step_AC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rk2_continuity_step_ETM (istep)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% perform the continuity step of the rk2 for ETM
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: istep
        integer, pointer :: thisPackCol, Npack
        logical :: isreset
        !%-----------------------------------------------------------------------------
        !%
        !if (crashYN) return
        !% Baseline continuity source:
        !% Compute net flowrates for channels, conduits 
        thisPackCol => col_elemP(ep_CC_H_ETM)
        Npack => npack_elemP(thisPackCol)
        if (Npack > 0) then
            call ll_continuity_netflowrate_CC (er_SourceContinuity, thisPackCol, Npack)
        end if

        ! print *, '----------aaaa '
        ! write(*,"(5f12.7)") elemR(iet(1),er_SourceContinuity)


        !% compute net flowrates for junction mains
        thisPackCol => col_elemP(ep_JM_ETM)
        Npack => npack_elemP(thisPackCol)
        if (Npack > 0) then
            call ll_continuity_netflowrate_JM (er_SourceContinuity, thisPackCol, Npack)
        end if

        ! print *, '------------bbbb  '
        ! write(*,"(5f12.7)") elemR(iet(1),er_SourceContinuity)

        !% Solve for volume in ETM step
        thisPackCol => col_elemP(ep_CCJM_H_ETM)
        Npack => npack_elemP(thisPackCol)
        if (Npack > 0) then
            call ll_continuity_volume_CCJM_ETM (er_Volume, thisPackCol, Npack, istep)
        end if

        ! print *, '------------cccc  '
        ! write(*,"(5f12.7)") elemR(iet(1),er_Volume)

        !% MOVED INTO SLOT_TOPLEVEL, called near top of geometry 
        ! !% compute Preissmann slot for conduits only if ETM solver is used
        ! if (setting%Solver%SolverSelect == ETM) then
        !     !% all the closed conduit elements
        !     thisPackCol => col_elemP(ep_CC_Closed_Elements)
        !     Npack => npack_elemP(thisPackCol)
        !     if (Npack > 0) then
        !         call ll_CC_slot_computation_ETM (thisPackCol, Npack)
        !     end if

        !     !% --- all the JM elements are closed elements 
        !     !%    (although SurchargeExtraDepth=0 is effectively open)
        !     thisPackCol => col_elemP(ep_JM)
        !     Npack => npack_elemP(thisPackCol)
        !     if (Npack > 0) then
        !         call ll_JM_slot_computation_ETM (thisPackCol, Npack)
        !     end if
        ! endif

        ! print *, '------------dddd  '
        ! write(*,"(5f12.7)") elemR(iet(1),er_Volume)



    end subroutine rk2_continuity_step_ETM
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rk2_continuity_step_AC (istep)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single RK2 continuity step for AC method
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: istep
        integer, pointer ::  thisMaskCol, thisPackCol, Npack
        logical :: isreset
        ! !%-----------------------------------------------------------------------------
        ! !%
        ! !if (crashYN) return
        ! !% Compute net flowrates for channels, conduits and special elements
        ! thisPackCol => col_elemP(ep_CC_AC)
        ! Npack => npack_elemP(thisPackCol)
        ! if (Npack > 0) then
        !     call ll_continuity_netflowrate_CC (er_SourceContinuity, thisPackCol, Npack)
        ! end if

        ! !% compute net flowrates for junction mains
        ! thisPackCol => col_elemP(ep_JM_AC)
        ! Npack => npack_elemP(thisPackCol)
        ! if (Npack > 0) then
        !     call ll_continuity_netflowrate_JM (er_SourceContinuity, thisPackCol, Npack)
        ! end if

        ! thisPackCol => col_elemP(ep_CCJM_H_AC_open)
        ! Npack => npack_elemP(thisPackCol)
        ! if (Npack > 0) then
        !     !% unique continuity source terms for AC open channel
        !     call ll_continuity_add_source_CCJM_AC_open (er_SourceContinuity, thisPackCol, Npack)
        !     !% unique continuity gamma terms for AC open channel
        !     call ll_continuity_add_gamma_CCJM_AC_open (er_GammaC, thisPackCol, Npack)
        !     !% solve for volume in AC open step
        !     call ll_continuity_volume_CCJM_AC_open (er_Volume, thisPackCol, Npack, istep)
        ! end if

        ! thisPackCol => col_elemP(ep_CCJM_H_ACsurcharged)
        ! Npack => npack_elemP(thisPackCol)
        ! if (Npack > 0) then
        !     !% unique continuity source terms for AC surcharged channel head
        !     call ll_continuity_add_source_CCJM_AC_surcharged (er_SourceContinuity, thisPackCol, Npack)
        !     !% solve for head in AC surcharged step
        !     call ll_continuity_head_CCJM_AC_surcharged (er_Head, thisPackCol, Npack, istep)
        ! end if

        ! !% adjust near-zero elements
        ! call adjust_limit_by_zerovalues (er_Volume, setting%ZeroValue%Volume, col_elemP(ep_CCJM_H_AC), .true.)

    end subroutine rk2_continuity_step_AC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rk2_momentum_step_ETM (istep)
        !%------------------------------------------------------------------
        !% Description:
        !% perform the momentum step of the rk2 for ETM (or ETM part of ETMAC)
        !%-------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: istep
            integer, pointer :: thisPackCol, Npack
            integer, pointer :: FMpackCol, nFMpack
            
            integer, parameter :: thisMethod = ETM
        !%-------------------------------------------------------------------
        !% Aliases
            thisPackCol => col_elemP(ep_CC_Q_ETM)
            Npack       => npack_elemP(thisPackCol)
        !%-------------------------------------------------------------------
        !%
        if (Npack > 0) then

            !% --- momentum K source terms for different methods for ETM
            call ll_momentum_Ksource_CC (er_Ksource, thisPackCol, Npack)
                !print *, '... Ksource :',elemR(1:3,er_Ksource)
                !  call util_CLprint (' after rk2 call ll_momentum_Ksource_CC')

            !% --- Common source for momentum on channels and conduits for ETM
            call ll_momentum_source_CC (er_SourceMomentum, thisPackCol, Npack)
                !print *, '... sM      :',elemR(1:3,er_SourceMomentum)
                !  call util_CLprint (' after rk2 call ll_momentum_source_CC')

            !% --- Common Gamma for momentum on channels and conduits for  ETM
            !%     Here for all channels and conduits, assuming CM roughness
            call ll_momentum_gammaCM_CC (er_GammaM, thisPackCol, Npack)
                !print *, '... gamma   :',elemR(1:3,er_GammaM)
                !  call util_CLprint (' after rk2 call ll_momentum_gamma_CC')

            !% --- handle force mains as Gamma terms
            !%     These overwrites the gamma from the CM roughness above
            if (setting%Solver%ForceMain%UseForceMainTF) then 
                !% --- surcharged Force main elements with Hazen-Williams roughness
                FMPackCol => col_elemP(ep_FM_HW_PSsurcharged)
                nFMpack   => npack_elemP(FMPackCol)
                if (nFMpack > 0) call ll_momentum_gammaFM_CC (er_GammaM, FMPackCol, nFMpack, HazenWilliams)

                !% --- surcharged Force Main elements with Darcy-Weisbach roughness
                FMPackCol => col_elemP(ep_FM_dw_PSsurcharged)
                nFMpack   => npack_elemP(FMPackCol)
                if (nFMpack > 0) call ll_momentum_gammaFM_CC (er_GammaM, FMPackCol, nFMpack, DarcyWeisbach)
            end if

            !% --- add minor loss term to gamma for all conduits
            call ll_minorloss_friction_gamma_CC (er_GammaM, thisPackCol, Npack)   

            !% --- Advance flowrate to n+1/2 for conduits and channels with ETM
            call ll_momentum_solve_CC (er_Velocity, thisPackCol, Npack, thisMethod, istep)
                !print *, '... vel     :',elemR(1:3,er_Velocity)
                !  call util_CLprint (' after rk2 call ll_momentum_solve_CC')

            !% --- velocity for ETM time march
            call ll_momentum_velocity_CC (er_Velocity, thisPackCol, Npack)
                !print *, '... vel     :',elemR(1:3,er_Velocity)
                !  call util_CLprint (' after rk2 call ll_momentum_velocity_CC')

            !% --- prevent backflow through flapgates
            call ll_enforce_flapgate_CC (er_Velocity, thisPackCol, Npack)
                !  call util_CLprint (' after rk2 call ll_enformce_flapgate_CC')

        end if

        !% --- update junction branches
        call ll_flowrate_and_velocity_JB(ETM,istep)

            ! print *, '   in rk2momentum 555'
            ! call util_CLprint (' after ll_flowrate_and_velocity_JB')

    end subroutine rk2_momentum_step_ETM
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rk2_momentum_step_AC (istep)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single hydrology step
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: istep
        integer, pointer    :: thisCol, Npack
        integer, parameter  :: thisMethod = AC
        !%-----------------------------------------------------------------------------
        !%
        !if (crashYN) return
        thisCol => col_elemP(ep_CC_Q_AC)
        Npack => npack_elemP(thisCol)
        if (Npack > 0) then
            !% momentum K source terms for different methods
            call ll_momentum_Ksource_CC (er_Ksource, thisCol, Npack)
            !% Source for momentum on channels and conduits
            call ll_momentum_source_CC (er_SourceMomentum, thisCol, Npack)
            !% Gamma for momentum on channels and conduits for Chezy-Manning roughness
            call ll_momentum_gammaCM_CC (er_GammaM, thisCol, Npack)
            !% additional momentum Gamma for AC time-march on channels and conduits
            call ll_momentum_add_gamma_CC_AC (er_GammaM, thisCol, Npack)
            !% additional momentum source terms for AC time-march on channels and conduits
            call ll_momentum_add_source_CC_AC (er_SourceMomentum, thisCol, Npack)
            !% AC elements advance flowrate to n+1(*) for conduits and channels
            call ll_momentum_solve_CC (er_Velocity, thisCol, Npack, thisMethod, istep)
        end if

        thisCol => col_elemP(ep_CC_Q_AC)
        Npack => npack_elemP(thisCol)
        if (Npack > 0) then
            !% velocity for AC time march
            call ll_momentum_velocity_CC (er_Velocity, thisCol,Npack)
        end if

    end subroutine rk2_momentum_step_AC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rk2_extrapolate_to_fullstep_ETM ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Finds ETM elements at time n+1/2 that are adjacent to fAC
        !% Makes temporary store of data for Q, H, V at n+1/2
        !% overwrites the Q, H, V location with an extrapolation to n+1.
        !%-----------------------------------------------------------------------------
        integer, pointer :: thisCol, Npack
        !%-----------------------------------------------------------------------------
        !if (crashYN) return

        !% brh20211212 hard stop until fixe is made
        print *, 'CODE ERROR problems with CCJB_eETM_i_fAC mask need to be fixed'
        stop 68795
        
        thisCol => col_elemP(ep_CCJB_eETM_i_fAC)
        Npack => npack_elemP(thisCol)


        if (Npack > 0) then
            !% temporary storage of n+1/2 data
            call ll_store_in_temporary (thisCol, Npack)

            !% extrapolation
            call ll_extrapolate_values (thisCol, Npack)

            !% update aux for extrapolated variables
    !        call update_auxiliary_variables_byPack (thisCol, Npack)
        end if

    end subroutine rk2_extrapolate_to_fullstep_ETM
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rk2_restore_to_midstep_ETM ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% reverses the effect of ETM_extrapolate_to_fullstep
        !% Finds ETM elements at time n+1/2 that are adjacent to fAC
        !% restsores data of Q, H, V at n+1/2
        !%-----------------------------------------------------------------------------
        integer, pointer ::  thisCol, Npack
        !%-----------------------------------------------------------------------------
        !%
        !if (crashYN) return
        thisCol => col_elemP(ep_CCJB_eETM_i_fAC)
        Npack   => npack_elemP(thisCol)

        if (Npack > 0) then
            !% temporary storage of n+1/2 data
            call ll_restore_from_temporary (thisCol, Npack)

            !% update aux for restored variables
            ! call update_auxiliary_variables_byPack (thisPackCol, Npack)
        end if

    end subroutine rk2_restore_to_midstep_ETM
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rk2_interpolate_to_halfstep_AC ()
        !%-----------------------------------------------------------------------------
    !     !% Description:
    !     !% Finds AC and Diag elements at time n+1/2 that are adjacent to fETM
    !     !% Makes temporary store of data for Q, H, V at n+1(*)
    !     !% overwrites the Q, H, V location with an interpolation to n+1/1.
    !     !%-----------------------------------------------------------------------------
    !     integer, pointer :: thisCol, Npack
    !     !%-----------------------------------------------------------------------------
    !     !if (crashYN) return
    !     !%
    !     thisCol => col_elemP( ep_CCJB_eAC_i_fETM)
    !     Npack => npack_elemP(thisCol)

    !     if (Npack > 0) then
    !         !% temporary storage of n+1 data
    !         call ll_store_in_temporary (thisCol, Npack)

    !         !% interpolation to half step
    !         call ll_interpolate_values (thisCol, Npack)

    !         !% update aux for interpolated variables
    !   !     call update_auxiliary_variables_byPack (thisPackCol, Npack)
    !     end if

    end subroutine rk2_interpolate_to_halfstep_AC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rk2_restore_to_fullstep_AC()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% reverses the effect of AC_interpolate_to_halfstep
        !% Finds AC elements at time n+1/2 that are adjacent to fETM
        !% restsores data of Q, H, V at n+1/2
        !%-----------------------------------------------------------------------------
        integer, pointer :: thisCol, Npack
    !     !%-----------------------------------------------------------------------------
    !     !if (crashYN) return
    !     thisCol = col_elemP(ep_CCJB_eAC_i_fETM)
    !     Npack => npack_elemP(thisCol)

    !     if (Npack > 0) then
    !         !% temporary storage of n+1 data
    !         call ll_restore_from_temporary (thisCol, Npack)

    !         !% update aux for restored data
    !  !       call update_auxiliary_variables_byPack (thisPackCol, Npack)
    !     end if

    end subroutine rk2_restore_to_fullstep_AC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rk2_store_conservative_fluxes (whichTM)
        !%------------------------------------------------------------------
        !% Description:
        !% store the intermediate face flow rates in the Rk2 which are
        !% the conservative flowrate over the time step
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: whichTM
            integer, pointer    :: NpackE, NpackF, NpackJ
            integer, pointer    :: thisColE,   thisColF, thisColJ, isbranch(:)
            integer, pointer    :: thisP(:), thisF(:), thisJ(:), fup(:), fdn(:)
            real(8), pointer    :: fQcons(:), fQ(:)
            integer :: ii
        !%------------------------------------------------------------------
        !% Preliminaries:
            select case (whichTM)
            case (ETM)
                thisColE  => col_elemP(ep_CC_ETM)
                thisColJ  => col_elemP(ep_JM_ETM)
                thisColF  => col_faceP(fp_Diag)
            case default
                print *, 'CODE ERROR: time march type unknown for # ', whichTM
                print *, 'which has key ',trim(reverseKey(whichTM))
                stop 398705
            end select
            
        !%------------------------------------------------------------------
        !% Aliases:
            NpackE => npack_elemP(thisColE)
            NpackF => npack_faceP(thisColF)
            NpackJ => npack_elemP(thisColJ)
            fup    => elemI(:,ei_Mface_uL)
            fdn    => elemI(:,ei_Mface_dL)
            fQcons => faceR(:,fr_Flowrate_Conservative)
            fQ     => faceR(:,fr_Flowrate)
            isbranch => elemSI(:,esi_JunctionBranch_Exists)
        !%------------------------------------------------------------------

        !% handle the flux faces of the time-marching elements
        if (NpackE > 0) then
            thisP  => elemP(1:NpackE, thisColE)
            fQcons(fup(thisP)) = fQ(fup(thisP))
            fQcons(fdn(thisP)) = fQ(fdn(thisP))
        end if        

        !% handle flux faces of junctions
        if (NpackJ > 0) then
            thisJ  => elemP(1:NpackJ,thisColJ)
            do ii=1,max_branch_per_node,2
                fQcons(fup(thisJ+ii  )) = fQ(fup(thisJ+ii  )) * real(isbranch(thisJ+ii  ),8)
                fQcons(fdn(thisJ+ii+1)) = fQ(fdn(thisJ+ii+1)) * real(isbranch(thisJ+ii+1),8)
            end do
        end if
    
        !% handle the flux faces of the diagnostic elements
        if (NpackF > 0) then
            thisF  => faceP(1:NpackF, thisColF)
            fQcons(thisF) = fQ(thisF)
        end if

        !%------------------------------------------------------------------
        !% Closing
        !%
    end subroutine rk2_store_conservative_fluxes
!%   
!%==========================================================================
    !%==========================================================================
!%
    subroutine rk2_dynamic_ManningsN (whichTM) 
        !%------------------------------------------------------------------
        !% Description:
        !% Updates the baseline Manning's n with a dynamic value
        !% adjustment
        !%
        !% DISABLED AS OF 20220817 -- HAS PROBLEMS WITH SMALL DEPTHS
        !%
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: whichTM
            integer, pointer    :: npack, thisColCC, thisColJM
            integer, pointer    :: BranchExists(:)
            integer, pointer    :: thisP(:), fUp(:), fDn(:), tM
            integer             :: ii, kk, tB(1), dpnorm_col
            real(8), pointer    :: dynamic_mn(:), mn(:), dp_norm(:)
            real(8), pointer    :: eHead(:), fHead_d(:), fHead_u(:)
            real(8), pointer    :: zBottom(:), volume(:), length(:)
            real(8), pointer    :: alpha, dt
            character(64) :: subroutine_name ='rk2_dynamic_ManningsN'
        !%------------------------------------------------------------------  
            if (.not. setting%Solver%ManningsN%useDynamicManningsN) return
        !%------------------------------------------------------------------  
        !% Aliases
            dpnorm_col = er_Temp01 !% not an alias!
            dt           => setting%Time%Hydraulics%Dt
            alpha        => setting%Solver%ManningsN%alpha
            dynamic_mn   => elemR(:,er_ManningsN_Dynamic)
            mn           => elemR(:,er_ManningsN)
            dp_norm      => elemR(:,dpnorm_col)
            eHead        => elemR(:,er_Head)
            length       => elemR(:,er_Length)
            zBottom      => elemR(:,er_Zbottom)
            volume       => elemR(:,er_Volume)
            fHead_d      => faceR(:,fr_Head_d)
            fHead_u      => faceR(:,fr_Head_u)
            fUp          => elemI(:,ei_Mface_uL)
            fDn          => elemI(:,ei_Mface_dL)
            BranchExists => elemSI(:,esi_JunctionBranch_Exists)
        !%------------------------------------------------------------------
        !% Preliminaries:
            select case (whichTM)
            case (ETM)
                thisColCC  => col_elemP(ep_CC_ETM)
                thisColJM  => col_elemP(ep_JM_ETM)
            case default
                print *, 'CODE ERROR: time march type not handled for # ', whichTM
                print *, 'which has key ',trim(reverseKey(whichTM))
                stop 398705
            end select    
        !%------------------------------------------------------------------  

        print *, 'DYNAMIC MANNINGS CANNOT BE USED'
        print *, 'change setting%Solver%ManningsN%useDynamicManningsN = .false.'    
        stop 20987344

        !% --- compute ManningsN for CC elements    
        npack      => npack_elemP(thisColCC)
        if (npack .ge.1)  then

            thisP      => elemP(1:npack,thisColCC)   
            !% --- the normalized pressure gradient scale
            ! dp_norm(thisP) = (  abs(fHead_d(fUp(thisP)) - eHead(thisP))   &
            !                   + abs(fHead_u(fDn(thisP)) - eHead(thisP)) ) &
            !                  / abs(eHead(thisP) - zBottom(thisP))

            dp_norm(thisp) = (abs(elemR(thisP,er_Velocity_N0) - elemR(thisP,er_Velocity_N1)) / setting%Time%Hydraulics%Dt)  &
            / ( &
            (elemR(thisP,er_ManningsN)**2) * (elemR(thisP,er_Velocity_N0)**2) * setting%Constant%gravity &
            / (elemR(thisP,er_HydRadius)**fourthirdsR) &
            )
        
            !print *, 'in ',trim(subroutine_name)
            !print *, dp_norm(33)
            ! print *, fHead_d(fUp(139)), eHead(139)
            ! print *, fHead_u(fDn(139)), eHead(139)
            ! print *, eHead(139),zBottom(139)
           
            call ll_get_dynamic_ManningsN (thisP, dpnorm_col) 

            !dynamic_mn(thisP) =  mn(thisP) &
            !    +  alpha *  (dt / ((length(thisP))**(onethirdR))) * (exp(dp_norm(thisP)) - oneR )        
        end if

        !% --- compute ManningsN for JB elements
        npack => npack_elemP(thisColJM)
 
        ! if (Npack > 0) then
        !     do ii=1,Npack
        !         tM => elemP(ii,thisColJM)  !% JM junction main ID
        !         !% --- handle the upstream branches
        !         do kk=1,max_branch_per_node,2
        !             tB(1) = tM + kk  !% JB branch ID
        !             if (BranchExists(tB(1))==1) then
        !                 !% --- normalized head difference is with upstream face
        !                 dp_norm(tB) = (abs(fHead_d(fUp(tB)) - eHead(tB))) &
        !                      / abs(eHead(tB) - zBottom(tB))
        !                 !% --- add the dynamic ManningsN    
        !                 call ll_get_dynamic_ManningsN (tB, dpnorm_col)       
        !             else
        !                 !% skip if not a valid branch
        !             end if
        !         end do
        !         do kk=2,max_branch_per_node,2
        !             tB(1) = tM + kk  !% JB branch ID
        !             if (BranchExists(tB(1))==1) then
        !                 !% --- normalized head difference is with downstream face
        !                 dp_norm(tB) = (abs(fHead_u(fDn(tB)) - eHead(tB))) &
        !                      / abs(eHead(tB) - zBottom(tB))
        !                 !% --- add the dynamic ManningsN  
        !                 call ll_get_dynamic_ManningsN (tB, dpnorm_col)    
        !             else
        !                 !% skip if not a valid branch
        !             end if
        !         end do
        !     end do
        ! end if

        ! print *, 'in ',trim(subroutine_name)
        ! print *, ' '
        ! print *, 'dynamic mn'
        ! print *, dynamic_mn(thisP)
        ! print *, ' '
        ! print *, 'dp norm'
        ! print *, dp_norm(thisP)

        ! print *, ' '
        ! print *, abs(fHead_d(fUp(thisP)) - eHead(thisP))
        ! print *, ' '
        ! print *, abs(fHead_u(fDn(thisP)) - eHead(thisP)) 
        ! print *, ' '
        ! print *, (eHead(thisP) - zBottom(thisP))
        ! print *, ' '
        ! print *, ' '

        ! print *, ' '
        ! print *, fHead_d(fUp(thisP))
        ! print *, ' '
        ! print *, eHead(thisP)
        ! print *, ' '
        ! print *, fHead_u(fDn(thisP))

        ! print *, ' '
        ! print *, 'mn '
        ! print *, mn(thisP)
        ! print *, ' '
        ! print *, 'volume '
        ! print *, volume(thisP)
        ! print *, ' '
        ! print *, 'volume**1/9'
        ! print *,  volume(thisP)**(oneninthR)
        ! print *, ' '
        ! print *, '1-e'
        ! print *, (oneR - exp(dp_norm(thisP)) )


        dp_norm(:) = nullvalueR  

       ! stop 398745

    end subroutine rk2_dynamic_ManningsN
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rk2_ForceMain_ManningsN ()
        !%------------------------------------------------------------------
        !% Description:
        !% Calls lower-level subroutines to set the equivalent Mannings n
        !% for face mains that have a free surface.
        !%------------------------------------------------------------------
        !% Declarations:
            integer, pointer :: thisPackCol, Npack, thisP(:)
        !%------------------------------------------------------------------
        !% Preliminaries:
            if (.not. setting%Solver%ForceMain%UseForceMainTF) return
        !%------------------------------------------------------------------
        !% Aliases:
        !%------------------------------------------------------------------
        !%
        !% --- for Hazen-Williams
        if (.not. setting%Solver%ForceMain%HazenWilliams_equivalent_isSetTF) then
            !% --- the HW equivalent only needs to be called once and the equivalent
            !%     Manning's n is stored for re-use
            thisPackCol => col_elemP(ep_FM_HW_all)
            Npack       => npack_elemP(thisPackCol)
            if (Npack > 0) then
                call ll_ForceMain_equivalent_manningsN (thisPackCol,Npack,HazenWilliams)
            endif
            !% --- only check this once. The ep_FM_HW_all is a static pack, so
            !%     if Npack = 0 on the first time through it will always be zero
            setting%Solver%ForceMain%HazenWilliams_equivalent_isSetTF = .true.
        else 
            !% --- skip if already set
        end if

        !% --- for Darcy-Weisbach
        !%     Equivalent Mannings n is function of friction factor, so needs 
        !%     to be re-computed at each time step for the non-surcharged FM elements
        thisPackCol => col_elemP(ep_FM_dw_PSnonSurcharged)
        Npack       => npack_elemP(thisPackCol)
        if (Npack > 0) then
            call ll_ForceMain_dw_friction (thisPackCol, Npack)
            call ll_ForceMain_equivalent_manningsN (thisPackCol,Npack,DarcyWeisbach)
        endif

        !%------------------------------------------------------------------
        !% Closing
        !%
    end subroutine rk2_ForceMain_ManningsN
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
        !% Closing
        !%
!%==========================================================================
!% END OF MODULE
!%==========================================================================
end module runge_kutta2

