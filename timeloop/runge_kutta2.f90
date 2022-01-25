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
    use utility, only: util_CLprint

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
            integer :: istep
            character(64) :: subroutine_name = 'rk2_toplevel_ETM'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (icrash) return
            if (setting%Debug%File%runge_kutta2) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------
        !% Aliases
        !%-----------------------------------------------------------------
        !% --- RK2 solution step -- single time advance step for CC and JM
        istep=1
        call rk2_step_ETM (istep)

        !% --- RK2 solution step -- update all non-diagnostic aux variables
        call update_auxiliary_variables (ETM)

        !% --- set the flagged zero and small depth cells (allow depth to change)
        !%     This does not reset the zero/small depth packing
        call adjust_zero_and_small_depth_elem (ETM, .false.)
     
        !% --- RK2 solution step  -- all face interpolation
        call face_interpolation(fp_all,ETM)

        !% --- set the zero and small depth fluxes
        call adjust_zero_and_small_depth_face (ETM, .true.)

        !% --- RK2 solution step  -- update diagnostic elements and faces
        call diagnostic_toplevel()

        !% --- RK2 solution step  -- make ad hoc adjustments
        call adjust_Vfilter (ETM) ! brh20220211 this is useful in lateral flow induced oscillations
        
        !% -- the conservative fluxes from N to N_1 are the values just before the second RK2 step
        call rk2_store_conservative_fluxes (ETM)

        !% --------------------------------------------------------------------------
        !% --- RK2 solution step -- RK2 second step for ETM 
        istep=2
        call rk2_step_ETM (istep)

        !% --- RK2 solution step -- update non-diagnostic auxiliary variables
        call update_auxiliary_variables(ETM)  

        !% --- set the flagged zero and small depth cells (allow depth to change)
        call adjust_zero_and_small_depth_elem (ETM, .false.)

        !% --- RK2 solution step -- update all faces
        call face_interpolation(fp_all,ETM)

        !% --- set the zero and small depth fluxes
        call adjust_zero_and_small_depth_face (ETM, .false.)
        
        !% --- RK2 solution step -- update diagnostic elements and faces
        call diagnostic_toplevel()
        
        !% --- RK2 solution step -- make ad hoc adjustments (V filter)
        call adjust_Vfilter (ETM)

        !% --- ensures that the Vfilter hasn't affected the zero/small depth cells        
        call adjust_zero_and_small_depth_elem (ETM, .true.)

        !% --- accumulate the volume overflow
        elemR(:,er_VolumeOverFlowTotal) = elemR(:,er_VolumeOverFlowTotal) + elemR(:,er_VolumeOverFlow)

   
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
        if (icrash) return
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
        if (icrash) return
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
        call face_interpolation(fp_all,ALLtm)
        

        !% step 5 -- update diagnostic elements and faces
        call diagnostic_toplevel ()

        !% step X -- make ad hoc adjustments
        call adjust_Vfilter (ALLtm)

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
            call face_interpolation (fp_AC,AC)
            

        end if

        !% step 7 -- update diagnostic elements and faces
        call diagnostic_toplevel()


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
            call face_interpolation (fp_all,ALLtm)
        end if

        !% step 9 -- update diagnostic elements and faces
        call diagnostic_toplevel

        !% step X -- make ad hoc adjustments
        call adjust_Vfilter (ALLtm)

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
        if (icrash) return
        !%
        !% perform the continuity step of the rk2 for ETM
        call rk2_continuity_step_ETM(istep)

        !% only adjust extremely small element volumes that have been introduced
        call adjust_limit_by_zerovalues (er_Volume, setting%ZeroValue%Volume/twentyR, col_elemP(ep_CCJM_H_ETM))

        !% perform the momentum step of the rk2 for ETM
        call rk2_momentum_step_ETM(istep)

       

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
        if (icrash) return
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
        if (icrash) return
        !% Baseline continuity source:
        !% Compute net flowrates for channels, conduits and special elements
        thisPackCol => col_elemP(ep_CC_H_ETM)
        Npack => npack_elemP(thisPackCol)
        if (Npack > 0) then
            call ll_continuity_netflowrate_CC (er_SourceContinuity, thisPackCol, Npack)
        end if

        !print *, '----------aaaa '
        !write(*,"(5f12.7)") elemR(ietmp,er_SourceContinuity)
        !write(*,"(5f12.7)") elemR(1:3,er_SourceContinuity)

        !% compute net flowrates for junction mains
        thisPackCol => col_elemP(ep_JM_ETM)
        Npack => npack_elemP(thisPackCol)
        if (Npack > 0) then
            call ll_continuity_netflowrate_JM (er_SourceContinuity, thisPackCol, Npack)
        end if

        !print *, '------------bbbb  '
        !write(*,"(5f12.7)") elemR(ietmp,er_SourceContinuity)
        !write(*,"(5f12.7)") elemR(1:3,er_SourceContinuity)

        !% Solve for volume in ETM step
        thisPackCol => col_elemP(ep_CCJM_H_ETM)
        Npack => npack_elemP(thisPackCol)
        if (Npack > 0) then
            call ll_continuity_volume_CCJM_ETM (er_Volume, thisPackCol, Npack, istep)
        end if

        !print *, '------------cccc  '
        !write(*,"(5f12.7)") elemR(ietmp,er_Volume)
        !write(*,"(5f12.7)") elemR(1:3,er_Volume)

        !% compute slot for conduits only if ETM solver is used
        if (setting%Solver%SolverSelect == ETM) then
            !% all the closed conduit elements
            thisPackCol => col_elemP(ep_Closed_Elements)
            Npack => npack_elemP(thisPackCol)
            if (Npack > 0) then
                call ll_slot_computation_ETM (thisPackCol, Npack)
            end if
        endif

        !print *, '------------dddd  '
        !write(*,"(5f12.7)") elemR(ietmp,er_Volume)
        !write(*,"(5f12.7)") elemR(1:3,er_Volume)

        !

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
        !%-----------------------------------------------------------------------------
        !%
        if (icrash) return
        !% Compute net flowrates for channels, conduits and special elements
        thisPackCol => col_elemP(ep_CC_AC)
        Npack => npack_elemP(thisPackCol)
        if (Npack > 0) then
            call ll_continuity_netflowrate_CC (er_SourceContinuity, thisPackCol, Npack)
        end if

        !% compute net flowrates for junction mains
        thisPackCol => col_elemP(ep_JM_AC)
        Npack => npack_elemP(thisPackCol)
        if (Npack > 0) then
            call ll_continuity_netflowrate_JM (er_SourceContinuity, thisPackCol, Npack)
        end if

        thisPackCol => col_elemP(ep_CCJM_H_AC_open)
        Npack => npack_elemP(thisPackCol)
        if (Npack > 0) then
            !% unique continuity source terms for AC open channel
            call ll_continuity_add_source_CCJM_AC_open (er_SourceContinuity, thisPackCol, Npack)
            !% unique continuity gamma terms for AC open channel
            call ll_continuity_add_gamma_CCJM_AC_open (er_GammaC, thisPackCol, Npack)
            !% solve for volume in AC open step
            call ll_continuity_volume_CCJM_AC_open (er_Volume, thisPackCol, Npack, istep)
        end if

        thisPackCol => col_elemP(ep_CCJM_H_AC_surcharged)
        Npack => npack_elemP(thisPackCol)
        if (Npack > 0) then
            !% unique continuity source terms for AC surcharged channel head
            call ll_continuity_add_source_CCJM_AC_surcharged (er_SourceContinuity, thisPackCol, Npack)
            !% solve for head in AC surcharged step
            call ll_continuity_head_CCJM_AC_surcharged (er_Head, thisPackCol, Npack, istep)
        end if

        !% adjust near-zero elements
        call adjust_limit_by_zerovalues (er_Volume, setting%ZeroValue%Volume, col_elemP(ep_CCJM_H_AC))

    end subroutine rk2_continuity_step_AC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rk2_momentum_step_ETM (istep)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% perform the momentum step of the rk2 for ETM (or ETM part of ETMAC)
        integer, intent(in) :: istep
        integer, pointer :: thisPackCol, Npack
        integer, parameter :: thisMethod = ETM
        !%-----------------------------------------------------------------------------
        thisPackCol => col_elemP(ep_CC_Q_ETM)
        Npack       => npack_elemP(thisPackCol)
        !%-----------------------------------------------------------------------------
        !%
        if (icrash) return
        if (Npack > 0) then
            !print *, '... vel    :',elemR(1:2,er_Velocity)

            !% momentum K source terms for different methods for ETM
            call ll_momentum_Ksource_CC (er_Ksource, thisPackCol, Npack)
            !print *, '... Ksource :',elemR(1:3,er_Ksource)

            !% Common source for momentum on channels and conduits for ETM
            call ll_momentum_source_CC (er_SourceMomentum, thisPackCol, Npack)
            !print *, '... sM      :',elemR(1:3,er_SourceMomentum)

            !% Common Gamma for momentum on channels and conduits for  ETM
            call ll_momentum_gamma_CC (er_GammaM, thisPackCol, Npack)
            !print *, '... gamma   :',elemR(1:3,er_GammaM)

            !% Advance flowrate to n+1/2 for conduits and channels with ETM
            call ll_momentum_solve_CC (er_Velocity, thisPackCol, Npack, thisMethod, istep)
            !print *, '... vel     :',elemR(1:3,er_Velocity)

            !% velocity for ETM time march
            call ll_momentum_velocity_CC (er_Velocity, thisPackCol, Npack)
            !print *, '... vel     :',elemR(1:3,er_Velocity)

        end if

        call ll_junction_branch_flowrate_and_velocity(ETM,istep)

       ! write(*,"(5f12.7)") elemR(1,er_Velocity)
       ! stop 835783

        ! !% junction branch momentum source
        ! if (setting%Junction%isDynamicYN) then
        !     call ll_momentum_source_JB (ETM, istep)
        ! end if

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
        if (icrash) return
        thisCol => col_elemP(ep_CC_Q_AC)
        Npack => npack_elemP(thisCol)
        if (Npack > 0) then
            !% momentum K source terms for different methods
            call ll_momentum_Ksource_CC (er_Ksource, thisCol, Npack)
            !% Source for momentum on channels and conduits
            call ll_momentum_source_CC (er_SourceMomentum, thisCol, Npack)
            !% Gamma for momentum on channels and conduits
            call ll_momentum_gamma_CC (er_GammaM, thisCol, Npack)
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
        if (icrash) return

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
        if (icrash) return
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
        !% Description:
        !% Finds AC and Diag elements at time n+1/2 that are adjacent to fETM
        !% Makes temporary store of data for Q, H, V at n+1(*)
        !% overwrites the Q, H, V location with an interpolation to n+1/1.
        !%-----------------------------------------------------------------------------
        integer, pointer :: thisCol, Npack
        !%-----------------------------------------------------------------------------
        if (icrash) return
        !%
        thisCol => col_elemP( ep_CCJB_eAC_i_fETM)
        Npack => npack_elemP(thisCol)

        if (Npack > 0) then
            !% temporary storage of n+1 data
            call ll_store_in_temporary (thisCol, Npack)

            !% interpolation to half step
            call ll_interpolate_values (thisCol, Npack)

            !% update aux for interpolated variables
      !     call update_auxiliary_variables_byPack (thisPackCol, Npack)
        end if

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
        !%-----------------------------------------------------------------------------
        if (icrash) return
        thisCol = col_elemP(ep_CCJB_eAC_i_fETM)
        Npack => npack_elemP(thisCol)

        if (Npack > 0) then
            !% temporary storage of n+1 data
            call ll_restore_from_temporary (thisCol, Npack)

            !% update aux for restored data
     !       call update_auxiliary_variables_byPack (thisPackCol, Npack)
        end if

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
                print *, 'CODE ERROR: incomplete code -- not developed yet for other time march'
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

