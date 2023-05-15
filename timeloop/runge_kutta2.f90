module runge_kutta2

    use define_globals
    use define_indexes
    use define_keys
    use define_settings, only: setting
    use update
    use face
    use geometry
    use geometry_lowlevel, only: llgeo_head_from_depth_pure
    use forcemain, only: forcemain_ManningsN
    use junction_elements
    use lowlevel_rk2
    use culvert_elements, only: culvert_toplevel
    use pack_mask_arrays
    use preissmann_slot
    use adjust
    use diagnostic_elements
    use utility, only: util_syncwrite
    use utility_crash
    use utility_unit_testing, only: util_utest_CLprint, util_utest_checkIsNan

    implicit none

    !%-----------------------------------------------------------------------------
    !% Description:
    !% Runge-Kutta time-marching method for both ETM and AC approaches with RK2
    !%-----------------------------------------------------------------------------

    private

    public :: rk2_toplevel_ETM
    !public :: rk2_toplevel_ETM_2
    !public :: rk2_toplevel_ETM_3
    !public :: rk2_toplevel_ETM_4
    !public :: rk2_toplevel_ETM_5
    !public :: rk2_toplevel_ETM_6
    !public :: rk2_toplevel_ETM_7
    !public :: rk2_toplevel_AC
    !public :: rk2_toplevel_ETMAC

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine rk2_toplevel_ETM ()
        !%------------------------------------------------------------------
        !% Description:
        !% single RK2 step for explicit time advance of SVE
        !%------------------------------------------------------------------
        !% Declarations:
            integer :: istep, ii, kk, mm
            integer, pointer :: Npack, thisP(:), fup(:), fdn(:), thisJB(:)
            integer, pointer :: JMar(:), fAdj(:)
            
            integer :: fr_hAdj, fr_aAdj

            real(8), pointer :: deltaH(:), energyQ(:), grav, dt
            real(8)          :: bsign
            real(8)          :: volume1, volume2, inflowVolume, outflowVolume

            logical :: isUpstreamFace, isUpstreamBranch

            integer :: epCCcol

            integer :: printJB = 75

            character(64) :: subroutine_name = 'rk2_toplevel_ETM'
            real(8) :: ConsSum, tval(1)
        !%------------------------------------------------------------------
        !% Preliminaries
        !% --- reset the overflow counter for this time level
            elemR(:,er_VolumeOverFlow) = zeroR     
            elemR(:,er_VolumeArtificialInflow) = zeroR   
        !%-----------------------------------------------------------------
            grav => setting%Constant%gravity
            dt   => setting%Time%Hydraulics%Dt
        !%-----------------------------------------------------------------

        !% --- total volume conservation setup
        volume1 = zeroR
        do ii=1,N_elem(1)
            if ((elemI(ii,ei_elementType) == JM) .or. (elemI(ii,ei_elementType) == CC)) then
                volume1 = volume1 + elemR(ii,er_Volume)
            end if
        end do
            
        !print *, setting%Time%Step, setting%Time%Hydraulics%Dt
      
          !  call util_utest_CLprint ('======= AAA  start of RK2 ==============================')    
            
        !%==================================    
        !% --- Diagnostic and junction adjustments before RK2
        istep = zeroI

        !% saz 20230507 commented out
        ! !% --- Diagnostic elements and faces
        ! if (N_diag > 0) then 
        !     !% --- update flowrates for aa diagnostic elements
        !     call diagnostic_by_type (ep_Diag, istep)  
        !         ! ! call util_utest_CLprint ('------- BBB  after diagnostic')

        !     !% --- push the diagnostic flowrate data to faces -- true is upstream, false is downstream
        !     call face_push_elemdata_to_face (ep_Diag, fr_Flowrate, er_Flowrate, elemR, .true.)
        !     call face_push_elemdata_to_face (ep_Diag, fr_Flowrate, er_Flowrate, elemR, .false.)
        !         ! ! call util_utest_CLprint ('------- CCC  after face_push_elemdata_to_face')
        ! end if

        !% AT THIS POINT: DIAGNOSTIC VALUES ENFORCED ON ELEMENTS AND FACES, BUT CONNECTED JB ARE INCONSISTENT

        !% --- Preliminary values for JM/JB elements
        !%     Note, this must be called even if no JM/JB on this image because 
        !%     the faces require synchronizing.
        call junction_preliminaries ()

            ! ! call util_utest_CLprint ('------- BBB after junction preliminaries')


        !%==================================  
        !% --- RK2 SOLUTION
        do istep = 1,2
            ! print *, ' '
            ! print *, istep, '= iSTEP /////////////////////////////////////// ISTEP = ',istep
            ! print *, ' '
        
            !% --- Half-timestep advance on CC for U and UVolume
            call rk2_step_ETM_CC (istep)  
                ! ! call util_utest_CLprint ('------- III  after rk2_step')

            !% --- Update all CC aux variables
            !%     Note, these updates CANNOT depend on face values
            !%     Through geometry, this sets Preissmann Slot variables
            call update_auxiliary_variables_CC (                  &
                ep_CC, ep_CC_Open_Elements, ep_CC_Closed_Elements, &
                .true., .false., dummyIdx)
                ! ! call util_utest_CLprint ('------- JJJ  after update_aux...CC step')

            !% --- zero and small depth adjustment for elements
            call adjust_element_toplevel (CC)
                ! ! call util_utest_CLprint ('------- KKK  after adjust element CC (before 2nd step junction)')

            !% --- JUNCTION 2nd STEP
            if (N_nJM > 0) then 
                if (istep == 1) then
                    !% --- ensure JB interpweights for Q are forced for JB dominance: QUESTION -- SHOULD THIS BE TRUE OR NOT?
                    Npack => npack_elemP(ep_JB)
                    if (Npack > 0) then 
                        thisP => elemP(1:Npack, ep_JB)
                        call update_interpweights_JB (thisP, Npack, .false.)
                    end if

                else if (istep == 2) then 
                    !% --- conservative storage advance for junction
                    call junction_second_step ()
                    ! ! call util_utest_CLprint ('------- PPP  after junction second step')
                end if
            end if

            !% --- interpolate all data to faces
            sync all
            call face_interpolation(fp_noBC_IorS, .true., .true., .true., .false., .true.) 
                ! ! call util_utest_CLprint ('------- OOO  after face interp')

            !% saz 20230507
            if (N_diag > 0) then 
                !% --- update flowrates for aa diagnostic elements
                call diagnostic_by_type (ep_Diag, istep)  
                    ! ! call util_utest_CLprint ('------- BBB  after diagnostic')
    
                !% --- push the diagnostic flowrate data to faces -- true is upstream, false is downstream
                call face_push_elemdata_to_face (ep_Diag, fr_Flowrate, er_Flowrate, elemR, .true.)
                call face_push_elemdata_to_face (ep_Diag, fr_Flowrate, er_Flowrate, elemR, .false.)
                    ! ! call util_utest_CLprint ('------- CCC  after face_push_elemdata_to_face')
            end if

            !% ==============================================================
            !% --- face sync (saz 20230507)
            !%     sync all the images first. then copy over the data between
            !%     shared-identical faces. then sync all images again
            sync all
            call face_shared_face_sync (fp_noBC_IorS)
            sync all
            !% 
            !% ==============================================================

            !% --- update various packs of zeroDepth faces
            call pack_CC_zeroDepth_interior_faces ()
            if (N_nJM > 0) then 
                call pack_JB_zeroDepth_interior_faces ()
            end if
            sync all

            call pack_CC_zeroDepth_shared_faces ()  
            if (N_nJM > 0) then
                call pack_JB_zeroDepth_shared_faces ()  
            end if

            call face_zeroDepth (fp_CC_downstream_is_zero_IorS, &
                fp_CC_upstream_is_zero_IorS,fp_CC_bothsides_are_zero_IorS)
                ! ! call util_utest_CLprint ('------- PPP.01 after face zerodepth ')

            if (N_nJM > 0) then

                ! print *, ' '
                ! print *, 'down ',faceP(1:npack_faceP(fp_JB_downstream_is_zero_IorS),fp_JB_downstream_is_zero_IorS )
                ! print *, ' '
                ! print *, 'up   ',faceP(1:npack_faceP(fp_JB_upstream_is_zero_IorS),fp_JB_upstream_is_zero_IorS )
                ! print *, ' '
                ! print *, 'both ',faceP(1:npack_faceP(fp_JB_bothsides_are_zero_IorS),fp_JB_bothsides_are_zero_IorS )
                ! print *, ' '


                !% --- set face geometry and flowrates where adjacent element is zero
                !%     only applies to faces with JB on one side
                call face_zeroDepth (fp_JB_downstream_is_zero_IorS, &
                    fp_JB_upstream_is_zero_IorS,fp_JB_bothsides_are_zero_IorS)
                    ! ! call util_utest_CLprint ('------- PPP.02 after face zerodepth ')
            end if                

            !% --- enforce open (1) closed (0) "setting" value from EPA-SWMM
            !%     for all CC and Diag elements (not allowed on junctions)
            call face_flowrate_for_openclosed_elem (ep_CCDiag)
                ! ! call util_utest_CLprint ('------- QQQ after openclosed setting (before 1st step junction solution)')

            !% QUESTION: DO WE NEED ANOTHER SYNC HERE? OR CAN face_flowrate_for_openclosedL_elem
            !% BE MOVED UPWARDS IN STEPPING SO THAT IT GETS SYNCED?

            !% --- JUNCTION FIRST STEP
            if (istep == 1) then 
                !% --- Junction first step RK estimate
                !%     Note that this must be called in every image, including
                !%     those that do not have junctions as it contains a sync
                call junction_first_step ()
            end if
                 ! ! call util_utest_CLprint ('------- VVV.01  before adjust Vfilter CC after junction first step')

            !% --- Filter flowrates to remove grid-scale checkerboard
            call adjust_Vfilter_CC ()

                ! ! call util_utest_CLprint ('------- VVV.02  after adjust Vfilter CC')

            if (istep == 1) then 
                !% -- fluxes at end of first RK2 step are the conservative fluxes enforced
                !%    in second step
                call rk2_store_conservative_fluxes (ALL) 
                    ! all util_utest_CLprint ('------- WWW  after  step 1 store conservative fluxes all')
            end if

                ! ! call util_utest_CLprint ('------- YYY end of RK step')

        end do

         ! ! call util_utest_CLprint ('========== ZZZ  end of RK2 ============================')

        !% --- overall volume conservation
        volume2 = zeroR
        do ii=1,N_elem(1)
            !if (elemR(ii,er_Depth) > setting%ZeroValue%Depth) then
            if ((elemI(ii,ei_elementType) == JM) .or. (elemI(ii,ei_elementType) == CC)) then
                volume2 = volume2 + elemR(ii,er_Volume)
            end if
        end do
        Npack => npack_faceP(fp_BCup)
        if (Npack > 0) then 
            thisP => faceP(1:Npack,fp_BCup)
            inflowVolume  = sum(faceR(thisP,fr_Flowrate)) * setting%Time%Hydraulics%Dt
        end if
        Npack => npack_faceP(fp_BCdn)
        if (Npack > 0) then 
            thisP => faceP(1:Npack,fp_BCdn)
            outflowVolume  = sum(faceR(thisP,fr_Flowrate)) * setting%Time%Hydraulics%Dt
        end if

        ! !print *, ' '
        ! VolumeErrorCumulative = VolumeErrorCumulative + volume2 - volume1 - inflowVolume + outflowVolume
        ! write(*,"(A, i6, 6e12.4)") 'VOLUMES ',setting%Time%Step, volume1, volume2, inflowVolume, outflowVolume, &
        ! volume2 - volume1 - inflowVolume + outflowVolume, VolumeErrorCumulative

        

        !print *, ' '
    end subroutine rk2_toplevel_ETM
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    subroutine rk2_step_ETM_CC (istep)
        !%------------------------------------------------------------------
        !% Performs rk2 step for volume and velocity for CC elements
        !% using ETM method
        !%------------------------------------------------------------------
            integer, intent(in) :: istep
            integer, pointer    :: thisPackCol, Npack
            integer, pointer    :: FMpackCol, nFMpack, thisP(:)

            integer, parameter :: thisMethod = ETM
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------
            
            ! ! ! ! ! ! ! ! !  ! ! call util_utest_CLprint ('-------  WWW 00  start of rk2 Step CC')

        !% --- CONTINUITY
        thisPackCol => col_elemP(ep_CC_H)
        Npack       => npack_elemP(thisPackCol)
        if (Npack > 0) then
            thisP => elemP(1:Npack,thisPackCol)

            !% --- Compute net flowrates for CC as source termo
            call ll_continuity_netflowrate_CC (er_SourceContinuity, thisPackCol, Npack)


                !! ! call util_utest_CLprint ('-------  WWW 01  after ll_continuity_netflowrate')

            !% --- Solve for new volume
            call ll_continuity_volume_ETM (er_Volume, thisPackCol, Npack, istep)

                !! ! call util_utest_CLprint ('-------  WWW 02  after ll_continuity_volume')


            !% --- adjust extremely small volumes that might be been introduced
            call adjust_limit_by_zerovalues &
                (er_Volume, setting%ZeroValue%Volume, thisP, .true.)

                !! ! call util_utest_CLprint ('-------  WWW 03  after adjust_limit_by_zerovalues')

        end if  



        !% --- MOMENTUM
        thisPackCol => col_elemP(ep_CC_Q)
        Npack       => npack_elemP(thisPackCol)
        if (Npack > 0) then

            !% --- momentum K source terms for different methods for ETM
            call ll_momentum_Ksource_CC (er_Ksource, thisPackCol, Npack)

              !  print *, 'Ksource ',elemR(6,er_KSource)

            !% --- Common source for momentum on channels and conduits for ETM
            call ll_momentum_source_CC (er_SourceMomentum, thisPackCol, Npack)

                !! ! call util_utest_CLprint ('-------  WWW 06  after ll_momentum_source_cc')
 
              ! print *, 'Source ',elemR(6,er_SourceMomentum)

            !% --- Common Gamma for momentum on channels and conduits for  ETM
            !%     Here for all channels and conduits, assuming CM roughness
            call ll_momentum_gammaCM_CC (er_GammaM, thisPackCol, Npack)

              !  print *, 'Gamma  ',elemR(6,er_GammaM)

            !% --- handle force mains as Gamma terms
            !%     These overwrite the gamma from the CM roughness above
            if (setting%Solver%ForceMain%AllowForceMainTF) then 

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

              !  print *, 'Gamma2 ',elemR(6,er_GammaM)

            !% --- Advance flowrate to n+1/2 for conduits and channels with ETM
            call ll_momentum_solve_CC (er_Velocity, thisPackCol, Npack, thisMethod, istep)

              !  print *, 'Vel 1   ',elemR(6,er_Velocity)

                !  ! ! call util_utest_CLprint ('-------  WWW 15  after ll_momentum_solve_cc')

            !% --- velocity for ETM time march
            call ll_momentum_velocity_CC (er_Velocity, thisPackCol, Npack)

              !  print *, 'Vel 2   ',elemR(6,er_Velocity)

            !% --- prevent backflow through flapgates
            call ll_enforce_flapgate_CC (er_Velocity, thisPackCol, Npack)

             !   print *, 'flap   ',elemR(6,er_Velocity)

            !% --- enforce zero velocity on elements that began as ZeroDepth
            call ll_enforce_zerodepth_velocity (er_Velocity, thisPackCol, Npack)

                ! print *, 'Zero   ',elemR(6,er_Velocity)

        end if
        
    end subroutine rk2_step_ETM_CC
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rk2_step_ETM (istep)
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Performs a single rk2 step for ETM
        !% continuity on CC and JM
        !% momentum on CC
        !% velocity/flowrate on JB
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: istep
        integer :: tmType, Npack
        !%-----------------------------------------------------------------------------

        print *, 'OBSOLETE'
        stop 2098734
        ! !%       
        !     ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !  ! ! call util_utest_CLprint ('before rk2 continuity step etm')

        ! !% perform the continuity step of the rk2 for ETM CC and JM
        ! call rk2_continuity_step_ETM(istep)

        !     ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !  ! ! call util_utest_CLprint ('after rk2 continuity step etm')

        ! !% only adjust extremely small element volumes that have been introduced
        ! Npack = npack_elemP(ep_CCJM_H_ETM)
        ! if (Npack > 0) then
        !     call adjust_limit_by_zerovalues &
        !         (er_Volume, setting%ZeroValue%Volume/twentyR, elemP(1:Npack,ep_CCJM_H_ETM), .true.)
        ! end if
        !     ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !  ! ! call util_utest_CLprint ('after rk2 call to adjust limit by zero')

        ! !% perform the momentum step of the rk2 for ETM
        ! call rk2_momentum_step_ETM(istep)

            ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !  ! ! call util_utest_CLprint (' after rk2 call to rk2_momentum_step_ETM')

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
        !% Description:  THIS WILL BECOM OBSOLETE WITH IMPLICIT JUNCTION
        !% perform the continuity step of the rk2 for ETM
        !%-----------------------------------------------------------------------------
        integer, intent(in) :: istep
        integer, pointer :: thisPackCol, Npack
        logical :: isreset
        !%-----------------------------------------------------------------------------
        !%
        ! if (setting%Time%Step > 37466) then 
            ! print *, 'in rk2_continuity_step at top'
            ! print *, elemR(50,er_SourceContinuity), elemR(50,er_Velocity)
        ! end if

        print *, 'OBSOLETE'
        stop 4098734

        ! !% Baseline continuity source:
        ! !% Compute net flowrates for channels, conduits 
        ! thisPackCol => col_elemP(ep_CC_H)
        ! Npack => npack_elemP(thisPackCol)
        ! if (Npack > 0) then
        !     call ll_continuity_netflowrate_CC (er_SourceContinuity, thisPackCol, Npack)
        ! end if


        !     ! print *, 'in rk2_continuity_step after netflowrate CC'
        !     ! print *, elemR(50,er_SourceContinuity), elemR(50,er_Velocity)
 


        ! !% compute net flowrates for junction mains
        ! thisPackCol => col_elemP(ep_JM)
        ! Npack => npack_elemP(thisPackCol)
        ! if (Npack > 0) then
        !     call ll_continuity_netflowrate_JM (er_SourceContinuity, thisPackCol, Npack)
        ! end if


        !     ! print *, 'in rk2_continuity_step after netflowrate JM'
        !     ! print *, elemR(50,er_SourceContinuity), elemR(50,er_Velocity)


        ! !% Solve for volume in ETM step
        ! thisPackCol => col_elemP(ep_CCJM_H_ETM)
        ! Npack => npack_elemP(thisPackCol)
        ! if (Npack > 0) then
        !     call ll_continuity_volume_ETM (er_Volume, thisPackCol, Npack, istep)
        ! end if

        !     ! print *, 'in rk2_continuity_step after volume CCJM'
        !     ! print *, elemR(50,er_Volume), elemR(50,er_Volume_N0)


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
        ! THIS WILL BECOME OBSOLETE WITH IMPLICIT JUNCTION
        !%------------------------------------------------------------------
        !% Description:
        !% perform the momentum step of the rk2 for ETM (or ETM part of ETMAC)
        !% Computed on CC elements and velocity/flowrate for JB
        !%-------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: istep
            integer, pointer :: thisPackCol, Npack
            integer, pointer :: FMpackCol, nFMpack
            
            integer, parameter :: thisMethod = ETM
        !%-------------------------------------------------------------------
        !% Aliases
            thisPackCol => col_elemP(ep_CC_Q)
            Npack       => npack_elemP(thisPackCol)
        !%-------------------------------------------------------------------
        !%
        if (Npack > 0) then

                ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !  ! ! call util_utest_CLprint (' start of rk2_momentum_step_ETM')

                ! print *, ' '
                ! print *, 'in rk2_momentum_step at start'
                ! print *, elemR(61,er_Ksource), elemR(61,er_Velocity), elemR(61,er_Flowrate)
   
            !% --- momentum K source terms for different methods for ETM
            call ll_momentum_Ksource_CC (er_Ksource, thisPackCol, Npack)


                ! print *, 'in rk2_momentum_step at A'
                ! print *, elemR(61,er_Ksource), elemR(61,er_Velocity), elemR(61,er_Flowrate)

                ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !  ! ! call util_utest_CLprint (' after rk2 call ll_momentum_Ksource_CC')


            !% --- Common source for momentum on channels and conduits for ETM
            call ll_momentum_source_CC (er_SourceMomentum, thisPackCol, Npack)

                ! print *, 'in rk2_momentum_step at B'
                ! print *, elemR(61,er_SourceMomentum), elemR(61,er_HydRadius), elemR(61,er_ManningsN), elemR(61,er_Flowrate)

                ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !  ! ! call util_utest_CLprint (' after rk2 call ll_momentum_source_CC')

            !% --- Common Gamma for momentum on channels and conduits for  ETM
            !%     Here for all channels and conduits, assuming CM roughness
            call ll_momentum_gammaCM_CC (er_GammaM, thisPackCol, Npack)

                ! print *, 'in rk2_momentum_step at C'
                ! print *, elemR(61,er_GammaM), elemR(61,er_Velocity), elemR(61,er_Flowrate)      
            
                ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !  ! ! call util_utest_CLprint (' after rk2 call ll_momentum_gammaCM_CC')

            !% --- handle force mains as Gamma terms
            !%     These overwrites the gamma from the CM roughness above
            if (setting%Solver%ForceMain%AllowForceMainTF) then 
                !% --- surcharged Force main elements with Hazen-Williams roughness
                FMPackCol => col_elemP(ep_FM_HW_PSsurcharged)
                nFMpack   => npack_elemP(FMPackCol)
                if (nFMpack > 0) call ll_momentum_gammaFM_CC (er_GammaM, FMPackCol, nFMpack, HazenWilliams)

                !% --- surcharged Force Main elements with Darcy-Weisbach roughness
                FMPackCol => col_elemP(ep_FM_dw_PSsurcharged)
                nFMpack   => npack_elemP(FMPackCol)
                if (nFMpack > 0) call ll_momentum_gammaFM_CC (er_GammaM, FMPackCol, nFMpack, DarcyWeisbach)
            end if

                ! print *, 'in rk2_momentum_step at D'
                ! print *, elemR(61,er_GammaM), elemR(61,er_Velocity), elemR(61,er_Flowrate)



            !% --- add minor loss term to gamma for all conduits
            call ll_minorloss_friction_gamma_CC (er_GammaM, thisPackCol, Npack)   

                ! print *, 'in rk2_momentum_step at E'
                ! print *, elemR(61,er_GammaM), elemR(61,er_Velocity), elemR(61,er_Flowrate)

                ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !  ! ! call util_utest_CLprint (' after rk2 call ll_minorloss_friction_gamma_CC')
  
            !% --- Advance flowrate to n+1/2 for conduits and channels with ETM
            call ll_momentum_solve_CC (er_Velocity, thisPackCol, Npack, thisMethod, istep)


                ! print *, 'in rk2_momentum_step at F'
                ! print *, elemR(61,er_GammaM), elemR(61,er_Velocity), elemR(61,er_Flowrate)
            
                ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !  ! ! call util_utest_CLprint (' after rk2 call ll_momentum_solve_CC')


            !% --- velocity for ETM time march
            call ll_momentum_velocity_CC (er_Velocity, thisPackCol, Npack)


                ! print *, 'in rk2_momentum_step at G'
                ! print *, elemR(61,er_GammaM), elemR(61,er_Velocity), elemR(61,er_Flowrate)
                ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !  ! ! call util_utest_CLprint (' after rk2 call ll_momentum_velocity_CC')


            !% --- prevent backflow through flapgates
            call ll_enforce_flapgate_CC (er_Velocity, thisPackCol, Npack)

                ! print *, 'in rk2_momentum_step at H'
                ! print *, elemR(61,er_GammaM), elemR(61,er_Velocity), elemR(61,er_Flowrate)

                ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !  ! ! call util_utest_CLprint (' after rk2 call ll_enforce_flapgate_CC')

        end if

        !% --- update junction branches
        ! select case (setting%Junction%Method)
        !     case (Implicit0)
        !         !% no action
        !     case (Explicit1)
        !         call ll_flowrate_and_velocity_JB(ETM,istep)
        !     case (Explicit2)
        !         call ll_flowrate_and_velocity_JB_2(ETM,istep)
        !     case default
        !         print *, 'CODE ERROR: unexpected case default'
        ! end select

            !     print *, 'in rk2_momentum_step at I'
            !     print *, elemR(61,er_GammaM), elemR(61,er_Velocity), elemR(61,er_Flowrate)

    
            ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !  ! ! call util_utest_CLprint (' after ll_flowrate_and_velocity_JB')

    end subroutine rk2_momentum_step_ETM
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rk2_momentum_step_AC (istep)
        !%-----------------------------------------------------------------------------
        ! !% Description:
        ! !% Performs a single hydrology step
        ! !%-----------------------------------------------------------------------------
        integer, intent(in) :: istep
        ! integer, pointer    :: thisCol, Npack
        ! integer, parameter  :: thisMethod = AC
        ! !%-----------------------------------------------------------------------------
        ! !%
        ! !if (crashYN) return
        ! thisCol => col_elemP(ep_CC_Q_AC)
        ! Npack => npack_elemP(thisCol)
        ! if (Npack > 0) then
        !     !% momentum K source terms for different methods
        !     call ll_momentum_Ksource_CC (er_Ksource, thisCol, Npack)
        !     !% Source for momentum on channels and conduits
        !     call ll_momentum_source_CC (er_SourceMomentum, thisCol, Npack)
        !     !% Gamma for momentum on channels and conduits for Chezy-Manning roughness
        !     call ll_momentum_gammaCM_CC (er_GammaM, thisCol, Npack)
        !     !% additional momentum Gamma for AC time-march on channels and conduits
        !     call ll_momentum_add_gamma_CC_AC (er_GammaM, thisCol, Npack)
        !     !% additional momentum source terms for AC time-march on channels and conduits
        !     call ll_momentum_add_source_CC_AC (er_SourceMomentum, thisCol, Npack)
        !     !% AC elements advance flowrate to n+1(*) for conduits and channels
        !     call ll_momentum_solve_CC (er_Velocity, thisCol, Npack, thisMethod, istep)
        ! end if

        ! thisCol => col_elemP(ep_CC_Q_AC)
        ! Npack => npack_elemP(thisCol)
        ! if (Npack > 0) then
        !     !% velocity for AC time march
        !     call ll_momentum_velocity_CC (er_Velocity, thisCol,Npack)
        ! end if

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
        
        ! thisCol => col_elemP(ep_CCJB_eETM_i_fAC)
        ! Npack => npack_elemP(thisCol)


        ! if (Npack > 0) then
        !     !% temporary storage of n+1/2 data
        !     call ll_store_in_temporary (thisCol, Npack)

        !     !% extrapolation
        !     call ll_extrapolate_values (thisCol, Npack)

        !     !% update aux for extrapolated variables
        ! !        call update_auxiliary_variables_byPack (thisCol, Npack)
        ! end if

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
        print *, 'OBSOLETE '
        stop 298734
        ! thisCol => col_elemP(ep_CCJB_eETM_i_fAC)
        ! Npack   => npack_elemP(thisCol)

        ! if (Npack > 0) then
        !     !% temporary storage of n+1/2 data
        !     call ll_restore_from_temporary (thisCol, Npack)

        !     !% update aux for restored variables
        !     ! call update_auxiliary_variables_byPack (thisPackCol, Npack)
        ! end if

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
    subroutine rk2_store_conservative_fluxes (faceset)
        !%------------------------------------------------------------------
        !% Description:
        !% store the intermediate face flow rates in the Rk2 which are
        !% the conservative flowrate over the time step
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: faceset !% --- either ALL or JB
            integer, pointer    :: npack, thisF(:)
                ! integer, intent(in) :: whichTM
                ! integer, pointer    :: NpackE, NpackF, NpackJ, nBarrel(:)
                ! integer, pointer    :: thisColE,   thisColF, thisColJ, isbranch(:)
                ! integer, pointer    :: thisP(:), thisF(:), thisJ(:), fup(:), fdn(:)
                ! real(8), pointer    :: fQcons(:), fQ(:)
                ! integer :: ii
        !%------------------------------------------------------------------
        !% Preliminaries:
            ! select case (whichTM)
            ! case (ETM)
            !     thisColE  => col_elemP(ep_CC_ETM)
            !     thisColJ  => col_elemP(ep_JM_ETM)
            !     thisColF  => col_faceP(fp_Diag_NotJB_interior) !% changed 20230317
            ! case default
            !     print *, 'CODE ERROR: time march type unknown for # ', whichTM
            !     print *, 'which has key ',trim(reverseKey(whichTM))
            !     stop 398705
            ! end select
            
        !%------------------------------------------------------------------
        !% Aliases:
            ! NpackE => npack_elemP(thisColE)
            ! NpackF => npack_faceP(thisColF)
            ! NpackJ => npack_elemP(thisColJ)
            ! nBarrel=> elemI(:,ei_barrels)
            ! fup    => elemI(:,ei_Mface_uL)
            ! fdn    => elemI(:,ei_Mface_dL)
            ! fQcons => faceR(:,fr_Flowrate_Conservative)
            ! fQ     => faceR(:,fr_Flowrate)
            ! isbranch => elemSI(:,esi_JunctionBranch_Exists)
        !%------------------------------------------------------------------

        select case (faceset)
        case (ALL)
            faceR(:,fr_Flowrate_Conservative) = faceR(:,fr_Flowrate)
        case (CCDiag)
            npack => npack_faceP(fp_notJB_all)
            if (npack > 0) then 
                thisF => faceP(1:npack,fp_notJB_all)
                faceR(thisF,fr_Flowrate_Conservative) = faceR(thisF,fr_Flowrate)
            end if
        case (JBDiag)
            npack => npack_faceP(fp_JBorDiag_all)
            if (npack > 0) then 
                thisF => faceP(1:npack,fp_JBorDiag_all)
                faceR(thisF,fr_Flowrate_Conservative) = faceR(thisF,fr_Flowrate)
            end if
        case (JB)
            npack => npack_faceP(fp_JB_all)
            if (npack > 0) then
                thisF => faceP(1:npack,fp_JB_all)
                faceR(thisF,fr_Flowrate_Conservative) = faceR(thisF,fr_Flowrate)
            end if
        case default
            print *, 'CODE ERROR: unexpected case default'
            call util_crashpoint(8802772)
        end select

        ! print *, 'fQcons 00 ',fQcons(elemI(32,ei_Mface_uL))
        ! !% handle the flux faces of the time-marching elements
        ! if (NpackE > 0) then
        !     thisP  => elemP(1:NpackE, thisColE)
        !     print *, ' '
        !     print *, 'thisP'
        !     print *, thisP
        !     print *, ' '
        !     fQcons(fup(thisP)) = fQ(fup(thisP))
        !     fQcons(fdn(thisP)) = fQ(fdn(thisP))
        ! end if        

        ! print *, 'fQcons AA ',fQcons(elemI(32,ei_Mface_uL))

        ! !% handle flux faces of junctions
        ! if (NpackJ > 0) then
        !     thisJ  => elemP(1:NpackJ,thisColJ)
        !     do ii=1,max_branch_per_node,2
        !         fQcons(fup(thisJ+ii  )) = fQ(fup(thisJ+ii  )) * real(isbranch(thisJ+ii  ),8) * real(nBarrel(thisP+ii  ),8)
        !         fQcons(fdn(thisJ+ii+1)) = fQ(fdn(thisJ+ii+1)) * real(isbranch(thisJ+ii+1),8) * real(nBarrel(thisP+ii+1),8)
        !     end do
        ! end if

        ! print *, 'fQcons BB ',fQcons(elemI(32,ei_Mface_uL))
    
        ! !% handle the flux faces of the diagnostic elements
        ! if (NpackF > 0) then
        !     thisF  => faceP(1:NpackF, thisColF)
        !     fQcons(thisF) = fQ(thisF)
        ! end if

        ! print *, 'fQcons CC ',fQcons(elemI(32,ei_Mface_uL))

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
            ! select case (whichTM)
            ! case (ETM)
                thisColCC  => col_elemP(ep_CC)
                thisColJM  => col_elemP(ep_JM)
            ! case default
            !     print *, 'CODE ERROR: time march type not handled for # ', whichTM
            !     print *, 'which has key ',trim(reverseKey(whichTM))
            !     stop 398705
            ! end select    
        !%------------------------------------------------------------------  

        print *, 'DYNAMIC MANNINGS CANNOT BE USED'
        print *, 'change setting%Solver%ManningsN%useDynamicManningsN = .false.'    
        stop 3076612

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

