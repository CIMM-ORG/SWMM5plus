module runge_kutta2
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Provides the single time step of hydraulic time march
    !%
    !% Methods:
    !% Runge-Kutta 2-step is used for conduits and channels with
    !% a junction diagnostic balance approach
    !%==========================================================================

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
    use rk2_lowlevel
    use culvert_elements, only: culvert_toplevel
    use pack_mask_arrays
    use preissmann_slot
    use adjust
    use diagnostic_elements
    use utility, only: util_syncwrite
    use utility_crash
    use utility_unit_testing, only: util_utest_CLprint, util_utest_checkIsNan

    implicit none

    private

    public :: rk2_toplevel

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine rk2_toplevel ()
        !%------------------------------------------------------------------
        !% Description:
        !% single RK2 step for explicit time advance of SVE
        !%------------------------------------------------------------------
        !% Declarations:
            integer          :: istep, ii, kk, mm
            integer          :: fr_hAdj, fr_aAdj
            integer, pointer :: Npack, thisP(:), fup(:), fdn(:), thisJB(:)
            integer, pointer :: JMar(:), fAdj(:)
            
            real(8), pointer :: deltaH(:), energyQ(:), grav, dt
            real(8)          :: bsign
            real(8)          :: volume1, volume2, inflowVolume, outflowVolume

            logical :: isUpstreamFace, isUpstreamBranch

            character(64) :: subroutine_name = 'rk2_toplevel_ETM'
        !%------------------------------------------------------------------
        !% Preliminaries
        !% --- reset the overflow counter for this time level
            elemR(:,er_VolumeOverFlow)         = zeroR     
            elemR(:,er_VolumeArtificialInflow) = zeroR   
        !%-----------------------------------------------------------------
        !% Aliases
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
                    
        !% --- istep is the RK substep counter, initially set to zero
        !%     for preliminaries
        istep = zeroI

        !% --- Preliminary values for JM/JB elements
        !%     Note, this must be called even if no JM/JB on this image because 
        !%     the faces require synchronizing.
        call junction_preliminaries ()

        !%==================================  
        !% --- RK2 SOLUTION
        do istep = 1,2
 
            !% --- Half-timestep advance on CC for U and UVolume
            call rk2_step_CC (istep)  

            !% --- Update all CC aux variables
            !%     Note, these updates CANNOT depend on face values
            !%     Through geometry, this sets Preissmann Slot variables
            call update_auxiliary_variables_CC (                   &
                ep_CC, ep_CC_Open_Elements, ep_CC_Closed_Elements, &
                .true., .false., dummyIdx)

            !% --- zero and small depth adjustment for elements
            call adjust_element_toplevel (CC)

            !% --- JUNCTION 1st Step setup, 2nd Step compute
            if (N_nJM > 0) then 
                if (istep == 1) then
                    !% --- ensure JB interpweights for Q are forced for JB dominance
                    Npack => npack_elemP(ep_JB)
                    if (Npack > 0) then 
                        thisP => elemP(1:Npack, ep_JB)
                        call update_interpweights_JB (thisP, Npack, .false.)
                    end if
                else if (istep == 2) then 
                    !% --- conservative storage advance for junction, second step
                    call junction_second_step ()
                end if
            end if

            !% --- interpolate all data to faces
            sync all
            call face_interpolation(fp_noBC_IorS, .true., .true., .true., .false., .true.) 

            if (N_diag > 0) then 
                !% --- update flowrates for aa diagnostic elements
                call diagnostic_by_type (ep_Diag, istep)  
                !% --- push the diagnostic flowrate data to faces -- true is upstream, false is downstream
                call face_push_elemdata_to_face (ep_Diag, fr_Flowrate, er_Flowrate, elemR, .true.)
                call face_push_elemdata_to_face (ep_Diag, fr_Flowrate, er_Flowrate, elemR, .false.)
            end if

            !% --- face sync
            !%     sync all the images first. then copy over the data between
            !%     shared-identical faces. then sync all images again
            sync all
            call face_shared_face_sync (fp_noBC_IorS)
            sync all

            !% --- update various packs of zeroDepth faces for changes in depths
            call pack_CC_zeroDepth_interior_faces ()
            if (N_nJM > 0) then 
                call pack_JB_zeroDepth_interior_faces ()
            end if

            !% --- transfer zero depth faces between images
            sync all 
            call pack_CC_zeroDepth_shared_faces ()  
            if (N_nJM > 0) then
                call pack_JB_zeroDepth_shared_faces ()  
            end if

            !% --- set adjacent values for zerodepth faces for CC 
            call face_zeroDepth (fp_CC_downstream_is_zero_IorS, &
                fp_CC_upstream_is_zero_IorS,fp_CC_bothsides_are_zero_IorS)

            if (N_nJM > 0) then
                !% --- set face geometry and flowrates where adjacent element is zero
                !%     only applies to faces with JB on one side
                call face_zeroDepth (fp_JB_downstream_is_zero_IorS, &
                    fp_JB_upstream_is_zero_IorS,fp_JB_bothsides_are_zero_IorS)
            end if                

            !% --- enforce open (1) closed (0) "setting" value from EPA-SWMM
            !%     for all CC and Diag elements (not allowed on junctions)
            call face_flowrate_for_openclosed_elem (ep_CCDiag)

            !% HACK QUESTION: DO WE NEED ANOTHER SYNC HERE AND DATA TRANSFER? 
            !% OR CAN face_flowrate_for_openclosed_elem
            !% BE MOVED UPWARDS IN STEPPING SO THAT IT GETS SYNCED?

            !% --- JUNCTION -- first step compute
            if (istep == 1) then 
                !% --- Junction first step RK estimate
                !%     Note that this must be called in every image, including
                !%     those that do not have junctions as it contains a sync
                call junction_first_step ()
            end if

            !% --- Filter flowrates to remove grid-scale checkerboard
            call adjust_Vfilter_CC ()

            if (istep == 1) then 
                !% -- fluxes at end of first RK2 step are the conservative fluxes enforced
                !%    in second step
                call rk2_store_conservative_fluxes (ALL) 
            end if

        end do

        !% RETAIN FOR DEBUGGING
        !% --- overall volume conservation
        ! volume2 = zeroR
        ! do ii=1,N_elem(1)
        !     if ((elemI(ii,ei_elementType) == JM) .or. (elemI(ii,ei_elementType) == CC)) then
        !         volume2 = volume2 + elemR(ii,er_Volume)
        !     end if
        ! end do
        ! Npack => npack_faceP(fp_BCup)
        ! if (Npack > 0) then 
        !     thisP => faceP(1:Npack,fp_BCup)
        !     inflowVolume  = sum(faceR(thisP,fr_Flowrate)) * setting%Time%Hydraulics%Dt
        ! end if
        ! Npack => npack_faceP(fp_BCdn)
        ! if (Npack > 0) then 
        !     thisP => faceP(1:Npack,fp_BCdn)
        !     outflowVolume  = sum(faceR(thisP,fr_Flowrate)) * setting%Time%Hydraulics%Dt
        ! end if
        ! print *, ' '
        ! VolumeErrorCumulative = VolumeErrorCumulative + volume2 - volume1 - inflowVolume + outflowVolume
        ! write(*,"(A, i6, 6e12.4)") 'VOLUMES ',setting%Time%Step, volume1, volume2, inflowVolume, outflowVolume, &
        ! volume2 - volume1 - inflowVolume + outflowVolume, VolumeErrorCumulative




    end subroutine rk2_toplevel
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    subroutine rk2_step_CC (istep)
        !%------------------------------------------------------------------
        !% Description:
        !% Performs rk2 step for volume and velocity for CC elements
        !% using ETM method
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: istep
            integer, pointer    :: thisPackCol, Npack
            integer, pointer    :: FMpackCol, nFMpack, thisP(:)
        !%------------------------------------------------------------------

        !% --- CONTINUITY
        thisPackCol => col_elemP(ep_CC_H)
        Npack       => npack_elemP(thisPackCol)
        if (Npack > 0) then
            thisP => elemP(1:Npack,thisPackCol)

            !% --- Compute net flowrates for CC as source termo
            call ll_continuity_netflowrate_CC (er_SourceContinuity, thisPackCol, Npack)

            !% --- Solve for new volume
            call ll_continuity_volume_CC (er_Volume, thisPackCol, Npack, istep)

            !% --- adjust extremely small volumes that might be been introduced
            call adjust_limit_by_zerovalues &
                (er_Volume, setting%ZeroValue%Volume, thisP, .true.)

        end if  

        !% --- MOMENTUM
        thisPackCol => col_elemP(ep_CC_Q)
        Npack       => npack_elemP(thisPackCol)
        if (Npack > 0) then

            !% --- momentum K source terms for different methods for ETM
            call ll_momentum_Ksource_CC (er_Ksource, thisPackCol, Npack)

            !% --- Common source for momentum on channels and conduits for ETM
            call ll_momentum_source_CC (er_SourceMomentum, thisPackCol, Npack)

            !% --- Common Gamma for momentum on channels and conduits for  ETM
            !%     Here for all channels and conduits, assuming CM roughness
            call ll_momentum_gammaCM_CC (er_GammaM, thisPackCol, Npack)

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

            !% --- Advance flowrate to n+1/2 for conduits and channels with ETM
            call ll_momentum_solve_CC (er_Velocity, thisPackCol, Npack, istep)

            !% --- velocity for ETM time march
            call ll_momentum_velocity_CC (er_Velocity, thisPackCol, Npack)

            !% --- prevent backflow through flapgates
            call ll_enforce_flapgate_CC (er_Velocity, thisPackCol, Npack)

            !% --- enforce zero velocity on elements that began as ZeroDepth
            call ll_enforce_zerodepth_velocity (er_Velocity, thisPackCol, Npack)

        end if
        
    end subroutine rk2_step_CC
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
        !%------------------------------------------------------------------            
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

    end subroutine rk2_store_conservative_fluxes
!%   
!%==========================================================================
!% END OF MODULE
!%==========================================================================
end module runge_kutta2

