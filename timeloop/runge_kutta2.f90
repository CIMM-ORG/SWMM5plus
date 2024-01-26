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
    use culvert_elements, only: culvert_toplevel  !% NOT WORKING AS OF 20230912
    use pack_mask_arrays
    use preissmann_slot
    use adjust
    use diagnostic_elements
    use air_entrapment
    use utility_crash
    use utility_unit_testing, only: util_utest_CLprint, util_utest_checkIsNan

    implicit none

    private

    public :: rk2_toplevel

    integer :: printIdx = 49
    integer :: stepcut = 73594
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
            integer          :: istep, ii, kk
            integer, pointer :: Npack, thisP(:)
            
            real(8), pointer :: grav, dt
            real(8)          :: volume1, volume2, inflowVolume, outflowVolume
            real(8)          :: totalvolume, sumlocaldiff, localcons
            
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

        !% --- debug total volume conservation
        if (setting%Debug%isGlobalVolumeBalance) then
            volume1 = zeroR
            Npack => npack_elemP(ep_CCJM)
            if (Npack > 0) then 
                thisP => elemP(1:Npack,ep_CCJM)
                volume1 = sum(elemR(thisP,er_Volume))
            endif
        end if
                    
        !% --- istep is the RK substep counter, initially set to zero
        !%     for preliminaries
        istep = zeroI

             !call util_utest_CLprint('AAAA start RK2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

        !% --- Preliminary values for JM/JB elements
        !%     Note, this must be called even if no JM/JB on this image because 
        !%     the faces require synchronizing.
        call junction_preliminaries ()

            !call util_utest_CLprint('BBBB after junction preliminaries')
        
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

                !call util_utest_CLprint('DDDD after update auxiliary CC')

            !% --- zero and small depth adjustment for elements
            call adjust_element_toplevel (CC)
            
                !call util_utest_CLprint('EEEE after adjust element toplevel CC')

            !% --- JUNCTION 1st Step setup, 2nd Step compute
            if (N_nJM > 0) then 
                if (istep == 1) then
                    !% --- update JB interpweights without forcing
                    Npack => npack_elemP(ep_JB)
                    if (Npack > 0) then 
                        thisP => elemP(1:Npack, ep_JB)
                        call update_interpweights_JB (thisP, Npack, .false.)

                        !call util_utest_CLprint('FFF after update interpweights JB')
                    end if
                else if (istep == 2) then 
                    !call util_utest_CLprint('TTTT before junction second step')
                    !% --- conservative storage advance for junction, second step
                    call junction_second_step ()
                    !call util_utest_CLprint('UUUU after junction second step')
                end if
            end if  


            !% --- interpolate all data to faces
            sync all
            call face_interpolation(fp_noBC_IorS, .true., .true., .true., .false., .true.) 

            !call util_utest_CLprint('GGG after face interpolation')

            if (N_diag > 0) then 
                !% --- update flowrates for aa diagnostic elements
                call diagnostic_by_type (ep_Diag, istep)  
                !% --- push the diagnostic flowrate data to faces -- true is upstream, false is downstream
                call face_push_elemdata_to_face (ep_Diag, fr_Flowrate, er_Flowrate, elemR, .true.)
                call face_push_elemdata_to_face (ep_Diag, fr_Flowrate, er_Flowrate, elemR, .false.)
            end if

            !call util_utest_CLprint('HHH after diagnostic')

            !% --- face sync
            !%     sync all the images first. then copy over the data between
            !%     shared-identical faces. then sync all images again
            sync all
            call face_shared_face_sync_single (fp_Diag_IorS,fr_Flowrate)
            sync all

            !% --- update face velocities after sync changes areas and flowrates
            call face_update_velocities (fp_Diag_IorS)

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

            !% --- face sync
            !%     sync all the images first. then copy over the data between
            !%     shared-identical faces. then sync all images again
            sync all
            call face_shared_face_sync (fp_noBC_IorS, [fr_flowrate,fr_Velocity_d,fr_Velocity_u])
            sync all

            !call util_utest_CLprint('PPPP before junction first step')

            !% --- JUNCTION -- first step compute
            if (istep == 1) then 
                !% --- Junction first step RK estimate
                !%     Note that this must be called in every image, including
                !%     those that do not have junctions as it contains a sync
                call junction_first_step ()

                !call util_utest_CLprint('QQQQ after junction first step')
            end if

            !% --- Filter flowrates to remove grid-scale checkerboard
            call adjust_Vfilter (istep)

            !call util_utest_CLprint('RRRR after V filter')

            if (istep == 1) then 
                !% -- fluxes at end of first RK2 step are the conservative fluxes enforced
                !%    in second step
                call rk2_store_conservative_fluxes (ALL) 

                !call util_utest_CLprint('SSSS end of RK2 first step')
            else 
                !%  --- no action 
            end if

            !% Air entrapment modeling
            if (setting%AirTracking%UseAirTrackingYN) then
                call air_entrapment_toplevel (istep)
            end if 

        end do

        !call util_utest_CLprint('ZZZZ end RK2')

        !% HACK --- this needs to be setup for multiple images and moved to the utility_debug
        if (setting%Debug%isGlobalVolumeBalance) then
            !% --- overall volume conservation
            !% --- initialization
            elemR(:,er_Temp01) = zeroR
            inflowVolume  = zeroR
            outflowVolume = zeroR
            sumlocaldiff  = zeroR
            !% --- get all in-line inflows
            Npack => npack_faceP(fp_BCup)
            if (Npack > 0) then 
                thisP => faceP(1:Npack,fp_BCup)
                inflowVolume  = inflowVolume + sum(faceR(thisP,fr_Flowrate)) * setting%Time%Hydraulics%Dt
            end if
            !% --- get all outfall outflows
            Npack => npack_faceP(fp_BCdn)
            if (Npack > 0) then 
                thisP => faceP(1:Npack,fp_BCdn)
                outflowVolume  = outflowVolume + sum(faceR(thisP,fr_Flowrate_Conservative)) * setting%Time%Hydraulics%Dt
            end if
            !% --- net inflows due to lateral (negatives are outflows)
            Npack => npack_elemP(ep_CCJM)
            thisP => elemP(1:Npack,ep_CCJM)
            volume2 = sum(elemR(thisP,er_Volume))
            inflowVolume = inflowVolume + sum(elemR(thisP,er_FlowrateLateral)) * setting%Time%Hydraulics%Dt
            !% --- compute local volume conservation for CC
            Npack => npack_elemP(ep_CC)
            thisP => elemP(1:Npack,ep_CC)
            elemR(thisP,er_Temp01) = setting%Time%Hydraulics%Dt               &
                * (  elemR(thisP,er_FlowrateLateral)                          &
                   + faceR(elemI(thisP,ei_Mface_uL),fr_Flowrate_Conservative) &
                   - faceR(elemI(thisP,ei_Mface_dL),fr_Flowrate_Conservative) &
                  )                                                           &
                - (elemR(thisP,er_Volume) - elemR(thisP,er_Volume_N0))
            sumlocaldiff = sumlocaldiff + sum(elemR(thisP,er_Temp01)) 
            !% --- compute local volume conservation for JM 
            !%     start with lateral inflows
            Npack =>   npack_elemP(ep_JM)
            thisP => elemP(1:Npack,ep_JM)
            elemR(thisP,er_Temp01) = elemR(thisP,er_FlowrateLateral) * setting%Time%Hydraulics%Dt &
                 - (elemR(thisP,er_Volume) - elemR(thisP,er_Volume_N0))
            !% --- next are branch flows, so we shift the packed array
            Npack =>   npack_elemP(ep_JB)
            thisP => elemP(1:Npack,ep_JB)
            !% --- accumulate flow volumes for upstream JB
            where ((elemSI(thisP,esi_JB_IsUpstream) .eq. oneI) .and. (elemSI(thisP,esi_JB_Exists) .eq. oneI))
                elemR  (elemSI(thisP,esi_JB_Main_Index),er_Temp01)     &
                = elemR(elemSI(thisP,esi_JB_Main_Index),er_Temp01)     &
                + faceR( elemI(thisP,ei_Mface_uL),fr_Flowrate_Conservative) * setting%Time%Hydraulics%Dt
            endwhere
            !% ---- accumulate flow volumes for downstream JB
            where ((elemSI(thisP,esi_JB_IsUpstream) .eq. zeroI) .and. (elemSI(thisP,esi_JB_Exists) .eq. oneI))
                elemR  (elemSI(thisP,esi_JB_Main_Index),er_Temp01)       &
                = elemR(elemSI(thisP,esi_JB_Main_Index),er_Temp01)     &
                + faceR( elemI(thisP,ei_Mface_dL),fr_Flowrate_Conservative) * setting%Time%Hydraulics%Dt
            endwhere
            !% --- create sum for JM
            Npack =>   npack_elemP(ep_JM)
            thisP => elemP(1:Npack,ep_JM)
            sumlocaldiff = sumlocaldiff + sum(elemR(thisP,er_Temp01)) 
            
            totalvolume  = max(volume1, volume2)
            if (totalvolume > oneR) then
                !% --- use normalized volume for large volumes
                setting%Debug%GlobalVolumeBalance = setting%Debug%GlobalVolumeBalance & 
                + (volume2 - volume1 - inflowVolume + outflowVolume)  / totalvolume
            else
                !% --- use raw values
                setting%Debug%GlobalVolumeBalance = setting%Debug%GlobalVolumeBalance & 
                    + (volume2 - volume1 - inflowVolume + outflowVolume)  
            end if

            ! !% DEBUG PRINTING DO NOT DELETE
            !  print *, setting%Time%Step, volume1, volume2, volume2 - volume1, &
            !       inflowVolume, outflowVolume,                                &
            !       volume2 - volume1 - inflowVolume + outflowVolume,           &
            !       setting%Debug%GlobalVolumeBalance 

            ! print *, 'Local conservation'
            ! localcons = zeroR
            ! do ii=1,N_elem(this_image())
            !     if ((elemI(ii,ei_elementType) .eq. CC) .or. (elemI(ii,ei_elementType) .eq. JM)) then
            !         localcons = localcons + elemR(ii,er_Temp01)
            !         !print *, ii, elemR(ii,er_Temp01)
            !         write(*,"(i4,12e12.4)") ii, elemR(ii,er_Volume), elemR(ii,er_Volume_N0), elemR(ii,er_Volume) -elemR(ii,er_Volume_N0)
            !     end if
                
            ! end do   
            ! print *, 'local cons ',localcons
        end if


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

            ! !call util_utest_CLprint('inside bbbb --------------------')

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

        ! !call util_utest_CLprint('inside cccc --------------------')
        
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
                print *, 'CODE ERROR unexpected case default'
                call util_crashpoint(8802772)
        end select

    end subroutine rk2_store_conservative_fluxes
!%   
!%==========================================================================
!% END OF MODULE
!%==========================================================================
end module runge_kutta2

