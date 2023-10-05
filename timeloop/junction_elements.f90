module junction_elements
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Computes junction head and flowrate
    !%
    !% Methods:
    !% balancing of dQ/dH for branches, overflow and storage
    !%==========================================================================
    use define_globals
    use define_keys
    use define_indexes
    use define_xsect_tables
    use define_settings, only: setting
    use adjust
    use diagnostic_elements, only: diagnostic_by_type
    use face
    use geometry
    ! use lowerlevel_junction
    use junction_lowlevel
    use geometry_lowlevel, only: llgeo_head_from_depth_pure
    use pack_mask_arrays
    use preissmann_slot
    use update
    use utility_crash, only: util_crashpoint

    use utility_unit_testing, only: util_utest_CLprint

    implicit none

    private
    
    public :: junction_preliminaries
    public :: junction_main_volume_advance
    public :: junction_first_step
    public :: junction_second_step

    integer :: printJM =30
    integer :: printJB =31
    integer :: stepCut = 52100

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine junction_preliminaries ()
        !%-----------------------------------------------------------------
        !% Description
        !%-----------------------------------------------------------------
        !% Provides computations for junctions prior to RK2 time stepping
        !% to (1) obtain consistency with diagnostic elements, (2) estimate
        !% the flowrate/velocity in the junction main, (3) use an energy
        !% equation argument for nominal outflow flowrates, and (4) compute
        !% the dQ/dH needed for the first step junction solution. 
        !%-----------------------------------------------------------------
        !% Declarations
            integer(kind=8), pointer :: step
        !%-----------------------------------------------------------------
            step => setting%Time%Step
        !%-----------------------------------------------------------------

        !% --- ensure Slot is correctly initialized
        if (N_nJM > 0) then 
            call lljunction_main_slotwidth (ep_JM)
        end if

        !% --- flowrate/velocity of the JM
        if (N_nJM > 0) then
            call lljunction_main_velocity (ep_JM)

                ! call util_utest_CLprint ('------- jjj.01  after lljunction_main_velocity')

            call lljunction_main_energyhead (ep_JM)

            ! call util_utest_CLprint ('------- jjj.02  after lljunction_main_energyhead')
        end if

        !% --- push inflows on CC upstream or downstream of JB elem to CC/JB face
        if (N_nJM > 0) then 
            call lljunction_push_inflows_from_CC_to_JB_face ()

            !% --- update the velocities for the new inflows.
            call face_update_velocities (fp_JB_IorS)

                ! call util_utest_CLprint ('------- jjj.03  after lljunction_push_inflowCC_flowrates_to_face')
        end if
           
        !% AT THIS POINT: JB-ADJACENT FACES NOW CONTAIN EITHER 
        !% (1) Diag fluxes or (2) CC inflows or (3) old data.
        !% However, the JB element fluxes are inconsistent with faces
        !% Need to resolve shared faces and then push face data to JB

        !% ==============================================================
        !% --- face sync 
        !%     sync all the images first. then copy over the data between
        !%     shared-identical faces. then sync all images again
        !%     This is needed so that we can pull face data to JB elements
        !%     when the image boundary is a JB face
        sync all
        call face_shared_face_sync (fp_JB_IorS,[fr_Flowrate,fr_Velocity_d,fr_Velocity_u])
        sync all
        !% ==============================================================

        !% --- ensure that all JB are consistent with adjacent face before the
        !%     energy equation is invoked for outflow
        if (N_nJM > 0) then 
            call face_pull_facedata_to_JBelem (ep_JM, fr_Flowrate,   er_Flowrate, .true.)
            call face_pull_facedata_to_JBelem (ep_JM, fr_Flowrate,   er_Flowrate, .false.)
            call face_pull_facedata_to_JBelem (ep_JM, fr_Velocity_d, er_Velocity, .true.)
            call face_pull_facedata_to_JBelem (ep_JM, fr_Velocity_u, er_Velocity, .false.)
        end if

        ! call util_utest_CLprint ('------- jjj.05 after face_pull_facedata')

        !% --- store junction-adjacent element data on face so that no-neighbor principal
        !%     is not violated.
        !%     QUESTION -- SHOULD THIS BE IN THE RK ITERATION FOR UPDATES?
        !%     ANSWER: NO, as long as the second step junction solution is NOT the backwards euler.
        !%     TO BE MOVED TO junction_branch_adjacent 
        if (N_nJM > 0) then 
            call lljunction_push_adjacent_CC_elemdata_to_face ()
        end if

        !% --- push JB adjacent diag data to faces
        if (npack_elemP(ep_Diag_JBadjacent) > 0) then
            call face_push_diag_adjacent_data_to_face (ep_Diag_JBadjacent)
        end if

        ! call util_utest_CLprint ('------- jjj.06 after lljunction_push_adjacent')
        !% ==============================================================
        !% --- face sync
        !%     sync all the images first. then copy over the data between
        !%     shared-identical faces. then sync all images again
        sync all
        call face_shared_face_sync (fp_JB_IorS,[fr_Head_Adjacent,fr_EnergyHead_Adjacent,         &
                                    fr_Topwidth_Adjacent,fr_Length_Adjacent, fr_Zcrest_Adjacent, &
                                    fr_Velocity_Adjacent,fr_Froude_Adjacent,fr_Depth_Adjacent,   &
                                    fr_dQdH_Adjacent])
        sync all
        !% 
        !% ==============================================================

 

        ! call util_utest_CLprint ('------- jjj.07 after face_push_diag')

         !% ==============================================================
        !% --- face sync
        !%     sync all the images first. then copy over the data between
        !%     shared-identical faces. then sync all images again
        ! sync all
        ! call face_shared_face_sync (fp_JB_IorS,[fr_Zcrest_Adjacent,fr_dQdH_Adjacent,fr_EnergyHead_Adjacent])
        ! sync all
        !% 
        !% ==============================================================

        !% --- JB ENERGY EQUATION compute flows/velocities on JB/CC outflow elements/faces from 
        !%     energy equation (does not affect JB with inflows or diagnostic adjacent )
        !%     TO BE MOVED TO junction_branch_element_flowrates?
        if (N_nJM > 0) then 
            call lljunction_branch_energy_outflow ()
            ! call util_utest_CLprint ('------- jjj.08 after lljunction_branch_energy_outflow')
        end if

        !% --- force the changed JB element flowrate and velocity alues to the faces for upstream (true)
        !%     and downstream (false) branches.
        call face_push_JBelem_to_face (ep_JM, fr_Flowrate,   er_Flowrate, .true.)
        call face_push_JBelem_to_face (ep_JM, fr_Flowrate,   er_Flowrate, .false.)
        call face_push_JBelem_to_face (ep_JM, fr_Velocity_d, er_Velocity, .true.)
        call face_push_JBelem_to_face (ep_JM, fr_Velocity_d, er_Velocity, .false.)
        call face_push_JBelem_to_face (ep_JM, fr_Velocity_u, er_Velocity, .true.)
        call face_push_JBelem_to_face (ep_JM, fr_Velocity_u, er_Velocity, .false.)
    
        !% AT THIS POINT: We now have JB elements and faces that are consistent. The Diag-adjacent and
        !% inflows have had the face values pushed to the JB, and JB/CC faces having the energy
        !% based flowrate/velocity from the element pushed to faces. Need to sync faces to
        !% ensure consistency across processors

        !% ==============================================================
        !% --- face sync
        !%     sync all the images first. then copy over the data between
        !%     shared-identical faces. then sync all images again
        sync all
        call face_shared_face_sync (fp_JB_IorS,[fr_Flowrate,fr_Velocity_u,fr_Velocity_d])
        sync all
        !% 
        !% ==============================================================

        !% --- store the junction dQdH used in Backwards Euler
        if (N_nJM > 0) then            
            call lljunction_branch_dQdH ()  
                ! call util_utest_CLprint ('------- jjj.09 after lljunction_branch_dQdH')
        end if

    end subroutine junction_preliminaries
!%
!%==========================================================================
!%==========================================================================
!%   
    subroutine junction_main_volume_advance (epCol,Npack)
        !%------------------------------------------------------------------
        !% Description
        !% computes a volume-conservative advance using the flux values
        !% from faceR(:,fr_FlowrateConservative)
        !% Note this begins from the N0 volume and uses the RK1 step values
        !% to provide the updates
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: epCol, Npack

            integer, pointer :: thisP(:)
            real(8), pointer :: dt

            integer :: mm, kk
        !%------------------------------------------------------------------
        !% Preliminaries
            if (Npack < 1) return
        !%------------------------------------------------------------------
        !% Aliases
            thisP => elemP(1:Npack,epCol)
            dt    => setting%Time%Hydraulics%Dt
        !%------------------------------------------------------------------
        !% cycle through the junctions
        do mm=1,Npack

            !% --- volume change due to lateral inflow
            elemR(thisp(mm),er_Volume)                        &
                = elemR(thisp(mm),er_Volume_N0)               &
                + elemR(thisP(mm),er_FlowrateLateral) * dt

            !% --- volume change due to overflow
            elemR(thisp(mm),er_Volume)                        &
                = elemR (thisp(mm),er_Volume)                 &
                + elemSR(thisP(mm),esr_JM_OverflowPondingRate) * dt      

            !% --- volume change due to overflow or ponding
            select case (elemSI(thisP(mm),esi_JM_OverflowType))
                case (NoOverflow)
                    !% no action 
                    
                case (OverflowWeir,OverflowOrifice)

                    ! if (setting%Time%Step > 54165)  then 
                    !     print *, 'Volume Overflow BBB', elemR(109,er_VolumeOverFlow)
                    ! end if

                    !% --- VolumeOverFlow > 0 is removing volume
                    ! elemR(thisp(mm),er_Volume) = elemR(thisp(mm),er_Volume) &
                    !         - elemR(thisP(mm),er_VolumeOverflow)

                    ! elemR(thisP(mm),er_VolumeOverFlowTotal) = elemR(thisP(mm),er_VolumeOverFlow) &
                    !     + elemR(thisP(mm),er_VolumeOverFlowTotal)

                    elemR(thisP(mm),er_VolumeOverFlowTotal)                  &
                        = elemR (thisP(mm),er_VolumeOverFlowTotal)           &
                        - elemSR(thisP(mm),esr_JM_OverflowPondingRate) * dt  


                case (PondedWeir,PondedOrifice)

                    !% --- VolumePonded > 0 is removing volume 
                    !%     VolumePonded < 0 is adding volume
                    ! elemR(thisp(mm),er_Volume) = elemR(thisp(mm),er_Volume) &
                    !         - elemR(thisP(mm),er_VolumePonded)

                    ! elemR(thisp(mm),er_Volume) = elemR(thisp(mm),er_Volume) &
                    !         + elemSR(thisP(mm),esr_JM_OverflowPondingRate) * dt

                    ! elemR(thisP(mm),er_VolumePondedTotal) = elemR(thisP(mm),er_VolumePonded) &
                    !     + elemR(thisP(mm),er_VolumePondedTotal)

                    elemR(thisP(mm),er_VolumePondedTotal)                   &
                        = elemR (thisP(mm),er_VolumePondedTotal)            &
                        - elemSR(thisP(mm),esr_JM_OverflowPondingRate) * dt  

                case default
                    print *, 'CODE ERROR unexpected case default'
                    call util_crashpoint(772223)
            end select

            !% --- volume change due to branch flows
            do kk=1,max_branch_per_node
                if (elemSI(thisP(mm)+kk,esi_JB_Exists) .ne. oneI) cycle
                if (mod(kk,2) == 0) then 
                    !% --- downstream branch
                    elemR(thisP(mm),er_Volume) = elemR(thisp(mm),er_Volume) &
                        - faceR(elemI(thisP(mm)+kk,ei_Mface_dL),fr_Flowrate_Conservative) * dt
                else 
                    !% --- upstream branch
                    elemR(thisP(mm),er_Volume) = elemR(thisp(mm),er_Volume) &
                        + faceR(elemI(thisP(mm)+kk,ei_Mface_uL),fr_Flowrate_Conservative) * dt
                end if
            end do

        end do

        ! print *, ' '
        ! print *, 'SINGLE POINT CHECK in junction_main_volume_advance'
        ! print *, -(elemR(116,er_Volume) - elemR(116,er_Volume_N0)) & 
        !          + elemR(116,er_FlowrateLateral)  * setting%Time%Hydraulics%Dt &
        !          + faceR(elemI(117,ei_Mface_uL),fr_Flowrate_Conservative) * setting%Time%Hydraulics%Dt &
        !          - faceR(elemI(118,ei_Mface_dL),fr_Flowrate_Conservative) * setting%Time%Hydraulics%Dt

    end subroutine junction_main_volume_advance
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine junction_first_step ()
        !%------------------------------------------------------------------
        !% Description:
        !% Provides first estimate of junction flows and storage
        !%------------------------------------------------------------------

        ! print *, '  '
        ! print *, 'AAAA ',elemR(101,er_Flowrate)
        ! print *, trim(reverseKey(elemI(101,ei_elementType)))
        ! print *, elemI(101,ei_Mface_dL)

        ! call util_utest_CLprint('aaaa top of junction 1st step -------------------------------------')
        ! print *, ' '

        if (N_nJM > 0) then 
            !% --- forces JB elem Q to faces (overriding the interpolation)
            !%     computes new JM Volume, Head, JB fluxes, JB DeltaQ
            !%     for conservation and dH change in head 
            !%     assigns new JB and JM aux values
            call junction_toplevel(1)

            ! call util_utest_CLprint('zzzz after junction toplevel -------------------------------------')
            ! print *, 'BBBB ',elemR(101,er_Flowrate)

        end if       
        !% ==============================================================
        !% --- face sync
        !%     sync all the images first. then copy over the data between
        !%     shared-identical faces. then sync all images again
        !%     This ensures faces between JB and adjacent element are
        !%     identical when face is shared between processors
        sync all
        call face_shared_face_sync (fp_JB_IorS,[fr_Flowrate,fr_DeltaQ])
        sync all
        !% 
        !% ==============================================================

        !% --- update velocities for new flowrates on JB faces
        call face_update_velocities (fp_JB_IorS)

        ! print *, 'DDDD ',elemR(101,er_Flowrate)

        !% brh20230927 not sure that this is needed
        ! if (N_nJM > 0) then 
        !     !% --- Adjust JB-adjacent CC elements using fr_DeltaQ flux changes
        !     !%     This fixes conservative flowrate, volume, velocity for upstream (true)
        !     !%     and downstream (false) branches. Note that flowrate is already
        !     !%     fixed in face_force_JBelem_to_face. Also calls update_auxiliary_data_CC
        !     !%     for associated geometry data updates
        !     call lljunction_CC_for_JBadjacent (ep_CC_UpstreamOfJunction,   1, .true.)
        !     call lljunction_CC_for_JBadjacent (ep_CC_DownstreamOfJunction, 1, .false.)

        !     !% --- NOTE we do not reset diagnostic faces because we cannot make them consistent
        !     !%     on both sides without violating the no-neighbor principle.

        ! end if

        ! print *, 'EEEE ',elemR(101,er_Flowrate)

        ! !% ==============================================================
        ! !% --- face sync
        ! !%     sync all the images first. then copy over the data between
        ! !%     shared-identical faces. then sync all images again
        ! !%     This ensures faces between JB and adjacent element are
        ! !%     identical when face is shared between processors
        ! sync all
        ! call face_shared_face_sync_single (fp_JB_IorS,fr_Flowrate_Conservative)
        ! sync all
        ! !% 
        ! !% ==============================================================

        !% --- these calls are outside of the if (N_nJM) statement to prevent any race conditions
        !% --- update various packs of zeroDepth faces
        call pack_JB_zeroDepth_interior_faces ()

        sync all
        call pack_JB_zeroDepth_shared_faces ()  !% HACK STUB ROUTINE NOT COMPLETE
        sync all

        !% --- set face geometry and flowrates where adjacent element is zero
        !%     only applies to faces with JB on one side
        call face_zeroDepth (fp_JB_downstream_is_zero_IorS, &
            fp_JB_upstream_is_zero_IorS,fp_JB_bothsides_are_zero_IorS)

            ! print *, 'FFFF',elemR(101,er_Flowrate)
        
    end subroutine junction_first_step
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine junction_second_step ()
        !%------------------------------------------------------------------
        !% Description
        !% Conservative flux and storage step for junctions
        !%------------------------------------------------------------------
            integer, pointer :: Npack, thisP(:)
            integer :: JMidx
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------

        Npack => npack_elemP(ep_JM)
        if (Npack > 0) then
            thisP => elemP(1:Npack,ep_JM)

            ! if (setting%Time%Step > stepCut) then
            !     JMidx = 30 
            !     print *, ''
            !     print *, 'in junction 2nd step before'
            !     print *, 'volumes ',elemR(JMidx,er_Volume), elemR(JMidx,er_Volume_N0)
            !     print *, 'Q vol up',elemR(JMidx+1,er_Flowrate) * setting%Time%Hydraulics%Dt, &
            !              faceR(elemI(JMidx+1,ei_Mface_uL),fr_Flowrate_Conservative) * setting%Time%Hydraulics%Dt
            !     print *, 'Q vol dn',elemR(JMidx+2,er_Flowrate) * setting%Time%Hydraulics%Dt, &
            !              faceR(elemI(JMidx+2,ei_Mface_dL),fr_Flowrate_Conservative) * setting%Time%Hydraulics%Dt
            !     print *, 'Qlat    ',elemR(JMidx,er_FlowrateLateral)  * setting%Time%Hydraulics%Dt
            !     print *, 'sum flow',&
            !        faceR(elemI(JMidx+1,ei_Mface_uL),fr_Flowrate_Conservative) * setting%Time%Hydraulics%Dt &
            !      - faceR(elemI(JMidx+2,ei_Mface_dL),fr_Flowrate_Conservative) * setting%Time%Hydraulics%Dt &
            !      + faceR(JMidx+1,er_FlowrateLateral) * setting%Time%Hydraulics%Dt
            !     print *, ' '
            ! end if

            ! if (setting%Time%Step > 52100) then 
            !     print *, ' '
            !     print *, 'in junction 2nd step before'
            !     print *, 'volume ',elemR(30,er_Volume), elemR(30,er_Volume_N0)
            !     print *, ' '
            ! end if


            !% --- new junction volume from conservative face fluxes
            call junction_main_volume_advance (ep_JM, Npack)

            ! if (setting%Time%Step > 52100) then 
            !     print *, ' '
            !     print *, 'in junction 2nd step'
            !     print *, 'volume ',elemR(30,er_Volume), elemR(30,er_Volume_N0)
            !     print *, ' '
            ! end if

            !% --- new junction plan area (non-surcharged functional, tabular storage only)
            call geo_plan_area_from_volume_JM (elemPGetm, npack_elemPGetm, col_elemPGetm)
            
            !% --- compute slots based on solved volume
            !%     includes the JunctionMain_Surcharge_Plan_Area
            call slot_JM (ep_JM, Npack)
            
            !% --- new junction depth 
            !%     NOTE: THIS USES storage_implied_depth_from_volume
            !%     that limits depth based on fulldepth, Slot is added back in slot_JM_adjustments
            call geo_depth_from_volume_JM (elemPGetm, npack_elemPGetm, col_elemPGetm)
        
            !% --- new JM head, ellDepth and area
            elemR(thisP,er_Head)     = llgeo_head_from_depth_pure (thisP,elemR(thisP,er_Depth))
            elemR(thisP,er_EllDepth) = elemR(thisP,er_Depth)
            elemR(thisP,er_Area)     = elemR(thisP,er_Depth) * sqrt(elemSR(thisP,esr_Storage_Plan_Area))

            !% --- add the Preissmann slot depths back to head 
            call slot_JM_adjustments (ep_JM, Npack)

            ! print *, ' '
            ! print *, 'SlotDepth after JM adjust ',elemR(printJM,er_SlotDepth)
            ! print *, 'head                      ',elemR(printJM,er_Head)
            ! print *, ' '

            !% --- adjust JM for small or zero depth
            call adjust_element_toplevel (JM)

            ! print *, 'JB head before'
            ! print *, elemR(169,er_Head), elemR(170,er_Head)
            ! print *, ' '

            !% --- assign JB values based on new JM head
            call geo_assign_JB_from_head (ep_JM) !% HACK  revise using ep_JB

            ! print *, 'JB head after'
            ! print *, elemR(169,er_Head), elemR(170,er_Head)
            ! print *, ' '

            !% --- Preissmann slot computations
            call slot_JB_computation (ep_JM)

            ! print *, 'SlotDepth after JB',elemR(169,er_SlotDepth),elemR(170,er_SlotDepth)
            ! print *, ' '

            !% --- adjust JB for small or zero depth
            call adjust_element_toplevel (JB)

            ! print *, 'Head end of 2nd step'
            ! print *, elemR(169,er_Head), elemR(168,er_Head), elemR(170,er_Head)
            ! print *, ' '
            
        end if

        !% --- auxiliary variables update
        !%     replaces update_auxiliary_variables_JMJB
        !%     NOTE: .true. in interpweight call forces Q weight on JB to minimum, which
        !%     means that face interpolation will have JB values
        !%     dominate over adjacent CC values, but will be
        !%     simple averaging with adjacent Diag Q values -- HERE WE USE TRUE IS THIS CORRECT?
        Npack => npack_elemP(ep_JB)
        if (Npack > 0) then 
            thisP => elemP(1:Npack, ep_JB)
            call update_Froude_number_element (thisP) 
            call update_wavespeed_element (thisP)
            call update_interpweights_JB (thisP, Npack, .false.)  !% 20230904brh change to FALSE
        end if

        !% --- wave speed, Froude number on JM
        Npack => npack_elemP(ep_JM)
        if (Npack > 0) then
            thisP => elemP(1:Npack, ep_JM)
            call update_wavespeed_element (thisP)
            call update_Froude_number_element (thisP) 
        end if

    end subroutine junction_second_step
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    subroutine junction_toplevel(istep)
        !%------------------------------------------------------------------
        !% Description
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: istep
            integer, pointer    :: Npack, thisP(:)
            integer             ::  mm
        !%------------------------------------------------------------------
        !% Preliminaries
            Npack => npack_elemP(ep_JM)   
            if (Npack < 1) return
        !%------------------------------------------------------------------
        !% Aliases
            thisP => elemP(1:Npack,ep_JM)
        !%------------------------------------------------------------------

        !% --- Consistency, store face values identical
        !% --- store face flowrate in JB for upstream (1) and downstream (2)
        !%     This is required because the face flowrates may have changed by 
        !%     interpolation on JB/CC branches
        call face_pull_facedata_to_JBelem (ep_JM, fr_Flowrate, er_Flowrate, .true.)
        call face_pull_facedata_to_JBelem (ep_JM, fr_Flowrate, er_Flowrate, .false.)

        ! do mm=1,Npack
        !     call lljunction_branch_getface (elemR(:,er_Flowrate),fr_Flowrate,thisP(mm),ei_Mface_uL,1)
        !     call lljunction_branch_getface (elemR(:,er_Flowrate),fr_Flowrate,thisP(mm),ei_Mface_dL,2)
        ! end do

        ! call util_utest_CLprint('bbb after branch getface -------------------------------------')

        !% --- compute the new junction element volume and head, JB flowrates
        !%     Does not change JB face values or JB values other than flowrate.
        call junction_calculation (thisP, Npack, istep)

        ! call util_utest_CLprint('ccc after junction calc -------------------------------------')

        call geo_assign_JB_from_head (ep_JM)

        !%  -- we need JB slot computations here
        call slot_JB_computation (ep_JM)

        ! call util_utest_CLprint('fff after slot JB -------------------------------------')

        !% --- force the JB element values to the faces for upstream (true)
        !%     and downstream (false) branches.
        !%     Forces elem flowrate, deltaQ
        call face_push_JBelem_to_face (ep_JM, fr_Flowrate, er_Flowrate, .true.)
        call face_push_JBelem_to_face (ep_JM, fr_Flowrate, er_Flowrate, .false.)
        call face_push_JBelem_to_face (ep_JM, fr_DeltaQ,   er_DeltaQ,   .true.)
        call face_push_JBelem_to_face (ep_JM, fr_DeltaQ,   er_DeltaQ,   .false.)

        !% --- make JB face head the same as JB element
        !%     note we do this to the fr_Head_u and then make both _u and _d 
        !%     identical
        call face_push_JBelem_to_face (ep_JM, fr_Head_u, er_Head, .false.)
        call face_push_JBelem_to_face (ep_JM, fr_Head_u, er_Head, .true.) 
        call face_make_up_dn_identical(fp_JB_IorS, fr_Head_u, fr_Head_d)

        call face_push_JBelem_to_face (ep_JM, fr_Area_u, er_Area, .false.)
        call face_push_JBelem_to_face (ep_JM, fr_Area_u, er_Area, .true.) 
        call face_make_up_dn_identical(fp_JB_IorS, fr_Area_u, fr_Area_d)

        !% --- auxiliary variables update
        !% --- wave speed, Froude number on JM
        Npack => npack_elemP(ep_JM)
        if (Npack > 0) then
            thisP => elemP(1:Npack, ep_JM)
            !% --- adjust JM for small or zero depth (may be redundant)
            call adjust_element_toplevel (JM)
            call update_wavespeed_element(thisP)
            call update_Froude_number_element (thisP) 
        end if

        ! call util_utest_CLprint('ggg after updates JM -------------------------------------')

        !% NOTE TRUE FORCES Q weight on JB to minimum, which
        !% means that face interpolation will have JB values
        !% dominate over adjacent CC values, but will be
        !% simple averaging with adjacen Diag Q values
        Npack => npack_elemP(ep_JB)
        if (Npack > 0) then 
            thisP => elemP(1:Npack, ep_JB)
            !% --- adjust JB for small or zero depth (may be redundant)
            call adjust_element_toplevel (JB)
            call update_Froude_number_element (thisP) 
            call update_wavespeed_element(thisP)
            call update_interpweights_JB (thisP, Npack, .false.)    !% 20230904brh change to FALSE
        end if

        ! call util_utest_CLprint('hhh after updates JB -------------------------------------')

    end subroutine junction_toplevel
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine junction_calculation (thisJM, Npack, istep)
        !%-----------------------------------------------------------------
        !% Description:
        !% Solves for the junction head and flows 
        !% Changes elemR(:,er_Volume), elemR(:,er_Flowrate), elem(:,er_Head)
        !% for Junction main and some storage rate data stored 
        !% for use in this routine
        !% ONLY CHANGES JB ELEMENT FLOWRATE DATA -- NOT HEAD
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: Npack, istep
            integer, dimension(Npack), intent(in) :: thisJM 

            integer :: mm, JMidx

            real(8), pointer :: Qstorage(:), Qoverflow(:), Qlateral(:)
            real(8), pointer :: MinHeadForOverflowPonding

            real(8) :: Qnet
            real(8) :: dQdHoverflow, dQdHstorage
            real(8) :: dH, resid, dHsave, dt1
            real(8), pointer :: dt !%, crk(:)

            real(8), dimension(2) :: Hbound

            logical :: isCrossingIntoSurcharge, isCrossingOutofSurcharge
            logical :: isCrossingIntoOverflowOrPonding, isCrossingOutofOverflowOrPonding
            logical :: isOverflow, isPonding, canoverflowOrPond
            

            real(8), parameter :: localEpsilon = 1.0d-6
        !%-----------------------------------------------------------------
        !% Aliases
            Qstorage    => elemSR(:,esr_JM_StorageRate) !% positive is increasing storage
            !% --- note that Qoverflow includes ponding rate
            Qoverflow   => elemSR(:,esr_JM_OverflowPondingRate) !% negative is outflow
            Qlateral    => elemR (:,er_FlowrateLateral) !% negative is outflow) 
            dt => setting%Time%Hydraulics%Dt 
        !%----------------------------------------------------------------- 

        do mm=1,Npack
            JMidx = thisJM(mm)
            ! isCrossingIntoSurcharge          = .false.
            ! isCrossingOutofSurcharge         = .false.
            ! isCrossingIntoOverflowOrPonding  = .false.
            ! isCrossingOutofOverflowOrPonding = .false.
            isOverflow                       = .false.
            isPonding                        = .false.
            MinHeadForOverflowPonding => elemSR(JMidx,esr_JM_MinHeadForOverflowPonding)
            
            if (elemSI(JMidx,esi_JM_OverflowType) == NoOverflow) then 
                canOverflowOrPond = .false.
            else 
                canOverflowOrPond = .true.
            end if

            ! if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then 
            !         print *, 'AAAA JB flowrate ',elemR(printJB,er_Flowrate)
            !         print *, 'plan area        ',elemSR(JMidx,esr_JM_Present_PlanArea)
            !         print *, 'slot width       ',elemR(JMidx,er_SlotWidth)
            !         print *, 'length           ',elemR(JMidx,er_Length)
            ! end if

            ! !% --- 20230913 -- presently sets large values
            ! call lljunction_main_head_bounds (JMidx, Hbound)
            
            ! !% --- convert Hbound to a deltaH bound
            ! Hbound = Hbound - elemR(JMidx,er_Head)

            !% --- set the present plan area
            elemSR(JMidx,esr_JM_Present_PlanArea) = lljunction_main_plan_area(JMidx)

            ! if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then 
            !     print *, 'A2222 JB flowrate ',elemR(printJB,er_Flowrate)
            !     print *, 'plan area        ',elemSR(JMidx,esr_JM_Present_PlanArea)
            ! end if

            !% --- set the overflow/ponding heads
            call lljunction_main_overflow_conditions (JMidx)

            ! if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then 
            ! !     print *, 'overflow depth: ',elemSR(JMidx,esr_JM_OverflowDepth)
            !     print *, 'BBBB'
            !         print *, 'dH           ',dH
            !         print *, 'Storage rate ',elemSR(JMidx,esr_JM_StorageRate)
            ! end if

            ! if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then 
            !     print *, 'BBBB JB flowrate ',elemR(printJB,er_Flowrate)
            ! end if
            !% --- get the net flowrate based on n data from all sources
            !%     including branches, overflow/ponding and lateral
            !%     if not overflow/ponding at PRESENT HEAD then logicals are returned as false
            !%     Overflow/ponding flow rates based on PRESENT Overflow/pond head diff,
            !%     which are zero if no
            call lljunction_main_netFlowrate &
                 (JMidx, Qnet, canOverflowOrPond, isOverflow, isPonding)

            !% --- fix flowrates if drying junction
            if ((elemR(JMidx,er_Volume_N0) + Qnet*dt) < setting%ZeroValue%Volume) then
                call lljunction_main_dryingfix (JMidx, Qnet)     
            end if


                !  if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then 
                !     print *, ''
                !     print *, 'in junction calculation after net flowrate'
                !     print *, 'net flowrate vol', Qnet  * setting%Time%Hydraulics%Dt
                !     print *, 'volumes ',elemR(JMidx,er_Volume), elemR(JMidx,er_Volume_N0)
                !     print *, 'Q vol up',elemR(JMidx+1,er_Flowrate) * setting%Time%Hydraulics%Dt
                !     print *, 'Q vol dn',elemR(JMidx+2,er_Flowrate) * setting%Time%Hydraulics%Dt
                !     print *, 'Qlat    ',elemR(JMidx,er_FlowrateLateral)  * setting%Time%Hydraulics%Dt
                !     print *, 'sum flow', &
                !        elemR(JMidx+1,er_Flowrate) * setting%Time%Hydraulics%Dt &
                !      - elemR(JMidx+2,er_Flowrate) * setting%Time%Hydraulics%Dt &
                !      + elemR(JMidx,er_FlowrateLateral) * setting%Time%Hydraulics%Dt
                !     print *, ' '
                ! end if

            ! if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then 
            !      print *, 'Overflow rate:  ',elemSR(JMidx,esr_JM_OverflowPondingRate)
            ! end if

                ! if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then 
                !     print *, 'CCCC'
                !     print *, 'dH           ',dH
                !     print *, 'Storage rate ',elemSR(JMidx,esr_JM_StorageRate)
                ! end if
                !  if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then 
                !     print *, 'CCCC JB flowrate ',elemR(printJB,er_Flowrate)
                ! end if     

            !% --- Compute dH
            call lljunction_main_dHcompute &
                (JMidx, dH, dQdHoverflow,  dQdHstorage, Qnet, Hbound, istep, &
                isOverflow, isPonding, .false., .false.)

                ! if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then 
                !     print *, 'DDDD JB flowrate ',elemR(printJB,er_Flowrate)
                ! end if

            !% --- Limit head increase for cases of in/out of surcharge or overflow
            !%     This ensures a 2-step process to cross the boundary.
            !% --- check if crossing in or out of overflow/ponding -- only applies
            !%     if NOT crossing in/out of surcharge
        ! dHsave = dH
        ! if (canOverflowOrPond) then 
        !     call lljunction_main_iscrossing_overflow_or_ponding  &
        !         ( JMidx,  dH, isCrossingIntoOverflowOrPonding, isCrossingOutofOverflowOrPonding)
        !     if (isCrossingIntoOverflowOrPonding .or. isCrossingOutofOverflowOrPonding) then 
        !         !% --- limit this dH to the distance from the present head to MinHeadForOverflowPonding
        !         !%     the second step is computed later
        !         dH = MinHeadForOverflowPonding - elemR(JMidx,er_Head) 
        !     end if       
        ! else
        !     !% --- no action, accept defaults
        ! end if
        
        ! if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then 
        !     print *, 'dhsave, dh 1: ',dHsave, dH
        ! end if

            !% --- check if crossing in or out of surcharge -- only applies if NOT overflow crossing
            !%     note that JM with zero for OverflowHeightAboveCrown are given false for
            !%     eYN_canSurcharge, so they will not use this algorithm
        ! if (elemYN(JMidx,eYN_canSurcharge)) then
        !     if ((.not. isCrossingIntoOverflowOrPonding ) .and. (.not. isCrossingOutofOverflowOrPonding)) then
        !         call lljunction_main_iscrossing_surcharge &
        !             ( JMidx,  dH, isCrossingIntoSurcharge, isCrossingOutofSurcharge)

        !             ! if (printJM == JMidx) then 
        !             !     print *, ' '
        !             !     print *, 'crossing ',isCrossingIntoSurcharge, isCrossingOutofSurcharge
        !             !     print *, ' '
        !             ! end if

        !         if (isCrossingIntoSurcharge .or. isCrossingOutofSurcharge) then 
        !             !% --- limit this dH to the distance from the present head to the crown
        !             dH = elemR(JMidx,er_Zcrown) - elemR(JMidx,er_Head) 
        !         end if
        !     else 
        !         !% --- no action
        !     end if
        ! else
        !     !% --- no action, accept defaults
        ! end if

            ! if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then 
            !     print *, 'DDDD'
            !     print *, 'dH           ',dH
            !     print *, 'Storage rate ',elemSR(JMidx,esr_JM_StorageRate)
            ! end if

            ! if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then 
            !     print *, 'Crosing in, Out ',isCrossingIntoSurcharge, ' ',isCrossingOutofSurcharge
            !     !print *, 'dhsave, dh 2: ',dHsave, dH
            ! end if

            ! if (printJM == JMidx) then 
            !     print *, ' '
            !     print *, 'dH, dHsave ',dH, dHsave
            !     print *, ' '
            ! end if

            !% --- update values (head, depth, volume, Qoverflow, DeltaQ, branch Q, QnetBranches
            !%     Qstorage) 
            !%     Note that if not a Crossing... then the dH is used for head, depth, volume.
            !%     When any Crossing... is true then head, depth, volume values are set to the Full values
            call lljunction_main_update_intermediate &
                (JMidx, istep, dH, dQdHoverflow, dQdHstorage, MinHeadForOverflowPonding, isOverflow, isPonding, &
                 isCrossingIntoOverflowOrPonding, isCrossingOutofOverflowOrPonding, &
                 isCrossingIntoSurcharge, isCrossingOutofSurcharge)


                !  if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then 
                !     print *, 'EEEE JB flowrate ',elemR(printJB,er_Flowrate)
                ! end if

                ! if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then 
                !     print *, 'EEEE'
                !     print *, 'dH           ',dH
                !     print *, 'Storage rate ',elemSR(JMidx,esr_JM_StorageRate)
                ! end if

                ! if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then 
                !     print *, 'Overflow rate:  ',elemSR(JMidx,esr_JM_OverflowPondingRate)
                ! end if
            !% --- second part: additional dH when crossing surcharge
            ! if ((isCrossingIntoSurcharge)         .or. (isCrossingOutofSurcharge) .or.            &
            !     (isCrossingIntoOverflowOrPonding) .or. (isCrossingOutofOverflowOrPonding)) then 

            !     !% --- Additional dH after crossing surcharge/free threshold
            !     dH = dHsave - dH
            !     call lljunction_main_dHcompute &
            !             (JMidx, dH, dQdHoverflow, dQdHstorage, Qnet, Hbound, istep, &
            !              isOverflow, isPonding, isCrossingIntoSurcharge, isCrossingOutofSurcharge)

            !     !% --- update values (head, depth, volume, deltaQ) based on dH
            !     !%     Use isCrossing... = .false. so that full values are
            !     !%     note used
            !     call lljunction_main_update_intermediate &
            !         (JMidx, istep, dH, dQdHoverflow, dQdHstorage, MinHeadForOverflow, isOverflow, isPonding, &
            !         .false., .false., .false., .false.)

            ! else
            !     !% --- no action
            ! end if
            
            ! if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then 
            !     print *, 'Volume Overflow CCC:  ',elemR(JMidx,er_VolumeOverFlow)
            ! end if

            ! call lljunction_main_update_final (JMidx, istep)

            ! if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then 
            !     print *, 'Volume Overflow DDD:  ',elemR(JMidx,er_VolumeOverFlow)
            ! end if

            !% -- THIS CAN BE COMMENTED AS THE RESIDUAL IS COMPUTED IN 
            !%    junction_mass_conservation
            ! resid = lljunction_conservation_residual (JMidx) 

            ! if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then 
            !     print *, ''
            !     print *, 'in junction calculation'
            !     print *, 'volumes ',elemR(JMidx,er_Volume), elemR(JMidx,er_Volume_N0)
            !     print *, 'Q vol up',elemR(JMidx+1,er_Flowrate) * setting%Time%Hydraulics%Dt
            !     print *, 'Q vol dn',elemR(JMidx+2,er_Flowrate) * setting%Time%Hydraulics%Dt
            !     print *, 'Qlat    ',elemR(JMidx,er_FlowrateLateral)  * setting%Time%Hydraulics%Dt
            !     print *, 'sum flow', &
            !        elemR(JMidx+1,er_Flowrate) * setting%Time%Hydraulics%Dt &
            !      - elemR(JMidx+2,er_Flowrate) * setting%Time%Hydraulics%Dt &
            !      + elemR(JMidx,er_FlowrateLateral) * setting%Time%Hydraulics%Dt
            !     print *, ' '
            ! end if

        end do

    end subroutine junction_calculation
!%
!%==========================================================================
!%==========================================================================
!% 
    subroutine junction_mass_conservation (thisColP, istep)
        !%-----------------------------------------------------------------
        !% Description:
        !% Finds residual of junctions and adjusts JB and face values
        !% to obtain mass conservation
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisColP, istep
            integer, pointer    :: thisJM(:), Npack, JMidx
            real(8), pointer    :: Qoverflow(:), Qstorage(:), Qlateral(:)

            integer :: mm
            real(8) :: resid !% residual in/out flowrate
            real(8) :: QnetIn, QnetOut, QnetBranches

            real(8), parameter :: local_epsilon = 1.0d-15
        !%-----------------------------------------------------------------
        !% Aliases
            Npack   => npack_elemP(thisColP)
            if (Npack < 1) return
            thisJM  => elemP(1:Npack,thisColP)
            Qoverflow => elemSR(:,esr_JM_OverflowPondingRate)
            Qstorage  => elemSR(:,esr_JM_StorageRate)
            Qlateral  => elemR (:,er_FlowrateLateral)
        !%-----------------------------------------------------------------

        do mm=1,Npack
            JMidx => thisJM(mm)

            ! !% --- compute junction residual
            resid = lljunction_conservation_residual (thisJM(mm)) 

            ! if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then 
            !     print *, 'Cons Resid before:', resid
            !     print *, 'vol - full volume:', elemR(printJM,er_Volume) - elemR(printJM,er_FullVolume)
            ! end if

            if (abs(resid) > local_epsilon) then 
                !% --- note that QnetIn > 0 and QnetOut < 0
                QnetIn  = lljunction_branch_Qnet (thisJM(mm),+oneI)
                QnetOut = lljunction_branch_Qnet (thisJM(mm),-oneI)
                if (Qoverflow(thisJM(mm)) > zeroR) then
                    QnetIn = QnetIn + Qoverflow(thisJM(mm))
                else
                    QnetOut = QnetOut + Qoverflow(thisJM(mm))
                end if

                ! if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then 
                !     print *, 'Qnet before', QnetIn, QnetOut, Qoverflow(printJM)
                ! end if

                !% --- ad hoc adjustment of elemR(:,er_Flowrate) of JB and Qoverflow of JM
                call lljunction_conservation_fix(thisJM(mm),resid, QnetIn, QnetOut)

                ! if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then 
                !     print *, 'Cons Resid after:', resid
                !     print *, 'vol - full volume:', elemR(printJM,er_Volume) - elemR(printJM,er_FullVolume)
                ! end if

                ! if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then 
                !     print *, 'Qnet after', QnetIn, QnetOut, Qoverflow(printJM)
                ! end if

                !% --- note that QnetIn, QnetOut, and elemSR(:,esr_JM_OverflowPondingRate) 
                !%     are all changed in lljunction_conservation_fix
               
                !% --- update net Q branches for residual changes
                QnetBranches = lljunction_main_sumBranches (thisJM(mm),er_Flowrate,elemR)

                !% --- volume conservative storage rate (not based on dH)
                Qstorage(thisJM(mm)) = QnetBranches + Qoverflow(thisJM(mm)) + Qlateral(thisJM(mm))

                !% --- update Volume, VolumeOverflow and JB face values
                call lljunction_main_update_Qdependent_values (thisJM(mm), istep)

                ! if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then 
                !     print *, 'Volume Oveflow', elemR(printJM,er_VolumeOverFlow)
                ! end if

                ! if ((setting%Time%Step > stepCut) .and. (JMidx == printJM)) then 
                !     print *, 'at end of junction mass conservation'
                !     print *, 'vol - full volume:', elemR(printJM,er_Volume) - elemR(printJM,er_FullVolume)
                ! end if

            end if 
        end do

    end subroutine junction_mass_conservation
!%    
!%========================================================================== 
!% END OF MODULE
!%==========================================================================
end module junction_elements