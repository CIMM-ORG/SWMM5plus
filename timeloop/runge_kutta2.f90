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
    public :: rk2_toplevel_ETM_2
    public :: rk2_toplevel_ETM_3
    public :: rk2_toplevel_ETM_4
    public :: rk2_toplevel_ETM_5
    public :: rk2_toplevel_ETM_6
    public :: rk2_toplevel_ETM_7
    public :: rk2_toplevel_AC
    public :: rk2_toplevel_ETMAC

    contains
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine rk2_toplevel_ETM_7 ()
        !%------------------------------------------------------------------
        !% Description:
        !% single RK2 step for explicit time advance of SVE
        !%------------------------------------------------------------------
        !% Declarations:
            integer :: istep, ii, kk, mm
            integer, pointer :: Npack, thisP(:), fup(:), fdn(:), thisJB(:)
            integer, pointer :: JMar(:), fAdj(:)
            
            integer :: fr_hAdj, fr_aAdj

            real(8), pointer :: deltaH(:), grav, dt
            real(8)          :: bsign
            real(8)          :: volume1, volume2, inflowVolume, outflowVolume

            logical :: isUpstreamFace, isUpstreamBranch

            integer :: epCCcol

            integer :: printJB = 75

            character(64) :: subroutine_name = 'rk2_toplevel_ETM_7'
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
            
            !print *,' '
            !print *, ' '
          ! call util_utest_CLprint ('======= AAA  start of RK2 ==============================')    
 
        !%==================================    
        !% --- Initial adjustments
        istep = zeroI

        !% --- update Diagnostic elements
        if (N_diag > 0) then 
            !% STEP A
            !% --- update flowrates for aa diagnostic elements
            call diagnostic_by_type (ep_Diag, istep)  

                ! ! call util_utest_CLprint ('------- BBB  after diagnostic')

            !% STEP B
            !% --- push the diagnostic flowrate data to faces -- true is upstream, false is downstream
            call face_push_elemdata_to_face (ep_Diag, fr_Flowrate, er_Flowrate, elemR, .true.)
            call face_push_elemdata_to_face (ep_Diag, fr_Flowrate, er_Flowrate, elemR, .false.)

                ! ! call util_utest_CLprint ('------- CCC  after face_push_elemdata_to_face')
        end if

        !% --- push inflows on CC upstream of JB to face
        if (N_nJM > 0) then 
            Npack => npack_elemP(ep_CC_UpstreamOfJunction)
            if (Npack > 0) then 
                thisP => elemP(1:Npack,ep_CC_UpstreamOfJunction)
                fAdj  => elemI(:,ei_Mface_dL)
                where (elemR(thisP,er_Flowrate) > zeroR)
                    faceR(fAdj(thisP),fr_Flowrate)   = elemR(thisP,er_Flowrate)
                    faceR(fAdj(thisP),fr_Velocity_u) = elemR(thisP,er_Flowrate) / faceR(fAdj(thisP),fr_Area_u)
                    faceR(fAdj(thisP),fr_Velocity_d) = elemR(thisP,er_Flowrate) / faceR(fAdj(thisP),fr_Area_d)
                    !% check if any shared face has been diverged
                    where (faceYN(fAdj(thisP),fYN_isSharedFace))
                        faceYN(fAdj(thisP),fYN_isSharedFaceDiverged) = .true.
                    end where
                endwhere
            end if
            Npack => npack_elemP(ep_CC_DownstreamOfJunction)
            if (Npack > 0) then 
                thisP => elemP(1:Npack,ep_CC_DownstreamOfJunction)
                fAdj  => elemI(:,ei_Mface_uL)
                where (elemR(thisP,er_Flowrate) < zeroR)
                    faceR(fAdj(thisP),fr_Flowrate) = elemR(thisP,er_Flowrate)
                    faceR(fAdj(thisP),fr_Velocity_u) = elemR(thisP,er_Flowrate) / faceR(fAdj(thisP),fr_Area_u)
                    faceR(fAdj(thisP),fr_Velocity_d) = elemR(thisP,er_Flowrate) / faceR(fAdj(thisP),fr_Area_d)
                    !% check if any shared face has been diverged
                    where (faceYN(fAdj(thisP),fYN_isSharedFace))
                        faceYN(fAdj(thisP),fYN_isSharedFaceDiverged) = .true.
                    end where
                endwhere
            end if
        end if
            ! call util_utest_CLprint ('------- CCC.1  after push inflows to JB faces')

        !% --- get the velocity in the junction main
        if (N_nJM > 0) then
            !% STEP E
            !% --- provide junction outlet JB velocity estimate
            call junction_main_velocity (ep_JM)

                ! call util_utest_CLprint ('------- DDD  after junction_main_velocity')
        end if

        !% --- store junction-adjacent element data on face
        !%     QUESTION -- SHOULD THIS BE IN THE RK ITERATION FOR UPDATES?
        !%     ANSWER: NO, as long as the second step junction solution is NOT the backwards euler.
        !%     TO BE MOVED TO junction_branch_adjacent 
        if (N_nJM > 0) then 
            do ii=1,2
                !% --- cycle over upstream and downstream faces of junction
                if (ii==1) then 
                    !% --- an upstream face is CC downstream of junction
                    isUpstreamFace = .true.
                    epCCcol = ep_CC_DownstreamOfJunction
                else
                    !% -- a downstream face is CC upstream of junction
                    isUpstreamFace = .false.
                    epCCcol = ep_CC_UpstreamOfJunction
                end if
                call face_push_elemdata_to_face (epCCcol, fr_Head_Adjacent,     er_Head,         elemR, isUpstreamface)
                call face_push_elemdata_to_face (epCCcol, fr_Topwidth_Adjacent, er_Topwidth,     elemR, isUpstreamface)
                call face_push_elemdata_to_face (epCCcol, fr_Length_Adjacent,   er_Length,       elemR, isUpstreamface)
                call face_push_elemdata_to_face (epCCcol, fr_Zcrest_Adjacent,   er_Zbottom,      elemR, isUpstreamface)
                call face_push_elemdata_to_face (epCCcol, fr_Velocity_Adjacent, er_Velocity,     elemR, isUpstreamface)
                call face_push_elemdata_to_face (epCCcol, fr_Froude_Adjacent,   er_FroudeNumber, elemR, isUpstreamface)
                call face_push_elemdata_to_face (epCCcol, fr_Depth_Adjacent,    er_Depth,        elemR, isUpstreamface)
            end do

                ! call util_utest_CLprint ('------- EEE  after face push of adjacent data')
        end if


        !% --- compute flows on JB elements from energy equation
        !% TO BE MOVED TO junction_branch_element_flowrates
        if (N_nJM > 0) then 
            !% --- array for the JM index
            JMar   => elemSI(:,esi_JunctionBranch_Main_Index)           
            !% --- cycle through nominal upstream and downstream JB
            do ii=1,2
                !% --- get the upstream or downstream JB elements
                if (ii==1) then 
                    !% --- upstream JB
                    isUpstreamBranch = .true.
                    bsign = oneR
                    Npack => npack_elemP(ep_JB_Upstream)
                    if (Npack > 0) then 
                        thisJB => elemP(1:Npack,ep_JB_Upstream)
                    end if
                    !% --- the JB-adjacent face is upstream
                    fAdj => elemI(:,ei_Mface_uL)
                    fr_hAdj = fr_Head_u !% -- use u for jump purposes? QUESTION
                    fr_aAdj = fr_Area_u !% QUESTION
                else
                    !% --- downstream JB
                    bsign = -oneR
                    isUpstreamBranch = .false.
                    Npack => npack_elemP(ep_JB_Downstream)
                    if (Npack > 0) then 
                        thisJB => elemP(1:Npack,ep_JB_Downstream)
                    end if
                    !% --- the JB-adjacent face is downstream
                    fAdj => elemI(:,ei_Mface_dL)
                    fr_hAdj = fr_Head_d !% --- use d for jump purposes? QUESTION
                    fr_aAdj = fr_Area_d !% QUESTION
                end if
                !% --- head difference from junction main to face
                deltaH => elemR(:,er_Temp01)
                deltaH = zeroR
                where (elemR(thisJB,er_Depth) > setting%ZeroValue%Depth)
                    deltaH(thisJB) =  elemR(JMar(thisJB),er_Head) - faceR(fAdj(thisJB),fr_hAdj)
                endwhere

                ! print *, ' '
                ! print *, 'thisJB ', thisJB
                ! print *, 'depth  ', elemR(thisJB,er_Depth)
                ! print *, 'flowr  ', elemR(thisJB,er_Flowrate)
                ! print *, 'deltaH ', deltaH(thisJB)
                ! print *, 'flo*b  ', elemR(thisJB,er_Flowrate) * bsign
                ! print *,' '


                !% --- Subcritical outflow junction (super overwrites below)
                where ((deltaH(thisJB) > 1.0d-8) .and. (elemR(thisJB,er_Flowrate) * bsign < -1.0d-8))
                        !% --- outflow from JB upstream 
                        !% --- or outflow from JB downstream
                        !%     NOTE: without junction main approach velocity
                        elemR(thisJB,er_Flowrate) = - bsign * faceR(fAdj(thisJB),fr_aAdj)                  &
                                * sqrt(                                                                     &
                                        (twoR * grav * deltaH(thisJB))                                      &
                                        /(oneR + elemSR(thisJB,esr_JunctionBranch_Kfactor))                 &
                                        )   
                endwhere 

                where ((deltaH(thisJB) > 1.0d-8) .and. (elemR(thisJB,er_Flowrate) * bsign > -1.0d-8))
                    !% --- inconsistent flow direction and pressure gradient   
                    elemR(thisJB,er_Flowrate) = zeroR
                endwhere 

                where  ((deltaH(thisJB) .le. 1.0d-8) .and. (deltaH(thisJB) .ge. -1.0d-8))
                    !% --- effectively zero pressure gradient
                    elemR(thisJB,er_Flowrate) = zeroR
                endwhere 

                where ((deltaH(thisJB) < -1.0d-8) .and. (elemR(thisJB,er_Flowrate) * bsign > -1.0d-8))
                    !% --- inflow using face flowrate
                    elemR(thisJB,er_Flowrate) =  faceR(fAdj(thisJB),fr_Flowrate) 
                endwhere

                where  ((deltaH(thisJB) < -1.0d-8) .and. (elemR(thisJB,er_Flowrate) * bsign .le. -1.0d-8))
                    !% --- inconsistent flow direction and pressure gradient   
                    elemR(thisJB,er_Flowrate) = zeroR
                endwhere
 

                ! do mm = 1,size(thisJB)
                !     if (thisJB(mm) == printJB) then 
                !         print *, JMar(thisJB(mm))
                !         print *, 'AAAA head     ',elemR(JMar(thisJB(mm)),er_Head), faceR(fAdj(thisJB(mm)),fr_hAdj)
                !         print *, 'AAAA dH       ',deltaH(thisJB(mm))
                !         print *, 'AAAA FLowrate ',thisJB(mm), elemR(thisJB(mm),er_Flowrate)
                !     end if
                ! end do
                ! print *, ' '
                ! print *, 'before ', elemR(7,er_Flowrate), elemR(5,er_FroudeNumber)

                !% THE SUPERCRITICAL CAUSED PROBLEMS FOR THIN LAYERS
                ! !% --- Supercritical junction (requires JM velocity and froude number)
                ! if (isUpstreamBranch) then
                !     where (elemR(JMar(thisJB),er_FroudeNumber) < -oneR)
                !         elemR(thisJB,er_Flowrate) = elemR(JMar(thisJB),er_Velocity) * faceR(fAdj(thisJB),fr_aAdj)
                !         !elemR(thisJB,er_Flowrate) = elemR(thisJB,er_Velocity) * faceR(fAdj(thisJB),fr_aAdj)
                !     endwhere     
                ! else
                !     where (elemR(JMar(thisJB),er_FroudeNumber) > oneR)
                !         elemR(thisJB,er_Flowrate) = elemR(JMar(thisJB),er_Velocity) * faceR(fAdj(thisJB),fr_aAdj)
                !         !elemR(thisJB,er_Flowrate) = elemR(thisJB,er_Velocity) * faceR(fAdj(thisJB),fr_aAdj)
                !     endwhere     
                ! end if  

                ! print *, 'after ',elemR(7,er_Flowrate)
                ! print *, ' '

                !% THIS CAUSED HUGE PROBLEMS
                ! ! !% --- volume limiter on outflow
                ! if (isUpstreamBranch) then 
                ! !     ! !% --- volume limit on outflow from upstream branch
                ! !     ! where ( elemR(JMar(thisJB),er_Volume) < -dt * elemR(thisJB,er_Flowrate) )
                ! !     !     elemR(thisJB,er_Flowrate) = - elemR(JMar(thisJB),er_Volume) / (dt * twoR)
                ! !     ! endwhere
                ! else 
                ! ! 20230430brh
                !     !% --- volume limit on outflow on downstream branch
                !     !where ( elemR(JMar(thisJB),er_Volume) < dt * elemR(thisJB,er_Flowrate) )
                !     !    elemR(thisJB,er_Flowrate) = elemR(JMar(thisJB),er_Volume) / (dt * twoR)
                !     do mm=1,size(thisJB)
                !         print *, thisJB(mm), elemR(thisJB(mm),er_Flowrate), elemR(JMar(thisJB(mm)),er_Volume) / (dt * twoR)
                !     end do
                !     !endwhere
                ! end if

                ! if (isUpstreamBranch) then 
                ! else 
                !     !do mm=1,size(thisJB)
                !         print *, 'vel to empty', 127, 125, elemR(125,er_Volume) / dt
                !     !end do
                ! end if

                !% --- update velocity
                where (elemR(thisJB,er_Area) > setting%ZeroValue%Area)
                    elemR(thisJB,er_Velocity) = elemR(thisJB,er_Flowrate) / elemR(thisJB,er_Area)
                elsewhere 
                    elemR(thisJB,er_Velocity) = zeroR
                endwhere

                !% --- apply velocity limiter
                where (abs(elemR(thisJB,er_Velocity)) > setting%Limiter%Velocity%Maximum)
                    elemR(thisJB,er_Velocity) = sign(setting%Limiter%Velocity%Maximum * 0.99d0, elemR(thisJB,er_Velocity))
                endwhere

                ! !% --- zero outflow from upstream branch if JB depth is zero
                ! where ((elemR(thisJB,er_Depth) .le. setting%ZeroValue%Depth)    .and. &
                !        (elemSI(thisJB,esi_JunctionBranch_IsUpstream) .eq. oneI) .and. &
                !        (elemR(thisJB,er_Flowrate) < zeroR) )

                !        elemR(thisJB,er_Flowrate) = zeroR 
                !        elemR(thisJB,er_Velocity) = zeroR

                ! endwhere

                ! !% --- zero outflow from downstream branch if JB depth is zero
                ! where ((elemR(thisJB,er_Depth) .le. setting%ZeroValue%Depth)    .and. &
                !        (elemSI(thisJB,esi_JunctionBranch_IsUpstream) .ne. oneI) .and. &
                !        (elemR(thisJB,er_Flowrate) > zeroR) )

                !        elemR(thisJB,er_Flowrate) = zeroR 
                !        elemR(thisJB,er_Velocity) = zeroR

                ! endwhere

                ! print *, ' '
                ! print *, 'in energy term ',elemR(7,er_Flowrate)
                ! print *, elemR(5,er_Volume) / elemR(7,er_Flowrate), setting%Time%Hydraulics%Dt
                ! print *, ' '

            end do

            Npack => npack_elemP(ep_JB)
            if (Npack > 0) then 
                thisP => elemP(1:Npack, ep_JB)
                call update_Froude_number_element (thisP) 
                call update_wavespeed_element (thisP)
                call update_interpweights_JB (thisP, Npack, .true.)
            end if

                ! call util_utest_CLprint ('------- FFF  after Energy-based JB flows')
        end if
    

        !% --- Push JB element flows to faces
        !%     Note that for inflows this doesn't change anything because the JB and faces
        !%     are still the same. This only changes tht outflows
        if (N_nJM > 0) then 
            !% --- push JB upstream branch Q to the upstream face
            call face_push_elemdata_to_face (ep_JB_Upstream,fr_Flowrate,   er_Flowrate, elemR,.true.)
            !% --- push JB downstream branch Q to the downstream face
            call face_push_elemdata_to_face (ep_JB_Downstream,fr_Flowrate, er_Flowrate, elemR,.false.)

                ! call util_utest_CLprint ('------- GGG after push JB flowrates to face')
        end if

        !% ==============================================================
        !% --- face sync (saz05022023)
        !%     sync all the images first. then copy over the data between
        !%     shared-identical faces. then sync all images again
        sync all

        call face_shared_face_sync (fp_noBC_IorS)

        sync all
        !% 
        !% ==============================================================
        
        !% --- store the junction dQdH used in Backwards Euler
        if (N_nJM > 0) then  
            
            Npack => npack_elemP(ep_JB)
            if (Npack > 0) then     
                thisP => elemP(1:Npack,ep_JB)

                !% --- dQdH for both CC-adjacent and Diag-adjacent JB 
                call junction_branch_dQdH (thisP, Npack, 1)  

                    ! call util_utest_CLprint ('------- HHH  after junction_branch_dQdH')
            end if
        end if

        !% --- ensure DeltaQ are zero 
        if (N_nJM > 0) then        
            faceR(:,fr_DeltaQ) = zeroR
            elemR(:,er_DeltaQ) = zeroR
        end if

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
            call update_auxiliary_variables_CC (                  &
                ep_CC, ep_CC_Open_Elements, ep_CC_Closed_Elements, &
                .true., .false., dummyIdx)
                
                ! call util_utest_CLprint ('------- JJJ  after update_aux...CC step')

            !% --- zero and small depth adjustment for elements
            call adjust_element_toplevel (CC)

                ! call util_utest_CLprint ('------- KKK  after adjust element CC')

            !% --- Conservative JM update in 2nd RK step
            !% MOVE THIS TO A JUNCTION SUBROUTINE
            if (N_nJM > 0) then 
                if (istep == 1) then
                    !% --- ensure JB interpweights for Q are forced for JB dominance: QUESTION -- SHOULD THIS BE TRUE OR NOT?
                    Npack => npack_elemP(ep_JB)
                    if (Npack > 0) then 
                        thisP => elemP(1:Npack, ep_JB)
                        call update_interpweights_JB (thisP, Npack, .true.)
                    end if

                else if (istep == 2)then 

                    !print *,'in junction '

                    Npack => npack_elemP(ep_JM)
                    if (Npack > 0) then
                        thisP => elemP(1:Npack,ep_JM)
                        !% --- new junction volume from conservative face fluxes
                        call junction_main_volume_advance (ep_JM, Npack)
                        !% --- new junction plan area
                        call geo_plan_area_from_volume_JM (elemPGetm, npack_elemPGetm, col_elemPGetm)
                        !print *, 'plan area ',elemSR(thisP(1),esr_Storage_Plan_Area)


                        !% --- new junction depth 
                        !% NOTE: THIS USES storage_implied_depth_from_volume
                        !% that limits depth based on fulldepth
                        call geo_depth_from_volume_JM (elemPGetm, npack_elemPGetm, col_elemPGetm)
                    
                        !% --- new JM head
                        elemR(thisP,er_Head) = llgeo_head_from_depth_pure (thisP,elemR(thisP,er_Depth))
                        elemR(thisP,er_EllDepth) = elemR(thisP,er_Depth)

                        !% NEED PREISSMANN SLOT STUFF HERE (OR BELOW?) TO GET CORRECT HEAD
                        !% QUESTION: SHOULD ELLDEPTH BE DEPTH OR HEAD - ZBOTTOM?

                        !% --- adjust JM for small or zero depth
                        call adjust_element_toplevel (JM)
                        !% --- assign JB values based on new JM head
                        call geo_assign_JB (ep_JM) !% HACK  revise using ep_JB
                        !% --- adjust JB for small or zero depth
                        call adjust_element_toplevel (JB)
                        

                        ! QUESTION? NEED THIS HERE? see  geometry_toplevel_JMJB
                        !% call slot_JB_computation (ep_JM)
            
                    end if

                    ! call util_utest_CLprint ('------- LLL  after JM 2nd step ')

                    !% --- auxiliary variables update
                    !%     replaces update_auxiliary_variables_JMJB
                    !% NOTE TRUE FORCES Q weight on JB to minimum, which
                    !% means that face interpolation will have JB values
                    !% dominate over adjacent CC values, but will be
                    !% simple averaging with adjacen Diag Q values -- HERE WE USE TRUE IS THIS CORRECT?
                    Npack => npack_elemP(ep_JB)
                    if (Npack > 0) then 
                        thisP => elemP(1:Npack, ep_JB)
                        call update_Froude_number_element (thisP) 
                        call update_wavespeed_element (thisP)
                        call update_interpweights_JB (thisP, Npack, .true.)
                    end if

                    !% --- wave speed, Froude number on JM
                    Npack => npack_elemP(ep_JM)
                    if (Npack > 0) then
                        thisP => elemP(1:Npack, ep_JM)
                        call update_wavespeed_element (thisP)
                        call update_Froude_number_element (thisP) 
                    end if
                
                        ! call util_utest_CLprint ('------- MMM  after update aux JMJB 2nd step ')

                    !% --- QUESTION -- IS THIS NEEDED HERE? 
                    !%     should this be outside of the IF/ENDIF for the istep=2?
                    if (N_diag > 0) then 
                            !% --- update flowrates for diagnostic elements that are not adjacent to JB
                            !call diagnostic_by_type (ep_Diag_notJBadjacent, istep)  
                            !20230423 test using all DIAG
                            call diagnostic_by_type (ep_Diag, istep) 

                            ! ! call util_utest_CLprint ('------- NNN  after update aux JMJB 2nd step ')
                    end if
                end if
            end if

            !% --- interpolate all data to faces
            sync all
            call face_interpolation(fp_noBC_IorS, .true., .true., .true., .false., .true.) 

                ! call util_utest_CLprint ('------- OOO  after face interp')

            !% --- reset the face zerodepth array 
            !%     note, this affects both interior shared faces, so this must be done
            !%     once before all the face_zeroDepth calls
            !%     HACK THIS ARRAY APPEARS TO BE UNUSED
           ! faceI(:,fi_zeroDepth) = zeroI  !% default to no zerodepth

             !% --- update various packs of zeroDepth faces 20230430brh
            call pack_CC_zeroDepth_interior_faces ()
            if (N_nJM > 0) then 
                call pack_JB_zeroDepth_interior_faces ()
            end if
            sync all
            call pack_CC_zeroDepth_shared_faces ()  !% HACK NEEDS REVIEW
            if (N_nJM > 0) then
                call pack_JB_zeroDepth_shared_faces ()  !% HACK STUB ROUTINE NOT COMPLETE
            end if

            ! print *, ' '
            ! print *, 'face pack JB'
            ! print *, faceP(1:npack_faceP(fp_JB_downstream_is_zero_IorS),fp_JB_downstream_is_zero_IorS)
            ! print *, ' '

            !% --- set face geometry and flowrates where adjacent element is zero
            !%     only applies to faces with CC on both sides

            ! print *, ' '
            ! print *, 'face with zero downstream'
            ! print *, faceP(1:npack_faceP(fp_CC_downstream_is_zero_IorS),fp_CC_downstream_is_zero_IorS)
            ! print *, ' '

            call face_zeroDepth (fp_CC_downstream_is_zero_IorS, &
                fp_CC_upstream_is_zero_IorS,fp_CC_bothsides_are_zero_IorS)

                ! call util_utest_CLprint ('------- PPP.01 after face zerodepth ')

            if (N_nJM > 0) then
                ! print *, ' '
                ! print *, faceP(1:npack_faceP(fp_JB_downstream_is_zero_IorS),fp_JB_downstream_is_zero_IorS)
                ! print *, ' '
                !% --- set face geometry and flowrates where adjacent element is zero
                !%     only applies to faces with JB on one side
                call face_zeroDepth (fp_JB_downstream_is_zero_IorS, &
                    fp_JB_upstream_is_zero_IorS,fp_JB_bothsides_are_zero_IorS)
            end if

                ! call util_utest_CLprint ('------- PPP.02 after face zerodepth ')

            !% --- enforce open (1) closed (0) "setting" value from EPA-SWMM
            !%     for all CC and Diag elements (not allowed on junctions)
            call face_flowrate_for_openclosed_elem (ep_CCDiag)

            !% HACK -- FOR SHARED WE NEED TO IDENTIFY ANY SHARED FACES THAT HAVE AN
            !% UPSTREAM elemR(:,er_Setting) = 0.0 and set their face Q and V to zero.

                ! call util_utest_CLprint ('------- QQQ after openclosed setting')

            !%=========================================
            !% BEGIN JUNCTION SOLUTION
            if (istep == 1) then 
                if (N_nJM > 0) then 

                    !% ---
                    !%     forces JB elem Q to faces (overriding the interpolation)
                    !%     computes new JM Volume, Head, JB fluxes, JB DeltaQ
                    !%     for conservation and dH change in head 
                    !%     assigns new JB and JM aux values
                    call junction_toplevel_5 (istep)

                        ! call util_utest_CLprint ('------- RRR  after junction toplevel')

                    !% --- force the JB element values to the faces for upstream (true)
                    !%     and downstream (false) branches.
                    !%     Forces elem  flowrate, deltaQ
                    !% HACK - INCLUDES DEPTH, ETC, BUT COMMENTED OUT FOR NOW
                    !% NOTE - THIS IS NOT in junction_toplevel because this
                    !% is the subroutine that requires face syncing afterwards
                    call face_force_JBadjacent_values (ep_JM, .true.)
                    call face_force_JBadjacent_values (ep_JM, .false.)

                        ! call util_utest_CLprint ('------- SSS  after face_force_JBadjacent')
                end if
                
                !sync all  !% cannot be in an if N_nJM statement
                !% ==============================================================
                !% --- face sync (saz05022023)
                !%     sync all the images first. then copy over the data between
                !%     shared-identical faces. then sync all images again
                sync all

                call face_shared_face_sync (fp_noBC_IorS)

                sync all
                !% 
                !% ==============================================================

                !%================================
                !% --- HACK need to force face data on all JB faces for image containing JB
                !%     element to the connected image. This is a direct transfer of face
                !%     Data without transferring ghost element data or interpolation
                !%================================

                if (N_nJM > 0) then 
                    !% --- Adjust JB-adjacent CC elements using fr_DeltaQ flux changes
                    !%     This fixes conservative flowrate, volume, velocity for upstream (true)
                    !%     and downstream (false) branches. Note that flowrate is already
                    !%     fixed in face_force_JBadjacent_values. Also calls update_auxiliary_data_CC
                    !%     for associated geometry data updates

    
                    call junction_CC_for_JBadjacent (ep_CC_UpstreamOfJunction,   istep, .true.)

                    call junction_CC_for_JBadjacent (ep_CC_DownstreamOfJunction, istep, .false.)
    
                        ! call util_utest_CLprint ('------- TTT  after junction_CC_for_JBadjacent')
    
                    ! if (N_diag > 0) then
                    ! CANNOT USE THIS BECAUSE YOU MIGHT HAVE JUNCTIONS ON BOTH SIDES OF
                    ! A DIAGNOSTIC WHICH MEANS YOU CANNOT RESET THE FAR FACE. INSTEAD WE
                    !% MUST CHANGE THE LOCAL JUNCTION TO MATCH THE DIAGNOSTIC.
                    !         !% --- set conservative fluxes for JB/Diag
                    !     !% --- Adjust JB-adjacent Diag elements and faces using new face values
                    !     !%     This fixes flowrate only. Note that this changes both JB/Diag
                    !     !%     faces and the non-JB face of the Diag element. This ensures that
                    !     !%     the JB face velocity from the junction solution is propagated through
                    !     !%     the entire Diag element to its opposite face.
                    !     call junction_Diag_for_JBadjacent (ep_Diag_UpstreamOfJunction,   istep, .true.)
                    !     call junction_Diag_for_JBadjacent (ep_Diag_DownstreamOfJunction, istep, .false.)
    
                    !         ! ! ! ! ! call util_utest_CLprint ('------- UUU  after junction_Diag_for_JBadjacent')
                    ! end if

                    !% --- update various packs of zeroDepth faces
                    call pack_JB_zeroDepth_interior_faces ()
                    sync all
                    call pack_JB_zeroDepth_shared_faces ()  !% HACK STUB ROUTINE NOT COMPLETE

                    !% --- set face geometry and flowrates where adjacent element is zero
                    !%     only applies to faces with JB on one side
                    call face_zeroDepth (fp_JB_downstream_is_zero_IorS, &
                        fp_JB_upstream_is_zero_IorS,fp_JB_bothsides_are_zero_IorS)

                        ! ! call util_utest_CLprint ('------- TTT.01  after face zerodepth for JB')

                end if
            end if
            !% END JUNCTION SOLUTION
            !%=========================================

                ! ! call util_utest_CLprint ('------- VVV.01  before adjust Vfilter CC')

            !% --- Filter flowrates to remove grid-scale checkerboard
            call adjust_Vfilter_CC ()

                ! ! call util_utest_CLprint ('------- VVV.02  after adjust Vfilter CC')

            if (istep == 1) then 
                !% -- the conservative fluxes from N to N_1 on CC are stored for CC and Diag
                !%    that are NOT adjacent to JB
                !call rk2_store_conservative_fluxes (CCDiag) 
                call rk2_store_conservative_fluxes (ALL) 

                    ! ! call util_utest_CLprint ('------- WWW  after  step 1 store conservative fluxes all')
            end if

                ! call util_utest_CLprint ('------- YYY end of RK step')

        end do

        ! call util_utest_CLprint ('========== ZZZ  end of RK2 ============================')

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
    end subroutine rk2_toplevel_ETM_7
!%
!%==========================================================================            
!%==========================================================================
!%
    subroutine rk2_toplevel_ETM_6
        !%------------------------------------------------------------------
        !% Description:
        !% single RK2 step for explicit time advance of SVE
        !%------------------------------------------------------------------
        !% Declarations:
            integer :: istep, ii, kk, mm
            integer, pointer :: Npack, thisP(:), fup(:), fdn(:)

            character(64) :: subroutine_name = 'rk2_toplevel_ETM_6'
            real(8) :: ConsSum, tval(1)
        !%------------------------------------------------------------------
        !% Preliminaries
        !% --- reset the overflow counter for this time level
            elemR(:,er_VolumeOverFlow) = zeroR     
            elemR(:,er_VolumeArtificialInflow) = zeroR   
        !%-----------------------------------------------------------------

            ! print *,' '
            ! print *, ' '
            ! ! ! ! call util_utest_CLprint ('======= AAA  start of RK2 ==============================')    
 
        !%==================================    
        !% --- Initial adjustments
        istep = zeroI

        !% --- setup for Diagnostic elements
        if (N_diag > 0) then 
            !% STEP A
            !% --- update flowrates for diagnostic elements that are not adjacent to JB
            !call diagnostic_by_type (ep_Diag_notJBadjacent, istep)  
            !20230423 test using all DIAG
            call diagnostic_by_type (ep_Diag, istep)  

                ! ! ! ! call util_utest_CLprint ('------- BBB  after diagnostic')

        !     !% STEP B
        !     !% --- push the flowrate data to faces -- true is upstream, false is downstream
        !     call face_push_elemdata_to_face (ep_Diag_notJBadjacent, fr_Flowrate, er_Flowrate, elemR, .true.)
        !     call face_push_elemdata_to_face (ep_Diag_notJBadjacent, fr_Flowrate, er_Flowrate, elemR, .false.)

        !         ! ! ! ! ! ! call util_utest_CLprint ('------- CCC  after face_push_elemdata_to_face')
        end if

       
        ! if ((setting%Junction%UseJunctionMomentumTF) .and. (N_nJM > 0)) then
        !     !% STEP E
        !     !% --- provide junction outlet JB velocity estimate
            ! call junction_main_velocity (ep_JM)

            !     ! ! ! ! call util_utest_CLprint ('------- DDD  after junction_main_velocity')
        ! end if

        !% STEP F: all diag-notJB are updated on element and faces

        !% --- JB element values to JB faces
        !%     and CC JBadjacent data to face adjacent storage.
        !%     Ensures consistency between diagnostic and JB
        !%     requires cross-image communication
        if (N_nJM > 0) then 
            !% STEPS G- S
            call junction_consistency (istep)

                ! ! ! ! call util_utest_CLprint ('------- EEE  after junction_consistency')
        end if

        !%======================================
        !% --- RK2 SOLUTION
        do istep = 1,2
            ! print *, ' '
            ! print *, istep, '= iSTEP /////////////////////////////////////// ISTEP = ',istep
            ! print *, ' '

            !% --- set the faceR(:,deltaQ) to zero as it is an accumulator for
            !%     changes both caused by face interpolation and changes caused
            !%     by junction solution
            if (N_nJM > 0) then 

                !% STEP rkA
                faceR(:,fr_DeltaQ) = zeroR
                elemR(:,er_DeltaQ) = zeroR


                !% --- store the junction dQdH used in Backwards Euler
                Npack => npack_elemP(ep_JB)
                if (Npack > 0) then     
                    thisP => elemP(1:Npack,ep_JB)

                    !% STEP rkB
                    !% --- dQdH for both CC-adjacent and Diag-adjacent JB 
                    call junction_branch_dQdH (thisP, Npack, istep)

                        !if (istep == 1) 
                        ! ! ! ! call util_utest_CLprint ('------- FFF  after junction_branch_dQdH')
                end if
            end if

            !% STEP rkC
            !% --- Half-timestep advance on CC for U and UVolume
            call rk2_step_ETM_CC (istep)  

                ! ! ! ! call util_utest_CLprint ('------- GGG  after rk2_step')

            !% STEP rkD
            !% --- Update all CC aux variables
            !%     Note, these updates CANNOT depend on face values
            call update_auxiliary_variables_CC (                  &
                ep_CC, ep_CC_Open_Elements, ep_CC_Closed_Elements, &
                .true., .false., dummyIdx)
                
                !if (istep == 1) 
                ! ! ! ! call util_utest_CLprint ('------- HHH  after update_aux...CC step')

            !% STEP rkE
            !% --- zero and small depth adjustment for elements
            call adjust_element_toplevel (CC)

                !if (istep == 1) 
                ! ! ! ! call util_utest_CLprint ('------- HHH.1  after adjust element CC')


            !% --- 2nd step advance on junctions
            if ((istep == 2) .and. (N_nJM > 0)) then 
                Npack => npack_elemP(ep_JM)
                thisP => elemP(1:Npack,ep_JM)
                do mm=1,Npack

                    elemR(thisp(mm),er_Volume) = elemR(thisp(mm),er_Volume_N0) &
                        + elemR(thisP(mm),er_FlowrateLateral) * setting%Time%Hydraulics%Dt

                    do kk=1,max_branch_per_node
                        if (elemSI(thisP(mm)+kk,esi_JunctionBranch_Exists) .ne. oneI) cycle
                        if (mod(kk,2) == 0) then 
                            !% --- downstream branch
                            elemR(thisP(mm),er_Volume) = elemR(thisp(mm),er_Volume) &
                                - faceR(elemI(thisP(mm)+kk,ei_Mface_dL),fr_Flowrate_Conservative) * setting%Time%Hydraulics%Dt
                        else 
                            !% --- upstream branch
                            elemR(thisP(mm),er_Volume) = elemR(thisp(mm),er_Volume) &
                                + faceR(elemI(thisP(mm)+kk,ei_Mface_uL),fr_Flowrate_Conservative) * setting%Time%Hydraulics%Dt
                        end if
                    end do

                    select case (elemSI(thisP(mm),esi_JunctionMain_Type))
                    case (NoStorage)
                        print *, 'CODE ERROR: NoStorage type not handled as of 20230424'
                        call util_crashpoint(55222238)
                    case (ImpliedStorage)
                        elemR(thisP(mm),er_Depth) = elemR(thisP(mm),er_Volume) &
                            / elemSR(thisP(mm),esr_Storage_Plan_Area)
                        elemR(thisp(mm),er_Head) = elemR(thisP(mm),er_Depth) + elemR(thisP(mm),er_Zbottom)
                    case (TabularStorage,FunctionalStorage)
                        print *, 'CODE ERROR: Tabular, functional storage types not handled as of 20230424'
                        call util_crashpoint(729382083)
                    case default
                        print *, 'CODE ERROR: unexpected case default'
                        call util_crashpoint(918732)
                    end select
                end do

                    ! ! ! ! call util_utest_CLprint ('------- HHH.2  after JM 2nd step ')

                call update_auxiliary_variables_JMJB ( .true.)

                    ! ! ! ! call util_utest_CLprint ('------- HHH.3  after update aux JMJB 2nd step ')

                ! if (N_diag > 0) then 
                !     !% STEP A
                !     !% --- update flowrates for diagnostic elements that are not adjacent to JB
                !     !call diagnostic_by_type (ep_Diag_notJBadjacent, istep)  
                !     !20230423 test using all DIAG
                !     call diagnostic_by_type (ep_Diag, istep) 
                ! end if
            end if

            if (N_nJM > 0) then 
                !% STEP rkF
                !% --- Sets JB/CC face interpweights are based on flow on both sides of face               
                Npack => npack_elemP(ep_JB_CC_Adjacent)  
                thisP => elemP(1:Npack,ep_JB_CC_Adjacent)
                ! print *, 'Npack JB CC adjacent ', Npack
                if (Npack > 0) then 
                    call update_interpweights_JB (thisP, Npack, .false.)
                end if
                    !if (istep == 1) 
                    ! ! ! ! call util_utest_CLprint ('------- JJJ  after update weights JB for JB/CC')

                ! !% STEP rkG
                ! !% --- Sets JB/Diag face interpweights for flow to minimum (which should be
                ! !%     same for diagnostic)     
                ! Npack => npack_elemP(ep_JB_Diag_Adjacent)
                ! thisP => elemP(1:Npack,ep_JB_Diag_Adjacent)
                ! if (Npack > 0) then 
                !     call update_interpweights_JB (thisP, Npack, .true.)
                ! end if

                !     if (istep == 1) ! ! ! ! ! ! call util_utest_CLprint ('------- KKK  after update weights JB for JB/Diag')
            end if

            !% STEP rkH
            !% --- interpolate all data to faces
            !% NOTE 20230419 THIS IS DIFFERENT THAN IN JUNCTION_TOPLEVEL_3
            sync all
            call face_interpolation(fp_noBC_IorS, .true., .true., .true., .false., .true.) 

                !if (istep == 1) 
                ! ! ! ! call util_utest_CLprint ('------- LLL  after face interp')


            !% STEP rkE.1 
            call pack_CC_zeroDepth_interior_faces ()
            sync all
            call pack_CC_zeroDepth_shared_faces ()

            call face_zeroDepth_geometry_interior(fp_CC_downstream_is_zero_IorS)
            call face_zeroDepth_geometry_interior(fp_CC_upstream_is_zero_IorS)
            call face_zeroDepth_geometry_interior(fp_CC_bothsides_are_zero_IorS)

                !if (istep == 1) 
                ! ! ! ! ! ! ! call util_utest_CLprint ('------- LLL.01  after face zerodepth ')

            call face_zeroDepth_geometry_shared(fp_CC_downstream_is_zero_IorS)
            call face_zeroDepth_geometry_shared(fp_CC_upstream_is_zero_IorS)
            call face_zeroDepth_geometry_shared(fp_CC_bothsides_are_zero_IorS)

                !if (istep == 1) 
                ! ! ! ! ! ! ! call util_utest_CLprint ('------- LLL.02  after face zerodepth ')

            call face_zeroDepth_flowrates_interior(fp_CC_downstream_is_zero_IorS)
            call face_zeroDepth_flowrates_interior(fp_CC_upstream_is_zero_IorS)
            call face_zeroDepth_flowrates_interior(fp_CC_bothsides_are_zero_IorS)

                !if (istep == 1) 
                ! ! ! ! ! ! ! call util_utest_CLprint ('------- LLL.03  after face zerodepth ')

            call face_zeroDepth_flowrates_shared(fp_CC_downstream_is_zero_IorS)
            call face_zeroDepth_flowrates_shared(fp_CC_upstream_is_zero_IorS)
            call face_zeroDepth_flowrates_shared(fp_CC_bothsides_are_zero_IorS)

                !if (istep == 1) 
                ! ! ! ! call util_utest_CLprint ('------- LLL.04  after face zerodepth ')

            call face_flowrate_for_openclosed_elem (ep_CCDiag)
            !% HACK -- FOR SHARED WE NEED TO IDENTIFY ANY SHARED FACES THAT HAVE AN
            !% UPSTREAM elemR(:,er_Setting) = 0.0 and set their face Q and V to zero.


               !if (istep == 1) 
                ! ! ! ! call util_utest_CLprint ('------- MMM  after face zerodepth ')
         
            !%=========================================
            !% BEGIN JUNCTION SOLUTION
            if ((istep == 1) .and. (N_nJM > 0)) then 
            ! if (N_nJM > 0) then 
                !% --- Backwards-Euler junction computation for Head and flowrate changes
                Npack => npack_elemP(ep_JM)   
                thisP => elemP(1:Npack,ep_JM)
                    !% STEPS J1 - J22
                    call junction_calculation_4 (thisP, Npack, istep)

                        ! ! ! ! call util_utest_CLprint ('------- NNN after junction calculation_4')

                    !% STEP J23
                    !% HACK --- PREISSMANN SLOT -- STEP J23.  May need something here or in junction_calculation_4

                    !% STEP J24
                    !% --- for cases where dH is limited, reset JB flowrates for mass conservation
                   call junction_mass_conservation (ep_JM, istep)

                        ! ! ! ! call util_utest_CLprint ('------- OOO after junction_mass_conservation ')

                    !% HACK --OLD CODE HAD UPDATE AUX_JMJB AND ADJUST ELEMENTS BEFORE MASS CONSERVATION

                    !% STEP J25
                    !% --- update auxiliary variable on JM and JB
                    !%     true indicates flowrate interp on all JB are set to minimum.
                    call update_auxiliary_variables_JMJB ( .true.)

                        ! ! ! ! call util_utest_CLprint ('------- PPP  after update_aux ..JMJB')

            end if   

            !% STEP J26
            !% --- adjust JB and JM for small or zero depth
            call adjust_element_toplevel (JB)
            call adjust_element_toplevel (JM)

                ! ! ! ! call util_utest_CLprint ('------- QQQ  after adjust element JB JM')

            if ((istep == 1) .and. (N_nJM > 0)) then 
            ! if (N_nJM > 0) then 
                !% STEP J27
                !% --- force the JB values to the faces for upstream (true)
                !%     and downstream (false) branches.
                !%     Forces elem  flowrate, deltaQ, and head to face
                !%     also computes new depths and areas for faces
                call face_force_JBadjacent_values (ep_JM, .true.)
                call face_force_JBadjacent_values (ep_JM, .false.)

                    ! ! ! ! call util_utest_CLprint ('------- RRR  after face_force_JBadjacent')

            end if

            sync all 

            !% STEP J28
            !%================================
            !% --- HACK need to force face data on all JB faces for image containing JB
            !%     element to the connected image. This is a direct transfer of face
            !%     Data without transferring ghost element data or interpolation
            !%================================


            if ((istep == 1) .and. (N_nJM > 0)) then
            ! if (N_nJM > 0) then
                !% STEP J29, J30
                !% --- Adjust JB-adjacent CC elements using new face values
                !%     This fixes flowrate, volume, and geometry for upstream (true)
                !%     and downstream (false) branches. Calls update_auxiliary_data_CC
                !%     for associated data

                call junction_CC_for_JBadjacent (ep_CC_UpstreamOfJunction,   istep, .true.)
                call junction_CC_for_JBadjacent (ep_CC_DownstreamOfJunction, istep, .false.)

                    !if (istep == 1) 
                    ! ! ! ! call util_utest_CLprint ('------- SSS  after junction_CC_for_JBadjacent')

                ! if (N_diag > 0) then
                ! CANNOT USE THIS BECAUSE YOU MIGHT HAVE JUNCTIONS ON BOTH SIDES OF
                ! A DIAGNOSTIC WHICH MEANS YOU CANNOT RESET THE FAR FACE. INSTEAD WE
                !% MUST CHANGE THE LOCAL JUNCTION TO MATCH THE DIAGNOSTIC.
                !         !% --- set conservative fluxes for JB/Diag
                !     !% --- Adjust JB-adjacent Diag elements and faces using new face values
                !     !%     This fixes flowrate only. Note that this changes both JB/Diag
                !     !%     faces and the non-JB face of the Diag element. This ensures that
                !     !%     the JB face velocity from the junction solution is propagated through
                !     !%     the entire Diag element to its opposite face.
                !     call junction_Diag_for_JBadjacent (ep_Diag_UpstreamOfJunction,   istep, .true.)
                !     call junction_Diag_for_JBadjacent (ep_Diag_DownstreamOfJunction, istep, .false.)

                !         ! ! ! ! ! call util_utest_CLprint ('------- TTT  after junction_Diag_for_JBadjacent')
                ! end if

                !% HACK OLD CODE HAD adjust_Face_toplevel () called here in junction_toplevel_3

                ! !% --- experiment with suppressing jump downstream of JB
                ! Npack => npack_elemP(ep_CC_DownstreamOfJunction)   
                ! thisP => elemP(1:Npack,ep_CC_DownstreamOfJunction)
                ! do ii=1,Npack 
                !     mm = thisP(ii)
                !     !% --- downstream face
                !     fdn => elemI(mm,ei_Mface_dL)
                !     if (faceI(fdn,fi_jump_type) == jump_from_upstream ) then 
                !         !% --- force elem area and velocity to be the face value
                !         elemR(mm,er_Area) = faceR(fdn,fr_Area_d)
                !         elemR(mm,er_Velocity) = elemR(mm,er_Flowrate) / elemR(mm,er_Area)
                !     end if
                ! end do



            end if
            !% END JUNCTION SOLUTION
            !%======================================

            !% HACK OLD CODE HAD CALL TO ALL face INTERPOLATION AFTER junction_toplevel_3, then 
            !% call to adjust-face_toplevel(fp_noBC_IorS)

            !% STEP rkJ
            sync all

            !% STEP rkK
            !%================================  
            !% HACK -- need to force flowrate data push from any fp_Diag face on the image 
            !% containing the diagnostic element to the associated face on the image that
            !% contains adjacent element (i.e., either CC or JB)
            !%================================

            ! if (N_diag > 0) then
            !     !% STEP rkL
            !     !% --- Adjust Velocities on all Diag faces. The real target for this
            !     !%     is the faces of a Diag where the Diag is adjacent ot a JB and its 
            !     !%     non-adjacent face has flowrate changed for mass conservation.
            !     call face_velocities(fp_Diag_all, .true.)

            !         if (istep == 1) ! ! ! ! ! ! call util_utest_CLprint ('------- VVV  after face velocities')

            ! end if

            !% 20230428 FACE FROUDE NUMBER
            Npack => npack_faceP(fp_noBC_IorS)
            thisP => faceP(1:Npack,fp_noBC_IorS)
            where (faceR(thisP,fr_Depth_u) > setting%ZeroValue%Depth)
                faceR(thisP,fr_FroudeNumber_u) = faceR(thisP,fr_Velocity_u) / sqrt(setting%Constant%gravity * faceR(thisP,fr_Area_u))
            endwhere
            where (faceR(thisP,fr_Depth_d) > setting%ZeroValue%Depth)
                faceR(thisP,fr_FroudeNumber_d) = faceR(thisP,fr_Velocity_d) / sqrt(setting%Constant%gravity * faceR(thisP,fr_Area_d))
            endwhere

                ! ! ! ! call util_utest_CLprint ('------- VVV.01  before adjust Vfilter CC')

            Npack => npack_elemP(ep_CC)    
            thisP => elemP(1:Npack,ep_CC)
            fup   => elemI(:,ei_Mface_uL)
            fdn   => elemI(:,ei_Mface_dL)
            where ( (faceR(fup(thisP),fr_FroudeNumber_d) < oneR) .and. &
                    (faceR(fdn(thisP),fr_FroudeNumber_u) < oneR) .and. &
                    (elemR(thisP,er_FroudeNumber)        > oneR) )
                elemR(thisP,er_Head)      = onehalfR * (faceR(fup(thisP),fr_Head_d)      + faceR(fdn(thisP),fr_Head_u))    
                ! elemR(thisP,er_Area)      = onehalfR * (faceR(fup(thisP),fr_Area_d)      + faceR(fdn(thisP),fr_Area_u))    
                !elemR(thisP,er_Depth)     = onehalfR * (faceR(fup(thisP),fr_Depth_d)     + faceR(fdn(thisP),fr_Depth_u))
                ! elemR(thisP,er_Velocity)  = onehalfR * (faceR(fup(thisP),fr_Velocity_d)  + faceR(fdn(thisP),fr_Velocity_u))    
            endwhere

            !% STEP rkM
            !% --- Filter flowrates to remove grid-scale checkerboard
            call adjust_Vfilter_CC () !% brhADDED 20230418


            
                !  ! ! ! ! call util_utest_CLprint ('------- VVV.02  after adjust Vfilter CC')


            !% HACK -- for better coding, change the rk2_store_conservative call to use the
            !% correct fp_... column as the argument

            if (istep == 1) then 

                !% STEP rkN
                !% -- the conservative fluxes from N to N_1 on CC are stored for CC and Diag
                !%    that are NOT adjacent to JB
                !call rk2_store_conservative_fluxes (CCDiag) 
                call rk2_store_conservative_fluxes (ALL) 

                    ! ! ! ! call util_utest_CLprint ('------- MMM  after  step 1 store conservative fluxes all')
            end if

            !% NOTE CONSERVATIVE FLUX FOR JB IS FIXED IN junction_CC_for_JBadjacent
            ! if ((N_nJM > 0) .and. (istep == 2)) then 
            !     !% STEP rkO
            !     !% --- store conservative fluxes on JB and Diag faces. Note that the CC/CC and
            !     !%     CC/Diag faces have flowrates set in the first RK2 step. Here we are 
            !     !%     getting the final correction since the volume in JM  at the end of the
            !     !%     junction correction is set from er_Volume_N0 in junction_update_Qdependent_values
            !     !%     and the volumes/flowrates in CC adjacent to JB and JBDiag have been adjusted
            !     !%     for the JM/JB solution.
            !     call rk2_store_conservative_fluxes (JBDiag) 

            !         ! ! ! ! ! call util_utest_CLprint ('------- UUU  after step 2 store conservative flux JBDiag')
            ! end if


            ! !% STEP rkP
            ! !% --- reset the overflow counter (we only save conservation in the 2nd step)
            ! if (istep == 1) then 
            !     elemR(:,er_VolumeOverFlow) = zeroR
            !     elemR(:,er_VolumeArtificialInflow) = zeroR 
            ! end if

            !if (istep == 1) 
            ! ! ! ! call util_utest_CLprint ('------- VVV end of RK step')


        end do

        ! ! ! ! call util_utest_CLprint ('========== ZZZ  end of RK2 ============================')

        ! tval = maxval(elemR(2:(N_elem(1)-1),er_Flowrate)) 
        ! if (tval(1) > oneR) then
        !     print *, 'Flowrate'
        !     do ii=1,N_elem(1)-1
        !         print *, ii, trim(reverseKey(elemI(ii,ei_elementType))), elemR(ii,er_Flowrate)
        !     end do
        !     stop 6098723
        ! end if

    end subroutine rk2_toplevel_ETM_6
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rk2_toplevel_ETM_5
        !%------------------------------------------------------------------
        !% Description:
        !% single RK2 step for explicit time advance of SVE
        !%------------------------------------------------------------------
        !% Declarations:
            integer :: istep, ii
            integer, pointer :: Npack, thisp(:)

            character(64) :: subroutine_name = 'rk2_toplevel_ETM_5'
            real(8) :: ConsSum
        !%------------------------------------------------------------------
        !% Preliminaries
        !% --- reset the overflow counter for this time level
            elemR(:,er_VolumeOverFlow) = zeroR     
            elemR(:,er_VolumeArtificialInflow) = zeroR   
        !%-----------------------------------------------------------------

            print *,' '
            print *, ' '
            ! ! ! ! ! call util_utest_CLprint ('======= AAA  start of RK2 ==============================')    
 
        !%==================================    
        !% --- Initial adjustments
        istep = zeroI

        !% --- setup for Diagnostic elements
        ! if (N_diag > 0) then 
        !     !% STEP A
        !     !% --- update flowrates for diagnostic elements that are not adjacent to JB
        !     call diagnostic_by_type (ep_Diag_notJBadjacent, istep)  

        !         ! ! ! ! ! ! call util_utest_CLprint ('------- BBB  after diagnostic')

        !     !% STEP B
        !     !% --- push the flowrate data to faces -- true is upstream, false is downstream
        !     call face_push_elemdata_to_face (ep_Diag_notJBadjacent, fr_Flowrate, er_Flowrate, elemR, .true.)
        !     call face_push_elemdata_to_face (ep_Diag_notJBadjacent, fr_Flowrate, er_Flowrate, elemR, .false.)

        !         ! ! ! ! ! ! call util_utest_CLprint ('------- CCC  after face_push_elemdata_to_face')
        ! end if

       
        ! if ((setting%Junction%UseJunctionMomentumTF) .and. (N_nJM > 0)) then
        !     !% STEP E
        !     !% --- provide junction outlet JB velocity estimate
        !     call junction_main_velocity (ep_JM)

        !         ! ! ! ! ! ! call util_utest_CLprint ('------- DDD  after junction_main_velocity')
        ! end if

        !% STEP F: all diag-notJB are updated on element and faces

        !% --- JB element values to JB faces
        !%     and CC JBadjacent data to face adjacent storage.
        !%     Ensures consistency between diagnostic and JB
        !%     requires cross-image communication
        if (N_nJM > 0) then 
            !% STEPS G- S
            call junction_consistency (istep)

                ! ! ! ! ! ! call util_utest_CLprint ('------- EEE  after junction_consistency')
        end if

        !%======================================
        !% --- RK2 SOLUTION
        do istep = 1,2
            ! print *, ' '
            ! print *, istep, '= iSTEP /////////////////////////////////////// ISTEP = ',istep
            ! print *, ' '

            !% --- set the faceR(:,deltaQ) to zero as it is an accumulator for
            !%     changes both caused by face interpolation and changes caused
            !%     by junction solution
            if (N_nJM > 0) then 

                !% STEP rkA
                faceR(:,fr_DeltaQ) = zeroR


                !% --- store the junction dQdH used in Backwards Euler
                Npack => npack_elemP(ep_JB)
                if (Npack > 0) then     
                    thisP => elemP(1:Npack,ep_JB)

                    !% STEP rkB
                    !% --- dQdH for both CC-adjacent and Diag-adjacent JB 
                    call junction_branch_dQdH (thisP, Npack, istep)

                        !if (istep == 1) 
                        ! ! ! ! ! ! call util_utest_CLprint ('------- FFF  after junction_branch_dQdH')
                end if
            end if

            !% STEP rkC
            !% --- Half-timestep advance on CC for U and UVolume
            call rk2_step_ETM_CC (istep)  

                ! ! ! ! ! ! call util_utest_CLprint ('------- GGG  after rk2_step')

            !% STEP rkD
            !% --- Update all CC aux variables
            !%     Note, these updates CANNOT depend on face values
            call update_auxiliary_variables_CC (                  &
                ep_CC, ep_CC_Open_Elements, ep_CC_Closed_Elements, &
                .true., .false., dummyIdx)
                
                !if (istep == 1) 
                ! ! ! ! ! ! call util_utest_CLprint ('------- HHH  after update_aux...CC step')

            !% STEP rkE
            !% --- zero and small depth adjustment for elements
            call adjust_element_toplevel (CC)

                !if (istep == 1) 
                ! ! ! ! ! ! call util_utest_CLprint ('------- HHH.1  after adjust element CC')

            if (N_nJM > 0) then 
                !% STEP rkF
                !% --- Sets JB/CC face interpweights are based on flow on both sides of face               
                Npack => npack_elemP(ep_JB_CC_Adjacent)  
                thisP => elemP(1:Npack,ep_JB_CC_Adjacent)
                ! print *, 'Npack JB CC adjacent ', Npack
                if (Npack > 0) then 
                    call update_interpweights_JB (thisP, Npack, .false.)
                end if

                    !if (istep == 1) 
                    ! ! ! ! ! ! call util_utest_CLprint ('------- JJJ  after update weights JB for JB/CC')

                ! !% STEP rkG
                ! !% --- Sets JB/Diag face interpweights for flow to minimum (which should be
                ! !%     same for diagnostic)     
                ! Npack => npack_elemP(ep_JB_Diag_Adjacent)
                ! thisP => elemP(1:Npack,ep_JB_Diag_Adjacent)
                ! if (Npack > 0) then 
                !     call update_interpweights_JB (thisP, Npack, .true.)
                ! end if

                !     if (istep == 1) ! ! ! ! ! ! call util_utest_CLprint ('------- KKK  after update weights JB for JB/Diag')
            end if

            !% STEP rkH
            !% --- interpolate all data to faces
            !% NOTE 20230419 THIS IS DIFFERENT THAN IN JUNCTION_TOPLEVEL_3
            sync all
            call face_interpolation(fp_noBC_IorS, .true., .true., .true., .false., .true.) 

                !if (istep == 1) 
                ! ! ! ! ! ! call util_utest_CLprint ('------- LLL  after face interp')


            !% STEP rkE.1 
            call pack_CC_zeroDepth_interior_faces ()
            sync all
            call pack_CC_zeroDepth_shared_faces ()

            call face_zeroDepth_geometry_interior(fp_CC_downstream_is_zero_IorS)
            call face_zeroDepth_geometry_interior(fp_CC_upstream_is_zero_IorS)
            call face_zeroDepth_geometry_interior(fp_CC_bothsides_are_zero_IorS)

                !if (istep == 1) 
                ! ! ! ! ! ! ! call util_utest_CLprint ('------- LLL.01  after face zerodepth ')

            call face_zeroDepth_geometry_shared(fp_CC_downstream_is_zero_IorS)
            call face_zeroDepth_geometry_shared(fp_CC_upstream_is_zero_IorS)
            call face_zeroDepth_geometry_shared(fp_CC_bothsides_are_zero_IorS)

                !if (istep == 1) 
                ! ! ! ! ! ! ! call util_utest_CLprint ('------- LLL.02  after face zerodepth ')

            call face_zeroDepth_flowrates_interior(fp_CC_downstream_is_zero_IorS)
            call face_zeroDepth_flowrates_interior(fp_CC_upstream_is_zero_IorS)
            call face_zeroDepth_flowrates_interior(fp_CC_bothsides_are_zero_IorS)

                !if (istep == 1) 
                ! ! ! ! ! ! ! call util_utest_CLprint ('------- LLL.03  after face zerodepth ')

            call face_zeroDepth_flowrates_shared(fp_CC_downstream_is_zero_IorS)
            call face_zeroDepth_flowrates_shared(fp_CC_upstream_is_zero_IorS)
            call face_zeroDepth_flowrates_shared(fp_CC_bothsides_are_zero_IorS)

                !if (istep == 1) 
                ! ! ! ! ! ! ! call util_utest_CLprint ('------- LLL.04  after face zerodepth ')

            call face_flowrate_for_openclosed_elem (ep_CCDiag)
            !% HACK -- FOR SHARED WE NEED TO IDENTIFY ANY SHARED FACES THAT HAVE AN
            !% UPSTREAM elemR(:,er_Setting) = 0.0 and set their face Q and V to zero.

               !if (istep == 1) 
            !    ! ! ! ! ! call util_utest_CLprint ('------- MMM  after face zerodepth ')

            !%=========================================
            !% BEGIN JUNCTION SOLUTION
            if (N_nJM > 0) then 
                !% --- Backwards-Euler junction computation for Head and flowrate changes
                Npack => npack_elemP(ep_JM)
                thisP => elemP(1:Npack,ep_JM)
                !% STEPS J1 - J22
                call junction_calculation_4 (thisP, Npack, istep)

                    ! ! ! ! ! ! call util_utest_CLprint ('------- NNN after junction calculation_4')

                !% STEP J23
                !% HACK --- PREISSMANN SLOT -- STEP J23.  May need something here or in junction_calculation_4

                !% HACK JUNCTION MASS CONSERVATION DELAYED UNTIL JUST BEFORE adjust_vflilter_CC IN OLD CODE

                !% HACK OLD CODE DID NOT HAVE A SECOND CALL TO JUNCITON MASS CONSERVATION
                !% WHICH WAS PROBABLY A PROBLEM!


                !% STEP J24
                !% --- for cases where dH is limited, reset JB flowrates for mass conservation
                call junction_mass_conservation (ep_JM, istep)

                    ! ! ! ! ! ! call util_utest_CLprint ('------- OOO after junction_mass_conservation ')

                !% HACK --OLD CODE HAD UPDATE AUX_JMJB AND ADJUST ELEMENTS BEFORE MASS CONSERVATION

                !% STEP J25
                !% --- update auxiliary variable on JM and JB
                !%     true indicates flowrate interp on all JB are set to minimum.
                call update_auxiliary_variables_JMJB ( .true.)

                    ! ! ! ! ! ! call util_utest_CLprint ('------- PPP  after update_aux ..JMJB')

                !% STEP J26
                !% --- adjust JB and JM for small or zero depth
                call adjust_element_toplevel (JB)
                call adjust_element_toplevel (JM)

                    ! ! ! ! ! ! call util_utest_CLprint ('------- QQQ  after adjust element JB JM')

                !% STEP J27
                !% --- force the JB values to the faces for upstream (true)
                !%     and downstream (false) branches.
                !%     Forces elem  flowrate, deltaQ, and head to face
                !%     also computes new depths and areas for faces
                call face_force_JBadjacent_values (ep_JM, .true.)
                call face_force_JBadjacent_values (ep_JM, .false.)

                    ! ! ! ! ! ! call util_utest_CLprint ('------- RRR  after face_force_JBadjacent')

            end if


            sync all 

            !% STEP J28
            !%================================
            !% --- HACK need to force face data on all JB faces for image containing JB
            !%     element to the connected image. This is a direct transfer of face
            !%     Data without transferring ghost element data or interpolatoin
            !%================================



            if (N_nJM > 0) then 
                !% STEP J29, J30
                !% --- Adjust JB-adjacent CC elements using new face values
                !%     This fixes flowrate, volume, and geometry for upstream (true)
                !%     and downstream (false) branches. Calls update_auxiliary_data_CC
                !%     for associated data

                call junction_CC_for_JBadjacent (ep_CC_UpstreamOfJunction,   istep, .true.)
                call junction_CC_for_JBadjacent (ep_CC_DownstreamOfJunction, istep, .false.)

                    !if (istep == 1) 
                    ! ! ! ! ! ! call util_utest_CLprint ('------- SSS  after junction_CC_for_JBadjacent')

                ! if (N_diag > 0) then
                !     !% --- Adjust JB-adjacent Diag elements and faces using new face values
                !     !%     This fixes flowrate only. Note that this changes both JB/Diag
                !     !%     faces and the non-JB face of the Diag element. This ensures that
                !     !%     the JB face velocity from the junction solution is propagated through
                !     !%     the entire Diag element to its opposite face.
                !     call junction_Diag_for_JBadjacent (ep_Diag_UpstreamOfJunction,   istep, .true.)
                !     call junction_Diag_for_JBadjacent (ep_Diag_DownstreamOfJunction, istep, .false.)

                !         if (istep == 1) ! ! ! ! ! ! call util_utest_CLprint ('------- TTT  after junction_Diag_for_JBadjacent')
                ! end if

                !% HACK OLD CODE HAD adjust_Face_toplevel () called here in junction_toplevel_3

            end if
            !% END JUNCTION SOLUTION
            !%======================================


            !% HACK OLD CODE HAD CALL TO ALL face INTERPOLATION AFTER junction_toplevel_3, then 
            !% call to adjust-face_toplevel(fp_noBC_IorS)


            !% STEP rkJ
            sync all

            !% STEP rkK
            !%================================  
            !% HACK -- need to force flowrate data push from any fp_Diag face on the image 
            !% containing the diagnostic element to the associated face on the image that
            !% contains adjacent element (i.e., either CC or JB)
            !%================================

            ! if (N_diag > 0) then
            !     !% STEP rkL
            !     !% --- Adjust Velocities on all Diag faces. The real target for this
            !     !%     is the faces of a Diag where the Diag is adjacent ot a JB and its 
            !     !%     non-adjacent face has flowrate changed for mass conservation.
            !     call face_velocities(fp_Diag_all, .true.)

            !         if (istep == 1) ! ! ! ! ! ! call util_utest_CLprint ('------- VVV  after face velocities')
            ! end if

            !% STEP rkM
            !% --- Filter flowrates to remove grid-scale checkerboard
            call adjust_Vfilter_CC () !% brhADDED 20230418


            !% HACK -- for better coding, change the rk2_store_conservative call to use the
            !% correct fp_... column as the argument

            ! if (istep == 1) then 

            !     !HACK THIS MAY NEED A DIFFERENT ARGUMENT TO GET BC
                
            !     !% STEP rkN
            !     !% -- the conservative fluxes from N to N_1 on CC are stored for CC and Diag
            !     !%    that are NOT adjacent to JB
            !     call rk2_store_conservative_fluxes (CCDiag) 

            !         if (istep == 1) ! ! ! ! ! ! call util_utest_CLprint ('------- MMM  after  step 1 store consertave fluxes CCDiag')
            ! end if

            ! if ((N_nJM > 0) .and. (istep == 2)) then 
            !     !% STEP rkO
            !     !% --- store conservative fluxes on JB and Diag faces. Note that the CC/CC and
            !     !%     CC/Diag faces have flowrates set in the first RK2 step. Here we are 
            !     !%     getting the final correction since the volume in JM  at the end of the
            !     !%     junction correction is set from er_Volume_N0 in junction_update_Qdependent_values
            !     !%     and the volumes/flowrates in CC adjacent to JB and JBDiag have been adjusted
            !     !%     for the JM/JB solution.
            !     call rk2_store_conservative_fluxes (JBDiag) 

            !         if (istep == 1) ! ! ! ! ! ! call util_utest_CLprint ('------- UUU  after step 2 store conservative flux JBDiag')
            ! end if


            ! !% STEP rkP
            ! !% --- reset the overflow counter (we only save conservation in the 2nd step)
            ! if (istep == 1) then 
            !     elemR(:,er_VolumeOverFlow) = zeroR
            !     elemR(:,er_VolumeArtificialInflow) = zeroR 
            ! end if

            !if (istep == 1) 
            ! ! ! ! ! ! call util_utest_CLprint ('------- VVV end of RK step')

        end do

        ! ! ! ! ! ! call util_utest_CLprint ('========== ZZZ  end of RK2 ============================')

    end subroutine rk2_toplevel_ETM_5
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rk2_toplevel_ETM_4 ()
        !%------------------------------------------------------------------
        !% Description:
        !% single RK2 step for explicit time advance of SVE
        !%------------------------------------------------------------------
        !% Declarations:
            integer :: istep, ii
            integer :: whichTM
            character(64) :: subroutine_name = 'rk2_toplevel_ETM_NEW'

            real(8) :: ConsSum
        !%------------------------------------------------------------------
        !% Preliminaries
            !% --- set the time-marching type for this routine
            whichTM = ETM
            !% --- reset the overflow counter for this time level
            elemR(:,er_VolumeOverFlow) = zeroR     
            elemR(:,er_VolumeArtificialInflow) = zeroR   
        !%-----------------------------------------------------------------
        !% Aliases
        !%-----------------------------------------------------------------
        !% --- first step of the RK2
        istep = oneI

            !print *, ' '
         ! ! ! ! ! ! call util_utest_CLprint ('======= AAA  start of RK2 ==============================')
            !print *, 'TEST20230327   GGG',elemR(22,er_Head), elemR(22,er_Zbottom)
            !stop 709873

        !% --- NOTE: Dynamic manning's n not included in this routine.

        !% --- compute Force Main Manning's N for force main elements
        !%     Note this is only computed at the start of the RK2 and the
        !%     value is held constant throughout the step.
        !%     20230128 -- modified pack so only CC elements, which should not matter
        if (setting%Solver%ForceMain%AllowForceMainTF) call forcemain_ManningsN ()   


        ! THOUGHT -- COMPUTE the JM velocity for outflow branches
        ! If this is larger than the JB velocity, then use this for the JB 
        ! flowrate interpolation.
        if (setting%Junction%UseJunctionMomentumTF) then
            call junction_main_velocity (ep_JM)
        end if


        !% --- store the junction face terms for time n 
        call junction_face_terms_linear (ep_JB, istep)    
        call util_crashstop(6228745)

            ! ! ! ! ! ! call util_utest_CLprint ('------- BBB  after junction_face_terms')

        !% --- Half-timestep advance on CC for U and UVolume
        call rk2_step_ETM_CC (istep)   
        call util_crashstop(7399287)

            ! ! ! ! ! ! call util_utest_CLprint ('------- CCC  after rk2_step_ETM_cc')
        
        !% --- Update all CC aux variables
        !%     Note, these updates CANNOT depend on face values
        call update_auxiliary_variables_CC (&
            ep_CC, ep_CC_Open_Elements, ep_CC_Closed_Elements, &
            .true., .false., dummyIdx)

            ! ! ! ! ! ! call util_utest_CLprint ('------- DDD  after update_auxiliary_variables_CC')
    
        call adjust_element_toplevel (CC)
        call util_crashstop (80298733)

                ! ! ! ! ! ! call util_utest_CLprint ('------- EEE after adjust element ... CC')

        !% --- Implicit junction for JB flowrate and head
        call junction_toplevel_3 (istep)
        call util_crashstop (11284298)
    
            ! ! ! ! ! ! call util_utest_CLprint ('------- FFF  after junction_toplevel')  
            
            ! print *, ' '
            ! print *,  'resid ',junction_conservation_residual(13)
            ! print *, ' '

        !% --- face interpolation (all faces)
        !%     because of forceJByn = .true. in last call to update_auxiliary_variables_JMJB 
        !%     this will favor JB element values at JB faces, which is the mass conservative
        !%     value
        sync all
        call face_interpolation(fp_noBC_IorS,.true.,.true.,.true.,.false.,.false.)
        call util_crashstop(72112342)

        ! print *, ' face flowrate debug'
        ! print *, elemR(78,er_Flowrate), faceR(73, fr_Flowrate), elemR(80,er_Flowrate)

            ! ! ! ! ! ! call util_utest_CLprint ('------- GGG  after face_interpolation')

            ! print *, ' '
            ! print *,  'resid ',junction_conservation_residual(13)
            ! print *, ' '

        call adjust_face_toplevel (fp_noBC_IorS)    
            
            ! ! ! ! ! ! call util_utest_CLprint ('------- HHH  after face ad hoc adjustments')

            ! print *, ' '
            ! print *,  'resid ',junction_conservation_residual(13)
            ! print *, ' '
        
        !% --- update diagnostic elements and faces
        call diagnostic_toplevel (ep_Diag, fp_Diag_IorS, istep)  
        call util_crashstop(982373)

            ! ! ! ! ! ! call util_utest_CLprint ('-------III after diagnostic_toplevel')

            ! print *, ' '
            ! print *,  'resid ',junction_conservation_residual(13)
            ! print *, ' '



        !% --- ensure mass conservation at junctions
        call junction_mass_conservation (ep_JM, istep)
        call util_crashstop(602984)

            ! ! ! ! ! ! call util_utest_CLprint ('-------JJJ after junction_mass_conservation')

        !% --- fix flowrates in diagnostic elements adjacent to JB
        ! call diagnostic_fix_JB_adjacent ()
        !call util_crashstop(998734)

            ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('-------JJJ01 after diagnostic_fix_JB_adjacent')

        !% --- RK2 solution step -- check culverts
        !%     HACK: Can this be included in diagnostic?
        !%     NEED CHECKING FOR JUNCTIONS
        ! call culvert_toplevel() NEED TO REWORK  20230317
        ! call util_crashstop(669753)
        ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- JJJ02  after culvert_toplevel')

        !% --- RK2 solution step  -- make ad hoc adjustments
        call adjust_Vfilter_CC () !% this is useful in reducing flow oscillations
        call util_crashstop(139871)

            ! ! ! ! ! ! call util_utest_CLprint ('------- JJJ03  after adjust_Vfilter_CC')
        
        !% -- the conservative fluxes from N to N_1 are the values just before the second RK2 step
        call rk2_store_conservative_fluxes (ALL) 

            ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- JJJ04  after rk2_store_conservative_fluxes')

        !% --- reset the overflow counter (we only save conservation in the 2nd step)
        elemR(:,er_VolumeOverFlow) = zeroR
        elemR(:,er_VolumeArtificialInflow) = zeroR 

        !% --------------------------------------------------------------------------
        !% --------------------------------------------------------------------------
        !% --- RK2 solution step -- RK2 second step for ETM 
        !%     CC advanced for continuity and momentum
        istep = twoI

        if (setting%Junction%UseJunctionMomentumTF) then 
            call junction_main_velocity (ep_JM)
        end if

        !% --- store the junction face terms for time n
        call junction_face_terms_linear (ep_JB, istep) 
        call util_crashstop(6281287)  

            ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- MMM  after junction_face_terms')

        !% -- Full-timestep advance on CC for U and UVolume
        call rk2_step_ETM_CC (istep)
        call util_crashstop(3298744)

            ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- NNN  after rk2_step_ETM_cc')

        !% --- update all CC aux variables
        !%     Note, these updates CANNOT depend on face values
        call update_auxiliary_variables_CC (&
            ep_CC, ep_CC_Open_Elements, ep_CC_Closed_Elements, &
            .true., .false., dummyIdx)
        call util_crashstop(2969983)

           ! ! ! ! ! ! ! call util_utest_CLprint ('------- OOO after update_auxiliary_variables_CC')

        call adjust_element_toplevel (CC)
        call util_crashstop(2508773)    

            ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- PPP  after adjust... CC')

        !% --- Implicit junction for JB flowrate and head
        call junction_toplevel_3 (istep)
        call util_crashstop (11287322)

            ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- RRR  after junction_toplevel')



        !% NOTE: modifications after this point only affect the non-conservative
        !% values at time n+1 and do not affect either the conservative fluxes from
        !% time n->n+1 or the volume at n+1

            ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('-------RRR022 after store conservative fluxes')


        sync all
        call face_interpolation(fp_noBC_IorS,.true.,.true.,.true.,.false.,.false.)
        call util_crashstop(72129873)

            ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- SSS  after face_interpolation')

       call adjust_face_toplevel (fp_noBC_IorS)
       call util_crashstop(5098723)

            ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- SSS02  after face ad hoc adjumstents')

        !% --- update diagnostic elements and faces
        call diagnostic_toplevel (ep_Diag, fp_Diag_IorS, istep)
        call util_crashstop(982332)

        !% -- the conservative fluxes from on JB are after the 2nd RK step
        call rk2_store_conservative_fluxes (JB) 


            ! ! if (setting%Time%Step > 5296) then
            !     print *, ' '
            !     print *, ' LOCAL CONSERVATION ======================================================='
            !     print *, 'before weir'
            !     ConsSum = zeroR
            !     do ii=1,50
            !         print *, ii, (elemR(ii,er_Volume) - elemR(ii,er_Volume_N0)) &
            !                 -(   faceR(elemI(ii,ei_Mface_uL),fr_Flowrate_Conservative)  &
            !                 - faceR(elemI(ii,ei_Mface_dL),fr_Flowrate_Conservative)  &
            !                 ) * setting%Time%Hydraulics%Dt
                
            !         if (.not. elemYN(ii,eYN_isZeroDepth)) then
            !             ConsSum = ConsSum + (elemR(ii,er_Volume) - elemR(ii,er_Volume_N0)) &
            !                     -(   faceR(elemI(ii,ei_Mface_uL),fr_Flowrate_Conservative)  &
            !                     - faceR(elemI(ii,ei_Mface_dL),fr_Flowrate_Conservative)  &
            !                     ) * setting%Time%Hydraulics%Dt
            !         end if
            !     end do
            !     print *, 'junctions'
            !     print *, 51, (elemR(51,er_Volume) - elemR(51,er_Volume_N0)) &
            !                 -(   faceR(elemI(52,ei_Mface_uL),fr_Flowrate_Conservative)  &
            !                    - faceR(elemI(53,ei_Mface_dL),fr_Flowrate_Conservative)  &
            !                 ) * setting%Time%Hydraulics%Dt
            !     if (elemR(51,er_Depth) > setting%ZeroValue%Depth) then            
            !         ConsSum = ConsSum + (elemR(51,er_Volume) - elemR(51,er_Volume_N0)) &
            !             -(   faceR(elemI(52,ei_Mface_uL),fr_Flowrate_Conservative)  &
            !                - faceR(elemI(53,ei_Mface_dL),fr_Flowrate_Conservative)  &
            !             ) * setting%Time%Hydraulics%Dt
            !     end if

                ! print *, 62, (elemR(62,er_Volume) - elemR(62,er_Volume_N0)) &
                !             -(   faceR(61,fr_Flowrate_Conservative)  &
                !             - faceR(62,fr_Flowrate_Conservative)  &
                !             ) * setting%Time%Hydraulics%Dt   
                ! if (elemR(62,er_Depth) > setting%ZeroValue%Depth) then            
                !     ConsSum = ConsSum +   (elemR(62,er_Volume) - elemR(62,er_Volume_N0)) &
                !         -(   faceR(61,fr_Flowrate_Conservative)  &
                !         - faceR(62,fr_Flowrate_Conservative)  &
                !         ) * setting%Time%Hydraulics%Dt 
                ! end if          
                                        
                ! print *, 'after weir'
                ! do ii=62,111
                !     print *, ii, (elemR(ii,er_Volume) - elemR(ii,er_Volume_N0)) &
                !                 -(   faceR(elemI(ii,ei_Mface_uL),fr_Flowrate_Conservative)  &
                !                 - faceR(elemI(ii,ei_Mface_dL),fr_Flowrate_Conservative)  &
                !                 ) * setting%Time%Hydraulics%Dt
                !     if (.not. elemYN(ii,eYN_isZeroDepth)) then            
                !     ConsSum = ConsSum + (elemR(ii,er_Volume) - elemR(ii,er_Volume_N0)) &
                !         -(   faceR(elemI(ii,ei_Mface_uL),fr_Flowrate_Conservative)  &
                !         - faceR(elemI(ii,ei_Mface_dL),fr_Flowrate_Conservative)  &
                !         ) * setting%Time%Hydraulics%Dt
                !     end if
                ! end do  
                
                ! print *, ' '
                ! print *, 'CONSERVATION SUM ',ConsSum
                ! print *, ' '

                ! if (abs(ConsSum) > 1e-13) then 
                !     print *, 'Stopping due to conservation'
                !     stop 5597822
                ! end if
            !end if

            ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- TTT  after diagnostic_toplevel')

        ! !% --- ensure mass conservation at junctions
        ! call junction_mass_conservation (ep_JM, istep)
        ! call util_crashstop(6209845)

        !     ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('-------UUU after junction_mass_conservation')

        ! !% --- fix flowrates in diagnostic elements adjacent to JB
        ! call diagnostic_fix_JB_adjacent ()
        ! call util_crashstop(998734)

        !     ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('-------UUU01 after diagnostic_fix_JB_adjacent')




        !% --- RK2 solution step -- check culverts
        !%     HACK: Can this be included in diagnostic?
        ! call culvert_toplevel()     NEED TO REWORK  20230317
        ! call util_crashstop(669743)

            !!! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- UUU02  after culvert_toplevel')

        !% --- RK2 solution step  -- make ad hoc adjustments
        call adjust_Vfilter_CC () !% this is useful in lateral flow induced oscillations
        call util_crashstop(13987)

            ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- VVV after adjust_Vfilter_CC')

        !% --- check zero and small depths that might be changed by Vfilter
        call adjust_element_toplevel (CC)

             ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- WWW after zero/small adjusts')

        !% --- accumulate the volume overflow
        elemR(:,er_VolumeOverFlowTotal)                                                    &
            = elemR(:,er_VolumeOverFlowTotal)         + elemR(:,er_VolumeOverFlow)

        elemR(:,er_VolumeArtificialInflowTotal)                                            &
             = elemR(:,er_VolumeArtificialInflowTotal)+ elemR(:,er_VolumeArtificialInflow)
! 

            ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- ZZZ end of RK2')

    end subroutine rk2_toplevel_ETM_4
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine rk2_toplevel_ETM_3 ()
        ! !%------------------------------------------------------------------
        ! !% Description:
        ! !% single RK2 step for explicit time advance of SVE
        ! !%------------------------------------------------------------------
        ! !% Declarations:
        !     integer :: istep, ii
        !     integer :: whichTM
        !     character(64) :: subroutine_name = 'rk2_toplevel_ETM_NEW'
        ! !%------------------------------------------------------------------
        ! !% Preliminaries
        !     !% --- set the time-marching type for this routine
        !     whichTM = ETM
        !     !% --- reset the overflow counter for this time level
        !     elemR(:,er_VolumeOverFlow) = zeroR     
        !     elemR(:,er_VolumeArtificialInflow) = zeroR   
        ! !%-----------------------------------------------------------------
        ! !% Aliases
        ! !%-----------------------------------------------------------------
        ! !% --- first step of the RK2
        ! istep = oneI

        !     print *, ' '
        !      !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('======= AAA  start of RK2 ==============================')

        ! !% --- NOTE: Dynamic manning's n not included in this routine.

        ! !% --- compute Force Main Manning's N for force main elements
        ! !%     Note this is only computed at the start of the RK2 and the
        ! !%     value is held constant throughout the step.
        ! !%     20230128 -- modified pack so only CC elements, which should not matter
        ! if (setting%Solver%ForceMain%AllowForceMainTF) call forcemain_ManningsN ()   

        ! !% --- store the junction face terms for time n
        ! call junction_face_terms (fp_JB_IorS, istep)    

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- BBB  after junction_face_terms')

        ! !% --- RK2 solution step -- single time advance step
        ! !%     CC advanced for continuity and momentum
        ! call rk2_step_ETM_CC (istep)   
        ! call util_crashstop(7399287)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- CCC  after rk2_step_ETM_cc')
        
        ! !% --- RK2 solution step -- update all CC aux variables
        ! !%     Note, these updates CANNOT depend on face values
        ! call update_auxiliary_variables_CC (whichTM)
        ! call util_crashstop(29687984)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- DDD  after update_auxiliary_variables_CC')

        ! !% --- identify zero depths (.true. is zero depth)
        ! call adjust_zero_or_small_depth_identify_NEW(CC,.true.)
        !     ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- DDD01  after adjusts')
        ! !% --- identify small depths (.false. is small depth)
        ! call adjust_zero_or_small_depth_identify_NEW(CC,.false.)
        !     ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- DDD02  after adjusts')
        ! !% --- create packed arrays
        ! call pack_small_and_zero_depth_elements (whichTM, CC)
        !     ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- DDD03  after adjusts')
        ! !% --- set the zero depths
        ! call adjust_zerodepth_element_values (whichTM, CC) 
        !     ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- DDD04  after adjusts')
        ! !% --- faces that have zero depth on one side or another
        ! call pack_CC_zeroDepth_interior_faces (fp_notJB_interior)
        !     ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- DDD05  after adjusts')
        ! !% --- apply limiters to fluxes and velocity
        ! !%     (.false. so that smalldepth fluxes are not set to zero)
        ! call adjust_smalldepth_element_fluxes_CC (whichTM, .false.)
        !     ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- DDD06  after adjusts')
        ! call adjust_limit_velocity_max_CC (whichTM) 

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- EEE  after adjusts')

        ! !% HACK --- rethink this -- should interp weights for junction Q always be to outside element?
        ! !% --- adjust interpweights to force adjacent element Q, psiL, and energyhead
        ! !%     onto the JB face
        ! call junction_force_Qinterpweights (ep_JM_ETM)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- FFF  after junction_force_Qinterpweights')

        ! !% --- face interpolation without JB
        ! sync all
        ! call face_interpolation(fp_notJB_interior,whichTM,.false.,.false.,.false.)
        ! call util_crashstop(72112342)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- GGG  after face_interpolation')

        ! !% --- flux adjustments (.true. so that conservative fluxes are altered)
        ! call adjust_smalldepth_face_fluxes_CC (.true.)
        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- GGG01  after face_adjustment')
        ! call adjust_zerodepth_face_fluxes_CC  (.true.)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- GGG02  after face_adjustment')

        ! !% --- RK2 solution step -- compute implicit junction
        ! call junction_toplevel_3 (whichTM, istep)
        ! call util_crashstop (11284298)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- HHH  after junction_toplevel')

        ! !% --- set JM and JB variables (resets the Q interpweights)
        ! call update_auxiliary_variables_JMJB(whichTM)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- III  after update_aux...JM')

        ! !% --- identify zero depths
        ! call adjust_zero_or_small_depth_identify_NEW(JB,.true.)
        !     ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('-----III01  after zero depth...JB')
        ! call adjust_zero_or_small_depth_identify_NEW(JM,.true.)
        !     ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- III02  after zero depth...JM')
        ! !% --- identify small depths
        ! call adjust_zero_or_small_depth_identify_NEW(JB,.false.)
        !     ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- III03  after small depth...JB')
        ! call adjust_zero_or_small_depth_identify_NEW(JM,.false.)
        !     ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- III4  after small depth...JM')
        ! !% --- packing
        ! call pack_CC_zeroDepth_interior_faces (fp_notJB_interior)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- JJJ  after adjust...JB JM')

        ! !% --- interpolation to get the non-flow values on faces of junction
        ! call face_interpolation(fp_JB_IorS,whichTM,.true.,.true.,.true.)  
        
        !          !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- KKK  after JB face interpolation ==========================================')

        ! !% --- update diagnostic elements and faces where NOT JB (ensures mass conservation)
        ! call diagnostic_toplevel (ep_Diag_notJBadjacent,fp_Diag_NotJB_interior,istep)  
        ! call util_crashstop(982373)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('-------LLL after diagnostic_toplevel')

        ! !% --- RK2 solution step -- check culverts
        ! !%     HACK: Can this be included in diagnostic?
        ! !%     NEED CHECKING FOR JUNCTIONS
        ! ! call culvert_toplevel() NEED TO REWORK  20230317
        ! ! call util_crashstop(669753)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- MMM  after culvert_toplevel')

        ! !% --- RK2 solution step  -- make ad hoc adjustments
        ! call adjust_Vfilter_CC (whichTM) !% this is useful in lateral flow induced oscillations
        ! call util_crashstop(139871)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- NNN  after adjust_Vfilter_CC')
        
        ! !% -- the conservative fluxes from N to N_1 are the values just before the second RK2 step
        ! call rk2_store_conservative_fluxes (whichTM)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- OOO  after rk2_store_conservative_fluxes')

        ! !% --- reset the overflow counter (we only save conservation in the 2nd step)
        ! elemR(:,er_VolumeOverFlow) = zeroR
        ! elemR(:,er_VolumeArtificialInflow) = zeroR

        ! !% --------------------------------------------------------------------------
        ! !% --------------------------------------------------------------------------
        ! !% --- RK2 solution step -- RK2 second step for ETM 
        ! !%     CC advanced for continuity and momentum
        ! istep = twoI

        ! !% --- store the junction face terms for time n
        ! call junction_face_terms (fp_JB_IorS, istep)    

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- PPP  after junction_face_terms')

        ! call rk2_step_ETM_CC (istep)
        ! call util_crashstop(3298744)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- QQQ  after rk2_step_ETM_cc')

        ! !% --- update all CC aux variables
        ! !%     Note, these updates CANNOT depend on face values
        ! call update_auxiliary_variables_CC (whichTM)
        ! call util_crashstop(2969983)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- RRR  after update_auxiliary_variables_CC')

        ! !% --- identify zero depths
        ! call adjust_zero_or_small_depth_identify_NEW(CC,.true.)
        ! !% --- identify small depths
        ! call adjust_zero_or_small_depth_identify_NEW(CC,.false.)
        ! !% --- create packed arrays
        ! call pack_small_and_zero_depth_elements (whichTM, CC)
        ! !% --- set the minimum geometry values where depth is zero depth
        ! call adjust_zerodepth_element_values (whichTM, CC) 
        ! !% --- faces that have zero depth on one side or another
        ! call pack_CC_zeroDepth_interior_faces (fp_notJB_interior)
        ! !% --- apply limiters to fluxes and velocity
        ! !%     .false. so that smalldepth fluxes are not set to zero
        ! call adjust_smalldepth_element_fluxes_CC (whichTM, .false.)
        ! call adjust_limit_velocity_max_CC (whichTM) 

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- SSS  after adjusts')

        ! !% --- update the JB faces to force adjacent element Q
        ! call junction_force_Qinterpweights (ep_JM_ETM)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- SSS01  after junction_force_Qinterpweights')

        ! !% --- face interpolation without JB
        ! sync all
        ! call face_interpolation(fp_notJB_interior,whichTM,.false.,.false.,.false.)
        ! call util_crashstop(72129873)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- TTT  after face_interpolation')

        ! !% --- flux adjusments with .false. so that conservative fluxes are not altered
        ! call adjust_smalldepth_face_fluxes_CC (whichTM,.false.)
        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- TTT01  after adjust smalldepth face flux')
        ! call adjust_zerodepth_face_fluxes_CC  (whichTM,.false.)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- TTT02  after adjust zerodepth face flux')

        ! !% --- RK2 solution step -- compute implicit junction
        ! call junction_toplevel_3 (whichTM, istep)
        ! call util_crashstop (112873)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- UUU  after junction_toplevel')

        ! !% --- set JM and JB variables (resets the Q interpweights)
        ! call update_auxiliary_variables_JMJB (whichTM)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- VVV after update_aux...JM')

        ! !% --- identify zero depths
        ! call adjust_zero_or_small_depth_identify_NEW(JB,.true.)
        ! call adjust_zero_or_small_depth_identify_NEW(JM,.true.)
        ! !% --- identify small depths
        ! call adjust_zero_or_small_depth_identify_NEW(JB,.false.)
        ! call adjust_zero_or_small_depth_identify_NEW(JM,.false.)
        ! call pack_CC_zeroDepth_interior_faces (fp_JB_IorS)

        !      !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- VVV01 after adjust...JM')

        ! !% --- interpolation to get the non-flow values on faces of junction
        ! !%     skipQ = true
        ! call face_interpolation(fp_JB_IorS,whichTM,.true.,.false.,.true.)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- WWW after face interp for JB')

        ! !% --- update diagnostic elements and faces
        !     !% -- IMPORTANT -- if diagnostic changes a face flowrate for a junction we need to
        ! !%    distribute the new volume residual over the junction branches WITHOUT changing
        ! !%    the head -- essentially an ad hoc fix to ensure mass conservation
        ! !%    There is likely a similar problem with culverts below.
        ! call diagnostic_toplevel (ep_Diag, fp_Diag_IorS, istep)
        ! call util_crashstop(982332)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- XXX  after diagnostic_toplevel')

        ! !% --- RK2 solution step -- check culverts
        ! !%     HACK: Can this be included in diagnostic?
        ! ! call culvert_toplevel()     NEED TO REWORK  20230317
        ! ! call util_crashstop(669743)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- YYY  after culvert_toplevel')

        ! !% --- RK2 solution step  -- make ad hoc adjustments
        ! call adjust_Vfilter_CC (whichTM) !% this is useful in lateral flow induced oscillations
        ! call util_crashstop(13987)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- ZZZ after adjust_Vfilter_CC')

        ! !% --- check zero and small depths that might be changed by Vfilter
        ! !% --- identify zero depths
        ! call adjust_zero_or_small_depth_identify_NEW(CC,.true.)
        ! !% --- identify small depths
        ! call adjust_zero_or_small_depth_identify_NEW(CC,.false.)
        ! !% --- create packed arrays
        ! call pack_small_and_zero_depth_elements (whichTM, CC)
        ! !% --- set the minimum geometry values where depth is zero depth
        ! call adjust_zerodepth_element_values (whichTM, CC) 
        ! !% --- HACK: NOT SURE IF THIS IS NEEDED FOR IMPLICIT APPROACH?
        ! call pack_CC_zeroDepth_interior_faces (fp_noBC_IorS)
        ! !% --- apply limiters to fluxes and velocity
        ! !%     (.false. so that smalldepth fluxes are not set to zero)
        ! call adjust_smalldepth_element_fluxes_CC (whichTM, .false.)
        ! call adjust_limit_velocity_max_CC (whichTM) 

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- ZZZ01 after zero/small adjusts')

        ! !% --- accumulate the volume overflow
        ! elemR(:,er_VolumeOverFlowTotal) = elemR(:,er_VolumeOverFlowTotal) + elemR(:,er_VolumeOverFlow)
        ! elemR(:,er_VolumeArtificialInflowTotal) = elemR(:,er_VolumeArtificialInflowTotal) + elemR(:,er_VolumeArtificialInflow)

        !    !  !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- ZZZ02 end of RK2')

    end subroutine rk2_toplevel_ETM_3
!%
!%==========================================================================
!%==========================================================================
!%    
    subroutine rk2_toplevel_ETM_2 ()
        ! !%------------------------------------------------------------------
        ! !% Description:
        ! !% single RK2 step for explicit time advance of SVE
        ! !%------------------------------------------------------------------
        ! !% Declarations:
        !     integer :: istep, ii
        !     integer :: whichTM
        !     character(64) :: subroutine_name = 'rk2_toplevel_ETM_NEW'
        ! !%------------------------------------------------------------------
        ! !% Preliminaries
        !     !% --- set the time-marching type for this routine
        !     whichTM = ETM
        !     !% --- reset the overflow counter for this time level
        !     elemR(:,er_VolumeOverFlow) = zeroR     
        !     elemR(:,er_VolumeArtificialInflow) = zeroR   
        ! !%-----------------------------------------------------------------
        ! !% Aliases
        ! !%-----------------------------------------------------------------

        !         ! print *, ' '
        !         !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('======= AAA  start of RK2 ==============================')

        ! !% --- NOTE: Dynamic manning's n not included in this routine.

        ! !% --- compute Force Main Manning's N for force main elements
        ! !%     Note this is only computed at the start of the RK2 and the
        ! !%     value is held constant throughout the step.
        ! !%     20230128 -- modified pack so only CC elements, which should not matter
        ! if (setting%Solver%ForceMain%AllowForceMainTF) call forcemain_ManningsN ()   

        ! !% --- RK2 solution step -- single time advance step
        ! !%     CC advanced for continuity and momentum
        ! istep=1
        ! call rk2_step_ETM_CC (istep)   
        ! call util_crashstop(739874)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- BBB  after rk2_step_ETM_cc')
        
        ! !% --- RK2 solution step -- update all CC aux variables
        ! !%     Note, these updates CANNOT depend on face values
        ! call update_auxiliary_variables_CC (whichTM)
        ! call util_crashstop(2968722)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- CCC  after update_auxiliary_variables_CC')

        ! !% --- identify zero depths (.true. is zero depth)
        ! call adjust_zero_or_small_depth_identify_NEW(CC,.true.)
        !     ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- CCC01  after adjusts')
        ! !% --- identify small depths (.false. is small depth)
        ! call adjust_zero_or_small_depth_identify_NEW(CC,.false.)
        !     ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- CCC02  after adjusts')
        ! !% --- create packed arrays
        ! call pack_small_and_zero_depth_elements (whichTM, CC)
        !     ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- CCC03  after adjusts')
        ! !% --- set the zero depths
        ! call adjust_zerodepth_element_values (whichTM, CC) 
        !     ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- CCC04  after adjusts')
        ! !% --- faces that have zero depth on one side or another
        ! call pack_CC_zeroDepth_interior_faces (fp_noBC_IorS)
        !     ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- CCC05  after adjusts')
        ! !% --- apply limiters to fluxes and velocity
        ! !%     (.false. so that smalldepth fluxes are not set to zero)
        ! call adjust_smalldepth_element_fluxes_CC (whichTM, .false.)
        !     ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- CCC06  after adjusts')
        ! call adjust_limit_velocity_max_CC (whichTM) 

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- DDD  after adjusts')


        ! !% --- update the psi CC elements adjacent to junction
        ! !%     must be done after zerodepth and smalldepth adjustments
        ! call update_element_psi_CC (ep_CC_UpstreamOfJunction)
        ! call update_element_psi_CC (ep_CC_DownstreamOfJunction)
        

        ! !% --- update the energy head on CC elements
        ! !%     must be done after zeredepth and smalldepth adjustments
        ! call update_element_energyHead_CC (ep_CC_ETM)

        ! !% --- adjust interpweights to force adjacent element Q, psiL, and energyhead
        ! !%     onto the JB face
        ! call junction_force_Qinterpweights (ep_JM_ETM)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- EEE  after junction_force_Qinterpweights')

        ! !% --- RK2 solution step  -- all face interpolation
        ! !%     junction JB faces get adjacent element Q values
        ! sync all
        ! call face_interpolation(fp_noBC_IorS,whichTM,.false.,.false.,.false.)
        ! call util_crashstop(72129873)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- FFF  after face_interpolation')

        ! !% --- flux adjustments (.true. so that conservative fluxes are altered)
        ! call adjust_smalldepth_face_fluxes_CC (whichTM,.true.)
        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- FFF01  after face_adjustment')
        ! call adjust_zerodepth_face_fluxes_CC  (whichTM,.true.)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- FFF02  after face_adjustment')

        ! !% --- RK2 solution step -- compute implicit junction
        ! call junction_toplevel (whichTM, istep)
        ! call util_crashstop (112875)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- GGG  after junction_toplevel')

        ! !% --- set JM and JB variables (resets the Q interpweights)
        ! call update_auxiliary_variables_JMJB (whichTM)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- HHH  after update_aux...JM')

        ! !% --- identify zero depths
        ! call adjust_zero_or_small_depth_identify_NEW(JB,.true.)
        !     ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('-------HHH01  after zero depth...JB')
        ! call adjust_zero_or_small_depth_identify_NEW(JM,.true.)
        !     ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- HHH02  after zero depth...JM')
        ! !% --- identify small depths
        ! call adjust_zero_or_small_depth_identify_NEW(JB,.false.)
        !     ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- HHH03  after small depth...JB')
        ! call adjust_zero_or_small_depth_identify_NEW(JM,.false.)
        !     ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- HHH04  after small depth...JM')
        ! !% --- packing
        ! call pack_CC_zeroDepth_interior_faces (fp_noBC_IorS)

        !     ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- HHH05  after adjust...JB JM')

        ! !% --- interpolation to get the non-flow values on faces of junction
        ! !%     skipQ = true
        !     ! print *, ' '
        !     ! print *, 'fp_JB_IorS', npack_faceP(fp_JB_IorS)
        !     ! print *, faceP(1:npack_faceP(fp_JB_IorS),fp_JB_IorS)
        !     ! print *, ' '

        ! call face_interpolation(fp_JB_IorS,whichTM,.true.,.true.,.true.)    

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- III  after JB face interpolation ==========================================')

        ! !% --- RK2 solution step  -- update diagnostic elements and faces
        ! !%     (.true. as this is RK first step)
        ! !%NOTE: WE DO NOT WANT THE DIAGNOSTIC TOPLEVEL CHANGING THE Q FLUX AT THE
        ! !% FACE WITH A JB --- rethink interpolation here. 
        ! !% -- IMPORTANT -- if diagnostic changes a face flowrate for a junction we need to
        ! !%    distribute the new volume residual over the junction branches WITHOUT changing
        ! !%    the head -- essentially an ad hoc fix to ensure mass conservation
        ! !%    There is likely a similar problem with culverts below.
        ! !call diagnostic_toplevel (.true.)  TEMPORARY CUTOUT 20230210 CAUSING PROBLEMS ON JUNCTIONS DOWNSTREAM FACES
        ! !call util_crashstop(982332)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- JJJ  after diagnostic_toplevel')

        ! !% --- RK2 solution step -- check culverts
        ! !%     HACK: Can this be included in diagnostic?
        ! call culvert_toplevel()
        ! call util_crashstop(669743)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- KKK  after culvert_toplevel')

        ! !% --- RK2 solution step  -- make ad hoc adjustments
        ! call adjust_Vfilter_CC (whichTM) !% this is useful in lateral flow induced oscillations
        ! call util_crashstop(13987)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- LLL  after adjust_Vfilter_CC')
        
        ! !% -- the conservative fluxes from N to N_1 are the values just before the second RK2 step
        ! call rk2_store_conservative_fluxes (whichTM)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- MMM  after rk2_store_conservative_fluxes')

        ! !% --- reset the overflow counter (we only save conservation in the 2nd step)
        ! elemR(:,er_VolumeOverFlow) = zeroR
        ! elemR(:,er_VolumeArtificialInflow) = zeroR

        ! !% --------------------------------------------------------------------------
        ! !% --------------------------------------------------------------------------
        ! !% --- RK2 solution step -- RK2 second step for ETM 
        ! !%     CC advanced for continuity and momentum
        ! istep=2
        ! call rk2_step_ETM_CC (istep)
        ! call util_crashstop(3298744)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- NNN  after rk2_step_ETM_cc')

        ! !% --- RK2 solution step -- update all CC aux variables
        ! !%     Note, these updates CANNOT depend on face values
        ! call update_auxiliary_variables_CC (whichTM)
        ! call util_crashstop(2968722)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- OOO  after update_auxiliary_variables_CC')

        ! !% --- identify zero depths
        ! call adjust_zero_or_small_depth_identify_NEW(CC,.true.)
        ! !% --- identify small depths
        ! call adjust_zero_or_small_depth_identify_NEW(CC,.false.)
        ! !% --- create packed arrays
        ! call pack_small_and_zero_depth_elements (whichTM, CC)
        ! !% --- set the minimum geometry values where depth is zero depth
        ! call adjust_zerodepth_element_values (whichTM, CC) 
        ! !% --- faces that have zero depth on one side or another
        ! call pack_CC_zeroDepth_interior_faces (fp_noBC_IorS)
        ! !% --- apply limiters to fluxes and velocity
        ! !%     .false. so that smalldepth fluxes are not set to zero
        ! call adjust_smalldepth_element_fluxes_CC (whichTM, .false.)
        ! call adjust_limit_velocity_max_CC (whichTM) 

        ! !% --- update the psi CC elements adjacent to junction
        ! !%     must be done after zerodepth and smalldepth adjustments
        ! call update_element_psi_CC (ep_CC_DownstreamOfJunction)
        ! call update_element_psi_CC (ep_CC_UpstreamOfJunction)

        ! !% --- update the energy head on CC elements
        ! !%     must be done after zeredepth and smalldepth adjustments
        ! call update_element_energyHead_CC (ep_CC_ETM)

        ! !% --- update the JB faces to force adjacent element Q, psiL, and energyhead
        ! call junction_force_Qinterpweights (ep_JM_ETM)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- PPP  after adjusts')

        ! !% --- RK2 solution step  -- all face interpolation
        ! !%     junction JB faces get adjacent element Q values
        ! sync all
        ! call face_interpolation(fp_noBC_IorS,whichTM,.false.,.false.,.false.)
        ! call util_crashstop(72129873)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- QQQ  after face_interpolation')

        ! !% --- flux adjusments with .false. so that conservative fluxes are not altered
        ! call adjust_smalldepth_face_fluxes_CC (whichTM,.false.)
        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- QQQ01  after adjust smalldepth face flux')
        ! call adjust_zerodepth_face_fluxes_CC  (whichTM,.false.)
        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- QQQ02  after adjust zerodepth face flux')

            

        ! !% --- RK2 solution step -- compute implicit junction
        ! call junction_toplevel (whichTM, istep)
        ! call util_crashstop (112873)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- RRR  after junction_toplevel')

        ! !% --- set JM and JB variables (resets the Q interpweights)
        ! call update_auxiliary_variables_JMJB (whichTM)

        ! !% --- identify zero depths
        ! call adjust_zero_or_small_depth_identify_NEW(JB,.true.)
        ! call adjust_zero_or_small_depth_identify_NEW(JM,.true.)
        ! !% --- identify small depths
        ! call adjust_zero_or_small_depth_identify_NEW(JB,.false.)
        ! call adjust_zero_or_small_depth_identify_NEW(JM,.false.)
        ! call pack_CC_zeroDepth_interior_faces (fp_noBC_IorS)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- SSS after update_aux...JM')

        ! !% --- interpolation to get the non-flow values on faces of junction
        ! !%     skipQ = true
        ! call face_interpolation(fp_JB_IorS,whichTM,.true.,.false.,.true.)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- TTT after face interp for JB')

        ! !% --- RK2 solution step  -- update diagnostic elements and faces
        ! !%     (.false. as this is RK second step)
        ! !%NOTE: WE DO NOT WANT THE DIAGNOSTIC TOPLEVEL CHANGING THE Q FLUX AT THE
        ! !% FACE WITH A JB --- rethink interpolation here.
        ! ! call diagnostic_toplevel (.false.)  TEMPORARY CUTOUT 20230210 CAUSING PROBLEMS ON JUNCTIONS DOWNSTREAM FACES
        ! ! call util_crashstop(982332)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- UUU  after diagnostic_toplevel')

        ! !% --- RK2 solution step -- check culverts
        ! !%     HACK: Can this be included in diagnostic?
        ! call culvert_toplevel()
        ! call util_crashstop(669743)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- VVV  after culvert_toplevel')

        ! !% --- RK2 solution step  -- make ad hoc adjustments
        ! call adjust_Vfilter_CC (whichTM) !% this is useful in lateral flow induced oscillations
        ! call util_crashstop(13987)

        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- WWW  after adjust_Vfilter_CC')

        ! !% --- check zero and small depths that might be changed by Vfilter
        ! !% --- identify zero depths
        ! call adjust_zero_or_small_depth_identify_NEW(CC,.true.)
        ! !% --- identify small depths
        ! call adjust_zero_or_small_depth_identify_NEW(CC,.false.)
        ! !% --- create packed arrays
        ! call pack_small_and_zero_depth_elements (whichTM, CC)
        ! !% --- set the minimum geometry values where depth is zero depth
        ! call adjust_zerodepth_element_values (whichTM, CC) 
        ! !% --- HACK: NOT SURE IF THIS IS NEEDED FOR IMPLICIT APPROACH?
        ! call pack_CC_zeroDepth_interior_faces (fp_noBC_IorS)
        ! !% --- apply limiters to fluxes and velocity
        ! !%     (.false. so that smalldepth fluxes are not set to zero)
        ! call adjust_smalldepth_element_fluxes_CC (whichTM, .false.)
        ! call adjust_limit_velocity_max_CC (whichTM) 

        ! !% --- accumulate the volume overflow
        ! elemR(:,er_VolumeOverFlowTotal) = elemR(:,er_VolumeOverFlowTotal) + elemR(:,er_VolumeOverFlow)
        ! elemR(:,er_VolumeArtificialInflowTotal) = elemR(:,er_VolumeArtificialInflowTotal) + elemR(:,er_VolumeArtificialInflow)
        


        !     !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('------- ZZZ end of RK2')

    end subroutine rk2_toplevel_ETM_2        
!%
!%==========================================================================
!%==========================================================================
!%    
    subroutine rk2_toplevel_ETM ()
        ! !%------------------------------------------------------------------
        ! !% Description:
        ! !% single RK2 step for explicit time advance of SVE
        ! !%------------------------------------------------------------------
        ! !% Declarations:
        !     integer :: istep, ii
        !     integer :: whichTM
        !     character(64) :: subroutine_name = 'rk2_toplevel_ETM'
        ! !%------------------------------------------------------------------
        ! !% Preliminaries
        !     if (setting%Debug%File%runge_kutta2) &
        !         write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !     !% --- set the time-marching type for this routine
        !     whichTM = ETM
        !     !% --- reset the overflow counter for this time level
        !     elemR(:,er_VolumeOverFlow) = zeroR     
        !     elemR(:,er_VolumeArtificialInflow) = zeroR   
        ! !%-----------------------------------------------------------------
        ! !% Aliases
        ! !%-----------------------------------------------------------------

        !     ! do ii=1,N_elem(1)
        !     !     print *, ii, elemR(ii,er_Depth)
        !     ! end do
        !     ! print *, ' '
        !     ! do ii=1, N_face(1)
        !     !     print *, ii, faceR(ii,fr_Depth_u), faceR(ii,fr_Depth_d)
        !     ! end do
        !     ! stop 498733

        !     !print *, ' '
        !     ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('======= AAA  start of RK2 ==============================')
        !     ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CheckIsNan ()

        !     !stop 409879
            

        ! !% --- compute the dynamic mannings N (DISABLED AS OF 20220817 brh)
        ! if (setting%Solver%ManningsN%useDynamicManningsN) then
        !     call rk2_dynamic_ManningsN (whichTM)
        ! else
        !     !% --- arguably not needed, but here to prevent bugs
        !     !%     if dynamic Mannings N is used in place of standard Mannings N
        !     elemR(:,er_ManningsN_Dynamic) = elemR(:,er_ManningsN)
        ! end if

        
        ! !% --- compute Force Main Manning's N for force main elements
        ! if (setting%Solver%ForceMain%AllowForceMainTF) call forcemain_ManningsN ()

        ! !% --- RK2 solution step -- single time advance step
        ! !%     CC advanced for continuity and momentum
        ! !%     JM advanced for continuity
        ! !%     JB determined from adjacent values.
        ! istep=1
        ! call rk2_step_ETM(istep)

        !     ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('BBB after volume/momentum step 1---------------------------')
   
        ! !% --- RK2 solution step -- update all non-diagnostic aux variables
        ! !%     Note, these updates CANNOT depend on face values
        ! call update_auxiliary_variables(whichTM)

        !     ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('CCC  after update aux step 1-----------------------')

        ! !% --- set the flagged zero and small depth cells (allow depth to change)
        ! !%     This does not reset the zero/small depth packing as we allow negative
        ! !%     volume values in first RK2 step (but not in final step)
        ! !%     TESTING MAKING THIS RESET 20221227 to solve problem with fr_Flowrate_Conservative
        ! call adjust_zero_and_small_depth_elem (whichTM, .true., .false.)
        ! call util_crashstop(340927)

        !     ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('DDD  after adjust zero/small elem-----------------')

        ! !% --- RK2 solution step  -- all face interpolation
        ! sync all
        ! call face_interpolation(fp_noBC_IorS,whichTM,.false.,.false.,.false.)

        !     ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('EEE  after face interpolation step 1---------------')

        !     ! if (setting%Time%Step == 15) then
        !     !     print *, ' '
        !     !     stop 733723
        !     !  end if

        ! !% --- set the zero and small depth fluxes
        ! !%     This resets the faces that are zero
        ! call adjust_zero_and_small_depth_face (whichTM, .true.)
        ! call util_crashstop(440223)

        !     ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('FFF  after zero/small face step 1-----------------')

        ! !% --- RK2 solution step  -- update diagnostic elements and faces
        ! call diagnostic_toplevel (.true.)
        ! call util_crashstop(402873)

        !     ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('GGG  after diagnostic step 1')

        ! !% --- RK2 solution step -- set the JB as a function of the face and JM values
        ! if (setting%Junction%Method .eq. Explicit2) then 
        !     call ll_alternate_JB (whichTM,istep)
        ! end if

        ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('GGG01  after after ll_alternate_JB')

        ! !% --- RK2 solution step -- check culverts
        ! call culvert_toplevel()
        ! call util_crashstop(669743)

        ! !% --- RK2 solution step  -- make ad hoc adjustments
        ! call adjust_Vfilter_CC (whichTM) ! brh20220211 this is useful in lateral flow induced oscillations
        ! call util_crashstop(13987)

        !     ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('HHH  after Vfilter step 1 -------------------------')
        
        ! !% -- the conservative fluxes from N to N_1 are the values just before the second RK2 step
        ! call rk2_store_conservative_fluxes (whichTM)

        !      ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('III  after consQ store step 1 ----------------------')

        ! !% --- reset the overflow counter (we only save conservation in the 2nd step)
        ! elemR(:,er_VolumeOverFlow) = zeroR
        ! elemR(:,er_VolumeArtificialInflow) = zeroR

        ! !% --------------------------------------------------------------------------
        ! !% --- RK2 solution step -- RK2 second step for ETM 
        ! !%     CC advanced for continuity and momentum
        ! !%     JM advanced for continuity
        ! !%     JB determined from adjacent values.

        ! istep=2
        ! call rk2_step_ETM (istep)
        
        !     ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('JJJ  after volume rk2 step 2 -----------------------')

        ! !% --- RK2 solution step -- update non-diagnostic auxiliary variables
        ! !%     Note, these updates CANNOT depend on face values
        ! call update_auxiliary_variables(whichTM)  

        !     ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('KKK  after update aux step 2 --------------------------')

        ! !% --- set the flagged zero and small depth cells (allow depth to change)
        ! !%     This DOES reset the packing 20221227brh
        ! !%     Reset is required so that we don't get negative volumes in final RK2 step
        ! call adjust_zero_and_small_depth_elem (whichTM, .true., .false.)
        ! call util_crashstop(12973)

        !     ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('LLL  after zero/small elem step 2 -------------------')

        ! !% --- RK2 solution step -- update all faces
        ! sync all
        ! call face_interpolation(fp_noBC_IorS,whichTM,.false.,.false.,.false.)

        !     ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('MMM  after face interp step 2 --------------------------')

        ! !% --- set the zero and small depth fluxes
        ! !%     ifixQcons = false to prevent conservation issues (cannot change Qcons after 2nd RK2 step)
        ! call adjust_zero_and_small_depth_face (whichTM, .false.)

        !      ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('NNN  after zero/small face step 2 ---------------------')

        ! !% --- RK2 solution step -- update diagnostic elements and faces
        ! call diagnostic_toplevel (.false.)
        ! call util_crashstop(662398)

        !     ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('OOO  after diagnostic step 2')

        ! !% --- RK2 solution step -- set the JB as a function of the face and JM values
        ! if (setting%Junction%Method .eq. Explicit2) then 
        !     call ll_alternate_JB (whichTM,istep)
        ! end if

        !     ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('OOO_01  after ll_alternae_JB step 2')

        ! !% --- RK2 solution step -- check culverts
        ! call culvert_toplevel()
        ! call util_crashstop(669742)        
        
        ! !% --- RK2 solution step -- make ad hoc adjustments (V filter)
        ! ! print *, 'vfilter before ',elemR(54,er_Flowrate)
        ! call adjust_Vfilter_CC (whichTM)
        ! call util_crashstop(449872)
        ! ! print *, 'vfilter after ',elemR(54,er_Flowrate)

        !     ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('PPP  after Vfilter step 2-----------------------------')

        ! !% --- ensures that the Vfilter hasn't affected the zero/small depth cells        
        ! call adjust_zero_and_small_depth_elem (whichTM, .true., .false.)
        ! call util_crashstop(64987)

        !     ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('QQQ  after zero/small elem step 2 (2nd time)')

        ! !% --- accumulate the volume overflow
        ! elemR(:,er_VolumeOverFlowTotal) = elemR(:,er_VolumeOverFlowTotal) + elemR(:,er_VolumeOverFlow)
        ! elemR(:,er_VolumeArtificialInflowTotal) = elemR(:,er_VolumeArtificialInflowTotal) + elemR(:,er_VolumeArtificialInflow)
        
        !     ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('ZZZ  after accumulate overflow step 2')
        !     ! print *, '==================================================='
        !     ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CheckIsNan ()

        !     ! if (setting%Time%Step == 25) then
        !     !    print *, ' '
        !     !    stop 6987588
        !     ! end if

        ! !%-----------------------------------------------------------------
        ! !% closing
        !     if (setting%Debug%File%runge_kutta2)  &
        !         write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
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
        ! !%-----------------------------------------------------------------------------
        ! !% Description:
        ! !%
        ! !%-----------------------------------------------------------------------------
        ! integer :: istep, faceMaskCol, thisCol
        ! integer, pointer :: Npack
        ! character(64) :: subroutine_name = 'rk2_toplevel_ETMAC'
        ! !%-----------------------------------------------------------------------------
        ! !if (crashYN) return
        ! if (setting%Debug%File%runge_kutta2) &
        !     write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        
        ! !% step 1 -- RK2 step 1 for ETM

        ! istep=1
        ! if (N_etm > 0) then
        !     call rk2_step_ETM (istep)
        ! end if

        ! !% step 2 -- RK2 step 1 for AC
        ! if (N_ac > 0) then
        !     !% step 2(a)
        !     call rk2_step_AC (istep)
        !     if (N_etm > 0) then
        !         !% step 2(b,c) create time n+1(1*) consistency for AC
        !         call rk2_extrapolate_to_fullstep_ETM()
        !     end if
        ! end if

        ! !% step 3 -- all aux variables for non-diagnostic
        ! call update_auxiliary_variables(ALLtm)

        ! !% step 4 -- all face interpolation
        ! sync all
        ! call face_interpolation(fp_noBC_IorS,ALLtm,.false.,.false.,.false.)
        
        ! !% step 5 -- update diagnostic elements and faces
        ! call diagnostic_toplevel (.true.)
        ! call util_crashstop(66234)

        ! !% step X -- make ad hoc adjustments
        ! call adjust_Vfilter_CC (ALLtm)
        ! call util_crashstop(53434)

        ! !% step 6 -- RK2 step 2 for AC
        ! istep=2
        ! if (N_ac > 0) then
        !     !% step 6(a)
        !     call rk2_step_AC (istep)
        !     if (N_etm > 0) then
        !         !% step 6(b)
        !         call rk2_restore_to_midstep_ETM()
        !         !% step 6(c,d)
        !         call rk2_interpolate_to_halfstep_AC()
        !     end if
        !     !% step 6(e)
        !     call update_auxiliary_variables (AC)

        !     !% step 6(f,g) -- update faces for AC elements
        !     !% HACK -- need to sync some of the processors before this, but cannot sync all because some have N_ac <1
        !     call face_interpolation (fp_AC,AC,.false.,.false.,.false.)
            

        ! end if

        ! !% step 7 -- update diagnostic elements and faces
        ! call diagnostic_toplevel(.false.)
        ! call util_crashstop(12293)


        ! !% step 8 -- RK2 step 2 for ETM
        ! if (N_etm > 0) then
        !     !% step 8(a)
        !     call rk2_step_ETM (istep)

        !     if (N_ac > 0) then
        !         !% step 8(b)
        !         call rk2_restore_to_fullstep_AC ()
        !     end if

        !     !% step 8(c)
        !     call update_auxiliary_variables(ETM)

        !     !% step 8(d,e) -- update all faces
        !     !% HACK -- need to sync some of the processors before this, but cannot sync all because some have N_ac <1
        !     call face_interpolation (fp_noBC_IorS,ALLtm,.false.,.false.,.false.)
        ! end if

        ! !% step 9 -- update diagnostic elements and faces
        ! call diagnostic_toplevel (.false.)
        ! call util_crashstop(23422)

        ! !% step X -- make ad hoc adjustments
        ! call adjust_Vfilter_CC (ALLtm)
        ! call util_crashstop(99287)

        ! if (setting%Debug%File%runge_kutta2)  &
        !     write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine rk2_toplevel_ETMAC
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
            
            ! ! ! ! ! ! ! call util_utest_CLprint ('-------  WWW 00  start of rk2 Step CC')

        !% --- CONTINUITY
        thisPackCol => col_elemP(ep_CC_H)
        Npack       => npack_elemP(thisPackCol)
        if (Npack > 0) then
            thisP => elemP(1:Npack,thisPackCol)

            !% --- Compute net flowrates for CC as source termo
            call ll_continuity_netflowrate_CC (er_SourceContinuity, thisPackCol, Npack)


                ! ! ! ! ! ! ! call util_utest_CLprint ('-------  WWW 01  after ll_continuity_netflowrate')

            !% --- Solve for new volume
            call ll_continuity_volume_ETM (er_Volume, thisPackCol, Npack, istep)

                ! ! ! ! ! ! ! call util_utest_CLprint ('-------  WWW 02  after ll_continuity_volume')


            !% --- adjust extremely small volumes that might be been introduced
            call adjust_limit_by_zerovalues &
                (er_Volume, setting%ZeroValue%Volume, thisP, .true.)

                ! ! ! ! ! ! ! call util_utest_CLprint ('-------  WWW 03  after adjust_limit_by_zerovalues')

        end if  



        !% --- MOMENTUM
        thisPackCol => col_elemP(ep_CC_Q)
        Npack       => npack_elemP(thisPackCol)
        if (Npack > 0) then

            !% --- momentum K source terms for different methods for ETM
            call ll_momentum_Ksource_CC (er_Ksource, thisPackCol, Npack)

                !print *, 'Ksource ',elemR(1,er_KSource)

            !% --- Common source for momentum on channels and conduits for ETM
            call ll_momentum_source_CC (er_SourceMomentum, thisPackCol, Npack)

                ! ! ! ! ! ! ! call util_utest_CLprint ('-------  WWW 06  after ll_momentum_source_cc')
 
               ! print *, 'Source ',elemR(1,er_SourceMomentum)

            !% --- Common Gamma for momentum on channels and conduits for  ETM
            !%     Here for all channels and conduits, assuming CM roughness
            call ll_momentum_gammaCM_CC (er_GammaM, thisPackCol, Npack)

                !print *, 'Gamma  ',elemR(1,er_GammaM)

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

                !print *, 'Gamma2 ',elemR(58,er_GammaM)

            !% --- Advance flowrate to n+1/2 for conduits and channels with ETM
            call ll_momentum_solve_CC (er_Velocity, thisPackCol, Npack, thisMethod, istep)

                !print *, 'Vel 1   ',elemR(1,er_Velocity)

                ! ! ! ! ! ! ! call util_utest_CLprint ('-------  WWW 15  after ll_momentum_solve_cc')

            !% --- velocity for ETM time march
            call ll_momentum_velocity_CC (er_Velocity, thisPackCol, Npack)

                !print *, 'Vel 2   ',elemR(1,er_Velocity)

            !% --- prevent backflow through flapgates
            call ll_enforce_flapgate_CC (er_Velocity, thisPackCol, Npack)

                !print *, 'flap   ',elemR(58,er_Velocity)

            !% --- enforce zero velocity on elements that began as ZeroDepth
            call ll_enforce_zerodepth_velocity (er_Velocity, thisPackCol, Npack)

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
        !     ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('before rk2 continuity step etm')

        ! !% perform the continuity step of the rk2 for ETM CC and JM
        ! call rk2_continuity_step_ETM(istep)

        !     ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('after rk2 continuity step etm')

        ! !% only adjust extremely small element volumes that have been introduced
        ! Npack = npack_elemP(ep_CCJM_H_ETM)
        ! if (Npack > 0) then
        !     call adjust_limit_by_zerovalues &
        !         (er_Volume, setting%ZeroValue%Volume/twentyR, elemP(1:Npack,ep_CCJM_H_ETM), .true.)
        ! end if
        !     ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint ('after rk2 call to adjust limit by zero')

        ! !% perform the momentum step of the rk2 for ETM
        ! call rk2_momentum_step_ETM(istep)

            ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint (' after rk2 call to rk2_momentum_step_ETM')

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

                ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint (' start of rk2_momentum_step_ETM')

                ! print *, ' '
                ! print *, 'in rk2_momentum_step at start'
                ! print *, elemR(61,er_Ksource), elemR(61,er_Velocity), elemR(61,er_Flowrate)
   
            !% --- momentum K source terms for different methods for ETM
            call ll_momentum_Ksource_CC (er_Ksource, thisPackCol, Npack)


                ! print *, 'in rk2_momentum_step at A'
                ! print *, elemR(61,er_Ksource), elemR(61,er_Velocity), elemR(61,er_Flowrate)

                ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint (' after rk2 call ll_momentum_Ksource_CC')


            !% --- Common source for momentum on channels and conduits for ETM
            call ll_momentum_source_CC (er_SourceMomentum, thisPackCol, Npack)

                ! print *, 'in rk2_momentum_step at B'
                ! print *, elemR(61,er_SourceMomentum), elemR(61,er_HydRadius), elemR(61,er_ManningsN), elemR(61,er_Flowrate)

                ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint (' after rk2 call ll_momentum_source_CC')

            !% --- Common Gamma for momentum on channels and conduits for  ETM
            !%     Here for all channels and conduits, assuming CM roughness
            call ll_momentum_gammaCM_CC (er_GammaM, thisPackCol, Npack)

                ! print *, 'in rk2_momentum_step at C'
                ! print *, elemR(61,er_GammaM), elemR(61,er_Velocity), elemR(61,er_Flowrate)      
            
                ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint (' after rk2 call ll_momentum_gammaCM_CC')

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

                ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint (' after rk2 call ll_minorloss_friction_gamma_CC')
  
            !% --- Advance flowrate to n+1/2 for conduits and channels with ETM
            call ll_momentum_solve_CC (er_Velocity, thisPackCol, Npack, thisMethod, istep)


                ! print *, 'in rk2_momentum_step at F'
                ! print *, elemR(61,er_GammaM), elemR(61,er_Velocity), elemR(61,er_Flowrate)
            
                ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint (' after rk2 call ll_momentum_solve_CC')


            !% --- velocity for ETM time march
            call ll_momentum_velocity_CC (er_Velocity, thisPackCol, Npack)


                ! print *, 'in rk2_momentum_step at G'
                ! print *, elemR(61,er_GammaM), elemR(61,er_Velocity), elemR(61,er_Flowrate)
                ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint (' after rk2 call ll_momentum_velocity_CC')


            !% --- prevent backflow through flapgates
            call ll_enforce_flapgate_CC (er_Velocity, thisPackCol, Npack)

                ! print *, 'in rk2_momentum_step at H'
                ! print *, elemR(61,er_GammaM), elemR(61,er_Velocity), elemR(61,er_Flowrate)

                ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint (' after rk2 call ll_enforce_flapgate_CC')

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

    
            ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! call util_utest_CLprint (' after ll_flowrate_and_velocity_JB')

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

