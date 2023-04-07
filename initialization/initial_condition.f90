!
!% module initial_condition
!% Provides setup of initial conditions for geometry and dynamics
!
!
!==========================================================================
!
!==========================================================================
!
module initial_condition

    use define_indexes
    use define_keys
    use define_globals
    use define_settings
    use define_xsect_tables
    use pack_mask_arrays
    use boundary_conditions
    use update
    use face
    use forcemain, only : forcemain_ManningsN
    use diagnostic_elements
    use geometry 
    use arch_conduit
    use circular_conduit
    use basket_handle_conduit
    use catenary_conduit
    use egg_shaped_conduit
    use gothic_conduit
    use horiz_ellipse_conduit
    use vert_ellipse_conduit
    use horse_shoe_conduit
    use semi_elliptical_conduit
    use semi_circular_conduit
    use filled_circular_conduit
    use geometry_lowlevel
    use irregular_channel, only: irregular_geometry_from_depth_singular
    ! use rectangular_channel, only: rectangular_area_from_depth
    ! use parabolic_channel, only: parabolic_area_from_depth
    ! use rectangular_conduit, only: rectangular_closed_area_from_depth
    ! use rectangular_triangular_conduit, only: rectangular_triangular_area_from_depth
    ! use trapezoidal_channel, only: trapezoidal_area_from_depth
    ! use triangular_channel, only: triangular_area_from_depth
    use storage_geometry
    use preissmann_slot, only: slot_initialize
    use adjust
    use xsect_tables
    use control_hydraulics, only: control_init_monitoring_and_action_from_EPASWMM
    use interface_, only: interface_get_nodef_attribute
    use utility_profiler
    use utility_allocate
    use utility_deallocate
    use utility_interpolate
    use utility_key_default
    use utility_crash, only: util_crashpoint
    !use utility, only: util_CLprint
    !use utility_unit_testing, only: util_utest_CLprint

    implicit none

    public :: init_IC_toplevel

    private

contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine init_IC_toplevel ()
        !%------------------------------------------------------------------
        !% Description:
        !% set up the initial conditions for all the elements
        !%------------------------------------------------------------------
        !% Declarations:
            integer          :: ii, iblank, whichTM
            integer, pointer :: whichSolver
            integer, allocatable :: tempP(:)  !% for debugging
            integer, pointer :: thisCol, npack, thisP(:)
            character(64)    :: subroutine_name = 'init_IC_toplevel'
        !%-------------------------------------------------------------------
        !% Preliminaries:
            if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
        !% Aliases
            whichSolver => setting%Solver%SolverSelect
            select case (whichSolver)
            case (ETM)
                whichTM = ETM
            case (AC)
                whichTM = AC
                print *, 'CONFIGURATION ERROR: AC solver is not functional'
                print *, 'use setting.Solver.Select = ETM'
                call util_crashpoint(7298744)
            case (ETM_AC)
                print *, 'CONFIGURATION ERROR: AC solver is not functional'
                print *, 'use setting.Solver.Select = ETM'
                call util_crashpoint(7298744)
            case default
                print *, 'CODE ERROR: solver type unknown for # ', whichSolver
                print *, 'which has key ',trim(reverseKey(whichSolver))
                call util_crashpoint(478383)
                !return
            end select      
        !%-------------------------------------------------------------------    
        !% --- default the TimeLastSet to 0.0 (elapsed) for all elements
        elemR(:,er_TimeLastSet) = zeroR

        !% --- initialize all the element Setting as 1
        !%     this is fully open for links, weirs, orifices, outlets, and on for pumps.
        elemR(1:size(elemR,1)-1,er_TargetSetting) = oneR
        elemR(1:size(elemR,1)-1,er_Setting)       = oneR

        !% --- initialize all the minor losses and seepage rates to zero
        elemR(:,er_KJunction_MinorLoss) = zeroR
        elemR(:,er_Kconduit_MinorLoss)  = zeroR
        elemR(:,er_SeepRate)            = zeroR

        !% --- initialize overflow
        elemR(:,er_VolumeOverFlow) = zeroR
        elemR(:,er_VolumeOverFlowTotal) = zeroR

        elemR(:,er_VolumeArtificialInflowTotal) = zeroR

        !% --- initialize barrels
        setting%Output%BarrelsExist = .false. !% will be set to true if barrels > 1 detected
        elemI(:,ei_barrels) = oneR

        !% --- initialize sedmient depths
        !%     Note: as of 20221006 only FilledCircular is allowed to have nonzero sediment depth
        !%     this corresponds to the "yBot" of the Filled Circular cross-section in EPA-SWMM
        elemR(:,er_SedimentDepth) = zeroR

            ! ! call util_utest_CLprint ('initial_condition near start')

        !% --- get data that can be extracted from links
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin init_IC_from_linkdata'
        call init_IC_from_linkdata ()

            ! ! call util_utest_CLprint ('initial_condition after IC_from_linkdata')

        !print *, 'TEST 98743987 ', elemR()

          !  stop 6398743

        sync all
        !% --- set up background geometry for weir, orifice, etc.
        !%     from adjacent elements
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin init_IC_diagnostic_geometry_from_adjacent'
        call init_IC_diagnostic_geometry_from_adjacent (.true.)

            ! ! call util_utest_CLprint ('initial_condition after diagnostic_geometry_from_adjacent')

        !% --- sync after all the link data has been extracted
        !%     the junction branch data is read in from the elemR arrays which
        !%     may need inter-image communication, thus the sync call is needed
        sync all
        
        !% --- get data that can be extracted from nodes
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin init_IC_for_nJM_from_nodedata'
        call init_IC_for_nJm_from_nodedata ()

            ! ! call util_utest_CLprint ('initial_condition afer IC_from_nodedata')

        !% --- second call for diagnostic that was next to JM/JB
        sync all
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin init_IC_diagnostic_geometry_from_adjacent'
        call init_IC_diagnostic_geometry_from_adjacent (.false.)

        !% --- identify special case diagnostic elements that have JB on either side
        !%     these have the face flowrates frozen in the junction residual computation
        call init_IC_diagnostic_JB_bounded ()

        !% --- error checking for nullvalues
        !%     keep for future debugging use.
        ! do ii=1,N_elem(1)
        !     if ((elemI(ii,ei_geometryType) == nullvalueI) .or. &
        !         (elemI(ii,ei_geometryType) == undefinedKey)) then

        !         if (    ((elemI(ii,ei_elementType) .eq. JB) .and.    &
        !                   elemSI(ii,esi_JunctionBranch_Exists)     ) &
        !             .or.                                             &
        !                 (elemI(ii,ei_elementType) .ne. JB) ) then
        !             print *, 'POSSIBLE PROBLEM IN link/node with nullvalue or undefined geometry Type'
        !             print *, 'ii ',ii, elemI(ii,ei_geometryType)
        !             print *, 'link id ',elemI(ii,ei_link_Gidx_BIPquick)
        !             print *, 'node id ',elemI(ii,ei_node_Gidx_BIPquick)
        !             if (elemI(ii,ei_link_Gidx_BIPquick) .ne. nullvalueI) then 
        !                 print *, trim(link%Names(elemI(ii,ei_link_Gidx_BIPquick))%str)
        !             end if
        !         end if
        !     end if
        ! end do

        !% --- set up the transect arrays
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin init_IC_elem_transect...'
        call init_IC_elem_transect_arrays ()
        call init_IC_elem_transect_geometry ()

            ! ! call util_utest_CLprint ('initial_condition after transect_arrays and _geometry')

        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin init_IC_error_check'
        call init_IC_error_check ()
       
        !% --- identify the small and zero depths (must be done before pack)
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin adjust small/zero depth'
        if (setting%SmallDepth%useMomentumCutoffYN) then
            call adjust_smalldepth_identify_all ()
        else
            !% --- if not using small depth algorithm, then cuttoff is the same as zero depth
            setting%SmallDepth%MomentumDepthCutoff = setting%ZeroValue%Depth
        end if
        call adjust_zerodepth_identify_all ()

            ! ! call util_utest_CLprint ('initial_condition after adjust_smalldepth and _zerodepth')

        !% ---zero out the lateral inflow column
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin init_set_zero_lateral_inflow'
        call init_IC_set_zero_lateral_inflow ()

            ! ! call util_utest_CLprint ('initial_condition after IC_set_zero_lateral_inflow')

        !% --- update time marching type
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin init_IC_solver_select '
        call init_IC_solver_select (whichSolver)

            ! ! call util_utest_CLprint ('initial_condition after IC_solver_select')
       
        !% --- set up all the static packs and masks
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin pack_mask arrays_all'
        call pack_mask_arrays_all ()
  
            ! ! call util_utest_CLprint ('initial_condition after pack_mask_arrays_all')

        !% --- initialize zerovalues for other than depth (must be done after pack)
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin IC_Zerovalues'
        call init_IC_ZeroValues_nondepth ()

            ! ! call util_utest_CLprint ('initial_condition after IC_ZeroValues_nondepth')

        !% --- set all the zero and small volumes
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin adjust small/zero depth 2'
        call adjust_zero_and_small_depth_elem (.true., .true.)

            ! ! call util_utest_CLprint ('initial_condition after adjust_zero_and_small_depth_elem')

        !% --- get the bottom slope
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin IC bottom slope'
        call init_IC_bottom_slope ()

            ! ! call util_utest_CLprint ('initial_condition after IC_bottom_slope')

        !     !% --- get beta (S0/n, used for section factor)
        !    ! if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin IC beta'
        !     call init_IC_beta ()

        !% --- set small volume values in elements
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin init_IC_set_SmallVolumes'
        call init_IC_set_SmallVolumes ()

            ! ! call util_utest_CLprint ('initial_condition after IC_set_SmallVolumes')

        !% --- initialize Preissmann slots
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin init_IC_slot'
        call init_IC_slot ()

            ! ! call util_utest_CLprint ('initial_condition after IC_slot')

        !% --- get the velocity and any other derived data
        !%     These are data needed before bc and aux variables are updated
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin init_IC_derived_data'
        call init_IC_derived_data()

            ! ! call util_utest_CLprint ('initial_condition after IC_derived_data')

        ! !% --- need geometry defined before BC processing  --- 20221024 brh
        ! call geometry_toplevel (whichTM)

        !% --- set the reference head (based on Zbottom values)
        !%     this must be called before bc_update() so that
        !%     the timeseries for head begins correctly
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin init_reference_head'
        call init_reference_head()

        !% --- remove the reference head values from arrays
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin init_subtract_reference_head'
        call init_subtract_reference_head()

            ! ! call util_utest_CLprint ('initial_condition after reference_head')

        !% --- create the packed set of nodes for BC
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin pack_nodes'
        call pack_nodes_BC()
        call util_allocate_bc()

        !% --- initialize Manning's n for forcemain elements
        if (setting%Solver%ForceMain%AllowForceMainTF) call forcemain_ManningsN ()

        !% --- initialize boundary conditions
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin init_bc'
        call init_bc()

            ! ! call util_utest_CLprint ('initial_condition after init_bc')

        !% --- setup the sectionfactor arrays needed for normal depth computation on outfall BC
        if ((setting%Output%Verbose) .and. (this_image() == 1))  print *, "begin init_uniformtable_array"
        call init_uniformtable_array()

        !% --- update the BC so that face interpolation works in update_aux...
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin bc_update'
        call bc_update()
        if (crashI==1) return

            ! ! call util_utest_CLprint ('initial_condition after bc_update')

        !% --- storing dummy values for branches that are invalid
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin branch dummy values'
        call init_IC_branch_dummy_values ()

        !% --- initialize branch values that need to be zero
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin branch zero values'
        call init_IC_branch_zero_values ()

            ! ! call util_utest_CLprint ('initial_condition after IC_branch_dummy_values')

        ! print *, 'TEST20230327 AAA'
        ! print *, elemR(8,er_head),  faceR(9, fr_Head_u), elemR(10,er_Head)
        ! print *, elemR(10,er_Head), faceR(11,fr_Head_u), elemR(11,er_Head)

        !% --- set all the auxiliary (dependent) variables
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin update_aux_variables CC'
        call update_auxiliary_variables_CC (&
            ep_CC_ETM, ep_CC_Open_Elements, ep_CC_Closed_Elements, &
            .true., .false., dummyIdx)

        ! print *, 'TEST20230327 BBBB'
        ! print *, elemR(8,er_head),  faceR(9, fr_Head_u), elemR(10,er_Head)
        ! print *, elemR(10,er_Head), faceR(11,fr_Head_u), elemR(11,er_Head)

        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin update_aux_variables JM'
        call update_auxiliary_variables_JMJB (.false.)

        ! print *, 'TEST20230327 CCCC'
        ! print *, elemR(8,er_head),  faceR(9, fr_Head_u), elemR(10,er_Head)
        ! print *, elemR(10,er_Head), faceR(11,fr_Head_u), elemR(11,er_Head)

            ! ! call util_utest_CLprint ('initial_condition after update_auxiliary_variables')

            !stop 5509873

        !% --- initialize old head 
        !%     HACK - make into a subroutine if more variables need initializing
        !%     after update_aux_var
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin setting old head'
        elemR(:,er_Head_N0) = elemR(:,er_Head)

        !% --- update diagnostic interpolation weights
        !%     (the interpolation weights of diagnostic elements
        !%     stays the same throughout the simulation. Thus, they
        !%     are only needed to be set at the top of the simulation)
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,  'begin init_IC_diagnostic_interpolation_weights'
        call init_IC_diagnostic_interpolation_weights()

        ! print *, 'TEST20230327 DDD'
        ! print *, elemR(8,er_head),  faceR(9, fr_Head_u), elemR(10,er_Head)
        ! print *, elemR(10,er_Head), faceR(11,fr_Head_u), elemR(11,er_Head)

            ! ! call util_utest_CLprint ('initial_condition after IC_diagnostic_interpolation_weights')

        !% --- set small values to diagnostic element interpolation sets
        !%     Needed so that junk values does not mess up the first interpolation
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin  init_IC_small_values_diagnostic_elements'
        call init_IC_small_values_diagnostic_elements

        ! print *, 'TEST20230327 EEE'
        ! print *, elemR(8,er_head),  faceR(9, fr_Head_u), elemR(10,er_Head)
        ! print *, elemR(10,er_Head), faceR(11,fr_Head_u), elemR(11,er_Head)

            ! ! call util_utest_CLprint ('initial_condition after IC_small_values_diagnostic')

        !% --- update faces
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin face_interpolation '
        call face_interpolation (fp_all,.true.,.true.,.true.,.false.,.false.)

        ! print *, 'TEST20230327 FFF'
        ! print *, elemR(8,er_head),  faceR(9, fr_Head_u), elemR(10,er_Head)
        ! print *, elemR(10,er_Head), faceR(11,fr_Head_u), elemR(11,er_Head)

            ! ! call util_utest_CLprint ('initial_condition after face_interpolation')

            !stop 5590873

        !% --- SET THE MONITOR AND ACTION POINTS FROM EPA-SWMM
        if ((setting%Output%Verbose) .and. (this_image() == 1))  print *, "begin controls init monitoring and action from EPSWMM"
        call control_init_monitoring_and_action_from_EPASWMM()

            ! ! call util_utest_CLprint ('initial_condition after control_init_monitoring...')

        !% --- update the initial condition in all diagnostic elements
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin diagnostic_toplevel'
        call diagnostic_toplevel (ep_Diag, fp_Diag, 2)

           ! ! call util_utest_CLprint ('initial_condition after diagnostic_toplevel')

        !% --- ensure that small and zero depth faces are correct
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *,'begin adjust small/zero depth 3'
        call adjust_zero_and_small_depth_face (.false.)

           ! ! call util_utest_CLprint ('initial_condition after adjust_zero_and_small_depth_face')

        !% ---populate er_ones columns with ones
        if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin init_IC_oneVectors'
        call init_IC_oneVectors ()

        ! ! call util_utest_CLprint ('initial_condition at end')

           !stop 666987

        ! print *, 'TEST20230327'
        ! print *, elemR(8,er_head),  faceR(9, fr_Head_u), elemR(10,er_Head)
        ! print *, elemR(10,er_Head), faceR(11,fr_Head_u), elemR(11,er_Head)
        ! stop 23908734


        ! print *, trim(reverseKey(elemI(iet(3),ei_geometryType)))
        ! print *, 'z bottom          ',elemR(iet(3),er_Zbottom)
        ! print *, 'zbreadtmax        ',elemR(iet(3),er_ZbreadthMax)
        ! print *, 'zcrown            ',elemR(iet(3),er_Zcrown)
        ! print *, 'DepthAtBreadthMax ',elemR(iet(3),er_DepthAtBreadthMax)
        ! print *, 'full depth        ',elemR(iet(3),er_FullDepth)
        ! print *, 'full topwidth     ',elemR(iet(3),er_FullTopWidth)
        ! stop 598734

        ! !% --- initialized the max face flowrate
        ! if ((setting%Output%Verbose) .and. (this_image() == 1)) print *, 'begin face_flowrate_max...'
        ! !call util_CLprint()
        ! call sleep(1)
        ! call face_flowrate_max_interior (fp_all)
        ! call face_flowrate_max_shared   (fp_all)

        ! print *, elemI(143,ei_elementType),trim(reverseKey(elemI(143,ei_elementType)))
        ! print *, 'pump ',elemR(143,er_Flowrate), elemR(143,er_Setting)
        ! stop 5098734

        !% Notes on initial conditions brh20211215
        !% dHdA is not initialized in channels except where timemarch is AC


        ! do ii=1,N_elem(this_image())
        !     if (elemI(ii,ei_node_Gidx_BIPquick) .ne. nullvalueI ) then
        !         if (elemI(ii,ei_elementType) .eq. JM) then
        !          print *, ii, elemI(ii,ei_node_Gidx_BIPquick), &
        !             trim(node%Names(elemI(ii,ei_node_Gidx_BIPquick))%str) &
        !             ,' ',trim(reverseKey(elemI(ii,ei_elementType)))
            
        !         end if
        !     end if
        ! end do

        ! do ii=1,N_elem(this_image())
        !     if (elemI(ii,ei_link_Gidx_BIPquick) .ne. nullvalueI ) then
        !         if (elemI(ii,ei_elementType) .eq. pump) then
        !             print *, ii, elemI(ii,ei_link_Gidx_BIPquick), &
        !             trim(link%Names(elemI(ii,ei_link_Gidx_BIPquick))%str) &
        !             ,' ',trim(reverseKey(elemI(ii,ei_elementType)))
        !         end if
        !     end if
        ! end do 
        ! stop 49734
      

   

        !%-------------------------------------------------------------------
        !% Closing
        if (setting%Debug%File%initial_condition) then
            print*, '----------------------------------------------------'
            print*, 'image = ', this_image()
            print*, '.....................elements.......................'
            print*, reversekey(elemI(:,ei_elementType)), 'element type'
            print*, reversekey(elemI(:,ei_geometryType)),'element geometry'
            print*, '-------------------Geometry Data--------------------'
            print*, elemR(:,er_Depth), 'depth'
            print*, elemR(:,er_Area), 'area'
            print*, elemR(:,er_Head), 'head'
            print*, elemR(:,er_Topwidth), 'topwidth'
            print*, elemR(:,er_EllDepth), 'L depth'
            !print*, elemR(:,er_HydDepth), 'hydraulic depth'
            print*, elemR(:,er_HydRadius), 'hydraulic radius'
            print*, elemR(:,er_Perimeter), 'wetted perimeter'
            print*, elemR(:,er_Volume),'volume'
            print*, '-------------------Dynamics Data--------------------'
            print*, elemR(:,er_Flowrate), 'flowrate'
            print*, elemR(:,er_Velocity), 'velocity'
            print*, elemR(:,er_FroudeNumber), 'froude Number'
            print*, elemR(:,er_InterpWeight_uQ), 'timescale Q up'
            print*, elemR(:,er_InterpWeight_dQ), 'timescale Q dn'
            print*, '..................faces..........................'
            print*, faceR(:,fr_Area_u), 'face area up'
            print*, faceR(:,fr_Area_d), 'face area dn'
            print*, faceR(:,fr_Head_u), 'face head up'
            print*, faceR(:,fr_Head_d), 'face head dn'
            print*, faceR(:,fr_Flowrate), 'face flowrate'
            !print*, faceR(:,fr_Topwidth_u), 'face topwidth up'
            !print*, faceR(:,fr_Topwidth_d), 'face topwidth dn'
            ! call execute_command_line('')
        end if

        if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine init_IC_toplevel
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    subroutine init_IC_from_linkdata ()
        !%------------------------------------------------------------------
        !% Description:
        !% get the initial depth, flowrate, and geometry data from links
        !%------------------------------------------------------------------
        !% Declarations:
            integer                                     :: ii, kk, pLink, npack
            integer, pointer                            :: thisLink, eIdx(:)
            integer, dimension(:), allocatable, target  :: packed_link_idx
            integer, dimension(:), allocatable, target  :: ePack
            integer           :: allocation_status, deallocation_status
            character(len=99) :: emsg
            character(64) :: subroutine_name = 'init_IC_from_linkdata'
        !%------------------------------------------------------------------
        !% Preliminaries
            !if (crashYN) return
            if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% pack all the link indexes in an image
        packed_link_idx = pack(link%I(:,li_idx), (link%I(:,li_P_image) == this_image()))

        !% find the number of links in an image
        pLink = size(packed_link_idx)

        ! !% allocate an array for packed comparison
        ! allocate(ePack(N_elem(this_image())), stat=allocation_status, errmsg=emsg)
        ! call util_allocate_check(allocation_status, emsg, 'ePack')

        !% cycle through the links in an image
        do ii = 1,pLink
             print *, ' '
             print *, 'link ',ii
             print *, ' '

            ! % necessary pointers
            thisLink    => packed_link_idx(ii)

                print *, 'calling get barrels'
            call init_IC_get_barrels_from_linkdata(thisLink)

                print *, 'calling get depth'
            !% --- note that this does NOT adjust depth for closed conduit crown height
            call init_IC_get_head_and_depth (thisLink)

                ! print *, 'aaa ',elemR(82,er_Depth), elemR(82,er_Head), elemYN(82,eYN_isPSsurcharged)

                print *, 'calling get flow and roughness'
            call init_IC_get_flow_and_roughness_from_linkdata (thisLink)

                ! print *, 'bbb ',elemR(82,er_Depth), elemR(82,er_Head), elemYN(82,eYN_isPSsurcharged)

               print *, 'calling get elemtype'
            call init_IC_get_elemtype_from_linkdata (thisLink)
            !    print *, thisLink, elemR(54,er_Depth), elemR(54,er_Volume)

            ! print *, 'ccc',elemR(82,er_Depth), elemR(82,er_Head), elemYN(82,eYN_isPSsurcharged)

                print *, 'calling get geometry'
            !% --- note this adjusts depth for closed conduit crown height and sets surcharge
            call init_IC_get_geometry_from_linkdata (thisLink)

            ! print *, 'ddd ',elemR(82,er_Depth), elemR(82,er_Head), elemYN(82,eYN_isPSsurcharged)

               print *, 'calling get flapgate'
            call init_IC_get_flapgate_from_linkdata (thisLink)

            ! print *, 'eee ',elemR(82,er_Depth), elemR(82,er_Head), elemYN(82,eYN_isPSsurcharged)

               print *, 'calling get forcemain'
            call init_IC_get_ForceMain_from_linkdata (thisLink)     

            ! print *, 'fff ',elemR(82,er_Depth), elemR(82,er_Head), elemYN(82,eYN_isPSsurcharged)

               print *, 'calling get culvert'
            call init_IC_get_culvert_from_linkdata(thisLink)
            !    print *, thisLink, elemR(54,er_Depth), elemR(54,er_Volume)

            ! print *, 'ggg ',elemR(82,er_Depth), elemR(82,er_Head), elemYN(82,eYN_isPSsurcharged)

            if ((setting%Output%Verbose) .and. (this_image() == 1)) then
                if (mod(ii,1000) == 0) then
                    print *, '... handling link ',ii
                end if
            end if

        end do

        !%------------------------------------------------------------------
        !% Closing
            !% deallocate the temporary array
            deallocate(packed_link_idx)
            ! deallocate(ePack, stat=deallocation_status, errmsg=emsg)
            ! call util_deallocate_check(deallocation_status, emsg, 'ePack')

            if (setting%Debug%File%initial_condition) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_from_linkdata
!
!==========================================================================
!%==========================================================================
!%
    subroutine init_IC_get_barrels_from_linkdata (thisLink)
        !%-----------------------------------------------------------------
        !% Description:
        !% Sets the number of barrels (default is one)
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in)  :: thisLink
            integer, pointer     :: firstE, lastE, fdn(:), fup(:), eBarrels(:)
            integer, pointer     :: fBarrels(:)
            character(64) :: subroutine_name = 'init_IC_get_barrels_from_linkdata'
        !%-----------------------------------------------------------------
        !% Preliminaries
        !%-----------------------------------------------------------------
        !% Aliases
            firstE      => link%I(thisLink,li_first_elem_idx)
            lastE       => link%I(thisLink,li_last_elem_idx)
            fdn         => elemI(:,ei_Mface_dL)
            fup         => elemI(:,ei_Mface_uL)
            eBarrels    => elemI(:,ei_barrels)
            fBarrels    => faceI(:,fi_barrels)
        !%-----------------------------------------------------------------
        
        !% --- for elements
        eBarrels(firstE:lastE) = link%I(thisLink,li_barrels)

        !print *, 'n barrels ', eBarrels(firstE)

        !% --- for faces
        fBarrels(fup(firstE))       = eBarrels(firstE)
        fBarrels(fdn(firstE:lastE)) = eBarrels(firstE:lastE)

        !% --- note that default for setting%Output%BarrelsExist is false, so
        !%     only need a single multi-barrel to make this true.
        if (any(eBarrels(firstE:lastE) > 1)) setting%Output%BarrelsExist = .true.

    end subroutine init_IC_get_barrels_from_linkdata
!%
!%==========================================================================
!==========================================================================
!
    subroutine init_IC_get_head_and_depth (thisLink)
        !%-----------------------------------------------------------------
        !% Description:
        !% get the initial depth data from links and nodes
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in)  :: thisLink
            integer              :: mm, ei_max, firstidx(1)
            integer, allocatable :: pElem(:)
            integer, pointer     :: LdepthType,  nUp, nDn
            integer, pointer     :: thisP(:)
            logical, pointer     :: hasFlapGate
            real(8), pointer     :: DepthUp, DepthDn
            real(8), pointer     :: zLinkUp, zLinkDn, Slope
            real(8), pointer     :: eDepth(:), eHead(:), eLength(:), eZbottom(:)
            real(8)              :: kappa,  headUp, headDn, linkLength, length2Here
            real(8)              :: dDelta, hDelta
            
            character(64) :: subroutine_name = 'init_IC_get_depth_from_linkdata'
        !%-----------------------------------------------------------------
        !% Preliminaries
            !if (crashYN) return
            if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------
        !% Aliases
            !% --- type of initial depth type
            LdepthType  => link%I(thisLink,li_InitialDepthType)
            !% --- upstream and downstream nodes
            nUp      => link%I(thisLink,li_Mnode_u)
            nDn      => link%I(thisLink,li_Mnode_d)
            !% --- flapgate on downstream node (e.g. outfall)
            hasFlapGate => node%YN(nDn,nYN_hasFlapGate)
            !% --- link upstream and downstream bottom elevation
            zLinkUp  => link%R(thisLink,lr_ZbottomUp)
            zLinkDn  => link%R(thisLink,lr_ZbottomDn)
            Slope    => link%R(thisLink,lr_Slope)
            !% --- depths upstream and downstream on link (not yet initialized)
            DepthUp     => link%R(thisLink,lr_InitialUpstreamDepth)
            DepthDn     => link%R(thisLink,lr_InitialDnstreamDepth)
            !% 
            eLength   => elemR(:,er_Length)
            eDepth    => elemR(:,er_Depth)
            eHead     => elemR(:,er_Head)
            eZbottom  => elemR(:,er_Zbottom)
        !%-----------------------------------------------------------------

            ! print *, ' '
            ! print *, 'in init_IC_get_depth: link ',trim(link%Names(thisLink)%str)

        !% --- Head upstream
        headUp = node%R(nUp,nr_Zbottom) + node%R(nUp,nr_InitialDepth)
        !% --- provisional head downstream
        headDn = node%R(nDn,nr_Zbottom) + node%R(nDn,nr_InitialDepth)

        ! print *, 'headUp ', headUp
        ! print *, 'headDn ', headDn


        !% --- set upstream link depths including effects of offsets
        !%     where head upstream is less than zbottom, depth is zero
        DepthUp = max(headUp - zLinkUp, zeroR)

        !print *, 'DepthUp ',DepthUp

            ! print *, ' =================================== '
            ! print *, 'thisLink ',thisLink, ' ',trim(link%Names(thisLink)%str)
            ! print *, 'node up  ',nUp, ' ',trim(node%Names(nUp)%str)
            ! print *, 'node dn  ',nDn, ' ',trim(node%Names(nDn)%str)
            ! print *, 'ZlinkUp  ',zLinkUp
            ! print *, 'zlinkDn  ',zLinkDn
            ! print *, 'DepthUp  ',node%R(nUp,nr_InitialDepth)
            ! print *, 'DepthDn  ',node%R(nDn,nr_InitialDepth)
            ! print *, 'HeadUp   ',headUp
            ! print *, 'HeadDn   ',headDn
            ! print *, 'Slope    ',Slope
            ! print *, ' node up z bottom  ', node%R(nUp,nr_Zbottom)
            ! print *, ' node dn z bototm  ', node%R(nDn,nr_Zbottom)
            ! print *, ' node up init depth',node%R(nUp,nr_InitialDepth)
            ! print *, ' node dn init depth',node%R(nDn,nr_InitialDepth)
            ! print *, ' '

    
        !% --- set downstream link depths including effects of offsets
        !%     where downstream head is less than zbottom, depth is zero
        DepthDn = max(headDn - zLinkDn, zeroR)

        !print *, 'Depth Dn ',DepthDn

        !% --- check for a downstream gate on the node
        !%     adjust depths and head as needed
        ! if (node%YN(nDn,nYN_hasFlapGate)) then
        !     !print *, 'has Flap Gate'
        !     if (DepthUp == zeroR) then
        !         !% --- if zero depth upstream, then downstream is also zero
        !         !%     and we switch to a uniform depth interpolation scheme
        !         !%     so that the entire link is dry
        !         DepthDn = zeroR
        !         headDn  = zLinkDn
        !         LdepthType = UniformDepth
        !     else
        !         !% --- for positive upstream depth
        !         !%     if upstream head is lower than downstream head
        !         !%     then flap gate is closed
        !         if (headUp < headDn) then
        !             !% --- closed flap gate
        !             !%     set downstream at the upstream head (ponding at gate)
        !             headDn  = headUp
        !             DepthDn = headDn - zLinkDn
        !             !% --- ensure the elements are handled by fixed head
        !             LdepthType = FixedHead
        !         else
        !             !% --- for upstream head > downstream head
        !             !%     ensure interpolation over link
        !             select case (LdepthType)
        !             case (UniformDepth, FixedHead)
        !                 LdepthType = LinearlyVaryingDepth
        !             case default 
        !                 !% --- continue with selected interpolation type
        !             end select
        !         end if
        !     end if
        ! end if

       ! print *, 'Depth Dn BBB ',DepthDn

        !% --- pack the elements for this link
        pElem = pack(elemI(:,ei_Lidx), (elemI(:,ei_link_Gidx_BIPquick) == thisLink))

        !print *, 'pElem ', pElem

        !% --- error checking, the upstream should be the first element in the pack
        firstidx = findloc(elemI(pElem,ei_link_Pos),1)
        if (firstidx(1) .ne. 1) then
            print *, 'CODE ERROR'
            print *, 'Possible problem in element ordering in a link'
            print *, 'error with link ',trim(link%Names(thisLink)%str)
            print *, elemI(pElem,ei_link_Pos)
            call util_crashpoint(55872)
        end if

        !% --- total depth delta
        ! dDelta = DepthUp - DepthDn
        hDelta = headUp - headDn

        ! print *, 'Depth Delta ',dDelta

        !% --- total length of all elements in link
        linkLength = sum(eLength(pElem))

        !% -- initialize length measure from the upper end of the link
        !%    to an interative element center
        length2Here = zeroR

        !% ---- set the initial depths and heads in each element
        if ((DepthUp > setting%ZeroValue%Depth) .and. &
            (DepthDn > setting%ZeroValue%Depth) ) then
            !% --- distribute head linearly along the link
            !% saz 01/11/2023
            !%     in this new method head values in elements are linearly interpolated
            !%     then the depths are recovered from those heads.
            do mm=1,size(pElem)
                !% --- use the length from upstream face to center of this element
                length2Here       = length2Here + onehalfR * eLength(pElem(mm))
                !% --- head by linear interpolation
                ! eHead(pElem(mm)) = headUp - Slope * length2Here
                eHead(pElem(mm)) = headUp - hDelta * length2Here / linkLength
                !% --- depth from head
                eDepth(pElem(mm)) = max(eHead(pElem(mm)) - eZbottom(pElem(mm)), zeroR) 
                !% --- add the remainder of this element to the length
                length2Here       = length2Here + onehalfR * eLength(pElem(mm))
            end do


        elseif ((DepthUp > setting%ZeroValue%Depth) .and. &
                (DepthDn .le. setting%ZeroValue%Depth)) then 
                  
            if (headUp .le. zLinkDn) then        
                !% --- conduit sloping upwards
                !%     upstream depth provides uniform head over the entire link  
                eHead(pElem) = headUp 
                eDepth(pElem) = eHead(pElem) - eZbottom(pElem)
            else 
                do mm=1,size(pElem)
                    !% --- use the length from upstream face to center of this element
                    length2Here       = length2Here + onehalfR * eLength(pElem(mm))
                    !% --- head by linear interpolation
                    ! eHead(pElem(mm)) = headUp - Slope * length2Here
                    eHead(pElem(mm)) = headUp - hDelta * length2Here / linkLength
                    !% --- depth from head
                    eDepth(pElem(mm)) = max(eHead(pElem(mm)) - eZbottom(pElem(mm)), zeroR) 
                    !% --- add the remainder of this element to the length
                    length2Here       = length2Here + onehalfR * eLength(pElem(mm))
                end do
                !% --- implies head upstream goes to zero depth downstream.
                !%     which is an inconsistent initial condition
                !%     This condition is allowed (needed for free overflow)
                !%     but could cause issues!

                ! print *, 'Inconsistent intial conditions at link '
                ! print *, trim(link%Names(thisLink)%str)
                ! print *, 'Depth at downstream node is zero and head'
                ! print *, 'at upstream node would imply flow from a '
                ! print *, 'non-zero depth to a zero depth, which is '
                ! print *, 'not a valid starting conditions. Please check'
                ! print *, 'the initial depths at the nodes connecting'
                ! print *, 'to this link'
                ! call util_crashpoint(7098234)
            end if



        elseif ((DepthDn > setting%ZeroValue%Depth) .and. &
                (DepthUp .le. setting%ZeroValue%Depth)) then
            !% --- downstream depth provides uniform head over the entire link
            if (HeadDn .le. zLinkUp) then 
                eHead(pElem) = HeadDn
                eDepth(pElem) = eHead(pElem) - eZbottom(pElem)
            else
                !% --- implied reverse gradient is not allowed
                print *, 'Inconsistent initial conditions at link '
                print *, trim(link%Names(thisLink)%str)
                print *, 'Depth at upstream node is zero and non-zero depth at'
                print *, 'downstream noded implies upstream node depth'
                print *, 'should be at least ',HeadDn - zLinkUp, ' meters'
                print *, 'or ',(HeadDn - zLinkUp)*3.28084d0, 'feet'
                print *, 'Please check the initial depths at the nodes connecting'
                print *, 'to this link.'
                call util_crashpoint(40187339)
            end if

        elseif ((DepthDn .le. setting%ZeroValue%Depth) .and. &
                (DepthUp .le. setting%ZeroValue%Depth)) then
            !% --- zero depths everywhere along the element.        
            eDepth(pElem) = setting%ZeroValue%Depth  * 0.99d0 
            eHead (pElem) = eZBottom(pElem) + setting%ZeroValue%Depth  * 0.99d0        
        else 
            print *, 'CODE ERROR: unexpected else.'
            print *, 'code should not have reached this point'
            call util_crashpoint(8898723)
        end if

        where(eDepth(pElem) .le. setting%ZeroValue%Depth)
            eDepth(pElem) = setting%ZeroValue%Depth * 0.99d0
            eHead(pElem)  = eZbottom(pElem) + setting%ZeroValue%Depth * 0.99d0
        endwhere

        !print *, 'in here in link ', elemR(15,er_Depth)
        !%------------------------------------------------------------------------------------
        !% saz 01/11/2023
        !% commented out to set up link depth based on node heads only 

        !% --- Single-length elements can only be UniformDepth or
        !%     FixedHead. Set interpolated schemes to UniformDepth
        ! if (size(pElem) == 1) then
        !     select case (LdepthType)
        !     case (LinearlyVaryingDepth, ExponentialDepth)
        !         LdepthType = UniformDepth
        !     case (default)
        !         !% no change
        !     end select
        ! end if

        ! print *, 'Depth Type = ',reverseKey(LdepthType)

        !% ---set the depths in link elements from links
        !%    Note these depths are the combination of water and sediment
        ! select case (LdepthType)

        !     case (UniformDepth)
        !         !% --- uniform depth uses the average of upstream and downstream depths
        !         eDepth(pElem) = onehalfR * (DepthUp + DepthDn)
        

        !     case (LinearlyVaryingDepth)
        !         !% --- linearly-varying depth distribution
        !         do mm=1,size(pElem)
        !             !% --- use the length from upstream face to center of this element
        !             length2Here       = length2Here + onehalfR * eLength(pElem(mm))
        !             !% --- depth by linear interpolation
        !             eDepth(pElem(mm)) = DepthUp - dDelta * length2Here /linkLength
        !             !% --- add the remainder of this element to the length
        !             length2Here       = length2Here + onehalfR * eLength(pElem(mm))
        !         end do

        !     case (ExponentialDepth)
        !         !% --- if the link has exponentially increasing or decreasing depth

        !         do mm=1,size(pElem)
        !             !% --- use the length from upstream face to center of this element
        !             length2Here       = length2Here + onehalfR * eLength(pElem(mm))
        !             !% --- normalized exponential decay
        !             kappa = - exp(oneR) * length2Here / linkLength
        !             !% --- depth by linear interpolation
        !             eDepth(pElem(mm)) = DepthDn + dDelta * exp(-kappa)
        !             !% --- add the remainder of this element to the length
        !             length2Here       = length2Here + onehalfR * eLength(pElem(mm))
        !         end do

        !     case (FixedHead)    
        !         !% --- set the downstream depth as a fixed head (ponding)
        !         !%     over all the elements in the link.
        !         eDepth(pElem) = max(headDn - eZbottom(pElem), zeroR)
            
        !     case default
        !         print *, 'In ', subroutine_name
        !         print *, 'CODE ERROR: unexpected initial depth type #', LdepthType,'  in link, ', thisLink
        !         print *, 'which has key ',trim(reverseKey(LdepthType)) 
        !         !stop 
        !         call util_crashpoint(83753)
        !         !return
        ! end select

        !% --- set zero values to zerodepth
        where (eDepth(pElem) < setting%ZeroValue%Depth)
            eDepth(pElem) = setting%ZeroValue%Depth * 0.99d0 
        endwhere

        ! print *, ' '
        ! print *, 'in init_IC_get_depth'
        ! print *, 'pElem ',pElem
        ! print *, 'head  ',eHead(pElem)
        ! print *, 'zbottom ', eZbottom(pElem)
        ! print *, 'depth ',eDepth(pElem)
        ! print *, ' '
    
        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_get_head_and_depth
!
!==========================================================================
!==========================================================================
!
    subroutine init_IC_get_flow_and_roughness_from_linkdata (thisLink)
        !%-----------------------------------------------------------------
        !% Description:
        !% get the initial flowrate and roughness data from links
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisLink
            integer :: ii
            integer, pointer :: firstelem, lastelem, tNode
            character(64) :: subroutine_name = 'init_IC_get_flow_and_roughness_from_linkdata'
        !--------------------------------------------------------------------------
            if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

            firstelem => link%I(thisLink,li_first_elem_idx)
            lastelem  => link%I(thisLink,li_last_elem_idx)
        !--------------------------------------------------------------------------   
        !% --- handle all the initial conditions that don't depend on geometry type
        where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
            elemR(:,er_Flowrate)           = link%R(thisLink,lr_FlowrateInitial) / link%I(thisLink,li_barrels)
            elemR(:,er_Flowrate_N0)        = link%R(thisLink,lr_FlowrateInitial) / link%I(thisLink,li_barrels)
            elemR(:,er_Flowrate_N1)        = link%R(thisLink,lr_FlowrateInitial) / link%I(thisLink,li_barrels)
            elemR(:,er_ManningsN)          = link%R(thisLink,lr_Roughness)
            !% --- distribute minor losses uniformly over all the elements in thi link
            elemR(:,er_Kconduit_MinorLoss) = link%R(thisLink,lr_Kconduit_MinorLoss) / (real(lastelem - firstelem + oneI,8))
            elemR(:,er_FlowrateLimit)      = link%R(thisLink,lr_FlowrateLimit)
            elemR(:,er_SeepRate)           = link%R(thisLink,lr_SeepRate)
            elemR(:,er_ManningsN_Dynamic)  = elemR(:,er_ManningsN)   
        endwhere

        !% --- assign minor losses
        !%     If connection not to an nJM junction, then add the entry/exit losses
        !%     to the Kconduit_MinorLoss
        !elemR(firstelem,er_Kentry_MinorLoss) = link%R(thisLink,lr_Kentry_MinorLoss)
        !elemR(lastelem ,er_Kexit_MinorLoss)  = link%R(thisLink,lr_Kexit_MinorLoss)
        tNode => link%I(thisLink,li_Mnode_u)
        if (node%I(tNode,ni_node_type) == nJM) then 
            !% --- this is an outlet from an nJM junction, which uses the entry minor
            !%     loss from the downstream link
            !%     HACK -- need a way to have different entry/exit for flow reversal
            elemR(firstelem,er_KJunction_MinorLoss) = link%R(thisLink,lr_Kentry_MinorLoss)
        else
            !% --- this is an nJ2 or a nBC junction, so the entry loss is added to the conduit loss
            elemR(firstelem,er_Kconduit_MinorLoss) = elemR(firstelem,er_Kconduit_MinorLoss) &
                                                   + link%R(thisLink,lr_Kentry_MinorLoss)
        end if

        tNode => link%I(thisLink,li_Mnode_d)
        if (node%I(tNode,ni_node_type) == nJM) then 
            !% --- this is an inlet to an nJM junction, which uses the exit minor
            !%     loss from the upstream link as the entrance loss to the junction
            !%     HACK -- need a way to have different entry/exit for flow reversal
            elemR(lastelem,er_KJunction_MinorLoss) = link%R(thisLink,lr_Kexit_MinorLoss)
        else
            !% --- this is an nJ2 or a nBC junction, so the exit loss is added to the conduit loss
            elemR(lastelem,er_Kconduit_MinorLoss) = elemR(lastelem,er_Kconduit_MinorLoss) &
                                                   + link%R(thisLink,lr_Kexit_MinorLoss)
        end if

        ! print *, ' '
        ! print *, 'in ',trim(subroutine_name)
        ! do ii=1,size(elemI,1)
        !     print *, ii, elemI(ii,ei_link_Gidx_BIPquick), thisLink
        ! end do

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_get_flow_and_roughness_from_linkdata
!
!==========================================================================
!==========================================================================
!
    subroutine init_IC_get_elemtype_from_linkdata (thisLink)
        !--------------------------------------------------------------------------
        !% get the geometry data from links
        !--------------------------------------------------------------------------

            integer, intent(in) :: thisLink
            integer, pointer    :: linkType
            integer             :: ii

            character(64) :: subroutine_name = 'init_IC_get_elemtype_from_linkdata'
        !--------------------------------------------------------------------------
            !if (crashYN) return
            if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% necessary pointers
        linkType      => link%I(thisLink,li_link_type)

        select case (linkType)

            case (lChannel)
                where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                    elemI(:,ei_elementType)     = CC
                    elemI(:,ei_HeqType)         = time_march
                    elemI(:,ei_QeqType)         = time_march
                    elemYN(:,eYN_canSurcharge)  = .false.
                endwhere


            case (lpipe)
                where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                    elemI(:,ei_elementType)     = CC
                    elemI(:,ei_HeqType)         = time_march
                    elemI(:,ei_QeqType)         = time_march
                    elemYN(:,eYN_canSurcharge)  =  .true.
                endwhere

            case (lweir)
                where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                    elemI(:,ei_elementType)         = weir
                    elemI(:,ei_QeqType)             = diagnostic
                    elemI(:,ei_HeqType)             = notused
                    elemYN(:,eYN_canSurcharge)      = link%YN(thisLink,lYN_weir_CanSurcharge)
                endwhere

            case (lOrifice)
                where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                    elemI(:,ei_elementType)            = orifice
                    elemI(:,ei_QeqType)                = diagnostic
                    elemI(:,ei_HeqType)                = notused
                    elemYN(:,eYN_canSurcharge)         = .true.
                endwhere

            case (lPump)
                where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                    elemI(:,ei_elementType)            = pump
                    elemI(:,ei_QeqType)                = diagnostic
                    elemI(:,ei_HeqType)                = notused
                    elemYN(:,eYN_canSurcharge)         = .true.
                endwhere

            case (lOutlet)
                where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                    elemI(:,ei_elementType)            = outlet
                    elemI(:,ei_QeqType)                = diagnostic
                    elemI(:,ei_HeqType)                = notused
                    elemYN(:,eYN_canSurcharge)         = .true.
                endwhere

            case default

                print *, 'in ', trim(subroutine_name)
                print *, 'CODE ERROR: unexpected link type, ', linkType,'  in the network'
                if ((linkType > 0) .and. (linkType < size(reverseKey))) then
                    print *, 'which has key number ',trim(reverseKey(linkType))
                else 
                    print *, 'key number is outside of allowed bounds.'
                end if 
                call util_crashpoint(65343)
        end select

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_get_elemtype_from_linkdata
!
!==========================================================================
!==========================================================================
!
    subroutine init_IC_get_flapgate_from_linkdata (thisLink)
        !%-----------------------------------------------------------------
        !% Description:
        !% Sets a flap gate (if it exists) to the last element in a link
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in)  :: thisLink
            logical, pointer     :: hasFlapGate
            integer, pointer     :: firstE, lastE
            
            character(64) :: subroutine_name = 'init_IC_get_flapgate_linkdata'
        !%-----------------------------------------------------------------
        !% Preliminaries
        !%-----------------------------------------------------------------
        !% Aliases
            hasFlapGate => link%YN(thisLink,lYN_hasFlapGate)
            firstE      => link%I(thisLink,li_first_elem_idx)
            lastE       => link%I(thisLink,li_last_elem_idx)
        !%-----------------------------------------------------------------
        !% --- initialize all conduit link flap gates to false
        elemYN(firstE:lastE,eYN_hasFlapGate) = .false.

        !% --- set any flap gate to the last element in the conduit
        if (hasFlapGate) then
            elemYN(lastE,eYN_hasFlapGate) = .true.
        end if
        
    end subroutine init_IC_get_flapgate_from_linkdata
!
!==========================================================================
!==========================================================================
!
    subroutine init_IC_get_ForceMain_from_linkdata (thisLink)
        !%-----------------------------------------------------------------
        !% Description:
        !% Sets the Force main coefficients
        !%-----------------------------------------------------------------
        !% Declarations:
        integer, intent(in)  :: thisLink
        integer, pointer     :: firstE, lastE, linkType, linkGeo
        
        character(64) :: subroutine_name = 'init_IC_get_ForceMain_from_linkdata'
        !%-----------------------------------------------------------------
        !% Preliminaries

        !%-----------------------------------------------------------------
        !% Aliases
        firstE      => link%I(thisLink,li_first_elem_idx)
        lastE       => link%I(thisLink,li_last_elem_idx)
        linkType    => link%I(thisLink,li_link_type)
        linkGeo     => link%I(thisLink,li_geometry)
        !%-----------------------------------------------------------------
        !%
        !% only call for pipes
        if (linkType .eq. lpipe) then
            !% --- if UseForceMain
            if (setting%Solver%ForceMain%AllowForceMainTF) then
                !% --- if FMallClosedConduits
                if (setting%Solver%ForceMain%FMallClosedConduitsTF) then
                    !% --- forcing all closed conduits to be Force Main
                    call init_IC_set_forcemain_elements (firstE, lastE, thisLink)
                else 
                    !% --- handle links designated as force main in SWMM input file
                    if (linkGeo .eq. lForce_main) then
                        call init_IC_set_forcemain_elements (firstE, lastE, thisLink)
                    else
                        !% --- not a force main
                        elemYN(firstE:lastE,eYN_isForceMain)      = .false.
                        elemSR(firstE:lastE,esr_Conduit_ForceMain_Coef)   = nullvalueR
                        elemSI(firstE:lastE,esi_Conduit_Forcemain_Method) = NotForceMain
                    end if
                end if
            else    
                !% --- if force mains are turned off
                elemYN(firstE:lastE,eYN_isForceMain)      = .false.
                elemSR(firstE:lastE,esr_Conduit_ForceMain_Coef)   = nullvalueR
                elemSI(firstE:lastE,esi_Conduit_Forcemain_Method) = NotForceMain
            end if
        else
            !% --- do not initialize FM storage for weir etc.
        end if

    end subroutine init_IC_get_ForceMain_from_linkdata
!    
!==========================================================================
!==========================================================================
!    
    subroutine init_IC_set_forcemain_elements (firstE, lastE, thisLink)
        !%-----------------------------------------------------------------
        !% Description:
        !% sets the force main conditions between the first element (firstE)
        !% and last element (lastE) of a link (thisLink)
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: firstE, lastE, thisLink
        !%-----------------------------------------------------------------
        !%-----------------------------------------------------------------
        !%
        !% --- set these elements to a force main    
        elemYN(firstE:lastE,eYN_isForceMain)      = .true.

        !% --- check as to whether SWMMinput methods or an overwrite of
        !%     the JSON file is used
        if (setting%Solver%ForceMain%UseSWMMinputMethodTF) then 
            !% --- using SWMM input method
            elemSR(firstE:lastE,esr_Conduit_ForceMain_Coef)   = link%R(thisLink,lr_ForceMain_Coef)
            elemSI(firstE:lastE,esi_Conduit_Forcemain_Method) = setting%SWMMinput%ForceMainEquationType
        else
            !% --- overwriting with default method from JSON file
            select case (setting%Solver%ForceMain%Default_method)
            case (HazenWilliams)
                elemSI(firstE:lastE,esi_Conduit_Forcemain_Method) = HazenWilliams
                elemSR(firstE:lastE,esr_Conduit_ForceMain_Coef)   = setting%Solver%ForceMain%Default_HazenWilliams_coef
            case (DarcyWeisbach)
                elemSI(firstE:lastE,esi_Conduit_Forcemain_Method) = DarcyWeisbach
                elemSR(firstE:lastE,esr_Conduit_ForceMain_Coef)   = setting%Solver%ForceMain%Default_DarcyWeisbach_roughness_mm
            case default 
                print *, 'CODE ERROR: unexpected case default'
                call util_crashpoint(7729873)
            end select
        end if

        !% --- error checking
        !%     Examines if roughness values for FM are consistent with what's expected for the
        !%     Hazen-Williams or Darcy-Weisbach approaches.
        if (setting%Solver%ForceMain%errorCheck_RoughnessTF) then 
            if (elemSI(firstE,esi_Conduit_Forcemain_Method) .eq. HazenWilliams) then 
                !% --- for Hazen Williams Force main
                if (elemSR(firstE,esr_Conduit_ForceMain_Coef) < 90.0) then
                    print *, 'POSSIBLE CONFIGURATION ERROR: Force Main Coefficients'
                    print *, 'The Hazen-Williams equation for Force Mains is invoked '
                    print *, '   however the HW roughness coefficient seems small for'
                    print *, '   an HW solution.' 
                    print *, 'At link name ',trim(link%Names(thisLink)%str)
                    print *, '  the HW roughness was ', elemSR(firstE,esr_Conduit_ForceMain_Coef)
                    print *, 'This might be because the roughness is for a Darcy-Weisbach'
                    print *, '  force main, in which case you need to change the FORCE_MAIN_EQUATION'
                    print *, '  in the SWMM input file.'
                    print *, 'If this coefficient (and all other small coefficients) are OK'
                    print *, '  then use setting.Solver.ForceMain.errorCheck_RoughnessTF = false'
                    print *, '  and re-run to pass this error check point.'
                    call util_crashpoint(509874)
                end if
            else
                !% --- for Darcy-Weisbach Force Main
                if (elemSR(firstE,esr_Conduit_ForceMain_Coef)*1000.d0 > 60) then 
                    print *, 'POSSIBLE CONFIGURATION ERROR: Force Main Coefficients'
                    print *, 'The Darcy-Weisbach equation for Force Mains is invoked '
                    print *, '   however the DW roughness coefficient seems large for'
                    print *, '   a DW solution.' 
                    print *, 'At link name ',trim(link%Names(thisLink)%str)
                    print *, '  the input DW roughness (in SI) was ', elemSR(firstE,esr_Conduit_ForceMain_Coef)*1000.d0, ' mm'
                    print *, 'This might be because the roughness is for a Hazen-Williams'
                    print *, '  force main, in which case you need to change the FORCE_MAIN_EQUATION'
                    print *, '  in the SWMM input file.'
                    print *, 'If this coefficient (and all other large coefficients) are OK'
                    print *, '  then use setting.Solver.ForceMain.errorCheck_RoughnessTF = false'
                    print *, '  and re-run to pass this error check point.'
                    call util_crashpoint(5098742)
                end if
            end if
        else
            !% -- no error checking
        end if

    end subroutine init_IC_set_forcemain_elements
!%    
!%==========================================================================  
!%==========================================================================
!%
    subroutine init_IC_get_culvert_from_linkdata (thisLink)
        !%-----------------------------------------------------------------
        !% Description:
        !% Sets up the culvert
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in)  :: thisLink
            integer, pointer     :: firstE, lastE, thisC
            integer              :: ii
            character(64) :: subroutine_name = 'init_IC_get_culvert_from_linkdata'
        !%-----------------------------------------------------------------
        !% Preliminaries
            !% --- if no culvert, return
            if (link%I(thisLink,li_culvertCode) == zeroI) return
        !%-----------------------------------------------------------------
        !% Aliases
            firstE      => link%I(thisLink,li_first_elem_idx)
            lastE       => link%I(thisLink,li_last_elem_idx)
        !%-----------------------------------------------------------------

        !% --- look for culvert
        if (link%I(thislink,li_culvertCode) == 0) then
            elemYN(firstE:lastE,eYN_isCulvert) = .false.
            !% --- elemSI is not initialized for non-culverts
            return
        elseif ((link%I(thislink,li_culvertCode) > 0 ) .and. &
                (link%I(thislink,li_culvertCode) <= NculvertTypes)) then
            elemYN(firstE:lastE,eYN_isCulvert) = .true.
        else
            print *, 'USER CONFIGURATION ERROR'
            print *, 'Culvert Code found with value of ',link%I(thisLink,li_culvertCode)
            print *, 'for link # ',thisLink
            print *, 'which is the link named ',trim(link%Names(thisLink)%str)
            print *, 'Allowable culvert codes are zero or greater and'
            print *, 'less than or equal to ',NculvertTypes
            call util_crashpoint(6628732)
        end if
   
        !% --- culverts are only defined on closed-conduit links
        if (link%I(thislink,li_link_type) == lpipe) then
            !% --- store the culvert code for all culvert elements
            elemSI(firstE:lastE,esi_Conduit_Culvert_Code) = link%I(thisLink,li_culvertCode)

            !% --- identify parts of the culvert
            if (firstE == lastE) then
                !% if only 1 element in link
                elemSI(firstE,esi_Conduit_Culvert_Part) = Culvert_InOut 
                elemSI(firstE,esi_Conduit_Culvert_OutletID) = firstE
            else 
                !% --- designate inlet
                elemSI(firstE,esi_Conduit_Culvert_Part) = Culvert_Inlet
                !% --- designate outlet
                elemSI(lastE ,esi_Conduit_Culvert_Part) = Culvert_Outlet
                !% --- designate interior barrel elements
                elemSI(firstE+1:lastE-1,esi_Conduit_Culvert_Part) = Culvert_Barrel
                !% --- store outlet location 
                elemSI(firstE,esi_Conduit_Culvert_OutletID) = lastE
            end if

            !% --- LOCAL STORE OF CULVERT VALUES:

            !% --- pointer for covenience
            thisC  => elemSI(firstE,esi_Conduit_Culvert_Code)

            !% --- convert the EquationForm real in the culvertValue to an integer
            if (culvertValue(thisC,1) == 1.d0) then 
                elemSI(firstE:lastE,esi_Conduit_Culvert_EquationForm) = oneI
            elseif (culvertValue(thisC,1) == 2.d0) then   
                elemSI(firstE:lastE,esi_Conduit_Culvert_EquationForm) = twoI
            else 
                print *, 'CODE ERROR: unexpected else'
                call util_crashpoint(739874)
            end if

            !% -- real data from culvertValue
            elemSR(firstE:lastE,esr_Conduit_Culvert_K)   = culvertValue(thisC,2)
            elemSR(firstE:lastE,esr_Conduit_Culvert_M)   = culvertValue(thisC,3)
            elemSR(firstE:lastE,esr_Conduit_Culvert_C)   = culvertValue(thisC,4)
            elemSR(firstE:lastE,esr_Conduit_Culvert_Y)   = culvertValue(thisC,5)
            elemSR(firstE:lastE,esr_Conduit_Culvert_SCF) = culvertValue(thisC,6)

        else 
            !% --- error: culvert not allowed for non-conduit elements
            print *, 'CONFIGURATION ERROR'
            print *, 'A culvert code has been found for a link that'
            print *, 'is not a closed conduit.  Only closed conduits'
            print *, 'can be culverts'
            print *, 'Problem for link # ',thisLink
            print *, 'Link Name ',trim(link%Names(thisLink)%str)
            call util_crashpoint(833287)
        end if
        
    end subroutine init_IC_get_culvert_from_linkdata
!%
!%==========================================================================    
!%==========================================================================
!%
    subroutine init_IC_get_geometry_from_linkdata (thisLink)
        !--------------------------------------------------------------------------
        !% get the geometry data from links
        !--------------------------------------------------------------------------

            integer, intent(in) :: thisLink
            integer, pointer    :: linkType
            character(64) :: subroutine_name = 'init_IC_get_geometry_from_linkdata'
        !--------------------------------------------------------------------------
            !if (crashYN) return
            if (setting%Debug%File%initial_condition) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% necessary pointers
        linkType      => link%I(thisLink,li_link_type)

        print *, 'in ',trim(subroutine_name)
        print *, thisLink, linkType
        print *, trim(reverseKey(linkType))

        select case (linkType)

            case (lChannel)
                
                !% get geometry data for channels
                call init_IC_get_channel_geometry (thisLink)

            case (lpipe)
                !% get geometry data for conduits
                call init_IC_get_conduit_geometry (thisLink)

            case (lweir)
                !% get geometry data for weirs
                call init_IC_get_weir_geometry (thisLink)

            case (lOrifice)
                !% get geometry data for orifices
                call init_IC_get_orifice_geometry (thisLink)

            case (lPump)
                !% get geometry data for pump
                call init_IC_get_pump_geometry (thisLink)

            case (lOutlet)
                !% get geometry data for link outlets
                ! print *, 'CODE ERROR:  an outlet link in the SWMM input file was found.'
                ! print *, 'This feature is not yet available in SWMM5+'
                ! call util_crashpoint(4409872)
                call init_IC_get_outlet_geometry (thisLink)

            case default

                print *, 'In ', subroutine_name
                print *, 'CODE ERROR: unexpected link type, ', linkType,'  in the network'
                print *, 'which has key ',trim(reverseKey(linkType))
                !stop 
                call util_crashpoint(99834)
                !return

        end select

        print *, 'at exit'         

        !call geometry_table_initialize ()

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_get_geometry_from_linkdata
!
!==========================================================================
!==========================================================================
!
    subroutine init_IC_get_channel_geometry (thisLink)
        !%-----------------------------------------------------------------
        !% Description:
        !% get the geometry data for open channel links
        !% and calculate element volumes
        !% Note that the "FullDepth" must be defined for open channels.    
        !%-------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisLink
            integer, pointer    :: geometryType, link_tidx, eIdx(:), thisP(:)
            integer, pointer    :: fUp(:), fDn(:)

            integer :: Npack, ii, mm

            !% for pure calls to scalars
            integer, dimension(1) :: Iarg
            real(8), dimension(1) :: Rarg
        

            real(8), pointer    :: depth(:)
            real(8), pointer    :: fullarea(:), fullperimeter(:)
            real(8), pointer    :: fulltopwidth(:), initialDepth(:)
            real(8), pointer    :: fullhydradius(:), fulldepth(:)

            character(64) :: subroutine_name = 'init_IC_get_channel_geometry'
        !%--------------------------------------------------------------------
        !% Preliminaries:
            if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%--------------------------------------------------------------------
        !% Aliases:

            !% --- pointer to geometry type
            geometryType => link%I(thisLink,li_geometry)

            !% --- pointer to element indexes
            eIdx         => elemI(:,ei_Lidx)

            !% --- pointer to face indexes
            fUp          => elemI(:,ei_Mface_uL)
            fDn          => elemI(:,ei_Mface_dL)

            !% --- pack the elements for this link in temporary array
            Npack = count(elemI(:,ei_link_Gidx_BIPquick) == thisLink)
            if (Npack < 1) return 
            elemI(1:Npack,ei_Temp01) = pack(eIdx,elemI(:,ei_link_Gidx_BIPquick) == thisLink)
            thisP => elemI(1:Npack,ei_Temp01)

            initialDepth => elemR(:,er_Temp01)
            depth        => elemR(:,er_Depth)
            fullarea     => elemR(:,er_FullArea)
            fulldepth    => elemR(:,er_FullDepth)
            fullhydradius=> elemR(:,er_FullHydRadius)
            fullperimeter=> elemR(:,er_FullPerimeter)
            fulltopwidth => elemR(:,er_FullTopwidth)
        !%--------------------------------------------------------------------

        !% --- temporarily store initial depth in temp array so that full depth
        !%     can replace it for computing full geometry with standard functions
        !%     We restore this to the regular depth before computing IC
        !elemR(thisP,er_Temp01) = elemR(thisP,er_Depth)
        initialDepth(thisP) = depth(thisP)
        

        ! print *, 'in ',trim(subroutine_name)
        ! print *, 'geometrytype ',geometryType,trim(reverseKey(geometryType))    

        select case (geometryType)

        case (lIrregular)
             
            ! print *, 'lIrregular'
            ! print *, 'thisLink ',thisLink

            !% --- transect index for this link
            link_tidx => link%I(thisLink,li_transect_idx)

            ! print *, 'link_tidx ',link_tidx

            ! print *, 'thisP ',thisP
            
            !% --- assign non-table transect data
            elemI(thisP,ei_geometryType)  = irregular
            elemI(thisP,ei_link_transect_idx)  = link_tidx

            ! print *, 'tr_widthMax ', tr_widthMax
            ! print *, 'link_tidx   ', link_tidx
            ! print *, size(link%transectR,1), size(link%transectR,2)

            ! print *, 'link%transectR(link_tidx,tr_widthMax)',link%transectR(link_tidx,tr_widthMax)
            
            !% --- independent data
            elemR(thisP,er_BreadthMax)          = link%transectR(link_tidx,tr_widthMax)

            elemR(thisP,er_AreaBelowBreadthMax) = link%transectR(link_tidx,tr_areaBelowBreadthMax)

            !% --- note, do not apply the full depth limiter function to transects!
            elemR(thisP,er_FullDepth)           = link%transectR(link_tidx,tr_depthFull)
 
            elemR(thisP,er_FullArea)            = link%transectR(link_tidx,tr_areaFull)

            elemR(thisP,er_FullTopwidth)        = link%transectR(link_tidx,tr_widthFull)

            elemR(thisP,er_ZbreadthMax)         = link%transectR(link_tidx,tr_depthAtBreadthMax) + elemR(thisP,er_Zbottom)

            elemR(thisP,er_FullHydRadius)       = link%transectR(link_tidx,tr_hydRadiusFull)

            !% --- full conditions
            elemR(thisP,er_FullPerimeter) = llgeo_perimeter_from_hydradius_and_area_pure &
                                                (thisP, fullhydradius(thisP), fullarea(thisP))  

            !% --- dependent data
            elemR(thisP,er_Zcrown)        = elemR(thisP,er_Zbottom)  + elemR(thisP,er_FullDepth)
            elemR(thisP,er_FullVolume)    = elemR(thisP,er_FullArea) * elemR(thisP,er_Length)
            
            !% ---NOTE the IC data for area, volume, etc cannot be initialized until the transect tables are setup, which is
            !%     delayed until after the JB are initialized.

        case (lParabolic)
            ! print *, 'lParabolic'
            
            elemI(thisP,ei_geometryType) = parabolic

            !% --- independent data
            elemSGR(thisP,esgr_Parabolic_Breadth)   = link%R(thisLink,lr_BreadthScale)
            elemSGR(thisP,esgr_Parabolic_Radius)    = elemSGR(thisP,esgr_Parabolic_Breadth) / twoR / sqrt(link%R(thisLink,lr_FullDepth))
            elemR(thisP,er_FullDepth)               = link%R(thisLink,lr_FullDepth)
            elemR(thisP,er_BreadthMax)              = link%R(thisLink,lr_BreadthScale)

            !% --- error checking
            if ((link%R(thisLink,lr_FullDepth) .le. zeroR) .or. &
                (link%R(thisLink,lr_BreadthScale) .le. zeroR)) then 
                print *, 'USER CONFIGURATION ERROR'
                print *, 'Parabolic cross section has zero specified for Full Height or Top Width'
                print *, 'Problem with link # ',thisLink
                print *, 'which is named ',trim(link%Names(thisLink)%str)
                call util_crashpoint(698704)
            end if
  
            !% --- full conditions
            elemR(thisP,er_FullArea)      = llgeo_parabolic_area_from_depth_pure &
                                                (thisP, fulldepth(thisP))

            elemR(thisP,er_FullPerimeter) = llgeo_parabolic_perimeter_from_depth_pure &
                                                (thisP, fulldepth(thisP))

            elemR(thisP,er_FullTopwidth)  = llgeo_parabolic_topwidth_from_depth_pure &
                                                (thisP, fulldepth(thisP))

            elemR(thisP,er_FullHydRadius) = llgeo_hydradius_from_area_and_perimeter_pure &
                                                (thisP, fullarea(thisP), fullperimeter(thisP))
            
            !% --- dependent  
            elemR(thisP,er_AreaBelowBreadthMax)     = elemR(thisP,er_FullArea) 
            elemR(thisP,er_ZbreadthMax)             = elemR(thisP,er_FullDepth) + elemR(thisP,er_Zbottom)
            elemR(thisP,er_Zcrown)                  = elemR(thisP,er_Zbottom)   + elemR(thisP,er_FullDepth)
            elemR(thisP,er_FullVolume)              = elemR(thisP,er_FullArea)  * elemR(thisP,er_Length)
             

            !% --- store IC data
            elemR(thisP,er_Perimeter)     = llgeo_parabolic_perimeter_from_depth_pure (thisP, depth(thisP))
            elemR(thisP,er_Topwidth)      = llgeo_parabolic_topwidth_from_depth_pure (thisP, depth(thisP))
            elemR(thisP,er_Area)          = llgeo_parabolic_area_from_depth_pure(thisP, depth(thisP))
            elemR(thisP,er_Area_N0)       = elemR(thisP,er_Area)
            elemR(thisP,er_Area_N1)       = elemR(thisP,er_Area)
            elemR(thisP,er_Volume)        = elemR(thisP,er_Area) * elemR(thisP,er_Length)
            elemR(thisP,er_Volume_N0)     = elemR(thisP,er_Volume)
            elemR(thisP,er_Volume_N1)     = elemR(thisP,er_Volume)

            where (elemR(thisP,er_Perimeter) > zeroR) 
                elemR(thisP,er_HydRadius) = elemR(thisP,er_Area) / elemR(thisP,er_Perimeter)
            elsewhere
                elemR(thisP,er_HydRadius) = zeroR
            endwhere

        case (lPower_function)
            ! print *, 'lPower_function'

            print *, 'CODE/CONFIGURATION ERROR: power function cross-sections not supported in SWMM5+'
            call util_crashpoint(4589723)

        case (lRectangular)
            ! print *, 'lRectangular'
            ! print *, 'thisP',thisP

            elemI(thisP,ei_geometryType) = rectangular

            !% --- independent data
            elemSGR(thisP,esgr_Rectangular_Breadth) = link%R(thisLink,lr_BreadthScale)
            elemR(thisP,er_Breadthmax)              = link%R(thisLink,lr_BreadthScale)
            elemR(thisP,er_FullDepth)               = init_IC_limited_fulldepth(link%R(thisLink,lr_FullDepth),thisLink)

            !% --- error checking
            if ((link%R(thisLink,lr_FullDepth) .le. zeroR) .or. &
                (link%R(thisLink,lr_BreadthScale) .le. zeroR)) then 
                print *, 'USER CONFIGURATION ERROR'
                print *, 'Rectangular open cross section has zero specified for Full Height or Top Width'
                print *, 'Problem with link # ',thisLink
                print *, 'which is named ',trim(link%Names(thisLink)%str)
                call util_crashpoint(6987041)
            end if

            !% --- custom functions using temporary store
            elemR(thisP,er_FullArea)      = llgeo_rectangular_area_from_depth_pure  &
                                            (thisP, fulldepth(thisP))

            elemR(thisP,er_FullPerimeter) = llgeo_rectangular_perimeter_from_depth_pure &
                                            (thisP, fulldepth(thisP))

            elemR(thisP,er_FullTopwidth)   = llgeo_rectangular_topwidth_from_depth_pure &
                                            (thisP, fulldepth(thisP))

            ! do ii=1,size(thisP)
            !     mm = thisP(ii)
            !     ! elemR(mm,er_FullArea)      = llgeo_rectangular_area_from_depth_pure &
            !     !                                     (nn, fulldepth(mm))

            !     elemR(mm,er_FullPerimeter) = llgeo_rectangular_closed_perimeter_from_depth_singular &
            !                                         (mm, fulldepth(mm))

            !     elemR(mm,er_FullTopwidth)  = llgeo_rectangular_closed_topwidth_from_depth_singular &
            !                                         (mm, fulldepth(mm))
            ! end do

            !% --- full conditions
            !elemR(thisP,er_FullHydDepth)  = llgeo_hyddepth_from_area_and_topwidth_pure &
            !                                    (thisP, fullarea(thisP), fulltopwidth(thisP))

            elemR(thisP,er_FullHydRadius) = llgeo_hydradius_from_area_and_perimeter_pure &
                                                (thisP, fullarea(thisP), fullperimeter(thisP))

            !elemR(thisP,er_FullEllDepth)       = llgeo_FullEll_pure(thisP)  

            !% --- dependent data
            elemR(thisP,er_BreadthMax)              = elemR(thisP,er_FullTopwidth)
            elemR(thisP,er_AreaBelowBreadthMax)     = elemR(thisP,er_FullArea)
            elemR(thisP,er_ZbreadthMax)             = elemR(thisP,er_FullDepth) + elemR(thisP,er_Zbottom)
            elemR(thisP,er_Zcrown)                  = elemR(thisP,er_Zbottom)   + elemR(thisP,er_FullDepth)
            elemR(thisP,er_FullVolume)              = elemR(thisP,er_FullArea)  * elemR(thisP,er_Length)       

            !% --- store IC data
            elemR(thisP,er_Perimeter)     = llgeo_rectangular_perimeter_from_depth_pure (thisP, depth(thisP))
            elemR(thisP,er_Topwidth)      = llgeo_rectangular_topwidth_from_depth_pure (thisP, depth(thisP))
            elemR(thisP,er_Area)          = llgeo_rectangular_area_from_depth_pure(thisP,depth(thisP))
            elemR(thisP,er_Area_N0)       = elemR(thisP,er_Area)
            elemR(thisP,er_Area_N1)       = elemR(thisP,er_Area)
            elemR(thisP,er_Volume)        = elemR(thisP,er_Area) * elemR(thisP,er_Length)
            elemR(thisP,er_Volume_N0)     = elemR(thisP,er_Volume)
            elemR(thisP,er_Volume_N1)     = elemR(thisP,er_Volume)

            where (elemR(thisP,er_Perimeter) > zeroR) 
                elemR(thisP,er_HydRadius) = elemR(thisP,er_Area) / elemR(thisP,er_Perimeter)
            elsewhere
                elemR(thisP,er_HydRadius) = zeroR
            endwhere

            ! print *, 'thisP ',thisP
            ! print *, 'full area ',elemR(thisP,er_FullArea)
            ! print *, 'full Depth ',elemR(thisP,er_FullDepth)
            ! print *, 'rectbreadth ',elemSGR(thisP,esgr_Rectangular_Breadth)
            ! stop 298374

        case (lTrapezoidal)
             !print *, 'lTrapezoidal'

            elemI(thisP,ei_geometryType) = trapezoidal

            !% --- independent data
            elemSGR(thisP,esgr_Trapezoidal_Breadth)    = link%R(thisLink,lr_BreadthScale)
            elemSGR(thisP,esgr_Trapezoidal_LeftSlope)  = link%R(thisLink,lr_LeftSlope)
            elemSGR(thisP,esgr_Trapezoidal_RightSlope) = link%R(thisLink,lr_RightSlope)
            elemR(thisP,er_FullDepth)                  = init_IC_limited_fulldepth(link%R(thisLink,lr_FullDepth),thisLink)

            ! print *, 'link trapezoidal breadth scale ',link%R(thisLink,lr_BreadthScale)
            ! print *, 'left slope ',link%R(thisLink,lr_LeftSlope)
            ! print *, 'rightslope ',link%R(thisLink,lr_RightSlope)
            ! print *, 'here in init_IC_get_channel_geometry'
            ! print *, thisP, elemR(thisP,er_FullDepth), elemSGR(thisP,esgr_Trapezoidal_Breadth)

            
            !% --- error checking
            if ((link%R(thisLink,lr_FullDepth)    .le. zeroR) .or. &
                (link%R(thisLink,lr_LeftSlope)    .le. zeroR) .or. &
                (link%R(thisLink,lr_RightSlope)   .le. zeroR) .or. &
                (link%R(thisLink,lr_BreadthScale) .le. zeroR)) then 
                print *, 'USER CONFIGURATION ERROR'
                print *, 'Trapezoidal open cross section has zero specified for Full Height,'
                print *, 'Base Width, Left Slope, or Right Slope. Note that a base width of '
                print *, 'zero should use a triangular cross section. Left/Right slopes of '
                print *, 'zero should be rectangular cross section.'
                print *, 'Problem with link # ',thisLink
                print *, 'which is named ',trim(link%Names(thisLink)%str)
                print *, 'FullDepth ',link%R(thisLink,lr_FullDepth) 
                print *, 'LeftSlope ',link%R(thisLink,lr_LeftSlope)
                print *, 'RightSlope',link%R(thisLink,lr_RightSlope)
                print *, 'Breadthscale ',link%R(thisLink,lr_BreadthScale)
                call util_crashpoint(6987042)
            end if

            !% --- full conditions
            elemR(thisP,er_FullArea)      = llgeo_trapezoidal_area_from_depth_pure &
                                                (thisP, fulldepth(thisP))

            elemR(thisP,er_FullPerimeter) = llgeo_trapezoidal_perimeter_from_depth_pure &
                                                (thisP, fulldepth(thisP))

            elemR(thisP,er_FullTopwidth)  = llgeo_trapezoidal_topwidth_from_depth_pure &
                                                (thisP, fulldepth(thisP))

            elemR(thisP,er_FullHydRadius) = llgeo_hydradius_from_area_and_perimeter_pure &
                                                (thisP, fullarea(thisP), fullperimeter(thisP))
            
            !% --- dependent data
            elemR(thisP,er_BreadthMax)              = elemR(thisP,er_FullTopwidth)
            elemR(thisP,er_AreaBelowBreadthMax)     = elemR(thisP,er_FullArea)
            elemR(thisP,er_ZbreadthMax)             = elemR(thisP,er_FullDepth) + elemR(thisP,er_Zbottom)
            elemR(thisP,er_Zcrown)                  = elemR(thisP,er_Zbottom)   + elemR(thisP,er_FullDepth)
            elemR(thisP,er_FullVolume)              = elemR(thisP,er_FullArea)  * elemR(thisP,er_Length)
            
        
            !% --- store IC data
            elemR(thisP,er_Perimeter)    = llgeo_trapezoidal_perimeter_from_depth_pure (thisP, depth(thisP))
            elemR(thisP,er_Topwidth)     = llgeo_trapezoidal_topwidth_from_depth_pure (thisP, depth(thisP))
            elemR(thisP,er_Area)         = llgeo_trapezoidal_area_from_depth_pure (thisP, depth(thisP))
            elemR(thisP,er_Area_N0)      = elemR(thisP,er_Area)
            elemR(thisP,er_Area_N1)      = elemR(thisP,er_Area)
            elemR(thisP,er_Volume)       = elemR(thisP,er_Area) * elemR(thisP,er_Length)
            elemR(thisP,er_Volume_N0)    = elemR(thisP,er_Volume)
            elemR(thisP,er_Volume_N1)    = elemR(thisP,er_Volume)     
            
            where (elemR(thisP,er_Perimeter) > zeroR) 
                elemR(thisP,er_HydRadius) = elemR(thisP,er_Area) / elemR(thisP,er_Perimeter)
            elsewhere
                elemR(thisP,er_HydRadius) = zeroR
            endwhere

            ! print *, 'trapezoidal link'
            ! print *, 'thisP ',thisP
            ! print *, 'this depth ',depth(thisP)
            ! print *, 'breadth ',elemSGR(thisP,esgr_Trapezoidal_Breadth)
            ! print *, 'thisTopwidth ',elemR(thisP,er_Topwidth)
            ! stop 59874


            
        case (lTriangular)
            ! print *, 'lTriangular'

            elemI(thisP,ei_geometryType) = triangular

            !% --- independent data
            elemSGR(thisP,esgr_Triangular_TopBreadth)  = link%R(thisLink,lr_BreadthScale)
            elemR(thisP,er_FullDepth)                  = init_IC_limited_fulldepth(link%R(thisLink,lr_FullDepth),thisLink)
            elemR(thisP,er_BreadthMax)                 = link%R(thisLink,lr_BreadthScale)
            elemSGR(thisP,esgr_Triangular_Slope)       = elemSGR(thisP,esgr_Triangular_TopBreadth) &
                                                         / (twoR * elemR(thisP,er_FullDepth))

            ! print *, 'triangular topbreadth ',   elemSGR(thisP,esgr_Triangular_TopBreadth), &
            ! 2.d0 * elemSGR(thisP,esgr_Triangular_Slope) * elemR(thisP,er_FullDepth)

            !% --- error checking
            if ((link%R(thisLink,lr_FullDepth)    .le. zeroR) .or. &
                (link%R(thisLink,lr_BreadthScale) .le. zeroR)) then 
                print *, 'USER CONFIGURATION ERROR'
                print *, 'Triangular open cross section has zero specified for Full Height or Top Width,'
                print *, 'Problem with link # ',thisLink
                print *, 'which is named ',trim(link%Names(thisLink)%str)
                call util_crashpoint(6987043)
            end if
                                                                                      

            !% --- full conditions
            elemR(thisP,er_FullArea)      = llgeo_triangular_area_from_depth_pure &
                                                (thisP, fulldepth(thisP))

            elemR(thisP,er_FullPerimeter) = llgeo_triangular_perimeter_from_depth_pure &
                                                (thisP, fulldepth(thisP))
                                                
            elemR(thisP,er_FullTopwidth)  = llgeo_triangular_topwidth_from_depth_pure &
                                                (thisP, fulldepth(thisP))
  
            !elemR(thisP,er_FullHydDepth)  = llgeo_hyddepth_from_area_and_topwidth_pure &
            !                                    (thisP, fullarea(thisP), fulltopwidth(thisP))

            elemR(thisP,er_FullHydRadius) = llgeo_hydradius_from_area_and_perimeter_pure &
                                                (thisP, fullarea(thisP), fullperimeter(thisP))

            !elemR(thisP,er_FullEllDepth)       = llgeo_FullEll_pure(thisP)
                                         
            !% --- dependent data
            
            elemR(thisP,er_AreaBelowBreadthMax)     = elemR(thisP,er_FullArea)
            elemR(thisP,er_ZbreadthMax)             = elemR(thisP,er_FullDepth) + elemR(thisP,er_Zbottom)
            elemR(thisP,er_Zcrown)                  = elemR(thisP,er_Zbottom)   + elemR(thisP,er_FullDepth)
            elemR(thisP,er_FullVolume)              = elemR(thisP,er_FullArea)  * elemR(thisP,er_Length)
            
            !% store IC data
            elemR(thisP,er_Perimeter)    = llgeo_triangular_perimeter_from_depth_pure (thisP, depth(thisP))
            elemR(thisP,er_Topwidth)     = llgeo_triangular_topwidth_from_depth_pure (thisP, depth(thisP))
            elemR(thisP,er_Area)         = llgeo_triangular_area_from_depth_pure(thisP, depth(thisP)) 
            elemR(thisP,er_Area_N0)      = elemR(thisP,er_Area)
            elemR(thisP,er_Area_N1)      = elemR(thisP,er_Area)
            elemR(thisP,er_Volume)       = elemR(thisP,er_Area) * elemR(thisP,er_Length)
            elemR(thisP,er_Volume_N0)    = elemR(thisP,er_Volume)
            elemR(thisP,er_Volume_N1)    = elemR(thisP,er_Volume)

            where (elemR(thisP,er_Perimeter) > zeroR) 
                elemR(thisP,er_HydRadius) = elemR(thisP,er_Area) / elemR(thisP,er_Perimeter)
            elsewhere
                elemR(thisP,er_HydRadius) = zeroR
            endwhere

        case default

            print *, 'In, ', subroutine_name
            print *, 'CODE ERROR -- geometry type unknown for # ',geometryType
            print *, 'which has key ',trim(reverseKey(geometryType))
            call util_crashpoint(98734)

        end select
    

        !% -- set the face values for the crown (full depth)
        faceR(fUp(thisP),fr_Zcrown_d) = faceR(fUp(thisP),fr_Zbottom) + elemR(thisP,er_FullDepth)
        faceR(fDn(thisP),fr_Zcrown_u) = faceR(fDn(thisP),fr_Zbottom) + elemR(thisP,er_FullDepth)

        !% --- reset the temporary space
        !%     Note, real must be first as int is used for thisP
        elemR(thisP,er_Temp01) = zeroR
        elemI(thisP,ei_Temp01) = nullvalueI

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine init_IC_get_channel_geometry
!
!==========================================================================
!==========================================================================
!
    subroutine init_IC_get_conduit_geometry (thisLink)
        !%-----------------------------------------------------------------
        !% Description:
        !% get the geometry data for closed conduits
        !% NOTE DepthBelowMaxBreadth is taken at the highest depth, or
        !% just slightly below that, where the max breadth occurs; e.g.
        !% for the basket handle the max breadth occurs between 0.2 and
        !% 0.28 of the normalized depth, so we use 0.27 so that 
        !% lookup tables fall in between two values with max breadth
        !%
        !% NOTE if geo_common_initialize is NOT called for a type of
        !% closed-conduit geometry, then the geometry MUST separately
        !% call slot_initialize
        !%-----------------------------------------------------------------
        !% Declarations
            integer :: ii, mm, Npack
            integer, intent(in) :: thisLink
            integer, pointer    :: geometryType, eIdx(:), thisP(:)
            integer, pointer    :: fUp(:), fDn(:)
            real(8), pointer    :: fullDepth(:), breadthMax(:), fullArea(:)
            real(8), pointer    :: depth(:), fullHydRadius(:)
            real(8), pointer    :: pi, topwidthDepth
            real(8)             :: bottomHydRadius, dummyA(1)
            character(64) :: subroutine_name = 'init_IC_get_conduit_geometry'
        !%-----------------------------------------------------------------
        !% Preliminaries
            !% --- pack the elements for this link in temporary array
            Npack = count(elemI(:,ei_link_Gidx_BIPquick) == thisLink)
            if (Npack < 1) return

            if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"   
        !%------------------------------------------------------------------
        !% Aliases
            !% pointer to geometry type
            geometryType => link%I(thisLink,li_geometry)
            pi => setting%Constant%pi

            !% --- pointer to element indexes
            eIdx         => elemI(:,ei_Lidx)

            !% --- pointer to face indexes
            fUp         => elemI(:,ei_Mface_uL)
            fDn         => elemI(:,ei_Mface_dL)

            !% --- pointer to full depth
            depth         => elemR(:,er_Depth)
            fullDepth     => elemR(:,er_FullDepth)
            breadthMax    => elemR(:,er_BreadthMax)
            fullArea      => elemR(:,er_FullArea)            
            fullHydRadius => elemR(:,er_FullHydRadius)
            !tempDepth     => elemR(:,er_Temp02)

            elemI(1:Npack,ei_Temp01) = pack(eIdx,elemI(:,ei_link_Gidx_BIPquick) == thisLink)
            thisP => elemI(1:Npack,ei_Temp01)
        !%------------------------------------------------------------------
        !% --- independent common data
        elemR(thisP,er_FullDepth)     = link%R(thisLink,lr_FullDepth)
        elemR(thisP,er_FullArea)      = link%R(thisLink,lr_FullArea)
        elemR(thisP,er_FullHydRadius) = link%R(thisLink,lr_FullHydRadius)

        !% -- set the face values for the crown
        faceR(fUp(thisP),fr_Zcrown_d) = faceR(fUp(thisP),fr_Zbottom) + elemR(thisP,er_FullDepth)
        faceR(fDn(thisP),fr_Zcrown_u) = faceR(fDn(thisP),fr_Zbottom) + elemR(thisP,er_FullDepth)

        ! print *, ' '
        ! print *, 'in ',trim(subroutine_name)
        ! print *, 'geometry type ',geometryType,trim(reverseKey(geometryType))
        ! print *, ' '

        select case (geometryType)

        case (lArch)  !% TABULAR

            elemI(thisP,ei_geometryType) = arch

            !% --- independent custom data
            elemR(thisP,er_BreadthMax)        = link%R(thisLink,lr_BreadthScale) 
            elemR(thisP,er_DepthAtBreadthMax) = 0.28d0 * elemR(thisP,er_FullDepth)

            !% --- error checking
            if ((link%R(thisLink,lr_FullDepth)    .le. zeroR) .or. &
                (link%R(thisLink,lr_BreadthScale) .le. zeroR)) then 
                print *, 'USER CONFIGURATION ERROR'
                print *, 'Arch cross section has zero specified for Full Height or Top Width,'
                print *, 'Problem with link # ',thisLink
                print *, 'which is named ',trim(link%Names(thisLink)%str)
                call util_crashpoint(6987044)
            end if
  
            call geo_common_initialize (thisP, arch, AArch, TArch, RArch, dummyA)

        case (lBasket_handle) !% TABULAR

            elemI(thisP,ei_geometryType) = basket_handle

            !% --- independent custom data
            elemR(thisP,er_BreadthMax)        = link%R(thisLink,lr_BreadthScale) 
            elemR(thisP,er_DepthAtBreadthMax)  = 0.27d0 * elemR(thisP,er_FullDepth)

            !% --- error checking
            if ((link%R(thisLink,lr_BreadthScale) .le. zeroR)) then 
                print *, 'USER CONFIGURATION ERROR'
                print *, 'BasketHandle cross section has zero specified for Full Height'
                print *, 'Problem with link # ',thisLink
                print *, 'which is named ',trim(link%Names(thisLink)%str)
                call util_crashpoint(6987045)
            end if

            call geo_common_initialize (thisP, basket_handle, ABasketHandle, TBasketHandle, RBasketHandle, dummyA)
    
        case (lCatenary)  !% TABULAR

            elemI(thisP,ei_geometryType) = catenary

            !% --- independent custom data
            elemR(thisP,er_BreadthMax)        = link%R(thisLink,lr_BreadthScale) 
            elemR(thisP,er_DepthAtBreadthMax) = 0.29d0 * elemR(thisP,er_FullDepth)

            !% --- error checking
            if ((link%R(thisLink,lr_FullDepth)    .le. zeroR) ) then 
                print *, 'USER CONFIGURATION ERROR'
                print *, 'Catenary cross section has zero specified for Full Height ,'
                print *, 'Problem with link # ',thisLink
                print *, 'which is named ',trim(link%Names(thisLink)%str)
                call util_crashpoint(6987046)
            end if

            call geo_common_initialize (thisP, catenary, ACatenary, TCatenary, dummyA, SCatenary)

        case (lCircular,lForce_main) !% TABULAR

            !% --- force mains are required to be circular pipe
            elemI(thisP,ei_geometryType)    = circular

            !% --- Get data for force main
            !%     NotethisP all force mains are circular pipes
            if (setting%Solver%ForceMain%AllowForceMainTF) then
                if (geometryType == lForce_main) then 
                    where (elemI(thisP,ei_link_Gidx_BIPquick) == thisLink)
                        elemYN(thisP,eYN_isForceMain)      = .TRUE.
                        elemSI(thisP,esi_Conduit_Forcemain_Method) = setting%SWMMinput%ForceMainEquationType
                        elemSR(thisP,esr_Conduit_ForceMain_Coef)   = link%R(thislink,lr_ForceMain_Coef)
                    endwhere
                endif
            else
                !% -- if global AllowForceMainTF is false, then
                !%    make sure FM from the SWMM input are set to
                !%    non-force-main.
                if (geometryType == lForce_main) then 
                    where (elemI(thisP,ei_link_Gidx_BIPquick) == thisLink)
                        elemSI(thisP,esi_Conduit_Forcemain_Method) = nullvalueI
                        elemYN(thisP,eYN_isForceMain)      = .FALSE.
                        elemSR(thisP,esr_Conduit_ForceMain_Coef)   = nullvalueR
                    endwhere
                end if
            end if
                        
            !% --- independent custom data
            elemSGR(thisP,esgr_Circular_Diameter) = link%R(thisLink,lr_BreadthScale)
            elemSGR(thisP,esgr_Circular_Radius)   = link%R(thisLink,lr_BreadthScale) / twoR
            elemR(thisP,er_BreadthMax)            = elemSGR(thisP,esgr_Circular_Diameter)
            elemR(thisP,er_DepthAtBreadthMax)     = 0.5d0 * elemR(thisP,er_FullDepth)

            !% --- error checking
            if ((link%R(thisLink,lr_BreadthScale)    .le. zeroR)) then 
                print *, 'USER CONFIGURATION ERROR'
                print *, 'Circular cross section has zero specified for Diameter ,'
                print *, 'Problem with link # ',thisLink
                print *, 'which is named ',trim(link%Names(thisLink)%str)
                call util_crashpoint(6987047)
            end if

                ! print *, 'calling geo_common_initialize for circular'
            call geo_common_initialize (thisP, circular, ACirc, TCirc, RCirc, dummyA)   
         
        case (lCustom)  !% TABULAR
            print *, 'CODE/CONFIGURATION ERROR: Custom conduit cross-sections not supported in SWMM5+'
            call util_crashpoint(7798723)

        case (lEggshaped)  !% TABULAR

            elemI(thisP,ei_geometryType) = eggshaped

            !% --- independent custom data
            elemR(thisP,er_BreadthMax)        = link%R(thisLink,lr_BreadthScale)
            elemR(thisP,er_DepthAtBreadthMax) = 0.64d0 * elemR(thisP,er_FullDepth)

            !% --- error checking
            if ((link%R(thisLink,lr_BreadthScale)    .le. zeroR)) then 
                print *, 'USER CONFIGURATION ERROR'
                print *, 'Eggshaped cross section has zero specified for FullHeight ,'
                print *, 'Problem with link # ',thisLink
                print *, 'which is named ',trim(link%Names(thisLink)%str)
                call util_crashpoint(6987048)
            end if

            call geo_common_initialize (thisP, eggshaped, AEgg, TEgg, REgg, dummyA)
        
        case (lFilled_circular)  !% ANALYTICAL
            !% --- note, Zbottom is always the bottom of the filled section (not the top of it)

            elemI(thisP,ei_geometryType) = filled_circular


            !% HACK -- THIS ASSUMES THAT INPUT VALUES OF FULLDEPTH IS FOR PIPE WITHOUT SEDIMENT

            !% --- independent data

            !% --- get the sediment depth
            elemR(thisP,er_SedimentDepth) = link%R(thisLink,lr_BottomDepth)

            !% --- error checking
            if ((link%R(thisLink,lr_BreadthScale)    .le. zeroR) .or. &
                (link%R(thisLink,lr_BottomDepth)      <  zeroR)) then 
                print *, 'USER CONFIGURATION ERROR'
                print *, 'Filled circular cross section has zero specified for FullHeight '
                print *, 'or less than zero for sediment depth,'
                print *, 'Problem with link # ',thisLink
                print *, 'which is named ',trim(link%Names(thisLink)%str)
                call util_crashpoint(6987049)
            end if

            !% --- reset the depth previously computed from nodes (without sediment)
            elemR(thisP,er_Depth) = elemR(thisP,er_Depth) - elemR(thisP,er_SedimentDepth)

            elemSGR(thisP,esgr_Filled_Circular_TotalPipeDiameter)    &
                 = link%R(thisLink,lr_FullDepth) + elemR(thisP,er_SedimentDepth)    !% HACK -- check what link full depth means

            elemSGR(thisP,esgR_Filled_Circular_TotalPipeArea)        &
                = (onefourthR * pi) * (elemSGR(thisP,esgr_Filled_Circular_TotalPipeDiameter)**2)

            elemSGR(thisP,esgR_Filled_Circular_TotalPipePerimeter)   &
                = pi * elemSGR(thisP,esgr_Filled_Circular_TotalPipeDiameter)    

            elemSGR(thisP,esgr_Filled_Circular_TotalPipeHydRadius)     &
                =   elemSGR(thisP,esgR_Filled_Circular_TotalPipeArea)  &
                  / elemSGR(thisP,esgR_Filled_Circular_TotalPipePerimeter)


            !% FOR INITIAL FILLED AREA CALCULATION, RESET THE FULL DEPTH TO TOTAL DIA OF THE PIPE
            elemR(thisP,er_BreadthMax) = elemSGR(thisP,esgr_Filled_Circular_TotalPipeDiameter)

            elemR(thisP,er_FullDepth)  =  elemSGR(thisP,esgr_Filled_Circular_TotalPipeDiameter) 

            do ii=1,size(thisP)
                mm = thisP(ii)
                if (elemR(mm,er_SedimentDepth) >= setting%ZeroValue%Depth) then

                    elemSGR(mm,esgr_Filled_Circular_bottomArea)               &
                        = llgeo_tabular_from_depth_singular                   &
                            (mm, elemR(mm,er_SedimentDepth), elemSGR(mm,esgR_Filled_Circular_TotalPipeArea),    &
                            setting%ZeroValue%Depth, zeroR, ACirc )

                    elemSGR(mm,esgr_Filled_Circular_bottomTopwidth)           &
                        = llgeo_tabular_from_depth_singular                   &
                            (mm, elemR(mm,er_SedimentDepth), breadthMax(mm),  &
                            setting%ZeroValue%Depth, zeroR, TCirc )

                    bottomHydRadius                                             &
                        = llgeo_tabular_from_depth_singular                     &
                            (mm, elemR(mm,er_SedimentDepth), fullHydRadius(mm), &
                            setting%ZeroValue%Depth, zeroR, RCirc )   

                    if (bottomHydRadius <= setting%ZeroValue%Depth) then
                        !% -- near zero hydraulic radius
                        elemSGR(mm,esgr_Filled_Circular_bottomPerimeter) = zeroR
                    else
                        elemSGR(mm,esgr_Filled_Circular_bottomPerimeter) &
                            = elemSGR(mm,esgr_Filled_Circular_bottomArea) / bottomHydRadius
                    end if
                else 
                    !% --- near zero sediment depths
                    !%     the setting%ZeroValues%... are not yet assigned.
                    elemSGR(mm,esgr_Filled_Circular_bottomArea)      = zeroR
                    elemSGR(mm,esgr_Filled_Circular_bottomTopwidth)  = zeroR
                    elemSGR(mm,esgr_Filled_Circular_bottomPerimeter) = zeroR
                end if


            end do

             !% --- for consistency in other uses, elemR store values for the flow section only
            elemR(thisP,er_FullDepth)     =  elemSGR(thisP,esgr_Filled_Circular_TotalPipeDiameter)   &
                                             - elemR(thisP,er_SedimentDepth)  

            elemR(thisP,er_FullArea)      =  elemSGR(thisP,esgR_Filled_Circular_TotalPipeArea)       &
                                           - elemSGR(thisP,esgr_Filled_Circular_bottomArea)

            elemR(thisP,er_FullPerimeter) =   elemSGR(thisP,esgr_Filled_Circular_TotalPipePerimeter) &
                                            - elemSGR(thisP,esgr_Filled_Circular_bottomPerimeter)    &
                                            + elemSGR(thisP,esgr_Filled_Circular_bottomTopwidth)

            elemR(thisP,er_FullHydRadius) = elemR(thisP,er_FullArea) / elemR(thisP,er_FullPerimeter) 

            !% --- Location of maximum breadth
            where (elemR(thisP,er_SedimentDepth)  &
                   > (onehalfR * elemSGR(thisP,esgr_Filled_Circular_TotalPipeDiameter)) )

                !% --- solid fill level is above the midpoint, then max breadth for flow is at top of solid fill
                elemR(thisp,er_DepthAtBreadthMax)   = zeroR 
                elemR(thisP,er_BreadthMax)          = elemSGR(thisP,esgr_Filled_Circular_bottomTopwidth)
                elemR(thisP,er_AreaBelowBreadthMax) = zeroR

            elsewhere
                !% --- solid fill level is below the midpoint, then max breadth for flow is at midpoint
                elemR(thisP,er_DepthAtBreadthMax)   =  onehalfR * elemSGR(thisP,esgr_Filled_Circular_TotalPipeDiameter) &
                                                       - elemR(thisP,er_SedimentDepth)

                elemR(thisP,er_BreadthMax)          = elemSGR(thisP,esgr_Filled_Circular_TotalPipeDiameter)

                elemR(thisP,er_AreaBelowBreadthMax) = onehalfR * elemSGR(thisP,esgR_Filled_Circular_TotalPipeArea)  &
                                                               - elemSGR(thisP,esgr_Filled_Circular_bottomArea) 
            endwhere

            call geo_common_initialize (thisP, filled_circular, dummyA, dummyA, dummyA, dummyA)
        
        case (lGothic) !% TABULAR

            elemI(thisP,ei_geometryType) = gothic

            !% --- independent custom data
            elemR(thisP,er_BreadthMax)        = link%R(thisLink,lr_BreadthScale) 
            elemR(thisP,er_DepthAtBreadthMax) = 0.49d0 * elemR(thisP,er_FullDepth)

            !% --- error checking
            if ((link%R(thisLink,lr_BreadthScale)    .le. zeroR)) then 
                print *, 'USER CONFIGURATION ERROR'
                print *, 'Gothic cross section has zero specified for FullHeight '
                print *, 'Problem with link # ',thisLink
                print *, 'which is named ',trim(link%Names(thisLink)%str)
                call util_crashpoint(698701)
            end if

            call geo_common_initialize (thisP, gothic, AGothic, TGothic, dummyA, SGothic)
          
        case (lHoriz_ellipse) !% TABULAR

            elemI(thisP,ei_geometryType) = horiz_ellipse

            !% --- independent custom data
            elemR(thisP,er_BreadthMax)        = link%R(thisLink,lr_BreadthScale) 
            elemR(thisP,er_DepthAtBreadthMax) = 0.5d0 * elemR(thisP,er_FullDepth)

            !% --- error checking
            if ((link%R(thisLink,lr_BreadthScale)    .le. zeroR) .or. &
                (link%R(thisLink,lr_FullDepth)       .le. zeroR) ) then 
                print *, 'USER CONFIGURATION ERROR'
                print *, 'Horiz Ellipse cross section has zero specified for FullHeight or Max width'
                print *, 'Problem with link # ',thisLink
                print *, 'which is named ',trim(link%Names(thisLink)%str)
                call util_crashpoint(698702)
            end if

            call geo_common_initialize (thisP, horiz_ellipse, AHorizEllip, THorizEllip, RHorizEllip, dummyA)

        case (lHorseshoe) !% TABULAR

            elemI(thisP,ei_geometryType) = horseshoe

            !% --- independent custom data
            elemR(thisP,er_BreadthMax)        = link%R(thisLink,lr_BreadthScale) 
            elemR(thisP,er_DepthAtBreadthMax) = 0.5d0 * elemR(thisP,er_FullDepth)

            !% --- error checking
            if ((link%R(thisLink,lr_BreadthScale)    .le. zeroR)) then 
                print *, 'USER CONFIGURATION ERROR'
                print *, 'Horseshoe cross section has zero specified for FullHeight '
                print *, 'Problem with link # ',thisLink
                print *, 'which is named ',trim(link%Names(thisLink)%str)
                call util_crashpoint(698703)
            end if

            call geo_common_initialize (thisP, horseshoe, AHorseShoe, THorseShoe, RHorseShoe, dummyA)

        case (lIrregular) !% ERROR
            print *, 'In ', trim(subroutine_name)
            print *, 'USER ERROR Irregular cross-section geometry not allowed for closed conduits (open-channel only) in SWMM5+'
            call util_crashpoint(4409874)    

        case (lMod_basket)  !% ANALYTICAL

            elemI(thisP,ei_geometryType)  = mod_basket

            !% --- independent custom data
            elemR(  thisP,eR_BreadthMax)             = link%R(thisLink,lr_BreadthScale)
            elemSGR(thisP,esgr_Mod_Basket_Rtop)      = link%R(thisLink,lr_BottomRadius)

            elemSGR(thisP,esgr_Mod_Basket_ThetaTop)  = twoR * asin( onehalfR * elemR(thisP,er_BreadthMax) &
                                                                    / elemSGR(thisP,esgr_Mod_Basket_Rtop) )

            elemSGR(thisP,esgr_Mod_Basket_Ytop)      = elemSGR(thisP,esgr_Mod_Basket_Rtop) &
                                                      * ( oneR - cos( elemSGR(thisP,esgr_Mod_Basket_ThetaTop) / twoR ) )

                
            elemSGR(thisP,esgr_Mod_Basket_Atop)      = onehalfR * ( elemSGR(thisP,esgr_Mod_Basket_Rtop)**2 )       &
                                                                * ( elemSGR(thisP,esgr_Mod_Basket_ThetaTop)        &
                                                                    - sin(elemSGR(thisP,esgr_Mod_Basket_ThetaTop)) &
                                                                  ) 

            elemR(thisP,er_DepthAtBreadthMax)        = elemR(thisP,er_FullDepth) - elemSGR(thisP,esgr_Mod_Basket_Ytop)                                                      

            !% --- error checking
            if ((link%R(thisLink,lr_BreadthScale)    .le. zeroR) .or. &
                (link%R(thisLink,lr_BottomRadius)    .le. zeroR) .or. &
                (link%R(thisLink,lr_FullDepth)       .le. zeroR)  ) then 
                print *, 'USER CONFIGURATION ERROR'
                print *, 'Mod Basket cross section has zero specified for FullHeight, Base width, '
                print *, 'or Top Radius'
                print *, 'Problem with link # ',thisLink
                print *, 'which is named ',trim(link%Names(thisLink)%str)
                call util_crashpoint(698704)
            end if

            call geo_common_initialize (thisP, mod_basket, dummyA, dummyA, dummyA, dummyA)

        case (lRectangular_closed)  !% ANALYTICAL

                elemI(thisP,ei_geometryType) = rectangular_closed

                !% --- independent data
                elemR(thisP,er_BreadthMax)              = link%R(thisLink,lr_BreadthScale)
                elemR(thisP,er_DepthAtBreadthMax)       = onehalfR * elemR(thisP,er_FullDepth)
                elemSGR(thisP,esgr_Rectangular_Breadth) = elemR(thisP,er_BreadthMax) 

                !% --- error checking
                if ((link%R(thisLink,lr_BreadthScale)    .le. zeroR) .or. &
                    (link%R(thisLink,lr_FullDepth)       .le. zeroR)  ) then 
                    print *, 'USER CONFIGURATION ERROR'
                    print *, 'Rectangular Closed cross section has zero specified for FullHeight or Top width, '
                    print *, 'Problem with link # ',thisLink
                    print *, 'which is named ',trim(link%Names(thisLink)%str)
                    call util_crashpoint(698705)
                end if   

                call geo_common_initialize (thisP, rectangular_closed, dummyA, dummyA, dummyA, dummyA)
                    
        case (lRect_round)

            elemI(thisP,ei_geometryType)                       = rect_round

            !% --- independent data
            elemSGR(thisP,esgr_Rectangular_Round_Ybot)   = link%R(thisLink,lr_BottomDepth)
            elemSGR(thisP,esgr_Rectangular_Round_Rbot)   = link%R(thisLink,lr_BottomRadius)
            elemR( thisP,er_BreadthMax)                  = link%R(thisLink,lr_BreadthScale)
            elemR( thisP,er_DepthAtBreadthMax)           = elemSGR(thisP,esgr_Rectangular_Round_Ybot) &
                                                         + elemSGR(thisP,esgr_Rectangular_Round_Rbot)

            !% --- error checking
            if ((link%R(thisLink,lr_BreadthScale)    .le. zeroR) .or. &
                (link%R(thisLink,lr_BottomRadius)    .le. zeroR) .or. &
                (link%R(thisLink,lr_FullDepth)       .le. zeroR)  ) then 
                print *, 'USER CONFIGURATION ERROR'
                print *, 'Rectangular Round cross section has zero specified for FullHeight or Top width, '
                print *, 'or bottom radius.'
                print *, 'Problem with link # ',thisLink
                print *, 'which is named ',trim(link%Names(thisLink)%str)
                call util_crashpoint(698706)
            end if  


            elemSGR(thisP,esgr_Rectangular_Round_ThetaBot)                        &
                    = twoR * asin(                                                &
                                    onehalfR * elemR(thisP,er_BreadthMax)         & 
                                    / elemSGR(thisP,esgr_Rectangular_Round_Rbot)  &
                                    )
            elemSGR(thisP,esgr_Rectangular_Round_Abot)                                   &
                    = onehalfR * (elemSGR(thisP,esgr_Rectangular_Round_Rbot)**2)    & 
                                * (                                                       &
                                    elemSGR(thisP,esgr_Rectangular_Round_ThetaBot)        &
                                    - sin(elemSGR(thisP,esgr_Rectangular_Round_ThetaBot)) &
                                    )

            call geo_common_initialize (thisP, rect_round, dummyA, dummyA, dummyA, dummyA)                

        case (lRect_triang) !% ANALYTICAL

            elemI(thisP,ei_geometryType) = rect_triang

            !% --- independent data
            elemSGR(thisP,esgr_Rectangular_Triangular_BottomDepth)  = link%R(thisLink,lr_BottomDepth)
            elemR(  thisP,er_BreadthMax)                            = link%R(thisLink,lr_BreadthScale)
            elemR(  thisP,er_DepthAtBreadthMax)                     = elemR(thisP,er_FullDepth)  

            !% --- error checking
            if ((link%R(thisLink,lr_BreadthScale)    .le. zeroR) .or. &
                (link%R(thisLink,lr_BottomDepth)    .le. zeroR) .or. &
                (link%R(thisLink,lr_FullDepth)       .le. zeroR)  ) then 
                print *, 'USER CONFIGURATION ERROR'
                print *, 'Rectangular triangular cross section has zero specified for FullHeight or Top width, '
                print *, 'or triangle height.'
                print *, 'Problem with link # ',thisLink
                print *, 'which is named ',trim(link%Names(thisLink)%str)
                call util_crashpoint(698706)
            end if  

            elemSGR(thisP,esgr_Rectangular_Triangular_BottomSlope)  &
                = elemR(thisP,er_BreadthMax)  / (twoR * elemSGR(thisP,esgr_Rectangular_Triangular_BottomDepth))

            elemSGR(thisP,esgr_Rectangular_Triangular_BottomArea)  &
                 = onehalfR * elemSGR(thisP,esgr_Rectangular_Triangular_BottomDepth) &
                            * elemR(thisP,er_BreadthMax)  
                            
            call geo_common_initialize (thisP, rect_triang, dummyA, dummyA, dummyA, dummyA)             
            
        case (lSemi_circular) !% TABULAR

            elemI(thisP,ei_geometryType) = semi_circular

            !% --- independent custom data
            elemR(thisP,er_BreadthMax)        = link%R(thisLink,lr_BreadthScale) 
            elemR(thisP,er_DepthAtBreadthMax) = 0.19d0 * elemR(thisP,er_FullDepth)

            if ((link%R(thisLink,lr_BreadthScale)    .le. zeroR)) then 
                print *, 'USER CONFIGURATION ERROR'
                print *, 'Semi Circular cross section has zero specified for FullHeight '
                print *, 'Problem with link # ',thisLink
                print *, 'which is named ',trim(link%Names(thisLink)%str)   
                call util_crashpoint(698706)
            end if

            call geo_common_initialize (thisP, semi_circular, ASemiCircular, TSemiCircular, dummyA, SSemiCircular)
        

        case (lSemi_elliptical) !% TABULAR

            elemI(thisP,ei_geometryType) = semi_elliptical

            !% --- independent custom data
            elemR(thisP,er_BreadthMax)        = link%R(thisLink,lr_BreadthScale) 
            elemR(thisP,er_DepthAtBreadthMax) = 0.24d0 * elemR(thisP,er_FullDepth)

            call geo_common_initialize (thisP, semi_elliptical, ASemiEllip, TSemiEllip, dummyA, SSemiEllip)
    
            !% --- error checking
            if ((link%R(thisLink,lr_BreadthScale)    .le. zeroR) .or. &
                (link%R(thisLink,lr_FullDepth)       .le. zeroR)) then 
                print *, 'USER CONFIGURATION ERROR'
                print *, 'Semi Elliptical cross section has zero specified for FullHeight or Max WIdth '
                print *, 'Problem with link # ',thisLink
                print *, 'which is named ',trim(link%Names(thisLink)%str)
                call util_crashpoint(698707)
            end if

        case (lVert_ellipse) !% TABULAR

            elemI(thisP,ei_geometryType) = vert_ellipse

            !% --- independent custom data
            elemR(thisP,er_BreadthMax)        = link%R(thisLink,lr_BreadthScale) 
            elemR(thisP,er_DepthAtBreadthMax) = 0.50d0 * elemR(thisP,er_FullDepth)

            if ((link%R(thisLink,lr_BreadthScale)    .le. zeroR) .or. &
                (link%R(thisLink,lr_FullDepth)       .le. zeroR)) then 
                print *, 'USER CONFIGURATION ERROR'
                print *, 'Vertical Elliptical cross section has zero specified for FullHeight or Max Width'
                print *, 'Problem with link # ',thisLink
                print *, 'which is named ',trim(link%Names(thisLink)%str)
                call util_crashpoint(698708)
            end if

            call geo_common_initialize (thisP, vert_ellipse, AVertEllip, TVertEllip, RVertEllip, dummyA)

        case default

            print *, 'In, ', trim(subroutine_name)
            print *, 'CODE ERROR geometry type unknown for # ', geometryType
            if ((geometryType > 0) .and. (geometryType < keys_lastplusone)) then
                print *, 'which has key ',trim(reverseKey(geometryType))
            else
                print *, 'which is not a valid geometry type index'
            end if
            call util_crashpoint(887344)

        end select

        !% --- reset temporary space
        elemI(:,ei_Temp01) = nullvalueI

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_get_conduit_geometry
!
!==========================================================================
!==========================================================================
!
    subroutine init_IC_get_weir_geometry (thisLink)
        !--------------------------------------------------------------------------
        !
        !% get the geometry and other data data for weir links
        !
        !--------------------------------------------------------------------------

            integer, intent(in) :: thisLink
            integer, pointer    :: specificWeirType
            integer, allocatable :: thisPack(:)
            integer :: ii

            character(64) :: subroutine_name = 'init_IC_get_weir_geometry'
        !--------------------------------------------------------------------------
            !if (crashYN) return
            if (setting%Debug%File%initial_condition) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% pointer to specific weir type
        specificWeirType => link%I(thisLink,li_link_sub_type)

        select case (specificWeirType)
            !% --- set up weir specific data
            case (lTrapezoidalWeir)

                where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                    !% integer data
                    elemSI(:,esi_Weir_SpecificType)          = trapezoidal_weir
                    elemSI(:,esi_Weir_GeometryType)          = trapezoidal
                    !% real data
                    elemSR(:,esr_Weir_FullDepth)             = link%R(thisLink,lr_FullDepth)  
                    elemSR(:,esr_Weir_EffectiveFullDepth)    = link%R(thisLink,lr_FullDepth)
                    !% the unit of weir discharge coefficient is (ft^0.5)/sec
                    !% thus converting it to (m^0.5)/sec by multiplying it with 0.5521
                    elemSR(:,esr_Weir_Rectangular)           = link%R(thisLink,lr_DischargeCoeff1) * 0.5521
                    elemSR(:,esr_Weir_Triangular)            = link%R(thisLink,lr_DischargeCoeff2) * 0.5521
                    elemSR(:,esr_Weir_TrapezoidalBreadth)    = link%R(thisLink,lr_BreadthScale)
                    elemSR(:,esr_Weir_TrapezoidalLeftSlope)  = link%R(thisLink,lr_SideSlope)
                    elemSR(:,esr_Weir_TrapezoidalRightSlope) = link%R(thisLink,lr_SideSlope)
                    elemSR(:,esr_Weir_FullArea)              = ( elemSR(:,esr_Weir_TrapezoidalBreadth) &
                                                                + onehalfR  &
                                                                * (  elemSR(:,esr_Weir_TrapezoidalLeftSlope) &
                                                                    + elemSR(:,esr_Weir_TrapezoidalRightSlope) &
                                                                    ) * elemSR(:,esr_Weir_FullDepth) &
                                                                ) * elemSR(:,esr_Weir_FullDepth)
                    elemSR(:,esr_Weir_Zcrest)                = elemR(:,er_Zbottom) + link%R(thisLink,lr_InletOffset)
                    elemSR(:,esr_Weir_Zcrown)                = elemSR(:,esr_Weir_Zcrest) + link%R(thisLink,lr_FullDepth)

                    !% --- default channel geometry (overwritten later by adjacent CC shape)
                    !%     assumes channel is rectangular 
                    elemI(:,ei_geometryType)            = rectangular
                    elemSGR(:,esgr_Rectangular_Breadth) = twoR * (  elemSR(:,esr_Weir_TrapezoidalBreadth) &
                                                               +    elemSR(:,esr_Weir_EffectiveFullDepth) &
                                                                *(  elemSR(:,esr_Weir_TrapezoidalLeftSlope) &
                                                                  + elemSR(:,esr_Weir_TrapezoidalRightSlope)) )
                    elemR(:,er_BreadthMax)              = elemSGR(:,esgr_Rectangular_Breadth)                                       
                    elemR(:,er_FullDepth)               = twoR * max(elemSR(:,esr_Weir_Zcrown) - elemR(:,er_Zbottom), elemSR(:,esr_Weir_FullDepth))  
                endwhere

            case (lSideFlowWeir)

                where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                    !% integer data
                    elemSI(:,esi_Weir_SpecificType)          = side_flow
                    elemSI(:,esi_Weir_GeometryType)          = rectangular
                    elemSI(:,esi_Weir_EndContractions)       = link%I(thisLink,li_weir_EndContractions)
                    !% real data
                    elemSR(:,esr_Weir_EffectiveFullDepth)    = link%R(thisLink,lr_FullDepth)
                    elemSR(:,esr_Weir_FullDepth)             = link%R(thisLink,lr_FullDepth) 
                    !% the unit of weir discharge coefficient is (ft^0.5)/sec
                    !% thus converting it to (m^0.5)/sec by multiplying it with 0.5521
                    elemSR(:,esr_Weir_Rectangular)           = link%R(thisLink,lr_DischargeCoeff1) * 0.5521
                    elemSR(:,esr_Weir_RectangularBreadth)    = link%R(thisLink,lr_BreadthScale)
                    elemSR(:,esr_Weir_FullArea)              = elemSR(:,esr_Weir_RectangularBreadth) * elemSR(:,esr_Weir_FullDepth)
                    elemSR(:,esr_Weir_Zcrest)                = elemR(:,er_Zbottom) + link%R(thisLink,lr_InletOffset)
                    elemSR(:,esr_Weir_Zcrown)                = elemSR(:,esr_Weir_Zcrest) + link%R(thisLink,lr_FullDepth)

                    !% --- default channel geometry (overwritten later by adjacent CC shape)
                    !%     assumes channel is rectangular and uses weir data
                    elemI(:,ei_geometryType)            = rectangular
                    elemSGR(:,esgr_Rectangular_Breadth) = twoR * elemSR(:,esr_Weir_RectangularBreadth) 
                    elemR(:,er_BreadthMax)              = elemSR(:,esr_Weir_RectangularBreadth)                                     
                    elemR(:,er_FullDepth)               = twoR * max(elemSR(:,esr_Weir_Zcrown) - elemR(:,er_Zbottom),elemR(:,er_FullDepth))
                endwhere

            case (lRoadWayWeir)

                where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)  
                    !% integer data
                    elemSI(:,esi_Weir_SpecificType)          = roadway_weir
                    elemSI(:,esi_Weir_GeometryType)          = rectangular
                    elemSI(:,esi_Weir_RoadSurface)           = link%I(thisLink,li_RoadSurface)
                    !% real data
                    elemSR(:,esr_Weir_EffectiveFullDepth)    = link%R(thisLink,lr_FullDepth)
                    elemSR(:,esr_Weir_FullDepth)             = link%R(thisLink,lr_FullDepth) 
                    !% the unit of weir discharge coefficient is (ft^0.5)/sec
                    !% thus converting it to (m^0.5)/sec by multiplying it with 0.5521
                    elemSR(:,esr_Weir_Rectangular)           = link%R(thisLink,lr_DischargeCoeff1) * 0.5521
                    elemSR(:,esr_Weir_RectangularBreadth)    = link%R(thisLink,lr_BreadthScale)
                    elemSR(:,esr_Weir_FullArea)              = elemSR(:,esr_Weir_RectangularBreadth) * elemSR(:,esr_Weir_FullDepth)
                    elemSR(:,esr_Weir_Zcrest)                = elemR(:,er_Zbottom) + link%R(thisLink,lr_InletOffset)
                    elemSR(:,esr_Weir_Zcrown)                = elemSR(:,esr_Weir_Zcrest) + link%R(thisLink,lr_FullDepth)
                    elemSR(:,esr_Wier_RoadWidth)             = link%R(thisLink,lr_RoadWidth)

                    !% --- default channel geometry (overwritten later by adjacent CC shape)
                    !%     assumes channel is rectangular and uses weir data
                    elemI(:,ei_geometryType)            = rectangular
                    elemSGR(:,esgr_Rectangular_Breadth) = twoR * elemSR(:,esr_Weir_RectangularBreadth) 
                    elemR(:,er_BreadthMax)              = elemSR(:,esr_Weir_RectangularBreadth)                                     
                    elemR(:,er_FullDepth)               = twoR * max(elemSR(:,esr_Weir_Zcrown) - elemR(:,er_Zbottom),elemR(:,er_FullDepth))
                end where

            case (lVnotchWeir)

                where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                    !% integer data
                    elemSI(:,esi_Weir_SpecificType)          = vnotch_weir
                    elemSI(:,esi_Weir_GeometryType)          = triangular
                    !% real data
                    elemSR(:,esr_Weir_EffectiveFullDepth)    = link%R(thisLink,lr_FullDepth)
                    elemSR(:,esr_Weir_FullDepth)             = link%R(thisLink,lr_FullDepth)
                    !% the unit of weir discharge coefficient is (ft^0.5)/sec
                    !% thus converting it to (m^0.5)/sec by multiplying it with 0.5521
                    elemSR(:,esr_Weir_Triangular)            = link%R(thisLink,lr_DischargeCoeff1) * 0.5521 
                    elemSR(:,esr_Weir_TriangularSideSlope)   = link%R(thisLink,lr_SideSlope)
                    elemSR(:,esr_Weir_FullArea)              = elemSR(:,esr_Weir_FullDepth) * elemSR(:, esr_Weir_FullDepth) &
                                                                * elemSR(:,esr_Weir_TriangularSideSlope) 
                    elemSR(:,esr_Weir_Zcrest)                = elemR(:,er_Zbottom) + link%R(thisLink,lr_InletOffset)
                    elemSR(:,esr_Weir_Zcrown)                = elemSR(:,esr_Weir_Zcrest) + link%R(thisLink,lr_FullDepth)

                    !% --- default channel geometry (overwritten later by adjacent CC shape)
                    !%     assumes channel is rectangular and twice the breadth of weir and
                    !%     used weir crown as the maximum overflow
                    elemI(:,ei_geometryType)            = rectangular
                    elemSGR(:,esgr_Rectangular_Breadth) = fourR * elemSR(:,esr_Weir_EffectiveFullDepth) &
                                                                * elemSR(:,esr_Weir_TriangularSideSlope)
                    elemR(:,er_BreadthMax)              = elemSGR(:,esgr_Rectangular_Breadth)                                       
                    elemR(:,er_FullDepth)               = twoR * max(elemSR(:,esr_Weir_Zcrown) - elemR(:,er_Zbottom),elemSR(:,esr_Weir_FullDepth)) 
                endwhere

            case (lTransverseWeir)

                where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                    !% integer data
                    elemSI(:,esi_Weir_SpecificType)          = transverse_weir
                    elemSI(:,esi_Weir_GeometryType)          = rectangular
                    elemSI(:,esi_Weir_EndContractions)       = link%I(thisLink,li_weir_EndContractions)
                    !% real data
                    elemSR(:,esr_Weir_EffectiveFullDepth)    = link%R(thisLink,lr_FullDepth)
                    elemSR(:,esr_Weir_FullDepth)             = link%R(thisLink,lr_FullDepth)
                    !% the unit of weir discharge coefficient is (ft^0.5)/sec
                    !% thus converting it to (m^0.5)/sec by multiplying it with 0.5521
                    elemSR(:,esr_Weir_Rectangular)           = link%R(thisLink,lr_DischargeCoeff1) * 0.5521
                    elemSR(:,esr_Weir_RectangularBreadth)    = link%R(thisLink,lr_BreadthScale)
                    elemSR(:,esr_Weir_FullArea)              = elemSR(:,esr_Weir_RectangularBreadth) * elemSR(:,esr_Weir_FullDepth)
                    elemSR(:,esr_Weir_Zcrest)                = elemR(:,er_Zbottom)  + link%R(thisLink,lr_InletOffset)
                    elemSR(:,esr_Weir_Zcrown)                = elemSR(:,esr_Weir_Zcrest) + link%R(thisLink,lr_FullDepth)

                    !% --- default channel geometry (overwritten later by adjacent CC shape)
                    !%     assumes channel is rectangular and uses weir data for channel
                    elemI(:,ei_geometryType)            = rectangular
                    elemSGR(:,esgr_Rectangular_Breadth) = twoR * elemSR(:,esr_Weir_RectangularBreadth) 
                    elemR(:,er_BreadthMax)              = elemSR(:,esr_Weir_RectangularBreadth)                                     
                    elemR(:,er_FullDepth)               = twoR* max(elemSR(:,esr_Weir_Zcrown) - elemR(:,er_Zbottom), elemSR(:,esr_Weir_FullDepth))
                endwhere

            case default

                print *, 'In ', trim(subroutine_name)
                print *, 'CODE ERROR: unknown weir type, ', specificWeirType,'  in network'
                print *, 'which has key ',trim(reverseKey(specificWeirType)) 
                call util_crashpoint(99834)
        end select

        !% --- set minimum crest height as 101% of the zero depth value for all weirs
        !%     this ensures that zero-height weir elements cannot cause flow for zerovalue depths
        thisPack = pack(elemI(:,ei_Lidx),(elemI(:,ei_link_Gidx_BIPquick) == thisLink) ) 
        do ii=1,size(thisPack)
            elemSR(thisPack(ii),esr_Weir_Zcrest) = &
                max( elemSR(thisPack(ii),esr_Weir_Zcrest), elemR(thisPack(ii),er_Zbottom) + setting%ZeroValue%Depth*oneOneHundredthR  )
        end do
        deallocate(thisPack)

        !% --- initialize a default rectangular channel as the background of the weir
        call init_IC_diagnostic_default_geometry (thisLink,weir)

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine init_IC_get_weir_geometry
!
!==========================================================================
!==========================================================================
!
    subroutine init_IC_get_orifice_geometry (thisLink)
        !%--------------------------------------------------------------------------
        !% Description:
        !% get the geometry and other data data for orifice links
        !%--------------------------------------------------------------------------
        !% Declarations
            integer, intent(in)  :: thisLink
            integer, pointer     :: specificOrificeType, OrificeGeometryType
            integer, allocatable :: thisPack(:)
            integer :: ii
            real(8), pointer     :: pi

            character(64) :: subroutine_name = 'init_IC_get_orifice_geometry'
        !--------------------------------------------------------------------------
            if (setting%Debug%File%initial_condition) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        pi => setting%Constant%pi
        !% pointer to specific orifice type
        specificOrificeType => link%I(thisLink,li_link_sub_type)   

        select case (specificOrificeType)
        !% copy orifice specific data
        case (lBottomOrifice)
            where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                !% integer data
                elemSI(:,esi_Orifice_SpecificType)      = bottom_orifice
            endwhere
        case (lSideOrifice)
            where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                !% integer data
                elemSI(:,esi_Orifice_SpecificType)       = side_orifice
            endwhere
        case default
            print *, 'In ', subroutine_name
            print *, 'CODE ERROR: unknown orifice type, ', specificOrificeType,'  in network'
            print *, 'which has key ',trim(reverseKey(specificOrificeType))
            !stop 
            call util_crashpoint(8863411)
            !return
        end select

        !% pointer to specific orifice geometry
        OrificeGeometryType => link%I(thisLink,li_geometry)

        select case (OrificeGeometryType)
            !% copy orifice specific geometry data
        case (lRectangular_closed)  !% brh20211219 added Rect_closed

            where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                !% integer data
                elemSI(:,esi_Orifice_GeometryType)       = rectangular_closed
                !% real data
                elemSR(:,esr_Orifice_FullDepth)          = link%R(thisLink,lr_FullDepth)
                elemSR(:,esr_Orifice_EffectiveFullDepth) = link%R(thisLink,lr_FullDepth)
                elemSR(:,esr_Orifice_DischargeCoeff)     = link%R(thisLink,lr_DischargeCoeff1)
                elemSR(:,esr_Orifice_Orate)              = link%R(thisLink,lr_DischargeCoeff2)
                elemSR(:,esr_Orifice_Zcrest)             = elemR(:,er_Zbottom) + link%R(thisLink,lr_InletOffset)
                elemSR(:,esr_Orifice_Zcrown)             = elemSR(:,eSr_Orifice_Zcrest) + link%R(thisLink,lr_FullDepth)
                elemSR(:,esr_Orifice_RectangularBreadth) = link%R(thisLink,lr_BreadthScale)
                elemSR(:,esr_Orifice_FullArea)           = elemSR(:,esr_Orifice_RectangularBreadth) * elemSR(:,esr_Orifice_FullDepth)
                elemSR(:,esr_Orifice_EffectiveFullArea)  = elemSR(:,esr_Orifice_RectangularBreadth) * elemSR(:,esr_Orifice_EffectiveFullDepth)    

                !% --- default channel geometry (overwritten later by adjacent CC shape)
                !%     assumes channel is rectangular and twice the breadth of orifice and
                !%     used orifice crown as the maximum overflow
                elemI(:,ei_geometryType)            = rectangular
                elemSGR(:,esgr_Rectangular_Breadth) = twoR * elemSR(:,esr_Orifice_RectangularBreadth)
                elemR(:,er_BreadthMax)              = elemSGR(:,esgr_Rectangular_Breadth)
                elemR(:,er_FullDepth)               = twoR * link%R(thisLink,lr_FullDepth)  
            end where

        case (lCircular)

            where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
                !% integer data
                elemSI(:,esi_Orifice_GeometryType)       = circular
                !% real data
                elemSR(:,esr_Orifice_FullDepth)          = link%R(thisLink,lr_FullDepth)
                elemSR(:,esr_Orifice_EffectiveFullDepth) = link%R(thisLink,lr_FullDepth)
                elemSR(:,esr_Orifice_FullArea)           = (pi / fourR) * elemSR(:,esr_Orifice_FullDepth) ** twoR
                elemSR(:,esr_Orifice_EffectiveFullArea)  = (pi / fourR) * elemSR(:,esr_Orifice_EffectiveFullDepth) ** twoR
                elemSR(:,esr_Orifice_DischargeCoeff)     = link%R(thisLink,lr_DischargeCoeff1)
                elemSR(:,esr_Orifice_Orate)              = link%R(thisLink,lr_DischargeCoeff2)
                elemSR(:,esr_Orifice_Zcrest)             = elemR(:,er_Zbottom) + link%R(thisLink,lr_InletOffset)
                elemSR(:,esr_Orifice_Zcrown)             = elemSR(:,esr_Orifice_Zcrest) + link%R(thisLink,lr_FullDepth)

                !% --- default channel geometry (overwritten later by adjacent CC shape)
                !%     assumes channel is rectangular and twice the breadth of orifice full depth and
                !%     and full depth is greater of orifice crown to bottom or OrificeFullDepth
                elemI(:,ei_geometryType)            = rectangular
                elemSGR(:,esgr_Rectangular_Breadth) = twoR * elemSR(:,esr_Orifice_FullDepth)
                elemR(:,er_BreadthMax)              = elemSGR(:,esgr_Rectangular_Breadth) 
                elemR(:,er_FullDepth)               = twoR * max(elemSR(:,esr_Orifice_Zcrown) - elemR(:,er_Zbottom),elemSR(:,esr_Orifice_FullDepth))
            end where

        case default
            print *, 'In ', subroutine_name
            print *, 'CODE ERROR: unknown orifice geometry type, ', OrificeGeometryType,'  in network'
            print *, 'which has key ',trim(reverseKey(OrificeGeometryType))
            !stop 
            call util_crashpoint(8345553)
            !return
        end select

        !% --- set minimum crest height as 101% of the zero depth value for all orifices
        !%     this ensures that zero-height orifice elements cannot cause flow for zerovalue depths
        thisPack = pack(elemI(:,ei_Lidx),(elemI(:,ei_link_Gidx_BIPquick) == thisLink) ) 
        do ii=1,size(thisPack)
            elemSR(thisPack(ii),esr_Orifice_Zcrest) = &
                max( elemSR(thisPack(ii),esr_Orifice_Zcrest), elemR(thisPack(ii),er_Zbottom) + setting%ZeroValue%Depth*1.01d0 )
        end do
        deallocate(thisPack)

        !% --- initialize a default rectangular channel as the background of the orifice
        call init_IC_diagnostic_default_geometry (thisLink, orifice)


        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_get_orifice_geometry
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_IC_get_pump_geometry (thisLink)
        !%-------------------------------------------------------------------
        !% Description:
        !% get the geometry for pumps
        !%-------------------------------------------------------------------
            integer             :: ii
            integer, intent(in) :: thisLink
            integer, pointer    :: specificPumpType, curveID, lastRow

            character(64) :: subroutine_name = 'init_IC_get_pump_geometry'
        !--------------------------------------------------------------------------
            if (setting%Debug%File%initial_condition) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% --- pointer to specific pump type
        specificPumpType => link%I(thisLink,li_link_sub_type)
        curveID          => link%I(thisLink,li_curve_id)
        lastRow          => curve(curveID)%NumRows

        !% --- cycle through the elements to find this link (HACK -- need to rewrite)
        do ii = 1,N_elem(this_image())
            if (elemI(ii,ei_link_Gidx_BIPquick) == thisLink) then  
                !% real data
                elemSR(ii,esr_Pump_yOn)     = link%R(thisLink,lr_yOn)
                elemSR(ii,esr_Pump_yOff)    = link%R(thisLink,lr_yOff)
                elemR(ii,er_Setting)        = link%R(thisLink,lr_initSetting)
                elemSI(ii,esi_Pump_IsControlled) = zeroI

                print *, '........pump initial condition setting ',elemR(ii,er_Setting)

                !% --- ensure pump ON depth is greater than zero.
                if (elemSR(ii,esr_Pump_yOn) == zeroR) then
                    elemSR(ii,esr_Pump_yOn) = setting%ZeroValue%Depth
                end if

                !% --- ensure pump OFF depth is greater than zero.
                if (elemSR(ii,esr_Pump_yOff) == zeroR) then
                    elemSR(ii,esr_Pump_yOff) = setting%ZeroValue%Depth
                end if

                ! print *, 'pump_yOn     ', elemSR(ii,esr_Pump_yOn)
                ! print *, 'pump_yOff    ', elemSR(ii,esr_Pump_yOff)
                ! print *, 'pump setting ',elemR(ii,er_Setting) 

                ! print *, 'link ',elemI(ii,ei_link_Gidx_BIPquick)
                ! print *, 'pump link name  ',trim(link%Names(elemI(ii,ei_link_Gidx_BIPquick))%str)

                if ((elemSR(ii,esr_Pump_yOff) > elemSR(ii,esr_Pump_yOn)) &
                    .and. (elemSR(ii,esr_Pump_yOff) > zeroR)) then 
                    print *, 'USER CONFIGURATION ERROR:'
                    print *, 'depth/head at which pumps shuts off is larger than'
                    print *, 'the depth/head at which pumps turns on, which provides'
                    print *, 'illogical behavior.'
                    print *, 'Pump at link # ',elemI(ii,ei_link_Gidx_BIPquick)
                    !print *, link%Names(1)%str
                    print *, 'pump link name  ',trim(link%Names(elemI(ii,ei_link_Gidx_BIPquick))%str)
                    call util_crashpoint(6439872)
                end if

                !% --- set nominal element length
                elemR(ii,er_Length)         = setting%Discretization%NominalElemLength

                if (curveID < zeroI) then
                    !% integer data
                    elemSI(ii,esi_Pump_SpecificType) = type_IdealPump 
                else
                    elemSI(ii,esi_Pump_CurveID) = curveID
                    elemSR(ii,esr_Pump_xMin)    = curve(curveID)%ValueArray(1,curve_pump_Xvar)
                    elemSR(ii,esr_Pump_xMax)    = curve(curveID)%ValueArray(lastRow,curve_pump_Xvar)
                    Curve(curveID)%ElemIdx      = ii
                    !% copy pump specific data
                    if (specificPumpType == lType1Pump) then
                        !% integer data
                        elemSI(ii,esi_Pump_SpecificType) = type1_Pump

                    else if (specificPumpType == lType2Pump) then
                        !% integer data
                        elemSI(ii,esi_Pump_SpecificType) = type2_Pump

                    else if (specificPumpType == lType3Pump) then
                        !% integer data
                        elemSI(ii,esi_Pump_SpecificType) = type3_Pump

                    else if (specificPumpType == lType4Pump) then
                        !% integer data
                        elemSI(ii,esi_Pump_SpecificType) = type4_Pump
                    else
                        print *, 'In ', subroutine_name
                        print *, 'CODE ERROR: unknown pump type, ', specificPumpType,'  in network'
                        print *, 'which has key ',trim(reverseKey(specificPumpType))
                        !stop 
                        call util_crashpoint(8863411)
                        !return
                    end if
                end if 
            end if
        end do
     
        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_get_pump_geometry
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_IC_get_outlet_geometry (thisLink)
        !%-----------------------------------------------------------------
        !% get the geometry and other data data for outlet links
        !% Note, these are uncommon -- and are NOT outfalls (which are nodes)
        !%-------------------------------------------------------------------
            integer             :: ii
            integer, intent(in) :: thisLink
            integer, pointer    :: specificOutletType, curveID, eIDx

            character(64) :: subroutine_name = 'init_IC_get_outlet_geometry'
        !%-------------------------------------------------------------------
            if (setting%Debug%File%initial_condition) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% pointer to specific outlet type
        specificOutletType => link%I(thisLink,li_link_sub_type)
        curveID            => link%I(thisLink,li_curve_id)

        do ii = 1,N_elem(this_image())
            if (elemI(ii,ei_link_Gidx_BIPquick) == thisLink) then

                !% real data
                elemSR(ii,esr_Outlet_Coefficient) = link%R(thisLink,lr_DischargeCoeff1)
                elemSR(ii,esr_Outlet_Exponent)    = link%R(thisLink,lr_DischargeCoeff2)
                elemSR(ii,esr_Outlet_Zcrest)      = elemR(ii,er_Zbottom) + link%R(thisLink,lr_InletOffset)

                !% --- set nominal element length
                elemR(ii,er_Length)         = setting%Discretization%NominalElemLength

                if ((specificOutletType == lNodeDepth) .and. (curveID == zeroI)) then
                    !% integer data
                    elemSI(ii,esi_Outlet_SpecificType)  = func_depth_outlet
                elseif ((specificOutletType == lNodeDepth) .and. (curveID /= zeroI)) then
                    !% integer data
                    elemSI(ii,esi_Outlet_SpecificType)  = tabl_depth_outlet
                    elemSI(ii,esi_Outlet_CurveID)       = curveID
                    Curve(curveID)%ElemIdx              = ii
                elseif ((specificOutletType == lNodeHead) .and. (curveID == zeroI)) then
                    !% integer data
                    elemSI(ii,esi_Outlet_SpecificType)  = func_head_outlet
                elseif ((specificOutletType == lNodeHead) .and. (curveID /= zeroI)) then
                    !% integer data
                    elemSI(ii,esi_Outlet_SpecificType)  = tabl_head_outlet
                    elemSI(ii,esi_Outlet_CurveID)       = curveID
                    Curve(curveID)%ElemIdx              = ii
                else
                    print*, 'In ', subroutine_name
                    print*, 'CODE ERROR: unknown outlet type, ', specificOutletType,'  in network'
                    print *, 'which has key ',trim(reverseKey(specificOutletType))
                    call util_crashpoint(82564)
                end if
            end if 
        end do

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_get_outlet_geometry
!%
!%=========================================================================
!%=========================================================================
!%
    subroutine init_IC_diagnostic_default_geometry (thisLink,thisType)
        !%-----------------------------------------------------------------
        !% Description:
        !% Provides default background channel geometry for diagnostic
        !% elements. This geometry will be overwritten by that of an
        !% adjacent cell that has valid channel/conduit geometry
        !% Assume Depth FullDepth, BreadthMax, and ZBottom  and
        !% egsr_Rectangular_Breadth already assigned
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisLink, thisType
            integer             :: geoType
            logical             :: canSurcharge
        !%-----------------------------------------------------------------

        !% --- set the geometry type for the input type
        !%     designed for flexibility, but presently requiring
        !%     rectangular_closed for any diagnostic default geometry
        !%     to prevent issues of overflow
        select case (thisType)
        case (weir)
            geoType = rectangular_closed
            canSurcharge = .true.
        case (orifice)
            geoType = rectangular_closed
            canSurcharge = .true.
        case default
            print *, 'CODE ERROR: unexpected case default'
            call util_crashpoint(5582366)
        end select    

        where (elemI(:,ei_link_Gidx_BIPquick) == thisLink)
            elemI(:,ei_geometryType)            = geoType
            elemR(:,er_Length)                  = setting%Discretization%NominalElemLength
            !elemR(:,er_FullHydDepth)            = elemR(:,er_FullDepth)
            elemR(:,er_FullPerimeter)           = elemR(:,er_BreadthMax) + twoR * elemR(:,er_FullDepth)
            elemR(:,er_ZbreadthMax)             = elemR(:,er_FullDepth) + elemR(:,er_Zbottom)  
            elemR(:,er_Zcrown)                  = elemR(:,er_Zbottom)   + elemR(:,er_FullDepth)
            elemR(:,er_FullArea)                = elemR(:,er_FullDepth) * elemR(:,er_BreadthMax)
            elemR(:,er_FullVolume)              = elemR(:,er_FullArea)  * elemR(:,er_Length)
            elemR(:,er_AreaBelowBreadthMax)     = elemR(:,er_FullArea)
           ! elemR(:,er_FullEllDepth)                 = elemR(:,er_FullDepth)
            
            !% store IC data
            elemR(:,er_Area)          = elemSGR(:,esgr_Rectangular_Breadth) * elemR(:,er_Depth)
            elemR(:,er_Area_N0)       = elemR(:,er_Area)
            elemR(:,er_Area_N1)       = elemR(:,er_Area)
            elemR(:,er_Volume)        = elemR(:,er_Area) * elemR(:,er_Length)
            elemR(:,er_Volume_N0)     = elemR(:,er_Volume)
            elemR(:,er_Volume_N1)     = elemR(:,er_Volume)
            elemR(:,er_EllDepth)      = elemR(:,er_Depth)
            !elemR(:,er_HydDepth)      = elemR(:,er_Depth)
            elemR(:,er_Perimeter)     = twoR * elemR(:,er_Depth) + elemR(:,er_BreadthMax)
            elemR(:,er_HydRadius)     = elemR(:,er_Area) / elemR(:,er_Perimeter)
            elemR(:,er_TopWidth)      = elemR(:,er_BreadthMax)
        endwhere
     
    end subroutine init_IC_diagnostic_default_geometry
!%
!%=========================================================================
!%=========================================================================
!%
    subroutine init_IC_diagnostic_geometry_from_adjacent (isFirstCall)
        !%-----------------------------------------------------------------
        !% Description:   20220511brh
        !% Provides the additional "background" geometry of
        !% diagnostic (pump, outlet only) elements based on its surroundings. This is the
        !% geometry of the channel/conduit in which the diagnostic element exists.  
        !% This ensures that a diagnostic element next to
        !% a JB branch has a valid geometry that can be used for the JB branch.
        !% THIS DOES NOT APPLY TO WEIRS OR ORIFICES, which get their background
        !% geometry from their weir/orifice information.
        !% PUMP -- if upstream element is CC, the pump takes on the
        !%   geometry of the upstream CC element. If the upstream element is
        !%   other than CC, then the pump takes on the geometry of the
        !%   downstream CC element. If the downstream element is also other than 
        !%   CC then an error is returned
        !% Outlet -- requires an upstream CC element
        !%
        !% Called initially for CC adjacent only, then for CC and JB when
        !% after JB have been updated in init_IC_for_nJm_from_nodedata
        !%-----------------------------------------------------------------
        !% Declarations
            logical, intent(in) :: isFirstCall !% true for first time through
            integer, dimension(:), allocatable, target :: packIdx
            integer, pointer :: Fidx, Aidx, thisP
            integer, pointer :: linkIdx
            integer :: ii, jj, Ci
            

            character(64) :: subroutine_name = 'init_IC_diagnostic_geometry_from_adjacent'
        !%-----------------------------------------------------------------
        !% Preliminaries:
            !% --- get the set of pumps, and outlets
            packIdx = pack(elemI(:,ei_Lidx), &
                    ((elemI(:,ei_elementType) .eq. pump) &
                    .or. &
                    (elemI(:,ei_elementType) .eq. outlet) ) )
        !%-----------------------------------------------------------------

        ! print *, ' '
        ! print *, '+++++++++++++++++++++++++++++++++++++++++++++++'
        ! print *, 'in ',trim(subroutine_name)

        !% --- cycle through to set geometry of diagnostic element
        !%     use the upstream geometry if it is CC
        do ii=1,size(packIdx)
            !% --- the present point
            thisP  => packIdx(ii)
        
                ! print *, '=========================='
                ! print *, 'thisP ',thisP
                ! print *, 'this elem type ',trim(reverseKey(elemI(thisP,ei_elementType)))
                ! print *, 'geometry type #',elemI(thisP,ei_geometryType)
                ! if (elemI(thisP,ei_geometryType) .ne. nullvalueI) then
                !     print *, 'geometry type name ',trim(reverseKey(elemI(thisP,ei_geometryType)))
                ! end if

            !% --- cycle if not a nullvalue geometry type
            if (elemI(thisP,ei_geometryType) .ne. undefinedKey) cycle 

            !% --- the link
            linkIdx => elemI(thisP,ei_link_Gidx_SWMM)
                ! print *, 'linkIdx  ' ,linkIdx
                ! print *, 'linkname ' ,trim(link%Names(linkIdx)%str)

            !% --- UPSTREAM ELEMENTS ----------------------------------------
            !% --- the upstream face
            Fidx => elemI(thisP,ei_Mface_uL)
                ! print *, 'up face ',Fidx

            !% --- identify the upstream element
            !%     which may be on a different image
            if (elemYN(thisP,eYN_isBoundary_up)) then
                Ci   =  faceI(Fidx,fi_Connected_image)
                Aidx => faceI(Fidx,fi_GhostElem_uL)
            else
                Ci   =  this_image()
                Aidx => faceI(Fidx,fi_Melem_uL)
            end if

                ! print *, 'Ci ',Ci
                ! print *, 'Up element  ',Aidx
                ! print *, 'Up elememtType ',trim(reverseKey(elemI(Aidx,ei_elementType)[Ci]))

                ! print *, 'is first call ',isFirstCall

            !% --- set geometry for thisP based on upstream elements where possible
            if (isFirstCall) then
                if (elemI(Aidx,ei_elementType)[Ci] == CC) then
                    call init_IC_set_implied_geometry (thisP, Aidx, Ci)
                else
                    !% --- if the upstream element is not CC, use the downstream element CC geometry
                    !%     for pumps, but fail for outlets
                    if (elemI(thisP,ei_elementType) == outlet) then
                        !% --- outlets are required to have upstream CC
                        print *, 'USER SYSTEM CONFIGURATION ERROR'
                        print *, 'An outlet requires at least one upstream link that is a'
                        print *, 'conduit or channel. This condition violated for'
                        print *, 'outlet with name ',trim(link%Names(linkIdx)%str)
                        call util_crashpoint(92873)
                    end if
                end if
            else 
                if ((elemI(Aidx,ei_elementType)[Ci] == CC) .or.        &
                    (elemI(Aidx,ei_elementType)[Ci] == JB)      ) then
                    call init_IC_set_implied_geometry (thisP, Aidx, Ci) 
                else
                    print *, 'NOT CC OR JB'
                    stop 590873 
                end if
            end if

            ! print *, 'new type ',trim(reverseKey(elemI(thisP,ei_geometryType)))
            ! print *, ' '

            !% --- Look downstream if this element still undefined
            if (elemI(thisP,ei_geometryType) .ne. undefinedKey) cycle 

            !% --- the downstream face
            Fidx => elemI(thisP,ei_Mface_dL)
                ! print *, 'dn face ',Fidx

            !% --- the downstream element
            !%     which may be on a different image
            if (elemYN(thisP,eYN_isBoundary_dn)) then
                Ci   =  faceI(Fidx,fi_Connected_image)
                Aidx => faceI(Fidx,fi_GhostElem_dL)
            else
                Ci   =  this_image()
                Aidx => faceI(Fidx,fi_Melem_dL)
            end if

                ! print *, 'Ci ',Ci
                ! print *, 'Dn element  ',Aidx
                ! print *, 'Dn elememtType ',trim(reverseKey(elemI(Aidx,ei_elementType)[Ci]))

                ! print *, 'downstream link     ',elemI(Aidx,ei_link_Gidx_BIPquick)
                ! if (elemI(Aidx,ei_link_Gidx_BIPquick) .ne. nullvalueI) then
                !     print *, 'downstream linkname ' ,trim(link%Names(elemI(Aidx,ei_link_Gidx_BIPquick))%str)
                ! endif
                ! print *, 'downstream node   ',elemI(Aidx,ei_node_Gidx_BIPquick)
                ! if (elemI(Aidx,ei_node_Gidx_BIPquick) .ne. nullvalueI) then
                !     print *, 'downstream nodename ',trim(node%Names(elemI(Aidx,ei_node_Gidx_BIPquick))%str)
                ! end if

            if (isFirstCall) then    
                !% --- the element type downstream
                if (elemI(Aidx,ei_elementType)[Ci] == CC) then
                    call init_IC_set_implied_geometry (thisP, Aidx, Ci) 
                else
                    ! if (elemI(Aidx,ei_elementType)[Ci] == JB) then
                    !     !% --- pump with both upstream and downstream not CC
                    !     !%     downstream is JB and upstream may be JB
                    !     !%     must wait to resolve geometry after JB assigned
                    !     !%     Assign nullvalueI to find this pump later.
                    !     elemI(thisP,ei_geometryType) = nullvalueI
                    ! else  
                    !     print *, ' '

                    !     !% --- pumps do not have default channel geometry, so they must
                    !     !%     have a CC element upstream or downstream.
                    !     print *, 'USER SYSTEM CONFIGURATION ERROR'
                    !     print *, 'A pump requires at least one upstream or downstream link that is a'
                    !     print *, 'conduit or channel or junction. This condition violated for'
                    !     print *, 'pump with name ',trim(link%Names(linkIdx)%str)
                    !     call util_crashpoint(2398789)
                    ! end if
                end if
            else 
                if ((elemI(Aidx,ei_elementType)[Ci] == CC) .or.        &
                    (elemI(Aidx,ei_elementType)[Ci] == JB)      ) then
                    call init_IC_set_implied_geometry (thisP, Aidx, Ci) 
                end if
            end if
        end do

        !%-----------------------------------------------------------------
        !% Closing:
            deallocate(packIdx)

    end subroutine init_IC_diagnostic_geometry_from_adjacent
!%
!%==========================================================================
!%==========================================================================
!% 
    subroutine init_IC_diagnostic_JB_bounded ()
        !%-----------------------------------------------------------------
        !% Description: 
        !% identifies the special case diagnostic elements that have JB
        !% on either side.
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, dimension(:), allocatable, target :: packIdx
            integer, pointer :: eIdx, fUp, fDn, AidxUp, AidxDn
            integer          :: ii, CiUp, CiDn
        !%-----------------------------------------------------------------
        !% Prelminiaries
            !% --- initialize all faceYN(:,fYN_isJB_QfrozenByDiag) to false. 
            !%     reset to true only if diagnostic bounded by two junctions
            !%     is found
            faceYN(:,fYN_isJB_QfrozenByDiag) = .false.

            !% --- initialize elemSI(:,esi_JunctionBranchCanModifyQ) to oneI
            !%     for all JB
            packIdx = pack(elemI(:,ei_Lidx), (elemI(:,ei_elementType) .eq. JB))
            if (size(packIdx) < 1) return !% no JB found, so not possible
            !% --- initialization to allowing modification
            elemSI(packIdx(1:size(packIdx)),esi_JunctionBranch_CanModifyQ) = oneI 

            deallocate(packIdx)
            
            !% --- get the set of weirs, orifices, and pumps (does not include outlet)
            packIdx = pack(elemI(:,ei_Lidx),              &
                ((elemI(:,ei_elementType) .eq. pump)      &
                 .or.                                     &
                 (elemI(:,ei_elementType) .eq. weir)      &
                 .or.                                     &
                 (elemI(:,ei_elementType) .eq. orifice)) )
        !%-----------------------------------------------------------------

        !% --- cycle through diagnostic elements        
        do ii=1,size(packIdx)
            !% --- element and face indexes on this image
            eIdx  => packIdx(ii)
            fUp   => elemI(eIdx,ei_Mface_uL)
            fDn   => elemI(eIdx,ei_Mface_dL)

            !% --- identify upstream element
            !%     which may be on a different image
            if (elemYN(eIdx,eYN_isBoundary_up)) then 
                CiUp   =  faceI(fUp,fi_Connected_image)
                AidxUp => faceI(fUp,fi_GhostElem_uL)
            else
                CiUp   =  this_image()
                AidxUp => faceI(fUp,fi_Melem_uL)
            end if

            !% --- identify downstream element
            !%     which may be on a different image
            if (elemYN(eIdx,eYN_isBoundary_dn)) then 
                CiDn   =  faceI(fDn,fi_Connected_image)
                AidxDn => faceI(fDn,fi_GhostElem_dL)
            else
                CiDn   =  this_image()
                AidxDn => faceI(fDn,fi_Melem_dL)
            end if

            if ((elemI(AidxUp,ei_elementType)[CiUp] == JB)   &
                .and.                                        &
                (elemI(AidxDn,ei_elementType)[CiDn] == JB) ) then
                !% --- diagnostic that requires special treatment
                !%     during junction computation
                faceYN(fUp,fYN_isJB_QfrozenByDiag) = .true.  
                faceYN(fDn,fYN_isJB_QfrozenByDiag) = .true.
                elemSI(eIdx,esi_JunctionBranch_CanModifyQ) = zeroI
            end if

        end do

        deallocate(packIdx)

    end subroutine init_IC_diagnostic_JB_bounded
!%
!%==========================================================================
!%==========================================================================
!% 
    subroutine init_IC_set_implied_geometry (thisP, Aidx, Ci)    
        !%-----------------------------------------------------------------
        !% Description
        !% Copies geometry from adjacent element Aidx in connected image Ci
        !% to thisP element
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisP, Aidx, Ci
            integer             :: thisCol(5), ii
        !%-----------------------------------------------------------------
        !% Preliminaries:
        !% --- set the column indexes of the fixed geometry data that are
        !%     independent of Z bottom and length that
        !%     are needed in a diagnostic element
            ii=1
            thisCol(ii) = er_AreaBelowBreadthMax; ii=ii+1
            thisCol(ii) = er_BreadthMax;          ii=ii+1
            !thisCol(ii) = er_FullEllDepth;             ii=ii+1
            thisCol(ii) = er_FullArea;            ii=ii+1
            thisCol(ii) = er_FullDepth;           ii=ii+1
            !thisCol(ii) = er_FullHydDepth;        ii=ii+1
            thisCol(ii) = er_FullPerimeter;       ii=ii+1
        !%-----------------------------------------------------------------
        
        !% --- if an adjacent element is a channel/conduit, use this for the background channel
        !$     geometry of the diagnostic element in which the weir/orifice/pump/outlet is embeded
        elemI(thisP,ei_geometryType) = elemI(Aidx,ei_geometryType)[Ci]
        elemR(thisP,thisCol)         = elemR(Aidx,thisCol)[Ci]

        !% --- initialize other consistent terms based on local length and zbottom
        elemR(thisP,er_FullVolume)   = elemR(thisP,er_FullArea) * elemR(thisP,er_Length)
        elemR(thisP,er_ZbreadthMax)  = elemR(thisP,er_Zbottom) &
                                        + elemR(Aidx,er_ZbreadthMax) - elemR(Aidx,er_Zbottom)
        elemR(thisP,er_Zcrown)       = elemR(thisP,er_Zbottom) &
                                             + elemR(Aidx,er_Zcrown) - elemR(Aidx,er_Zbottom)
        !% --- copy special geometry
        call init_IC_diagnostic_special_geometry (thisP, Aidx, Ci)

    end subroutine init_IC_set_implied_geometry
!%
!%==========================================================================
!%==========================================================================
!% 
    subroutine init_IC_diagnostic_special_geometry (thisP, Aidx, Ci)
        !%-----------------------------------------------------------------
        !% Description:
        !% Copies the special fixed geoemtry (depends on element geometry type)
        !% from the adjacent cell (Aidx) to this cell (thisP) where
        !% Aidx is on the connected image (Ci). This is used to get the
        !% geometry for a JB junction branch
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: thisP, Aidx, Ci
            character(64) :: subroutine_name = 'init_IC_diagnostic_special_geometry'
        !%-----------------------------------------------------------------
        !%-----------------------------------------------------------------
        !% --- copy over special geometry data depending on geometry type
        select case (elemI(thisP,ei_geometryType))
            case (rectangular)
                elemSGR(thisP,esgr_Rectangular_Breadth)    = elemSGR(Aidx,esgr_Rectangular_Breadth)[Ci]
            case (trapezoidal)
                elemSGR(thisP,esgr_Trapezoidal_Breadth)    = elemSGR(Aidx,esgr_Trapezoidal_Breadth)[Ci]
                elemSGR(thisP,esgr_Trapezoidal_LeftSlope)  = elemSGR(Aidx,esgr_Trapezoidal_LeftSlope)[Ci]
                elemSGR(thisP,esgr_Trapezoidal_RightSlope) = elemSGR(Aidx,esgr_Trapezoidal_RightSlope)[Ci]
            case (triangular)
                elemSGR(thisP,esgr_Triangular_TopBreadth)  = elemSGR(Aidx,esgr_Triangular_TopBreadth)[Ci]
                elemSGR(thisP,esgr_Triangular_Slope)       = elemSGR(Aidx,esgr_Triangular_Slope)[Ci] 
            case (irregular)
                elemI(thisP,ei_link_transect_idx)          = elemI(Aidx,ei_link_transect_idx)[Ci]
            case (rectangular_closed)
                elemSGR(thisP,esgr_Rectangular_Breadth)    = elemSGR(Aidx,esgr_Rectangular_Breadth)[Ci]
            case (circular)
                elemSGR(thisP,esgr_Circular_Diameter)      = elemSGR(Aidx,esgr_Circular_Diameter)[Ci]
                elemSGR(thisP,esgr_Circular_Radius)        = elemSGR(Aidx,esgr_Circular_Radius)[Ci]
            case default
                print *, 'CODE ERROR unexpected geometry'
                print *, 'ei_geometryType index # ',elemI(thisP,ei_geometryType)
                print *, 'which represents ',reverseKey(elemI(thisP,ei_geometryType))
                print *, 'is not handled in subroutine ',trim(subroutine_name)
                call util_crashpoint(99376)
        end select
                
    end subroutine init_IC_diagnostic_special_geometry
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_IC_for_nJm_from_nodedata ()
        !--------------------------------------------------------------------------
        !% get the initial depth, and geometry data from nJm nodes
        !--------------------------------------------------------------------------

            integer                       :: ii, image, pJunction
            integer, pointer              :: thisJunctionNode
            integer, allocatable, target  :: packed_nJm_idx(:)

            character(64) :: subroutine_name = 'init_IC_for_nJm_from_nodedata'
        !--------------------------------------------------------------------------
            !if (crashYN) return
            if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% Setting the local image value
        image = this_image()

        !% pack all the link indexes in an image
        packed_nJm_idx = pack( node%I(:,ni_idx), &
                             ((node%I(:,ni_P_image) == image) .and. &
                              (node%I(:,ni_node_type) == nJm) ) )

        ! print *, 'packed ',packed_nJm_idx   
        ! print *, node%I(:,ni_node_type)
   

        !% find the number of links in an image
        pJunction = size(packed_nJm_idx)

        ! print *, 'pJunction ',pJunction

        !% cycle through the links in an image
        do ii = 1,pJunction
                ! print *, ' '
                ! print *, 'pJunction  in nJM from nodedata ',ii, packed_nJm_idx(ii)

            !% necessary pointers
            thisJunctionNode => packed_nJm_idx(ii)
            call init_IC_get_junction_data (thisJunctionNode)
        end do

        !% deallocate the temporary array
        deallocate(packed_nJm_idx)

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_for_nJm_from_nodedata
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_IC_test_nJ2_data ()

        integer :: ii

        do ii=1,N_node
            print *, ii
            print *, node%I(ii,ni_node_type), reverseKey(node%I(ii,ni_node_type))
            print *, node%I(ii,ni_N_link_u), node%I(ii,ni_N_link_d)
            print *, 'curve ID      ',node%I(ii,ni_curve_ID)
            print *, 'assigned      ',node%I(ii,ni_assigned)
            print *, 'elem idx      ',node%I(ii,ni_elem_idx)
            print *, 'face idx      ',node%I(ii,ni_face_idx)
            print *, 'Z bottom      ',node%R(ii,nr_Zbottom)
            print *, 'init depth    ',node%R(ii,nr_InitialDepth)
            print *, 'full depth    ',node%R(ii,nr_FullDepth)
        end do

        print *, 'up element ', faceI(7,fi_Melem_uL)
        print *, 'up element ', faceI(13,fi_Melem_uL)

        stop 509873

    end subroutine init_IC_test_nJ2_data
!%
!%==========================================================================
!%==========================================================================
!
    subroutine init_IC_get_junction_data (thisJunctionNode)        
        !%-----------------------------------------------------------------
        !% Description:
        !% get data for the multi branch junction elements
        !%-----------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: thisJunctionNode

            integer              :: ii, jj, kk, JMidx, JBidx, Aidx, Ci
            integer, pointer     :: BranchIdx, JBgeometryType, JmType, curveID, NumRows
            integer, pointer     :: Fidx, F2idx
            integer              :: nbranches
            real(8), allocatable :: integrated_volume(:)
            real(8)              :: LupMax, LdnMax
            real(8)              :: aa,bb, lowZ
            logical              :: isupstream

            character(64) :: subroutine_name = 'init_IC_get_junction_data'
        !%--------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%initial_condition) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%................................................................
        !% Junction main
        !%................................................................
        !% find the first element ID associated with that nJm
        !% masked on the global node number for this node.

        ! print *, 'thisJunctionNode ',thisJunctionNode
        ! print *, ' '

        ! print *, 'Lidx'
        ! print *, elemI(:,ei_Lidx)

        ! print *, ' '
        ! print *, elemI(:,ei_node_Gidx_SWMM)

        JMidx = minval(elemI(:,ei_Lidx), elemI(:,ei_node_Gidx_SWMM) == thisJunctionNode)

        ! print *, ' '
        ! do ii=1,N_elem(1)
        !     print *, ii, elemI(ii,ei_Lidx), elemI(ii,ei_node_Gidx_SWMM)
        ! end do

        ! print *, '    ============================================='
        ! print *, 'JMidx ',JMidx, ' ',node%YN(thisJunctionNode,nYN_has_storage)

        !% the first element index is a junction main
        elemI(JMidx,ei_elementType)  = JM
        elemI(JMidx,ei_HeqType)      = time_march
        elemI(JMidx,ei_QeqType)      = notused

        !% set the type of junction main
        if (node%YN(thisJunctionNode,nYN_has_storage)) then
            if (node%I(thisJunctionNode,ni_curve_ID) .eq. 0) then
                !% --- functional storage
                elemSI(JMidx,esi_JunctionMain_Type)   = FunctionalStorage
                elemSR(JMidx,esr_Storage_Constant)    = node%R(thisJunctionNode,nr_StorageConstant)
                elemSR(JMidx,esr_Storage_Coefficient) = node%R(thisJunctionNode,nr_StorageCoeff)
                elemSR(JMidx,esr_Storage_Exponent)    = node%R(thisJunctionNode,nr_StorageExponent)                    
            else
                !% --- tabular storage
                elemSI(JMidx,esi_JunctionMain_Type) = TabularStorage
                elemSI(JMidx,esi_JunctionMain_Curve_ID) = node%I(thisJunctionNode,ni_curve_ID)
            end if
            !% --- common data
            elemSR(JMidx,esr_Storage_FractionEvap)= node%R(thisJunctionNode,nr_StorageFevap)
        else
            !%-----------------------------------------------------------------------
            !% Junction main with implied or no storage
            !%-----------------------------------------------------------------------
            if (setting%Junction%ForceStorage) then 
                !% --- implied storage
                elemSI(JMidx,esi_JunctionMain_Type)    = ImpliedStorage
                setting%Junction%SurfaceArea_Minimum   = setting%SWMMinput%SurfaceArea_Minimum
            else 
                !% --- no storage
                elemSI(JMidx,esi_JunctionMain_Type)    = NoStorage
                setting%Junction%SurfaceArea_Minimum   = zeroR
            end if
            elemI (JMidx,ei_geometryType)          = rectangular
            elemSR(JMidx,esr_Storage_FractionEvap) = zeroR  !% --- no evap from implied storage junction

            ! !% --- check for low Zbottom in an implied storage
            ! !%     Require junction bottom to be at or above the minimum branch
            ! !%     Zbottom with implied storage
            ! lowZ = huge(oneR)
            ! do kk=1,max_branch_per_node
            !     if (elemSI(JMidx+kk,esi_JunctionBranch_Exists) == oneI) then 
            !         lowZ = min(lowZ, elemR(JMidx+kk,er_Zbottom))
            !     end if
            ! end do
            ! elemR(JMidx,er_Zbottom) = max(elemR(JMidx,er_Zbottom), lowZ)
        end if

        ! print *, 'type ',elemSI(JMidx,esi_JunctionMain_Type),ImpliedStorage

        !% --- junction main depth and head from initial conditions
        elemR(JMidx,er_Depth)     = node%R(thisJunctionNode,nr_InitialDepth)
        elemR(JMidx,er_Head)      = elemR(JMidx,er_Depth) + elemR(JMidx,er_Zbottom)
        elemR(JMidx,er_FullDepth) = node%R(thisJunctionNode,nr_FullDepth)
        elemR(JMidx,er_Zcrown)    = elemR(JMidx,er_FullDepth) + elemR(JMidx,er_Zbottom)
        
        !% --- overflow volume accumulator
        elemR(JMidx,er_VolumeOverFlowTotal) = zeroR

        elemR(JMidx,er_VolumeArtificialInflowTotal) = zeroR

        !% --- ponded area is stored in elemSR array
        if (setting%SWMMinput%AllowPonding) then
            elemSR(JMidx,esr_JunctionMain_PondedArea) = node%R(thisJunctionNode,nr_PondedArea)
        else
            elemSR(JMidx,esr_JunctionMain_PondedArea) = zeroR
        end if
        elemR(JMidx,er_VolumePonded) = zeroR

        !% --- all JM "can" surcharge, but are limited by their SurchargeExtraDepth
        !%     which might be 0 (effectively preventing any surcharge)
        elemYN(JMidx,eYN_canSurcharge) = .true.

        !% --- check for initialization of surcharge extra depth
        if (node%R(thisJunctionNode,nr_SurchargeExtraDepth) == nullvalueR) then 
            print *, 'ERROR: Surcharge Extra Depth at a junction not initialized'
            print *, 'This should not happen! Likely problem forinitialization code'
            call util_crashpoint(8838723)
        end if

        !% --- Set the extra head above the crown for maximum surcharge at Junction
        elemSR(JMidx,esr_JunctionMain_SurchargeExtraDepth)      &
            = node%R(thisJunctionNode,nr_SurchargeExtraDepth)

        !% --- Set the overflow and surcharge conditions
        !% --- check for infinite extra depth 
        !%     if infinite marker (1000) is used, then no oveflow allowed
        !%     applies to both 1000 m and 1000 ft as input.
        if  ( ( (elemSR(JMidx,esr_JunctionMain_SurchargeExtraDepth)                &
                .le. 1.001d0 * setting%Junction%InfiniteExtraDepthValue)           &
                .and.                                                              &
                (elemSR(JMidx,esr_JunctionMain_SurchargeExtraDepth)                &
                .ge. 0.999d0 * setting%Junction%InfiniteExtraDepthValue)           &
                )                                                                  &
            .or.                                                                   &
                ( (elemSR(JMidx,esr_JunctionMain_SurchargeExtraDepth)              &
                .le. 1.001d0 * setting%Junction%InfiniteExtraDepthValue*0.3048d0)  & 
                .and.                                                              &
                (elemSR(JMidx,esr_JunctionMain_SurchargeExtraDepth)                &
                .ge. 0.999d0 * setting%Junction%InfiniteExtraDepthValue*0.3048d0)  & 
                )                                                                  &
            ) then 
            !% --- set type to NoOverflow 
            elemSI(JMidx,esi_JunctionMain_OverflowType) = NoOverflow    
        else
            !% --- not infinite depth
            if (elemSR(JMidx,esr_JunctionMain_SurchargeExtraDepth) .eq. zeroR) then 
                !% --- open channel storage
                if (elemSR(JMidx,esr_JunctionMain_PondedArea) == zeroR) then
                    !% --- use the overflow weir algorithm
                    elemSI(JMidx,esi_JunctionMain_OverflowType) = OverflowWeir
                else
                    !% --- use ponded overflow algorithm
                    elemSI(JMidx,esi_JunctionMain_OverflowType) = Ponded 
                end if
            else 
                !% --- closed conduit overflow
                if (elemSR(JMidx,esr_JunctionMain_PondedArea) == zeroR) then
                    !% --- use oveflow orifice
                    elemSI(JMidx,esi_JunctionMain_OverflowType) = OverflowOrifice
                    !% --- HACK using default orifice length and height for overflow
                    elemSR(JMidx,esr_JunctionMain_OverflowOrifice_Length) = setting%Junction%OverflowOrificeLength
                    elemSR(JMidx,esr_Junctionmain_OverflowOrifice_Height) = setting%Junction%OverflowOrificeHeight
                else
                    !% --- use ponded overflow
                    elemSI(JMidx,esi_JunctionMain_OverflowType) = Ponded 
                end if
            end if
        end if

        !else 
            !elemSR(JMidx,esr_JunctionMain_SurchargeExtraDepth) = zeroR

            !% --- THE StorageSurchargeExtraDepth IS NOT YET IMPLEMENTED!
            !%     It requires computation of elemR(:,er_FullArea) in initial_conditions
            !%     for tabular and functional storage
            !% --- for storage nodes without extra surcharge depth, check
            !%     for a default surcharge depth.
            !%     NOTE: this code is needed because EPA-SWMM does not allow
            !%     storage nodes to surcharge, hence it does not provide an
            !%     input of Extra Surcharge Depth. SWMM5+ allows 
            !%     a storage node to surcharge, but there needs to be a 
            !%     future extension to provide the value for each storage node.  
            !%     As an intermediate step we use a default "StorageSurchargeExtraDepth" 
            !%     a the max surcharge at a Storage node.
            ! if ((node%YN(thisJunctionNode,nYN_has_storage)) .and. &
            !     (setting%Junction%StorageSurchargeExtraDepth > zeroR) ) then
            !     !% --- use the default extra surcharge depth    
            !     elemSR(JMidx,esr_JunctionMain_SurchargeExtraDepth) &
            !         = setting%Junction%StorageSurchargeExtraDepth 
            ! else
            !     !% --- if NOT a storage node or if IS a storage node
            !     !%     but the setting%Junction%StorageSurchargeExtraDepth = 0.0
            !     !%     then the junction cannot surcharge.
            !     elemSR(JMidx,esr_JunctionMain_SurchargeExtraDepth) = zeroR
            ! end if
        !end if    

        !% JM elements are not solved for momentum.
        elemR(JMidx,er_Flowrate)     = zeroR
        elemR(JMidx,er_Velocity)     = zeroR

        !% JM elements always have a single barrel
        elemI(JMidx,ei_barrels)      = oneR

        !% wave speed is the gravity wave speed for the depth
        elemR(JMidx,er_WaveSpeed)    = sqrt(setting%constant%gravity * elemR(JMidx,er_Depth))
        elemR(JMidx,er_FroudeNumber) = zeroR

        !% --- self index
        elemSI(JMidx,esi_JunctionBranch_Main_Index ) = JMidx

        ! print *, 'JM values ',elemR(JMidx,er_Head), elemR(JMidx,er_Head) - 198.12d0
        ! print *, 'Depth     ', elemR(JMidx,er_Depth)
        ! print *, 'zbottom   ', elemR(JMidx,er_Zbottom)
   
        !%................................................................
        !% Junction Branches
        !%................................................................
        !% loopthrough all the branches
        !% HACK -- replace much of this with a call to the standard geometry after all other IC have
        !% been done. The branch depth should be based on the upstream or downstream depth of the
        !% adjacent element.
        do ii = 1,max_branch_per_node

            ! print *, '================================',trim(subroutine_name)
            ! print *, ii, 'in branch for JMidx ',JMidx
            ! print *, 'node ',elemI(JMidx,ei_node_Gidx_SWMM)
            ! print *, 'node name ',trim(node%Names(elemI(JMidx,ei_node_Gidx_SWMM))%str)

            !% 20220406brh Rewritten to use adjacent element geometry initialization where possible.

            !% --- find the element id of junction branches
            JBidx = JMidx + ii

            ! print *, 'JBidx ',JBidx
            
            !% --- main index associated with branch
            !elemI(JBidx,ei_main_idx_for_branch) = JMidx
            elemSI(JBidx,esi_JunctionBranch_Main_Index) = JMidx

            !% --- set the adjacent element id (Aidx) where elemI and elemR data can be extracted
            !%     note that elemSGR etc are not yet assigned
            if (mod(ii,2) == zeroI) then
                Fidx => elemI(JBidx,ei_MFace_dL)
                !% --- even are downstream branches
                isupstream = .false.
                elemSI(JBidx,esi_JunctionBranch_IsUpstream) = zeroI
                if (elemYN(JBidx,eYN_isBoundary_dn)) then
                    Ci   = faceI(Fidx,fi_Connected_image)
                    Aidx = faceI(Fidx,fi_GhostElem_dL)
                else
                    Ci   = this_image()
                    Aidx = faceI(Fidx,fi_Melem_dL)
                end if
            else
                !% --- odd are upstream branches
                isupstream = .true.
                elemSI(JBidx,esi_JunctionBranch_IsUpstream) = oneI
                Fidx => elemI(JBidx,ei_MFace_uL)
                if (elemYN(JBidx,eYN_isBoundary_up)) then
                    Ci   = faceI(Fidx,fi_Connected_image)
                    Aidx = faceI(Fidx,fi_GhostElem_uL)
                else
                    Ci   = this_image()
                    Aidx = faceI(Fidx,fi_Melem_uL)
                end if
            end if

            !% --- set the junction branch element type
            elemI(JBidx,ei_elementType) = JB

            ! print *, 'exist ', elemSI(JBidx,esi_JunctionBranch_Exists)

            !% --- set the geometry for existing branches
            !%     Note that elemSI(,...Exists) is set in init_network_handle_nJm
            if (.not. elemSI(JBidx,esi_JunctionBranch_Exists) == oneI) cycle

            BranchIdx      => elemSI(JBidx,esi_JunctionBranch_Link_Connection)
            JBgeometryType => link%I(BranchIdx,li_geometry)

            ! print *, 'linkgeo', link%I(BranchIdx,li_geometry), reverseKey(link%I(BranchIdx,li_geometry))

            !% --- set the JB to time_march for use with splitting between AC
            !%     and ETM in rk2_extrapolate_to_fullstep_ETM, rk2_restore_to_midstep_ETM
            !%     rk2_interpolate_to_halfstep_AC, k2_restore_to_fullstep_AC

            !% TESTING REMOVAL 20220720 brh
            elemI(JBidx,ei_HeqType) = notused !% time_march
            elemI(JBidx,ei_QeqType) = notused !%time_march

            !% ---Junction branch k-factor 
            !%    If the user does not input the K-factor for junction branches entrance/exit loses then
            !%    use default from setting
            if (node%R(thisJunctionNode,nr_JunctionBranch_Kfactor) .ne. nullvalueR) then
                elemSR(JBidx,esr_JunctionBranch_Kfactor) = node%R(thisJunctionNode,nr_JunctionBranch_Kfactor)
            else
                elemSR(JBidx,esr_JunctionBranch_Kfactor) = setting%Junction%kFactor
            end if

            !% --- set the initial head and to the same as the junction main
            elemR(JBidx,er_Head)    = elemR(JMidx,er_Head)
            !% --- set the depth consistent with JB bottom
            elemR(JBidx,er_Depth)   = elemR(JBidx,er_Head) - elemR(JBidx,er_Zbottom)
            !% --- check for dry
            if (elemR(JBidx,er_Head) < elemR(JBidx,er_Zbottom)) then
                elemR(JBidx,er_Head) = elemR(JBidx,er_Zbottom)
                elemR(JBidx,er_Depth) = zeroR
            end if

            elemR(JBidx,er_VolumeOverFlow) = zeroR
            elemR(JBidx,er_VolumeOverFlowTotal) = zeroR

            elemR(JBidx,er_VolumeArtificialInflowTotal) = zeroR

            !% --- JB elements initialized for momentum
            elemR(JBidx,er_Flowrate)     = elemR(Aidx,er_Flowrate)[Ci] !% flowrate of adjacent element
            elemR(JBidx,er_WaveSpeed)    = sqrt(setting%constant%gravity * elemR(JBidx,er_Depth))
            elemR(JBidx,er_FroudeNumber) = zeroR

            !% --- Set the geometry from the adjacent elements
            !%     Must evaluate across connected images
            !%     JB inherits geometry type from connected branch
            elemI(JBidx,ei_geometryType)        = elemI(Aidx,ei_geometryType)[Ci]

            !% --- handle nullvalue geometry (can occur when adjacent element is diagnostic)
            !%     Looks for the next link upstream. If it is a channel or
            !%     conduit then its geometry can be assigned to the JB.
            if (elemI(Aidx,ei_geometryType)[Ci] == undefinedKey) then 
                call init_IC_JB_nullvalue_geometry &
                    (Aidx, Ci, thisJunctionNode, JBidx, isupstream)
            end if

            !% --- branch has same number of barrels as the connected element
            elemI(JBidx,ei_barrels)             = elemI(Aidx,ei_barrels)[Ci]

            !% --- Set the face flowrates and barrels such that it does not blowup

            if (elemI(JBidx, ei_Mface_uL) /= nullvalueI) then
                !print *, elemI(JBidx, ei_Mface_uL), faceI(elemI(JBidx, ei_Mface_uL),fi_barrels)
                faceR(elemI(JBidx, ei_Mface_uL),fr_flowrate) = elemR(JBidx,er_Flowrate) 
                faceI(elemI(JBidx, ei_Mface_uL),fi_barrels)  = elemI(JBidx,ei_barrels) 
            else if (elemI(JBidx, ei_Mface_dL) /= nullvalueI) then
                !print *, elemI(JBidx, ei_Mface_dL), faceI(elemI(JBidx, ei_Mface_dL),fi_barrels)
                faceR(elemI(JBidx, ei_Mface_dL),fr_flowrate) = elemR(JBidx,er_Flowrate)
                faceI(elemI(JBidx, ei_Mface_dL),fi_barrels)  = elemI(JBidx,ei_barrels)  
            else 
                print *, 'JBidx null face both down and up ',JBidx
                stop 650987
            end if

            !% --- check if the connected element is CC
            if (elemI(Aidx,ei_elementType) == CC) then 
                elemSI(JBidx,esi_JunctionBranch_CC_adjacent) = oneI
            else
                elemSI(JBidx,esi_JunctionBranch_CC_adjacent) = zeroI
            end if

            !% --- check if the connected element is Diagnostic weir, pump or orifice
            if ((elemI(Aidx,ei_elementType) == weir)    .or. &
                (elemI(Aidx,ei_elementTYpe) == orifice) .or. &
                (elemI(Aidx,ei_elementType) == pump)           ) then 
                elemSI(JBidx,esi_JunctionBranch_Diag_adjacent) = oneI
            else
                elemSI(JBidx,esi_JunctionBranch_Diag_adjacent) = zeroI
            end if

            !% --- Ability to surcharge is set by JM
            !%     Note that JB (if surcharged) isn't subject to the max surcharge depth 
            !%     of its JM. That is, a JB, if allowed to surcharge can surcharge to any
            !%     level, but typically won't be much about the JM since the JM head
            !%     drives the JB head.
            !%     Note that this might be perceived as a logic problem: a branch 
            !%     inherits geometry of the adjacent element,
            !%     which allows "surcharge" to exist on a branch that is considered
            !%     an open channel. This occurs when a channel is draining into
            !%     a closed junction. In this case we think of the JB as
            !%     having the flow characteristics of the adjacent channel, but
            !%     the head is inherited from the JM. Thus, a JB can have open
            !%     channel flow characteristics but a head based on the associated
            !%     closed JM.
            if (elemSR(JMidx,esr_JunctionMain_SurchargeExtraDepth) > zeroR) then 
                !% --- where JM is allowed to surcharge
                elemYN(JBidx,eYN_canSurcharge) = .true.
            else 
                !% --- where JM surcharge is limited to zero
                elemYN(JBidx,eYN_canSurcharge) = .false.
            end if

            !% --- OLD APPROACH (20220922) using surcharge based on connected element
            !elemYN(JBidx,eYN_canSurcharge)      = elemYN(Aidx,eYN_canSurcharge)[Ci]


            ! print *, 'JBidx',JBidx, elemI(JBidx,ei_geometryType), trim(reverseKey(elemI(JBidx,ei_geometryType)))
            ! print *, 'Aidx ',Aidx,  elemI(Aidx,ei_geometryType),  trim(reverseKey(elemI(Aidx,ei_geometryType)))

            select case  (elemI(JBidx,ei_geometryType))

            case (rectangular, trapezoidal, parabolic, triangular, rect_triang, rect_round, rectangular_closed, &
                    filled_circular, arch, semi_circular, circular, semi_elliptical, catenary, basket_handle,   &
                    horseshoe, gothic, eggshaped, horiz_ellipse, vert_ellipse, mod_basket, irregular)
                !% --- Copy all the geometry specific data from the adjacent element cell
                !%     Note that because irregular transect tables are not yet initialized, the
                !%     Area and Volume here will be junk for an irregular cross-section and will need to be
                !%     reset after transect tables are initialized. This occurs because we have
                !%     to cycle through all the CC, JM/JB before we can set the element transect
                !%     tables.
                elemR(JBidx,er_Area)                = elemR(Aidx,er_Area)[Ci]
                elemR(JBidx,er_AreaBelowBreadthMax) = elemR(Aidx,er_AreaBelowBreadthMax)[Ci]
                elemR(JBidx,er_BreadthMax)          = elemR(Aidx,er_BreadthMax)[Ci]
                elemR(JBidx,er_FullArea)            = elemR(Aidx,er_FullArea)[Ci]
                elemR(JBidx,er_FullDepth)           = elemR(Aidx,er_FullDepth)[Ci]
                !elemR(JBidx,er_FullHydDepth)        = elemR(Aidx,er_FullHydDepth)[Ci]
                elemR(JBidx,er_FullHydRadius)       = elemR(Aidx,er_FullHydRadius)[Ci]
                elemR(JBidx,er_FullPerimeter)       = elemR(Aidx,er_FullPerimeter)[Ci]
                elemR(JBidx,er_FullTopwidth)        = elemR(Aidx,er_FullTopwidth)[Ci]
                !% --- reference the Zbreadth max to the local bottom
                elemR(JBidx,er_ZbreadthMax)         = (elemR(Aidx,er_ZbreadthMax)[Ci] - elemR(Aidx,er_Zbottom)[Ci]) + elemR(JBidx,er_Zbottom)
                !% --- reference the Zcrown to the local bottom
                elemR(JBidx,er_Zcrown)              = (elemR(Aidx,er_Zcrown)[Ci] - elemR(Aidx,er_Zbottom)[Ci]) + elemR(JBidx,er_Zbottom)         
                elemR(JBidx,er_ManningsN)           = elemR(Aidx,er_ManningsN)[Ci]
                elemR(JBidx,er_ManningsN_Dynamic)   = elemR(Aidx,er_ManningsN)[Ci]
                elemI(JBidx,ei_link_transect_idx)   = elemI(Aidx,ei_link_transect_idx)[Ci]
                !% copy the entire row of the elemSGR array
                elemSGR(JBidx,:)                    = elemSGR(Aidx,:)[Ci]

            case (undefinedKey)
                print *, 'in ',trim(subroutine_name)
                print *, 'CODE ERROR: undefinedKey for ei_geometryType for junction'
                print *, 'at JBidx ',JBidx
                print * , ' '
                call util_crashpoint (23374)
                !return

            case default
                print *, 'in ',trim(subroutine_name)
                print *, 'CODE ERROR: unknown geometry type ',elemI(JBidx,ei_geometryType)
                print *, 'which has key ',trim(reverseKey(elemI(JBidx,ei_geometryType)))
                call util_crashpoint (4473)
                !return
            end select

            !% set the velocity
            if (elemR(JBidx,er_Area) .gt. setting%ZeroValue%Area) then ! BRHbugfix 20210813
                elemR(JBidx,er_Velocity) = elemR(JBidx,er_Flowrate) / elemR(JBidx,er_Area)
            else
                elemR(JBidx,er_Velocity) = zeroR
            end if

            !% Common geometry that do not depend on cross-section
            elemR(JBidx,er_Length)       = setting%Discretization%NominalElemLength
            elemR(JBidx,er_Area_N0)      = elemR(JBidx,er_Area)
            elemR(JBidx,er_Area_N1)      = elemR(JBidx,er_Area)
            elemR(JBidx,er_FullVolume)   = elemR(JBidx,er_FullArea)  * elemR(JBidx,er_Length)
            elemR(JBidx,er_Volume)       = elemR(JBidx,er_Area) * elemR(JBidx,er_Length)
            elemR(JBidx,er_Volume_N0)    = elemR(JBidx,er_Volume)
            elemR(JBidx,er_Volume_N1)    = elemR(JBidx,er_Volume)
            

            !print *, 'in initial', JBidx, elemR(JBidx,er_Length)

        end do
        
        !% --- set a JM length based on longest branches (20220711brh)
        LupMax = elemR(JMidx+1,er_Length) * real(elemSI(JMidx+1,esi_JunctionBranch_Exists),8)                              
        do ii=2,max_up_branch_per_node
            JBidx = JMidx + 2*ii - 1
            LupMax = max(elemR(JBidx,er_Length) * real(elemSI(JBidx,esi_JunctionBranch_Exists),8), LupMax)
        end do    
        LdnMax = elemR(JMidx+2,er_Length) * real(elemSI(JMidx+2,esi_JunctionBranch_Exists),8)  
        do ii=2,max_dn_branch_per_node
            JBidx = JMidx + 2*ii
            LdnMax = max(elemR(JBidx,er_Length) * real(elemSI(JBidx,esi_JunctionBranch_Exists),8), LdnMax)    
        end do
        elemR(JMidx,er_Length) = LupMax + LdnMax   
        
        !% --- get junction main geometry based on type
        JmType => elemSI(JMidx,esi_JunctionMain_Type)

        select case (JmType)

            case (ImpliedStorage)
                elemSR(JMidx,esr_Storage_Plan_Area) = setting%Junction%SurfaceArea_Minimum
                elemR (JMidx,er_FullVolume)         = elemSR(JMidx,esr_Storage_Plan_Area) * elemR(JMidx, er_FullDepth)
                elemR (JMidx,er_FullArea)           = elemSR(JMidx,esr_Storage_Plan_Area)
                elemR (JMidx,er_BreadthMax)         = sqrt(elemSR(JMidx,esr_Storage_Plan_Area) )
                elemR (JMidx,er_Length)             = sqrt(elemSR(JMidx,esr_Storage_Plan_Area) )

            case (NoStorage)
                elemSR(JMidx,esr_Storage_Plan_Area)      = zeroR
                elemR (JMidx,er_FullVolume)              = zeroR
                elemR (JMidx,er_FullArea)                = zeroR  
                elemR (JMidx,er_BreadthMax)              = zeroR     
                

                !% REMOVING 20230326
                    ! do ii=1,max_branch_per_node
                    !     JBidx = JMidx + ii
                    !     if (.not. elemSI(JBidx,esi_JunctionBranch_Exists) == oneI) cycle

                    !     BranchIdx      => elemSI(JBidx,esi_JunctionBranch_Link_Connection)
                    !     JBgeometryType => elemI(JBidx,ei_geometryType)

                    !     !% -- get the full area by summation of full area branches
                    !     elemR(JMidx,er_FullArea) = elemR(JMidx,er_FullArea) + init_IC_get_branch_fullarea(JBidx)

                    !     !% OBSOLETE: HACK -- needs to be removed after new junction proven 20230208
                    !     !% --- compute the storage plan area by accumulation from branches
                    !     select case (JBgeometryType)
                    !         case (rectangular,rectangular_closed)
                    !             elemSR(JMidx,esr_Storage_Plan_Area) = elemSR(JMidx,esr_Storage_Plan_Area)  &
                    !             +(real(elemSI( JBidx,esi_JunctionBranch_Exists),8)                       &
                    !                 * elemR(  JBidx,er_Length)                                          &
                    !                 * elemSGR(JBidx,esgr_Rectangular_Breadth) )

                    !         case (trapezoidal)
                    !             elemSR(JMidx,esr_Storage_Plan_Area) = elemSR(JMidx,esr_Storage_Plan_Area)  &
                    !                     +(real(elemSI( JBidx,esi_JunctionBranch_Exists),8)                  &
                    !                         * elemR(  JBidx,er_Length)                                     &
                    !                         * ( elemSGR(JBidx,esgr_Trapezoidal_Breadth)                    &
                    !                             +elemSGR(JBidx,esgr_Trapezoidal_LeftSlope)                  &
                    !                             +elemSGR(JBidx,esgr_Trapezoidal_RightSlope)))

                    !         case (triangular)
                    !             elemSR(JMidx,esr_Storage_Plan_Area) = elemSR(JMidx,esr_Storage_Plan_Area)  &
                    !                         +(real(elemSI( JBidx,esi_JunctionBranch_Exists),8)               &
                    !                             * elemR(  JBidx,er_Length)                                   &
                    !                             * (elemSGR(JBidx,esgr_Triangular_TopBreadth)/twoR) )
                    !         case (parabolic)
                    !             elemSR(JMidx,esr_Storage_Plan_Area) = elemSR(JMidx,esr_Storage_Plan_Area)  &
                    !             +(real(elemSI( JBidx,esi_JunctionBranch_Exists),8)                       &
                    !                 * elemR(  JBidx,er_Length)                                          &
                    !                 * elemSGR(JBidx,esgr_Parabolic_Breadth) )

                    !         case (rect_triang)
                    !             elemSR(JMidx,esr_Storage_Plan_Area) = elemSR(JMidx,esr_Storage_Plan_Area)  &
                    !                         +(real(elemSI(JBidx,esi_JunctionBranch_Exists),8)               &
                    !                             * elemR(  JBidx,er_Length)                                   &
                    !                             * (elemR( JBidx,er_BreadthMax)/twoR) )

                    !         case (basket_handle)
                    !             elemSR(JMidx,esr_Storage_Plan_Area) = elemSR(JMidx,esr_Storage_Plan_Area)  &
                    !                         +(real(elemSI( JBidx,esi_JunctionBranch_Exists),8)               &
                    !                             * elemR(  JBidx,er_Length)                                   &
                    !                             * (elemR(JBidx,er_BreadthMax)/twoR) )
                            
                    !         case (arch)
                    !             elemSR(JMidx,esr_Storage_Plan_Area) = elemSR(JMidx,esr_Storage_Plan_Area)  &
                    !                         +(real(elemSI( JBidx,esi_JunctionBranch_Exists),8)               &
                    !                             * elemR(  JBidx,er_Length)                                   &
                    !                             * (elemR(JBidx,er_BreadthMax)/twoR) )

                    !         case (horiz_ellipse)
                    !             elemSR(JMidx,esr_Storage_Plan_Area) = elemSR(JMidx,esr_Storage_Plan_Area)  &
                    !                         +(real(elemSI( JBidx,esi_JunctionBranch_Exists),8)               &
                    !                             * elemR(  JBidx,er_Length)                                   &
                    !                             * (elemR(JBidx,er_BreadthMax)/twoR) )
                            
                    !         case (vert_ellipse)
                    !             elemSR(JMidx,esr_Storage_Plan_Area) = elemSR(JMidx,esr_Storage_Plan_Area)  &
                    !                         +(real(elemSI( JBidx,esi_JunctionBranch_Exists),8)               &
                    !                             * elemR(  JBidx,er_Length)                                   &
                    !                             * (elemR(JBidx,er_BreadthMax)/twoR) )

                    !         case (eggshaped)
                    !             elemSR(JMidx,esr_Storage_Plan_Area) = elemSR(JMidx,esr_Storage_Plan_Area)  &
                    !                         +(real(elemSI( JBidx,esi_JunctionBranch_Exists),8)               &
                    !                             * elemR(  JBidx,er_Length)                                   &
                    !                             * (elemR(JBidx,er_BreadthMax)/twoR) )

                    !         case (horseshoe)
                    !             elemSR(JMidx,esr_Storage_Plan_Area) = elemSR(JMidx,esr_Storage_Plan_Area)  &
                    !                         +(real(elemSI( JBidx,esi_JunctionBranch_Exists),8)               &
                    !                             * elemR(  JBidx,er_Length)                                   &
                    !                             * (elemR(JBidx,er_BreadthMax)/twoR) )
                    !         case (catenary)
                    !             elemSR(JMidx,esr_Storage_Plan_Area) = elemSR(JMidx,esr_Storage_Plan_Area)  &
                    !                         +(real(elemSI( JBidx,esi_JunctionBranch_Exists),8)               &
                    !                             * elemR(  JBidx,er_Length)                                   &
                    !                             * (elemR(JBidx,er_BreadthMax)/twoR) )

                    !         case (gothic)
                    !             elemSR(JMidx,esr_Storage_Plan_Area) = elemSR(JMidx,esr_Storage_Plan_Area)  &
                    !                         +(real(elemSI( JBidx,esi_JunctionBranch_Exists),8)               &
                    !                             * elemR(  JBidx,er_Length)                                   &
                    !                             * (elemR(JBidx,er_BreadthMax)/twoR) )
                            
                    !         case (semi_elliptical)
                    !             elemSR(JMidx,esr_Storage_Plan_Area) = elemSR(JMidx,esr_Storage_Plan_Area)  &
                    !                         +(real(elemSI( JBidx,esi_JunctionBranch_Exists),8)               &
                    !                             * elemR(  JBidx,er_Length)                                   &
                    !                             * (elemR(JBidx,er_BreadthMax)/twoR) )
                                                
                    !         case (circular,force_main)
                    !             elemSR(JMidx,esr_Storage_Plan_Area) = elemSR(JMidx,esr_Storage_Plan_Area)  &
                    !             +(real(elemSI( JBidx,esi_JunctionBranch_Exists),8)                       &
                    !                 * elemR(  JBidx,er_Length)                                          &
                    !                 * (elemR(JBidx,er_BreadthMax)/twoR) )

                    !         case (semi_circular)
                    !             elemSR(JMidx,esr_Storage_Plan_Area) = elemSR(JMidx,esr_Storage_Plan_Area)  &
                    !                         +(real(elemSI( JBidx,esi_JunctionBranch_Exists),8)               &
                    !                             * elemR(  JBidx,er_Length)                                   &
                    !                             * (elemR(JBidx,er_BreadthMax)/twoR) )
                            
                    !         case (filled_circular)
                    !             elemSR(JMidx,esr_Storage_Plan_Area) = elemSR(JMidx,esr_Storage_Plan_Area)  &
                    !             +(real(elemSI( JBidx,esi_JunctionBranch_Exists),8)                          &
                    !                 * elemR(  JBidx,er_Length)                                             &
                    !                 * elemR(  JBidx,er_BreadthMax) )
                            
                    !         case (mod_basket)
                    !             elemSR(JMidx,esr_Storage_Plan_Area) = elemSR(JMidx,esr_Storage_Plan_Area)  &
                    !             +(real(elemSI( JBidx,esi_JunctionBranch_Exists),8)                          &
                    !                 * elemR(  JBidx,er_Length)                                             &
                    !                 * elemR(  JBidx,er_BreadthMax) )
                            
                    !         case (rect_round)
                    !             elemSR(JMidx,esr_Storage_Plan_Area) = elemSR(JMidx,esr_Storage_Plan_Area)  &
                    !                         +(real(elemSI( JBidx,esi_JunctionBranch_Exists),8)               &
                    !                             * elemR(  JBidx,er_Length)                                   &
                    !                             * (elemR(JBidx,er_BreadthMax)/twoR) )

                    !         case (irregular)    
                    !             !% --- for irregular geometry, we use the average breadth for the area below
                    !             !%     the maximum breadth as the characteristic width of the branch
                    !             elemSR(JMidx,esr_Storage_Plan_Area) = elemSR(JMidx,esr_Storage_Plan_Area)  &
                    !             +(real(elemSI( JBidx,esi_JunctionBranch_Exists),8)                       &
                    !                 * elemR(  JBidx,er_Length)                                          &
                    !                 * elemR(  JBidx,er_AreaBelowBreadthMax) / elemR(JBidx,er_ZbreadthMax))

                    !         case default
                    !             print *, 'In, ', subroutine_name
                    !             print *, 'CODE ERROR: geometry type unknown for # ', JBgeometryType
                    !             print *, 'which has key ',trim(reverseKey(JBgeometryType))
                    !             !stop 
                    !             call util_crashpoint(308974)
                    !             !return
                    !     end select
                    ! end do



            case (FunctionalStorage)     
                !% --- create a storage curve from the function
                call storage_create_curve_from_function (JMidx)
                !% --- get this curveID
                CurveID => elemSI(JMidx,esi_JunctionMain_Curve_ID)
                elemR(JMidx,er_FullVolume) = maxval(curve(CurveID)%ValueArray(:,curve_storage_volume))
                !% --- Computing the Full (cross-sectional) Area for storage as
                !%     Area = FullDepth * AverageTopwidth
                !%          = FullDepth * sqrt(Average(PlanArea))
                !%          = FullDepth * sqrt(FullVolume/FullDepth)
                !%          = sqrt(FullDepth * Full Volume) 
                !%     which would be simply sqrt(FullVolume * FullDepth)
                !%     This would allow the preissmann slot to be set for surcharged
                !%      on storage elements.
                elemR(JMidx,er_FullArea)   = sqrt( elemR(JMidx,er_FullVolume) * elemR(JMidx,er_FullDepth) )
                !% --- max breadth approximated as sqrt of max planar area
                elemR(JMidx,er_BreadthMax) = sqrt(maxval(curve(CurveID)%ValueArray(:,curve_storage_area)))

                elemR(JMidx,er_Length) = sqrt(elemR(JMidx,er_FullArea))

            case (TabularStorage)
                CurveID => elemSI(JMidx,esi_JunctionMain_Curve_ID)
                !% --- set the element index for the curve
                Curve(CurveID)%ElemIdx = JMidx

                !% SWMM5+ needs a volume vs depth relationship thus Trapezoidal rule is used
                !% to get to integrate the area vs depth curve
                call storage_create_integrated_volume_curve (CurveID)
                
                !% --- set full values based on curve
                elemR(JMidx,er_FullVolume) = maxval(curve(CurveID)%ValueArray(:,curve_storage_volume))
                !% --- see note in Functional Storage
                elemR(JMidx,er_FullArea)   = sqrt( elemR(JMidx,er_FullVolume) * elemR(JMidx,er_FullDepth) )
                !% --- max breadth approximated as sqrt of max planar area
                elemR(JMidx,er_BreadthMax) = sqrt(maxval(curve(CurveID)%ValueArray(:,curve_storage_area)))

                elemR(JMidx,er_Length) = sqrt(elemR(JMidx,er_FullArea))

            case default
                !% IMPORTANT -- if any other new type is defined, make sure that
                !% subroutine geo_depth_from_volume is updated
                print *, 'In, ', subroutine_name
                print *, 'CODE ERROR junction main type unknown for # ', JmType
                print *, 'which has key ',trim(reverseKey(JmType))
                !stop 
                call util_crashpoint(54895)
                !return

        end select

        !% -- initial conditions volume
        elemR(JMidx,er_Volume)     = storage_volume_from_depth_singular (JMidx,elemR(JMidx,er_Depth))  
        elemR(JMidx,er_Volume_N0)  = elemR(JMidx,er_Volume)
        elemR(JMidx,er_Volume_N1)  = elemR(JMidx,er_Volume)

        !% ---initial conditions plan storage areas
        select case (JmType)
            case (NoStorage, ImpliedStorage)
                !% -- already specified as constant
            case (FunctionalStorage,TabularStorage)
                call util_curve_lookup_singular(CurveID, er_Volume, er_Temp01, curve_storage_volume, &
                                                curve_storage_area, 1)
                elemSR(JMidx,esr_Storage_Plan_Area) = elemR(JMidx,er_Temp01)                          
            case default 
                print *, 'CODE ERROR: Unexpected case default'
                call util_crashpoint(6098734)
        end select


        
        !print *, 'JM ',JMidx, elemR(JMidx,er_BreadthMax)
        !stop 2397840

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_get_junction_data   
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_IC_JB_nullvalue_geometry  &
         (Aidx, Ci, thisJunctionNode, JBidx, isupstream)
        !%------------------------------------------------------------------
        !% Description: 
        !% handles cases where JB is adjacent to a diagnostic element
        !% without inherently-defined geometry
        !% Returns the Aidx and Ci of an element whose geometry can be used
        !% for inferring geometry of JB
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(inout) :: Aidx !% adjacent element index
            integer, intent(inout) :: Ci   !% adjacent element connected image 
            integer, intent(in)    :: thisJunctionNode !% node being handled
            integer, intent(in)    :: JBidx !% junction branch being handled
            logical, intent(in)    :: isupstream !% if JB is an upstream branch
            integer                :: adjLink, nextNode, farLink

            character(64)  :: subroutine_name = 'init_IC_JB_nullvalue_geometry'
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------      
            ! print *, 'in ',trim(subroutine_name)

        !% --- define the adjacent link
        adjLink = elemI(Aidx,ei_link_Gidx_BIPquick)[Ci]
        
            ! print *, 'adjLink ',adjLink

        !% --- DOWNSTREAM INFERENCE -----------------------------------
        if (.not. isupstream) then 
            !% --- get the next downstream node
            nextnode = link%I(adjLink,li_Mnode_d)

                ! print *, 'next node dn',nextnode
                ! print *, 'node name    ',trim(node%Names(nextnode)%str)

            !% --- check if only one link connected downstream
            if (node%I(nextnode,ni_N_link_d) == 1) then 
                farLink = node%I(nextnode,ni_Mlink_d1)

                    ! print *, 'far link dn',farLink
                    ! print *, 'link name   ', trim(link%Names(farLink)%str)
                    ! print *, 'link type   ', link%I(farLink,li_link_type)
                    ! print *, 'reverseKey  ',reverseKey(link%I(farLink,li_link_type))

                !% --- check if link type can be used to infer geometry    
                if ((link%I(farLink,li_link_type) .eq. lPipe) .or. &
                    (link%I(farLink,li_link_type) .eq. lChannel)) then 
                    !% --- set the connected image and adjacent element to
                    !%     the far link to use for JB geometry
                    Ci   = link%I(farLink,li_P_image)
                    Aidx = link%I(farLink,li_last_elem_idx)
                    elemI(JBidx,ei_geometryType) = elemI(Aidx,ei_geometryType)[Ci]   
                else
                    !% --- far link cannot be used because wrong type
                    !%     set to null
                    Ci   = nullvalueI
                    Aidx = nullvalueI
                end if
            else 
                !% --- far link cannot be used because more than 1 connection
                !%     set to null
                Ci   = nullvalueI
                Aidx = nullvalueI
            end if 

            !% --- in case a pipe/channel not found downstream of adjacent link
            if ((Ci == nullvalueI) .or. (Aidx == nullvalueI)) then
                !% -- check for a single link upstream that could be used
                !%    to assign geometry. Only applicable if there is
                !%    only 1 upstream link, otherwise we cannot infer a
                !%    geometry.
                if (node%I(thisJunctionNode,ni_N_link_u) == 1) then
                    !% --- get the upstream link
                    farLink = node%I(thisJunctionNode,ni_Mlink_u1)

                    !% --- check if link type can be used to infer geometry  
                    if ((link%I(farLink,li_link_type) .eq. lPipe) .or. &
                        (link%I(farLink,li_link_type) .eq. lChannel)) then 
                        !% --- set the connected image and adjacent element to
                        !%     the far link to use for JB geometry
                        Ci   = link%I(farLink,li_P_image)
                        Aidx = link%I(farLink,li_last_elem_idx)
                        elemI(JBidx,ei_geometryType) = elemI(Aidx,ei_geometryType)[Ci]  
                    else 
                        !% no change, Ci=nullvalueI
                    end if
                else 
                    !% no change, Ci=nullvalueI
                end if
            else 
                !% no change, Ci and Aidx have been found    
            end if

        !% --- UPSTREAM INFERENCE ---------------------------------------
        else
            !% --- get the next upstream node
            nextnode = link%I(adjLink,li_Mnode_u)

                ! print *, 'next node up ',nextnode
                ! print *, 'node name    ',trim(node%Names(nextnode)%str)

            !% --- check if only one link connected upstream
            if (node%I(nextnode,ni_N_link_u) == 1) then 
                farLink = node%I(nextnode,ni_Mlink_u1)

                    ! print *, 'far link up ',farLink
                    ! print *, 'link name   ', trim(link%Names(farLink)%str)
                    ! print *, 'link type   ', link%I(farLink,li_link_type)
                    ! print *, 'reverseKey  ',reverseKey(link%I(farLink,li_link_type))

                !% --- check if link type can be used to infer geometry  
                if ((link%I(farLink,li_link_type) .eq. lPipe) .or. &
                    (link%I(farLink,li_link_type) .eq. lChannel)) then 
                    !% --- set the connected image and adjacent element to
                    !%     the far link to use for JB geometry
                    Ci   = link%I(farLink,li_P_image)
                    Aidx = link%I(farLink,li_last_elem_idx)
                    elemI(JBidx,ei_geometryType) = elemI(Aidx,ei_geometryType)[Ci]   
                else
                    !% --- far link cannot be used because wrong type
                    !%     set to null
                    Ci   = nullvalueI
                    Aidx = nullvalueI
                end if
            else 
                !% --- far link cannot be used becuase there are multiple links
                !%     set to null
                Ci   = nullvalueI
                Aidx = nullvalueI
            end if

            !% --- in case a pipe/channel not found upstream of adjacent link
            if ((Ci == nullvalueI) .or. (Aidx == nullvalueI)) then
                !% -- check for a single link downstream that could be used
                !%    to assign geometry. Only applicable if there is
                !%    only 1 downstream link, otherwise we cannot infer a
                !%    geometry.
                if (node%I(thisJunctionNode,ni_N_link_d) == 1) then
                    farLink = node%I(thisJunctionNode,ni_Mlink_d1)

                    !% --- check if link type can be used to infer geometry  
                    if ((link%I(farLink,li_link_type) .eq. lPipe) .or. &
                        (link%I(farLink,li_link_type) .eq. lChannel)) then 
                        !% --- set the connected image and adjacent element to
                        !%     the far link to use for JB geometry
                        Ci   = link%I(farLink,li_P_image)
                        Aidx = link%I(farLink,li_last_elem_idx)
                        elemI(JBidx,ei_geometryType) = elemI(Aidx,ei_geometryType)[Ci]  
                    else 
                        !% no change, Ci=nullvalueI  
                    end if
                else 
                    !% no change, Ci=nullvalueI  
                end if
            else 
                !% no change, Ci and Aidx have been found
            end if
        end if


        !% --- check for error remaining:
        if ((Ci == nullvalueI) .or. (Aidx == nullvalueI)) then
            print *, 'USER SYSTEM CONFIGURATION ERROR'
            print *, 'located at node index ',thisJunctionNode,' named: ',trim(node%Names(thisJunctionNode)%str)
            if (isupstream) then 
                print *, 'with the upstream link index   ',adjLink,' named: ',trim(link%Names(adjLink)%str)
            else 
                print *, 'with the downstream link index ',adjLink,' named: ',trim(link%Names(adjLink)%str)
            end if
            print *, 'PROBLEM: Cannot define geometry of the junction branch.'
            if ((node%I(thisJunctionNode,ni_N_link_d) + node%I(thisJunctionNode,ni_N_link_u)) == 2) then
                print *, 'SWMM5+ requires either a channel/conduit link connected upstream/downstream '
                print *, 'of this link or a channel/conduit link on the opposite side of the node'
                print *, '(e.g., the downstream side if this is an upstream link on the node).'

            else
                print *, 'SWMM5+ requires a channel/conduit link connected upstream/downstream to this link.'
            end if
            print *, 'This configuration is required to set implied geometry of junction branches'
            call util_crashpoint(6798723)
        end if

    end subroutine init_IC_JB_nullvalue_geometry
!%
!%==========================================================================
!%==========================================================================
!%    
    real(8) function init_IC_get_branch_fullarea (JBidx) result(outvalue)  
        !%------------------------------------------------------------------
        !% Description
        !% gets the full area for a branch if it exists
        !%------------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: JBidx
        !%------------------------------------------------------------------
            outvalue = (real(elemSI( JBidx,esi_JunctionBranch_Exists),8) &
                       * elemR(  JBidx,er_FullArea)) 

    end function init_IC_get_branch_fullarea
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_IC_elem_transect_arrays ()
        !%------------------------------------------------------------------
        !% Description:
        !% initializes the transect tables for elements (as opposed to the
        !% link transect tables, which are already stored). Note that this
        !% must be done after both CC and JM/JB element geometry are 
        !% initialized so that we have a count of all the elements with
        !% irregular geometry.
        !% This presumes that the ei_link_transect_idx has already been 
        !% assigned
        !%------------------------------------------------------------------
        !% Declarations:
            integer, pointer :: geometryType(:)
            integer, pointer :: elemTransectIdx(:), linkTransectIdx(:)
            integer          :: ii, thisTransectIdx

            character(64) :: subroutine_name = 'init_IC_elem_transect_arrays'
        !%------------------------------------------------------------------
        !% Aliases:
            geometryType    => elemI(:,ei_geometryType)
            elemTransectIdx => elemI(:,ei_transect_idx)
            linkTransectIdx => elemI(:,ei_link_transect_idx)
        !%------------------------------------------------------------------
        !% Preliminaries:
            !% --- count the total number of elements with irregular cross-sections
            N_transect = count(elemI(:,ei_geometryType) .eq. irregular)
            if (N_transect .le. 0) return
        !%------------------------------------------------------------------

        !% --- allocate the element transect arrays
        call util_allocate_element_transect ()

        thisTransectIdx = 0
        !% --- cycle through the elements
        do ii=1,N_elem(this_image())

            if (geometryType(ii) .ne. irregular) cycle
            !% --- increment the element transect index
            thisTransectIdx = thisTransectIdx +1
            !% --- store the element transect index
            elemTransectIdx(ii) = thisTransectIdx

            !% --- store the link transect data for each element
            transectTableDepthR(thisTransectIdx,:,:) = link%transectTableDepthR(linkTransectIdx(ii),:,:)
            transectTableAreaR (thisTransectIdx,:,:) = link%transectTableAreaR (linkTransectIdx(ii),:,:)
            transectI(thisTransectIdx,:)             = link%transectI          (linkTransectIdx(ii),:)
            transectR(thisTransectIdx,:)             = link%transectR          (linkTransectIdx(ii),:)
            transectID(thisTransectIdx)              = link%transectID         (linkTransectIdx(ii))%str
        end do
        
        !% --- HACK: we should develop an approach to allow smoothing of irregular geometry
        !%     tables where 2 links connect and there is an abrupt change of geometry that isn't
        !%     really supported by the data. This should be an optional algorithm that 
        !%     looks at the size of the adjacent links and smooths some fraction of the
        !%     elements on either side of the transition. Note that this should only bee
        !%     applied across a nj2 node (face) that joins exactly 2 links without storage.

        !%------------------------------------------------------------------
    end subroutine init_IC_elem_transect_arrays
!%
!%==========================================================================
!%==========================================================================
!%    
    subroutine init_IC_elem_transect_geometry ()
        !%------------------------------------------------------------------
        !% Description:
        !% initializes the time=0 area, volume etc. that depend on the 
        !% initial depth for irregular cross-sections. 
        !% This is delayed from all the other geometry initialization because
        !% the transect tables must be setup prior to computing initial
        !% geometry
        !%------------------------------------------------------------------
        !% Declarations
            integer, allocatable :: tpack(:)
            integer :: npack, ii
            real(8), pointer :: area(:), area0(:), area1(:), fullarea(:)
            real(8), pointer :: depth(:), fulldepth(:), length(:), ellDepth(:)
            real(8), pointer :: hydDepth(:), hydRadius(:),  fullHydRadius(:)
            real(8), pointer :: topwidth(:), fullTopWidth(:), perimeter(:)
            real(8), pointer :: volume(:), volume0(:), volume1(:)
            real(8), pointer :: thisTable(:,:)
        !%------------------------------------------------------------------
        !% Aliases:
            area          => elemR(:,er_Area)
            area0         => elemR(:,er_Area_N0)
            area1         => elemR(:,er_Area_N1)
            fullarea      => elemR(:,er_FullArea)
            depth         => elemR(:,er_Depth)
            fulldepth     => elemR(:,er_FullDepth)
            ellDepth      => elemR(:,er_EllDepth)
            hydRadius     => elemR(:,er_HydRadius)
            fullHydRadius => elemR(:,er_FullHydRadius)
            length        => elemR(:,er_Length)
            perimeter     => elemR(:,er_Perimeter)
            topwidth      => elemR(:,er_Topwidth)
            fullTopWidth  => elemR(:,er_FullTopWidth)
            volume        => elemR(:,er_Volume)
            volume0       => elemR(:,er_Volume_N0)
            volume1       => elemR(:,er_Volume_N1)
        !%------------------------------------------------------------------
        !% Preliminaries
            npack = count(elemI(:,ei_geometryType) .eq. irregular)
            if (npack .eq. 0) return
        !%------------------------------------------------------------------
        !% --- get the packed irregular element list
        tpack = pack(elemI(:,ei_Lidx), elemI(:,ei_geometryType) .eq. irregular)

        !% --- temporary store of the normalized depth
        depth(tpack) = depth(tpack) / fulldepth(tpack)

        !% --- get first guess at normalized area from normalized depth
        !%     NOTE that the table that produces depth from area is used
        !%     in the time-march and is NOT exactly invertible with the
        !%     area from depth table. Thus, we have a two-step process
        !%     for initial conditions to get the area that provides the
        !%     area in the area-to-depth table for the initial depth.
        thisTable => transectTableDepthR(:,:,tt_area)
        call xsect_table_lookup_array (area, depth, thisTable, tpack) 

        !% --- Find the area associated with the initial depth in the
        !%     depth-from-area lookup table
        thisTable => transectTableAreaR(:,:,tt_depth)
        do ii=1,size(tpack)
            !% --- get the area that matches the depth
            area(tpack(ii)) = xsect_find_x_matching_y (depth(tpack(ii)), thisTable(elemI(tpack(ii),ei_transect_idx),:) )
        end do
        !% --- Recompute the depth from area
        call xsect_table_lookup_array (depth, area, thisTable, tpack)

        !% --- get physical area
        area(tpack) = area(tpack) * fullarea(tpack)

        !% --- get normalized hydraulic radius from normalized depth
        thisTable => transectTableDepthR(:,:,tt_hydradius)
        call xsect_table_lookup_array (hydRadius, depth, thisTable, tpack) 

        !% --- get physical hydraulic radius
        hydRadius(tpack) = hydRadius(tpack) * fullHydRadius(tpack)

        !% --- get normalized top width from normalized depth
        thisTable => transectTableDepthR(:,:,tt_width)
        call xsect_table_lookup_array (topwidth, depth, thisTable, tpack) 

        !% --- get physical topwidth
        topwidth(tpack) = topwidth(tpack) * fullTopWidth(tpack)

        !% --- restore the physical depth
        depth(tpack) = depth(tpack) * fulldepth(tpack)



        !% --- derived data
        area0(tpack)     = area(tpack)
        area1(tpack)     = area(tpack)
        !hydDepth(tpack)  = area(tpack) / topwidth(tpack)
        !ell(tpack)       = hydDepth(tpack)  !% HACK -- assumes x-sect is continually increasing in width with depth
        perimeter(tpack) = area(tpack) / hydRadius(tpack)
        volume(tpack)    = area(tpack) * length(tpack)
        volume0(tpack)   = volume(tpack)
        volume1(tpack)   = volume(tpack)


        
        !%------------------------------------------------------------------
        !% Closing:
            deallocate(tpack)

    end subroutine init_IC_elem_transect_geometry
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_IC_derived_data ()
        !%------------------------------------------------------------------
        !% Description:
        !% Initial conditions for data derived from data already read from
        !% the input file
        !%------------------------------------------------------------------
        !% Declarations
            logical, pointer :: isSmallVol(:)
            integer, pointer :: npack, thisP(:)
            real(8), pointer :: area(:), flowrate(:), velocity(:)
            integer          :: thisCol
        !%------------------------------------------------------------------
        !% Preliminaries
        !%------------------------------------------------------------------
        !% Aliases
            if (setting%SmallDepth%useMomentumCutoffYN) then
                thisCol    = ep_CC_NOTsmalldepth
            else
                thisCol    = ep_CC_NOTzerodepth
            end if
            npack      => npack_elemP(thisCol)
            thisP      => elemP(1:npack,thisCol)
            area       => elemR(:,er_Area)
            flowrate   => elemR(:,er_Flowrate)
            velocity   => elemR(:,er_Velocity)
        !%------------------------------------------------------------------
        
        elemR(:,er_Velocity) = zeroR
        elemR(:,er_GammaM) = zeroR
        elemR(:,er_GammaC) = zeroR
        faceR(:,fr_GammaM) = zeroR

        if (npack < 1) return
        velocity(thisP) = flowrate(thisP) / area(thisP)

        where (velocity(thisP) > setting%Limiter%Velocity%Maximum)
           ! velocity(thisP) = zeroR
            velocity(thisP) = 0.99d0
        end where
        
        elemR(:,er_Velocity_N0) = velocity
        elemR(:,er_Velocity_N1) = velocity

        !%------------------------------------------------------------------
        !% Closing
    end subroutine init_IC_derived_data
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_IC_solver_select (whichSolver)
        !%------------------------------------------------------------------
        !% Desscription
        !% select the solver based on depth for all the elements
        !%------------------------------------------------------------------
        !% Declarations
            integer, intent(in) :: whichSolver
            character(64)       :: subroutine_name = 'init_IC_solver_select'
        !%------------------------------------------------------------------
        !% Preliminaries:
            !if (crashYN) return
            if (setting%Debug%File%initial_condition) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        select case (whichSolver)
        case (ETM)
            where ( (elemI(:,ei_HeqType) == time_march) .or. &
                    (elemI(:,ei_QeqType) == time_march) )
                elemI(:,ei_tmType) = ETM
            endwhere
        case (AC)
            print*, 'In, ', subroutine_name
            print*, 'AC solver is not handeled at this moment'
            !stop 
            call util_crashpoint(83974)
            !return
        case (ETM_AC)
            print*, 'In, ', subroutine_name
            print*, 'ETM-AC solver is not handeled at this moment'
            !stop 
            call util_crashpoint(2975)
            !return
        case default
            print *, 'In, ', subroutine_name
            print *, 'CODE ERROR: unknown solver, ', whichSolver
            print *, 'which has key ',trim(reverseKey(whichSolver))
            !stop 
            call util_crashpoint(81878)
            !return
        end select

        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%initial_condition) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_solver_select
!
!==========================================================================
!==========================================================================
!
    subroutine init_IC_small_values_diagnostic_elements ()
        !--------------------------------------------------------------------------
        !
        !% set the volume, area, head, other geometry, and flow to zero values
        !% for the diagnostic elements so no error is induced in the primary
        !% face update
        !
        !--------------------------------------------------------------------------

            character(64)       :: subroutine_name = 'init_IC_small_values_diagnostic_elements'

        !--------------------------------------------------------------------------
            !if (crashYN) return
            if (setting%Debug%File%initial_condition) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        where ( (elemI(:,ei_QeqType) == diagnostic) .or. (elemI(:,ei_HeqType) == diagnostic))
            !% HACK: settings%ZeroValues should be used here
            !% when the code is finalized
            elemR(:,er_Area)     = setting%ZeroValue%Area  ! 1.0e-6
            elemR(:,er_Topwidth) = setting%ZeroValue%Topwidth !1.0e-6
            !elemR(:,er_HydDepth) = setting%ZeroValue%Depth  ! 1.0e-6
            elemR(:,er_EllDepth) = setting%ZeroValue%Depth  ! 1.0e-6
            elemR(:,er_Flowrate) = zeroR
            elemR(:,er_Head)     = setting%ZeroValue%Depth + elemR(:,er_Zbottom) !1.0e-6
        endwhere

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_small_values_diagnostic_elements
!
!==========================================================================
!==========================================================================
!
    subroutine init_IC_diagnostic_interpolation_weights ()
        !--------------------------------------------------------------------------
        !
        !% set the interpolation weights for diagnostic elements
        !
        !--------------------------------------------------------------------------

            character(64)       :: subroutine_name = 'init_IC_diagnostic_interpolation_weights'

            integer, pointer ::  Npack, thisP(:), tM
            integer :: ii, kk, tB
        !--------------------------------------------------------------------------
            !if (crashYN) return
            if (setting%Debug%File%initial_condition) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"


        !% Q-diagnostic elements will have minimum interp weights for Q
        !% and maximum interp values for G and H
        !% Theses serve to force the Q value of the diagnostic element to the faces
        !% the G and H values are obtained from adjacent elements.
        where (elemI(:,ei_QeqType) == diagnostic)
            elemR(:,er_InterpWeight_uQ) = setting%Limiter%InterpWeight%Minimum
            elemR(:,er_InterpWeight_dQ) = setting%Limiter%InterpWeight%Minimum
            elemR(:,er_InterpWeight_uG) = setting%Limiter%InterpWeight%Maximum
            elemR(:,er_InterpWeight_dG) = setting%Limiter%InterpWeight%Maximum
            elemR(:,er_InterpWeight_uH) = setting%Limiter%InterpWeight%Maximum
            elemR(:,er_InterpWeight_dH) = setting%Limiter%InterpWeight%Maximum
            elemR(:,er_InterpWeight_uP) = setting%Limiter%InterpWeight%Maximum
            elemR(:,er_InterpWeight_dP) = setting%Limiter%InterpWeight%Maximum
        endwhere

        !% H-diagnostic elements will have minimum interp weights for H and G
        !% and maximum interp weights for Q
        !% These serve to force the G, and H values of the diagnostic element to the faces
        !% and the Q value is obtained from adjacent elements
        where (elemI(:,ei_HeqType) == diagnostic)
            elemR(:,er_InterpWeight_uQ) = setting%Limiter%InterpWeight%Maximum
            elemR(:,er_InterpWeight_dQ) = setting%Limiter%InterpWeight%Maximum
            elemR(:,er_InterpWeight_uG) = setting%Limiter%InterpWeight%Minimum
            elemR(:,er_InterpWeight_dG) = setting%Limiter%InterpWeight%Minimum
            elemR(:,er_InterpWeight_uH) = setting%Limiter%InterpWeight%Minimum
            elemR(:,er_InterpWeight_dH) = setting%Limiter%InterpWeight%Minimum
            elemR(:,er_InterpWeight_uP) = setting%Limiter%InterpWeight%Maximum
            elemR(:,er_InterpWeight_dP) = setting%Limiter%InterpWeight%Maximum
        endwhere

        ! !% brh 20220204 -- ccommenting this approach
        ! !% Branch elements have invariant interpolation weights so are computed here
        ! !% These are designed so that the face of a JB gets the flowrate from the
        ! !% adjacent CC conduit or channel, but the geometry and head are from the JB.
        ! Npack => npack_elemP(ep_JM_ALLtm)
        ! if (Npack > 0) then
        !     thisP  => elemP(1:Npack,ep_JM_ALLtm)
        !     do ii=1,Npack
        !         tM => thisP(ii) !% junction main ID
        !         do kk=1,max_branch_per_node
        !             tB = tM + kk !% junction branch ID
        !             elemR(tB,er_InterpWeight_uQ) = setting%Limiter%InterpWeight%Maximum
        !             elemR(tB,er_InterpWeight_dQ) = setting%Limiter%InterpWeight%Maximum
        !             elemR(tB,er_InterpWeight_uG) = setting%Limiter%InterpWeight%Minimum
        !             elemR(tB,er_InterpWeight_dG) = setting%Limiter%InterpWeight%Minimum
        !             elemR(tB,er_InterpWeight_uH) = setting%Limiter%InterpWeight%Minimum
        !             elemR(tB,er_InterpWeight_dH) = setting%Limiter%InterpWeight%Minimum
        !         end do
        !     end do
        ! end if

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_diagnostic_interpolation_weights

!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_IC_branch_dummy_values ()    
        !%------------------------------------------------------------------
        !% Description:
        !% assigns dummy values that non-zero and not excessivley large 
        !% to the velocity, depth, area, etc. of non-valid junction branches
        !% This allows these branches to be used in array computations 
        !% without causing either divide by zero or overflow/underflow.
        !%------------------------------------------------------------------
            integer, pointer :: npack, thisP(:), BranchExists(:)
            integer :: thisCol, ii
        !%------------------------------------------------------------------
            thisCol = ep_JM
            npack   => npack_elemP(thisCol)
            if (npack < 1) return
            thisP         => elemP(1:npack,thisCol)
            BranchExists  => elemSI(:,esi_JunctionBranch_Exists)
        !%------------------------------------------------------------------

        do ii=1,max_branch_per_node
            where (BranchExists(thisP+ii) .ne. oneI)
                elemR(thisP+ii,er_Area) = zeroR
                elemR(thisP+ii,er_Depth) = zeroR
                elemR(thisP+ii,er_Head) = zeroR
                elemR(thisP+ii,er_Velocity) = zeroR
                elemR(thisP+ii,er_Volume) = zeroR
            end where
        end do

    end subroutine init_IC_branch_dummy_values
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_IC_branch_zero_values ()
        !%------------------------------------------------------------------
        !% Description:
        !% assigns zero to _JB values as IC.
        !%------------------------------------------------------------------
            integer, pointer :: npack, thisP(:)
            integer :: kk
        !%------------------------------------------------------------------
            npack   => npack_elemP(ep_JM)
            if (npack < 1) return
            thisP     => elemP(1:npack,ep_JM)
        !%------------------------------------------------------------------
        ! do kk = 1,max_branch_per_node
        !     elemR(thisP+kk,er_2B_psiL)    = zeroR
        !     elemR(thisP+kk,er_EnergyHead) = zeroR
        ! end do

    end subroutine init_IC_branch_zero_values
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_IC_set_SmallVolumes ()
        !%------------------------------------------------------------------
        !% set the small volume values in elements that are used with
        !% the small depth cutoff to change the time-marching solution
        !% The ratio of the (fixed)small volume to the actual volume is 
        !% used to blend the CM small-volume solution with the SVE solution.
        !%
        !% NOTE: tpack is a 2D array although only 1D is used. This is
        !% to ensure that it is compatible with the elemPGx arrays used
        !% in llgeo_tabular_depth_from_column
        !%------------------------------------------------------------------
        !% Declarations:
            character(64)       :: subroutine_name = 'init_IC_set_SmallVolumes'
            !logical, pointer    :: useSmallVolumes
            real(8), pointer    :: MomentumDepthCutoff, smallVolume(:), length(:)
            real(8), pointer    :: theta(:), radius(:),  area(:)
            !real(8), pointer    :: trapB(:), trapL(:), trapR(:),rectB(:)
            real(8), pointer    :: depth(:)
            real(8), pointer    :: tempDepth(:), tempArea(:), Atable(:)
            integer, pointer    :: geoType(:), tPack(:,:), eIdx(:) 
            integer             :: npack, ii, indx, kk
            integer, dimension(11) :: tabXsectType
        !%------------------------------------------------------------------
        !% Preliminaries
            !if (crashYN) return
            if (setting%Debug%File%initial_condition) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% Aliases
            !useSmallVolumes  => setting%SmallDepth%UseSmallVolumesYN
            MomentumDepthCutoff      => setting%SmallDepth%MomentumDepthCutoff
            geoType          => elemI(:,ei_geometryType)
            !% --- note, tPack is 2D, but 2nd column is dummy
            tPack            => elemI(:,ei_Temp01:ei_Temp02)
            eIdx             => elemI(:,ei_Lidx)
            smallVolume      => elemR(:,er_SmallVolume)
            length           => elemR(:,er_Length)
            area             => elemR(:,er_Area)
            depth            => elemR(:,er_Depth)
            tempDepth        => elemR(:,er_Temp01)
            tempArea         => elemR(:,er_Temp02)
            !rectB            => elemSGR(:,esgr_Rectangular_Breadth)
            radius           => elemSGR(:,esgr_Circular_Radius)
            !trapL            => elemSGR(:,esgr_Trapezoidal_LeftSlope)
            !trapR            => elemSGR(:,esgr_Trapezoidal_RightSlope)
            !trapB            => elemSGR(:,esgr_Trapezoidal_Breadth)
            theta            => elemR(:,er_Temp03)
        !%------------------------------------------------------------------
        !% More preliminaries
            elemR(:,er_SmallVolume) = zeroR

            !% error checking for circular pipes
            theta = zeroR ! temporary use of theta space for comparison, this isn't theta!
            where (geoType == circular)
                theta = radius - MomentumDepthCutoff
            end where
            if (any(theta < zeroR)) then
                print *, 'FATAL ERROR'
                print *, 'Small Volume depth cutoff ',MomentumDepthCutoff
                print *, 'is larger or equal to radius of the smallest pipe '
                print *, 'Small Volume depth cutoff must be smaller than the radius.'
                !stop 
                call util_crashpoint(398705)
                !return
            end if
        !%------------------------------------------------------------------

        tabXsectType = (/ arch, basket_handle, catenary, circular, eggshaped, gothic, &
            horiz_ellipse, horseshoe, semi_circular, semi_elliptical, vert_ellipse /)

        !% --- temporarily store depth and area. 
        !%     Replace depth with the cutoff depth so that we can use standard 
        !%     depth-to-area-functions
        !%     Area is overwritten in tabular_are_from_depth needed for conduits.
        !%     These must be reversed at end of subroutine
        tempDepth = depth
        depth     = MomentumDepthCutoff
        tempArea  = area

        !% --- cycle through tabulated cross-sections
        do ii=1,size(tabXsectType)
            select case (tabXsectType(ii))
            case (arch)
                Atable => AArch
            case (basket_handle)
                Atable => ABasketHandle
            case (catenary)
                Atable => ACatenary
            case (circular)
                Atable => ACirc
            case (eggshaped)
                Atable => AEgg
            case (gothic)
                Atable => AGothic
            case (horiz_ellipse)
                Atable => AHorizEllip
            case (horseshoe)
                Atable => AHorseShoe
            case (semi_circular)
                Atable => ASemiCircular
            case (semi_elliptical)
                Atable => ASemiEllip
            case (vert_ellipse)
                Atable => AVertEllip
            end select

            tPack(:,1) = zeroI
            npack = count(geoType == arch)
            if (npack > 0) then
                tPack(1:npack,1) = pack(eIdx,geoType == tabXsectType(ii))
                !% --- get area associated with small volume
                call llgeo_tabular_area_from_depth(tpack, Npack,1, Atable, zeroR)
                !% --- small volume = A * length
                smallvolume(tPack(1:npack,1)) = area(tPack(1:npack,1)) &
                                                * length(tpack(1:npack,1))
            end if
        end do

        !% --- parabolic channel
        tPack(:,1) = zeroI
        npack = count(geoType == parabolic)
        if (npack > 0) then
            tPack(1:npack,1) = pack(eIdx,geoType == parabolic)
            smallvolume(tPack(1:npack,1)) = llgeo_parabolic_area_from_depth_pure      &
                                            (tPack(1:npack,1),depth(tPack(1:npack,1)) ) &
                                            * length(tPack(1:npack,1))
        end if

        !% --- power function channel
        tPack(:,1) = zeroI
        npack = count(geoType == power_function)
        if (npack > 0) then
            print *, 'error power function not completed'
            call util_crashpoint(6698723)
            tPack(1:npack,1) = pack(eIdx,geoType == power_function)
            !smallvolume(tPack(1:npack)) = llgeo_powerfunction_area_from_depth_pure(tPack(1:npack)) * length(tPack(1:npack))
        end if

        !% --- rectangular channel
        tPack(:,1) = zeroI
        npack = count(geoType == rectangular)
        if (npack > 0) then
            tPack(1:npack,1) = pack(eIdx,geoType == rectangular)
            smallvolume(tPack(1:npack,1)) = llgeo_rectangular_area_from_depth_pure &
                                            (tPack(1:npack,1),depth(tpack(1:npack,1))) &
                                            * length(tPack(1:npack,1))
        end if

        !% --- trapezoidal channel 
        tPack(:,1) = zeroI
        npack = count(geoType == trapezoidal)
        if (npack > 0) then
            tPack(1:npack,1) = pack(eIdx,geoType == trapezoidal)
            smallvolume(tPack(1:npack,1)) = llgeo_trapezoidal_area_from_depth_pure &
                                            (tPack(1:npack,1),depth(tpack(1:npack,1)))  &
                                            * length(tPack(1:npack,1))
        end if  

        !% --- triangular channel 
        tPack(:,1) = zeroI
        npack = count(geoType == triangular)
        if (npack > 0) then
            tPack(1:npack,1) = pack(eIdx,geoType == triangular)
            smallvolume(tPack(1:npack,1)) = llgeo_triangular_area_from_depth_pure &
                                            (tPack(1:npack,1),depth(tpack(1:npack,1)))  &
                                            * length(tPack(1:npack,1))
        end if 

        !% --- irregular channel 
        tPack(:,1) = zeroI
        npack = count(geoType == irregular)
        if (npack > 0) then
            tPack(1:npack,1) = pack(eIdx,geoType == irregular)
            do kk=1,npack
                smallvolume(tPack(kk,1)) = irregular_geometry_from_depth_singular &
                    (tpack(kk,1),tt_area, depth(tpack(kk,1)), elemR(tpack(kk,1),er_FullArea), setting%ZeroValue%Area)
            end do
        end if 

        !% ---- CLOSED CONDUITS

        !% ---  filled circular conduit
        tPack(:,1) = zeroI
        npack = count(geoType == filled_circular)
        if (npack > 0) then
            tPack(1:npack,1) = pack(eIdx,geoType == filled_circular)
            do kk=1,npack 
                smallvolume(tpack(kk,1)) = llgeo_filled_circular_area_from_depth_singular &
                                            (tpack(kk,1), depth(tpack(kk,1)), setting%ZeroValue%Area)
            end do
        end if
    
        !% ---  Modified basket conduit
        tPack(:,1) = zeroI
        npack = count(geoType == mod_basket)
        if (npack > 0) then
            tPack(1:npack,1) = pack(eIdx,geoType == mod_basket)
            do kk=1,npack 
                smallvolume(tpack(kk,1)) = llgeo_mod_basket_area_from_depth_singular &
                                            (tpack(kk,1), depth(tpack(kk,1)), setting%ZeroValue%Area)
            end do
        end if

        !% --- rectangular closed conduit
        tPack(:,1) = zeroI
        npack = count(geoType == rectangular_closed)
        if (npack > 0) then
            tPack(1:npack,1) = pack(eIdx,geoType == rectangular_closed)
            do kk=1,npack 
                smallvolume(tpack(kk,1)) = llgeo_rectangular_closed_area_from_depth_singular &
                                            (tpack(kk,1), depth(tpack(kk,1)), setting%ZeroValue%Area)
            end do
        end if

        !% ---  rectangular round conduit
        tPack(:,1) = zeroI
        npack = count(geoType == rect_round)
        if (npack > 0) then
            tPack(1:npack,1) = pack(eIdx,geoType == rect_round)
            do kk=1,npack
                smallvolume(tpack(kk,1)) = llgeo_rect_round_area_from_depth_singular &
                                            (tpack(kk,1), depth(tpack(kk,1)), setting%ZeroValue%Area)
            end do
        end if

        !% ---  rectangular triang conduit
        tPack(:,1) = zeroI
        npack = count(geoType == rect_triang)
        if (npack > 0) then
            tPack(1:npack,1) = pack(eIdx,geoType == rect_triang)
            do kk=1,npack 
                smallvolume(tpack(kk,1)) = llgeo_rectangular_triangular_area_from_depth_singular &
                                            (tpack(kk,1), depth(tpack(kk,1)), setting%ZeroValue%Area)
            end do
        end if
    
        !% ---  custom conduit
        tPack(:,1) = zeroI
        npack = count(geoType == custom)
        if (npack > 0) then
            tPack(1:npack,1) = pack(eIdx,geoType == custom)
            print *, 'PROBLEM IN SMALL VOLUMES FOR custom CONDUIT'
            call util_crashpoint(6298738)
        end if



        !% restore the initial condition depth to the depth and area vectors
        depth = tempDepth
        area  = tempArea
 
        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%initial_condition) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_set_SmallVolumes
!
!==========================================================================
!==========================================================================
!
    subroutine init_IC_set_zero_lateral_inflow ()
        !%-----------------------------------------------------------------
        !% set all the lateral inflows to zero before start of a simulation
        !%-----------------------------------------------------------------

        elemR(:,er_FlowrateLateral) = zeroR

    end subroutine init_IC_set_zero_lateral_inflow
!
!==========================================================================
!==========================================================================
!
    subroutine init_IC_oneVectors ()
        !%-----------------------------------------------------------------
        !% set up a vector of real ones (useful in sign functions)
        !%-----------------------------------------------------------------

        elemR(:,er_ones) = oneR

    end subroutine init_IC_oneVectors
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_IC_slot ()
        !%-----------------------------------------------------------------
        !% Description:
        !% initialize Preissmann Slot
        !% get the geometry data for conduit links and calculate element volumes
        !%-----------------------------------------------------------------
        !% Declarations:
            integer :: ii
            integer, pointer    :: SlotMethod
            real(8), pointer    :: TargetPCelerity, grav, Alpha
            character(64) :: subroutine_name = 'init_IC_slot'
        !%-----------------------------------------------------------------
            if (setting%Debug%File%initial_condition) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% pointer to geometry type
        SlotMethod          => setting%Solver%PreissmannSlot%Method
        TargetPCelerity     => setting%Solver%PreissmannSlot%TargetCelerity
        Alpha               => setting%Solver%PreissmannSlot%Alpha
        grav                => setting%Constant%gravity

        !% initialize slots
        elemR(1:size(elemR,1)-1,er_SlotVolume)            = zeroR
        elemR(1:size(elemR,1)-1,er_SlotArea)              = zeroR
        elemR(1:size(elemR,1)-1,er_SlotWidth)             = zeroR
        elemR(1:size(elemR,1)-1,er_dSlotArea)             = zeroR
        elemR(1:size(elemR,1)-1,er_dSlotDepth)            = zeroR
        elemR(1:size(elemR,1)-1,er_dSlotVolume)           = zeroR
        elemR(1:size(elemR,1)-1,er_SlotVolume_N0)         = zeroR
        elemR(1:size(elemR,1)-1,er_Preissmann_Celerity)   = zeroR
        elemR(1:size(elemR,1)-1,er_Surcharge_Time)        = zeroR      
        elemR(1:size(elemR,1)-1,er_SlotDepth_N0)          = elemR(1:size(elemR,1)-1,er_SlotDepth)    
        elemR(1:size(elemR,1)-1,er_Preissmann_Number_initial) = TargetPCelerity / (Alpha * sqrt(grav &
                                                              * elemR(1:size(elemR,1)-1,er_FullDepth)))
    
        !% only calculate slots for ETM time-march
        if (setting%Solver%SolverSelect == ETM) then

            !% --- initialization where starting condition is surcharged
            where (elemR(:,er_Head) > elemR(:,er_Zcrown))
                elemYN(:,eYN_isPSsurcharged) = .true.
                elemR (:,er_SlotDepth)      = elemR(:,er_Head) - elemR(:,er_Zcrown)
            endwhere

            !% --- initialize PS dependent variables
            select case (SlotMethod)

            case (StaticSlot)

                elemR(1:size(elemR,1)-1,er_Preissmann_Number) = oneR

                where (elemYN(:,eYN_isPSsurcharged))
                    elemR(:,er_Preissmann_Celerity) = TargetPCelerity / elemR(:,er_Preissmann_Number)
                    elemR(:,er_SlotWidth)           = (grav * elemR(:,er_FullArea)) / (elemR(:,er_Preissmann_Celerity)**2.0)
                    elemR(:,er_SlotArea)            = elemR(:,er_SlotDepth) * elemR(:,er_SlotWidth) 
                    elemR(:,er_SlotVolume)          = elemR(:,er_SlotArea) * elemR(:,er_Length)
                    !% --- add slot volume to total volume (which was set to full volume)
                    elemR(:,er_Volume)              = elemR(:,er_Volume) + elemR(:,er_SlotVolume)
                end where

            case (DynamicSlot)

                elemR(1:size(elemR,1)-1,er_Preissmann_Number)     = TargetPCelerity / (Alpha * sqrt(grav * elemR(1:size(elemR,1)-1,er_FullDepth)))
                elemR(1:size(elemR,1)-1,er_Preissmann_Number_N0)  = elemR(1:size(elemR,1)-1,er_Preissmann_Number)
                where (elemYN(:,eYN_isPSsurcharged))
                    elemR(:,er_Preissmann_Celerity) = TargetPCelerity / elemR(:,er_Preissmann_Number)
                    elemR(:,er_SlotWidth)           = (grav * elemR(:,er_FullArea)) / (elemR(:,er_Preissmann_Celerity)**2.0)
                    elemR(:,er_SlotArea)            = elemR(:,er_SlotDepth) * elemR(:,er_SlotWidth)
                    elemR(:,er_SlotVolume)          = elemR(:,er_SlotArea) * elemR(:,er_Length)
                    !% --- add slot volume to total volume (which was set to full volume)
                    elemR(:,er_Volume)              = elemR(:,er_Volume) + elemR(:,er_SlotVolume)
                end where

            case default
                !% should not reach this stage
                print*, 'In ', subroutine_name
                print *, 'CODE ERROR Slot Method type unknown for # ', SlotMethod
                print *, 'which has key ',trim(reverseKey(SlotMethod))
                stop 38756
            end select
        end if

        if (setting%Debug%File%initial_condition) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_IC_slot
!
!==========================================================================
!==========================================================================
!
    subroutine init_reference_head ()
        !%------------------------------------------------------------------
        !% Description:
        !% computes the reference head 
        !%------------------------------------------------------------------
        !% Declarations:
            integer, pointer :: Npack, thisP(:)
            integer :: er_set(5), esr_set(8), fr_set(3)
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------   
        !% Aliases
            !% --- use only the time-marching elements to set reference head
            Npack => npack_elemP(ep_ALLtm)
            thisP => elemP(1:Npack,ep_ALLtm)
        !%------------------------------------------------------------------  

        !% --- Get the reference head  
        if (Npack > zeroI) then      
            if (setting%Solver%SubtractReferenceHead) then
                setting%Solver%ReferenceHead  = minval(elemR(thisP,er_Zbottom))
                setting%Solver%ReferenceHead  &
                    = real(floor(setting%Solver%ReferenceHead),8) 
            else
                setting%Solver%ReferenceHead = zeroR
            end if
            setting%Solver%MaxZbottom     = maxval(elemR(thisP,er_Zbottom))
            setting%Solver%MinZbottom     = minval(elemR(thisP,er_Zbottom))
            setting%Solver%AverageZbottom =    sum(elemR(thisP,er_Zbottom)) / Npack
        end if
        !% set the min and max over all the processors.
        call co_min(setting%Solver%ReferenceHead)
        call co_max(setting%Solver%MaxZbottom)
        call co_min(setting%Solver%MinZbottom)
        call co_max(setting%Solver%AverageZbottom)

        !% --- bug checking
        ! if (this_image() == 1) then
        !     write(*,*) 'reference head   ',setting%Solver%ReferenceHead
        !     write(*,*) 'max z bottom     ',setting%Solver%MaxZbottom
        !     write(*,*) 'average z bottom ',setting%Solver%AverageZbottom
        !     write(*,*) 'min z bottom     ',setting%Solver%MinZbottom
        ! end if
  
        !%------------------------------------------------------------------    
        !% Closing 
    end subroutine init_reference_head
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_subtract_reference_head ()
        !%------------------------------------------------------------------
        !% Description:
        !% removes reference head from Z values in all arrays
        !% except for BC, which is done in bc_fetch()
        !%------------------------------------------------------------------
        !% Declarations:
            integer, pointer :: Npack, thisP(:)
            integer :: er_set(5), esr_set(8), fr_set(3)
        !%------------------------------------------------------------------
        !% Preliminaries:
        !%------------------------------------------------------------------   
        !% Aliases
            !% --- use only the time-marching elements to set reference head
            Npack => npack_elemP(ep_ALLtm)
            thisP => elemP(1:Npack,ep_ALLtm)
        !%------------------------------------------------------------------   
        !% list of indexes using Z reference
                ! er_Head
                ! er_Head_N0
                ! er_Zbottom
                ! er_ZbreadthMax
                ! er_Zcrown
                ! esr_Weir_NominalDownstreamHead
                ! esr_Weir_Zcrown
                ! esr_Weir_Zcrest
                ! esr_Orifice_NominalDownstreamHead
                ! esr_Orifice_Zcrown
                ! esr_Orifce_Zcrest
                ! esr_Outlet_NominalDownstreamHead
                ! esr_Outlet_Zcrest
                ! fr_Head_u
                ! fr_Head_d
                ! fr_Zbottom

        !% --- subtract the reference head from elemR
        er_set = (/er_Head,        &
                  er_Head_N0,     &
                  er_Zbottom,      &
                  er_ZbreadthMax, &
                  er_Zcrown/)
        where (elemR(:,er_set) .ne. nullValueR)        
            elemR(:,er_set) = elemR(:,er_set) - setting%Solver%ReferenceHead
        endwhere

        ! !% --- subtract the reference head from elemSR
        ! esr_set = (/esr_Weir_NominalDownstreamHead,    &
        !            esr_Weir_Zcrown,                   &
        !            esr_Weir_Zcrest,                   &
        !            esr_Orifice_NominalDownstreamHead, &
        !            esr_Orifice_Zcrown,                &
        !            esr_Orifice_Zcrest,                &
        !            esr_Outlet_NominalDownstreamHead,  &
        !            esr_Outlet_Zcrest/)

        ! where (elemSR(:,esr_set) .ne. nullValueR)        
        !        elemSR(:,esr_set) = elemSR(:,esr_set) - setting%Solver%ReferenceHead
        ! endwhere
        
        !% substract the reference head from weir elemSR
        where (elemI(:,ei_elementType) == weir)      
               elemSR(:,esr_Weir_NominalDownstreamHead) = elemSR(:,esr_Weir_NominalDownstreamHead) &
                                                - setting%Solver%ReferenceHead
               elemSR(:,esr_Weir_Zcrown) = elemSR(:,esr_Weir_Zcrown) - setting%Solver%ReferenceHead
               elemSR(:,esr_Weir_Zcrest) = elemSR(:,esr_Weir_Zcrest) - setting%Solver%ReferenceHead
        endwhere

        !% substract the reference head from orifice elemSR
        where (elemI(:,ei_elementType) == orifice)      
               elemSR(:,esr_Orifice_NominalDownstreamHead) = elemSR(:,esr_Orifice_NominalDownstreamHead) &
                                                - setting%Solver%ReferenceHead
               elemSR(:,esr_Orifice_Zcrown) = elemSR(:,esr_Orifice_Zcrown) - setting%Solver%ReferenceHead
               elemSR(:,esr_Orifice_Zcrest) = elemSR(:,esr_Orifice_Zcrest) - setting%Solver%ReferenceHead
        endwhere

        !% substract the reference head from outlet elemSR
        where (elemI(:,ei_elementType) == outlet)      
               elemSR(:,esr_Outlet_NominalDownstreamHead) = elemSR(:,esr_Outlet_NominalDownstreamHead) &
                                                - setting%Solver%ReferenceHead
               elemSR(:,esr_Outlet_Zcrest) = elemSR(:,esr_Outlet_Zcrest) - setting%Solver%ReferenceHead
        endwhere

        !% --- subtract the refence head from faceR
        fr_set = (/fr_Head_u, &
                  fr_Head_d, &
                  fr_Zbottom/)
        where (faceR(:,fr_set) .ne. nullValueR)        
               faceR(:,fr_set) = faceR(:,fr_set) - setting%Solver%ReferenceHead
        endwhere          

        !%------------------------------------------------------------------    
        !% Closing 
    end subroutine init_subtract_reference_head
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_bc()
        !%-----------------------------------------------------------------------------
        !%
        !% Description:
        !%    Initializes boundary connditions
        !%
        !% Notes:
        !%    The structures are general enough to support 3 types of BCs:
        !%
        !%    BCup: updstream boundary condition which can be inflow or head BC
        !%    BCdn: downstream boundary condition which can be inflow or head BC
        !%    BClat: lateral inflow coming into and nJ2 or nJm node.
        !%
        !%    However, the code only supports inflow BCs for BCup and BClat,
        !%    and head BCs for BCdn, mimimcking EPA-SWMM 5.13 functionalities.
        !%    Further developments allowing other types of inflow and head BCs,
        !%    should store the respective BC in either the BC%inflowX or the
        !%    BC%headX arrays defining the corresponding type of BC (i.e., BCup,
        !%    BCdn, and BClat) in the BC%xI(:,bi_category) column.
        !%
        !%-----------------------------------------------------------------------------
        integer :: ii, nidx, ntype, counter_bc_er, outfallType
        integer :: SWMMtseriesIdx, SWMMbasepatType
        character(64) :: subroutine_name = "init_bc"
        !%-----------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%initialization)  &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        
        if (setting%Profile%useYN) call util_profiler_start (pfc_init_bc)

        !% --- set the key values to undefinedKey
        call util_key_default_bc()

        !% --- Set all to null with fetch to 1 and upper index to 0
        if (N_flowBC > 0) then
            BC%flowI = nullvalueI
            BC%flowR = nullvalueR
            BC%flowTimeseries = nullValueR
            BC%flowR(:, br_timeInterval) = abs(nullvalueR)  !% ensure positive
            BC%flowI(:,bi_fetch) = 1
            BC%flowI(:,bi_TS_upper_idx) = 0  !% latest position of upper bound in flow table
        end if
        if (N_headBC > 0) then
            BC%headI = nullvalueI
            BC%headTimeseries = nullValueR
            BC%headR(:, br_timeInterval) = abs(nullvalueR)  !% ensure positive
            BC%headI(:,bi_fetch) = 1
            BC%headI(:,bi_TS_upper_idx) = 0
        end if

        !% --- Initialize Inflow BCs
        if (N_flowBC > 0) then
            do ii = 1, N_flowBC
                nidx  = node%P%have_flowBC(ii)
                ntype = node%I(nidx, ni_node_type)

                ! print *, ' '
                ! print *, 'in ',trim(subroutine_name)
                ! print *, ii, nidx
                ! print *, ntype, reverseKey(ntype)
                ! print *, 'ext inflow ',node%YN(nidx, nYN_has_extInflow)
                ! print *, 'dwf Inflow ', node%YN(nidx, nYN_has_dwfInflow)

                !% Handle Inflow BCs (BCup and BClat only)
                if (node%YN(nidx, nYN_has_extInflow) .or. node%YN(nidx, nYN_has_dwfInflow)) then

                    BC%flowI(ii, bi_node_idx) = nidx
                    BC%flowI(ii, bi_idx)      = ii
                    BC%flowYN(ii,bYN_read_input_series) = .true.
                    BC%flowI(ii, bi_face_idx) = node%I(nidx,ni_face_idx)
                    BC%flowI(ii, bi_elem_idx) = node%I(nidx,ni_elem_idx)

                    ! print *, ''
                    ! print *, 'in ',trim(subroutine_name)
                    ! print *, 'node idx  ',nidx
                    ! print *, 'node name ',trim(node%Names(nidx)%str)
                    ! print *, 'node type ',trim(reverseKey(node%I(nidx,ni_node_type)))

                    !% --- assign category and face index
                    select case (ntype)
                    case (nJm)
                        !% --- standard junction
                        BC%flowI(ii, bi_category) = BClat
                        BC%flowI(ii, bi_face_idx) = nullvalueI
                        BC%flowI(ii, bi_elem_idx) = node%I(nidx,ni_elem_idx)
                    case (nJ1)
                        !% --- dead end without BCup
                        BC%flowI(ii, bi_category) = BClat
                        print *, 'CODE NEEDS TESTING: BClat inflow for dead-end nJ1 node has not been tested'
                        call util_crashpoint(5586688)
                    case (nJ2) 
                        !% --- face node (no storage) with lateral inflow into adjacent element
                        BC%flowI(ii, bi_category) = BClat
                        !BC%flowI(ii, bi_elem_idx) = node%I(nidx, ni_elemface_idx) !% elem idx OBSOLETE
                        ! print *, 'CODE NEEDS TESTING: BClat inflow for nJ2 node has not been tested'
                        ! print *, 'BC flow index ',ii
                        ! call util_crashpoint(7783723)
                    case (nBCup)
                        BC%flowI(ii, bi_category) = BCup
                        !BC%flowI(ii, bi_face_idx) = node%I(nidx, ni_elemface_idx) !% face idx OBSOLETE
                    case default
                        print *, "Error, BC type can't be an inflow BC for node " // node%Names(nidx)%str
                        !stop 
                        call util_crashpoint(739845)
                        !return
                    end select

                    !% HACK Pattern needs checking --- the following may be wrong! brh20211221
                    !% check whether there is a pattern (-1 is no pattern) for this inflow
                    SWMMbasepatType = &
                        interface_get_nodef_attribute(nidx, api_nodef_extInflow_basePat_type)
                    !% brh20211216 should be _type?    
                    !rm nbasepat = &
                    !rm    interface_get_nodef_attribute(nidx, api_nodef_extInflow_basePat)
                    
                    !% check whether there is a time series 
                    !% (-1 is none, >0 is index, API_NULL_VALUE_I is error, which crashes API)
                    SWMMtseriesIdx = &
                        interface_get_nodef_attribute(nidx, api_nodef_extInflow_tSeries)

                    !% BC does not have fixed value if its associated with dwfInflow
                    !% or if extInflow has tseries or pattern
                    BC%flowI(ii, bi_subcategory) = BCQ_tseries
                    if (.not. node%YN(nidx, nYN_has_dwfInflow)) then !% extInflow only
                        !% brh 20211216 modified the following for basepatType == -1 rather than basepat /= -1
                        if ((SWMMtseriesIdx == -1) .and. (SWMMbasepatType == -1)) then
                            BC%flowI(ii, bi_subcategory) = BCQ_fixed
                        end if
                    end if

                else
                    print *, "CODE ERROR: unexpected else."
                    print *, "Only nodes with extInflow or dwfInflow can have inflow BC"
                    call util_crashpoint(826549)

                end if
            end do
        end if

        !     print *, ' '
        !     print *, ' in ',trim(subroutine_name)
        !     do ii=1,N_flowBC
        !         print *, ii, node%P%have_flowBC(ii)
        !         print *, 'node type ', node%I(nidx, ni_node_type), trim(reverseKey(node%I(nidx, ni_node_type)))
        !         print *, BC%flowI(ii,bi_idx),BC%flowI(ii,bi_node_idx)
        !         print *, 'elem, face : ',BC%flowI(ii,bi_elem_idx), BC%flowI(ii,bi_face_idx)
        !         print *, 'face up of elem ', elemI(BC%flowI(ii,bi_elem_idx),ei_Mface_uL)
        !         print *, 'elem dn of face ', faceI(BC%flowI(ii,bi_face_idx),fi_Melem_dL)
        !     end do
        !     print *, ' '
        !     print *,'dummy idx is ',dummyIdx
        !     do ii=1,N_elem(1)
        !         print *, ' '
        !         print *, faceI(elemI(ii,ei_Mface_uL),fi_Melem_uL)
        !         print *, elemI(ii,ei_Mface_uL), ii, elemI(ii,ei_Mface_dL)
        !         print *, faceI(elemI(ii,ei_Mface_dL),fi_Melem_dL)
        !     end do


        !% --- Initialize Head BCs
        if (N_headBC > 0) then
            do ii = 1, N_headBC
                nidx =  node%P%have_headBC(ii)
                ntype = node%I(nidx, ni_node_type)

                BC%headI(ii, bi_idx) = ii
                BC%headI(ii, bi_node_idx) = nidx
                BC%headI(ii, bi_face_idx) = node%I(nidx, ni_face_idx) 
                BC%headI(ii, bi_elem_idx) = node%I(nidx, ni_elem_idx)

                select case (ntype)
                case (nBCdn)
                    BC%headI(ii, bi_category) = BCdn
                case default
                    print *, "CONFIGURATION OR CODE ERROR: a head boundary condition is "
                    print *, "designated on something other than an nBCdn node, which is not allowed"
                    print *, "node index is ",nidx
                    print *, "node name is  ", trim(node%Names(nidx)%str) 
                    if (ntype < (keys_lastplusone-1)) then
                        print *, "node type is ",reverseKey(ntype)
                    else
                        print *, "node type # is invalid: ",ntype
                    end if
                    call util_crashpoint(57635)
                    !return
                end select

                !% --- get the outfall type
                outfallType = int(interface_get_nodef_attribute(nidx, api_nodef_outfall_type))
                select case (outfallType)
                case (API_FREE_OUTFALL)
                    !% debug test 20220725brh
                    BC%headI(ii, bi_subcategory) = BCH_free
                    BC%headYN(ii, bYN_read_input_series) = .false.

                case (API_NORMAL_OUTFALL)
                    !% debug tested 20220729brh
                    BC%headI(ii, bi_subcategory) = BCH_normal
                    BC%headYN(ii, bYN_read_input_series) = .false.

                case (API_FIXED_OUTFALL) 
                    !% debug tested 20220729brh
                    BC%headI(ii, bi_subcategory) = BCH_fixed
                    BC%headYN(ii, bYN_read_input_series) = .false.

                case (API_TIDAL_OUTFALL)
                    !% debug tested 20220729brh
                    BC%headI(ii, bi_subcategory) = BCH_tidal
                    BC%headYN(ii, bYN_read_input_series) = .true.

                case (API_TIMESERIES_OUTFALL)
                    !% debug tested 2020729brh
                    BC%headI(ii, bi_subcategory) = BCH_tseries
                    BC%headYN(ii, bYN_read_input_series) = .true.

                case default
                    print *, 'CODE ERROR: unexpected case default'
                    call util_crashpoint(33875)
                end select

                !% --- check for a flap gate
                if (interface_get_nodef_attribute(nidx, api_nodef_hasFlapGate) == oneR) then
                    BC%headYN(ii,bYN_hasFlapGate) = .true.
                else
                    BC%headYN(ii,bYN_hasFlapGate) = .false.
                endif

            end do
        end if

        ! print *, ' '
        ! print *, ' in ',trim(subroutine_name)
        ! do ii=1,N_headBC
        !     print *, ii, node%P%have_headBC(ii)
        !     print *, 'node type ', node%I(nidx, ni_node_type), trim(reverseKey(node%I(nidx, ni_node_type)))
        !     print *, BC%headI(ii,bi_idx),BC%headI(ii,bi_node_idx)
        !     print *, 'elem, face : ',BC%headI(ii,bi_elem_idx), BC%headI(ii,bi_face_idx)
        !     print *, 'face dn of elem ', elemI(BC%headI(ii,bi_elem_idx),ei_Mface_dL)
        !     print *, 'elem up of face ', faceI(BC%headI(ii,bi_face_idx),fi_Melem_uL)
        ! end do
        ! print *, ' '
        ! print *,'dummy idx is ',dummyIdx
        ! do ii=1,N_elem(1)
        !     print *, ' '
        !     print *, faceI(elemI(ii,ei_Mface_uL),fi_Melem_uL)
        !     print *, elemI(ii,ei_Mface_uL), ii, elemI(ii,ei_Mface_dL)
        !     print *, faceI(elemI(ii,ei_Mface_dL),fi_Melem_dL)
        ! end do

    
        call bc_step()
        if (crashI==1) return

        call pack_data_BC()

        if (setting%Profile%useYN) call util_profiler_stop (pfc_init_bc)

        if (setting%Debug%File%initialization)  &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine init_bc
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_uniformtable_array ()
        !%------------------------------------------------------------------
        !% Description:
        !% initializes sectionfactor arrays (depth = f(sectionFactor))
        !% for computing normal depth
        !% 20220726 -- only stored for elements upstream of an outfall
        !%------------------------------------------------------------------
            integer       :: ii,  lastUT_idx  , kk    
            real(8) :: thisdepth, thiswidth, Qdelta
            character(64) :: subroutine_name = 'init_uniformtable_array'
        !%------------------------------------------------------------------

        ! print *, 'in ',trim(subroutine_name)
        
        call util_allocate_uniformtable_array()

        lastUT_idx = 0  !% last used index to uniform table

        !% --- set up uniform tables for section factor and critical flow for head BC locations
        call init_BChead_uniformtable (lastUT_idx)

        ! print *, 'Qcrit max ',uniformTableR(1,utr_QcritMax) , uniformTableR(1,utr_AreaMax) * sqrt(setting%Constant%gravity * uniformTableR(1,utr_DepthMax) )
        ! print *, 'Area max  ',uniformTableR(1,utr_AreaMax)
        ! print *, 'DepthMax  ',uniformTableR(1,utr_DepthMax)

        ! stop 50987

        !% THIS IS WHERE WE WOULD INSERT ANY OTHER UNIFORM TABLE INITIATIONS
        !% NEW DATA STARTs FROM lastUT_idx+1

        !% --- fill of values for each location
        do ii = 1,size(uniformTableDataR,1)

               !print *, 'TEMP COMMENT OUT OF SECTION FACTOR'

            !% --- uniformly-distributed section factor
            call init_uniformtabledata_Uvalue(ii,utr_SFmax, utd_SF_uniform)

            !% -- uniformly-distributed critical flow
            call init_uniformtabledata_Uvalue(ii,utr_QcritMax, utd_Qcrit_uniform)
   
            !% --- nonuniform values mapping from section factors
               ! print *, 'sectionfactors by depth ----------------'
            call init_uniformtabledata_nonUvalue (ii, utd_SF_depth_nonuniform, utd_SF_uniform)
               ! print *, 'sectionfactors by area ----------------'
            call init_uniformtabledata_nonUvalue (ii, utd_SF_area_nonuniform,  utd_SF_uniform)
   
            !% --- nonuniform values mapping from critical flow
               ! print *, 'Qcritical by depth ------------------'
            call init_uniformtabledata_nonUvalue (ii, utd_Qcrit_depth_nonuniform, utd_Qcrit_uniform)
               ! print *, 'Qcritical by area ------------------'
            call init_uniformtabledata_nonUvalue (ii, utd_Qcrit_area_nonuniform,  utd_Qcrit_uniform)

            ! print *, ' '
            ! print *, 'depths for uniformly-distributed section factors'
            ! do kk=1,size(uniformTableDataR,2)
            !     print *, ii, uniformTableDataR(ii,kk,utd_SF_depth_nonuniform), &
            !     uniformTableDataR(ii,kk,utd_SF_depth_nonuniform) * uniformTableR(ii,utr_DepthMax)
            ! end do

            ! print *, ' '
            ! print *, 'implied sF from depth and actual sf for rectangular for simple case'
            ! thiswidth = elemR(1,er_BreadthMax)
            ! do kk=1,size(uniformTableDataR,2)
            !     thisdepth = uniformTableDataR(ii,kk,utd_SF_depth_nonuniform) * uniformTableR(ii,utr_DepthMax)
            !     print *, ii, ((thiswidth*thisdepth)**(5.d0/3.d0)) &
            !      / ((thiswidth + twoR*thisdepth)**(2.d0/3.d0)), &
            !      uniformtableDataR(ii,kk,utd_SF_uniform)* uniformTableR(ii,utr_SFmax)
            ! end do

            !print *, ' '
           ! print *, 'implied Qcrit from depth and Qcrit from table for rectangular simple case'
            !thiswidth = elemR(1,er_BreadthMax)
            ! do kk=1,size(uniformTableDataR,2)
            !     thisdepth = uniformTableDataR(ii,kk,utd_Qcrit_depth_nonuniform) * uniformTableR(ii,utr_DepthMax)
            !     print *, ii, (thisdepth * thiswidth) * sqrt(9.81d0 * thisdepth), &
            !     uniformtableDataR(ii,kk,utd_Qcrit_uniform)* uniformTableR(ii,utr_QcritMax)
            ! end do

            ! print *, ' '
            ! Qdelta = uniformTableR(ii,utr_QcritMax) / real(size(uniformTableDataR,2)-1,8)
            ! print *, 'Qmax, Qdelta, size'
            ! print *, uniformTableR(ii,utr_QcritMax), Qdelta, size(uniformTableDataR,2)-1,8
            ! do kk=1,size(uniformTableDataR,2)
            !     print *, kk, uniformTableDataR(ii,kk,utd_Qcrit_Depth_nonuniform)
            ! end do

            !stop 34987
        end do

        ! print *, ' '
        ! print *, 'Qmax ',uniformTableR(1,utr_QcritMax)
        ! do kk=1,size(uniformTableDataR, 2)
        !     !print *, kk, uniformTableDataR(1,kk,utd_Qcrit_depth_nonuniform) * uniformTableR(1,utr_DepthMax), &
        !     !             uniformTableDataR(1,kk,utd_Qcrit_Area_nonuniform)  * uniformTableR(1,utr_AreaMax)
        !     print *, uniformTableDataR(1,kk,utd_Qcrit_Area_nonuniform)  * uniformTableR(1,utr_AreaMax) &
        !                 * sqrt(setting%Constant%gravity &
        !                        * uniformTableDataR(1,kk,utd_Qcrit_depth_nonuniform) * uniformTableR(1,utr_DepthMax) ) &
        !                 , (uniformTableR(1,utr_QcritMax)  / real((N_uniformTableData_items-1),8)) * real(kk-1,8)
        ! end do

        !  stop 5559783

    end subroutine init_uniformtable_array    
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_BChead_uniformtable (UT_idx)
        !%------------------------------------------------------------------ 
        !% Description
        !% Get the maximum values (no necessarily at full!) that are used
        !% for normalizing the uniform tables that are lookup by SF or Depth
        !% UT_idx is the last uniform table index used, which is incremented
        !% as more table data is added
        !% NOTE: see init_uniformtable_array for the actual table values
        !%------------------------------------------------------------------ 
        !% Declarations
            integer, intent (inout) :: UT_idx
            integer, pointer        :: eIdx
            integer                 :: ii, jj
            real(8), pointer        :: grav
            real(8)                 :: sf, qcrit, thisDepth, deltaD, depthTol
            character(64)           :: subroutine_name = 'init_BChead_uniformtable'
        !%------------------------------------------------------------------ 
        !% Aliases
            grav         => setting%Constant%gravity
        !%------------------------------------------------------------------ 
        !% --- return if there are no head BC
        if (N_headBC < 1) return

        do ii = 1,N_headBC

            !% --- increment over the last-used uniform table index
            UT_idx = UT_idx + 1

            !% --- the element index for the element upstream of the BC
            eIdx => BC%headI(ii, bi_elem_idx)

            !% --- store indexes
            uniformTableI(UT_idx,uti_idx)        = UT_idx  !% self store
            uniformTableI(UT_idx,uti_elem_idx)   = eIdx    !% element lcoation
            uniformTableI(UT_idx,uti_BChead_idx) = ii      !% BC head index
            BC%headI     (ii    ,bi_UTidx)       = UT_idx  !% ensure BC head knows the UT index

            !% --- store the maximum depths and areas for the location
            uniformTableR(UT_idx,utr_DepthMax) =  elemR(eIdx,er_FullDepth)
            uniformTableR(UT_idx,utr_AreaMax)  =  elemR(eIdx,er_FullArea)

            ! print *, 'FULL AREA ',uniformTableR(UT_idx,utr_AreaMax)
            ! stop 3409872

            !% --- maximum Qcrit flow is where Fr = 1 or Q = A sqrt(gH)
            uniformTableR(UT_idx,utr_QcritMax)  &
                =   uniformTableR(UT_idx,utr_AreaMax) &
                    * sqrt(grav * uniformTableR(UT_idx,utr_DepthMax))

            !% --- get max value of SectionFactor by stepping through cross-section
            !%     this allows us to deal with slight non-monotonic behavior in nearly full conduits
            thisDepth = zeroR
            deltaD = uniformTableR(UT_idx,utr_DepthMax) / onethousandR
            uniformTableR(UT_idx,utr_SFmax)    = zeroR
            
            !% --- include a depth tolerance to prevent round-off from
            !%     creating a step larger than the max depth
            depthTol = deltaD / tenR
            jj=0
            !% --- cycle through all the depths to find the maximum section factor
            do while (thisdepth .le. (uniformTableR(UT_idx,utr_DepthMax)-depthTol))
                thisDepth = thisDepth + deltaD

                !% --- section factor at this depth
                sf = geo_sectionfactor_from_depth_singular (eIdx, thisDepth, setting%ZeroValue%Area, setting%ZeroValue%Depth)

                !% --- check if this is the max sf thus far
                uniformTableR(UT_idx,utr_SFmax)    = max(uniformTableR(UT_idx,utr_SFmax),sf)
            end do
        end do

    end subroutine init_BChead_uniformtable
!%
!%==========================================================================
!%==========================================================================
!%  
    subroutine init_uniformtabledata_nonUvalue ( &
        UT_idx,     &  ! index of the uniform table
        utd_nonU,   &  ! slice in uniformTableDataR where nonuniform data are stored
        utd_uniform &  ! slice in uniformTableDataR where corresponding uniform data are stored
        )    
        !%------------------------------------------------------------------ 
        !% Description
        !% initializes a non-uniform value in the uniformTableDataR array
        !%------------------------------------------------------------------ 
            integer, intent (in) :: UT_idx, utd_nonU, utd_uniform
            integer              :: Utype, NUtype, jj, utr_max
            integer, pointer     :: eIdx
            real(8), pointer     ::  grav
            real(8)  :: thisUvalue, deltaDepth, deltaUvalue, errorU
            real(8)  :: testUvalue, testDepth, testArea, testPerimeter
            !real(8)  :: interpUvalue
            real(8)  :: oldtestUvalue, oldtestDepth, oldtestArea, oldtestPerimeter
            real(8)  :: thisDepth, thisArea, thisPerimeter
            real(8), parameter :: uTol = 1.d-3
            logical :: isIncreasing
            character(64) :: subroutine_name = 'init_uniformtabledata_nonUvalue'
        !%------------------------------------------------------------------ 
        !% Aliases
            eIdx => uniformTableI(UT_idx,uti_elem_idx)  ! element index
            grav => setting%Constant%gravity
        !%------------------------------------------------------------------ 
        !% --- set the type for the nonuniform data
        !%     must be consistent with type of max data
        !%     must be consistent with a utd_... index,
        select case (utd_nonU)
            case (utd_SF_depth_nonuniform, utd_Qcrit_depth_nonuniform)
                NUtype = DepthData
            case (utd_SF_area_nonuniform, utd_Qcrit_area_nonuniform)
                NUtype = AreaData
            case default
                print *, 'CODE ERROR: unexpected case default'
                call util_crashpoint(6629873)
        end select
    
        !% set the type for the uniform data -- must be a utd_... index
        select case (utd_uniform)
            case (utd_SF_uniform)
                Utype = SectionFactorData
                utr_max = utr_SFmax
            case (utd_Qcrit_uniform)
                Utype = QcriticalData
                utr_max = utr_QcritMax
            case default
                print *, 'CODE ERROR: unexpected case default'
                call util_crashpoint(3609433)
        end select

        !% --- get the uniform data delta
        deltaUvalue = uniformTableR(UT_idx,utr_max) /  real((N_uniformTableData_items-1),8)

            print *, 'deltaUvalue ',deltaUvalue

        !print *,'Umax ',uniformTableR(UT_idx,utr_max)
        !print *, 'by depth ',geo_sectionfactor_from_depth_singular (eIdx,20.d0)
        !print *, 'UT_idx ',UT_idx
 
        
        !% --- Get delta step for stepping through the non-uniform computation
        !%     looking for at least 3 digits of precision in cycling through nonuniform
        !%     values
        !%     Note: We ALWAYS step through in depth
        deltaDepth = uniformTableR(UT_idx,utr_DepthMax) / real(1000*(N_uniformTableData_items-1),8)
        if (deltaDepth < setting%Eps%Machine) then
            print *, 'CONFIGURATION OR CODE ERROR: too small of a depth step in ',trim(subroutine_name)
        end if

            print *, 'deltaDepth ',deltaDepth

        testUvalue    = zeroR
        testDepth     = zeroR
        testArea      = zeroR
        testPerimeter = zeroR

        !% --- initialization: store all zeros for the first table items
        uniformTableDataR(UT_idx,1,utd_nonU) = zeroR

        !% --- retain zeros as the first table items, so start at column 2.
        do jj = 2, N_uniformTableData_items
            !% --- increment to the next value of the uniform data (unnormalize)
            thisUvalue = uniformTableDataR(UT_idx,jj,utd_uniform) * uniformTableR(UT_idx,utr_max)

            ! print *, ' '
            ! print *, 'jj, thisQcrit ',jj, thisUvalue

            !% --- iterate to find depth that provides uniform value just below and
            !%     just above the target (thisUvalue)
            isIncreasing = .true.
            do while ((testUvalue < thisUvalue) &
                     .and. (testDepth + deltaDepth .le. elemR(eIdx,er_FullDepth)) &
                     .and. isIncreasing)
                !% --- store the previous (low) guess
                oldtestUvalue    = testUvalue
                oldtestDepth     = testDepth
                oldtestArea      = testArea
                oldtestPerimeter = testPerimeter
                !% --- increment the test depth
                testDepth     = testDepth + deltaDepth
                testArea      = geo_area_from_depth_singular (eIdx, testDepth, zeroR)
                !% --- compute values for incremented depth
                select case (Utype)
                case (SectionFactorData)
                    testUvalue    = geo_sectionfactor_from_depth_singular (eIdx, testDepth, setting%ZeroValue%Area, setting%ZeroValue%Depth)
                case (QcriticalData)
                    testUvalue    = geo_Qcritical_from_depth_singular (eIdx, testDepth, zeroR)
                    !print *, 'test Qcrit ',testUvalue
                    !print *,  '   depth, depthfromArea ',testDepth, testArea / 1.5d0
                case default
                    print *, 'CODE ERROR: unexpected case default'
                    call util_crashpoint(608723)
                end select

                !print *, 'old,new Uvalue ',oldtestUvalue, testUvalue
                !% --- for monotonic, exit will be when testUvalue >= thisUvalue
                !% --- as soon as non-monotonic is found, the remainder of the
                !%     array uses the final depth value
                if (oldtestUvalue > testUvalue) isIncreasing = .false.
                ! print *, 'isIncreasing',isIncreasing
                ! print *, 'test value: ',testUvalue, thisUvalue
                ! print *, 'test depth: ',testDepth + deltaDepth, elemR(eIdx,er_FullDepth)
            end do

            !%--- get the best estimate of the value of the Depth at thisUvalue
            if (testUvalue .eq. thisUvalue) then
                thisDepth = testDepth
                thisArea  = testArea
            elseif (testUvalue < thisUvalue) then
                !% --- exited on depth exceeding max or non-monotonic, so use last values
                thisDepth  =  testDepth
                thisArea   =  testArea
            else
                !% --- interpolate across the two available values that bracket thisUvalue
                thisDepth  = oldtestDepth  +        deltaDepth            *  (thisUvalue - oldtestUvalue) / deltaUvalue
                thisArea   = oldtestArea   + (testArea  - oldtestArea)    *  (thisUvalue - oldtestUvalue) / deltaUvalue
                !print *, 'old, this, test Uvalue ',oldtestUvalue, thisUvalue, testUvalue
                !print *, 'ratio ',(thisUvalue - oldtestUvalue) / deltaUvalue
            endif

            !print *, 'thisDepth Qcrit ',thisDepth, thisArea * sqrt(setting%Constant%gravity * thisDepth)

            !% --- store the table data (normalized)   
            select case (NUtype)
            case (DepthData) 
                uniformTableDataR(UT_idx,jj,utd_nonU) = thisDepth / uniformTableR(UT_idx,utr_DepthMax)
            case (AreaData)
                uniformTableDataR(UT_idx,jj,utd_nonU) = thisArea  / uniformTableR(UT_idx,utr_AreaMax)
            case default
                print *, 'CODE ERROR: unexpected case default'
                call util_crashpoint(2398542)
            end select

            !% --- final check for this item
            select case (Utype)
            case (SectionFactorData)
                ! print*, '**************************************'
                ! print*, thisDepth, 'thisDepth'
                testUvalue    = geo_sectionfactor_from_depth_singular (eIdx, thisDepth, setting%ZeroValue%Area, setting%ZeroValue%Depth)
            case (QcriticalData)
                testUvalue    = geo_Qcritical_from_depth_singular (eIdx, thisDepth, ZeroR)
            case default
                print *, 'CODE ERROR: unexpected case default'
            end select
            !% --- relative error
            errorU = abs((thisUvalue - testUvalue) / uniformTableR(UT_idx,utr_max))

            ! print *, ' '
            ! print *, oldtestDepth, thisDepth, testDepth
            ! print *, oldtestUvalue, thisUvalue, testUvalue
            ! print *, 'max value ',uniformTableR(UT_idx,utr_max)
            ! print *, 'error ',errorU
            ! print *, ' '

            if (errorU > uTol) then
                print *, 'CODE ERROR in geometry processing for uniform table.'
                print *, 'tolerance setting is ',uTol
                print *, 'relative error is ',errorU
                call util_crashpoint(69873)
            end if

            ! if (jj > 50) then
            !     stop 298734
            ! end if
        end do
          

    end subroutine init_uniformtabledata_nonUvalue
!%
!%==========================================================================
!%==========================================================================
!%    
    subroutine init_uniformtabledata_Uvalue ( &
         UT_idx,    &  ! index of the uniform table
         utr_max,   &  ! column in uniformTableR where max uniform value is stored
         utd_uniform & ! slice in uniformTableDataR where uniform data are stored
        )
        !%------------------------------------------------------------------ 
        !% Description:
        !% computes and stores a normalized uniform data set in uniformTableDataR
        !% Note that if the minimum of the data is not equal to zero, the data
        !% is offset by the minimum so that the normalized uniform data always
        !% is from zero to one.
        !%------------------------------------------------------------------ 
        !% Declarations
            integer, intent(in) :: UT_idx, utr_max, utd_uniform
            real(8), pointer :: uniformMax
            real(8)          :: thisValue, normDelta
            integer          :: jj
            character(64)    :: subroutine_name = 'init_uniformtabledata_Uvalue'
        !%------------------------------------------------------------------ 

        ! print *, 'in ',trim(subroutine_name)
        ! print *, 'UT_idx ',UT_idx
        ! print *, 'utr_max ',utr_max

        !% --- maximum and mininum values of the uniform data
        uniformMax => uniformTableR(UT_idx,utr_max)

        ! print *, '----- uniformMax ',uniformMax

        !% --- step sizes in the uniform table
        normDelta = uniformMax / real(N_uniformTableData_items-1,8)

        ! print *, '----- normDelta ',normDelta

        !% --- store the zero as starting point for normalized table
        uniformTableDataR(UT_idx,1,utd_uniform) = zeroR
        thisValue = zeroR

        !% --- retain zeros as the first table items, so start at column 2.
        do jj = 2, N_uniformTableData_items
              !% --- increment to the next value of the uniform data
            thisValue = thisValue + normDelta
            !% --- store the table data (normalized)    
            uniformTableDataR(UT_idx,jj,utd_uniform) = thisValue / uniformMax     
            
        !    print *, '----- ',jj, thisValue,  uniformTableDataR(UT_idx,jj,utd_uniform)

        end do

       ! stop 29873
        
    end subroutine init_uniformtabledata_Uvalue
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_IC_bottom_slope ()
        !%------------------------------------------------------------------ 
        !% Description:
        !% computes the bottom slope of all channel and conduit elements
        !%------------------------------------------------------------------
            integer, pointer :: npack, thisP(:), fup(:), fdn(:), Fidx
            integer          :: thisCol, mm, ii, JBidx, Aidx, Ci
            real(8), pointer :: slope(:), length(:), fZbottom(:)
        !%------------------------------------------------------------------
        !% Aliases
            thisCol = ep_CC_ALLtm
            npack   => npack_elemP(thisCol)
            if (npack < 1) return
            thisP   => elemP(1:npack,thisCol)
            fup     => elemI(:,ei_Mface_uL)
            fdn     => elemI(:,ei_Mface_dL)
            slope   => elemR(:,er_BottomSlope)
            length  => elemR(:,er_Length)
            fZbottom => faceR(:,fr_Zbottom)
        !%------------------------------------------------------------------
        
        slope(thisP) =  (fZbottom(fup(thisP)) - fZbottom(fdn(thisP))) / length(thisP)

        !% --- check for slopes that are too small
        where (abs(slope(thisP)) < setting%ZeroValue%Slope)
            slope(thisP) = sign(setting%ZeroValue%Slope,slope(thisP))
        endwhere

        !% initialize bottom slope for JB -- BRUTE FORCE HACK
        do mm=1,N_elem(this_image())
            if (elemI(mm,ei_elementType) == JM) then
                !% -- upstream branches
                do ii=1,max_branch_per_node,2
                    JBidx = mm+ii
                    if (elemSI(JBidx,esi_JunctionBranch_Exists) == oneI) then 
                        Fidx => elemI(JBidx,ei_MFace_uL)
                        if (elemYN(JBidx,eYN_isBoundary_up)) then
                            Ci   = faceI(Fidx,fi_Connected_image)
                            Aidx = faceI(Fidx,fi_GhostElem_uL)
                        else
                            Ci   = this_image()
                            Aidx = faceI(Fidx,fi_Melem_uL)
                        end if
                        !% --- branch inherits slope of adjacent branch
                        elemR(JBidx,er_BottomSlope) = elemR(Aidx,er_BottomSlope)[Ci]  
                    else 
                        !% no action
                    end if
                end do
                !% --- downstream branches
                do ii=2,max_branch_per_node,2
                    JBidx = mm+ii
                    if (elemSI(JBidx,esi_JunctionBranch_Exists) == oneI) then 
                        Fidx => elemI(JBidx,ei_MFace_dL)
                        if (elemYN(JBidx,eYN_isBoundary_dn)) then
                            Ci   = faceI(Fidx,fi_Connected_image)
                            Aidx = faceI(Fidx,fi_GhostElem_dL)
                        else
                            Ci   = this_image()
                            Aidx = faceI(Fidx,fi_Melem_dL)
                        end if
                        !% --- branch inherits slope of adjacent branch
                        elemR(JBidx,er_BottomSlope) = elemR(Aidx,er_BottomSlope)[Ci]  
                    else 
                        !% no action
                    end if
                end do
            end if
        end do


        !%------------------------------------------------------------------
    end subroutine init_IC_bottom_slope    
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_IC_ZeroValues_nondepth ()
        !%------------------------------------------------------------------
        !% Description:
        !% ensures consistent initialization of zero values. 
        !% The ZeroValue%Depth must already be set
        !%------------------------------------------------------------------
        !% Declarations
            real(8), pointer :: area0, topwidth0, volume0, depth0, slope0, lengthNominal
            integer, pointer :: Npack, thisP, allP(:)
            integer, pointer :: elemPGx(:,:), npack_elemPGx(:), col_elemPGx(:)
            integer :: ii
            real(8) :: volumeIncrease, volume0a
        !%------------------------------------------------------------------
        !% Aliases
            area0     => setting%ZeroValue%Area
            topwidth0 => setting%ZeroValue%Topwidth
            volume0   => setting%ZeroValue%Volume
            depth0    => setting%ZeroValue%Depth  !% already set
            slope0    => setting%ZeroValue%Slope
            lengthNominal => setting%Discretization%NominalElemLength

            !% --- used for computing depth by type
            ! select case (whichTM)
            ! case (ETM)
                elemPGx                => elemPGetm(:,:)
                npack_elemPGx          => npack_elemPGetm(:)
                col_elemPGx            => col_elemPGetm(:)
            ! case default 
            !     print *, 'CODE ERROR, additional TM methods needed'
            ! end select
        !%------------------------------------------------------------------
        if (.not. setting%ZeroValue%UseZeroValues) return

        !% --- depth zero is used as set by json file
        if (depth0 .le. onethousandR * setting%Eps%Machine) then
            print *, 'USER ERROR: setting.ZeroValue.Depth is too small '
            print *, 'selected value is   ',depth0
            print *, 'minimum required is ', onethousandR * setting%Eps%Machine
            call util_crashpoint(798523)
            return
        end if

        !% --- slope zero is used as set by json file
        if (slope0 .le. onethousandR *setting%Eps%Machine) then
            print *, 'USER ERROR: setting.ZeroValue.Slope is too small '
            print *, 'selected value is   ',slope0
            print *, 'minimum required is ', onethousandR * setting%Eps%Machine
            call util_crashpoint(7985237)
            return
        end if

        !% --- cycle through to set ZeroValues consistent with depth
        !%     use the set of all time-marching elements
        Npack => npack_elemP(ep_ALLtm)
        if (Npack > 0) then
            !thisP => elemP(1:Npack,ep_ALLtm)

            !% --- temporary store of initial depth and replace with zero depth
            elemR(:,er_Temp04) = elemR(:,er_Depth)
            elemR(:,er_Depth)  = depth0 * 0.99d0

            do ii=1,Npack
                thisP => elemP(ii,ep_ALLtm)
                select case (elemI(thisP,ei_elementType))
                case (CC)
                    !% temporary store a values for zero depth
                    elemR(thisP,er_Temp01) = geo_topwidth_from_depth_singular (thisP, depth0, zeroR)
                    elemR(thisP,er_Temp02) = geo_area_from_depth_singular     (thisP, depth0, zeroR) 
                    !% volume is area * length
                    elemR(thisP,er_Temp03) = elemR(thisP,er_Temp02) * elemR(thisP,er_Length)
                case (JM)
                    !% topwidth and area are ignored for JM
                    elemR(thisP,er_Temp01) = abs(nullvalueR)
                    elemR(thisP,er_Temp02) = abs(nullvalueR)
                    elemR(thisP,er_Temp03) = storage_volume_from_depth_singular(thisP,depth0)
                case default
                    print *, 'CODE ERROR: unexpected case default'
                    print *, 'element type not handeled for type # ',elemI(thisP,ei_elementType)
                    print *, 'at element index ',thisP
                    print *, trim(reverseKey(elemI(thisP,ei_elementType)))
                    call util_crashpoint(6629873)
                end select


                ! print *, thisP, trim(reverseKey(elemI(thisP,ei_elementType))), &
                !  elemI(thisP,ei_geometryType), trim(reverseKey(elemI(thisP,ei_geometryType))), &
                !   elemR(thisP,er_Temp01)
                !  print *, ' ' 
                !  print *, thisP
                !  print *, trim(reverseKey(elemI(thisP,ei_elementType)))
                !  print *, 'depth0    = ',depth0
                !  print *, 'topwidth = ',elemR(thisP,er_Temp01)
                !  print *, 'area     = ',elemR(thisP,er_Temp02)
                !  print *, 'volume   = ',elemR(thisP,er_Temp03)

                            
            end do

            !% --- get the minimum values, use 1/2 to ensure
            !%     that a zerovalue for depth will have a larger
            !%     value of topwidth, area, and volume than the
            !%     zerovalues of the respective terms
            allP => elemP(1:Npack,ep_ALLtm)
            topwidth0 = minval( elemR(allP,er_Temp01)) * onehalfR
            area0     = minval( elemR(allP,er_Temp02)) * onehalfR
            volume0   = minval( elemR(allP,er_Temp03)) * onehalfR

            !% Ensure zero values are not too small
            if (topwidth0 .le. setting%Eps%Machine) then
                topwidth0 = onethousandR * setting%Eps%Machine
            end if

            if (area0 .le. setting%Eps%Machine) then
                area0 = onethousandR * setting%Eps%Machine
            end if

            if (volume0 .le. setting%Eps%Machine) then
                volume0 = onethousandR * setting%Eps%Machine
            end if

            !% --- checking scale consistency
            topwidth0 = max(topwidth0, area0 / depth0)
            volume0   = min(volume0, area0 * setting%Discretization%NominalElemLength)


            !% --- reset temporary arrays used above
            elemR(:,er_Temp01) = nullvalueR
            elemR(:,er_Temp02) = nullvalueR
            elemR(:,er_Temp03) = nullvalueR
            
            !% --- temporary store of original volume and
            !%     overwrite with volume0
            elemR(:,er_Temp03) = elemR(:,er_Volume)
            elemR(:,er_Volume) = volume0

            !% --- temporary store of original 
            elemR(:,er_Temp02) = elemR(:,er_Area)
        
            !% --- check the predicted depth0 from volume0--------------------------------
            !%     Goal is to ensure that D = f(V) returns D < D0 when V = V0
            !% --- store the base level volume0
            volume0a = volume0
            !% --- get the depth predicted from volume0 -- output stored in elemR(:,er_Depth)    
            call geo_depth_from_volume_by_type_allCC (elemPGetm, npack_elemPGetm, col_elemPGetm)

            !% --- cycle through elements to ensure that depth0 obtained from volume0
            !%     is smaller than the volume obtained by from depth0
            !%     If volume0 returns a depth larger than depth0, then reset volume0
            do ii = 1,N_elem(1)
                !print *, ii, elemR(ii,er_Depth),  depth0 - elemR(ii,er_Depth), elemR(ii,er_Volume)
                if (elemR(ii,er_Depth) > depth0) then
                    !% --- depth for volume0 is larger than depth0
                    volumeIncrease = (elemR(ii,er_Depth) - depth0) * elemR(ii,er_Topwidth) * elemR(ii,er_Length)
                    volume0 = min(volume0, min(volume0a - volumeIncrease, onehundredR*setting%Eps%Machine) )
                end if
            end do

            !% --- reset the depth from depth0 to IC value
            elemR(:,er_Depth) = elemR(:,er_Temp04)

            !% --- reset the volume from volume0 to IC value
            elemR(:,er_Volume) = elemR(:,er_Temp03)

            elemR(:,er_Temp03) = nullvalueR
            elemR(:,er_Temp04) = nullvalueR
            !  print *, ' '
            !  print *, 'depth0   ',depth0
            !  print *, 'topwidth0',topwidth0
            !  print *, 'area0    ',area0
            !  print *, 'volume0  ',volume0   

            ! !% --- compute the topwidths for zero depth
            ! !%     temporary store initial condition topwidth
            ! elemR(:,er_Temp01) = elemR(:,er_Topwidth)
            ! !% --- get the topwidth at zero depth using packed geometry arrays
            ! call geo_topwidth_from_depth (elemPGalltm, npack_elemPGalltm, col_elemPGalltm)
            ! !% --- use the minimum topwidth at zero depth as the smallest topwidth
            ! topwidth0 = minval(elemR(thisP,er_Topwidth))             
            ! !% --- return initial condition values to topwidth 
            ! elemR(:,er_Topwidth) = elemR(:,er_Temp01)
            ! !% --- return initial condition values to depth
            ! elemR(:,er_Depth)    = elemR(:,er_Temp02)
          
            ! !% OLD the zero topwidth is 5% of the max breadth        
            ! !OLD topwidth = minval(elemR(thisP,er_BreadthMax)) / twentyR

            ! !% the zerovalue area is 50% of the product of zerovalue depth and topwidth
            ! area0 = onehalfR * topwidth0 * depth0

            ! !% the zero value volume uses 5% of the volume at minimum depth
            ! volume0 = area0 * minval(elemR(thisP,er_Length)) / twentyR

            ! ! print *, topwidth, area, depth, volume, minval(elemR(thisP,er_Length))
        else
            print *, 'unexpected error -- no time-marching elements found '
            !stop 
            call util_crashpoint(398733)
            !return
        end if

        if (depth0 < setting%Eps%Machine) then
            print *, depth0
            print *, 'error, setting%ZeroValue%Depth is too small'
            !stop 
            call util_crashpoint(39870951)
            !return
        end if

        if (topwidth0 < setting%Eps%Machine) then
            print *, topwidth0
            print *, 'error, setting%ZeroValue%TopWidth is too small'
            !stop 
            call util_crashpoint(39870952)
            !return
        end if

        if (area0 < setting%Eps%Machine) then
            print *, area0
            print *, 'error, setting%ZeroValue%Area is too small'
            !stop 
            call util_crashpoint(93764)
            !return
        end if

        if (volume0 < setting%Eps%Machine) then
            print *, volume0
            print *, 'error, setting%ZeroValue%Volume is too small'
            !stop 
            call util_crashpoint(77395)
            !return
        end if

        !%------------------------------------------------------------------
        !% Closing
    end subroutine init_IC_ZeroValues_nondepth
!%
!%==========================================================================
!%==========================================================================
!%
    real(8) function init_IC_limited_fulldepth (thisDepth, thisLink) result(outDepth)
        !%------------------------------------------------------------------
        !% Description:
        !% Checks the input full depth and applies limiter (if needed)
        !% Note this should only be called for open-channel geometries
        !% This should NOT be applied to transect geometries.
        !%------------------------------------------------------------------  
        !% Declarations:
            integer, intent(in) :: thisLink
            real(8), intent(in) :: thisDepth
        !%------------------------------------------------------------------  

        if (setting%Link%OpenChannelLimitDepthYN) then 
            !% --- limit the output Depth
            outDepth = min(thisDepth, setting%Link%OpenChannelFullDepth)
        else
            if (thisDepth .eq. nullvalueR) then 
                print *, 'Unexpected initialization error: '
                print *, 'The maximum depth in link # ',thisLink
                print *, 'is set to the nullvalueR ',nullvalueR
                print *, 'Problem in SWMM link name ',trim(link%Names(thisLink)%str)
                call util_crashpoint(6698723)
            else
                outDepth = thisDepth
            end if
        end if

    end function init_IC_limited_fulldepth
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine init_IC_error_check ()
        !%------------------------------------------------------------------
        !% Description
        !% Configuration error checking
        !%------------------------------------------------------------------
        !% Declarations:
            integer            :: ii 
            integer, pointer   :: fUp, eUp, JMidx
        !%------------------------------------------------------------------
    
        do ii=1,N_elem(this_image())
            !% --- check that type 1 pump has upstream nJm that does NOT have implied storage
            if (elemI(ii,ei_elementType) == pump) then 
                ! print *, 'pump element ',ii
                ! print *, 'link # ',elemI(ii,ei_link_Gidx_BIPquick)
                ! print *, 'link name ',trim(link%Names(elemI(ii,ei_link_Gidx_BIPquick))%str)
                if (elemSI(ii,esi_Pump_SpecificType) == type1_Pump) then 
                    fUp => elemI(ii,ei_Mface_uL)
                    ! print *, 'face up ',fUp
                    eUp => faceI(fUp,fi_Melem_uL)
                    ! print *, 'elem up ',eUp
                    ! print *, 'elem up type # ',elemI(eUp,ei_elementType)
                    ! print *, 'elem up type name ',trim(reversekey(elemI(eUp,ei_elementType)))
                    if (elemI(eUp,ei_elementType) .ne. JB) then 
                        print *, 'CODE ERROR: upstream of a type1 pump should be JB'
                        call util_crashpoint(4429873)
                    else
                        ! print *, elemSI(eUp,esi_JunctionBranch_Main_Index)
                        JMidx => elemSI(eUp,esi_JunctionBranch_Main_Index)
                        if (elemSI(JMidx,esi_JunctionMain_Type) == NoStorage) then 
                            print *, 'USER CONFIGURATION ERROR'
                            print *, 'NoStorage found for Pump Type 1 node.'
                            print *, 'Change node to tabular storage or functional storage, or '
                            print *, 'use setting%Junction%ForceStorage == true to get implied storage'
                            print *, 'link number ',elemI(ii,ei_link_Gidx_SWMM)
                            print *, 'link name   ',trim(link%Names(elemI(ii,ei_link_Gidx_SWMM))%str)
                            call util_crashpoint(788734)
                        else
                            !% upstream element of pump has defined storage
                        end if
                    end if
                end if
            end if
        end do

    end subroutine  init_IC_error_check
!%
!%==========================================================================    
!% END MODULE
!%==========================================================================
!%
end module initial_condition
