module utility_allocate

    use define_globals
    use define_keys
    use define_indexes
    use define_settings, only: setting
    use define_globals
    use interface
    use utility

    ! use utility, only: utility_check_allocation

    implicit none

!-----------------------------------------------------------------------------
!
! Description:
!   This has all the large array allocation. The only other allocation
!   should be in boundary conditions (module bc)
!   All the top-level storage should be allocated in this module
!
!-----------------------------------------------------------------------------

    private

    ! utility_allocate constants
    integer           :: allocation_status
    character(len=99) ::              emsg

    ! public members
    public :: util_allocate_linknode
    public :: util_allocate_subcatch
    public :: util_allocate_partitioning_arrays
    public :: util_allocate_elemX_faceX
    public :: util_allocate_columns
    public :: util_allocate_bc
    public :: util_allocate_profiler
    public :: util_allocate_outputML_elemtypes
    public :: util_allocate_outputML_facetypes
    public :: util_allocate_outputML_storage
    public :: util_allocate_outputML_times
    public :: util_allocate_outputML_filenames
    public :: util_allocate_curves
    public :: util_allocate_curve_entries
    public :: util_allocate_check


contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine util_allocate_linknode()
        !%------------------------------------------------------------------
        !% Description:
        !%   Allocates the link and node storage used for the coarse representation
        !%   of the network connectivity
        !%
        !% Method:
        !%   The tables node%I, link%I, node%R, link%R, node%YN, link%YN, are allocated
        !%   These are defined in globals.f08). Every time memory is allocated, the
        !%   util_allocate_check functionality (from utility.f08) is used to
        !%   determine wheter or not there was an error during the allocation.
        !-------------------------------------------------------------------
        !% Declarations
            character(64) :: subroutine_name = 'util_allocate_linknode'
            integer       :: additional_rows = 0
            integer       :: ii, obj_name_len
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (icrash) return
            if (setting%Debug%File%utility_allocate) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
        !% If BIPquick is being used for Partitioning, include additional rows to the link-node arrays
        if (setting%Partitioning%PartitioningMethod == BQuick) then
            additional_rows = num_images() - 1
        end if

        allocate(node%I(N_node + additional_rows, Ncol_nodeI), stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'node%I')
        node%I(:,:) = nullvalueI

        allocate(link%I(SWMM_N_link + additional_rows, Ncol_linkI), stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'link%I')
        link%I(:,:) = nullvalueI

        allocate(node%R(N_node + additional_rows, Ncol_nodeR), stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'node%R')
        node%R(:,:) = nullvalueR

        allocate(link%R(SWMM_N_link + additional_rows, Ncol_linkR), stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'link%R')
        link%R(:,:) = nullvalueR

        allocate(node%YN(N_node + additional_rows, Ncol_nodeYN), stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'node%YN')
        node%YN(:,:) = nullvalueL

        allocate(link%YN(SWMM_N_link + additional_rows, Ncol_linkYN), stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'link%YN')
        link%YN(:,:) = nullvalueL

        !% |
        !% | Only names of objects present in EPA-SWMM are stored
        !% v

        !% --- allocate storage for node names
        allocate(node%Names(N_node), stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'node%Names')

        !% --- get the length of the node name and allocate node%Names(:)%str to the correct size
        do ii = 1, N_node
            obj_name_len = interface_get_obj_name_len(ii, API_NODE)
            allocate(character(obj_name_len) :: node%Names(ii)%str, stat=allocation_status, errmsg=emsg)
            call util_allocate_check(allocation_status, emsg, 'character(obj_name_len) :: node%Names(ii)%str')
        end do

        allocate(link%Names(N_link), stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'link%Names')

        do ii = 1, N_link
            obj_name_len = interface_get_obj_name_len(ii, API_LINK)
            allocate(character(obj_name_len) :: link%Names(ii)%str, stat=allocation_status, errmsg=emsg)
            call util_allocate_check(allocation_status, emsg, 'character(obj_name_len) :: link%Names(ii)%str')
        end do

        !% allocate link_node_output_idx
        allocate(node_output_idx(N_node + additional_rows),stat=allocation_status,errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'node_output_idx')

        allocate(link_output_idx(SWMM_N_link + additional_rows), stat=allocation_status,errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'link_output_idx')

        !%-------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%utility_allocate) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine util_allocate_linknode
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_allocate_subcatch()
        !%------------------------------------------------------------------
        !% Description:
        !%   Allocates the subcatchment storage 
        !%
        !% Method:
        !%   Every time memory is allocated, the
        !%   util_allocate_check functionality (from utility.f90) is used to
        !%   determine wheter or not there was an error during the allocation.
        !-------------------------------------------------------------------
        !% Declarations
            character(64) :: subroutine_name = 'util_allocate_subcatch'
            integer       :: ii, obj_name_len
        !%-------------------------------------------------------------------
        !% Preliminaries
            if (icrash) return
            if (setting%Debug%File%utility_allocate) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------

        !% subcatchR
        allocate(subcatchR(SWMM_N_subcatch, Ncol_subcatchR), stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'subcatchR')
        subcatchR(:,:) = nullvalueR

        !% subcatchI
        allocate(subcatchI(SWMM_N_subcatch, Ncol_subcatchI), stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'subcatchI')
        subcatchR(:,:) = nullvalueI

        !% subcatchYN
        allocate(subcatchYN(SWMM_N_subcatch, Ncol_subcatchYN), stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'subcatchYN')
        subcatchR(:,:) = nullvalueL

        !%-------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%utility_allocate) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine util_allocate_subcatch
!%
!%==========================================================================
!%==========================================================================
!%    
    subroutine util_allocate_partitioning_arrays()
        if (icrash) return
        allocate(adjacent_links(max_branch_per_node))
        allocate(elem_per_image(num_images()))
        allocate(image_full(num_images()))

        !% If BIPquick is being used for Partitioning, allocate additional arrays
        if (setting%Partitioning%PartitioningMethod == BQuick) then
            call util_count_node_types(N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2)

            allocate(B_nodeI(size(node%I,1), max_us_branch_per_node))
            allocate(B_nodeR(size(node%R,1), twoI))
            allocate(B_roots(N_nBCdn))
            allocate(totalweight_visited_nodes(size(node%I, oneI)))
            allocate(partitioned_nodes(size(node%I, oneI)))
            allocate(partitioned_links(size(link%I, oneI)))
            allocate(weight_range(size(link%I, oneI), twoI))
            allocate(accounted_for_links(size(link%I, oneI)))
            allocate(phantom_link_tracker(size(link%I, oneI)))
        end if
    end subroutine util_allocate_partitioning_arrays
!
!==========================================================================
!==========================================================================
!!% MOVED TO UTILITY DEALLOCATE
    ! subroutine util_deallocate_partitioning_arrays()

    !     if (allocated(adjacent_links)) deallocate(adjacent_links)
    !     if (allocated(elem_per_image)) deallocate(elem_per_image)
    !     if (allocated(image_full)) deallocate(image_full)

    !     !% If BIPquick is being used for Partitioning, allocate additional arrays
    !     if (setting%Partitioning%PartitioningMethod == BQuick) then
    !         deallocate(B_nodeI)
    !         deallocate(B_nodeR)
    !         deallocate(totalweight_visited_nodes)
    !         deallocate(partitioned_nodes)
    !         deallocate(partitioned_links)
    !         deallocate(weight_range)
    !         deallocate(accounted_for_links)
    !         deallocate(phantom_link_tracker)
    !     end if

    ! end subroutine util_deallocate_partitioning_arrays
!
!==========================================================================
!==========================================================================
!
    subroutine util_allocate_elemX_faceX ()
        ! the max_caf_elem and max_caf_face are the maximum length of the coarray
        ! across all employed images
        ! ==========================
        ! This will be excuted at parallel level
        ! ==========================
        integer :: ii
        integer, pointer :: ncol
        character(64) :: subroutine_name = 'util_allocate_elemX_faceX'

        !-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%utility_allocate) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !==== elem allocation ====
        ncol => Ncol_elemR ! the maxmiumu number of columns
        allocate(elemR(max_caf_elem_N+N_dummy_elem, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'elemR')
        elemR(:,:) = nullvalueR

        ncol => Ncol_elemI
        allocate(elemI(max_caf_elem_N+N_dummy_elem, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'elemI')
        elemI(:,:) = nullvalueI

        ncol => Ncol_elemYN
        allocate(elemYN(max_caf_elem_N+N_dummy_elem, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'elemYN')
        elemYN(:,:) = nullvalueL

        ncol => Ncol_elemP
        allocate(elemP(max_caf_elem_N+N_dummy_elem, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'elemP')
        elemP(:,:) = nullvalueI

        ncol => Ncol_elemPGalltm
        allocate(elemPGalltm(max_caf_elem_N+N_dummy_elem, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'elemPGalltm')
        elemPGalltm(:,:) = nullvalueI

        ncol => Ncol_elemPGetm
        allocate(elemPGetm(max_caf_elem_N+N_dummy_elem, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'elemPGetm')
        elemPGetm(:,:) = nullvalueI

        ncol => Ncol_elemPGac
        allocate(elemPGac(max_caf_elem_N+N_dummy_elem, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'elemPGac')
        elemPGac(:,:) = nullvalueI

        ncol => Ncol_elemSI
        allocate(elemSI(max_caf_elem_N+N_dummy_elem, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'elemSI')
        elemSI(:,:) = nullvalueI

        ncol => Ncol_elemSR
        allocate(elemSR(max_caf_elem_N+N_dummy_elem, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'elemSR')
        elemSR(:,:) = nullvalueR

        ncol => Ncol_elemSGR
        allocate(elemSGR(max_caf_elem_N+N_dummy_elem, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'elemSGR')
        elemSGR(:,:) = nullvalueR

        !==== face allocation ====
        ncol => Ncol_faceR
        allocate(faceR(max_caf_face_N, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'faceR')
        faceR(:,:) = nullvalueR

        ncol=> Ncol_faceI
        allocate(faceI(max_caf_face_N, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'faceI')
        faceI(:,:) = nullvalueI

        ncol=> Ncol_faceYN
        allocate(faceYN(max_caf_face_N, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'faceYN')
        faceYN(:,:) = nullvalueL

        ncol=> Ncol_faceP
        allocate(faceP(max_caf_face_N, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'Ncol_faceP')
        faceP(:,:) = nullvalueI

        allocate(facePS(max_caf_face_N, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'facePS')
        facePS(:,:) = nullvalueI

        ncol=> Ncol_faceM
        allocate(faceM(max_caf_face_N, ncol)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'faceM')
        faceM(:,:) = nullvalueL

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine util_allocate_elemX_faceX
!
!==========================================================================
!==========================================================================
!
    subroutine util_allocate_columns()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   All the enumerated variables can not be used as pointers. Thus the
        !   variables are stored in col_elemX(:) arrays that is a target
        !
        !-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'util_allocate_columns'
        !-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%utility_allocate) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% allocation of the col_elemX and npack_elemX
        call util_allocate_col_elemI
        call util_allocate_col_elemP
        call util_allocate_col_elemPGalltm
        call util_allocate_col_elemPGac
        call util_allocate_col_elemPGetm
        call util_allocate_col_elemR
        call util_allocate_col_elemSI
        call util_allocate_col_elemSR
        call util_allocate_col_elemSGR
        call util_allocate_col_elemWDI
        call util_allocate_col_elemWDR
        call util_allocate_col_elemYN
        call util_allocate_col_faceI
        call util_allocate_col_faceM
        call util_allocate_col_faceP
        call util_allocate_col_facePS
        call util_allocate_col_faceR
        call util_allocate_col_faceYN

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_allocate_columns
!
!==========================================================================
!==========================================================================
!
    subroutine util_allocate_bc()
        !-----------------------------------------------------------------------------
        ! allocate storage for boundary conditions.
        !-----------------------------------------------------------------------------
        character(64)      :: subroutine_name = 'util_allocate_bc'
        integer            :: ii, allocation_status, bc_node
        character(len=99)  :: emsg
        !-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%utility_allocate) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        if (setting%BC%TimeSlotsStored < 2) then
            print *, "Error: the number of slots has to be greater than 2"
            stop
        end if

        if (N_headBC > 0) then
            allocate(BC%headI(N_headBC, N_headI), stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg, 'BC%headI')
            BC%headI(:,:) = nullvalueI

            allocate(BC%headYN(N_headBC, N_headYN), stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg, 'BC%headYN')
            BC%headYN(:,:) = nullvalueL

            allocate(BC%headR_timeseries(N_headBC, setting%BC%TimeSlotsStored, N_headR), stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg, 'BC%headR_timeseries')
            BC%headR_timeseries(:,:,:) = nullvalueR

            allocate(BC%headIdx(N_headBC), stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg, 'BC%headIdx')
            BC%headIdx(:) = nullvalueI

            allocate(BC%headRI(N_headBC), stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg, 'BC%headRI')
            BC%headRI(:) = nullvalueR

        end if

        if (N_flowBC > 0) then
            allocate(BC%flowI(N_flowBC, N_flowI), stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg, 'BC%flowI')
            BC%flowI(:,:) = nullvalueI

            allocate(BC%flowYN(N_flowBC, N_flowYN), stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg, 'BC%flowYN')
            BC%flowYN(:,:) = nullvalueL

            allocate(BC%flowR_timeseries(N_flowBC, setting%BC%TimeSlotsStored, N_flowR), stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg, 'BC%flowR_timeseries')
            BC%flowR_timeseries(:,:,:) = nullvalueR

            allocate(BC%flowIdx(N_flowBC), stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg, 'BC%flowIdx')
            BC%flowIdx(:) = nullvalueI

            allocate(BC%flowRI(N_flowBC), stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg, 'BC%flowRI(N_flowBC)')
            BC%flowRI(:) = nullvalueR
            
        end if

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_allocate_bc
!
!==========================================================================
!==========================================================================
!
    subroutine util_allocate_profiler ()
        !-----------------------------------------------------------------------------
        character(64)      :: subroutine_name = 'util_allocate_profiler'
        integer            :: ii, allocation_status, bc_node
        character(len=99)  :: emsg
        !-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        if (setting%Profile%useYN) then
            !% allocate profiler data
            allocate(profiler_data(Nrow_pf,Ncol_pf), stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg, 'profiler_data')
            profiler_data(:,:) = zeroR

            !% allocate storage of profiled procedure name
            allocate(profiler_procedure_name(Ncol_pf), stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg, 'profiler_procedure_name')

            !% allocate storage of profiled procedure level (1 = upper, 2 = middle, 3 = lower)
            allocate(profiler_procedure_level(Ncol_pf), stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg, 'profiler_procedure_level')

            !% assign profile procedure name and level
            profiler_procedure_name(:) = 'name unassigned in code'
            profiler_procedure_level(:) = 0

            profiler_procedure_name(pfc_initialize_all) = 'initialize_all'
            profiler_procedure_level(pfc_initialize_all) = 1

            profiler_procedure_name(pfc_init_partitioning) = 'init_partitioning'
            profiler_procedure_level(pfc_init_partitioning) = 1

            profiler_procedure_name(pfc_init_network_define_toplevel) = 'init_network_define_toplevel'
            profiler_procedure_level(pfc_init_network_define_toplevel) = 1

            profiler_procedure_name(pfc_init_bc) = 'init_bc'
            profiler_procedure_level(pfc_init_bc) = 3

            profiler_procedure_name(pfc_init_IC_setup) = 'init_IC_setup'
            profiler_procedure_level(pfc_init_IC_setup) = 1

            profiler_procedure_name(pfc_init_IC_from_linkdata) = 'init_IC_from_linkdata'
            profiler_procedure_level(pfc_init_IC_from_linkdata) = 2

            profiler_procedure_name(pfc_init_IC_get_depth_from_linkdata) = 'init_IC_get_depth_from_linkdata'
            profiler_procedure_level(pfc_init_IC_get_depth_from_linkdata) = 3

            profiler_procedure_name(pfc_init_IC_get_flow_roughness_from_linkdata) = 'init_IC_get_flow_roughness_from_linkdata'
            profiler_procedure_level(pfc_init_IC_get_flow_roughness_from_linkdata) = 3

            profiler_procedure_name(pfc_init_IC_get_elemtype_from_linkdata) = 'init_IC_get_elemtype_from_linkdata'
            profiler_procedure_level(pfc_init_IC_get_elemtype_from_linkdata) = 3

            profiler_procedure_name(pfc_init_IC_get_geometry_from_linkdata) = 'init_IC_get_geometry_from_linkdata'
            profiler_procedure_level(pfc_init_IC_get_geometry_from_linkdata) = 2

            profiler_procedure_name(pfc_init_IC_get_channel_geometry) = 'init_IC_get_channel_geometry'
            profiler_procedure_level(pfc_init_IC_get_channel_geometry) = 3

            profiler_procedure_name(pfc_init_IC_get_conduit_geometry) = 'init_IC_get_conduit_geometry'
            profiler_procedure_level(pfc_init_IC_get_conduit_geometry) = 3

            profiler_procedure_name(pfc_init_IC_get_weir_geometry) = 'init_IC_get_weir_geometry'
            profiler_procedure_level(pfc_init_IC_get_weir_geometry) = 3

            profiler_procedure_name(pfc_init_IC_get_orifice_geometry) = 'init_IC_get_orifice_geometry'
            profiler_procedure_level(pfc_init_IC_get_orifice_geometry) = 3

            profiler_procedure_name(pfc_geo_assign_JB) = 'geo_assign_JB'
            profiler_procedure_level(pfc_geo_assign_JB) = 3

            profiler_procedure_name(pfc_init_IC_get_channel_conduit_velocity) = 'init_IC_get_channel_conduit_velocity'
            profiler_procedure_level(pfc_init_IC_get_channel_conduit_velocity) = 3

            profiler_procedure_name(pfc_init_IC_from_nodedata) = 'init_IC_from_nodedata'
            profiler_procedure_level(pfc_init_IC_from_nodedata) = 2

            profiler_procedure_name(pfc_init_IC_get_junction_data) = 'init_IC_get_junction_data'
            profiler_procedure_level(pfc_init_IC_get_junction_data) = 3

            profiler_procedure_name(pfc_update_auxiliary_variables) = 'update_auxiliary_variables'
            profiler_procedure_level(pfc_update_auxiliary_variables) = 2

            profiler_procedure_name(pfc_init_IC_set_SmallVolumes) = 'init_IC_set_SmallVolumes'
            profiler_procedure_level(pfc_init_IC_set_SmallVolumes) = 3

            profiler_procedure_name(pfc_init_IC_diagnostic_interpolation_weights) = 'init_IC_diagnostic_interpolation_weights'
            profiler_procedure_level(pfc_init_IC_diagnostic_interpolation_weights) = 3

            profiler_procedure_name(pfc_face_interpolation) = 'face_interpolation'
            profiler_procedure_level(pfc_face_interpolation) =  2

            profiler_procedure_name(pfc_diagnostic_toplevel) = 'diagnostic_toplevel'
            profiler_procedure_level(pfc_diagnostic_toplevel) = 2
        end if

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine util_allocate_profiler
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_allocate_outputML_elemtypes ()
        !%-----------------------------------------------------------------------------
        !% allocates the output type arrays for the element data that are output
        !%-----------------------------------------------------------------------------
        integer            :: allocation_status
        character(len=99)  :: emsg
        character(64)       :: subroutine_name = 'util_allocate_outputML_elemtypes'
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !% --- don't execute if no output elements
        if (.not. setting%Output%OutputElementsExist) return

        !% --- bug check
        if (N_OutTypeElem < 1) then
            write(*,"(A)") 'ERROR (code) the N_OutTypeElem is less than 1 (i.e. no output types selected)'
                write(*,"(A)") '... which should have caused Output.OutputElementsExist = .false.'
                write(*,"(A,i8)") '... setting%Output%N_OutTypeElem      ', N_OutTypeElem
                write(*,"(A,i8)") '... etting%Output%OutputElementsExist ', setting%Output%OutputElementsExist
                stop
        end if

        !% --- allocate the output types for elements
        allocate(output_types_elemR(N_OutTypeElem), stat=allocation_status, errmsg=emsg)
        call util_allocate_check (allocation_status, emsg, 'output_types_elemR')
        output_types_elemR(:) = nullvalueI

        !% --- allocate the output type processing for elements
        allocate(output_typeProcessing_elemR(N_OutTypeElem), stat=allocation_status, errmsg=emsg)
        call util_allocate_check (allocation_status, emsg, 'output_typeProcessing_elemR')
        output_typeProcessing_elemR(:) = nullvalueI


        !% --- allocate the output typeNames
        allocate(output_typeNames_elemR(N_OutTypeElem), stat=allocation_status, errmsg=emsg)
        call util_allocate_check (allocation_status, emsg, 'output_typeNames_elemR')
        output_typeNames_elemR(:) = ""

        !% --- allocate the output typeUnits
        allocate(output_typeUnits_elemR(N_OutTypeElem), stat=allocation_status, errmsg=emsg)
        call util_allocate_check (allocation_status, emsg, 'output_typeUnits_elemR')
        output_typeUnits_elemR(:) = ""

        !% --- allocate the output typeNames for  elements + time
        allocate(output_typeNames_withTime_elemR(N_OutTypeElem+1), stat=allocation_status, errmsg=emsg)
        call util_allocate_check (allocation_status, emsg, 'output_typeNames_withTime_elemR')
        output_typeNames_withTime_elemR(:) = ""

        !% --- allocate the output typeUnits for  elements + time
        allocate(output_typeUnits_withTime_elemR(N_OutTypeElem+1), stat=allocation_status, errmsg=emsg)
        call util_allocate_check (allocation_status, emsg, 'output_typeUnits_withTime_elemR')
        output_typeUnits_withTime_elemR(:) = ""

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine util_allocate_outputML_elemtypes
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_allocate_outputML_facetypes ()
        !%-----------------------------------------------------------------------------
        !% allocates the output type arrays for the face data that are output
        !%-----------------------------------------------------------------------------
        integer            :: allocation_status
        character(len=99)  :: emsg
        character(64)       :: subroutine_name = 'util_allocate_outputML_facetypes'
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% --- don't execute if no output faces
        if (.not. setting%Output%OutputFacesExist) return

        !% --- bug check
        if (N_OutTypeFace < 1) then
            write(*,"(A)") 'ERROR (code) the N_OutTypeFace is less than 1 (i.e. no output types selected)'
                write(*,"(A)") '... which should have caused Output.OutputElementsExist = .false.'
                write(*,"(A,i8)") '... setting%Output%N_OutTypeFace      ', N_OutTypeFace
                write(*,"(A,i8)") '... etting%Output%OutputElementsExist ', setting%Output%OutputFacesExist
                stop
        end if

        !% --- allocate the output types for faces
        allocate(output_types_faceR(N_OutTypeFace), stat=allocation_status, errmsg=emsg)
        call util_allocate_check (allocation_status, emsg, 'output_types_faceR')
        output_types_faceR(:) = nullvalueI

        !% --- allocate the output type processing for faces
        allocate(output_typeProcessing_faceR(N_OutTypeFace), stat=allocation_status, errmsg=emsg)
        call util_allocate_check (allocation_status, emsg, 'output_typeProcessing_faceR')
        output_typeProcessing_faceR(:) = nullvalueI

        !% --- allocate the output typeNames for faces
        allocate(output_typeNames_faceR(N_OutTypeFace), stat=allocation_status, errmsg=emsg)
        call util_allocate_check (allocation_status, emsg, 'output_typeNames_faceR')
        output_typeNames_faceR(:) = ""

        !% --- allocate the output typeUnits for faces
        allocate(output_typeUnits_faceR(N_OutTypeFace), stat=allocation_status, errmsg=emsg)
        call util_allocate_check (allocation_status, emsg, 'output_typeUnits_faceR')
        output_typeUnits_faceR(:) = ""

        !% --- allocate the output typenames for faces + time
        allocate(output_typeNames_withTime_faceR(N_OutTypeFace+1), stat=allocation_status, errmsg=emsg)
        call util_allocate_check (allocation_status, emsg, 'output_typeNames_withTime_faceR')
        output_typeNames_withTime_faceR(:) = ""

        !% --- allocate the output typeUnits for faces + time
        allocate(output_typeUnits_withTime_faceR(N_OutTypeFace+1), stat=allocation_status, errmsg=emsg)
        call util_allocate_check (allocation_status, emsg, 'output_typeUnits_withTime_faceR')
        output_typeUnits_withTime_faceR(:) = ""

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine util_allocate_outputML_facetypes
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_allocate_outputML_times ()
        !%-----------------------------------------------------------------------------
        !% allocates the output time array for multi-level output
        !%-----------------------------------------------------------------------------
        integer            :: allocation_status
        integer, pointer   :: nLevel
        character(len=99)  :: emsg
        character(64)      :: subroutine_name = 'utility_allocate_outputtype'
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% --- don't do this is output is suppressed
        if (setting%Output%Report%suppress_MultiLevel_Output) return

        nLevel => setting%Output%StoredLevels

        !% --- bug check
        if (nLevel < 1) then
            write(*,"(A)") 'ERROR (code) the Output.StoredLevels is less than 1...'
            write(*,"(A)") '... which should have caused Output.Report.suppress_Multilevel_Output = .true.'
            write(*,"(A,i8)") '... setting%Output%StoredLevels              ', setting%Output%StoredLevels
            write(*,"(A,i8)") '... etting%Output%Report%suppress_MultiLevel_Output ', setting%Output%Report%suppress_MultiLevel_Output
            stop
        end if

        allocate(output_times(nLevel), stat=allocation_status, errmsg=emsg)
        call util_allocate_check (allocation_status, emsg, 'output_times')
        output_times(:) = nullvalueR

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine util_allocate_outputML_times
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_allocate_outputML_filenames ()
        !%-----------------------------------------------------------------------------
        !% allocates the output filename array for multi-level output
        !%-----------------------------------------------------------------------------
        integer, pointer   :: nLevel
        integer            :: allocation_status
        character(len=99)  :: emsg
        character(64)      :: subroutine_name = 'utility_allocate_outputML_filenames'
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% --- don't do this is output is suppressed
        if (setting%Output%Report%suppress_MultiLevel_Output) return

        nLevel => setting%Output%StoredFileNames

        !% --- bug check
        if (nLevel < 1) then
            write(*,"(A)") 'ERROR (code) the Output.StoredLevels is less than 1...'
            write(*,"(A)") '... which should have caused Output.Report.suppress_Multilevel_Output = .true.'
            write(*,"(A,i8)") '... setting%Output%StoredLevels              ', setting%Output%StoredLevels
            write(*,"(A,i8)") '... etting%Output%Report%suppress_MultiLevel_Output ', setting%Output%Report%suppress_MultiLevel_Output
            stop
        end if

        allocate(output_binary_filenames(nLevel), stat=allocation_status, errmsg=emsg)
        call util_allocate_check (allocation_status, emsg, 'output_binary_filenames')
        output_binary_filenames(:) = "null"

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine util_allocate_outputML_filenames
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_allocate_outputML_storage ()
        !%-----------------------------------------------------------------------------
        !% Description:
        !% Creates multi-time-level storage for output data
        !%-----------------------------------------------------------------------------
        integer, pointer  :: nLevel, nType, nMaxLevel
        integer           :: nElem, nFace, nTotal
        integer           :: allocation_status
        character(len=99) :: emsg
        character(64)     :: subroutine_name = ' util_allocate_outputML_storage'
        !%-----------------------------------------------------------------------------
        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% --- don't do this is output is suppressed
        if (setting%Output%Report%suppress_MultiLevel_Output) return

        !% --- shorthand for the stored levels
        nLevel => setting%Output%StoredLevels

        !% --- get the total number of time levels for the report
        !% --- increase by 2 for start and end files
        if (setting%Output%Report%TimeInterval > zeroR) then
            setting%Output%MaxExpectedLevels = 2 + &
                ceiling((setting%Time%End - setting%Output%Report%StartTime) &
                      /setting%Output%Report%TimeInterval)
        else
            write (*,"(A)") 'The Report Interval (SWMM inp file REPORT_STEP) is less than zero.'
            write (*,"(A)") 'This will suppress multi-level output files'
            setting%Output%Report%suppress_MultiLevel_Output = .true.
            return
        end if

        !% --- check and adjust stored level output so as not to waste memory
        if ( setting%Output%StoredLevels > setting%Output%MaxExpectedLevels+ 2) then
            write (*,"(A)") 'Changing Output.StoredLevels to Output.MaxExpectedLevels...'
            write (*,"(A)") '... based on requested SWMM report step. Old and new stored levels are...'
            write (*,"(2i8)") setting%Output%StoredLevels, setting%Output%MaxExpectedLevels
            setting%Output%StoredLevels = setting%Output%MaxExpectedLevels
        end if

        !% --- bug check
        if (nLevel < 1) then
            write(*,"(A)") 'ERROR (code) the Output.StoredLevels is less than 1...'
            write(*,"(A)") '... which should have caused Output.Report.suppress_Multilevel_Output = .true.'
            write(*,"(A,i8)") '... setting%Output%StoredLevels              ', setting%Output%StoredLevels
            write(*,"(A,i8)") '... etting%Output%Report%suppress_MultiLevel_Output=', setting%Output%Report%suppress_MultiLevel_Output
            stop
        end if

        if (setting%Output%OutputElementsExist) then
            !% --- shortand for the output types
            nType => N_OutTypeElem
            !print*, nType, 'nType in image = ', this_image()
            !% --- bug check
            if (nType < 1) then
                write(*,"(A)") 'ERROR (code) the N_OutTypeElem is less than 1 (i.e. no output types selected)'
                write(*,"(A)") '... which should have caused Output.OutputElementsExist = .false.'
                write(*,"(A,i8)") '... setting%Output%N_OutTypeElem       ', N_OutTypeElem
                write(*,"(A,L)") '... etting%Output%sOutputElementsExist =', setting%Output%OutputElementsExist
                stop
            end if

            !% --- set the maximum number of output elements in any image
            nElem = maxval(N_OutElem(:))

            !% --- bug check
            if (nElem < 1) then
                write(*,"(A)") 'ERROR (code) the maximum number of elements output from an image, ...'
                write(*,"(A)") '... maxval(N_OutElem(images)), is less than 1...'
                write(*,"(A)") '... which should have caused Output.OutputElementsExist = .false.'
                write(*,"(A,L)") '... setting%Output%OutputElementsExist =', setting%Output%OutputElementsExist
                write(*,"(A)") '... full listing of N_OutElem(:)...'
                write(*,"(I8)") N_OutElem(:)
                stop
            end if

            !% allocate the multi-level element storage for each image
            allocate(elemOutR(nElem, nType,nLevel)[*], stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg, 'elemOutR')
            elemOutR(:,:,:) = nullvalueR

            !% ----------------------------
            !% --- Allocate the Element combined storage from all images for stored time steps
            !% ----------------------------

            !% --- get the total number of output elements on all images
            nTotal = sum(N_OutElem(:))

            !% allocate the full network multi-level output array to one processor
            if (this_image() == 1) then
                !% --- get space for combined element data
                allocate(OutElemDataR(nTotal,nType,nLevel), stat=allocation_status, errmsg=emsg)
                call util_allocate_check(allocation_status, emsg, 'OutElemDataR')
                OutElemDataR(:,:,:) = nullvalueR

                !brh rm !% --- get space for the indexes to the elemR() etc array
                !brh rm allocate(OutElemGidx(nTotal), stat=allocation_status, errmsg=emsg)
                !brh rm call util_allocate_check(allocation_status, emsg, 'OutElemGidx')
                !brh rm OutElemGidx(:) = nullvalueI

                !% --- get space for integer data for OutElemFixedI
                allocate(OutElemFixedI(nTotal,Ncol_oefi), stat=allocation_status, errmsg=emsg)
                call util_allocate_check(allocation_status, emsg, 'OutElemFixedI')
                OutElemFixedI(:,:) = nullvalueI
            end if

        end if

        if (setting%Output%OutputFacesExist) then

            !% --- shortand for the output types
            nType => N_OutTypeFace

            !% --- bug check
            if (nType < 1) then
                write(*,"(A)") 'ERROR (code) the N_OutTypeFace is less than 1 (i.e. no output types selected)'
                write(*,"(A)") '... which should have caused Output.OutputFacseExist = .false.'
                write(*,"(A,i8)") '... setting%Output%N_OutTypeNFaces     ', N_OutTypeFace
                write(*,"(A,L)") '... setting%Output%OutputFacesExist =', setting%Output%OutputFacesExist
                stop
            end if


            !% --- set the maximum number of output faces in any image
            nFace = maxval(N_OutFace(:))

            !% --- bug check
            if (nFace < 1) then
                write(*,"(A)") 'ERROR (code) the maximum number of faces output from an image, ...'
                write(*,"(A)") '... maxval(N_OutFace(images)), is less than 1...'
                write(*,"(A)") '... which should have caused Output.OutputFacesExist = .false.'
                write(*,"(A,L)") '... setting%Output%OutputFacesExist =', setting%Output%OutputFacesExist
                write(*,"(A)") '... full listing of N_OutElem(:)...'
                write(*,"(I8)") N_OutElem(:)
                stop
            end if

            !% allocate the multi-level element storage for each image
            allocate(faceOutR(nFace,nType,nLevel)[*], stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg, 'faceOutR')
            faceOutR(:,:,:) = nullvalueR

            !% ----------------------------
            !% --- Allocate the Face combined storage from all images for stored time steps
            !% ----------------------------

            !% --- get the total number of output elements on all images
            nTotal = sum(N_OutFace(:))

            !% allocate the full network multi-level output array to one processor
            if (this_image() == 1) then
                !% --- get space for combined element data
                allocate(OutFaceDataR(nTotal,nType,nLevel), stat=allocation_status, errmsg=emsg)
                call util_allocate_check(allocation_status, emsg, 'OutFaceDataR')
                OutFaceDataR(:,:,:) = nullvalueR

                !brh rm !% --- get space for the indexes to the elemR() etc array
                !brh rm  allocate(OutFaceGidx(nTotal), stat=allocation_status, errmsg=emsg)
                !brh rm  call util_allocate_check(allocation_status, emsg, 'OutFaceGidx')
                !brh rm  OutFaceGidx(:) = nullvalueI

                !% --- get space for integer data for OutElemFixedI
                allocate(OutFaceFixedI(nTotal,Ncol_offi), stat=allocation_status, errmsg=emsg)
                call util_allocate_check(allocation_status, emsg, 'OutFaceFixedI')
                OutFaceFixedI(:,:) = nullvalueI
            end if

        end if

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine util_allocate_outputML_storage
!%
!%==========================================================================
!% PRIVATE
!%==========================================================================
!%
    subroutine util_allocate_col_elemI()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemI is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated ei_... array_index parameter
        !
        !-----------------------------------------------------------------------------
        integer, pointer    :: ncol
        integer             :: ii, jj
        character(64)       :: subroutine_name = 'util_allocate_col_elemI'
        !-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%utility_allocate) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        ncol => Ncol_elemI

        !% allocate an array for storing the column
        allocate( col_elemI(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'col_elemI')

        !% this array can be used as a pointer target in defining masks
        col_elemI(:) = [(ii,ii=1,ncol)]

        !%--------------------------------------------------------------
        !% the code below is a quick print check to see if
        !% the coarray have been set up properly
        ! if (this_image() == 1) then
        !     do jj = 1, num_images()
        !         print*, jj, 'image no'
        !         print*, col_elemI(:)[jj], 'col_elemI(:)[jj]'
        !     end do
        ! end if
        ! print*, 'press return to continue'
        ! read(*,*)
        !%--------------------------------------------------------------

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine util_allocate_col_elemI
!
!==========================================================================
!==========================================================================
!
    subroutine util_allocate_col_elemP()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemP is a vector of the columns in the elemP arrays
        !   that correspond to the enumerated ep_... array_index parameter
        !
        !   the npack_elemP(:) vector contains the number of packed elements
        !   for a given column.
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_elemP'

        !-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% define the maximum number of columns as
        ncol => Ncol_elemP

        !% allocate an array for storing the size of each packed type
        allocate(npack_elemP(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'npack_elemP')

        !% allocate an array for storing the column of each packed type
        allocate( col_elemP(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'col_elemP')

        !% this array can be used as a pointer target in defining masks
        col_elemP(:) = [(ii,ii=1,ncol)]

        !% zero the number of packed items (to be defined in the packing)
        npack_elemP(:) = 0

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_allocate_col_elemP
!
!==========================================================================
!==========================================================================
!
    subroutine util_allocate_col_elemPGalltm()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemPGalltm is a vector of the columns in the elemPGalltm array
        !   that correspond to the enumerated epg_... array_index parameters
        !
        !   the npack_elemPGalltm(:) vector contains the number of packed elements
        !   for a given column.
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_elemPGalltm'

        !-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%utility_allocate) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% define the maximum number of columns as
        ncol => Ncol_elemPGalltm !% whatever the last item in the enumerator

        !% allocate an array for storing the size of each packed type
        allocate( npack_elemPGalltm(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'npack_elemPGalltm')

        !% allocate an array for storing the enum type of each column
        allocate( col_elemPGalltm(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'col_elemPGalltm')

        !% this array can be used as a pointer target in defining masks
        col_elemPGalltm(:) = [(ii,ii=1,ncol)]

        !% zero the number of packed items (to be defined in the packing)
        npack_elemPGalltm(:) = 0

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_allocate_col_elemPGalltm
!
!==========================================================================
!==========================================================================
!
    subroutine util_allocate_col_elemPGac()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemPGac is a vector of the columns in the elemPGac arrays
        !   that correspond to the enumerated epg_... array_index parameters
        !
        !   the npack_elemPGac(:) vector contains the number of packed elements
        !   for a given column.
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_elemPGac'

        !-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%utility_allocate) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% define the maximum number of columns as
        ncol => Ncol_elemPGac!% whatever the last item in the enumerator

        !% allocate an array for storing the size of each packed type
        allocate( npack_elemPGac(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'npack_elemPGac')

        !% allocate an array for storing the enum type of each column
        allocate( col_elemPGac(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'col_elemPGac')

        !% this array can be used as a pointer target in defining masks
        col_elemPGac(:) = [(ii,ii=1,ncol)]

        !% zero the number of packed items (to be defined in the packing)
        npack_elemPGac(:) = 0

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_allocate_col_elemPGac
!
!==========================================================================
!==========================================================================
!
    subroutine util_allocate_col_elemPGetm()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemPGetm is a vector of the columns in the elemPGetm arrays
        !   that correspond to the enumerated epg_... array_index parameters
        !
        !   the npack_elemPGetm(:) vector contains the number of packed elements
        !   for a given column.
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_elemPGetm'

        !-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%utility_allocate) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% define the maximum number of columns as
        ncol => Ncol_elemPGetm   !% whatever the last item in the enumerator

        !% allocate an array for storing the size of each packed type
        allocate( npack_elemPGetm(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'npack_elemPGetm')

        !% allocate an array for storing the enum type of each column
        allocate( col_elemPGetm(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'col_elemPGetm')

        !% this array can be used as a pointer target in defining masks
        col_elemPGetm(:) = [(ii,ii=1,ncol)]

        !% zero the number of packed items (to be defined in the packing)
        npack_elemPGetm(:) = 0

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine util_allocate_col_elemPGetm
!
!==========================================================================
!==========================================================================
!
    subroutine util_allocate_col_elemR()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemR is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated er_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_elemR'

        !-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%utility_allocate) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% define the maximum number of columns as
        ncol => Ncol_elemR

        !% allocate an array for storing the column
        allocate( col_elemR(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'col_elemR')

        !% this array can be used as a pointer target in defining masks
        col_elemR(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_allocate_col_elemR
!
!==========================================================================
!==========================================================================
!
    subroutine util_allocate_col_elemSI()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemSI is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated esi_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_elemSI'

        !-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%utility_allocate) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% define the maximum number of columns as
        ncol => Ncol_elemSI

        !% allocate an array for storing the column
        allocate( col_elemSI(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'col_elemSI')

        !% this array can be used as a pointer target in defining masks
        col_elemSI(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_allocate_col_elemSI
!
!==========================================================================
!==========================================================================
!
    subroutine util_allocate_col_elemSR()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemSR is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated esr_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_elemSR'

        !-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%utility_allocate) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% define the maximum number of columns as
        ncol => Ncol_elemSR

        !% allocate an array for storing the column
        allocate( col_elemSR(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'col_elemSR')

        !% this array can be used as a pointer target in defining masks
        col_elemSR(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_allocate_col_elemSR
!
!==========================================================================
!==========================================================================
!
    subroutine util_allocate_col_elemSGR()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemSGR is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated esgr_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_elemSGR'

        !-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%utility_allocate) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% define the maximum number of columns as
        ncol => Ncol_elemSGR

        !% allocate an array for storing the column
        allocate( col_elemSGR(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'col_elemSGR')

        !% this array can be used as a pointer target in defining masks
        col_elemSGR(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_allocate_col_elemSGR
!
!==========================================================================
!==========================================================================
!
    subroutine util_allocate_col_elemWDI()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemWDI is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated ewdi_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_elemWDI'

        !-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%utility_allocate) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% define the maximum number of columns as
        ncol => Ncol_elemWDI

        !% allocate an array for storing the column
        allocate( col_elemWDI(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'col_elemWDI')

        !% this array can be used as a pointer target in defining masks
        col_elemWDI(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_allocate_col_elemWDI
!
!==========================================================================
!==========================================================================
!
    subroutine util_allocate_col_elemWDR()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemWDR is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated ewdr_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_elemWDI'

        !-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%utility_allocate) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% define the maximum number of columns as
        ncol => Ncol_elemWDR

        !% allocate an array for storing the column
        allocate( col_elemWDR(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'col_elemWDR')

        !% this array can be used as a pointer target in defining masks
        col_elemWDR(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_allocate_col_elemWDR
!
!==========================================================================
!==========================================================================
!
    subroutine util_allocate_col_elemYN()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_elemYN is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated eYN_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_elemYN'

        !-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%utility_allocate) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% define the maximum number of columns as
        ncol => Ncol_elemYN

        !% allocate an array for storing the column
        allocate( col_elemYN(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'col_elemYN')

        !% this array can be used as a pointer target in defining masks
        col_elemYN(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_allocate_col_elemYN
!
!==========================================================================
!==========================================================================
!
    subroutine util_allocate_col_faceI()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_faceI is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated fi_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_faceI'

        !-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%utility_allocate) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% define the maximum number of columns as
        ncol => Ncol_faceI

        !% allocate an array for storing the column
        allocate( col_faceI(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'col_faceI')

        !% this array can be used as a pointer target in defining masks
        col_faceI(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_allocate_col_faceI
!
!==========================================================================
!==========================================================================
!
    subroutine util_allocate_col_faceM()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_faceM is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated fM_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_faceM'

        !-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%utility_allocate) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% define the maximum number of columns as
        ncol => Ncol_faceM

        !% allocate an array for storing the column
        allocate( col_faceM(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'col_faceM')

        !% this array can be used as a pointer target in defining masks
        col_faceM(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_allocate_col_faceM
!
!==========================================================================
!==========================================================================
!
    subroutine util_allocate_col_faceP()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   packed arrays for faces
        !   the col_faceP is a vector of the columns in the faceP arrays
        !   that correspond to the enumerated fp_... array_index parameters
        !
        !   the npack_faceP(:) vector contains the number of packed elements
        !   for a given column.
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_faceP'

        !-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%utility_allocate) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% define the maximum number of columns as
        ncol => Ncol_faceP

        !% allocate an array for storing the size of each packed type
        allocate( npack_faceP(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'npack_faceP')

        !% allocate an array for storing the column of each packed type
        allocate( col_faceP(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'col_faceP')

        !% this array can be used as a pointer target in defining masks
        col_faceP(:) = [(ii,ii=1,ncol)]

        !% zero the number of packed items (to be defined in the packing)
        npack_faceP(:) = 0

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_allocate_col_faceP
!
!==========================================================================
!==========================================================================
!
    subroutine util_allocate_col_facePS()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   packed arrays for the shared (internal boundary)faces
        !   the col_facePS is a vector of the columns in the facePS arrays
        !   that correspond to the enumerated fp_... array_index parameters
        !   col_facePS has the same number of columns as col_faceP because
        !   all the packs for internal faces are needed for shared faces as
        !   well.
        !
        !   the npack_facePS(:) vector contains the number of packed elements
        !   for a given column.
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_facePS'

        !-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%utility_allocate) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% define the maximum number of columns as
        ncol => Ncol_faceP

        !% allocate an array for storing the size of each packed type
        allocate( npack_facePS(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'npack_facePS')

        !% allocate an array for storing the column of each packed type
        allocate( col_facePS(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'col_facePS')

        !% this array can be used as a pointer target in defining masks
        col_facePS(:) = [(ii,ii=1,ncol)]

        !% zero the number of packed items (to be defined in the packing)
        npack_facePS(:) = 0

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_allocate_col_facePS
!
!==========================================================================
!==========================================================================
!
    subroutine util_allocate_col_faceR()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_faceR is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated fr_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_faceR'

        !-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%utility_allocate) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% define the maximum number of columns as
        ncol => Ncol_faceR  !global

        !% allocate an array of column indexes that can be used as targets of pointers
        allocate( col_faceR(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'col_faceR')

        !% this array can be used as a pointer target in defining masks
        col_faceR(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_allocate_col_faceR
!
!==========================================================================
!==========================================================================
!
    subroutine util_allocate_col_faceYN()
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   the col_faceYN is a vector of the columns in the faceR arrays
        !   that correspond to the enumerated fYN_... array_index parameter
        !
        !-----------------------------------------------------------------------------

        integer, pointer    :: ncol
        integer             :: ii
        character(64)       :: subroutine_name = 'util_allocate_col_faceYN'

        !-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%utility_allocate) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% define the maximum number of columns as
        ncol => Ncol_faceYN  !global

        !% allocate an array of column indexes that can be used as targets of pointers
        allocate( col_faceYN(ncol)[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'col_faceYN')

        !% this array can be used as a pointer target in defining masks
        col_faceYN(:) = [(ii,ii=1,ncol)]

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_allocate_col_faceYN
!
!==========================================================================
!==========================================================================
!
    subroutine util_allocate_curves ()
        !-----------------------------------------------------------------------------
        !
        !
        !-----------------------------------------------------------------------------

        character(64)       :: subroutine_name = 'util_allocate_curves'

        !-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%utility_allocate) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% allocate curves
        allocate( curve(N_curve), stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'curve')

        curve%ID = nullvalueI
        curve%Type = nullvalueI
        curve%RefersTo = nullvalueI
        curve%NumRows = nullvalueI
        curve%ElemIdx = nullvalueI
        curve%FaceIdx = nullvalueI

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_allocate_curves
!
!==========================================================================
!==========================================================================
!
    subroutine util_allocate_curve_entries (curve_idx, num_entries)
        !-----------------------------------------------------------------------------
        !
        !
        !-----------------------------------------------------------------------------

        integer, intent(in) :: curve_idx, num_entries
        character(64)       :: subroutine_name = 'util_allocate_curve_entries'

        !-----------------------------------------------------------------------------
        if (icrash) return
        if (setting%Debug%File%utility_allocate) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% allocate the value array of curve
        allocate( curve(curve_idx)%ValueArray(num_entries,Ncol_curve), stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'curve')

        curve(curve_idx)%ValueArray = nullvalueR

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_allocate_curve_entries
!
!==========================================================================
!==========================================================================
!
    subroutine util_allocate_check(allocation_status, emsg, locationstring)
        !-----------------------------------------------------------------------------
        !
        ! Description:
        !   Checks allocation status and stops if there is an error
        !
        !-----------------------------------------------------------------------------

            integer,           intent(in   ) :: allocation_status
            character(len=*),  intent(in   ) :: emsg
            character(len=*),   intent(in   ) :: locationstring !% unique identifier of location

            character(64):: subroutine_name = 'util_allocate_check'

        !-----------------------------------------------------------------------------
            if (icrash) return
            if (setting%Debug%File%utility) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

            if (allocation_status > 0) then
                print *, allocation_status
                print *, trim(emsg)
                print *, 'variable trying to allocate = ',trim(locationstring)
                stop
            end if

            if (setting%Debug%File%utility) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_allocate_check

end module utility_allocate
