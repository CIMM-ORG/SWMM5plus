module utility_allocate

    use define_globals
    use define_keys
    use define_indexes
    use define_settings, only: setting
    use interface_
    use utility
    use utility_crash, only: util_crashpoint

    ! use utility, only: utility_check_allocation

    implicit none

!-----------------------------------------------------------------------------
! Description:
!   This has all the large array allocation. The only other allocation
!   should be in boundary conditions (module bc)
!   All the top-level storage should be allocated in this module
!-----------------------------------------------------------------------------

    private

    ! utility_allocate constants
    integer           :: allocation_status
    character(len=99) ::              emsg

    ! public members
    public :: util_allocate_scalar_for_images
    public :: util_allocate_secondary_coarrays
    public :: util_allocate_linknode
    public :: util_allocate_monitor_points
    public :: util_allocate_action_points
    public :: util_allocate_link_transect
    public :: util_allocate_element_transect
    public :: util_allocate_uniformtable_array
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
    public :: util_allocate_boundary_ghost_elem_array
    public :: util_allocate_temporary_arrays
    public :: util_allocate_output_profiles
    


contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
    subroutine util_allocate_scalar_for_images ()    
        !%------------------------------------------------------------------
        !% Description
        !% allocates coarrays that appear as scalars on an image
        !% These are NOT coarrays, but are indexed off the present image
        !%------------------------------------------------------------------

        allocate(N_elem(num_images()))
        allocate(N_face(num_images()))
        allocate(N_unique_face(num_images()))

    end subroutine util_allocate_scalar_for_images
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_allocate_secondary_coarrays ()
        !%------------------------------------------------------------------
        !% Description:
        !% allocates and initializes secondary coarrays (i.e., not the
        !% primary coarrays of elemX, faceX, etc.)
        !%------------------------------------------------------------------

        ! !% --- allocate a corray to store the number of control/monitoring
        ! !%     points on each image
        ! allocate(N_ConMon(num_images())[*],stat=allocation_status, errmsg= emsg)
        ! call util_allocate_check (allocation_status, emsg, 'N_ConMon')
        ! N_ConMon = zeroI

        !% --- allocate a coarray for the number of output elements on each image
        !%     the array is first used in outputML_size_OutElem_by_image
        allocate(N_OutElem(num_images())[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'N_OutElem')
        N_OutElem = zeroI
        !% --- similar coarray for number of faces that will be output for each image
        allocate(N_OutFace(num_images())[*], stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'N_OutFace')
        N_OutFace = zeroI

        !% --- allocate coarrays for the setting%Output%... controls
        allocate(setting%Output%ElementsExist_byImage(num_images())[*], &
                stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'setting%Output%ElementsExist_byImage')
        setting%Output%ElementsExist_byImage = .false.

        allocate(setting%Output%FacesExist_byImage(num_images())[*], &
                stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'setting%Output%FacesExist_byImage')
        setting%Output%FacesExist_byImage = .false.

        !% --- allocate coarrays for communicating monitor image data to action image
        !%     for controls
        allocate(monitorPassR(Ncol_MonitoringPointR)[*], &
                stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'monitorPassR')
        setting%Output%FacesExist_byImage = .false.

        !%------------------------------------------------------------------ 
    end subroutine util_allocate_secondary_coarrays
!%    
!%==========================================================================
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
            !if (crashYN) return
            if (setting%Debug%File%utility_allocate) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
        !% If BIPquick is being used for Partitioning, include additional rows to the link-node arrays
        if (setting%Partitioning%PartitioningMethod == BQuick) then
            additional_rows = num_images() - 1
        end if

        allocate(node%I(N_node + additional_rows, Ncol_nodeI)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'node%I')
        node%I(:,:) = nullvalueI

        allocate(link%I(setting%SWMMinput%N_link + additional_rows, Ncol_linkI)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'link%I')
        link%I(:,:) = nullvalueI

        allocate(node%R(N_node + additional_rows, Ncol_nodeR)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'node%R')
        node%R(:,:) = nullvalueR

        allocate(link%R(setting%SWMMinput%N_link + additional_rows, Ncol_linkR)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'link%R')
        link%R(:,:) = nullvalueR

        allocate(node%YN(N_node + additional_rows, Ncol_nodeYN)[*], stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'node%YN')
        node%YN(:,:) = nullvalueL

        allocate(link%YN(setting%SWMMinput%N_link + additional_rows, Ncol_linkYN)[*], stat=allocation_status, errmsg=emsg)
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

        allocate(link_output_idx(setting%SWMMinput%N_link + additional_rows), stat=allocation_status,errmsg=emsg)
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
    subroutine util_allocate_monitor_points()
        !%------------------------------------------------------------------
        !% Description
        !% allocates storage space for monitoring elements
        !% Note these are NOT coarrays
        !%------------------------------------------------------------------
        !% Declarations:
            integer, pointer :: nelem, ncol
        !%------------------------------------------------------------------
        !% Aliases
            nelem => N_MonitorPoint  !% number of control/monitoring points
        !%------------------------------------------------------------------

        allocate (monitorI(nelem,Ncol_MonitoringPointI), stat=allocation_status, errmsg=emsg)   
        call util_allocate_check(allocation_status, emsg, 'monitorI') 
        monitorI(:,:) = nullvalueI

        allocate (monitorR(nelem,Ncol_MonitoringPointR), stat=allocation_status, errmsg=emsg)   
        call util_allocate_check(allocation_status, emsg, 'monitorR') 
        monitorR(:,:) = nullvalueR

        ! allocate (monitorYN(nelem,Ncol_MonitoringPointYN), stat=allocation_status, errmsg=emsg)   
        ! call util_allocate_check(allocation_status, emsg, 'monitorYN') 
        ! monitorYN(:,:) = .false. 
        
        ! !% --- correspondence between conmonR columns and elemR columns
        ! allocate(cmR_eR_col(Ncol_ConMonR), stat=allocation_status, errmsg=emsg)
        ! call util_allocate_check(allocation_status, emsg, 'cmR_eR_col')
        ! !% --- initialize the packed array of columns that link the columns
        ! !%     in the conmonR to the elemR array. These must match the cmr_...
        ! !%     definitions in define_indexes
        ! cmR_eR_col(:) = (/ er_Depth, er_Flowrate, er_Head, er_Velocity, er_Volume /)

        ! !% --- integer data
        ! ncol  => Ncol_conmonI !% the maximum number of columns
        ! allocate(conmonI(nelem, ncol), stat=allocation_status, errmsg=emsg)
        ! call util_allocate_check(allocation_status, emsg, 'conmonI')
        ! conmonI(:,:) = nullValueI

        ! !% --- real data
        ! ncol  => Ncol_conmonR !% the maximum number of columns
        ! allocate(conmonR(nelem, ncol), stat=allocation_status, errmsg=emsg)
        ! call util_allocate_check(allocation_status, emsg, 'conmonR')
        ! conmonR(:,:) = nullValueR

        ! !% --- logical data
        ! ncol  => Ncol_conmonYN !% the maximum number of columns
        ! allocate(conmonYN(nelem, ncol), stat=allocation_status, errmsg=emsg)
        ! call util_allocate_check(allocation_status, emsg, 'conmonYN')
        ! conmonYN(:,:) = .false.

    end subroutine util_allocate_monitor_points
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_allocate_action_points ()
        !%------------------------------------------------------------------
        !% Description
        !% allocates storage space for monitoring elements
        !% Note these are NOT coarrays
        !%------------------------------------------------------------------
        !% Declarations:
        integer, pointer :: nelem, ncol
        !%------------------------------------------------------------------
        !% Aliases
            nelem => N_ActionPoint  !% number of control/monitoring points
        !%------------------------------------------------------------------

        allocate (actionI(nelem,Ncol_ActionPointI), stat=allocation_status, errmsg=emsg)   
        call util_allocate_check(allocation_status, emsg, 'actionI') 
        actionI(:,:) = nullvalueI

        allocate (actionR(nelem,Ncol_ActionPointR), stat=allocation_status, errmsg=emsg)   
        call util_allocate_check(allocation_status, emsg, 'actionR') 
        actionR(:,:) = nullvalueR

    end subroutine util_allocate_action_points
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_allocate_link_transect ()

        integer :: ii

        !% 2D transectI
        allocate(link%transectI(setting%SWMMinput%N_link_transect,Ncol_transectI),stat=allocation_status,errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'transectI')
        link%transectI(:,:) = nullvalueI

        !% 2D transectR
        allocate(link%transectR(setting%SWMMinput%N_link_transect,Ncol_transectR),stat=allocation_status,errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'transectR')
        link%transectR(:,:) = nullvalueR

        !% 3D transectTableDepthR -- data by uniform depth depth
        allocate(link%transectTableDepthR(setting%SWMMinput%N_link_transect,N_transect_depth_items,Ncol_transectTable),stat=allocation_status,errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'transectTableDepthR')
        link%transectTableDepthR(:,:,:) = nullvalueR
        
        !% 3D transectTableAreaA -- data by uniform area discretization
        allocate(link%transectTableAreaR(setting%SWMMinput%N_link_transect,N_transect_area_items,Ncol_transectTable),stat=allocation_status,errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'transectTableAreaR')
        link%transectTableAreaR(:,:,:) = nullvalueR

        !% 1D transectID
        allocate(link%transectID(setting%SWMMinput%N_link_transect),stat=allocation_status,errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'transectID')
        
        !% --- HACK -- using fixed len=16 for transectID. Should be changed
        !%     similar to allocation of link%Names to get length of ID
        do ii = 1,setting%SWMMinput%N_link_transect
            allocate(character(16) :: link%transectID(ii)%str, stat=allocation_status, errmsg=emsg)
            call util_allocate_check(allocation_status, emsg, 'character(32) :: link%transectID%str')
            link%transectID(ii)%str = '0'
        end do


    end subroutine util_allocate_link_transect    
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_allocate_element_transect ()

        !% 2D transectI
        allocate(transectI(N_transect,Ncol_transectI),stat=allocation_status,errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'transectI')
        transectI(:,:) = nullvalueI

        !% 2D transectR
        allocate(transectR(N_transect,Ncol_transectR),stat=allocation_status,errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'transectR')
        transectR(:,:) = nullvalueR

        !% 3D transectTableDepthR -- data by uniform depth depth
        allocate(transectTableDepthR(N_transect,N_transect_depth_items,Ncol_transectTable),stat=allocation_status,errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'transectTableDepthR')
        transectTableDepthR(:,:,:) = nullvalueR
        
        !% 3D transectTableAreaA -- data by uniform area discretization
        allocate(transectTableAreaR(N_transect,N_transect_area_items,Ncol_transectTable),stat=allocation_status,errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'transectTableAreaR')
        transectTableAreaR(:,:,:) = nullvalueR

        !% 1D transectID
        allocate(transectID(N_transect),stat=allocation_status,errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'transectID')
        transectID(:) = ""

    end subroutine util_allocate_element_transect    
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_allocate_uniformtable_array ()
        !%------------------------------------------------------------------
        !% Description:
        !% Allocates space for the lookup table for geometry = f(sectionfactor), etc.
        !%------------------------------------------------------------------
        !%------------------------------------------------------------------
        !% --- only headBC locations at this time 20220726brh
        !%     HACK -- to be expanded later to include transects and other geometry tables
        N_uniformTable_locations = N_headBC

        !% --- storage for integer data
        allocate(uniformTableI &
                (N_uniformTable_locations,Ncol_uniformTableI), &
                stat=allocation_status,errmsg=emsg)

        call util_allocate_check(allocation_status, emsg, 'uniformTableI')
        uniformTableI(:,:) = nullvalueI

        !% --- storage real data
        allocate(uniformTableR &
                (N_uniformTable_locations,Ncol_uniformTableR), &
                stat=allocation_status,errmsg=emsg)
                
        call util_allocate_check(allocation_status, emsg, 'uniformTableR')
        uniformTableR(:,:) = nullvalueR

        !% ---3D table indexed by uniformly-divided section factor
        allocate(uniformTableDataR &
                 (N_uniformTable_locations, N_uniformTableData_items, Ncol_uniformTableDataR), &
                  stat=allocation_status,errmsg=emsg)

        call util_allocate_check(allocation_status, emsg, 'uniformTableDataR')
        uniformTableDataR(:,:,:) = nullvalueR

    end subroutine util_allocate_uniformtable_array
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
            !if (crashYN) return
            if (setting%Debug%File%utility_allocate) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------

        !% subcatchR
        allocate(subcatchR(setting%SWMMinput%N_subcatch, Ncol_subcatchR), stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'subcatchR')
        subcatchR(:,:) = nullvalueR

        !% subcatchI
        allocate(subcatchI(setting%SWMMinput%N_subcatch, Ncol_subcatchI), stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'subcatchI')
        subcatchR(:,:) = nullvalueI

        !% subcatchYN
        allocate(subcatchYN(setting%SWMMinput%N_subcatch, Ncol_subcatchYN), stat=allocation_status, errmsg=emsg)
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
        !if (crashYN) return
        allocate(adjacent_links(max_branch_per_node))
        allocate(elem_per_image(num_images()))
        allocate(image_full(num_images()))

        !% If BIPquick is being used for Partitioning, allocate additional arrays
        if (setting%Partitioning%PartitioningMethod == BQuick) then
            call util_count_node_types(N_nBCup, N_nBCdn, N_nJm, N_nStorage, N_nJ2, N_nJ1)

            allocate(B_nodeI(size(node%I,1), max_up_branch_per_node))
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
        !if (crashYN) return
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

        ! ncol=> Ncol_faceM
        ! allocate(faceM(max_caf_face_N, ncol)[*], stat=allocation_status, errmsg=emsg)
        ! call util_allocate_check(allocation_status, emsg, 'faceM')
        ! faceM(:,:) = nullvalueL

        !% HOLD FOR LATER IMPLEMENTATION 20220930
        ! !% ---3D geometry table for fast look up by element
        ! allocate(geometryTableR &
        !          (max_caf_elem_N+N_dummy_elem, Ncol2_GeometryTableR, Ncol3_GeometryTableR), &
        !           stat=allocation_status,errmsg=emsg)
        ! geometryTableR(:,:,:) = nullvalueR

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine util_allocate_elemX_faceX
!
!==========================================================================
!==========================================================================
!
    subroutine util_allocate_columns()
        !%-----------------------------------------------------------------
        !% Description:
        !%   All the enumerated variables can not be used as pointers. Thus the
        !%   column variables are stored in col_elemX(:) arrays that are
        !%   allowed to be targets of pointers.
        !%-----------------------------------------------------------------
            character(64) :: subroutine_name = 'util_allocate_columns'
        !%-----------------------------------------------------------------
        if (setting%Debug%File%utility_allocate) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        !% allocation of the col_elemX and npack_elemX
        call util_allocate_col_elemI ()
        call util_allocate_col_elemP ()
        call util_allocate_col_elemPGalltm ()
        call util_allocate_col_elemPGac ()
        call util_allocate_col_elemPGetm ()
        call util_allocate_col_elemR ()
        call util_allocate_col_elemSI ()
        call util_allocate_col_elemSR ()
        call util_allocate_col_elemSGR ()
        !call util_allocate_col_elemWDI
        !call util_allocate_col_elemWDR
        call util_allocate_col_elemYN ()
        call util_allocate_col_faceI ()
        !call util_allocate_col_faceM
        call util_allocate_col_faceP ()
        call util_allocate_col_facePS ()
        call util_allocate_col_faceR ()
        call util_allocate_col_faceYN ()
        ! call util_allocate_col_conmonI ()
        ! call util_allocate_col_conmonR ()
        ! call util_allocate_col_conmonYN ()

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

            allocate(BC%headR(N_headBC, N_headR), stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg, 'BC%headR')
            BC%headR(:,:) = nullvalueR

            allocate(BC%headTimeseries(N_headBC, setting%BC%TimeSlotsStored, N_headR_TS), stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg, 'BC%headTimeseries')
            BC%headTimeseries(:,:,:) = nullvalueR

            !allocate(BC%headI(:,bi_TS_upper_idx)(N_headBC), stat=allocation_status, errmsg=emsg)
            !call util_allocate_check (allocation_status, emsg, 'BC%headI(:,bi_TS_upper_idx)')
            !BC%headI(:,bi_TS_upper_idx)(:) = nullvalueI

            !allocate(BC%head_value(N_headBC), stat=allocation_status, errmsg=emsg)
            !call util_allocate_check (allocation_status, emsg, 'BC%head_value')
            !BC%head_value(:) = nullvalueR

            !allocate(BC%hasFlapGateYN(N_headBC), stat=allocation_status, errmsg=emsg)
            !call util_allocate_check (allocation_status, emsg, 'BC%hasFlapGateYN')
            !BC%hasFlapGateYN = .false.

        end if

        if (N_flowBC > 0) then
            allocate(BC%flowI(N_flowBC, N_flowI), stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg, 'BC%flowI')
            BC%flowI(:,:) = nullvalueI

            allocate(BC%flowYN(N_flowBC, N_flowYN), stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg, 'BC%flowYN')
            BC%flowYN(:,:) = nullvalueL

            allocate(BC%flowR(N_flowBC, N_flowR), stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg, 'BC%flowR')
            BC%flowR(:,:) = nullvalueR

            allocate(BC%flowTimeseries(N_flowBC, setting%BC%TimeSlotsStored, N_flowR_TS), stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg, 'BC%flowTimeseries')
            BC%flowTimeseries(:,:,:) = nullvalueR
            
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
        !if (crashYN) return
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
        !%------------------------------------------------------------------
        !% Description:
        !% allocates the output type arrays for the element data that are output
        !%------------------------------------------------------------------
        !% Declarations:
            integer            :: allocation_status
            character(len=99)  :: emsg
            character(64)       :: subroutine_name = 'util_allocate_outputML_elemtypes'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (setting%Output%Report%suppress_MultiLevel_Output) return
            if (.not. setting%Output%ElementsExist_global) return
            if (setting%Debug%File%utility_allocate) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% --- bug check
        if (N_OutTypeElem < 1) then
            write(*,"(A)") 'ERROR (code) the N_OutTypeElem is less than 1 (i.e. no output types selected)'
            write(*,"(A)") '... which should have caused Output.ElementsExist_glboal = .false.'
            write(*,"(A,i8)") '... setting%Output%N_OutTypeElem      ', N_OutTypeElem
            write(*,"(A,i8)") '... etting%Output%ElementsExist_glboal ', setting%Output%ElementsExist_global
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

        !% --- allocate the output type logical for whether this output is multiplied by number of barrels
        allocate(output_typeMultiplyByBarrels_elemR(N_OutTypeElem), stat=allocation_status, errmsg=emsg)
        call util_allocate_check (allocation_status, emsg, 'output_typeMultiplyByBarrels_elemR')
        output_typeMultiplyByBarrels_elemR(:) = zeroI

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

        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%utility_allocate) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine util_allocate_outputML_elemtypes
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_allocate_outputML_facetypes ()
        !%------------------------------------------------------------------
        !% Description:
        !% allocates the output type arrays for the face data that are output
        !%------------------------------------------------------------------
        !% Declarations:
            integer            :: allocation_status
            character(len=99)  :: emsg
            character(64)       :: subroutine_name = 'util_allocate_outputML_facetypes'
        !%-------------------------------------------------------------------
        !% Preliminaries:
            if (setting%Output%Report%suppress_MultiLevel_Output) return
            if (.not. setting%Output%FacesExist_global) return
            if (setting%Debug%File%utility_allocate) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-------------------------------------------------------------------
        !% --- bug check
        if (N_OutTypeFace < 1) then
            write(*,"(A)") 'ERROR (code) the N_OutTypeFace is less than 1 (i.e. no output types selected)'
                write(*,"(A)") '... which should have caused Output.ElementsExist_global = .false.'
                write(*,"(A,i8)") '... setting%Output%N_OutTypeFace      ', N_OutTypeFace
                write(*,"(A,i8)") '... etting%Output%ElementsExist_byImage ', setting%Output%FacesExist_global
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

        !% --- allocate the output logical for whether this output is multiplied by number of barrels
        allocate(output_typeMultiplyByBarrels_faceR(N_OutTypeFace), stat=allocation_status, errmsg=emsg)
        call util_allocate_check (allocation_status, emsg, 'output_typeMultiplyByBarrels_faceR')
        output_typeMultiplyByBarrels_faceR(:) = zeroI

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

        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%utility_allocate) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine util_allocate_outputML_facetypes
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_allocate_outputML_times ()
        !%------------------------------------------------------------------
        !% Description:
        !% allocates the output time array for multi-level output
        !%------------------------------------------------------------------
        !% Declarations:
            integer            :: allocation_status
            integer, pointer   :: nLevel
            character(len=99)  :: emsg
            character(64)      :: subroutine_name = 'utility_allocate_outputtype'
        !%------------------------------------------------------------------
        !% Preliminaries:
            if (setting%Output%Report%suppress_MultiLevel_Output) return
            if (setting%Debug%File%utility_allocate) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
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

        !%------------------------------------------------------------------
        !% Closing
            if (setting%Debug%File%utility_allocate) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine util_allocate_outputML_times
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_allocate_outputML_filenames ()
        !%------------------------------------------------------------------
        !% Description:
        !% allocates the output filename array for multi-level output
        !%------------------------------------------------------------------
        !% Declarations:
            integer, pointer   :: nLevel
            integer            :: allocation_status
            character(len=99)  :: emsg
            character(64)      :: subroutine_name = 'utility_allocate_outputML_filenames'
        !%------------------------------------------------------------------
        !% Preliminaries:
            if (setting%Output%Report%suppress_MultiLevel_Output) return
            if (setting%Debug%File%utility_allocate) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
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

        !%------------------------------------------------------------------
        !% Closing:
            if (setting%Debug%File%utility_allocate) &
                write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine util_allocate_outputML_filenames
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_allocate_outputML_storage ()
        !%-------------------------------------------------------------------
        !% Description:
        !% Creates multi-time-level storage for output data
        !%-------------------------------------------------------------------
            integer, pointer  :: nLevel, nTypeElem, nTypeFace, nMaxLevel
            integer           :: nTotalElem, nTotalFace, nMaxElem, nMaxFace
            integer           :: allocation_status
            character(len=99) :: emsg
            character(64)     :: subroutine_name = ' util_allocate_outputML_storage'
        !%--------------------------------------------------------------------
            if (setting%Output%Report%suppress_MultiLevel_Output) return
            if (setting%Debug%File%utility_allocate) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%--------------------------------------------------------------------
        !% --- shorthand for the stored levels
        nLevel         => setting%Output%StoredLevels

        !% --- get the total number of time levels for the report
        !% --- increase by 2 for start and end files
        if (setting%Output%Report%TimeInterval > zeroR) then
            setting%Output%MaxExpectedLevels = 2 + &
                ceiling((setting%Time%End - setting%Output%Report%StartTime) &
                      /setting%Output%Report%TimeInterval)
        else
            !% NOTE -- the suppress option needs to always be invoked consistently
            !% across all processor images.
            write (*,"(A)") 'The Report Interval (SWMM inp file REPORT_STEP) is less than zero.'
            write (*,"(A)") 'This will suppress multi-level output files'
            setting%Output%Report%suppress_MultiLevel_Output = .true.
            return
        end if

        !% --- check and adjust stored level output so as not to waste memory
        if ( setting%Output%StoredLevels > setting%Output%MaxExpectedLevels+ 2) then
            if (this_image() == 1) then
                write (*,"(A,i5)") ' ... changing output levels stored before writing; originally: ',nLevel 
                write (*,"(A,i5)") '                Now using the max expected output time levels: ',setting%Output%MaxExpectedLevels+2
            end if
            !% --- reset the nLevel global setting
            nLevel = setting%Output%MaxExpectedLevels + 2
        end if

        !% --- bug check
        if (nLevel < 1) then
            write(*,"(A)") 'ERROR (code) the Output.StoredLevels is less than 1...'
            write(*,"(A)") '... which should have caused Output.Report.suppress_Multilevel_Output = .true.'
            write(*,"(A,i8)") '... setting%Output%StoredLevels              ', setting%Output%StoredLevels
            write(*,"(A,i8)") '... etting%Output%Report%suppress_MultiLevel_Output=', setting%Output%Report%suppress_MultiLevel_Output
            stop
        end if

        !% --- Allocation for Local Element Output
        if (setting%Output%ElementsExist_global) then
            !% --- shortand for the output types
            nTypeElem => N_OutTypeElem
            !print*, nType, 'nType in image = ', this_image()

            !% --- bug check
            if (nTypeElem < 1) then
                write(*,"(A)") 'ERROR (code) the N_OutTypeElem is less than 1 (i.e. no output types selected)'
                write(*,"(A)") '... which requires Output.ElementsExist_global = .false.'
                write(*,"(A,i8)") '... setting%Output%N_OutTypeElem       ', N_OutTypeElem
                write(*,"(A,L)") '... setting%Output%sElementsExist_global =', setting%Output%ElementsExist_global
                stop 984754
            end if

            !% --- set the maximum number of output elements in any image
            nMaxElem = maxval(N_OutElem(:))

            !% --- bug check
            if (nMaxElem < 1) then
                write(*,"(A)") 'ERROR (code) the maximum number of elements output from an image, ...'
                write(*,"(A)") '... maxval(N_OutElem(images)), is less than 1...'
                write(*,"(A)") '... which should have caused Output.ElementsExist_global = .false.'
                write(*,"(A,L)") '... setting%Output%ElementsExist_byImage =', setting%Output%ElementsExist_global
                write(*,"(A)") '... full listing of N_OutElem(:)...'
                write(*,"(I8)") N_OutElem(:)
                stop 387053
            end if

            if (int(nMaxElem,8) * int(nTypeElem,8) * int(nLevel,8) > setting%Output%MemoryStorageMax) then
                print *, 'CONFIGURATION ERROR: the output stored is probably too large'
                print *, '  based on the value in setting.Output.MemoryStorageMax.'
                print *, 'The number of time levels stored before writing (nLevel) is ',nLevel
                print *, '  which is set in setting.Output.StoredLevels of JSON file. '
                print *, 'There are ',nMaxElem,' elements and ',nTypeElem, 'element types'
                print *, '  such that nLevel * nMaxElem * nTypeElem = ',int(nMaxElem,8) * int(nTypeElem,8) * int(nLevel,8)
                print *, 'Suggest reducing setting.Output.StoredLevels below ', &
                            setting%Output%MemoryStorageMax / (int(nMaxElem,8) * int(nTypeElem,8) )
                call util_crashpoint(5598723)
            end if


            !% --- allocate the local  multi-level element storage for each image
            !%    note that EVERY image allocates this array based on the maximum 
            !%    number of elements on ANY image, even if the local image has zero
            !%    output elements
            allocate(elemOutR(nMaxElem, nTypeElem, nLevel)[*], stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg, 'elemOutR')
            elemOutR(:,:,:) = nullvalueR

            !% ----------------------------
            !% --- Allocate the Element combined storage from all images for stored time steps
            !% ----------------------------
            !% --- get the total number of output elements on all images
            nTotalElem = sum(N_OutElem(:))
            !% allocate the full network multi-level output array to one processor
            if (this_image() == 1) then
                !% --- get space for combined element data on image(1)
                allocate(OutElemDataR(nTotalElem,nTypeElem,nLevel), stat=allocation_status, errmsg=emsg)
                call util_allocate_check(allocation_status, emsg, 'OutElemDataR')
                OutElemDataR(:,:,:) = nullvalueR

                !% --- get space for integer data for OutElemFixedI on image(1)
                allocate(OutElemFixedI(nTotalElem,Ncol_oefi), stat=allocation_status, errmsg=emsg)
                call util_allocate_check(allocation_status, emsg, 'OutElemFixedI')
                OutElemFixedI(:,:) = nullvalueI

                !% --- allocation needed for dummy space element conversion on image(1)
                allocate(thisElementOut(nMaxElem), stat=allocation_status, errmsg=emsg)
                call util_allocate_check(allocation_status, emsg, 'thisElementOut')
                thisElementOut(:) = nullvalueI
            end if
        end if

        !% --- Allocation for Faces output
        if (setting%Output%FacesExist_global) then

            !% --- shortand for the output types
            nTypeFace => N_OutTypeFace

            !% --- bug check
            if (nTypeFace < 1) then
                write(*,"(A)") 'ERROR (code) the N_OutTypeFace is less than 1 (i.e. no output types selected)'
                write(*,"(A)") '... which should have caused Output.OutputFacseExist_global = .false.'
                write(*,"(A,i8)") '... setting%Output%N_OutTypeNFaces     ', N_OutTypeFace
                write(*,"(A,L)") '... setting%Output%FacesExist_global =', setting%Output%FacesExist_global
                stop 883753
            end if

            !% --- set the maximum number of output faces in any image
            nMaxFace = maxval(N_OutFace(:))

            !% --- bug check
            if (nMaxFace < 1) then
                write(*,"(A)") 'ERROR (code) the maximum number of faces output from an image, ...'
                write(*,"(A)") '... maxval(N_OutFace(images)), is less than 1...'
                write(*,"(A)") '... which should have caused Output.FacesExist_global = .false.'
                write(*,"(A,L)") '... setting%Output%FacesExist_byImage =', setting%Output%FacesExist_global
                write(*,"(A)") '... full listing of N_OutElem(:)...'
                write(*,"(I8)") N_OutElem(:)
                stop 448723
            end if

            !% allocate the multi-level element storage for each image
            !%    note that EVERY image allocates this array based on the maximum 
            !%    number of faces on ANY image, even if the local image has zero
            !%    output faces
            allocate(faceOutR(nMaxFace,nTypeFace,nLevel)[*], stat=allocation_status, errmsg=emsg)
            call util_allocate_check (allocation_status, emsg, 'faceOutR')
            faceOutR(:,:,:) = nullvalueR

            !% ----------------------------
            !% --- Allocate the Face combined storage from all images for stored time steps
            !% ----------------------------
            !% --- get the total number of output elements on all images
            nTotalFace = sum(N_OutFace(:))

            !% allocate the full network multi-level output array to one processor
            if (this_image() == 1) then
                !% --- get space for combined element data
                allocate(OutFaceDataR(nTotalFace,nTypeFace,nLevel), stat=allocation_status, errmsg=emsg)
                call util_allocate_check(allocation_status, emsg, 'OutFaceDataR')
                OutFaceDataR(:,:,:) = nullvalueR

                !% --- get space for integer data for OutElemFixedI
                allocate(OutFaceFixedI(nTotalFace,Ncol_offi), stat=allocation_status, errmsg=emsg)
                call util_allocate_check(allocation_status, emsg, 'OutFaceFixedI')
                OutFaceFixedI(:,:) = nullvalueI

                !% --- allocation needed for dummy space face conversion on image(1)
                allocate(thisFaceOut(nMaxFace), stat=allocation_status, errmsg=emsg)
                call util_allocate_check(allocation_status, emsg, 'thisFaceOut')
                thisFaceOut(:) = nullvalueI

                

            end if
        end if

        if (     (setting%Output%ElementsExist_global)    &
            .or. (setting%Output%FacesExist_global)    ) then
            !% --- allocate space for storing the number of elements in each link (including non-output links)
            !% --- note this MUST be the size of the setting%SWMMinput%N_link so that we can later pack indexes
            allocate(SWMMlink_num_elements(setting%SWMMinput%N_link), stat=allocation_status, errmsg=emsg)
            call util_allocate_check(allocation_status, emsg, 'SWMMlink_num_elements')
            SWMMlink_num_elements(:) = 0

            !% --- allocate space for storing the number of elements in each node (including non-output nodes)
            !% --- note this MUST be the size of the setting%SWMMinput%N_nodeso that we can later pack indexes
            allocate(SWMMnode_num_elements(setting%SWMMinput%N_node), stat=allocation_status, errmsg=emsg)
            call util_allocate_check(allocation_status, emsg, 'SWMMnode_num_elements')
            SWMMnode_num_elements(:) = 0

            !% --- allocate space for storing the number of faces in each node (including non-output node)
            !% --- this should be one for all faces. It is allocate to be able to create common output routines
            allocate(SWMMnode_num_faces(setting%SWMMinput%N_node), stat=allocation_status, errmsg=emsg)
            call util_allocate_check(allocation_status, emsg, 'SWMMnode_num_faces')
            SWMMnode_num_faces(:) = 0

        end if

        !%-------------------------------------------------------------------
        !% Closing
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
        !if (crashYN) return
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
        !if (crashYN) return
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
        !if (crashYN) return
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
        !if (crashYN) return
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
        !if (crashYN) return
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
        !if (crashYN) return
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
        !if (crashYN) return
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
        !if (crashYN) return
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
        !if (crashYN) return
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
!  OBSOLETE 20220616
    ! subroutine util_allocate_col_elemWDI()
    !     !-----------------------------------------------------------------------------
    !     !
    !     ! Description:
    !     !   the col_elemWDI is a vector of the columns in the faceR arrays
    !     !   that correspond to the enumerated ewdi_... array_index parameter
    !     !
    !     !-----------------------------------------------------------------------------

    !     integer, pointer    :: ncol
    !     integer             :: ii
    !     character(64)       :: subroutine_name = 'util_allocate_col_elemWDI'

    !     !-----------------------------------------------------------------------------
    !     !if (crashYN) return
    !     if (setting%Debug%File%utility_allocate) &
    !         write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    !     !% define the maximum number of columns as
    !     ncol => Ncol_elemWDI

    !     !% allocate an array for storing the column
    !     allocate( col_elemWDI(ncol)[*], stat=allocation_status, errmsg= emsg)
    !     call util_allocate_check (allocation_status, emsg, 'col_elemWDI')

    !     !% this array can be used as a pointer target in defining masks
    !     col_elemWDI(:) = [(ii,ii=1,ncol)]

    !     if (setting%Debug%File%utility_allocate) &
    !     write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    ! end subroutine util_allocate_col_elemWDI
!
!==========================================================================
!==========================================================================
!
    ! subroutine util_allocate_col_elemWDR()
    !% OBSOLETE 20220616
    !     !-----------------------------------------------------------------------------
    !     !
    !     ! Description:
    !     !   the col_elemWDR is a vector of the columns in the faceR arrays
    !     !   that correspond to the enumerated ewdr_... array_index parameter
    !     !
    !     !-----------------------------------------------------------------------------

    !     integer, pointer    :: ncol
    !     integer             :: ii
    !     character(64)       :: subroutine_name = 'util_allocate_col_elemWDI'

    !     !-----------------------------------------------------------------------------
    !     !if (crashYN) return
    !     if (setting%Debug%File%utility_allocate) &
    !         write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    !     !% define the maximum number of columns as
    !     ncol => Ncol_elemWDR

    !     !% allocate an array for storing the column
    !     allocate( col_elemWDR(ncol)[*], stat=allocation_status, errmsg= emsg)
    !     call util_allocate_check (allocation_status, emsg, 'col_elemWDR')

    !     !% this array can be used as a pointer target in defining masks
    !     col_elemWDR(:) = [(ii,ii=1,ncol)]

    !     if (setting%Debug%File%utility_allocate) &
    !     write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    ! end subroutine util_allocate_col_elemWDR
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
        !if (crashYN) return
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
        !if (crashYN) return
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
    ! subroutine util_allocate_col_faceM()
    !     !-----------------------------------------------------------------------------
    !     !
    !     ! Description:
    !     !   the col_faceM is a vector of the columns in the faceR arrays
    !     !   that correspond to the enumerated fM_... array_index parameter
    !     !
    !     !-----------------------------------------------------------------------------

    !     integer, pointer    :: ncol
    !     integer             :: ii
    !     character(64)       :: subroutine_name = 'util_allocate_col_faceM'

    !     !-----------------------------------------------------------------------------
    !     !if (crashYN) return
    !     if (setting%Debug%File%utility_allocate) &
    !         write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    !     !% define the maximum number of columns as
    !     ncol => Ncol_faceM

    !     !% allocate an array for storing the column
    !     allocate( col_faceM(ncol)[*], stat=allocation_status, errmsg= emsg)
    !     call util_allocate_check (allocation_status, emsg, 'col_faceM')

    !     !% this array can be used as a pointer target in defining masks
    !     col_faceM(:) = [(ii,ii=1,ncol)]

    !     if (setting%Debug%File%utility_allocate) &
    !     write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    ! end subroutine util_allocate_col_faceM
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
        !if (crashYN) return
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
        !if (crashYN) return
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
        !if (crashYN) return
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
        !if (crashYN) return
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
!%
!%==========================================================================
!%==========================================================================
!%
!     subroutine util_allocate_col_conmonI()
!         !%------------------------------------------------------------------
!         !% Description
!         !% allocates vector of columns in the conmonI array that correspond
!         !% to the enumerated indexes
!         !%------------------------------------------------------------------
!         !% Declarations
!             integer, pointer :: ncol
!             integer          :: ii
!         !%------------------------------------------------------------------
!         !% define the maximum number of columns as
!         ncol => N_ConMon(this_image())

!         !% allocate an array of column indexes that can be used as targets of pointers
!         allocate( col_conmonI(ncol)[*], stat=allocation_status, errmsg= emsg)
!         call util_allocate_check (allocation_status, emsg, 'col_conmonI')

!         !% this array can be used as a pointer target in defining masks
!         col_conmonI(:) = [(ii,ii=1,ncol)]

!     end subroutine util_allocate_col_conmonI
! !%
! !%==========================================================================
! !%==========================================================================
! !%
!     subroutine util_allocate_col_conmonR()
!         !%------------------------------------------------------------------
!         !% Description
!         !% allocates vector of columns in the conmonR array that correspond
!         !% to the enumerated indexes
!         !%------------------------------------------------------------------
!         !% Declarations
!             integer, pointer :: ncol
!             integer          :: ii
!         !%------------------------------------------------------------------
!         !% define the maximum number of columns as
!         ncol => N_ConMon(this_image())

!         !% allocate an array of column indexes that can be used as targets of pointers
!         allocate( col_conmonR(ncol)[*], stat=allocation_status, errmsg= emsg)
!         call util_allocate_check (allocation_status, emsg, 'col_conmonR')

!         !% this array can be used as a pointer target in defining masks
!         col_conmonR(:) = [(ii,ii=1,ncol)]

!     end subroutine util_allocate_col_conmonR
! !%
! !%==========================================================================
! !%==========================================================================
! !%
!     subroutine util_allocate_col_conmonYN()
!         !%------------------------------------------------------------------
!         !% Description
!         !% allocates vector of columns in the conmonYN array that correspond
!         !% to the enumerated indexes
!         !%------------------------------------------------------------------
!         !% Declarations
!             integer, pointer :: ncol
!             integer          :: ii
!         !%------------------------------------------------------------------
!         !% define the maximum number of columns as
!         ncol => N_ConMon(this_image())

!         !% allocate an array of column indexes that can be used as targets of pointers
!         allocate( col_conmonYN(ncol)[*], stat=allocation_status, errmsg= emsg)
!         call util_allocate_check (allocation_status, emsg, 'col_conmonYN')

!         !% this array can be used as a pointer target in defining masks
!         col_conmonYN(:) = [(ii,ii=1,ncol)]

!     end subroutine util_allocate_col_conmonYN
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_allocate_curves ()
        !-----------------------------------------------------------------------------
        !
        !
        !-----------------------------------------------------------------------------

        character(64)       :: subroutine_name = 'util_allocate_curves'

        !-----------------------------------------------------------------------------
        !if (crashYN) return
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
    subroutine util_allocate_boundary_ghost_elem_array ()
        ! allocate array to copy ghost element data from connected images
        integer :: ii, nrow, n_elemB
        integer, pointer :: ncol, Nfaces, Nelems
        character(64) :: subroutine_name = 'util_allocate_boundary_ghost_elem_array'

        !-----------------------------------------------------------------------------
        !if (crashYN) return
        if (setting%Debug%File%utility_allocate) &
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"


        Nelems => N_elem(this_image())
        nrow = count(elemYN(1:Nelems,eYN_isBoundary_up)) + count(elemYN(1:Nelems,eYN_isBoundary_dn))
        ncol => Ncol_elemBGR

        allocate(elemGR(nrow, ncol), stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'elemGR')
        elemGR(:,:) = nullvalueR

        allocate(elemB[*], stat=allocation_status, errmsg=emsg)

        allocate(elemB%R(nrow, ncol), stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'elemB%R')
        elemB%R(:,:) = nullvalueR

        ncol => Ncol_elemBGI
        allocate(elemB%I(nrow, ncol), stat=allocation_status, errmsg=emsg)
        call util_allocate_check(allocation_status, emsg, 'elemB%I')
        elemB%I(:,:) = nullvalueI

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine util_allocate_boundary_ghost_elem_array
!
!==========================================================================
!==========================================================================
!
    subroutine util_allocate_curve_entries (curve_idx, num_entries)
        !%-----------------------------------------------------------------
        !% Description:
        !% allocates the curve table of curve_indx for the number of values
        !% expected (num_entries)
        !%-----------------------------------------------------------------
        !% Declarations:
            integer, intent(in) :: curve_idx, num_entries
            character(64)       :: subroutine_name = 'util_allocate_curve_entries'
        !%-----------------------------------------------------------------
        !% Preliminaries
            !if (crashYN) return
            if (setting%Debug%File%utility_allocate) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%-----------------------------------------------------------------
                
        !% allocate the value array of curve of the given curve_idx
        allocate( curve(curve_idx)%ValueArray(num_entries,Ncol_curve), stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'curve entries')

        curve(curve_idx)%ValueArray = nullvalueR

        if (setting%Debug%File%utility_allocate) &
        write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

    end subroutine util_allocate_curve_entries
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine util_allocate_temporary_arrays ()
        !%-----------------------------------------------------------------
        !% Description:
        !% allocates miscellaneous temporary arrays that don't fit standard
        !% sizes
        !%-----------------------------------------------------------------
        !% Preliminaries
            !if (crashYN) return
        !%-----------------------------------------------------------------

        !% array of integers of same size as number of BCup face
        if (size(BC%P%BCup) > 0) then
            allocate( temp_BCupI(size(BC%P%BCup),N_tempBCupI), stat=allocation_status, errmsg= emsg)
            call util_allocate_check (allocation_status, emsg, 'temp_BCupI')

            allocate( temp_BCupR(size(BC%P%BCup),N_tempBCupR), stat=allocation_status, errmsg= emsg)
            call util_allocate_check (allocation_status, emsg, 'temp_BCupR')
        else
            !% if no BCup faces, then use size 1 to prevent seg fault
            allocate( temp_BCupI(1,N_tempBCupI), stat=allocation_status, errmsg= emsg)
            call util_allocate_check (allocation_status, emsg, 'temp_BCupI')

            allocate( temp_BCupR(1,N_tempBCupR), stat=allocation_status, errmsg= emsg)
            call util_allocate_check (allocation_status, emsg, 'temp_BCupR')
        end if

        temp_BCupI(:,:) = nullvalueI
        temp_BCupR(:,:) = nullvalueR

    end subroutine util_allocate_temporary_arrays


    subroutine util_allocate_output_profiles()

        allocate(output_profile_ids(max_profiles_N,max_links_profile_N),stat=allocation_status, errmsg= emsg)
        call util_allocate_check (allocation_status, emsg, 'output_profile_ids')


    end subroutine util_allocate_output_profiles





!%    
!%==========================================================================
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
            !if (crashYN) return
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
!%    
!%==========================================================================
!% END OF MODULE
!%==========================================================================
!%
end module utility_allocate
