module define_globals
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Global parameters and variables for SWMM5+
    !%
    !% Methods:
    !% These are global parameters and variables that are NOT associated with
    !% the 2D data storage arrays or the setting structure.
    !%
    !% The globals for data storage are found in modules array_index and data_keys
    !% The globals for setting are found in module setting_definition
    !%
    !% NOTE: don't use setting%... structure in define_globals to prevent linking problems
    !%==========================================================================
    use define_types
    use define_api_keys

    implicit none

    public

    !% reverseKeys is an array to find the string name of a key number in define_keys
    !% used for debugging
    character(len=32), allocatable :: reverseKey(:)
    
    integer(kind=8) :: irecCount = 0

    character (len=256) ::  outstring[*]

    real(8) :: Qcumulative = 0.d0
    real(8) :: VolumeErrorCumulative = 0.d0
!%
!%===================================================================================
!% VARIABLES
!%===================================================================================
    !% Developer note -- most variables are in setting%..., the only variables that should
    !% be in globals are ones that are internal to the code (no user setting) and always
    !% referred to with a short name
    
    !% --- internal flag for crash condition
    integer :: crashI = 0 

    !% set to true if present time step is a spin-up time step
    logical :: inSpinUpYN

!%===================================================================================
!% ARRAYS
!%===================================================================================

    !% --- Number of maximum branches for a junction
    !% note that these must always be matched with the same
    !% number of up and down nodes. 
    !% NOTE: element branch indexes in (e.g.) elemR arrays
    !%  are always sequential after the Junction Main index, 
    !% i.e. JBidx = JMidx + n, where odd n are upstream branches and even n are
    !% dnstream. IMPORTANT: this convention is different from column numbers used for
    !% node%I(:,ni_Mlink_u1), etc.
    !% The ADDBRANCH tag can be searched on to find places where branch changes are
    !% important
    integer, parameter :: max_up_branch_per_node = 5   !% ADDBRANCH 
    integer, parameter :: max_dn_branch_per_node = 5   !% ADDBRANCH
    integer, parameter :: max_branch_per_node = 10     !% ADDBRANCH

    !% --- branchsign is lobal is used for junction branches (JB)
    !%    for upstream (+1.d0) and downstream (-1.d0)
    real(8) :: branchsign(max_branch_per_node)

    !% --- Main Arrays - Allocated in utility_allocate.f08
    !%     links are the building blocks from SWMM link-node formulation
    type(LinkArray), target :: link

    !% --- nodes are the building blocks from the SWMM link-node formulation
    type(NodeArray), target :: node

    !% --- transect data (irregular cross-sections)
    integer, allocatable, target :: transectI(:,:)
    real(8), allocatable, target :: transectR(:,:)
    real(8), allocatable, target :: transectTableDepthR(:,:,:) ! transect table by depth
    real(8), allocatable, target :: transectTableAreaR(:,:,:)  ! transect table by area
    character(len=64), allocatable, target :: transectID(:)

    !% --- section factor lookup data (for limited number of elements)
    integer, allocatable, target :: uniformTableI(:,:)
    real(8), allocatable, target :: uniformTableR(:,:)
    real(8), allocatable, target :: uniformTableDataR(:,:,:) 

    !% --- boundary elements array
    type(BoundaryElemArray), allocatable :: elemB[:]

    !% --- Boundary conditions
    type(BCArray), target :: BC

    !% --- Curve typesnode_overflow                  
    type(curveType), dimension(:), allocatable, target :: curve

    !%  --- columns of element and face arrays
    integer, allocatable, target :: col_elemI(:)[:]                                !% columns of elemI array
    integer, allocatable, target :: col_elemP(:)[:],       npack_elemP(:)[:]       !% columns and number of packs for elemP array
    !integer, allocatable, target :: col_elemPGalltm(:)[:], npack_elemPGalltm(:)[:] !% columns and number of packs for elemPG array for all tm
    !integer, allocatable, target :: col_elemPGac(:)[:],    npack_elemPGac(:)[:]    !% columns and number of packs for elemPG array for ac tm
    integer, allocatable, target :: col_elemPGetm(:)[:],   npack_elemPGetm(:)[:]   !% columns and number of packs for elemPG array for etm
    integer, allocatable, target :: col_elemR(:)[:]                                !% columns of elemR array
    integer, allocatable, target :: col_elemSI(:)[:]                               !% columns of elemSI array
    integer, allocatable, target :: col_elemSR(:)[:]                               !% columns of elemSR array
    integer, allocatable, target :: col_elemSGR(:)[:]                              !% columns of elemSGR array
    integer, allocatable, target :: col_elemYN(:)[:]                               !% columns of elemYN array
    integer, allocatable, target :: col_faceI(:)[:]                                !% columns of faceI array
    integer, allocatable, target :: col_faceM(:)[:]                                !% columns of faceM array
    integer, allocatable, target :: col_faceP(:)[:],       npack_faceP(:)[:]       !% columns and number of packs for faceP array
    integer, allocatable, target :: col_facePS(:)[:],      npack_facePS(:)[:]      !% columns and number of packs for facePS array
    integer, allocatable, target :: col_faceR(:)[:]                                !% columns of faceR array
    integer, allocatable, target :: col_faceYN(:)[:]                               !% columns of faceYN array

    !%  --- number of dummy row in elem arrays (do not change)
    integer, parameter :: N_dummy_elem = 1

    !%  --- vector of number of elements and faces across images
    integer, allocatable, target :: N_elem(:)
    integer, allocatable, target :: N_face(:)
    integer, allocatable, target :: N_unique_face(:)

    !% --- output elements on each image
    integer, allocatable, target :: N_OutElem(:)[:]

    !% --- output faces on each image
    integer, allocatable, target :: N_OutFace(:)[:]

    !% --- number of control action and monitoring points in system (identical on each image)
    integer,  target :: N_ActionPoint
    integer,  target :: N_MonitorPoint

    !% --- number of profiles for output
    integer :: N_Profiles
    integer :: N_links_in_profile           !% number of links (maximum) in any profile
    integer :: N_items_in_profile           !% number of items (maximum) in any profile (i.e., link,node,link...)
    integer :: N_links_on_profile_line = 10 !% number of links on any line of profile definition
    integer :: N_char_of_profile_name  = 64 !% number of characters allowed in a profile name

    !%--- number of culverts on an image
    !integer, allocatable, target :: N_culvert(:)

    !%  --- elements in coarray
    real(8), allocatable, target :: elemR(:,:)[:]       !% coarray for elements
    integer, allocatable, target :: elemI(:,:)[:]       !% coarray for element Interger
    logical, allocatable, target :: elemYN(:,:)[:]      !% coarray for element logical
    integer, allocatable, target :: elemP(:,:)[:]       !% coarray for element pack array
    real(8), allocatable, target :: elemOutR(:,:,:)[:]  !% coarray for packed, multi-level output storage (index,type,level)

    integer, allocatable, target :: elemPGetm(:,:)[:]   !% coarray for element pack geometry array
    integer, allocatable, target :: elemSI(:,:)[:]      !% coarray for special element Integer
    real(8), allocatable, target :: elemSR(:,:)[:]      !% coarray for special elemen Real
    real(8), allocatable, target :: elemSGR(:,:)[:]     !% coarray for special element geometry Real
    real(8), allocatable, target :: elemGR(:,:)         !% array to copy required ghost element data (not a coarray)

    logical, allocatable, target :: elemIsNan(:,:)    !% used for checking values for NaN during runtime

    integer, allocatable, target :: actionI(:,:)        !% action data for EPA-SWMM controls
    real(8), allocatable, target :: actionR(:,:)        !% actionR MIGHT NOT BE NEEDED

    integer, allocatable, target :: monitorI(:,:)       !% monitor data for EPA-SWMM controls
    real(8), allocatable, target :: monitorR(:,:)  

    real(8), allocatable, target :: monitorPassR(:)[:]  !% monitor data to be passed between coarrays

    !% --- Culvert parameters
    integer, parameter :: NculvertTypes = 57
    integer, parameter :: NculvertParams = 6
    real(8), dimension(NculvertTypes,NculvertParams) :: culvertValue

    !%  --- faces in coarray
    real(8), allocatable, target :: faceR(:,:)[:]       !% coarray for faces real data
    integer, allocatable, target :: faceI(:,:)[:]       !% coarray for faces integer data
    logical, allocatable, target :: faceYN(:,:)[:]      !% coarray for faces logical data
    integer, allocatable, target :: faceP(:,:)[:]       !% coarray for faces pack array
    integer, allocatable, target :: facePS(:,:)[:]      !% coarray for shared faces pack array
    real(8), allocatable, target :: faceOutR(:,:,:)[:]  !% coarray for packed, multi-level output storage (index,type,level)

    logical, allocatable, target :: faceIsNan(:,:)    !% used for checking values for NaN during runtime

    !% --- subcatchments -- NOT coarray
    real(8), allocatable, target :: subcatchR(:,:)[:]   !% subcatchment real data
    integer, allocatable, target :: subcatchI(:,:)      !% subcatchment integer data
    logical, allocatable, target :: subcatchYN(:,:)     !% subcatchment logical data

    !% --- BIPquick Arrays - (De)Allocated in BIPquick.f90
    integer, allocatable :: B_nodeI(:,:)
    real(8), allocatable :: B_nodeR(:,:)
    integer, allocatable :: B_roots(:)
    real(8), allocatable :: weight_range(:,:)
    logical, allocatable :: totalweight_visited_nodes(:)
    logical, allocatable :: partitioned_nodes(:)
    logical, allocatable :: partitioned_links(:)
    logical, allocatable :: accounted_for_links(:)
    integer, allocatable :: phantom_link_tracker(:)

    !% --- Partitioning module Allocatables - Allocated and Deallocated in execute_partitioning.f08
    integer, allocatable :: adjacent_links(:)
    integer, allocatable :: elem_per_image(:)
    logical, allocatable :: image_full(:)

    !% --- Storage for entrapped air in links
    integer, allocatable, target :: sc_link_Idx(:,:)        !% link indexes of super conduits
    integer, allocatable, target :: links_per_sc(:)         !% save the total number of links in a super conduit
    integer, allocatable, target :: conduitElemMapsI(:,:,:) !% 3-d array to store element and face maps to a corresponding conduits
    integer, allocatable, target :: AirI(:,:,:)             !% 3-d array to store integer of entrapped air pockets
    real(8), allocatable, target :: AirR(:,:,:)             !% 3-d array to store reals of entrapped air pockets 
    logical, allocatable, target :: AirYN(:,:,:)            !% 3-d array to store logicals of entrapped pockets   

    !% --- link and node output_idx
    integer, allocatable :: link_output_idx(:)
    integer, allocatable :: node_output_idx(:)

    integer, allocatable, target           :: output_static_types_Link(:)
    integer, allocatable, target           :: output_static_typeProcessing_Link(:)
    integer, allocatable, target           :: output_static_typeMultiplyByBarrels_Link(:)
    character(len=64), allocatable, target :: output_static_typeNames_Link(:)
    character(len=16), allocatable         :: output_static_typeUnits_Link(:)
    real(8), allocatable                   :: output_static_link(:,:)[:]

    integer, allocatable, target           :: output_static_types_Node(:)
    integer, allocatable, target           :: output_static_typeProcessing_Node(:)
    integer, allocatable, target           :: output_static_typeMultiplyByBarrels_Node(:)
    character(len=64), allocatable, target :: output_static_typeNames_Node(:)
    character(len=16), allocatable         :: output_static_typeUnits_Node(:)
    real(8), allocatable                   :: output_static_node(:,:)[:]

    !% --- element output types
    integer, allocatable, target           :: output_types_elemR(:)
    integer, allocatable, target           :: output_typeProcessing_elemR(:)
    integer, allocatable, target           :: output_typeMultiplyByBarrels_elemR(:)
    character(len=64), allocatable, target :: output_typeNames_elemR(:)
    character(len=16), allocatable         :: output_typeUnits_elemR(:)
    character(len=64), allocatable, target :: output_typeNames_withTime_elemR(:)
    character(len=15), allocatable, target :: output_typeUnits_withTime_elemR(:)

    integer, allocatable, target           :: output_static_types_elemR(:)
    integer, allocatable, target           :: output_static_typeProcessing_elemR(:)
    integer, allocatable, target           :: output_static_typeMultiplyByBarrels_elemR(:)
    character(len=64), allocatable, target :: output_static_typeNames_elemR(:)
    character(len=16), allocatable         :: output_static_typeUnits_elemR(:)
    real(8), allocatable                   :: output_static_elem(:,:)[:]

    !% --- face output types
    integer, allocatable, target           :: output_types_faceR(:)
    integer, allocatable, target           :: output_typeProcessing_faceR(:)
    integer, allocatable, target           :: output_typeMultiplyByBarrels_faceR(:)
    character(len=64), allocatable, target :: output_typeNames_faceR(:)
    character(len=16), allocatable         :: output_typeUnits_faceR(:)
    character(len=64), allocatable, target :: output_typeNames_withTime_faceR(:)
    character(len=16), allocatable         :: output_typeUnits_withTime_faceR(:)

    !% --- output times
    real(8), allocatable :: output_times(:)

    !% --- character length used in HDF5 output
    integer, parameter :: stringLength_HDF5 = 150

    !% --- Profile_Outputs
    integer, allocatable, target :: output_profile_ids(:,:)
    integer, allocatable, target :: output_profile_link_idx(:,:)
    integer, allocatable, target :: Number_of_links_per_profile(:)
    character(len=stringLength_HDF5), allocatable, target :: output_profile_names(:)
    character(len=stringLength_HDF5), allocatable, target :: output_profile_link_names(:,:)
    character(len=stringLength_HDF5), allocatable, target :: output_profile_node_names(:,:)

    !% --- filenames for output binaries
    character(len=256), allocatable, target :: output_binary_filenames(:)
    character(len=256), allocatable, target :: output_binary_filenames_all(:)

    !% --- storage for all real data on all output elements for limited time levels
    !%     (outelement, type, time-level) real data -- not coarray
    real(8), allocatable, target :: OutElemDataR(:,:,:)

    !% --- storage for fixed integer data needed in output elements -- not coarray
    !%     (outelement, dataindex)
    integer, allocatable, target :: OutElemFixedI(:,:)

    !% --- storage for all real data on all output faces for limited time levels
    !%     (outelement, type, time-level) real data -- not coarray
    real(8), allocatable, target :: OutFaceDataR(:,:,:)

    !% --- storage for fixed integer data needed in output faces -- not coarray
    !%     (outface, dataindex)
    integer, allocatable, target :: OutFaceFixedI(:,:)

    integer, allocatable, target :: thisElementOut(:), thisFaceOut(:)

    integer, allocatable :: SWMMlink_num_elements(:) !% number of elements in each output link
    integer, allocatable :: SWMMnode_num_elements(:) !% number of elements in each output node
    integer, allocatable :: SWMMnode_num_faces(:)    !% number of faces in each output node

    !% --- maximum string length in link%Names or node%Names
    !%     value set in util_allocate_linknode()
    integer :: max_names_string_length = 0

    !% --- Temporary arrays that don't fit in any of the standard array structures
    integer, allocatable, target :: temp_BCupI(:,:)  !% number of BCup faces.
    real(8), allocatable, target :: temp_BCupR(:,:)
    integer :: N_tempBCupI = 1
    integer :: N_tempBCupR = 4
    
    !% --- Profiling Timer
    type(wall_clk) :: timer

    !% --- profiling storage
    real(8), allocatable :: profiler_data(:,:)
    character (len=64), allocatable :: profiler_procedure_name(:)
    integer, allocatable :: profiler_procedure_level(:)

    !% --- counter for transects
    integer :: lastTransectIdx = 0

    !% --- dummy index for elements/faces that do not exist
    !%     value set in init_network_set_dummy_elem()
    integer               :: dummyIdx
    integer, dimension(1) :: dummyUnitArray

    !% ---- dummy targets
    integer, target :: dummyTargetI(1) 
    real(8), target :: dummyTargetR(1)

!% ===================================================================================
!% CONSTANTS
!% ===================================================================================

    !% --- note that nullvalueI > 0 is required. A large value MUST be used
    !%     as this will be a limit on the largest network size
    integer, parameter :: nullvalueI = 998877  
    real(8), parameter :: nullvalueR = 9.98877d16
    logical, parameter :: nullvalueL = .false.
    real(8), parameter :: zeroR = 0.d0
    real(8), parameter :: oneR = 1.d0
    real(8), parameter :: twoR = 2.d0
    real(8), parameter :: threeR = 3.d0
    real(8), parameter :: fourR = 4.d0
    real(8), parameter :: fiveR = 5.d0
    real(8), parameter :: sixR = 6.d0
    real(8), parameter :: sevenR = 7.d0
    real(8), parameter :: eightR = 8.d0
    real(8), parameter :: nineR = 9.d0
    real(8), parameter :: tenR = 10.d0
    real(8), parameter :: sixteenR = 16.d0
    real(8), parameter :: twentyR = 20.d0
    real(8), parameter :: twentyfourR = 24.d0
    real(8), parameter :: sixtyR = 60.d0
    real(8), parameter :: fourtyeightR = 48.d0
    real(8), parameter :: ninetyR = 90.d0
    real(8), parameter :: onehundredR = 100.d0
    real(8), parameter :: fivehundredR = 500.d0
    real(8), parameter :: onethousandR = 1000.d0
    real(8), parameter :: pi = 4.d0*datan(1.d0)

    real(8), parameter :: oneOneThousandthR = oneR / onethousandR
    real(8), parameter :: oneOneHundredthR = oneR / onehundredR
    real(8), parameter :: onetenthR = oneR / tenR
    real(8), parameter :: oneninthR = oneR / nineR
    real(8), parameter :: oneeighthR = oneR / eightR
    real(8), parameter :: onesixthR = oneR / sixR
    real(8), parameter :: onefifthR = oneR / fiveR
    real(8), parameter :: onefourthR = oneR / fourR
    real(8), parameter :: onethirdR = oneR / threeR
    real(8), parameter :: onehalfR = oneR / twoR
    real(8), parameter :: twothirdR = twoR / threeR
    real(8), parameter :: threefourthR = threeR / fourR
    real(8), parameter :: threehalfR = threeR / twoR
    real(8), parameter :: fourthirdsR = fourR / threeR

    real(8), parameter :: seconds_per_minute = 60.d0
    real(8), parameter :: seconds_per_hour = 3600.d0
    real(8), parameter :: seconds_per_day  = 86400.d0

    real(8), parameter :: meters_per_ft = 0.3048d0

    integer, parameter :: dummyI = 1  !% used when dummy argument isneeded
    integer, parameter :: zeroI  = 0
    integer, parameter :: oneI   = 1
    integer, parameter :: twoI   = 2
    integer, parameter :: threeI = 3
    integer, parameter :: fourI  = 4
    integer, parameter :: fiveI  = 5
    integer, parameter :: sixI   = 6

    !% --- Number of objects
    integer :: N_transect            ! # of irregular cross-section transects for elements
    integer :: N_transect_depth_items
    integer :: N_transect_area_items
    integer :: N_uniformTableData_items = 51  !% 51 is consistent with tables of EPA-SWMM
    integer :: N_uniformTable_locations
    integer :: N_link
    integer :: N_conduit
    integer :: N_super_conduit
    integer :: N_node
    integer :: N_headBCnode !% number of head BC defined on nodes
    integer :: N_flowBCnode !% number of flow BC defined on nodes (may be pushed to links)
    integer :: N_flowBCLink !% number of flow BC defined on links
    integer :: N_flowBCelem !% number of flow BC on elements (includes link-elem and node-elem)
    integer :: N_nBCup
    integer :: N_nBCdn
    integer :: N_nBClat
    integer :: N_nJm
    integer :: N_nStorage
    integer :: N_nJ2
    integer :: N_nJ1
    integer :: N_diag
    integer :: N_ac
    integer :: N_etm
    integer :: N_link_output
    integer :: N_node_output
    integer :: N_monitor_types = 7 !% # of data types transferred from monitorR in monitorPassR
    integer :: N_Total_Curves !% sum of swmm input curves and additional storage curves
    integer :: N_subcatch_runon
    integer, target :: N_OutTypeElem
    integer, target :: N_Out_static_TypeElem
    integer, target :: N_Out_static_TypeLink
    integer, target :: N_Out_static_TypeNode
    integer, target :: N_OutTypeFace

    !% --- Coarray variables
    integer :: max_caf_elem_N    ! size of all elem array in coarray
    integer :: max_caf_face_N    ! size of all face array in coarray

    !% --- assign status parameters for nodes
    integer, parameter :: nUnassigned = 998877
    integer, parameter :: nAssigned   = 1
    integer, parameter :: nDeferred   = -1

    !% --- assign status parameters for links
    integer, parameter :: lUnassigned = 998877
    integer, parameter :: lAssigned   = 1
    integer, parameter :: lDeferred   = -1

    !% --- Constants for Junction
    integer, parameter :: J_elem_add = max_branch_per_node+1 ! Supplement elements for junction 
    integer, parameter :: J_face_add = max_branch_per_node   ! Supplement faces for junction

    !% --- default number of elements for different node types
    integer, parameter :: N_elem_nJ2      = 0 ! 2-link nodes are assigned to a single face
    integer, parameter :: N_elem_nJm      = max_branch_per_node+1 ! M-link nodes are assigned a maximum of 7 elements !%brh20211219
    integer, parameter :: N_elem_nStorage = 1 ! Storage nodes are assigned to 1 element
    integer, parameter :: N_elem_nBCdn    = 0 ! Downstream BC nodes are assigned to 0 element (only a face)
    integer, parameter :: N_elem_nBCup    = 0 ! Upstream BC nodes are assigned to 0 element (only a face)
    integer, parameter :: N_elem_nJ1      = 0 ! Upstream non-BC nodes are assigned to 0 elements (only a face)

    !% --- default for edge and non-edge node
    !% --- HACK: EdgeNode and NonEdgeNode are identified as 1 and 0 in the partition code
    !%           storing in the globals as parameter for now
    integer, parameter :: nonEdgeNode  = 0    ! Upstream BC nodes are assigned to 1 element
    integer, parameter :: EdgeNode     = 1    ! Edge node of a partition

    ! --- default for number of functional storage
    integer :: N_FunctionalStorage = 0

    !% --- maximum cfl across the network. Used for reporting 
    real(8) :: cfl_max

    !% --- datetime related variables
    !%     generally use an "epoch time" starting from 01/01/1900
    integer, parameter :: datedelta = 693594 !% --- days from 01/01/0000 to 01/01/1900
    integer, parameter :: secsperday = 86400
    integer :: dayspermonth(12,2) = &
        reshape((/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, & ! normal years
        31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/), (/12,2/)) ! leap years
!%==========================================================================
!% End of module
!%==========================================================================
end module define_globals
