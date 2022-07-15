! module define_globals
!
! These are global parameters and variables that are not associated with
! the 2D data storage arrays or the setting structure.
!
! The globals for data storage are found in modules array_index and data_keys
! The globals for setting are found in module setting_definition
!
!==========================================================================
module define_globals

    use define_types
    use define_api_keys

    implicit none

    public

    !% reverseKeys is an array to find the string name of a key number in define_keys
    !% used for debugging
    character(len=32), allocatable :: reverseKey(:)
    

    ! integer :: iet(10) = (/45,46,48,47,49,115,116,  113, 114,50/)
    ! integer :: ift(3) = (/47,48,            49/)
      
    !integer :: iet(2) = (/193,1 /)
    !integer :: ift(3) = (/205,1, 2/)

   !  integer :: iet(2) = (/192,216/)
   !  integer :: ift(3) = (/203,204,217/)

   !integer :: iet(10) = (/45, 46, 47, 48, 49, 50, 113, 114, 115, 116/)
   ! integer :: ift(3) = (/205,1,2/)

    !integer :: ift(4) = (/13,1,2,11/)

    integer :: ift(6) = (/15,1,2,11,12,13/)

    !integer :: iet(7) = (/13,15,24,2,1,3,12 /) !% orifice
    !integer :: iet(7) = (/13,17,25,4,1,3,12 /) !% weir


    integer :: iet(9) = (/15,17,26,2,1,3,12,13,14 /) !% orifice
    !integer :: iet(7) = (/13,17,25,4,1,3,12 /) !% weir

    integer(kind=8) :: irecCount = 0

    character (len=256) ::  outstring[*]


!% ===================================================================================
!% VARIABLES
!% ===================================================================================
!% Developer note -- most variables are in setting%..., the only variables that should
!% be in globals are ones that are internal to the code (no user setting) and always
!% referred to with a short name


    !logical :: crashYN = .false. !% error condition
    integer :: crashI = 0 

    !% set to true if present time step is a spin-up time step
    logical :: inSpinUpYN

!% ===================================================================================
!% ARRAYS
!% ===================================================================================


    !% Number of maximum branches for a junction
    !% note that these must always be matched with the same
    !% number of up and down nodes. 
    integer, parameter :: max_up_branch_per_node = 5   !% ADDBRANCH
    integer, parameter :: max_dn_branch_per_node = 5   !% ADDBRANCH
    integer, parameter :: max_branch_per_node = 10      !% ADDBRANCH

    real(8) :: branchsign(max_branch_per_node)

    !% Main Arrays - Allocated in utility_allocate.f08
    !%  links are the building blocks from SWMM link-node formulation
    type(LinkArray), target :: link

    !%  nodes are the building blocks from the SWMM link-node formulation
    type(NodeArray), target :: node

    integer, allocatable, target :: transectI(:,:)
    real(8), allocatable, target :: transectR(:,:)
    real(8), allocatable, target :: transectTableDepthR(:,:,:) ! transect table by depth
    real(8), allocatable, target :: transectTableAreaR(:,:,:)  ! transect table by area
    character(len=64), allocatable, target :: transectID(:)

    !% boundary elements array
    type(BoundaryElemArray), allocatable :: elemB[:]

    !% Boundary conditions
    type(BCArray), target :: BC

    !!% Curve typesnode_overflow                    // 28
    type(curveType), dimension(:), allocatable, target :: curve

    !%  columns of element and face arrays
    integer, allocatable, target :: col_elemI(:)[:]                                !% columns of elemI array
    integer, allocatable, target :: col_elemP(:)[:],       npack_elemP(:)[:]       !% columns and number of packs for elemP array
    integer, allocatable, target :: col_elemPGalltm(:)[:], npack_elemPGalltm(:)[:] !% columns and number of packs for elemPG array for all tm
    integer, allocatable, target :: col_elemPGac(:)[:],    npack_elemPGac(:)[:]    !% columns and number of packs for elemPG array for ac tm
    integer, allocatable, target :: col_elemPGetm(:)[:],   npack_elemPGetm(:)[:]   !% columns and number of packs for elemPG array for etm
    integer, allocatable, target :: col_elemR(:)[:]                                !% columns of elemR array
    integer, allocatable, target :: col_elemSI(:)[:]                               !% columns of elemSI array
    integer, allocatable, target :: col_elemSR(:)[:]                               !% columns of elemSR array
    integer, allocatable, target :: col_elemSGR(:)[:]                              !% columns of elemSGR array
    !integer, allocatable, target :: col_elemWDI(:)[:]                              !% columns of elemWDI array
    !integer, allocatable, target :: col_elemWDR(:)[:]                              !% columns of elemWDR array
    integer, allocatable, target :: col_elemYN(:)[:]                               !% columns of elemYN array
    integer, allocatable, target :: col_faceI(:)[:]                                !% columns of faceI array
    integer, allocatable, target :: col_faceM(:)[:]                                !% columns of faceM array
    integer, allocatable, target :: col_faceP(:)[:],       npack_faceP(:)[:]       !% columns and number of packs for faceP array
    integer, allocatable, target :: col_facePS(:)[:],      npack_facePS(:)[:]      !% columns and number of packs for facePS array
    integer, allocatable, target :: col_faceR(:)[:]                                !% columns of faceR array
    integer, allocatable, target :: col_faceYN(:)[:]                               !% columns of faceYN array
    ! integer, allocatable, target :: col_conmonI(:)[:]                              !% columns of conmomI array
    ! integer, allocatable, target :: col_conmonR(:)[:]                              !% columns of conmonR array
    ! integer, allocatable, target :: col_conmonYN(:)[:]                             !% columns of conmonYN array

    !%  number of dummy row in elem arrays (do not change)
    integer, parameter :: N_dummy_elem = 1

    !%  vector of number of elements and faces across images
    integer, allocatable, target :: N_elem(:)
    integer, allocatable, target :: N_face(:)
    integer, allocatable, target :: N_unique_face(:)
    !% output elements on each image
    integer, allocatable, target :: N_OutElem(:)[:]
    !% output nodes on each image
    integer, allocatable, target :: N_OutFace(:)[:]
    !% number of control action and monitoring points in system (identical on each image)
    integer,  target :: N_ActionPoint
    integer,  target :: N_MonitorPoint

    !% --- packed array of columns from elemR used in conmonR
    integer, allocatable, target :: cmR_eR_col(:)

    !%  elems in coarray
    real(8), allocatable, target :: elemR(:,:)[:]       !% coarray for elements
    integer, allocatable, target :: elemI(:,:)[:]       !% coarray for element Interger
    logical, allocatable, target :: elemYN(:,:)[:]      !% coarray for element logical
    integer, allocatable, target :: elemP(:,:)[:]       !% coarray for element pack array
    real(8), allocatable, target :: elemOutR(:,:,:)[:]  !% coarray for packed, multi-level output storage (index,type,level)

    integer, allocatable, target :: elemPGalltm(:,:)[:] !% coarray for element pack geometry array
    integer, allocatable, target :: elemPGac(:,:)[:]    !% coarray for element pack geometry array
    integer, allocatable, target :: elemPGetm(:,:)[:]   !% coarray for element pack geometry array
    integer, allocatable, target :: elemSI(:,:)[:]      !% coarray for special element Integer
    real(8), allocatable, target :: elemSR(:,:)[:]      !% coarray for special elemen Real
    real(8), allocatable, target :: elemSGR(:,:)[:]     !% coarray for special element geometry Real
    real(8), allocatable, target :: elemGR(:,:)         !% array to copy required ghost element data (not a coarray)

    ! integer, allocatable, target :: conmonI(:,:)        !% element data for control/monitoring points (not a coarray)
    ! real(8), allocatable, target :: conmonR(:,:)        !% real data for control/monitoring points (not a coarray)
    ! logical, allocatable, target :: conmonYN(:,:)       !% logical data for control/monitoring points (not a coarray)

    integer, allocatable, target :: actionI(:,:)        !% action data for EPA-SWMM controls
    real(8), allocatable, target :: actionR(:,:)        !% actionR MIGHT NOT BE NEEDED
    !logical, allocatable, target :: actionYN(:,:)

    integer, allocatable, target :: monitorI(:,:)       !% monitor data for EPA-SWMM controls
    real(8), allocatable, target :: monitorR(:,:)
    !logical, allocatable, target :: monitorYN(:,:)   

    real(8), allocatable, target :: monitorPassR(:)[:]  !% monitor data to be passed between coarrays

    !%  faces in coarray
    real(8), allocatable, target :: faceR(:,:)[:]       !% coarray for faces real data
    integer, allocatable, target :: faceI(:,:)[:]       !% coarray for faces integer data
    logical, allocatable, target :: faceYN(:,:)[:]      !% coarray for faces logical data
    integer, allocatable, target :: faceP(:,:)[:]       !% coarray for faces pack array
    integer, allocatable, target :: facePS(:,:)[:]      !% coarray for shared faces pack array
    !logical, allocatable, target :: faceM(:,:)[:]       !% coarray for faces mask array
    real(8), allocatable, target :: faceOutR(:,:,:)[:]  !% coarray for packed, multi-level output storage (index,type,level)

    !% subcatchments -- NOT coarray
    real(8), allocatable, target :: subcatchR(:,:)      !% subcatchment real data
    integer, allocatable, target :: subcatchI(:,:)      !% subcatchment integer data
    logical, allocatable, target :: subcatchYN(:,:)     !% subcatchment logical data

    !% BIPquick Arrays - (De)Allocated in BIPquick.f08
    integer, allocatable :: B_nodeI(:,:)
    real(8), allocatable :: B_nodeR(:,:)
    integer, allocatable :: B_roots(:)
    real(8), allocatable :: weight_range(:,:)
    logical, allocatable :: totalweight_visited_nodes(:)
    logical, allocatable :: partitioned_nodes(:)
    logical, allocatable :: partitioned_links(:)
    logical, allocatable :: accounted_for_links(:)
    integer, allocatable :: phantom_link_tracker(:)

    !% Partitioning module Allocatables - Allocated and Deallocated in execute_partitioning.f08
    integer, allocatable :: adjacent_links(:)
    integer, allocatable :: elem_per_image(:)
    logical, allocatable :: image_full(:)

    !%link and node output_idx
    integer, allocatable :: link_output_idx(:)
    integer, allocatable :: node_output_idx(:)

    !% element output types
    integer, allocatable, target           :: output_types_elemR(:)
    integer, allocatable, target           :: output_typeProcessing_elemR(:)
    character(len=64), allocatable, target :: output_typeNames_elemR(:)
    character(len=16), allocatable         :: output_typeUnits_elemR(:)
    character(len=64), allocatable, target :: output_typeNames_withTime_elemR(:)
    character(len=15), allocatable, target :: output_typeUnits_withTime_elemR(:)


    !% face output types
    integer, allocatable, target           :: output_types_faceR(:)
    integer, allocatable, target           :: output_typeProcessing_faceR(:)
    character(len=64), allocatable, target :: output_typeNames_faceR(:)
    character(len=16), allocatable         :: output_typeUnits_faceR(:)
    character(len=64), allocatable, target :: output_typeNames_withTime_faceR(:)
    character(len=16), allocatable         :: output_typeUnits_withTime_faceR(:)

    !% output times
    real(8), allocatable :: output_times(:)
    ! filenames for output binaries
    character(len=256), allocatable, target :: output_binary_filenames(:)
    character(len=256), allocatable, target :: output_binary_filenames_all(:)

    !% storage of global element index for output -- not coarray
    !brh rm integer, allocatable, target :: OutElemGidx(:)

    !% storage for all real data on all output elements for limited time levels
    !% (outelement, type, time-level) real data -- not coarray
    real(8), allocatable, target :: OutElemDataR(:,:,:)
    !% storage for fixed integer data needed in output elements -- not coarray
    !% (outelement, dataindex)
    integer, allocatable, target :: OutElemFixedI(:,:)

    !% storage of global face index for output -- not coarray
    !brh rm integer, allocatable, target :: OutFaceGidx(:)

    !% storage for all real data on all output faces for limited time levels
    !% (outelement, type, time-level) real data -- not coarray
    real(8), allocatable, target :: OutFaceDataR(:,:,:)
    !% storage for fixed integer data needed in output faces -- not coarray
    !% (outface, dataindex)
    integer, allocatable, target :: OutFaceFixedI(:,:)

    integer, allocatable, target :: thisElementOut(:), thisFaceOut(:)

    integer, allocatable :: SWMMlink_num_elements(:) !% number of elements in each output link
    integer, allocatable :: SWMMnode_num_elements(:) !% number of elements in each output link
    integer, allocatable :: SWMMnode_num_faces(:)    !% number of faces in each output node

    !% Temporary arrays that don't fit in any of the standard array structures
    !% brh 20220401
    integer, allocatable, target :: temp_BCupI(:,:)  !% number of BCup faces.
    real(8), allocatable, target :: temp_BCupR(:,:)
    integer :: N_tempBCupI = 1
    integer :: N_tempBCupR = 4
    
    !% Profiling Timer
    type(wall_clk) :: timer

    !% profiling storage
    real(8), allocatable :: profiler_data(:,:)
    character (len=64), allocatable :: profiler_procedure_name(:)
    integer, allocatable :: profiler_procedure_level(:)

    !% counter for transects
    integer :: lastTransectIdx = 0

!% ===================================================================================
!% CONSTANTS
!% ===================================================================================

    !% note that nullvalueI < 0 is required
    integer, parameter :: nullvalueI = 998877        !% note this places limit on largest network!
    real(8), parameter :: nullvalueR = 9.98877d16
    logical, parameter :: nullvalueL = .false.
    real(8), parameter :: negoneR = -1.0
    real(8), parameter :: zeroR = 0.d0
    real(8), parameter :: oneR = 1.d0
    real(8), parameter :: twoR = 2.d0
    real(8), parameter :: threeR = 3.d0
    real(8), parameter :: fourR = 4.d0
    real(8), parameter :: fiveR = 5.d0
    real(8), parameter :: sixR = 6.d0
    real(8), parameter :: eightR = 8.d0
    real(8), parameter :: tenR = 10.d0
    real(8), parameter :: twentyR = 20.d0
    real(8), parameter :: twentyfourR = 24.d0
    real(8), parameter :: sixtyR = 60.d0
    real(8), parameter :: onehundredR = 100.d0
    real(8), parameter :: onethousandR = 1000.d0
    real(8), parameter :: pi = 4.d0*datan(1.d0)

    real(8), parameter :: oneOneThounsandthR = oneR / onethousandR
    real(8), parameter :: oneOneHundredthR = oneR / onehundredR
    real(8), parameter :: onetenthR = oneR / tenR
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

    integer            :: dummyI
    integer, parameter :: zeroI = 0
    integer, parameter :: oneI = 1
    integer, parameter :: twoI = 2
    integer, parameter :: threeI = 3
    integer, parameter :: fourI = 4
    integer, parameter :: fiveI = 5
    integer, parameter :: sixI = 6

    !% Number of objects
    !integer :: SWMM_N_subcatch
    !integer :: SWMM_N_link
    !integer :: SWMM_N_node
    !integer :: SWMM_N_pollutant
    !integer :: SWMM_N_control
    !integer :: SWMM_N_divider
    !integer :: setting%SWMMinput%N_link_transect  ! # of irregular cross-section transects defined for links
    !integer :: SWMM_N_transect_depth_items
    integer :: N_transect            ! # of irregular cross-section transects for elements
    integer :: N_transect_depth_items
    integer :: N_transect_area_items
    integer :: N_link
    integer :: N_node
    integer :: N_headBC
    integer :: N_flowBC
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
    !integer :: SWMM_N_curve
    integer :: N_Curve
    integer, target :: N_OutTypeElem
    integer, target :: N_OutTypeFace

    !% Number of API parameters  ! REVISED APPROACH 20220422brh
    !% brh20211207s
    !rm integer, parameter :: N_api_node_attributes = api_node_overflow
    !integer, parameter :: N_api_nodef_attributes = api_nodef_totalEnd
    !rm integer, parameter :: N_api_link_attributes = api_linkf_conduit_length
    !integer, parameter :: N_api_linkf_attributes = api_linkf_commonBreak-1
    !% brh20211207e
    !integer, parameter :: N_api_linkf_type_attributes  = api_linkf_typeBreak  - api_nodef_totalEnd -1
    !integer, parameter :: N_api_linkf_xsect_attributes = api_linkf_totalEnd   - api_linkf_typeBreak -1
    !integer, parameter :: N_api_total_linkf_attributes = N_api_linkf_attributes + N_api_linkf_type_attributes &
    !                                                   + N_api_linkf_xsect_attributes
    !integer, parameter :: N_api_total_table_attributes = api_table_refers_to

    !% Coarray variables
    integer :: max_caf_elem_N    ! size of all elem array in coarray
    integer :: max_caf_face_N    ! size of all face array in coarray

    !% assign status parameters for nodes
    integer, parameter :: nUnassigned = 998877
    integer, parameter :: nAssigned   = 1
    integer, parameter :: nDeferred   = -1

    ! assign status parameters for links
    integer, parameter :: lUnassigned = 998877
    integer, parameter :: lAssigned   = 1
    integer, parameter :: lDeferred   = -1

    ! Constants for Junction
    integer, parameter :: J_elem_add = max_branch_per_node+1 ! Supplement elements for junction !%brh20211219
    integer, parameter :: J_face_add = max_branch_per_node   ! Supplement faces for junction !%brh20211219

    ! default number of elements for different node types
    integer, parameter :: N_elem_nJ2      = 0 ! 2-link nodes are assigned to a single face
    integer, parameter :: N_elem_nJm      = max_branch_per_node+1 ! M-link nodes are assigned a maximum of 7 elements !%brh20211219
    integer, parameter :: N_elem_nStorage = 1 ! Storage nodes are assigned to 1 element
    integer, parameter :: N_elem_nBCdn    = 0 ! Downstream BC nodes are assigned to 0 element (only a face)
    integer, parameter :: N_elem_nBCup    = 0 ! Upstream BC nodes are assigned to 0 element (only a face)
    integer, parameter :: N_elem_nJ1      = 0 ! Upstream non-BC nodes are assigned to 0 elements (only a face)

    ! useful shortcuts
    !% NOTE: don't use setting%... structure in define_globals to prevent linking problems
    !rm 20210629 griano integer, parameter :: debuglevelall = 0 ! set to 1 to get print of subroutine calls

    integer, parameter   :: NoAdjust       = 1   !% no link length adjustment has done
    integer, parameter   :: OneSideAdjust  = 2   !% one sided link length adjustment has done
    integer, parameter   :: BothSideAdjust = 3   !% both sided link length adjustment has done
    integer, parameter   :: DiagAdjust     = 4   !% if the connected link is an diagnostic element

    ! default for edge and non-edge node
    integer, parameter :: nonEdgeNode = 0 ! Upstream BC nodes are assigned to 1 element
    integer, parameter :: EdgeNode    = 1 ! Edge node of a partition

    ! ! defaults from initial depth type in links
    ! integer, parameter :: Uniform           = 1
    ! integer, parameter :: LinearlyVarying   = 2
    ! integer, parameter :: ExponentialDecay  = 3

    ! default for number of functional storage
    integer :: N_FunctionalStorage = 0

    !% maximum cfl across the network. Used for reporting 
    real(8) :: cfl_max
    real(8) :: cfl_max_CC
    real(8) :: cfl_max_JBJM

    !% datetime related variables
    integer, parameter :: datedelta = 693594
    integer, parameter :: secsperday = 86400
    integer :: dayspermonth(12,2) = &
        reshape((/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, & ! normal years
        31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/), (/12,2/)) ! leap years


end module define_globals
