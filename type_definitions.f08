! module type_definitions
!
! These are derived type definitions that are used in globals, setting, and
! elsewhere.
!
! Note that each of these is intended to be dependent only on derived
! types within this module.
!
! Types that are defined outside of here or setting_definition should
! be confined to the module in which they are defined.
!
!==========================================================================
module type_definitions

    implicit none

    type string
        character(len=:), allocatable :: str
    end type string

    !%  diagnostic%Volume
    type diagnosticVolumeType
        integer  :: Step
        real     :: Time
        real     :: Volume
        real     :: VolumeChange
        real     :: NetInflowVolume
        real     :: InflowRate
        real     :: OutflowRate
        real     :: ConservationThisStep ! + is artificial source, - is sink
        real     :: ConservationTotal
    end type diagnosticVolumeType

    type diagnosticType
        type(diagnosticVolumeType)  :: Volume
    end type diagnosticType

    type real_array
        integer :: max_size = 0
        integer :: len = 0
        real, allocatable :: array(:)
    end type real_array

    !% boundary condition data
    type bcType
        integer :: Idx
        integer :: NodeID
        integer :: FaceID
        integer :: ElemGhostID
        integer :: ElemInsideID
        integer :: Updn      ! bc_updn_...  (0 = upstream,  1 = downstream)
        integer :: Category  ! bc_category_... (0 = elevation, 1 = inflowrate)
        real, allocatable :: TimeArray
        real, allocatable :: ValueArray
        real    :: ThisValue
        real    :: ThisTime
        real    :: ThisFlowrate
    end type bcType

    !% output file location
    type outputfileType
        integer        :: UnitNumber = 0
        character(256) :: FileName   = 'dummy.txt'
        character(256) :: FolderName = 'dummyFolder'
        character(256) :: FolderPath = './'
        character(32)  :: FileStatus = 'new'
        character(512) :: WriteName  = ''
        logical        :: IsOpen     = .false.
    end type outputfileType

    !% specific information for a debug file
    type debugfileType
        type (outputfileType) :: FileInfo
        character(32)         :: ArrayName
        integer               :: ColumnIndex
    end type debugfileType

    !% specific information for a threaded output file
    type threadedfileType
        type (outputfileType) :: FileInfo
        character(32)         :: DataName
    end type threadedfileType

    type integer_array
        integer :: max_size = 0
        integer :: len = 0
        integer, allocatable :: array(:)
    end type integer_array

    ! TABLE OBJECT
    type real_table
        integer :: table_type
        integer :: dim
        integer, allocatable :: tsize(:)
        type(real_array), allocatable :: data(:)
    end type real_table

    ! TIME SERIES OBJECT
    type tseries
        integer :: current = 1
        type(real_table) :: table
    end type tseries

    ! PATTERN OBJECT
    type pattern
        integer :: ptype
        integer :: count
        real, dimension(24) :: factor
    end type pattern

    ! EXTERNAL INFLOW OBJECT
    type totalInflow
        integer :: node_id ! index to element thar receives inflow
        type(real_array) :: xy(2)
        ! t_series*sfactor + base_pat*baseline
        integer :: ext_t_series ! time_series
        integer :: ext_base_pat ! pattern
        real :: ext_baseline ! constant baseline value
        real :: ext_sfactor ! time series scaling factor
        ! ---------------------------------------------------------
        real :: dwf_avgValue ! average inflow value
        integer :: dwf_monthly_pattern
        integer :: dwf_daily_pattern
        integer :: dwf_hourly_pattern
        integer :: dwf_weekly_pattern
    end type

    type graph_node
        integer :: node_id
        type(integer_array) :: neighbors
        type(integer_array) :: link_id
        type(real_array) :: neighbor_flows
    end type graph_node

    type graph
        integer :: num_vertices
        type(graph_node), allocatable, dimension(:) :: g ! graph linked lists
        integer, allocatable, dimension(:) :: in_degree ! list with in-degrees of node
    end type graph

    ! --- File Handling
    type steady_state_record
        character(len=52) :: id_time
        real(8) :: flowrate
        real(8) :: wet_area
        real(8) :: depth
        real(8) :: froude
    end type steady_state_record
    !==========================================================================
    ! END OF MODULE type_definitions
    !==========================================================================
end module type_definitions
