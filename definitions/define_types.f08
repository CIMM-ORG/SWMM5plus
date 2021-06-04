! module define_types
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
module define_types

    implicit none

    type string
        character(len=:), allocatable :: str
    end type string

    !%  control data
    type controlType
        integer :: Idx
        integer :: LinkId
        integer :: ElemId
        real(8), dimension(:), allocatable :: TimeArray
        real(8), dimension(:), allocatable :: HeightArray
        real(8), dimension(:), allocatable :: AreaArray
        real(8)    :: HeightStart
        real(8)    :: HeightNow
        real(8)    :: AreaNow
        real(8)    :: AreaPrior
        real(8)    :: FullDepth
        real(8)    :: GateTimeChange1
        real(8)    :: GateTimeChange2
        real(8)    :: GateHeightChange1
        real(8)    :: GateHeightChange2
        real(8)    :: HeightMinimum
        real(8)    :: GateSpeed
        logical :: CanMove
        logical :: MovedThisStep
    end type controlType

    !%  diagnostic%Volume
    type diagnosticVolumeType
        integer  :: Step
        real(8)     :: Time
        real(8)     :: Volume
        real(8)     :: VolumeChange
        real(8)     :: NetInflowVolume
        real(8)     :: InflowRate
        real(8)     :: OutflowRate
        real(8)     :: ConservationThisStep ! + is artificial source, - is sink
        real(8)     :: ConservationTotal
    end type diagnosticVolumeType

    type diagnosticType
        type(diagnosticVolumeType)  :: Volume
    end type diagnosticType

    type real_array
        integer :: max_size = 0
        integer :: len = 0
        real(8), allocatable :: array(:)
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
        real(8), allocatable :: TimeArray(:)
        real(8), allocatable :: ValueArray(:)
        real(8)    :: ThisValue
        real(8)    :: ThisTime
        real(8)    :: ThisFlowrate
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

    ! PATTERN OBJECT
    type pattern
        integer :: ptype
        integer :: count
        real(8), dimension(24) :: factor
    end type pattern

    ! EXTERNAL INFLOW OBJECT
    type totalInflow
        integer :: node_id ! index to element thar receives inflow
        type(real_table) :: xy
        ! t_series*sfactor + base_pat*baseline
        integer :: ext_t_series = -1 ! time_series
        integer :: ext_base_pat = -1 ! pattern
        real(8) :: ext_baseline = 0! constant baseline value
        real(8) :: ext_sfactor = 0! time series scaling factor
        ! ---------------------------------------------------------
        real(8) :: dwf_avgValue = 0 ! average inflow value
        integer :: dwf_monthly_pattern = -1
        integer :: dwf_daily_pattern = -1
        integer :: dwf_hourly_pattern = -1
        integer :: dwf_weekend_pattern = -1
    end type totalInflow

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
end module define_types