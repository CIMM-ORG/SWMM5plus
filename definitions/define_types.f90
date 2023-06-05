module define_types
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Defines types used globally in SWMM5+
    !%
    !% Methods
    !% These are derived type definitions that are used in globals, setting, and
    !% elsewhere.
    !
    !% Note that each of these is intended to be dependent only on derived
    !% types within this module.
    !
    !% Types that are defined outside of here or setting_definition should
    !% be confined to the module in which they are defined.
    !%==========================================================================
    implicit none

    type string
        character(len=:), allocatable :: str
    end type string

    !% curve data
    type curveType
        integer :: ID
        integer :: Type
        integer :: RefersTo
        integer :: NumRows
        integer :: ElemIdx
        integer :: FaceIdx
        real(8), dimension(:,:), allocatable :: ValueArray
    end type

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

    type steady_state_record
        character(len=52) :: id_time
        real(8) :: flowrate
        real(8) :: wet_area
        real(8) :: depth
        real(8) :: froude
    end type steady_state_record

    !% ==============================================================
    !% Arrays
    !% ==============================================================

    type NodePack
        integer,      allocatable :: have_flowBC(:)
        integer,      allocatable :: have_headBC(:)
        integer,      allocatable :: have_output(:)
    end type NodePack

    type LinkPack
        integer, allocatable :: have_output(:)
    end type LinkPack

    type NodeArray
        integer,      allocatable :: I(:,:)[:]   !% integer data for nodes
        real(8),      allocatable :: R(:,:)[:]   !% real data for nodes
        logical,      allocatable :: YN(:,:)[:]  !% logical data for nodes
        type(string), allocatable :: Names(:) !% names for nodes retrieved from EPA-SWMM
        type(NodePack)            :: P        !% packs for nodes
    end type NodeArray

    type LinkArray
        integer,      allocatable :: I(:,:)[:]   !% integer data for links
        real(8),      allocatable :: R(:,:)[:]   !% real data for links
        logical,      allocatable :: YN(:,:)[:]  !% logical data for links
        integer,      allocatable :: transectI(:,:) !% integer data for irregular transects
        real(8),      allocatable :: transectR(:,:) !% real data for irregular transects
        real(8),      allocatable :: transectTableDepthR(:,:,:) !% table of real values by uniform depth distribution
        real(8),      allocatable :: transectTableAreaR(:,:,:)  !% table of real values by uniform area distribution 
        type(string), allocatable :: transectID(:) !% string transect ID from EPA-SWMM input
        type(string), allocatable :: Names(:) !% names for links retrieved from EPA-SWMM
        type(LinkPack)            :: P        !% pack for links
    end type LinkArray

    type :: BoundaryElemArray
        real(8),      allocatable :: R(:,:)   !% real data for elemB
        integer,      allocatable :: I(:,:)   !% integer data for elemB
    end type BoundaryElemArray

    type BCPack
        !% --- initialization of the the pack array for BC interpolation is later
        integer, allocatable :: BCup(:)
        integer, allocatable :: BClat(:)
        integer, allocatable :: BCdn(:)
    end type BCPack

    type BCArray
        integer,     allocatable :: flowI(:,:)              !% integer data for inflow BCs
        real(8),     allocatable :: flowR(:,:)
        logical,     allocatable :: flowYN(:,:)             !% logical data for inflow BCs
        real(8),     allocatable :: flowTimeseries(:,:,:) !% time series data for inflow BC
        integer,     allocatable :: headI(:,:)              !% integer data for elevation BCs
        real(8),     allocatable :: headR(:,:)
        logical,     allocatable :: headYN(:,:)             !% logical data for head BCs
        real(8),     allocatable :: headTimeseries(:,:,:) !% time series data for elevation BC
        type(BCPack)             :: P                       !% packs of boundary conditions
    end type BCArray

    !% Dynamic Array
    type f_array
        integer :: max_size = 0
        integer :: len = 0
        real, allocatable :: arr(:)
    end type f_array

    type job
        type(f_array) :: time_stamps
        real :: end_time
        integer :: id = -1
    end type job
    
    !% Profiling types
    type wall_clk
       type(job), allocatable :: jobs(:)
       integer :: max_num_jobs = 0
       integer :: num_jobs = 0
    end type wall_clk
 
end module define_types
