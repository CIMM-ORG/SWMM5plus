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
    end type NodePack

    type NodeArray
        integer,      allocatable :: I(:,:)   !% integer data for nodes
        real(8),      allocatable :: R(:,:)   !% real data for nodes
        logical,      allocatable :: YN(:,:)  !% logical data for nodes
        type(string), allocatable :: Names(:) !% names for nodes retrieved from EPA-SWMM
        type(NodePack)            :: P        !% packs for nodes
    end type NodeArray

    type LinkArray
        integer,      allocatable :: I(:,:)   !% integer data for links
        real(8),      allocatable :: R(:,:)   !% real data for links
        logical,      allocatable :: YN(:,:)  !% logical data for links
        type(string), allocatable :: Names(:) !% names for links retrieved from EPA-SWMM
    end type LinkArray

    type BCPack
        !% We initialze the pack array for BC interpolation use later
        integer, allocatable :: BCup(:)
        integer, allocatable :: BClat(:)
        integer, allocatable :: BCdn(:)
    end type BCPack

    type BCArray
        integer,     allocatable :: flowI(:,:)              !% integer data for inflow BCs
        real(8),     allocatable :: flowR_timeseries(:,:,:) !% time series data for inflow BC
        integer,     allocatable :: flowIdx(:)              !% indexes of current entry in flowR_timeseries
        real(8),     allocatable :: flowRI(:)               !% values of interpolated inflows at current time
        integer,     allocatable :: headI(:,:)              !% integer data for elevation BCs
        real(8),     allocatable :: headR_timeseries(:,:,:) !% time series data for elevation BC
        integer,     allocatable :: headIdx(:)              !% indexes of current entry in headR_timeseries
        real(8),     allocatable :: headRI(:)               !% values of interpolated heads at current time
        type(BCPack)             :: P                       !% packs of boundary conditions
    end type BCArray
end module define_types