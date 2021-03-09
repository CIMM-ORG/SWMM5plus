! module type defintions
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

    public

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

    !% boundary condition data
    type bcType
        integer :: Idx
        integer :: NodeID
        integer :: FaceID
        integer :: ElemGhostID
        integer :: ElemInsideID
        integer :: Updn      ! bc_updn_...  (0 = upstream,  1 = downstream)
        integer :: Category  ! bc_category_... (0 = elevation, 1 = inflowrate)
        real(8), dimension(:), allocatable :: TimeArray
        real(8), dimension(:), allocatable :: ValueArray
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


    !==========================================================================
    ! END OF MODULE type_definitions
    !==========================================================================
end module type_definitions
