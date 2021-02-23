!
! 2019-10 @EhsanMadadi added keys for widt-depth cases.
!==========================================================================
module setting_definition
    !
    ! defines the setting%.... structure that contains the controls for the
    ! simulation
    !

    use type_definitions

    implicit none

    public

    !% THIRD LEVEL TYPES ----------------------------------------------

    !% setting%Method%AdjustVshapedFlowrate
    type adjustVshapedFlowrateType
        logical ::  Apply = .true.
        real(8)    ::  Coef  = 1.0
    end type adjustVshapedFlowrateType

    !% setting%Method%AdjustWidthDepth
    type adjustWidthDepthType
        logical :: Apply               = .true.
        logical :: AddDownstreamBuffer = .false.
        real(8)   :: AdjustFractionMax   = 0.05
        real(8)   :: AdjustFraction      = 0.01
        real(8)   :: SmallWidth          = 1e-8
        real(8)   :: DownstreamMinLength = 0.0
        real(8)   :: depthMaxExpected    = 5.0
        real(8)   :: angleMinimum        = 0.1
        real(8)   :: areaMaximum         = 2.0
        real(8)   :: cellSizeTarget      = 10.0
    endtype adjustWidthDepthType

    !% setting%Limiter%BC
    type BClimiterType
        logical :: UseInflowFroudeNumberLimiter = .true.
        real(8)   :: FroudeInflowMaximum = 1.5 ! max value of Fr at inflow
    endtype BClimiterType

    !%  setting%Limiter%flowrate
    type flowrateType
        logical :: UseFaceVolumeTransport = .true.
        real(8)   :: FaceVolumeTransport = 0.5 ! Fraction of usptream volume that can be transported in on time step
    endtype flowrateType

    !%  setting%Limiter%Timescale
    type timescaleType
        real(8)   :: Maximum      = 1e6
        real(8)   :: Minimum      = 1e-6
    endtype timescaleType

    !%  setting%Limiter%Velocity
    type velocityType
        logical :: UseLimitMax  = .true.
        real(8)   :: Maximum      = 20 ! m/s
    endtype velocityType

    !% SECOND LEVEL TYPES ----------------------------------------------

    !%  setting%Constant
    type constantType
        real(8) :: gravity = 9.81  ! m^2/s
    end type constantType

    !%  setting%DebugOut
    type debugoutType
        logical :: elem2R  = .false.
        logical :: elem2I  = .false.
        logical :: elem2YN = .false.
        logical :: elemMR  = .false.
        logical :: elemMI  = .false.
        logical :: elemMYN = .false.
        logical :: faceR   = .false.
        logical :: faceI   = .false.
        logical :: faceYN  = .false.
        logical :: nodeR   = .false.
        logical :: nodeI   = .false.
        logical :: nodeYN  = .false.
        logical :: linkR   = .false.
        logical :: linkI   = .false.
        logical :: linkYN  = .false.
        integer :: DisplayInterval = 1
        logical :: SuppressAllFiles  = .false.
        logical :: SuppressTimeStep  = .false.
        logical :: SuppressTimeValue = .false.
        logical :: SuppressNdat      = .false.
        character(len=64)  :: FolderName  = 'debugoutputA'
        character(len=64)  :: FileName    = 'debug'
        character(len=256) :: FolderPath    = './'
    end type

    !%  setting%Eps
    type epsilonType
        real(8)   :: FroudeJump                     = 0.1 ! +- small non-dimensional range for hyd jump discrimination
        real(8)   :: InflowDepthIncreaseFroudeLimit = 0.1 ! Fractional increase in depth under froude limitation
    end type epsilonType

    !%  setting%Limiter
    type limiterType
        type(BClimiterType) :: BC
        type(flowrateType)  :: Flowrate
        type(velocityType)  :: Velocity
        type(timescaleType) :: Timescale
    endtype limiterType

    !%  setting%Method
    type methodType
        type(adjustVshapedFlowrateType) :: AdjustVshapedFlowrate
        type(adjustWidthDepthType)      :: AdjustWidthDepth
    end type methodType

    !%  setting%OutputThreadedLink
    type outputThreadedLinkType
        integer               :: DisplayInterval = 1
        logical               :: SuppressAllFiles  = .false.
        logical               :: UseThisOutput  = .false.
        logical               :: area           = .false.
        logical               :: flowrate       = .false.
        logical               :: velocity       = .false.
        logical               :: eta            = .false.
        logical               :: depth          = .false.
        character(len=64)     :: FolderName  = 'OutputThreaded'
        character(len=64)     :: FileName    = 'out'
        character(len=256)    :: FolderPath    = './'
    end type outputThreadedLinkType

    !%  setting%SmallVolume
    type smallvolumeType
        ! Dont using small volumes for weir case. Needed to be changed later
        logical ::  UseSmallVolumes = .false. ! YN to determine if smallvolume adjustments used
        real(8)   ::  DepthCutoff      = 0.01  ! m Determines where small volumes begin
        real(8)   ::  ManningsN        = 0.01
        real(8)   ::  MinimumTopwidth  = 0.5   ! m   Minimum value used for smallvolume reset
        real(8)   ::  MinimumArea      = 0.005  ! m^2
        real(8)   ::  MinimumPerimeter = 0.52  ! m
        real(8)   ::  MinimumHydRadius = 0.009  ! m
    end type smallvolumeType

    !%  setting%Step
    type stepType
        integer :: First
        integer :: Final
        integer :: Current
        integer :: FromRestart
    end type stepType

    !%  setting%TestCase
    type testcaseType
        logical       :: UseTestCase
        character(64) :: TestName
    end type testcaseType

    !%  setting%Time
    type timeType
        character(14)   ::  DateTimeStamp
        real(8) :: Dt  ! s
        real(8) :: StartTime = 0.0
        real(8) :: EndTime   = 0.0
        real(8) :: ThisTime
        real(8) :: NextTime
    end type timeType

    !%  setting%ZeroValue
    type zerovalueType
        logical :: UseZeroValues = .true.
        real(8)   :: Area         = 1.0e-7  ! m^2
        real(8)   :: Depth        = 1.0e-4  ! m
        real(8)   :: Flowrate     = 0.0     ! m^3/s
        real(8)   :: Topwidth     = 1.0e-4  ! m
        real(8)   :: Velocity     = 0.0     ! m/s
        real(8)   :: Volume       = 1.0e-7  ! m^3 !%%%%%%%%%%%%%%%%%%I changed it from 1.0e-6
    end type zerovalueType

    !%  setting%Weir
    type WeirType
        real(8)   :: WeirDischargeCoeff   = 1.4 ! m^3/s
        real(8)   :: EndContraction       = 0.0    ! Number of End Contraction(0, 1, 2)
        real(8)   :: WeirHeight           = 1.5  ! Vertical Height of Weir Opening m
        real(8)   :: WeirWidth            = 3.0
        real(8)   :: WeirSideSlope        = 1.0
        real(8)   :: WeirInletOffset      = 1.0
    end type WeirType

    !% FIRST LEVEL TYPE  ----------------------------------------------
    type settingType
        integer :: dummy
        type(constantType)          :: Constant     ! constants
        type(debugoutType)          :: DebugOut     ! control of debougout files
        type(epsilonType)           :: Eps          ! epsilons used to provide bandwidth for comparisons
        type(limiterType)           :: Limiter      ! maximum and minimum limiters
        type(methodType)            :: Method       ! controls over simulation methods
        type(outputThreadedLinkType):: OutputThreadedLink ! controls output for threaded link
        type(smallvolumeType)       :: SmallVolume   ! controls for small volumes
        type(stepType)              :: Step         ! controls over simulation time stepping
        type(testcaseType)          :: TestCase     ! custom setup for test cases
        type(timeType)              :: Time         ! controls of time step
        type(zerovalueType)         :: ZeroValue    ! finite values to represent small or negative values
        type(WeirType)              :: Weir         ! This contains the weir settings (only for initial implimentation purpose)
    end type settingType



    type(settingType), target :: setting


    integer, private :: debuglevel = 0

contains
    !==========================================================================
    !==========================================================================
    subroutine setting_default
        !
        !   This subroutine can be used in conjunction with the test_case to setup
        !   a simulation or set of simulations
        !
        !--------------------------------------------------------------------------

        setting%Debugout%DisplayInterval = 10

        setting%Time%dt = 5
        setting%Step%First = 1
        setting%Step%Final = 40000
        setting%Step%Current = 1
        setting%Step%FromRestart = 1

        setting%Time%StartTime = 0.0
        setting%Time%EndTime = setting%Time%StartTime  &
            + setting%Time%dt * (setting%Step%Final - setting%Step%First + 1)

    end subroutine setting_default
    !
    !==========================================================================
    ! END OF MODULE
    !==========================================================================
end module setting_definition
