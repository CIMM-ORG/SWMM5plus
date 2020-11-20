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
        real    ::  Coef  = 0.1
    end type adjustVshapedFlowrateType

    !% setting%Method%AdjustWidthDepth
    type adjustWidthDepthType
        logical :: Apply               = .true.
        logical :: AddDownstreamBuffer = .false.
        real    :: AdjustFractionMax   = 0.05
        real    :: AdjustFraction      = 0.01
        real    :: SmallWidth          = 1e-8
        real    :: DownstreamMinLength = 0.0
        real    :: depthMaxExpected    = 5.0
        real    :: angleMinimum        = 0.1
        real    :: areaMaximum         = 2.0
        real    :: cellSizeTarget      = 10.0
    endtype adjustWidthDepthType

    !% setting%Limiter%BC
    type BClimiterType
        logical :: UseInflowFroudeNumberLimiter = .true.
        real    :: FroudeInflowMaximum = 1.5 ! max value of Fr at inflow
    endtype BClimiterType

    !%  setting%Limiter%flowrate
    type flowrateType
        logical :: UseFaceVolumeTransport = .true.
        real    :: FaceVolumeTransport = 0.5 ! Fraction of usptream volume that can be transported in on time step
    endtype flowrateType

    !%  setting%Limiter%Timescale
    type timescaleType
        real    :: Maximum      = 1e6
        real    :: Minimum      = 1e-6
    endtype timescaleType

    !%  setting%Limiter%Velocity
    type velocityType
        logical :: UseLimitMax  = .true.
        real    :: Maximum      = 20.0 ! m/s
    endtype velocityType
        
    !%  setting%DefaultAC%Switch
    type switchType
        real    :: Depth        = 0.9   ! switch to AC solver if depth/depthMax >= 0.9
        real    :: Area         = 0.9   ! switch to AC solver if area/areaMax >= 0.9
        real    :: Buffer       = 0.05  ! 5% buffer for the switch
    end type switchType


    !%  setting%DefaultAC%Anomaly
    type anomalyType
        logical :: DensityCorrection = .false.  ! density anomaly correction to handle residual (or non-convergence) of the AC
        real    :: DnsityLowcutoff   = 1e-10    ! Note: setting the lowcutoff to zero can cause small truncation error velocities 
        real    :: DensityHighcutoff = 0.1      ! to cause waves that buildup over time into a sloshing flow. The lowcutoff avoid this
        real    :: OpenPipFactor     = 1.0      ! fraction of residual that is corrected - generally should be 1.0
        real    :: FullPipeFactor    = 1.0      ! fraction of residual that is corrected - generally should be 1.0
    end type anomalyType

    !%  setting%DefaultAC%Iter
    type iterType
        integer :: Max         = 100         ! cutoff for AC iterations without convergence
        integer :: Min         = 3           ! reset by code so that itermin * dtau >= dt
        integer :: Firststep   = 100         ! allows more iterations in first time step
    end type itertype

    !%  setting%DefaultAC%CFL
    type cflType
        real    :: CFLmax      = 2.0        ! maximum cfl for the AC dtau -- may be higher than 1.0)
        real    :: CFLsmall    = 0.05       ! small maximum cfl when reset to larger dtau
    end type cflType

    !%  setting%DefaultAC%Celerity
    type celerityType
        real    :: RC          = 1.0        ! celerity ratio of AC wave speed to gravity wave speed (1.0 works)
    end type celerityType

    !%  setting%DefaultAC%Convergence
    type convergenceType
        real    :: Hrelative   = 1e-2       ! AC convergence relative change in L2 norm when cutting off AC
        real    :: Qrelative   = 1e-2       ! AC convergence relative change in L2 norm when cutting off AC
        real    :: Habsolute   = 1e-5       ! AC convergence change in absolute L2 norm when cutting off AC
        real    :: Qabsolute   = 1e-5       ! AC convergence change in absolute L2 norm when cutting off AC
    end type convergenceType

    !%  setting%DefaultAC%dtauFactor
    type dtaufactorType
        real :: dtdtau      = 1.0 / 3.0     ! ratio of dt/dtau (changes with flow)
        real :: dtdtauMax   = 1.0 / 3.0     ! ensure the minimum number of iterations.
    end type dtaufactorType

    !% SECOND LEVEL TYPES ----------------------------------------------

    !%  setting%Constant
    type constantType
        real  :: gravity = 9.81  ! m^2/s
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
        real    :: FroudeJump                     = 0.1 ! +- small non-dimensional range for hyd jump discrimination
        real    :: InflowDepthIncreaseFroudeLimit = 0.1 ! Fractional increase in depth under froude limitation
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
        logical ::  UseSmallVolumes = .true. ! YN to determine if smallvolume adjustments used
        real    ::  DepthCutoff      = 0.01  ! m Determines where small volumes begin
        real    ::  ManningsN        = 0.01
        real    ::  MinimumTopwidth  = 0.5   ! m   Minimum value used for smallvolume reset
        real    ::  MinimumArea      = 0.005  ! m^2
        real    ::  MinimumPerimeter = 0.52  ! m
        real    ::  MinimumHydRadius = 0.009  ! m
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
        real :: Dt  ! s
        real :: StartTime = 0.0
        real :: EndTime   = 0.0
        real :: ThisTime
        real :: NextTime
    end type timeType

    !%  setting%ZeroValue
    type zerovalueType
        logical :: UseZeroValues = .true.
        real    :: Area         = 1.0e-6  ! m^2
        real    :: Depth        = 1.0e-4  ! m
        real    :: Flowrate     = 0.0     ! m^3/s
        real    :: Topwidth     = 1.0e-4  ! m
        real    :: Velocity     = 0.0     ! m/s
        real    :: Volume       = 1.0e-6  ! m^3 !%%%%%%%%%%%%%%%%%%I changed it from 1.0e-6
    end type zerovalueType

    !%  setting%DefaultAC
    type defaultACType
        character(len=64)       :: Tsource     = 'T00' ! 'T10', 'T20'
        character(len=64)       :: TimeStencil = 'backwards3'   ! 'CN' 
        logical                 :: PrintConvergence = .false.
        real                    :: dtau  
        type(dtaufactorType)    :: dtauFactor
        type(convergenceType)   :: Convergence
        type(celerityType)      :: Celerity
        type(cflType)           :: CLF 
        type(iterType)          :: Iter
        type(anomalyType)       :: Anomaly
        type(switchType)        :: Switch
    end type defaultACType

    !%  setting%BCondition
    type bconditionType
        logical                 :: InflowRampup     = .true.
        real                    :: InflowRampupTime = 100.0
        real                    :: flowrateIC       = 0.01
    end type bconditionType 

    !%  settin%HydJump
    type hydjumpType
        character(len=64)       :: FaceHydJumpMethod = 'momentum_match' !'extrapolate_surface' ! 'energy_limit', 'momentum_match', 'rectangular', 'simple'
    end type hydjumpType

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
        type(stepType)              :: Step          ! controls over simulation time stepping
        type(testcaseType)          :: TestCase      ! custom setup for test cases
        type(timeType)              :: Time          ! controls of time step
        type(zerovalueType)         :: ZeroValue     ! finite values to represent small or negative values
        type(defaultACType)         :: DefaultAC     ! This contains the default settings for atrificial compressibility
        type(bconditionType)        :: BCondition    ! This contains rampup boundary condition
        type(hydjumpType)           :: HydJump       ! Hydraulic jump method type
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

        setting%Debugout%DisplayInterval = 1

        setting%Time%dt = 0.005
        setting%Step%First = 1
        setting%Step%Final = 8000
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
