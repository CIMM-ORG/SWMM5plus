!==========================================================================
 module setting_definition
! 
! defines the setting%.... structure that contains the controls for the 
! simulation
!
    implicit none
    
    public

!% THIRD LEVEL TYPES ----------------------------------------------

    !% setting%Method%AdjustVshapedFlowrate
    type adjustVshapedFlowrateType
        logical ::  Apply = .true.
        real    ::  Coef  = 0.1
    end type adjustVshapedFlowrateType

    !%  setting%Limiter%flowrate   
    type flowrateType
        real    :: FaceVolumeTransport = 0.5 ! Fraction of usptream volume that can be transported in on time step
    endtype flowrateType

    !%  setting%Limiter%Timescale
    type timescaleType
        real    :: Maximum      = 1e6
        real    :: Minimum      = 1e-6
    endtype timescaleType
    
    !%  setting%Limiter%Velocity    
    type velocityType
        real    :: Maximum      = 10 ! m/s
    endtype velocityType

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
        character(len=64) :: FolderName  = 'debugoutputA'
        character(len=64) :: FileName    = 'debug'
        character(len=256) :: FolderPath    = './'
    end type

    !%  setting%Eps
    type epsilonType
        real    :: FroudeJump                     = 0.1 ! +- small non-dimensional range for hyd jump discrimination  
        real    :: InflowDepthIncreaseFroudeLimit = 0.1 ! Fractional increase in depth under froude limitation
    end type epsilonType  

    !%  setting%Limiter    
    type limiterType
        type(flowrateType   :: Flowrate
        type(velocityType)  :: Velocity   
        type(timescaleType) :: Timescale
    endtype limiterType
    
    !%  setting%Method
    type methodType
        type(adjustVshapedFlowrateType) :: AdjustVshapedFlowrate
    end type methodType
  
    !%  setting%SmallVolume
    type smallvolumeType
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
        real    :: Area         = 1.0e-7  ! m^2
        real    :: Depth        = 1.0e-4  ! m
        real    :: Flowrate     = 0.0     ! m^3/s
        real    :: Topwidth     = 1.0e-4  ! m
        real    :: Velocity     = 0.0     ! m/s
        real    :: Volume       = 1.0e-6  ! m^3
    end type zerovalueType
    
    

!% FIRST LEVEL TYPE  ----------------------------------------------      
    type settingType
        integer :: dummy
        type(constantType)    :: Constant     ! constants
        type(debugoutType)    :: DebugOut     ! control of debougout files
        type(epsilonType)     :: Eps          ! epsilons used to provide bandwidth for comparisons
        type(limiterType)     :: Limiter      ! maximum and minimum limiters
        type(methodType)      :: Method       ! controls over simulation methods
        type(smallvolumeType) :: SmallVolume   ! controls for small volumes
        type(stepType)        :: Step         ! controls over simulation time stepping
        type(testcaseType)    :: TestCase     ! custom setup for test cases
        type(timeType)        :: Time         ! controls of time step
        type(zerovalueType)   :: ZeroValue    ! finite values to represent small or negative values
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

 setting%Debugout%DisplayInterval = 1000

 setting%Time%dt = 50.0
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