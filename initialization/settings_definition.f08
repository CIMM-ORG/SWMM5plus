module settings_definition

    use data_keys
    use json_module

    implicit none
    public

    ! Third Level Types

    ! -
    ! --
    ! ---

    ! setting%ACmethod%Anomaly
    type ACmethodAnomalyType
        ! density anomaly correction to handle residual (or non-convergence) of the AC
        logical :: UseDensityCorrection = .false.
        ! Note: setting the lowcutoff to zero can cause small truncation error velocities
        ! to cause waves that buildup over time into a sloshing flow. The lowcutoff avoid this
        real(8) :: DensityLowCutoff = 1e-10
        real(8) :: DensityHighCutoff = 0.1
        ! fraction of residual that is corrected - generally should be 1.0
        real(8) :: OpenPipeFactor = 1.0
        ! fraction of residual that is corrected - generally should be 1.0
        real(8) :: FullPipeFactor = 1.0
    end type ACmethodAnomalyType

    ! setting%ACmethod%CFL
    type ACmethodCFLType
        ! maximum cfl for the AC dtau -- may be higher than 1.0)
        real(8) :: CFLmax = 2.0
        ! small maximum cfl when reset to larger dtau
        real(8) :: CFLsmall = 0.05
    end type ACmethodCFLType

    ! setting%ACmethod%Celerity
    type ACmethodCelerityType
        ! celerity ratio of AC wave speed to gravity wave speed (1.0 works)
        real(8) :: RC = 1.0
    end type ACmethodCelerityType

    ! setting%ACmethod%Convergence
    type ACmethodConvergenceType
        ! AC convergence relative change in L2 norm when cutting off AC
        real(8) :: Hrelative = 1e-2
        ! AC convergence relative change in L2 norm when cutting off AC
        real(8) :: Qrelative = 1e-2
        ! AC convergence change in absolute L2 norm when cutting off AC
        real(8) :: Habsolute = 1e-5
        ! AC convergence change in absolute L2 norm when cutting off AC
        real(8) :: Qabsolute = 1e-5
    end type ACmethodConvergenceType

    ! setting%ACmethod%Iter
    type ACmethodIterType
        ! cutoff for AC iterations without convergence
        integer :: Max = 100
        ! reset by code so that itermin * dtau >  dt
        integer :: Min = 3
        ! allows more iterations in first time step
        integer :: Firststep = 100
    end type ACmethodIterType

    ! setting%ACmethod%Switch
    type ACmethodSwitchType
        ! switch to AC solver if depth/depthMax >  0.9
        real(8) :: Depth = 0.9
        ! switch to AC solver if area/areaMax >  0.9
        real(8) :: Area = 0.9
        ! 5% buffer for the switch
        real(8) :: Buffer = 0.05
    end type ACmethodSwitchType

    ! setting%Adjust%Flowrate
    type AdjustFlowrateType
        logical :: Apply = .true.
        real(8) :: Coef = 1.0
        integer :: Approach = vshape
    end type AdjustFlowrateType

    ! setting%Adjust%Head
    type AdjustHeadType
        logical :: Apply = .true.
        real(8) :: Coef = 1.0
        integer :: Approach = vshape
    end type AdjustHeadType

    ! setting%Adjust%WidthDepth
    type AdjustWidthDepthType
        logical :: Apply = .true.
    endtype AdjustWidthDepthType

    ! setting%Limiter%BC
    type LimiterBCType
        logical :: UseInflowLimiter = .true.
        integer :: Approach = FroudeNumber
        ! max value of Fr at inflow
        real(8) :: FroudeInflowMaximum = 1.5
    endtype LimiterBCType

    ! setting%Limiter%Flowrate
    type LimiterFlowrateType
        logical :: UseFaceVolumeTransport = .true.
        ! Fraction of usptream volume that can be
        ! transported in on time step
        real(8) :: FaceVolumeTransport = 0.5
    endtype LimiterFlowrateType

    ! setting%Limiter%Timescale
    type LimiterTimescaleType
        real(8) :: Maximum = 1e6
        real(8) :: Minimum = 1e-6
    endtype LimiterTimescaleType

    ! setting%Limiter%Velocity
    type LimiterVelocityType
        logical :: UseLimitMax = .true.
        real(8) :: Maximum = 10.0 ! m/s
    endtype LimiterVelocityType

    ! -
    ! --
    ! ---

    ! Second Level Types

    ! -
    ! --

    ! setting%ACmethodType
    type ACmethodType
        type(ACmethodAnomalyType) :: Anomaly
        type(ACmethodCFLType) :: CFL
        type(ACmethodCelerityType) :: Celerity
        type(ACmethodConvergenceType) :: Convergence
        type(ACmethodIterType) :: Iter
        type(ACmethodSwitchType) :: Switch
    end type ACmethodType

    ! setting%Adjust
    type AdjustType
        type(AdjustFlowrateType) :: Flowrate
        type(AdjustHeadType) :: Head
        type(AdjustWidthDepthType) :: WidthDepth
    end type AdjustType

    ! setting%Constant
    type ConstantType
        real(8) :: gravity = 9.81 ! m^2/s
    end type ConstantType

    ! setting%DebugOut
    type DebugOutType
        logical :: elem2R = .false.
        logical :: elem2I = .false.
        logical :: elem2YN = .false.
        logical :: elemMR = .false.
        logical :: elemMI = .false.
        logical :: elemMYN = .false.
        logical :: faceR = .false.
        logical :: faceI = .false.
        logical :: faceYN = .false.
        logical :: nodeR = .false.
        logical :: nodeI = .false.
        logical :: nodeYN = .false.
        logical :: linkR = .false.
        logical :: linkI = .false.
        logical :: linkYN = .false.
        integer :: DisplayInterval = 1
        logical :: SuppressAllFiles = .false.
        logical :: SuppressTimeStep = .false.
        logical :: SuppressTimeValue = .false.
        logical :: SuppressNdat = .false.
        character(len=12) :: FolderName = 'debugoutputA'
        character(len=5) :: FileName = 'debug'
        character(len=256) :: FolderPath = './'
    end type DebugOutType

    ! setting%Eps
    type EpsilonType
        ! +- small non-dimensional range for hyd jump discrimination
        real(8) :: FroudeJump = 0.1
        ! Fractional increase in depth under froude limitation
        real(8) :: InflowDepthIncreaseFroudeLimit = 0.1
    end type EpsilonType

    ! setting%Limiter
    type LimiterType
        type(LimiterBCType) :: BC
        type(LimiterFlowrateType) :: Flowrate
        type(LimiterTimescaleType) :: Timescale
        type(LimiterVelocityType) :: Velocity
    endtype LimiterType

    ! setting%SmallVolume
    type SmallVolumeType
        ! Dont using small volumes for weir case.
        ! Needed to be changed later SmallVolumeType
        logical :: UseSmallVolumes = .true.
        real(8) :: DepthCutoff = 0.01 ! m
        real(8) :: ManningsN = 0.01
        real(8) :: MinimumArea = 0.005 ! m^2
        real(8) :: MinimumHydRadius = 0.009 ! m
        real(8) :: MinimumPerimeter = 0.52 ! m
        real(8) :: MinimumTopwidth = 0.5 ! m
    end type SmallVolumeType

    ! setting%Solver
    type SolverType
        real(8), dimension(2) :: crk2 = \(0.5, 1.0\)
        integer :: MomentumSourceM = T00
        logical :: PreissmanSlot = .true.
        integer :: SolverSelect = SVE
    end type SolverType

    ! setting%Step
    type StepType
        integer :: First
        integer :: Current
        integer :: Final
    end type StepType

    ! setting%Time
    type TimeType
        character(14) :: DateTimeStamp
        real(8) :: Dt ! s
        real(8) :: StartTime
        real(8) :: NextTime
        real(8) :: EndTime
    end type TimeType

    ! setting%ZeroValue
    type ZerovalueType
        logical :: UseZeroValues = .true.
        real(8) :: Area = 1.0e-6 ! m^2
        real(8) :: Depth = 1.0e-4 ! m
        real(8) :: Topwidth = 1.0e-4 ! m
        real(8) :: Volume = 1.0e-6 ! m^3
    end type ZerovalueType

    ! -
    ! --

    ! First Level Type (setting)
    type settingType
        type(ACmethodType) :: ACmethod
        type(AdjustType) :: Adjust
        type(ConstantType) :: Constant ! Constants
        type(DebugOutType) :: DebugOut ! control of debougout files
        type(EpsilonType) :: Eps ! epsilons used to provide bandwidth for comparisons
        type(LimiterType) :: Limiter ! maximum and minimum limiters
        type(SmallVolumeType) :: SmallVolume ! controls for small volumes
        type(SolverType) :: Solver ! switch for solver
        type(StepType) :: Step ! controls over simulation time stepping
        type(TimeType) :: Time ! controls of time step
        type(ZeroValueType) :: ZeroValue ! finite values to represent small or negative values
    end type settingType

contains

    subroutine load_settings(fpath)
        character(len=:), intent(in) :: fpath
        type(json_file) :: json
    end subroutine load_settings
end module settings_definition