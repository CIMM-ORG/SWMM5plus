module define_settings

    use json_module
    use define_keys
    use define_globals
    use utility_string, only: util_lower_case

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

    !% setting%ACmethod%ImplicitCoef
    !% coefficients of the implicit stencil
    type ACmethodImplicitCoefType
        real(8) :: a1 = +threehalfR
        real(8) :: a2 = -twoR
        real(8) :: a3 = +onehalfR
    end type ACmethodImplicitCoefType

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
        logical :: Apply = .false.
        real(8) :: Coef = 1.0
        integer :: Approach = vshape
    end type AdjustHeadType

    ! setting%Adjust%WidthDepth
    type AdjustWidthDepthType
        logical :: Apply = .true.
    end type AdjustWidthDepthType

    ! setting%Limiter%BC
    type LimiterBCType
        logical :: UseInflowLimiter = .true.
        integer :: Approach = FroudeNumber
        ! max value of Fr at inflow
        real(8) :: FroudeInflowMaximum = 1.5
    end type LimiterBCType

    ! setting%Limiter%Flowrate
    type LimiterFlowrateType
        logical :: UseFaceVolumeTransport = .true.
        ! Fraction of usptream volume that can be
        ! transported in on time step
        real(8) :: FaceVolumeTransport = 0.5
    end type LimiterFlowrateType

    ! setting%Limiter%InterpWeight
    type LimiterInterpWeightType
        real(8) :: Maximum = 1e6
        real(8) :: Minimum = 1e-6
    end type LimiterInterpWeightType

    ! setting%Limiter%Velocity
    type LimiterVelocityType
        logical :: UseLimitMax = .true.
        real(8) :: Maximum = 10.0 ! m/s
    end type LimiterVelocityType

    ! setting%Limiter%ArraySize
    type LimiterArraySizeType
        integer :: TemporalInflows = 10
        integer :: TotallInflows = 50
    end type LimiterArraySizeType

    ! setting%Time%Hydrology or Hydraulics
    type HydrologyHydraulicsTimeType
        real(8) :: Dt = 100  !% time step (seconds)
        real(8) :: timeNow = 0.0 !% the time at the start of the present step
        real(8) :: timeNext = 0.0 !% the time af the end of the present step
        real(8) :: timeFinal = 1.0  !% the time this loop ends
        integer :: stepNow = 0 !% the present time step
        integer :: stepNext = 1 !% the next time step
        integer :: stepFinal = 1 !% the maximum number of steps allowed in this loop
    end type HydrologyHydraulicsTimeType  

    type WeirConstantType
        real(8) :: WeirExponent
        real(8) :: WeirContractionFactor
        real(8) :: SideFlowWeirCrestExponent
        real(8) :: VillemonteCorrectionExponent
    endtype WeirConstantType

    ! setting%Debug%File
    type DebugFileYNType
        logical :: adjust           = .false.
        logical :: define_globals   = .false.
        logical :: define_indexes   = .false.
        logical :: define_keys      = .false.
        logical :: define_settings  = .false.
        logical :: define_types     = .false.
        logical :: diagnostic_elements = .false.
        logical :: discretization   = .false.
        logical :: face             = .false.
        logical :: finalization     = .false.
        logical :: geometry         = .false.
        logical :: initial_condition = .false.
        logical :: initialization   = .false.
        logical :: jump             = .false.
        logical :: lowlevel_rk2     = .false.
        logical :: network_define   = .false.
        logical :: orifice_elements = .false.
        logical :: rectangular_channel = .false.
        logical :: runge_kutta2     = .false.
        logical :: pack_mask_arrays = .false.
        logical :: partitioning     = .false.
        logical :: pump_elements    = .false.
        logical :: interface        = .false.
        logical :: timeloop         = .false.
        logical :: update           = .false.
        logical :: utility          = .false.
        logical :: utility_allocate = .false.
        logical :: utility_deallocate = .false.
        logical :: utility_array    = .false.
        logical :: utility_datetime = .false.
        logical :: utility_string   = .false.
        logical :: weir_elements    = .false.
    end type DebugFileYNType

    ! setting%Debug%FileGroup
    type DebugFileGroupYNType
        logical :: all              = .false.
        logical :: definitions      = .false.
        logical :: finalization     = .false.
        logical :: geometry         = .false.
        logical :: initialization   = .false.
        logical :: interface        = .false.
        logical :: timeloop         = .false.
        logical :: utility          = .false.
    end type DebugFileGroupYNType

    ! -
    ! --
    ! ---

    ! Second Level Types

    ! -
    ! --

    ! setting%ACmethodType
    type ACmethodType
        real(8) :: dtau = 1.0
        type(ACmethodAnomalyType) :: Anomaly
        type(ACmethodImplicitCoefType) :: ImplicitCoef
        type(ACmethodCFLType) :: CFL
        type(ACmethodCelerityType) :: Celerity
        type(ACmethodConvergenceType) :: Convergence
        type(ACmethodIterType) :: Iter
        type(ACmethodSwitchType) :: Switch
    end type ACmethodType

    ! setting%Adjust
    type AdjustType
        type(AdjustFlowrateType)   :: Flowrate
        type(AdjustHeadType)       :: Head
        type(AdjustWidthDepthType) :: WidthDepth
    end type AdjustType

    ! setting%BC
    type BCPropertiesType
        integer :: BCSlots = 10
    end type BCPropertiesType

    ! setting%Constant
    type ConstantType
        real(8) :: gravity = 9.81 ! m^2/s
    end type ConstantType

    !% setting%Discretization
    type DiscretizationType
        real(8) :: NominalElemLength  = 10.0
        real(8) :: LinkShortingFactor = 0.33
    end type DiscretizationType

    ! setting%Eps
    type EpsilonType
        ! +- small non-dimensional range for hyd jump discrimination
        real(8) :: FroudeJump = 0.1
        ! Fractional increase in depth under froude limitation
        real(8) :: InflowDepthIncreaseFroudeLimit = 0.1
    end type EpsilonType

    ! setting%Limiter
    type LimiterType
        type(LimiterBCType)           :: BC
        type(LimiterFlowrateType)     :: Flowrate
        type(LimiterInterpWeightType) :: InterpWeight
        type(LimiterVelocityType)     :: Velocity
        type(LimiterArraySizeType)    :: ArraySize
    end type LimiterType

    type OrificeType
        real(8) :: SharpCrestedWeirCoefficient
        real(8) :: TransverseWeirExponent
        real(8) :: VillemonteCorrectionExponent
    end type OrificeType

    !% setting%Partitioning
    type PartitioningType
        integer :: PartitioningMethod = BLink
    endtype PartitioningType

    !% setting%Simulation
    type SimulationType
        logical :: useHydrology = .true.
        logical :: useHydraulics = .false.
    end type SimulationType

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
        real(8), dimension(2) :: crk2 = [0.5, 1.0]
        integer :: MomentumSourceMethod = T00
        logical :: PreissmanSlot = .true.
        integer :: SolverSelect = ETM_AC
        real(8) :: SwitchFractionDn = 0.8
        real(8) :: SwitchFractionUp = 0.9
    end type SolverType

    !% REMOVED 20210607 brh -- rolled into setting%Time
    ! ! setting%Step
    ! type StepType
    !     integer :: First
    !     integer :: Current
    !     integer :: Final
    ! end type StepType
    !%  setting%TestCase

    type TestCaseType
        logical       :: UseTestCase = .false.
        character(64) :: TestName
    end type TestCaseType

    ! setting%Time
    type TimeType
        character(14) :: DateTimeStamp
        ! real(8) :: Dt ! s !% REMOVED 20210607 brh -- rolled into setting%Time
        real(8) :: StartTime
        ! real(8) :: NextTime !% REMOVED 20210607 brh -- rolled into setting%Time
        real(8) :: EndTime
        type(HydrologyHydraulicsTimeType) :: Hydrology
        type(HydrologyHydraulicsTimeType) :: Hydraulics     
    end type TimeType

    type WeirType
        type(WeirConstantType) :: Transverse
        type(WeirConstantType) :: SideFlow
        type(WeirConstantType) :: VNotch
        type(WeirConstantType) :: Trapezoidal
    end type WeirType

    !% setting%VariableDT
    type VariableDTType
        logical :: Apply = .true.
        real(8) :: CFL_hi_max = 0.7
        real(8) :: CFL_target = 0.45
        real(8) :: CFL_lo_max = 0.3
        real(8) :: decreaseFactor = 0.8
        real(8) :: increaseFactor = 1.2
        integer :: NstepsForCheck = 10
        integer :: LastCheckStep = 0
    end type VariableDTType

    ! setting%ZeroValue
    ! Note that zerovalue.depth = zerovalue.area/zerovalue.topwidth
    type ZerovalueType
        logical :: UseZeroValues = .true.
        real(8) :: Area = 1.0e-6 ! m^2
        real(8) :: Depth = 1.0e-4 ! m
        real(8) :: Topwidth = 1.0e-2 ! m
        real(8) :: Volume = 1.0e-4 ! m^3
    end type ZerovalueType

    !% setting%Debug
    type DebugType
        logical :: Tests = .false.
        type(DebugFileYNType) :: File
        type(DebugFileGroupYNType) :: FileGroup
    end type DebugType

    !% setting%Paths
    type PathType
        character(len=256) :: project ! project path
        character(len=256) :: setting = "" ! path to settings JSON file
        character(len=256) :: inp = "" ! path to SWMM input (.inp) file
        character(len=256) :: rpt ! path to SWMM report (.rpt) file
        character(len=256) :: out ! path to SWMM output (.out) file
    end type PathType

    ! -
    ! --

    ! First Level Type (setting)
    type settingType
        type(ACmethodType)       :: ACmethod
        type(AdjustType)         :: Adjust
        type(BCPropertiesType)   :: BC
        type(ConstantType)       :: Constant ! Constants
        type(DiscretizationType) :: Discretization
        type(EpsilonType)        :: Eps ! epsilons used to provide bandwidth for comparisons
        type(LimiterType)        :: Limiter ! maximum and minimum limiters
        type(OrificeType)        :: Orifice 
        type(PartitioningType)   :: Partitioning
        type(SimulationType)     :: Simulation    
        type(SmallVolumeType)    :: SmallVolume ! controls for small volumes
        type(SolverType)         :: Solver ! switch for solver
    !rm 20210607 brh    type(StepType)           :: Step ! controls over simulation time stepping
        type(TimeType)           :: Time ! controls of time step
        type(VariableDTType)     :: VariableDT
        type(WeirType)           :: Weir
        type(ZeroValueType)      :: ZeroValue ! finite values to represent small or negative values
        type(TestCaseType)       :: TestCase
        type(PathType)           :: Paths
        type(DebugType)          :: Debug
        logical                  :: Verbose
        logical                  :: Warning = .true.
    end type settingType

    type(settingType), target :: setting

contains

    subroutine def_load_settings(fpath)
    !-----------------------------------------------------------------------------
    !
    ! Description:
    !
    !
    ! Method:
    !
    !
    !-----------------------------------------------------------------------------
        character(len=254), intent(in) :: fpath
        character(kind=json_CK, len=:), allocatable :: c
        real(8) :: real_value
        integer :: integer_value
        logical :: logical_value
        logical :: found
        type(json_file) :: json

        character(64) :: subroutine_name = 'def_load_settings'

        if (setting%Debug%File%define_settings) print *, '*** enter ', subroutine_name

        ! ---  Define paths
        call getcwd(setting%Paths%project)
        setting%Paths%setting = trim(setting%Paths%project) // '/definitions/settings.json'

        call json%initialize()
        call json%load(filename = fpath)

        ! Load ACmethod Settings
        call json%get('ACmethod.dtau', real_value, found)
        setting%ACmethod%dtau = real_value
        if (.not. found) stop 10
        call json%get('ACmethod.Anomaly.DensityLowCutoff', real_value, found)
        setting%ACmethod%Anomaly%DensityLowCutoff = real_value
        if (.not. found) stop 20
        call json%get('ACmethod.Anomaly.FullPipeFactor', real_value, found)
        setting%ACmethod%Anomaly%FullPipeFactor = real_value
        if (.not. found) stop 30
        call json%get('ACmethod.Anomaly.OpenPipeFactor', real_value, found)
        setting%ACmethod%Anomaly%OpenPipeFactor = real_value
        if (.not. found) stop 40
        call json%get('ACmethod.Anomaly.UseDensityCorrection', logical_value, found)
        setting%ACmethod%Anomaly%UseDensityCorrection = logical_value
        if (.not. found) stop 50
        call json%get('ACmethod.Anomaly.DensityHighCutoff', real_value, found)
        setting%ACmethod%Anomaly%DensityHighCutoff = real_value
        if (.not. found) stop 60

        ! Load implicit stencil coefficients
        call json%get('ACmethod.ImplicitCoef.a1', real_value, found)
        setting%ACmethod%ImplicitCoef%a1 = real_value
        if (.not. found) stop 61
        call json%get('ACmethod.ImplicitCoef.a2', real_value, found)
        setting%ACmethod%ImplicitCoef%a2 = real_value
        if (.not. found) stop 62
        call json%get('ACmethod.ImplicitCoef.a3', real_value, found)
        setting%ACmethod%ImplicitCoef%a3 = real_value
        if (.not. found) stop 63

        ! Load CFL Settings
        call json%get('ACmethod.CFL.CFLmax', real_value, found)
        setting%ACmethod%CFL%CFLmax = real_value
        if (.not. found) stop 70
        call json%get('ACmethod.CFL.CFLmax', real_value, found)
        setting%ACmethod%CFL%CFLmax = real_value
        if (.not. found) stop 80

        ! Load Celerity Settings
        call json%get('ACmethod.Celerity.RC', real_value, found)
        setting%ACmethod%Celerity%RC = real_value
        if (.not. found) stop 90

        ! Load Convergence Settings
        call json%get('ACmethod.Convergence.Habsolute', real_value, found)
        setting%ACmethod%Convergence%Habsolute = real_value
        if (.not. found) stop 100
        call json%get('ACmethod.Convergence.Hrelative', real_value, found)
        setting%ACmethod%Convergence%Hrelative = real_value
        if (.not. found) stop 110
        call json%get('ACmethod.Convergence.Qabsolute', real_value, found)
        setting%ACmethod%Convergence%Qabsolute = real_value
        if (.not. found) stop 120
        call json%get('ACmethod.Convergence.Qrelative', real_value, found)
        setting%ACmethod%Convergence%Qrelative = real_value
        if (.not. found) stop 130

        ! Load Iter Settings
        call json%get('ACmethod.Iter.Firststep', integer_value, found)
        setting%ACmethod%Iter%Firststep = integer_value
        if (.not. found) stop 140
        call json%get('ACmethod.Iter.Max', integer_value, found)
        setting%ACmethod%Iter%Max = integer_value
        if (.not. found) stop 150
        call json%get('ACmethod.Iter.Min', integer_value, found)
        setting%ACmethod%Iter%Min = integer_value
        if (.not. found) stop 160

        ! Load Switch Settings
        call json%get('ACmethod.Switch.Area', real_value, found)
        setting%ACmethod%Switch%Area = real_value
        if (.not. found) stop 170
        call json%get('ACmethod.Switch.Buffer', real_value, found)
        setting%ACmethod%Switch%Buffer = real_value
        if (.not. found) stop 180
        call json%get('ACmethod.Switch.Depth', real_value, found)
        setting%ACmethod%Switch%Depth = real_value
        if (.not. found) stop 190

        ! Load Adjust Settings
        call json%get('Adjust.Flowrate.Apply', logical_value, found)
        setting%Adjust%Flowrate%Apply = logical_value
        if (.not. found) stop 200
        call json%get('Adjust.Flowrate.approach', c, found)
        call util_lower_case(c)
        if (c == 'vshape') then
            setting%Adjust%Flowrate%approach = vshape
        else
            print *, "Error, Adjust.Flowrate.approach not compatible"
            stop 210
        end if
        call json%get('Adjust.Flowrate.Coef', real_value, found)
        setting%Adjust%Flowrate%Coef = real_value
        if (.not. found) stop 220

        call json%get('Adjust.Head.Apply', logical_value, found)
        setting%Adjust%Head%Apply = logical_value
        if (.not. found) stop 230
        call json%get('Adjust.Head.approach', c, found)
        call util_lower_case(c)
        if (c == 'vshape_surcharge_only') then
            setting%Adjust%Head%approach = vshape_surcharge_only
        else
            print *, "Error, Adjust.Head.approach not compatible"
            stop 240
        end if
        call json%get('Adjust.Head.Coef', real_value, found)
        setting%Adjust%Head%Coef = real_value
        if (.not. found) stop 250

        ! Load BC Settings
        call json%get('BC.BCSlots', real_value, found)
        setting%BC%BCslots = real_value
        if (.not. found) stop 260

        ! Load Constant Settings
        call json%get('Constant.gravity', real_value, found)
        setting%Constant%gravity = real_value
        if (.not. found) stop 270

        ! For element length adjustment
        call json%get("Discretization.NominalElemLength", real_value, found)
        setting%Discretization%NominalElemLength = real_value
        if (.not. found) stop 280
        call json%get("Discretization.LinkShortingFactor", real_value, found)
        setting%Discretization%LinkShortingFactor = real_value
        if (.not. found) stop 290

        ! Load Eps Settings
        call json%get('Eps.FroudeJump', real_value, found)
        setting%Eps%FroudeJump = real_value
        if (.not. found) stop 300
        call json%get('Eps.InflowDepthIncreaseFroudeLimit', real_value, found)
        setting%Eps%InflowDepthIncreaseFroudeLimit = real_value
        if (.not. found) stop 310

        ! Load Limiter Settings
        call json%get('Limiter.BC.approach', c, found)
        call util_lower_case(c)
        if (c == 'froudenumber') then
            setting%Limiter%BC%approach = FroudeNumber
        else
            print *, "Error, Limiter.BC.approach not compatible"
            stop
            stop
        end if
        if (.not. found) stop 320
        call json%get('Limiter.BC.FroudeInflowMaximum', real_value, found)
        setting%Limiter%BC%FroudeInflowMaximum = real_value
        if (.not. found) stop 330
        call json%get('Limiter.BC.UseInflowLimiter', logical_value, found)
        setting%Limiter%BC%UseInflowLimiter = logical_value
        if (.not. found) stop 340

        call json%get('Limiter.Flowrate.FaceVolumeTransport', real_value, found)
        setting%Limiter%Flowrate%FaceVolumeTransport = real_value
        if (.not. found) stop 350
        call json%get('Limiter.Flowrate.UseFaceVolumeTransport', logical_value, found)
        setting%Limiter%Flowrate%UseFaceVolumeTransport = logical_value
        if (.not. found) stop 360

        call json%get('Limiter.InterpWeight.Maximum', real_value, found)
        setting%Limiter%InterpWeight%Maximum = real_value
        if (.not. found) stop 370
        call json%get('Limiter.InterpWeight.Minimum', real_value, found)
        setting%Limiter%InterpWeight%Minimum = real_value
        if (.not. found) stop 380

        call json%get('Limiter.Velocity.Maximum', real_value, found)
        setting%Limiter%Velocity%Maximum = real_value
        if (.not. found) stop 390
        call json%get('Limiter.Velocity.UseLimitMax', logical_value, found)
        setting%Limiter%Velocity%UseLimitMax = logical_value
        if (.not. found) stop 400


        ! Load BIPQuick settings
        call json%get('Partitioning.PartitioningMethod', c, found)
        call util_lower_case(c)
        if (c == 'default') then
            setting%Partitioning%PartitioningMethod = Default
        else if (c == 'bquick') then
            setting%Partitioning%PartitioningMethod = BQuick
        else if (c == 'random') then
            setting%Partitioning%PartitioningMethod = Random
        else if (c == 'blink') then
            setting%Partitioning%PartitioningMethod = BLink
        else
            print *, "Error, the setting '" // trim(c) // "' is not supported for PartitioningMethod"
            stop 410
        end if
        if (.not. found) stop 420

        call json%get("Simulation.useHydrology", logical_value, found)
        setting%Simulation%useHydrology = logical_value
        if (.not. found) stop 430
        call json%get("Simulation.useHydraulics", logical_value, found)
        setting%Simulation%useHydraulics = logical_value
        if (.not. found) stop 440

        ! Load SmallVolume Settings
        call json%get('SmallVolume.DepthCutoff', real_value, found)
        setting%SmallVolume%DepthCutoff = real_value
        if (.not. found) stop 450
        call json%get('SmallVolume.ManningsN', real_value, found)
        setting%SmallVolume%ManningsN = real_value
        if (.not. found) stop 460
        call json%get('SmallVolume.MinimumArea', real_value, found)
        setting%SmallVolume%MinimumArea = real_value
        if (.not. found) stop 470
        call json%get('SmallVolume.MinimumHydRadius', real_value, found)
        setting%SmallVolume%MinimumHydRadius = real_value
        if (.not. found) stop 480
        call json%get('SmallVolume.MinimumPerimeter', real_value, found)
        setting%SmallVolume%MinimumPerimeter = real_value
        if (.not. found) stop 490
        call json%get('SmallVolume.MinimumTopwidth', real_value, found)
        setting%SmallVolume%MinimumTopwidth = real_value
        if (.not. found) stop 500
        call json%get('SmallVolume.UseSmallVolumes', logical_value, found)
        setting%SmallVolume%UseSmallVolumes = logical_value
        if (.not. found) stop 510

        ! Load Solver Settings
        call json%get('Solver.crk2', real_value, found)
        setting%Solver%crk2 = real_value
        if (.not. found) stop 520
        call json%get('Solver.MomentumSourceMethod', c, found)
        call util_lower_case(c)
        if (c == 't00') then
            setting%Solver%MomentumSourceMethod = T00
        else if (c == 't10') then
            setting%Solver%MomentumSourceMethod = T10
        else if (c == 't20') then
            setting%Solver%MomentumSourceMethod = T20
        else
            print *, "Error, the setting '" // trim(c) // "' is not supported for MomentumSourceMethod"
            stop 530
        end if
        if (.not. found) stop 540
        call json%get('Solver.PreissmanSlot', logical_value, found)
        setting%Solver%PreissmanSlot = logical_value
        if (.not. found) stop 550
        call json%get('Solver.SolverSelect', c, found)
        call util_lower_case(c)
        if (c == 'etm') then
            setting%Solver%SolverSelect = ETM
        else if (c == 'etm_ac') then
            setting%Solver%SolverSelect = ETM_AC
        else if (c == 'ac') then
            setting%Solver%SolverSelect = AC
        else
            print *, "Error, the setting '" // trim(c) // "' is not supported for SolverSelect"
            stop 560
        end if
        if (.not. found) stop 570
        call json%get('Solver.SwitchFractionDn', real_value, found)
        setting%Solver%SwitchFractionDn = real_value
        if (.not. found) stop 580
        call json%get('Solver.SwitchFractionUp', real_value, found)
        setting%Solver%SwitchFractionUp = real_value
        if (.not. found) stop 590       

        ! Load Step Settings
        !rm 20210607 brh call json%get('Step.Current', real_value, found)
        !rm 20210607 brh setting%Step%Current = real_value
        !rm 20210607 brh if (.not. found) stop 48
        !rm 20210607 brh call json%get('Step.Final', real_value, found)
        !rm 20210607 brh setting%Step%Final = real_value
        !rm 20210607 brh if (.not. found) stop 49
        !rm 20210607 brh call json%get('Step.First', real_value, found)
        !rm 20210607 brh setting%Step%First = real_value
        !rm 20210607 brh if (.not. found) stop 50

        ! Load Time Settings
        call json%get('Time.DateTimeStamp', c, found)
        setting%Time%DateTimeStamp = c
        if (.not. found) stop 600
        !rm 20210607 brh call json%get('Time.Dt', real_value, found)
        !rm 20210607 brh setting%Time%Dt = real_value
        !rm 20210607 brh if (.not. found) stop 52
        call json%get('Time.EndTime', real_value, found)
        setting%Time%EndTime = real_value
        !rm 20210607 brh if (.not. found) stop 49
        !rm 20210607 brh call json%get('Time.NextTime', real_value, found)
        !rm 20210607 brh setting%Time%NextTime = real_value
        if (.not. found) stop 610
        call json%get('Time.StartTime', real_value, found)
        setting%Time%StartTime = real_value
        if (.not. found) stop 620
        call json%get('Time.Hydrology.Dt', real_value, found)
        setting%Time%Hydrology%Dt = real_value
        if (.not. found) stop 630
        call json%get('Time.Hydrology.timeNow', real_value, found)
        setting%Time%Hydrology%timeNow = real_value
        if (.not. found) stop 640
        call json%get('Time.Hydrology.timeNext', real_value, found)
        setting%Time%Hydrology%timeNext = real_value
        if (.not. found) stop 650
        call json%get('Time.Hydrology.timeFinal', real_value, found)
        setting%Time%Hydrology%timeFinal = real_value
        if (.not. found) stop 660   
        call json%get('Time.Hydrology.stepNow', integer_value, found)
        setting%Time%Hydrology%stepNow = integer_value
        if (.not. found) stop 670
        call json%get('Time.Hydrology.stepNext', integer_value, found)
        setting%Time%Hydrology%stepNext = integer_value
        if (.not. found) stop 680      
        call json%get('Time.Hydrology.stepFinal', integer_value, found)
        setting%Time%Hydrology%stepFinal = integer_value
        if (.not. found) stop 690   
        call json%get('Time.Hydraulics.Dt', real_value, found)
        setting%Time%Hydraulics%Dt = real_value
        if (.not. found) stop 700
        call json%get('Time.Hydraulics.timeNow', real_value, found)
        setting%Time%Hydraulics%timeNow = real_value
        if (.not. found) stop 710
        call json%get('Time.Hydraulics.timeNext', real_value, found)
        setting%Time%Hydraulics%timeNext = real_value
        if (.not. found) stop 720
        call json%get('Time.Hydraulics.timeFinal', real_value, found)
        setting%Time%Hydraulics%timeFinal = real_value
        if (.not. found) stop 730 
        call json%get('Time.Hydraulics.stepNow', integer_value, found)
        setting%Time%Hydraulics%stepNow = integer_value
        if (.not. found) stop 740
        call json%get('Time.Hydraulics.stepNext', integer_value, found)
        setting%Time%Hydraulics%stepNext = integer_value
        if (.not. found) stop 750    
        call json%get('Time.Hydraulics.stepFinal', integer_value, found)
        setting%Time%Hydraulics%stepFinal = integer_value
        if (.not. found) stop 760  

        call json%get("Weir.Transverse.WeirExponent", real_value, found)
        setting%Weir%Transverse%WeirExponent = real_value
        if (.not. found) stop 7601
        call json%get("Weir.Transverse.WeirContractionFactor", real_value, found)
        setting%Weir%Transverse%WeirContractionFactor = real_value
        if (.not. found) stop 7602
        call json%get("Weir.Transverse.SideFlowWeirCrestExponent", real_value, found)
        setting%Weir%Transverse%SideFlowWeirCrestExponent = real_value
        if (.not. found) stop 7603
        call json%get("Weir.Transverse.VillemonteCorrectionExponent", real_value, found)
        setting%Weir%Transverse%VillemonteCorrectionExponent = real_value
        if (.not. found) stop 7604

        call json%get("Weir.SideFlow.WeirExponent", real_value, found)
        setting%Weir%SideFlow%WeirExponent = real_value
        if (.not. found) stop 7611
        call json%get("Weir.SideFlow.WeirContractionFactor", real_value, found)
        setting%Weir%SideFlow%WeirContractionFactor = real_value
        if (.not. found) stop 7612
        call json%get("Weir.SideFlow.SideFlowWeirCrestExponent", real_value, found)
        setting%Weir%SideFlow%SideFlowWeirCrestExponent = real_value
        if (.not. found) stop 7613
        call json%get("Weir.SideFlow.VillemonteCorrectionExponent", real_value, found)
        setting%Weir%SideFlow%VillemonteCorrectionExponent = real_value
        if (.not. found) stop 7614

        call json%get("Weir.VNotch.WeirExponent", real_value, found)
        setting%Weir%VNotch%WeirExponent = real_value
        if (.not. found) stop 7621
        call json%get("Weir.VNotch.WeirContractionFactor", real_value, found)
        setting%Weir%VNotch%WeirContractionFactor = real_value
        if (.not. found) stop 7622
        call json%get("Weir.VNotch.SideFlowWeirCrestExponent", real_value, found)
        setting%Weir%VNotch%SideFlowWeirCrestExponent = real_value
        if (.not. found) stop 7623
        call json%get("Weir.VNotch.VillemonteCorrectionExponent", real_value, found)
        setting%Weir%VNotch%VillemonteCorrectionExponent = real_value
        if (.not. found) stop 7624


        call json%get("Weir.Trapezoidal.WeirExponent", real_value, found)
        setting%Weir%Trapezoidal%WeirExponent = real_value
        if (.not. found) stop 7631
        call json%get("Weir.Trapezoidal.WeirContractionFactor", real_value, found)
        setting%Weir%Trapezoidal%WeirContractionFactor = real_value
        if (.not. found) stop 7632
        call json%get("Weir.Trapezoidal.SideFlowWeirCrestExponent", real_value, found)
        setting%Weir%Trapezoidal%SideFlowWeirCrestExponent = real_value
        if (.not. found) stop 7633
        call json%get("Weir.Trapezoidal.VillemonteCorrectionExponent", real_value, found)
        setting%Weir%Trapezoidal%VillemonteCorrectionExponent = real_value
        if (.not. found) stop 7634

        !% load variable time step settings
        call json%get("VariableDT.Apply", logical_value, found)
        setting%VariableDT%Apply = logical_value
        if (.not. found) stop 770
        call json%get("VariableDT.CFL_hi_max", real_value, found)
        setting%VariableDT%CFL_hi_max = real_value
        if (.not. found) stop 780
        call json%get("VariableDT.CFL_target", real_value, found)
        setting%VariableDT%CFL_target = real_value
        if (.not. found) stop 790
        call json%get("VariableDT.CFL_lo_max", real_value, found)
        setting%VariableDT%CFL_lo_max = real_value
        if (.not. found) stop 800
        call json%get("VariableDT.decreaseFactor", real_value, found)
        setting%VariableDT%decreaseFactor = real_value
        if (.not. found) stop 810
        call json%get("VariableDT.increaseFactor", real_value, found)
        setting%VariableDT%increaseFactor = real_value
        if (.not. found) stop 820
        call json%get("VariableDT.NstepsForCheck", integer_value, found)
        setting%VariableDT%NstepsForCheck = integer_value
        if (.not. found) stop 830

        ! Load ZeroValue Settings
        call json%get('ZeroValue.Area', real_value, found)
        setting%ZeroValue%Area = real_value
        if (.not. found) stop 840
        call json%get('ZeroValue.Depth', real_value, found)
        setting%ZeroValue%Depth = real_value
        if (.not. found) stop 850
        call json%get('ZeroValue.Topwidth', real_value, found)
        setting%ZeroValue%Topwidth = real_value
        if (.not. found) stop 860
        call json%get('ZeroValue.UseZeroValues', logical_value, found)
        setting%ZeroValue%UseZeroValues = logical_value
        if (.not. found) stop 870
        call json%get('ZeroValue.Volume', real_value, found)
        setting%ZeroValue%Volume = real_value
        if (.not. found) stop 880

        ! Load Debug Settings
        call json%get('Debug.File.adjust', logical_value, found)
        setting%Debug%File%adjust = logical_value
        if (.not. found) stop 885
        call json%get('Debug.File.define_globals', logical_value, found)
        setting%Debug%File%define_globals = logical_value
        if (.not. found) stop 890
        call json%get('Debug.File.define_indexes', logical_value, found)
        setting%Debug%File%define_indexes = logical_value
        if (.not. found) stop 900
        call json%get('Debug.File.define_keys', logical_value, found)
        setting%Debug%File%define_keys = logical_value
        if (.not. found) stop 910
        call json%get('Debug.File.define_settings', logical_value, found)
        setting%Debug%File%define_settings = logical_value
        if (.not. found) stop 920
        call json%get('Debug.File.define_types', logical_value, found)
        setting%Debug%File%define_types = logical_value
        if (.not. found) stop 930
        call json%get('Debug.File.diagnostic_elements', logical_value, found)
        setting%Debug%File%diagnostic_elements = logical_value
        if (.not. found) stop 935
        call json%get('Debug.File.discretization', logical_value, found)
        setting%Debug%File%discretization = logical_value
        if (.not. found) stop 940
        call json%get('Debug.File.face', logical_value, found)
        setting%Debug%File%face = logical_value
        if (.not. found) stop 942
        call json%get('Debug.File.geometry', logical_value, found)
        setting%Debug%File%geometry = logical_value
        if (.not. found) stop 943        
        call json%get('Debug.File.interface', logical_value, found)
        setting%Debug%File%interface = logical_value
        if (.not. found) stop 945
        call json%get('Debug.File.initial_condition', logical_value, found)
        setting%Debug%File%initial_condition = logical_value
        if (.not. found) stop 945
        call json%get('Debug.File.initialization', logical_value, found)
        setting%Debug%File%initialization = logical_value
        if (.not. found) stop 950
        call json%get('Debug.File.jump', logical_value, found)
        setting%Debug%File%jump = logical_value
        if (.not. found) stop 953
        call json%get('Debug.File.lowlevel_rk2', logical_value, found)
        setting%Debug%File%lowlevel_rk2 = logical_value
        if (.not. found) stop 955
        call json%get('Debug.File.network_define', logical_value, found)
        setting%Debug%File%network_define = logical_value
        if (.not. found) stop 960
        call json%get('Debug.File.orifice_elements', logical_value, found)
        setting%Debug%File%orifice_elements = logical_value
        if (.not. found) stop 965
        call json%get('Debug.File.pack_mask_arrays', logical_value, found)
        setting%Debug%File%pack_mask_arrays = logical_value
        if (.not. found) stop 965
        call json%get('Debug.File.partitioning', logical_value, found)
        setting%Debug%File%partitioning = logical_value
        if (.not. found) stop 970
        call json%get('Debug.File.pump_elements', logical_value, found)
        setting%Debug%File%pump_elements = logical_value
        if (.not. found) stop 975
        call json%get('Debug.File.rectangular_channel', logical_value, found)
        setting%Debug%File%rectangular_channel = logical_value
        if (.not. found) stop 978
        call json%get('Debug.File.runge_kutta2', logical_value, found)
        setting%Debug%File%runge_kutta2 = logical_value
        if (.not. found) stop 980
        call json%get('Debug.File.timeloop', logical_value, found)
        setting%Debug%File%timeloop = logical_value
        if (.not. found) stop 982
        call json%get('Debug.File.update', logical_value, found)
        setting%Debug%File%update = logical_value
        if (.not. found) stop 985
        call json%get('Debug.File.utility_allocate', logical_value, found)
        setting%Debug%File%utility_allocate = logical_value
        if (.not. found) stop 990
        call json%get('Debug.File.utility_deallocate', logical_value, found)
        setting%Debug%File%utility_deallocate = logical_value
        if (.not. found) stop 1000
        call json%get('Debug.File.utility_array', logical_value, found)
        setting%Debug%File%utility_array = logical_value
        if (.not. found) stop 1010
        call json%get('Debug.File.utility_datetime', logical_value, found)
        setting%Debug%File%utility_datetime = logical_value
        if (.not. found) stop 1020
        call json%get('Debug.File.utility_string', logical_value, found)
        setting%Debug%File%utility_string = logical_value
        if (.not. found) stop 1030
        call json%get('Debug.File.utility', logical_value, found)
        setting%Debug%File%utility = logical_value
        call json%get('Debug.File.weir_elements', logical_value, found)
        setting%Debug%File%weir_elements = logical_value
        if (.not. found) stop 1035
        if (.not. found) stop 1040
        call json%get('Debug.FileGroup.all', logical_value, found)
        setting%Debug%FileGroup%all = logical_value
        if (.not. found) stop 1050
        call json%get('Debug.FileGroup.definitions', logical_value, found)
        setting%Debug%FileGroup%definitions = logical_value
        if (.not. found) stop 1060
        call json%get('Debug.FileGroup.finalization', logical_value, found)
        setting%Debug%FileGroup%finalization = logical_value
        if (.not. found) stop 1070
        call json%get('Debug.FileGroup.geometry', logical_value, found)
        setting%Debug%FileGroup%geometry = logical_value
        if (.not. found) stop 1075
        call json%get('Debug.FileGroup.initialization', logical_value, found)
        setting%Debug%FileGroup%initialization = logical_value
        if (.not. found) stop 1080
        call json%get('Debug.FileGroup.interface', logical_value, found)
        setting%Debug%FileGroup%interface = logical_value
        if (.not. found) stop 1090
        call json%get('Debug.FileGroup.timeloop', logical_value, found)
        setting%Debug%FileGroup%timeloop = logical_value
        if (.not. found) stop 1095
        call json%get('Debug.FileGroup.utility', logical_value, found)
        setting%Debug%FileGroup%interface = logical_value
        if (.not. found) stop 1100
        call def_update_debug_options()

        ! Load verbose or non-verbose run
        call json%get('Verbose', logical_value, found)
        setting%Verbose = logical_value
        if (.not. found) stop 1200

        call json%destroy()
        if (json%failed()) stop 1210

        if (setting%Debug%File%define_settings) print *, '*** leave ', subroutine_name
    end subroutine def_load_settings

    subroutine def_update_debug_options()
        if (setting%Debug%FileGroup%all) then
            setting%Debug%FileGroup%definitions = .true.
            setting%Debug%FileGroup%finalization = .true.
            setting%Debug%FileGroup%geometry = .true.
            setting%Debug%FileGroup%initialization = .true.
            setting%Debug%FileGroup%interface = .true.
            setting%Debug%FileGroup%timeloop  = .true.
            setting%Debug%FileGroup%utility = .true.
        end if
        if (setting%Debug%FileGroup%definitions) then
            setting%Debug%File%define_globals = .true.
            setting%Debug%File%define_indexes = .true.
            setting%Debug%File%define_keys = .true.
            setting%Debug%File%define_settings = .true.
            setting%Debug%File%define_types = .true.
        end if
        if (setting%Debug%FileGroup%finalization) then
            setting%Debug%File%finalization = .true.
        end if        
        if (setting%Debug%FileGroup%geometry) then
            setting%Debug%File%geometry = .true.
            setting%Debug%File%rectangular_channel = .true.
        end if
        if (setting%Debug%FileGroup%initialization) then
            setting%Debug%File%discretization = .true.
            setting%Debug%File%initial_condition = .true.
            setting%Debug%File%initialization = .true.
            setting%Debug%File%network_define = .true.
            setting%Debug%File%partitioning = .true.
            setting%Debug%File%pack_mask_arrays = .true.
        end if
        if (setting%Debug%FileGroup%interface) then
            setting%Debug%File%interface = .true.
        end if
        if (setting%Debug%FileGroup%timeloop) then
            setting%Debug%File%adjust = .true.
            setting%Debug%File%diagnostic_elements = .true.
            setting%Debug%File%face = .true.
            setting%Debug%File%jump = .true.
            setting%Debug%File%lowlevel_rk2 = .true.
            setting%Debug%File%orifice_elements = .true.
            setting%Debug%File%pump_elements = .true.
            setting%Debug%File%runge_kutta2 = .true.
            setting%Debug%File%timeloop = .true.
            setting%Debug%File%update = .true.
            setting%Debug%File%weir_elements = .true.
        endif    
        if (setting%Debug%FileGroup%utility) then
            setting%Debug%File%utility_allocate = .true.
            setting%Debug%File%utility_deallocate = .true.
            setting%Debug%File%utility_array = .true.
            setting%Debug%File%utility_datetime = .true.
            setting%Debug%File%utility_string = .true.
            setting%Debug%File%utility = .true.
        end if
    end subroutine def_update_debug_options
end module define_settings
