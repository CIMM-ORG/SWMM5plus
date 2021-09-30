module define_settings

    use json_module
    use define_keys
    use define_globals
    use utility_string, only: util_lower_case

    implicit none
    public

    !% -------------------------------------------------------------------------------
    !% Notes:
    !%    Settings are stored in the setting object which is public.
    !%    Settings can be defined in three ways:
    !%        1) via flags
    !%        2) via JSON file
    !%        3) via default values
    !%
    !%    If the settings are not defined via JSON file, default values
    !%    should be defined in the declaration of the derived type
    !%    associated to the setting. Settings flags (e.g., verbose -v)
    !%    overwrite any value defined via default or JSON values.
    !% -------------------------------------------------------------------------------

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

    ! setting%Output%CommandLine
    type CommandLineType
        logical :: quiet = .false.
        integer :: interval = 10
    end type CommandLineType

    ! setting%Limiter%BC
    type LimiterBCType
        logical :: UseInflowLimiter = .true.
        integer :: Approach = FroudeNumber
        ! max value of Fr at inflow
        real(8) :: FroudeInflowMaximum = 1.5
    end type LimiterBCType

    ! setting%Limiter%Channel
    type LimiterChannelType
        real(8) :: LargeDepthFactor = 10.0
    end type LimiterChannelType

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

    type LimiterDtType
        logical :: UseLimitMin = .true.
        real(8) :: Minimum     = 1e-6
    end type LimiterDtType

    type WeirConstantType
        real(8) :: WeirExponent
        real(8) :: WeirContractionFactor
        real(8) :: SideFlowWeirCrestExponent
        real(8) :: VillemonteCorrectionExponent
    endtype WeirConstantType

    ! setting%Debug%File
    type DebugFileYNType
        logical :: adjust              = .false.
        logical :: boundary_conditions = .false.
        logical :: c_library           = .false.
        logical :: define_globals      = .false.
        logical :: define_indexes      = .false.
        logical :: define_keys         = .false.
        logical :: define_settings     = .false.
        logical :: define_types        = .false.
        logical :: diagnostic_elements = .false.
        logical :: BIPquick            = .false.
        logical :: discretization      = .false.
        logical :: face                = .false.
        logical :: finalization        = .false.
        logical :: geometry            = .false.
        logical :: initial_condition   = .false.
        logical :: initialization      = .false.
        logical :: jump                = .false.
        logical :: lowlevel_rk2        = .false.
        logical :: network_define      = .false.
        logical :: orifice_elements    = .false.
        logical :: rectangular_channel = .false.
        logical :: trapezoidal_channel = .false.
        logical :: runge_kutta2        = .false.
        logical :: pack_mask_arrays    = .false.
        logical :: partitioning        = .false.
        logical :: pump_elements       = .false.
        logical :: interface           = .false.
        logical :: timeloop            = .false.
        logical :: update              = .false.
        logical :: utility             = .false.
        logical :: utility_allocate    = .false.
        logical :: utility_deallocate  = .false.
        logical :: utility_array       = .false.
        logical :: utility_datetime    = .false.
        logical :: utility_interpolate = .false.
        logical :: utility_output      = .false.
        logical :: utility_string      = .false.
        logical :: weir_elements       = .false.
        logical :: output              = .false.
    end type DebugFileYNType

    ! setting%Debug%FileGroup
    type DebugFileGroupYNType
        logical :: all              = .false.
        logical :: definitions      = .false.
        logical :: finalization     = .false.
        logical :: geometry         = .false.
        logical :: initialization   = .false.
        logical :: interface        = .false.
        logical :: output           = .false.
        logical :: timeloop         = .false.
        logical :: utility          = .false.
    end type DebugFileGroupYNType

    type TimeStepType
        real(8) :: Dt
        integer :: Step
    end type TimeStepType

    type RealTimeType
        integer :: EpochStartSeconds = 0
        integer :: EpochTimeLoopStartSeconds = 0
        integer :: EpochNowSeconds  = 0
    end type RealTimeType

    type CPUTimeType
        real (8) :: EpochStartSeconds = 0.0
        real (8) :: EpochNowSeconds  = 0.0
        real (8) :: EpochFinishSeconds  = 0.0
    end type CPUTimeType   

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
        integer :: slots = 10
        logical :: disableInterpolation = .false.
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

    !% setting%FaceInterp
    type FaceInterpType
        integer :: DownJBFaceInterp = Static
    end type FaceInterpType

    ! setting%Junction
    type JunctionType
        real(8) :: kFactor  = 0.0   !% default entrance/exit losses at junction branch (use 0.0 as needs debugging)
        real(8) :: HeadCoef = 1.0   !% junction branch head coef for diagnostic junction (must be > 0)
        real(8) :: CFLlimit = 0.5   !% limiter on CFL to control dynamic junction
        logical :: isDynamic = .true.
    end type JunctionType

    ! setting%Limiter
    type LimiterType
        type(LimiterBCType)           :: BC
        type(LimiterChannelType)      :: Channel
        type(LimiterFlowrateType)     :: Flowrate
        type(LimiterInterpWeightType) :: InterpWeight
        type(LimiterVelocityType)     :: Velocity
        type(LimiterArraySizeType)    :: ArraySize
        type(LimiterDtType)           :: Dt
    end type LimiterType

    type LinkType
        integer ::        DefaultInitDepthType = LinearlyVarying ! uniform, linear, exponential
        ! HACK - TODO file for properties of specific links
        character(512) :: PropertiesFile
    end type LinkType

    !% HACK - not defined in settings.json
    type OrificeType
        real(8) :: SharpCrestedWeirCoefficient
        real(8) :: TransverseWeirExponent
        real(8) :: VillemonteCorrectionExponent
    end type OrificeType

    !% setting%Partitioning
    type PartitioningType
        integer :: PartitioningMethod = BLink
    endtype PartitioningType

    !% setting%PreissmannSlot
    type PreissmannSlotType
        integer :: PreissmannSlotMethod = VariableSlot
        real(8) :: CelerityFactor = 1.0
    end type PreissmannSlotType

    !% setting%Simulation
    type SimulationType
        logical :: useHydrology = .true.
        logical :: useHydraulics = .true.
        integer :: useSWMMC = 0
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

    type TestCaseType
        logical       :: UseTestCase = .false.
        character(64) :: TestName
    end type TestCaseType

    ! setting%Time
    type TimeType
        type(TimeStepType) :: Hydraulics
        type(TimeStepType) :: Hydrology
        logical            :: matchHydrologyStep
        character(14)      :: DateTimeStamp
        integer            :: Step
        real(8)            :: Dt
        real(8)            :: DtTol
        real(8)            :: Start
        real(8)            :: Now
        real(8)            :: End
        real(8)            :: StartEpoch
        real(8)            :: EndEpoch
        type(RealTimeType) :: Real
        type(CPUTimeType)  :: CPU
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
        logical :: Input
        logical :: Output
    end type DebugType

    !% setting%Profile
    type ProfileType
        logical :: YN = .false.
        !logical :: Tests = .false.
        !type(ProfileFileYNType) :: File
        !type(ProfileFileGroupYNType) :: FileGroup
        !logical :: Input
        !logical :: Output
    end type ProfileType

    !% setting%Paths
    type PathType
        character(len=256) :: project ! project path
        character(len=256) :: setting = "" ! path to settings JSON file
        character(len=256) :: inp = "" ! path to SWMM input (.inp) file
        character(len=256) :: rpt ! path to SWMM report (.rpt) file
        character(len=256) :: out ! path to SWMM output (.out) file
    end type PathType

    !% setting%Output
    type OutputType
        logical :: report
        real(8) :: reportStartTime
        real(8) :: reportDt
        integer :: reportStep
        integer :: Slots = 20
        character(len=256) :: nodes_file = "node_input.csv"
        character(len=256) :: links_file = "link_input.csv"
        type(CommandLineType) :: CommandLine
    end type OutputType


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
        type(FaceInterpType)     :: FaceInterp ! Temporary: setting for face interpolation in downstream JB
        type(JunctionType)       :: Junction
        type(LimiterType)        :: Limiter ! maximum and minimum limiters
        type(LinkType)           :: Link
        type(OrificeType)        :: Orifice
        type(PartitioningType)   :: Partitioning
        type(PreissmannSlotType) :: PreissmannSlot
        type(SimulationType)     :: Simulation
        type(SmallVolumeType)    :: SmallVolume ! controls for small volumes
        type(SolverType)         :: Solver ! switch for solver
        type(TimeType)           :: Time ! controls of time step
        type(VariableDTType)     :: VariableDT
        type(WeirType)           :: Weir
        type(ZeroValueType)      :: ZeroValue ! finite values to represent small or negative values
        type(TestCaseType)       :: TestCase
        type(PathType)           :: Paths
        type(DebugType)          :: Debug
        type(ProfileType)        :: Profile
        type(OutputType)         :: Output
        logical                  :: Verbose
        logical                  :: Warning = .true.
    end type settingType

    type(settingType), target :: setting

contains

    subroutine def_load_settings()
    !%-----------------------------------------------------------------------------
    !% Description:
    !%    Loads setting values from external JSON file.
    !%-----------------------------------------------------------------------------
        character(kind=json_CK, len=:), allocatable :: c
        real(8) :: real_value
        integer :: integer_value
        logical :: logical_value
        logical :: found
        type(json_file) :: json

        character(64) :: subroutine_name = 'def_load_settings'

        if (setting%Debug%File%define_settings) &
            write(*,"(A,i5,A)") '*** enter ' // subroutine_name // " [Processor ", this_image(), "]"

        call json%initialize()
        call json%load(filename = trim(setting%paths%setting))

        ! Load ACmethod Settings
        call json%get('ACmethod.dtau', real_value, found)
        setting%ACmethod%dtau = real_value
        if (.not. found) stop "Error - setting " // 'ACmethod.dtau not found'
        call json%get('ACmethod.Anomaly.DensityLowCutoff', real_value, found)
        setting%ACmethod%Anomaly%DensityLowCutoff = real_value
        if (.not. found) stop "Error - setting " // 'ACmethod.Anomaly.DensityLowCutoff not found'
        call json%get('ACmethod.Anomaly.FullPipeFactor', real_value, found)
        setting%ACmethod%Anomaly%FullPipeFactor = real_value
        if (.not. found) stop "Error - setting " // 'ACmethod.Anomaly.FullPipeFactor not found'
        call json%get('ACmethod.Anomaly.OpenPipeFactor', real_value, found)
        setting%ACmethod%Anomaly%OpenPipeFactor = real_value
        if (.not. found) stop "Error - setting " // 'ACmethod.Anomaly.OpenPipeFactor not found'
        call json%get('ACmethod.Anomaly.UseDensityCorrection', logical_value, found)
        setting%ACmethod%Anomaly%UseDensityCorrection = logical_value
        if (.not. found) stop "Error - setting " // 'ACmethod.Anomaly.UseDensityCorrection not found'
        call json%get('ACmethod.Anomaly.DensityHighCutoff', real_value, found)
        setting%ACmethod%Anomaly%DensityHighCutoff = real_value
        if (.not. found) stop "Error - setting " // 'ACmethod.Anomaly.DensityHighCutoff not found'

        ! Load implicit stencil coefficients
        call json%get('ACmethod.ImplicitCoef.a1', real_value, found)
        setting%ACmethod%ImplicitCoef%a1 = real_value
        if (.not. found) stop "Error - setting " // 'ACmethod.ImplicitCoef.a1 not found'
        call json%get('ACmethod.ImplicitCoef.a2', real_value, found)
        setting%ACmethod%ImplicitCoef%a2 = real_value
        if (.not. found) stop "Error - setting " // 'ACmethod.ImplicitCoef.a2 not found'
        call json%get('ACmethod.ImplicitCoef.a3', real_value, found)
        setting%ACmethod%ImplicitCoef%a3 = real_value
        if (.not. found) stop "Error - setting " // 'ACmethod.ImplicitCoef.a3 not found'

        ! Load CFL Settings
        call json%get('ACmethod.CFL.CFLmax', real_value, found)
        setting%ACmethod%CFL%CFLmax = real_value
        if (.not. found) stop "Error - setting " // 'ACmethod.CFL.CFLmax not found'
        call json%get('ACmethod.CFL.CFLmax', real_value, found)
        setting%ACmethod%CFL%CFLmax = real_value
        if (.not. found) stop "Error - setting " // 'ACmethod.CFL.CFLmax not found'

        ! Load Celerity Settings
        call json%get('ACmethod.Celerity.RC', real_value, found)
        setting%ACmethod%Celerity%RC = real_value
        if (.not. found) stop "Error - setting " // 'ACmethod.Celerity.RC not found'

        ! Load Convergence Settings
        call json%get('ACmethod.Convergence.Habsolute', real_value, found)
        setting%ACmethod%Convergence%Habsolute = real_value
        if (.not. found) stop "Error - setting " // 'ACmethod.Convergence.Habsolute not found'
        call json%get('ACmethod.Convergence.Hrelative', real_value, found)
        setting%ACmethod%Convergence%Hrelative = real_value
        if (.not. found) stop "Error - setting " // 'ACmethod.Convergence.Hrelative not found'
        call json%get('ACmethod.Convergence.Qabsolute', real_value, found)
        setting%ACmethod%Convergence%Qabsolute = real_value
        if (.not. found) stop "Error - setting " // 'ACmethod.Convergence.Qabsolute not found'
        call json%get('ACmethod.Convergence.Qrelative', real_value, found)
        setting%ACmethod%Convergence%Qrelative = real_value
        if (.not. found) stop "Error - setting " // 'ACmethod.Convergence.Qrelative not found'

        ! Load Iter Settings
        call json%get('ACmethod.Iter.Firststep', integer_value, found)
        setting%ACmethod%Iter%Firststep = integer_value
        if (.not. found) stop "Error - setting " // 'ACmethod.Iter.Firststep not found'
        call json%get('ACmethod.Iter.Max', integer_value, found)
        setting%ACmethod%Iter%Max = integer_value
        if (.not. found) stop "Error - setting " // 'ACmethod.Iter.Max not found'
        call json%get('ACmethod.Iter.Min', integer_value, found)
        setting%ACmethod%Iter%Min = integer_value
        if (.not. found) stop "Error - setting " // 'ACmethod.Iter.Min not found'

        ! Load Switch Settings
        call json%get('ACmethod.Switch.Area', real_value, found)
        setting%ACmethod%Switch%Area = real_value
        if (.not. found) stop "Error - setting " // 'ACmethod.Switch.Area not found'
        call json%get('ACmethod.Switch.Buffer', real_value, found)
        setting%ACmethod%Switch%Buffer = real_value
        if (.not. found) stop "Error - setting " // 'ACmethod.Switch.Buffer not found'
        call json%get('ACmethod.Switch.Depth', real_value, found)
        setting%ACmethod%Switch%Depth = real_value
        if (.not. found) stop "Error - setting " // 'ACmethod.Switch.Depth not found'

        ! Load Adjust Settings
        call json%get('Adjust.Flowrate.Apply', logical_value, found)
        setting%Adjust%Flowrate%Apply = logical_value
        if (.not. found) stop "Error - setting " // 'Adjust.Flowrate.Apply not found'
        call json%get('Adjust.Flowrate.approach', c, found)
        call util_lower_case(c)
        if (c == 'vshape') then
            setting%Adjust%Flowrate%approach = vshape
        else
            stop "Error, Adjust.Flowrate.approach not compatible"
        end if
        call json%get('Adjust.Flowrate.Coef', real_value, found)
        setting%Adjust%Flowrate%Coef = real_value
        if (.not. found) stop "Error - setting " // 'Adjust.Flowrate.Coef not found'

        call json%get('Adjust.Head.Apply', logical_value, found)
        setting%Adjust%Head%Apply = logical_value
        if (.not. found) stop "Error - setting " // 'Adjust.Head.Apply not found'
        call json%get('Adjust.Head.approach', c, found)
        call util_lower_case(c)
        if (c == 'vshape_surcharge_only') then
            setting%Adjust%Head%approach = vshape_surcharge_only
        else
            stop "Error, Adjust.Head.approach not compatible"
        end if
        call json%get('Adjust.Head.Coef', real_value, found)
        setting%Adjust%Head%Coef = real_value
        if (.not. found) stop "Error - setting " // 'Adjust.Head.Coef not found'

        ! Load BC Settings
        call json%get('BC.slots', real_value, found)
        setting%BC%slots = real_value
        if (.not. found) stop "Error - setting " // 'BC.slots not found'
        call json%get('BC.disableInterpolation', logical_value, found)
        setting%BC%disableInterpolation = logical_value
        if (.not. found) stop "Error - setting " // 'BC.disableInterpolation not found'

        ! Load Constant Settings
        call json%get('Constant.gravity', real_value, found)
        setting%Constant%gravity = real_value
        if (.not. found) stop "Error - setting " // 'Constant.gravity not found'

        ! For element length adjustment
        call json%get('Discretization.NominalElemLength', real_value, found)
        setting%Discretization%NominalElemLength = real_value
        if (.not. found) stop "Error - setting " // 'Discretization.NominalElemLength not found'
        call json%get('Discretization.LinkShortingFactor', real_value, found)
        setting%Discretization%LinkShortingFactor = real_value
        if (.not. found) stop "Error - setting " // 'Discretization.LinkShortingFactor not found'

        ! Load Eps Settings
        call json%get('Eps.FroudeJump', real_value, found)
        setting%Eps%FroudeJump = real_value
        if (.not. found) stop "Error - setting " // 'Eps.FroudeJump not found'
        call json%get('Eps.InflowDepthIncreaseFroudeLimit', real_value, found)
        setting%Eps%InflowDepthIncreaseFroudeLimit = real_value
        if (.not. found) stop "Error - setting " // 'Eps.InflowDepthIncreaseFroudeLimit not found'

        ! Load FaceInterp Settings
        call json%get('FaceInterp.DownJBFaceInterp', c, found)
        call util_lower_case(c)
        if (c == 'static') then
            setting%FaceInterp%DownJBFaceInterp = static
        else if (c == 'dynamic') then
            setting%FaceInterp%DownJBFaceInterp = dynamic
        else
            print *, "Error, the setting '" // trim(c) // "' is not supported for DownJBFaceInterp"
            stop "in " // subroutine_name
        end if

        ! Load Junction Settings
        call json%get('Junction.kFactor', real_value, found)
        setting%Junction%kFactor = real_value
        if (.not. found) stop "Error - setting " // 'Junction.kFactor not found'

        call json%get('Junction.HeadCoef', real_value, found)
        setting%Junction%HeadCoef = real_value
        if (.not. found) stop "Error - setting " // 'Junction.HeeadCoef not found'

        call json%get('Junction.CFLlimit', real_value, found)
        setting%Junction%CFLlimit = real_value
        if (.not. found) stop "Error - setting " // 'Junction.CFLlimit not found'

        call json%get('Junction.isDynamic', logical_value, found)
        setting%Junction%isDynamic = logical_value
        if (.not. found) stop "Error - setting " // 'Limiter.Junction.isDynamic not found'

        ! Load Limiter Settings
        call json%get('Limiter.BC.approach', c, found)
        call util_lower_case(c)
        if (c == 'froudenumber') then
            setting%Limiter%BC%approach = FroudeNumber
        else
            stop "Error, Limiter.BC.approach not compatible"
        end if
        if (.not. found) stop "Error - setting " // 'Limiter.BC.approach not found'
        call json%get('Limiter.BC.FroudeInflowMaximum', real_value, found)
        setting%Limiter%BC%FroudeInflowMaximum = real_value
        if (.not. found) stop "Error - setting " // 'Limiter.BC.FroudeInflowMaximum not found'
        call json%get('Limiter.BC.UseInflowLimiter', logical_value, found)
        setting%Limiter%BC%UseInflowLimiter = logical_value
        if (.not. found) stop "Error - setting " // 'Limiter.BC.UseInflowLimiter not found'

        call json%get('Limiter.Channel.LargeDepthFactor', real_value, found)
        setting%Limiter%Channel%LargeDepthFactor = real_value
        if (.not. found) stop "Error - setting " // 'Limiter.Channel.LargeDepthFactor not found'

        call json%get('Limiter.Flowrate.FaceVolumeTransport', real_value, found)
        setting%Limiter%Flowrate%FaceVolumeTransport = real_value
        if (.not. found) stop "Error - setting " // 'Limiter.Flowrate.FaceVolumeTransport not found'
        call json%get('Limiter.Flowrate.UseFaceVolumeTransport', logical_value, found)
        setting%Limiter%Flowrate%UseFaceVolumeTransport = logical_value
        if (.not. found) stop "Error - setting " // 'Limiter.Flowrate.UseFaceVolumeTransport not found'

        call json%get('Limiter.InterpWeight.Maximum', real_value, found)
        setting%Limiter%InterpWeight%Maximum = real_value
        if (.not. found) stop "Error - setting " // 'Limiter.InterpWeight.Maximum not found'
        call json%get('Limiter.InterpWeight.Minimum', real_value, found)
        setting%Limiter%InterpWeight%Minimum = real_value
        if (.not. found) stop "Error - setting " // 'Limiter.InterpWeight.Minimum not found'

        call json%get('Limiter.Velocity.Maximum', real_value, found)
        setting%Limiter%Velocity%Maximum = real_value
        if (.not. found) stop "Error - setting " // 'Limiter.Velocity.Maximum not found'
        call json%get('Limiter.Velocity.UseLimitMax', logical_value, found)
        setting%Limiter%Velocity%UseLimitMax = logical_value
        if (.not. found) stop "Error - setting " // 'Limiter.Velocity.UseLimitMax not found'

        call json%get('Limiter.Dt.Minimum', real_value, found)
        setting%Limiter%Dt%Minimum = real_value
        if (.not. found) stop "Error - setting " // 'Limiter.Dt.Minimum not found'
        call json%get('Limiter.Dt.UseLimitMin', logical_value, found)
        setting%Limiter%Dt%UseLimitMin = logical_value
        if (.not. found) stop "Error - setting " // 'Limiter.Dt.UseLimitMin not found'

        ! Load Link settings
        call json%get('Link.DefaultInitDepthType', c, found)
        call util_lower_case(c)
        if (c == 'linear') then
            setting%Link%DefaultInitDepthType = LinearlyVarying
        else if (c == 'uniform') then
            setting%Link%DefaultInitDepthType = Uniform
        else if (c == 'exponential') then
            setting%Link%DefaultInitDepthType = ExponentialDecay
        else
            print *, "Error, the setting '" // trim(c) // "' is not supported for Link.DefaultInitDepthType"
            stop "in " // subroutine_name
        end if
        if (.not. found) stop "Error - setting " // 'Link.DefaultInitDepthType not found'

        call json%get('Link.PropertiesFile', c, found)
        setting%Link%PropertiesFile = c
        if (.not. found) stop "Error - setting " // 'Link.PropertiesFile not found'

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
            print *, "Error, the setting '" // trim(c) // "' is not supported for Partitioning.PartitioningMethod"
            stop 420
        end if
        if (.not. found) stop "Error - setting " // 'Partitioning.PartitioningMethod not found'

        ! Load PreissmanSlot settings
        call json%get('PreissmannSlot.PreissmannSlotMethod', c, found)
        call util_lower_case(c)
        if (c == 'staticslot') then
            setting%PreissmannSlot%PreissmannSlotMethod = StaticSlot
        else if (c == 'variableslot') then
            setting%PreissmannSlot%PreissmannSlotMethod = VariableSlot
        else
            print *, "Error, the setting '" // trim(c) // "' is not supported for PreissmannSlot.PreissmannSlotMethod"
            stop 470
        end if
        if (.not. found) stop "Error - setting " // 'PreissmannSlot.PreissmannSlotMethod not found'
        call json%get('PreissmannSlot.CelerityFactor', real_value, found)
        setting%PreissmannSlot%CelerityFactor = real_value
        if (.not. found) stop "Error - setting " // 'PreissmannSlot.CelerityFactor not found'

        call json%get('Simulation.useHydrology', logical_value, found)
        setting%Simulation%useHydrology = logical_value
        if (.not. found) stop "Error - setting " // 'Simulation.useHydrology not found'
        call json%get('Simulation.useHydraulics', logical_value, found)
        setting%Simulation%useHydraulics = logical_value
        if (.not. found) stop "Error - setting " // 'Simulation.useHydraulics not found'
        call json%get('Simulation.useSWMMC', logical_value, found)
        if (logical_value) then
            setting%Simulation%useSWMMC = 1
        else
            setting%Simulation%useSWMMC = 0
        end if
        if (.not. found) stop "Error - setting " // 'Simulation.useSWMMC not found'

        ! Load SmallVolume Settings
        call json%get('SmallVolume.DepthCutoff', real_value, found)
        setting%SmallVolume%DepthCutoff = real_value
        if (.not. found) stop "Error - setting " // 'SmallVolume.DepthCutoff not found'
        call json%get('SmallVolume.ManningsN', real_value, found)
        setting%SmallVolume%ManningsN = real_value
        if (.not. found) stop "Error - setting " // 'SmallVolume.ManningsN not found'
        call json%get('SmallVolume.MinimumArea', real_value, found)
        setting%SmallVolume%MinimumArea = real_value
        if (.not. found) stop "Error - setting " // 'SmallVolume.MinimumArea not found'
        call json%get('SmallVolume.MinimumHydRadius', real_value, found)
        setting%SmallVolume%MinimumHydRadius = real_value
        if (.not. found) stop "Error - setting " // 'SmallVolume.MinimumHydRadius not found'
        call json%get('SmallVolume.MinimumPerimeter', real_value, found)
        setting%SmallVolume%MinimumPerimeter = real_value
        if (.not. found) stop "Error - setting " // 'SmallVolume.MinimumPerimeter not found'
        call json%get('SmallVolume.MinimumTopwidth', real_value, found)
        setting%SmallVolume%MinimumTopwidth = real_value
        if (.not. found) stop "Error - setting " // 'SmallVolume.MinimumTopwidth not found'
        call json%get('SmallVolume.UseSmallVolumes', logical_value, found)
        setting%SmallVolume%UseSmallVolumes = logical_value
        if (.not. found) stop "Error - setting " // 'SmallVolume.UseSmallVolumes not found'

        ! Load Solver Settings
        ! call json%get('Solver.crk2', real_value, found)
        ! setting%Solver%crk2 = [real_value, real_value]
        ! if (.not. found) stop "Error - setting " // 'Solver.crk2 not found'
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
            stop 540
        end if
        if (.not. found) stop "Error - setting " // 'Solver.MomentumSourceMethod not found'
        call json%get('Solver.PreissmanSlot', logical_value, found)
        setting%Solver%PreissmanSlot = logical_value
        if (.not. found) stop "Error - setting " // 'Solver.PreissmanSlot not found'
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
            stop 570
        end if
        if (.not. found) stop "Error - setting " // 'Solver.SolverSelect not found'
        call json%get('Solver.SwitchFractionDn', real_value, found)
        setting%Solver%SwitchFractionDn = real_value
        if (.not. found) stop "Error - setting " // 'Solver.SwitchFractionDn not found'
        call json%get('Solver.SwitchFractionUp', real_value, found)
        setting%Solver%SwitchFractionUp = real_value
        if (.not. found) stop "Error - setting " // 'Solver.SwitchFractionUp not found'

        ! Load Step Settings
        !rm 20210607 brh call json%get('Step.Current', real_value, found)
        !rm 20210607 brh setting%Step%Current = real_value
        !rm 20210607 brh if (.not. found) stop "Error - setting " // 'Step.Current not found'
        !rm 20210607 brh call json%get('Step.Final', real_value, found)
        !rm 20210607 brh setting%Step%Final = real_value
        !rm 20210607 brh if (.not. found) stop "Error - setting " // 'Step.Final not found'
        !rm 20210607 brh call json%get('Step.First', real_value, found)
        !rm 20210607 brh setting%Step%First = real_value
        !rm 20210607 brh if (.not. found) stop "Error - setting " // 'Step.First not found'

        ! Load Time Settings
        call json%get('Time.matchHydrologyStep', logical_value, found)
        setting%Time%matchHydrologyStep = logical_value
        if (.not. found) stop "Error - setting " // 'Time.matchHydrologyStep not found'
        call json%get('Time.Start', real_value, found)
        setting%Time%Start = real_value
        if (.not. found) stop "Error - setting " // 'Time.Start not found'
        call json%get('Time.Now', real_value, found)
        setting%Time%Now = real_value
        if (.not. found) stop "Error - setting " // 'Time.Now not found'
        call json%get('Time.End', real_value, found)
        setting%Time%End = real_value
        if (.not. found) stop "Error - setting " // 'Time.End not found'
        call json%get('Time.Dt', real_value, found)
        setting%Time%Dt = real_value
        if (.not. found) stop "Error - setting " // 'Time.Dt not found'
        call json%get('Time.Step', integer_value, found)
        setting%Time%Step = integer_value
        if (.not. found) stop "Error - setting " // 'Time.Step not found'
        call json%get('Time.Hydraulics.Dt', real_value, found)
        setting%Time%Hydraulics%Dt = real_value
        if (.not. found) stop "Error - setting " // 'Time.Hydraulics.Dt not found'
        call json%get('Time.Hydraulics.Step', integer_value, found)
        setting%Time%Hydraulics%Step = integer_value
        if (.not. found) stop "Error - setting " // 'Time.Hydraulics.Step not found'
        call json%get('Time.Hydrology.Dt', real_value, found)
        setting%Time%Hydrology%Dt = real_value
        if (.not. found) stop "Error - setting " // 'Time.Hydrology.Dt not found'
        call json%get('Time.Hydrology.Step', integer_value, found)
        setting%Time%Hydrology%Step = integer_value
        if (.not. found) stop "Error - setting " // 'Time.Hydrology.Step not found'
        call json%get('Time.DtTol', real_value, found)
        setting%Time%DtTol = real_value
        if (.not. found) stop "Error - setting " // 'Time.Hydrology.DtTol not found'
        call json%get('Time.DateTimeStamp', c, found)
        setting%Time%DateTimeStamp = c
        if (.not. found) stop "Error - setting " // 'Time.DateTimeStamp not found'

        ! NOTE: these are NOT initialized because we set the times before the json file is read
        !call json%get('Time.Real.EpochStartSeconds', integer_value, found)
        !setting%Time%Real%EpochStartSeconds = integer_value
        !if (.not. found) stop "Error - setting " // 'Time.Real.EpochStartSeconds not found'
        !call json%get('Time.Real.EpochNowSeconds', integer_value, found)
        !setting%Time%Real%EpochNowSeconds = integer_value
        !if (.not. found) stop "Error - setting " // 'Time.Real.EpochNowSeconds not found'
        !call json%get('Time.CPU.EpochStartSeconds', integer_value, found)

        ! NOTE: these are NOT initialized because we set the times before the json file is read
        !setting%Time%CPU%EpochStartSeconds = real_value
        !if (.not. found) stop "Error - setting " // 'Time.CPU.EpochStartSeconds not found'
        !call json%get('Time.CPU.EpochNowSeconds', real_value, found)
        !setting%Time%CPU%EpochNowSeconds = real_value
        !if (.not. found) stop "Error - setting " // 'Time.CPU.EpochNowSeconds not found'
        !call json%get('Time.CPU.EpochFinishSeconds', real_value, found)
        !setting%Time%CPU%EpochFinishSeconds = real_value
        !if (.not. found) stop "Error - setting " // 'Time.CPU.EpochFinishSeconds not found'

        ! Transverse Weir settings
        call json%get('Weir.Transverse.WeirExponent', real_value, found)
        setting%Weir%Transverse%WeirExponent = real_value
        if (.not. found) stop "Error - setting " // 'Weir.Transverse.WeirExponent not found'
        call json%get('Weir.Transverse.WeirContractionFactor', real_value, found)
        setting%Weir%Transverse%WeirContractionFactor = real_value
        if (.not. found) stop "Error - setting " // 'Weir.Transverse.WeirContractionFactor not found'
        call json%get('Weir.Transverse.SideFlowWeirCrestExponent', real_value, found)
        setting%Weir%Transverse%SideFlowWeirCrestExponent = real_value
        if (.not. found) stop "Error - setting " // 'Weir.Transverse.SideFlowWeirCrestExponent not found'
        call json%get('Weir.Transverse.VillemonteCorrectionExponent', real_value, found)
        setting%Weir%Transverse%VillemonteCorrectionExponent = real_value
        if (.not. found) stop "Error - setting " // 'Weir.Transverse.VillemonteCorrectionExponent not found'

        ! Sideflow Weir settings
        call json%get('Weir.SideFlow.WeirExponent', real_value, found)
        setting%Weir%SideFlow%WeirExponent = real_value
        if (.not. found) stop "Error - setting " // 'Weir.SideFlow.WeirExponent not found'
        call json%get('Weir.SideFlow.WeirContractionFactor', real_value, found)
        setting%Weir%SideFlow%WeirContractionFactor = real_value
        if (.not. found) stop "Error - setting " // 'Weir.SideFlow.WeirContractionFactor not found'
        call json%get('Weir.SideFlow.SideFlowWeirCrestExponent', real_value, found)
        setting%Weir%SideFlow%SideFlowWeirCrestExponent = real_value
        if (.not. found) stop "Error - setting " // 'Weir.SideFlow.SideFlowWeirCrestExponent not found'
        call json%get('Weir.SideFlow.VillemonteCorrectionExponent', real_value, found)
        setting%Weir%SideFlow%VillemonteCorrectionExponent = real_value
        if (.not. found) stop "Error - setting " // 'Weir.SideFlow.VillemonteCorrectionExponent not found'

        call json%get('Weir.VNotch.WeirExponent', real_value, found)
        setting%Weir%VNotch%WeirExponent = real_value
        if (.not. found) stop "Error - setting " // 'Weir.VNotch.WeirExponent not found'
        call json%get('Weir.VNotch.WeirContractionFactor', real_value, found)
        setting%Weir%VNotch%WeirContractionFactor = real_value
        if (.not. found) stop "Error - setting " // 'Weir.VNotch.WeirContractionFactor not found'
        call json%get('Weir.VNotch.SideFlowWeirCrestExponent', real_value, found)
        setting%Weir%VNotch%SideFlowWeirCrestExponent = real_value
        if (.not. found) stop "Error - setting " // 'Weir.VNotch.SideFlowWeirCrestExponent not found'
        call json%get('Weir.VNotch.VillemonteCorrectionExponent', real_value, found)
        setting%Weir%VNotch%VillemonteCorrectionExponent = real_value
        if (.not. found) stop "Error - setting " // 'Weir.VNotch.VillemonteCorrectionExponent not found'


        call json%get('Weir.Trapezoidal.WeirExponent', real_value, found)
        setting%Weir%Trapezoidal%WeirExponent = real_value
        if (.not. found) stop "Error - setting " // 'Weir.Trapezoidal.WeirExponent not found'
        call json%get('Weir.Trapezoidal.WeirContractionFactor', real_value, found)
        setting%Weir%Trapezoidal%WeirContractionFactor = real_value
        if (.not. found) stop "Error - setting " // 'Weir.Trapezoidal.WeirContractionFactor not found'
        call json%get('Weir.Trapezoidal.SideFlowWeirCrestExponent', real_value, found)
        setting%Weir%Trapezoidal%SideFlowWeirCrestExponent = real_value
        if (.not. found) stop "Error - setting " // 'Weir.Trapezoidal.SideFlowWeirCrestExponent not found'
        call json%get('Weir.Trapezoidal.VillemonteCorrectionExponent', real_value, found)
        setting%Weir%Trapezoidal%VillemonteCorrectionExponent = real_value
        if (.not. found) stop "Error - setting " // 'Weir.Trapezoidal.VillemonteCorrectionExponent not found'

        !% load variable time step settings
        call json%get('VariableDT.Apply', logical_value, found)
        setting%VariableDT%Apply = logical_value
        if (.not. found) stop "Error - setting " // 'VariableDT.Apply not found'
        call json%get('VariableDT.CFL_hi_max', real_value, found)
        setting%VariableDT%CFL_hi_max = real_value
        if (.not. found) stop "Error - setting " // 'VariableDT.CFL_hi_max not found'
        call json%get('VariableDT.CFL_target', real_value, found)
        setting%VariableDT%CFL_target = real_value
        if (.not. found) stop "Error - setting " // 'VariableDT.CFL_target not found'
        call json%get('VariableDT.CFL_lo_max', real_value, found)
        setting%VariableDT%CFL_lo_max = real_value
        if (.not. found) stop "Error - setting " // 'VariableDT.CFL_lo_max not found'
        call json%get('VariableDT.decreaseFactor', real_value, found)
        setting%VariableDT%decreaseFactor = real_value
        if (.not. found) stop "Error - setting " // 'VariableDT.decreaseFactor not found'
        call json%get('VariableDT.increaseFactor', real_value, found)
        setting%VariableDT%increaseFactor = real_value
        if (.not. found) stop "Error - setting " // 'VariableDT.increaseFactor not found'
        call json%get('VariableDT.NstepsForCheck', integer_value, found)
        setting%VariableDT%NstepsForCheck = integer_value
        if (.not. found) stop "Error - setting " // 'VariableDT.NstepsForCheck not found'

        ! Load ZeroValue Settings
        call json%get('ZeroValue.Area', real_value, found)
        setting%ZeroValue%Area = real_value
        if (.not. found) stop "Error - setting " // 'ZeroValue.Area not found'
        call json%get('ZeroValue.Depth', real_value, found)
        setting%ZeroValue%Depth = real_value
        if (.not. found) stop "Error - setting " // 'ZeroValue.Depth not found'
        call json%get('ZeroValue.Topwidth', real_value, found)
        setting%ZeroValue%Topwidth = real_value
        if (.not. found) stop "Error - setting " // 'ZeroValue.Topwidth not found'
        call json%get('ZeroValue.UseZeroValues', logical_value, found)
        setting%ZeroValue%UseZeroValues = logical_value
        if (.not. found) stop "Error - setting " // 'ZeroValue.UseZeroValues not found'
        call json%get('ZeroValue.Volume', real_value, found)
        setting%ZeroValue%Volume = real_value
        if (.not. found) stop "Error - setting " // 'ZeroValue.Volume not found'

        ! Load Output settings
        call json%get('Output.report', logical_value, found)
        setting%Output%report = logical_value
        if (.not. found) stop "Error - setting " // 'Output.report not found'
        call json%get('Output.reportStartTime', real_value, found)
        setting%Output%reportStartTime = real_value
        if (.not. found) stop "Error - setting " // 'Output.reportStartTime not found'
        call json%get('Output.reportDt', real_value, found)
        setting%Output%reportDt = real_value
        if (.not. found) stop "Error - setting " // 'Output.reportDt not found'
        call json%get('Output.reportStep', integer_value, found)
        setting%Output%reportStep = integer_value
        if (.not. found) stop "Error - setting " // 'Output.reportStep not found'
        call json%get('Output.slots', integer_value, found)
        setting%Output%slots = integer_value
        if (.not. found) stop "Error - setting " // 'Output.slots not found'
        call json%get('Output.links_file', c, found)
        setting%Output%links_file = c
        if (.not. found) stop "Error - setting " // 'Output.links_file not found'
        call json%get('Output.nodes_file', c, found)
        setting%Output%nodes_file = c
        if (.not. found) stop "Error - setting " // 'Output.nodes_file not found'
        call json%get('Output.CommandLine.quiet', logical_value, found)
        setting%Output%CommandLine%quiet = logical_value
        if (.not. found) stop "Error - setting " // 'Output.CommandLine.quiet not found'
        call json%get('Output.CommandLine.interval', integer_value, found)
        setting%Output%CommandLine%interval = integer_value
        if (.not. found) stop "Error - setting " // 'Output.CommandLine.interval not found'

        ! Load verbose or non-verbose run
        call json%get('Verbose', logical_value, found)
        setting%Verbose = logical_value
        if (.not. found) stop "Error - setting " // 'Verbose not found'

        ! Load Debug Settings
        call json%get('Debug.Tests', logical_value, found)
        setting%Debug%Tests = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.Tests not found'
        call json%get('Debug.File.adjust', logical_value, found)
        setting%Debug%File%adjust = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.adjust not found'
        call json%get('Debug.File.boundary_conditions', logical_value, found)
        setting%Debug%File%boundary_conditions = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.boundary_conditions not found'
        call json%get('Debug.File.c_library', logical_value, found)
        setting%Debug%File%c_library = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.c_library not found'
        call json%get('Debug.File.define_globals', logical_value, found)
        setting%Debug%File%define_globals = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.define_globals not found'
        call json%get('Debug.File.define_indexes', logical_value, found)
        setting%Debug%File%define_indexes = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.define_indexes not found'
        call json%get('Debug.File.define_keys', logical_value, found)
        setting%Debug%File%define_keys = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.define_keys not found'
        call json%get('Debug.File.define_settings', logical_value, found)
        setting%Debug%File%define_settings = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.define_settings not found'
        call json%get('Debug.File.define_types', logical_value, found)
        setting%Debug%File%define_types = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.define_types not found'
        call json%get('Debug.File.diagnostic_elements', logical_value, found)
        setting%Debug%File%diagnostic_elements = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.diagnostic_elements not found'
        call json%get('Debug.File.discretization', logical_value, found)
        setting%Debug%File%discretization = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.discretization not found'
        call json%get('Debug.File.face', logical_value, found)
        setting%Debug%File%face = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.face not found'
        call json%get('Debug.File.geometry', logical_value, found)
        setting%Debug%File%geometry = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.geometry not found'
        call json%get('Debug.File.interface', logical_value, found)
        setting%Debug%File%interface = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.interface not found'
        call json%get('Debug.File.initial_condition', logical_value, found)
        setting%Debug%File%initial_condition = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.initial_condition not found'
        call json%get('Debug.File.initialization', logical_value, found)
        setting%Debug%File%initialization = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.initialization not found'
        call json%get('Debug.File.jump', logical_value, found)
        setting%Debug%File%jump = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.jump not found'
        call json%get('Debug.File.lowlevel_rk2', logical_value, found)
        setting%Debug%File%lowlevel_rk2 = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.lowlevel_rk2 not found'
        call json%get('Debug.File.network_define', logical_value, found)
        setting%Debug%File%network_define = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.network_define not found'
        call json%get('Debug.File.orifice_elements', logical_value, found)
        setting%Debug%File%orifice_elements = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.orifice_elements not found'
        call json%get('Debug.File.pack_mask_arrays', logical_value, found)
        setting%Debug%File%pack_mask_arrays = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.pack_mask_arrays not found'
        call json%get('Debug.File.partitioning', logical_value, found)
        setting%Debug%File%partitioning = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.partitioning not found'
        call json%get('Debug.File.pump_elements', logical_value, found)
        setting%Debug%File%pump_elements = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.pump_elements not found'
        call json%get('Debug.File.rectangular_channel', logical_value, found)
        setting%Debug%File%rectangular_channel = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.rectangular_channel not found'
        call json%get('Debug.File.trapezoidal_channel', logical_value, found)
        setting%Debug%File%trapezoidal_channel = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.trapezoidal_channel not found'
        call json%get('Debug.File.runge_kutta2', logical_value, found)
        setting%Debug%File%runge_kutta2 = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.runge_kutta2 not found'
        call json%get('Debug.File.timeloop', logical_value, found)
        setting%Debug%File%timeloop = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.timeloop not found'
        call json%get('Debug.File.update', logical_value, found)
        setting%Debug%File%update = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.update not found'
        call json%get('Debug.File.utility_allocate', logical_value, found)
        setting%Debug%File%utility_allocate = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.utility_allocate not found'
        call json%get('Debug.File.utility_deallocate', logical_value, found)
        setting%Debug%File%utility_deallocate = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.utility_deallocate not found'
        call json%get('Debug.File.utility_array', logical_value, found)
        setting%Debug%File%utility_array = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.utility_array not found'
        call json%get('Debug.File.utility_datetime', logical_value, found)
        setting%Debug%File%utility_datetime = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.utility_datetime not found'
        call json%get('Debug.File.utility_interpolate', logical_value, found)
        setting%Debug%File%utility_interpolate = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.utility_interpolate not found'
        call json%get('Debug.File.utility_output', logical_value, found)
        setting%Debug%File%utility_output = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.utility_output not found'
        call json%get('Debug.File.utility_string', logical_value, found)
        setting%Debug%File%utility_string = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.utility_string not found'
        call json%get('Debug.File.utility', logical_value, found)
        setting%Debug%File%utility = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.utility not found'
        call json%get('Debug.File.weir_elements', logical_value, found)
        setting%Debug%File%weir_elements = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.weir_elements not found'
        call json%get('Debug.File.output', logical_value, found)
        setting%Debug%File%output = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.File.output not found'
        call json%get('Debug.FileGroup.all', logical_value, found)
        setting%Debug%FileGroup%all = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.FileGroup.all not found'
        call json%get('Debug.FileGroup.definitions', logical_value, found)
        setting%Debug%FileGroup%definitions = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.FileGroup.definitions not found'
        call json%get('Debug.FileGroup.finalization', logical_value, found)
        setting%Debug%FileGroup%finalization = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.FileGroup.finalization not found'
        call json%get('Debug.FileGroup.geometry', logical_value, found)
        setting%Debug%FileGroup%geometry = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.FileGroup.geometry not found'
        call json%get('Debug.FileGroup.initialization', logical_value, found)
        setting%Debug%FileGroup%initialization = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.FileGroup.initialization not found'
        call json%get('Debug.FileGroup.interface', logical_value, found)
        setting%Debug%FileGroup%interface = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.FileGroup.interface not found'
        call json%get('Debug.FileGroup.output', logical_value, found)
        setting%Debug%FileGroup%output = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.FileGroup.output not found'
        call json%get('Debug.FileGroup.timeloop', logical_value, found)
        setting%Debug%FileGroup%timeloop = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.FileGroup.timeloop not found'
        call json%get('Debug.FileGroup.utility', logical_value, found)
        setting%Debug%FileGroup%utility = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.FileGroup.utility not found'
        call json%get('Debug.Input', logical_value, found)
        setting%Debug%Input = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.Input not found'
        call json%get('Debug.Output', logical_value, found)
        setting%Debug%Output = logical_value
        if (.not. found) stop "Error - setting " // 'Debug.Output not found'

        call def_update_debug_options()



        !% Load Profile Settings
        call json%get('Profile.YN', logical_value, found)
        setting%Profile%YN = logical_value
        if (.not. found) stop "Error - setting " // 'Profile.YN not found'

        ! call json%get('Profile.Tests', logical_value, found)
        ! setting%Profile%Tests = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.Tests not found'
        ! call json%get('Profile.File.adjust', logical_value, found)
        ! setting%Profile%File%adjust = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.adjust not found'
        ! call json%get('Profile.File.BIPquick', logical_value, found)
        ! setting%Profile%File%BIPquick = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.BIPquick not found'
        ! call json%get('Profile.File.boundary_conditions', logical_value, found)
        ! setting%Profile%File%boundary_conditions = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.boundary_conditions not found'
        ! call json%get('Profile.File.c_library', logical_value, found)
        ! setting%Profile%File%c_library = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.c_library not found'
        ! call json%get('Profile.File.define_globals', logical_value, found)
        ! setting%Profile%File%define_globals = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.define_globals not found'
        ! call json%get('Profile.File.define_indexes', logical_value, found)
        ! setting%Profile%File%define_indexes = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.define_indexes not found'
        ! call json%get('Profile.File.define_keys', logical_value, found)
        ! setting%Profile%File%define_keys = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.define_keys not found'
        ! call json%get('Profile.File.define_settings', logical_value, found)
        ! setting%Profile%File%define_settings = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.define_settings not found'
        ! call json%get('Profile.File.define_types', logical_value, found)
        ! setting%Profile%File%define_types = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.define_types not found'
        ! call json%get('Profile.File.diagnostic_elements', logical_value, found)
        ! setting%Profile%File%diagnostic_elements = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.diagnostic_elements not found'
        ! call json%get('Profile.File.discretization', logical_value, found)
        ! setting%Profile%File%discretization = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.discretization not found'
        ! call json%get('Profile.File.face', logical_value, found)
        ! setting%Profile%File%face = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.face not found'
        ! call json%get('Profile.File.geometry', logical_value, found)
        ! setting%Profile%File%geometry = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.geometry not found'
        ! call json%get('Profile.File.interface', logical_value, found)
        ! setting%Profile%File%interface = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.interface not found'
        ! call json%get('Profile.File.initial_condition', logical_value, found)
        ! setting%Profile%File%initial_condition = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.initial_condition not found'
        ! call json%get('Profile.File.initialization', logical_value, found)
        ! setting%Profile%File%initialization = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.initialization not found'
        ! call json%get('Profile.File.jump', logical_value, found)
        ! setting%Profile%File%jump = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.jump not found'
        ! call json%get('Profile.File.lowlevel_rk2', logical_value, found)
        ! setting%Profile%File%lowlevel_rk2 = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.lowlevel_rk2 not found'
        ! call json%get('Profile.File.network_define', logical_value, found)
        ! setting%Profile%File%network_define = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.network_define not found'
        ! call json%get('Profile.File.orifice_elements', logical_value, found)
        ! setting%Profile%File%orifice_elements = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.orifice_elements not found'
        ! call json%get('Profile.File.pack_mask_arrays', logical_value, found)
        ! setting%Profile%File%pack_mask_arrays = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.pack_mask_arrays not found'
        ! call json%get('Profile.File.partitioning', logical_value, found)
        ! setting%Profile%File%partitioning = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.partitioning not found'
        ! call json%get('Profile.File.pump_elements', logical_value, found)
        ! setting%Profile%File%pump_elements = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.pump_elements not found'
        ! call json%get('Profile.File.rectangular_channel', logical_value, found)
        ! setting%Profile%File%rectangular_channel = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.rectangular_channel not found'
        ! call json%get('Profile.File.trapezoidal_channel', logical_value, found)
        ! setting%Profile%File%trapezoidal_channel = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.trapezoidal_channel not found'
        ! call json%get('Profile.File.runge_kutta2', logical_value, found)
        ! setting%Profile%File%runge_kutta2 = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.runge_kutta2 not found'
        ! call json%get('Profile.File.timeloop', logical_value, found)
        ! setting%Profile%File%timeloop = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.timeloop not found'
        ! call json%get('Profile.File.update', logical_value, found)
        ! setting%Profile%File%update = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.update not found'
        ! call json%get('Profile.File.utility_allocate', logical_value, found)
        ! setting%Profile%File%utility_allocate = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.utility_allocate not found'
        ! call json%get('Profile.File.utility_deallocate', logical_value, found)
        ! setting%Profile%File%utility_deallocate = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.utility_deallocate not found'
        ! call json%get('Profile.File.utility_array', logical_value, found)
        ! setting%Profile%File%utility_array = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.utility_array not found'
        ! call json%get('Profile.File.utility_datetime', logical_value, found)
        ! setting%Profile%File%utility_datetime = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.utility_datetime not found'
        ! call json%get('Profile.File.utility_interpolate', logical_value, found)
        ! setting%Profile%File%utility_interpolate = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.utility_interpolate not found'
        ! call json%get('Profile.File.utility_output', logical_value, found)
        ! setting%Profile%File%utility_output = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.utility_output not found'
        ! call json%get('Profile.File.utility_string', logical_value, found)
        ! setting%Profile%File%utility_string = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.utility_string not found'
        ! call json%get('Profile.File.utility', logical_value, found)
        ! setting%Profile%File%utility = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.utility not found'
        ! call json%get('Profile.File.weir_elements', logical_value, found)
        ! setting%Profile%File%weir_elements = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.weir_elements not found'
        ! call json%get('Profile.File.output', logical_value, found)
        ! setting%Profile%File%output = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.File.output not found'
        ! call json%get('Profile.FileGroup.all', logical_value, found)
        ! setting%Profile%FileGroup%all = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.FileGroup.all not found'
        ! call json%get('Profile.FileGroup.definitions', logical_value, found)
        ! setting%Profile%FileGroup%definitions = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.FileGroup.definitions not found'
        ! call json%get('Profile.FileGroup.finalization', logical_value, found)
        ! setting%Profile%FileGroup%finalization = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.FileGroup.finalization not found'
        ! call json%get('Profile.FileGroup.geometry', logical_value, found)
        ! setting%Profile%FileGroup%geometry = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.FileGroup.geometry not found'
        ! call json%get('Profile.FileGroup.initialization', logical_value, found)
        ! setting%Profile%FileGroup%initialization = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.FileGroup.initialization not found'
        ! call json%get('Profile.FileGroup.interface', logical_value, found)
        ! setting%Profile%FileGroup%interface = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.FileGroup.interface not found'
        ! call json%get('Profile.FileGroup.output', logical_value, found)
        ! setting%Profile%FileGroup%output = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.FileGroup.output not found'
        ! call json%get('Profile.FileGroup.timeloop', logical_value, found)
        ! setting%Profile%FileGroup%timeloop = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.FileGroup.timeloop not found'
        ! call json%get('Profile.FileGroup.utility', logical_value, found)
        ! setting%Profile%FileGroup%utility = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.FileGroup.utility not found'
        ! call json%get('Profile.Input', logical_value, found)
        ! setting%Profile%Input = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.Input not found'
        ! call json%get('Profile.Output', logical_value, found)
        ! setting%Profile%Output = logical_value
        ! if (.not. found) stop "Error - setting " // 'Profile.Output not found'

        ! call def_update_profile_options()

        call json%destroy()
        if (json%failed()) stop "JSON failed to destroy"

        if (setting%Debug%File%define_settings) &
            write(*,"(A,i5,A)") '*** leave ' // subroutine_name // " [Processor ", this_image(), "]"
    end subroutine def_load_settings

    subroutine def_update_debug_options()
        if (setting%Debug%FileGroup%all) then
            setting%Debug%FileGroup%definitions = .true.
            setting%Debug%FileGroup%finalization = .true.
            setting%Debug%FileGroup%geometry = .true.
            setting%Debug%FileGroup%initialization = .true.
            setting%Debug%FileGroup%interface = .true.
            setting%Debug%FileGroup%output  = .true.
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
            setting%Debug%File%trapezoidal_channel = .true.
        end if
        if (setting%Debug%FileGroup%initialization) then
            setting%Debug%File%discretization = .true.
            setting%Debug%File%initial_condition = .true.
            setting%Debug%File%initialization = .true.
            setting%Debug%File%network_define = .true.
            setting%Debug%File%partitioning = .true.
            setting%Debug%File%BIPquick = .true.
            setting%Debug%File%pack_mask_arrays = .true.
        end if
        if (setting%Debug%FileGroup%interface) then
            setting%Debug%File%c_library = .true.
            setting%Debug%File%interface = .true.
        end if
        if (setting%Debug%FileGroup%output) then
            setting%Debug%File%output = .true.
        end if
        if (setting%Debug%FileGroup%timeloop) then
            setting%Debug%File%adjust = .true.
            setting%Debug%File%boundary_conditions = .true.
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
        end if
        if (setting%Debug%FileGroup%utility) then
            setting%Debug%File%utility_allocate = .true.
            setting%Debug%File%utility_deallocate = .true.
            setting%Debug%File%utility_array = .true.
            setting%Debug%File%utility_datetime = .true.
            setting%Debug%File%utility_interpolate = .true.
            setting%Debug%File%utility_output = .true.
            setting%Debug%File%utility_string = .true.
            setting%Debug%File%utility = .true.
        end if
    end subroutine def_update_debug_options

    ! subroutine def_update_profile_options()
    !     if (setting%Profile%FileGroup%all) then
    !         setting%Profile%FileGroup%definitions = .true.
    !         setting%Profile%FileGroup%finalization = .true.
    !         setting%Profile%FileGroup%geometry = .true.
    !         setting%Profile%FileGroup%initialization = .true.
    !         setting%Profile%FileGroup%interface = .true.
    !         setting%Profile%FileGroup%timeloop  = .true.
    !         setting%Profile%FileGroup%utility = .true.
    !     end if
    !     if (setting%Profile%FileGroup%definitions) then
    !         setting%Profile%File%define_globals = .true.
    !         setting%Profile%File%define_indexes = .true.
    !         setting%Profile%File%define_keys = .true.
    !         setting%Profile%File%define_settings = .true.
    !         setting%Profile%File%define_types = .true.
    !     end if
    !     if (setting%Profile%FileGroup%finalization) then
    !         setting%Profile%File%finalization = .true.
    !     end if
    !     if (setting%Profile%FileGroup%geometry) then
    !         setting%Profile%File%geometry = .true.
    !         setting%Profile%File%rectangular_channel = .true.
    !         setting%Profile%File%trapezoidal_channel = .true.
    !     end if
    !     if (setting%Profile%FileGroup%initialization) then
    !         setting%Profile%File%discretization = .true.
    !         setting%Profile%File%initial_condition = .true.
    !         setting%Profile%File%initialization = .true.
    !         setting%Profile%File%network_define = .true.
    !         setting%Profile%File%partitioning = .true.
    !         setting%Profile%File%BIPquick = .true.
    !         setting%Profile%File%pack_mask_arrays = .true.
    !     end if
    !     if (setting%Profile%FileGroup%interface) then
    !         setting%Profile%File%c_library = .true.
    !         setting%Profile%File%interface = .true.
    !     end if
    !     if (setting%Profile%FileGroup%timeloop) then
    !         setting%Profile%File%adjust = .true.
    !         setting%Profile%File%boundary_conditions = .true.
    !         setting%Profile%File%diagnostic_elements = .true.
    !         setting%Profile%File%face = .true.
    !         setting%Profile%File%jump = .true.
    !         setting%Profile%File%lowlevel_rk2 = .true.
    !         setting%Profile%File%orifice_elements = .true.
    !         setting%Profile%File%pump_elements = .true.
    !         setting%Profile%File%runge_kutta2 = .true.
    !         setting%Profile%File%timeloop = .true.
    !         setting%Profile%File%update = .true.
    !         setting%Profile%File%weir_elements = .true.
    !     endif
    !     if (setting%Profile%FileGroup%utility) then
    !         setting%Profile%File%utility_allocate = .true.
    !         setting%Profile%File%utility_deallocate = .true.
    !         setting%Profile%File%utility_array = .true.
    !         setting%Profile%File%utility_datetime = .true.
    !         setting%Profile%File%utility_interpolate = .true.
    !         setting%Profile%File%utility_output = .true.
    !         setting%Profile%File%utility_string = .true.
    !         setting%Profile%File%utility = .true.
    !     end if
    ! end subroutine def_update_profile_options
end module define_settings
