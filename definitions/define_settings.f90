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
!%
!%==========================================================================
!% PUBLIC TYPES
!%==========================================================================
!%
    !% ---------------------------------------------------------------
    !% Third Level Types

    !% setting%ACmethod%Anomaly
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

    !% setting%ACmethod%CFL
    type ACmethodCFLType
        ! maximum cfl for the AC dtau -- may be higher than 1.0)
        real(8) :: CFLmax = 2.0
        ! small maximum cfl when reset to larger dtau
        real(8) :: CFLsmall = 0.05
    end type ACmethodCFLType

    !% setting%ACmethod%Celerity
    type ACmethodCelerityType
        ! celerity ratio of AC wave speed to gravity wave speed (1.0 works)
        real(8) :: RC = 1.0
    end type ACmethodCelerityType

    !% setting%ACmethod%Convergence
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

    !% setting%ACmethod%Iter
    type ACmethodIterType
        ! cutoff for AC iterations without convergence
        integer :: Max = 100
        ! reset by code so that itermin * dtau >  dt
        integer :: Min = 3
        ! allows more iterations in first time step
        integer :: Firststep = 100
    end type ACmethodIterType

    !% setting%ACmethod%Switch
    type ACmethodSwitchType
        ! switch to AC solver if depth/depthMax >  0.9
        real(8) :: Depth = 0.9
        ! switch to AC solver if area/areaMax >  0.9
        real(8) :: Area = 0.9
        ! 5% buffer for the switch
        real(8) :: Buffer = 0.05
    end type ACmethodSwitchType

    !% setting%Adjust%Flowrate
    type AdjustFlowrateType
        logical :: Apply = .true.
        real(8) :: Coef = 1.0
        integer :: Approach = vshape
    end type AdjustFlowrateType

    !% setting%Adjust%Head
    type AdjustHeadType
        logical :: Apply = .false.
        real(8) :: Coef = 1.0
        integer :: Approach = vshape
    end type AdjustHeadType

    !% setting%Adjust%WidthDepth
    type AdjustWidthDepthType
        logical :: Apply = .true.
    end type AdjustWidthDepthType

    !% setting%Output%CommandLine
    type CommandLineType
        logical :: quiet = .false.
        integer :: interval = 10
    end type CommandLineType

    !% setting%Time%CPU
    type CPUTimeType
        real (8) :: EpochStartSeconds = 0.0
        real (8) :: EpochNowSeconds  = 0.0
        real (8) :: EpochFinishSeconds  = 0.0
    end type CPUTimeType

    !% setting%Output%DataLink
    type DataOutType
        logical :: isAreaOut         = .true.
        logical :: isDepthOut        = .true.
        logical :: isFlowrateOut     = .true.
        logical :: isFroudeNumberOut = .false.
        logical :: isHeadOut         = .true.
        logical :: isHydRadiusOut    = .false.
        logical :: isPerimeterOut    = .false.
        logical :: isSlotWidthOut    = .false.
        logical :: isSlotDepthOut    = .false.
        logical :: isTopWidthOut     = .false.
        logical :: isVelocityOut     = .true.
        logical :: isVolumeOut       = .true.
        logical :: isWaveSpeedOut    = .false.
    end type DataOutType

    !% setting%Limiter%BC
    type LimiterBCType
        logical :: UseInflowLimiter = .true.
        integer :: Approach = FroudeNumber
        ! max value of Fr at inflow
        real(8) :: FroudeInflowMaximum = 1.5
    end type LimiterBCType

    !% setting%Limiter%Channel
    type LimiterChannelType
        real(8) :: LargeDepthFactor = 10.0
    end type LimiterChannelType

    !% setting%Limiter%Flowrate
    type LimiterFlowrateType
        logical :: UseFaceVolumeTransport = .true.
        ! Fraction of usptream volume that can be
        ! transported in on time step
        real(8) :: FaceVolumeTransport = 0.5
    end type LimiterFlowrateType

    !% setting%Limiter%InterpWeight
    type LimiterInterpWeightType
        real(8) :: Maximum = 1e6
        real(8) :: Minimum = 1e-6
    end type LimiterInterpWeightType

    !% setting%Limiter%Velocity
    type LimiterVelocityType
        logical :: UseLimitMax = .true.
        real(8) :: Maximum = 10.0 ! m/s
    end type LimiterVelocityType

    !% setting%Limiter%ArraySize
    type LimiterArraySizeType
        integer :: TemporalInflows = 10
        integer :: TotallInflows = 50
    end type LimiterArraySizeType

    !% setting%Limiter%Dt
    type LimiterDtType
        logical :: UseLimitMin = .true.
        real(8) :: Minimum     = 1e-6
    end type LimiterDtType

    !% setting%Time%Real
    type RealTimeType
        integer :: EpochStartSeconds = 0
        integer :: EpochTimeLoopStartSeconds = 0
        integer :: EpochNowSeconds  = 0
    end type RealTimeType

    !% setting%Time% ...Hydraulics, Hydrology
    type TimeStepType
        real(8) :: Dt
        integer :: Step
    end type TimeStepType

    !% setting%File%UnitNumber
    type UnitNumberType
        integer :: inp_file
        integer :: out_file
        integer :: rpt_file
        integer :: setting_file
        integer :: links_input_file
        integer :: nodes_input_file
        integer :: debug_setup_linkR_file
        integer :: debug_setup_linkI_file
        integer :: debug_setup_nodeR_file
        integer :: debug_setup_nodeI_file
        integer :: debug_setup_nodeYN_file
        !integer :: outputML_combined_file
        integer :: outputML_filename_file
        integer :: outputML_control_file

    !     integer :: debug_output_linkR_file
    !     integer :: debug_output_linkI_file
    !     integer :: debug_output_nodeR_file
    !     integer :: debug_output_nodeI_file
    !     integer :: debug_output_nodeYN_file
    !     integer :: debug_output_elemR_file
    !     integer :: debug_output_faceR_file
    !    integer :: swmm5_output_linkR_file
    !    integer :: swmm5_output_linkI_file
    !    integer :: swmm5_output_nodeR_file
    !    integer :: swmm5_output_nodeI_file
    end type

    !% setting%Weir% ...Transverse, SideFlow, VNotch, Trapezoidal
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

    !% ---------------------------------------------------------------
    ! Second Level Types

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

    ! setting%CaseName
    type CaseNameType
        character(256) ::    Long
        character(16)  ::    Short
        character(31)  ::    withTimeStamp
    end type CaseNameType

    ! setting%Constant
    type ConstantType
        real(8) :: gravity = 9.81 ! m^2/s
    end type ConstantType

    !% setting%Debug
    type DebugType
        logical :: Tests = .false.
        type(DebugFileYNType) :: File
        type(DebugFileGroupYNType) :: FileGroup
        logical :: Setup
        logical :: Output
    end type DebugType

    !% setting%Discretization
    type DiscretizationType
        real(8) :: NominalElemLength   = 10.0
        integer :: MinElemLengthMethod = ElemLengthAdjust
        real(8) :: MinElemLengthFactor = 0.50
        real(8) :: LinkShortingFactor  = 0.33
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

    !% setting%File
    type FileType
        !% standard files and folders
        character(len=256)   :: base_folder = ""
        character(len=256)   :: library_folder = ""
        character(len=256)   :: project_folder = "" ! project path
        character(len=256)   :: output_folder= "" !
        character(len=256)   :: output_timestamp_subfolder = ""
        character(len=256)   :: output_temp_subfolder = "temp"
        character(len=256)   :: setting_file = "" ! path to settings JSON file
        character(len=256)   :: input_kernel = "" ! main part of input file name
        character(len=256)   :: output_kernel= "" ! main part ouf output file name
        character(len=256)   :: inp_file = "" ! path to SWMM input (.inp) file
        character(len=256)   :: rpt_file = "" ! path to SWMM report (.rpt) file
        character(len=256)   :: out_file = "" ! path to SWMM output (.out) file
        logical              :: force_folder_creation = .false.
        integer              :: last_unit = 1000 !% starting point for assigning unit numbers
        type(UnitNumberType) :: UnitNumber

        !% for multi-level output
        character(len=256) :: outputML_Link_kernel = "link"
        character(len=256) :: outputML_Node_kernel = "node"
        character(len=256) :: outputML_combinedfile_kernel = "combined_output" !% filenames that combine across images
        character(len=256) :: outputML_filename_file = 'output_filenames.txt'  !% list of interim filenames used in combined output
        character(len=256) :: outputML_control_file = 'control.unf'  !% global parameters that control the processing of combined output
        integer            :: outputML_Ncombined_file_written = 0
        integer            :: outputML_total_timelevels_written = 0

        !% for csv dump output
        character(len=256) :: debug_setup_link_folder = ""
        character(len=256) :: debug_setup_node_folder = ""
        character(len=256) :: debug_output_link_folder = ""
        character(len=256) :: debug_output_node_folder = ""
        character(len=256) :: debug_output_elemR_folder = ""
        character(len=256) :: debug_output_faceR_folder = ""
        character(len=256) :: debug_output_summary_folder = ""
        character(len=256) :: swmm5_output_link_folder = ""
        character(len=256) :: swmm5_output_node_folder = ""
        character(len=256) :: links_input_file = "links_input.csv"
        character(len=256) :: nodes_input_file = "nodes_input.csv"
        character(len=256) :: debug_setup_linkR_file = "linkR.csv"
        character(len=256) :: debug_setup_linkI_file = "linkI.csv"
        character(len=256) :: debug_setup_nodeR_file = "nodeR.csv"
        character(len=256) :: debug_setup_nodeI_file = "nodeI.csv"
        character(len=256) :: debug_setup_nodeYN_file = "nodeYN.csv"
        logical :: links_input_file_exist = .false.
        logical :: nodes_input_file_exist = .false.
    end type FileType

    ! setting%Junction
    type JunctionType
        real(8) :: kFactor      = 0.0   !% default entrance/exit losses at junction branch (use 0.0 as needs debugging)
        real(8) :: HeadCoef     = 1.0   !% junction branch head coef for diagnostic junction (must be > 0)
        real(8) :: CFLlimit     = 0.5   !% limiter on CFL to control dynamic junction
        integer :: FunStorageN  = 10    !% number of curve entries for functional storage
        logical :: isDynamic    = .true.
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

    !% setting%Orifice
    type OrificeType
        real(8) :: SharpCrestedWeirCoefficient
        real(8) :: TransverseWeirExponent
        real(8) :: VillemonteCorrectionExponent
    end type OrificeType

    !% setting%Output
    type OutputType
        logical :: report
        real(8) :: reportStartTime
        real(8) :: reportDt
        integer :: reportStep
        integer :: reportTimeUnits = InHours
        integer :: LastLevel = 0
        integer :: MaxExpectedLevels = 0
        integer :: StoredLevels = 100
        logical :: OutputElementsExist = .false.
        logical :: OutputFacesExist = .false.
        integer :: StoredFileNames = 2
        logical :: UseFileNameFile = .false.
        !integer :: Slots = 20 !% remove?
        integer :: max_links_csv = 100
        integer :: max_nodes_csv = 100
        logical :: print_links_csv = .false.
        logical :: print_nodes_csv = .false.
        logical :: suppress_MultiLevel_Output = .false.
        logical :: Verbose = .true.
        logical :: Warning = .true.
        type(CommandLineType) :: CommandLine
        type(DataOutType) :: DataOut
    end type OutputType

    !% setting%Partitioning
    type PartitioningType
        integer :: PartitioningMethod = BLink
    endtype PartitioningType

    !% setting%PreissmannSlot
    type PreissmannSlotType
        integer :: PreissmannSlotMethod = VariableSlot
        real(8) :: CelerityFactor = 1.0
    end type PreissmannSlotType

    !% setting%Profile
    type ProfileType
        logical :: YN = .false.
        !logical :: Tests = .false.
        !type(ProfileFileYNType) :: File
        !type(ProfileFileGroupYNType) :: FileGroup
        !logical :: Input
        !logical :: Output
    end type ProfileType

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

    !% setting%Weir
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

    !% ---------------------------------------------------------------
    !% First Level Type (setting)

    type settingType
        type(ACmethodType)       :: ACmethod
        type(AdjustType)         :: Adjust
        type(BCPropertiesType)   :: BC
        type(CaseNameType)       :: CaseName  ! name of case
        type(ConstantType)       :: Constant ! Constants
        type(DebugType)          :: Debug
        type(DiscretizationType) :: Discretization
        type(EpsilonType)        :: Eps ! epsilons used to provide bandwidth for comparisons
        type(FaceInterpType)     :: FaceInterp ! Temporary: setting for face interpolation in downstream JB
        type(FileType)           :: File
        type(JunctionType)       :: Junction
        type(LimiterType)        :: Limiter ! maximum and minimum limiters
        type(LinkType)           :: Link
        type(OrificeType)        :: Orifice
        type(OutputType)         :: Output
        type(PartitioningType)   :: Partitioning
        type(PreissmannSlotType) :: PreissmannSlot
        type(ProfileType)        :: Profile
        type(SimulationType)     :: Simulation
        type(SmallVolumeType)    :: SmallVolume ! controls for small volumes
        type(SolverType)         :: Solver ! switch for solver
        type(TestCaseType)       :: TestCase
        type(TimeType)           :: Time ! controls of time step
        type(VariableDTType)     :: VariableDT
        type(WeirType)           :: Weir
        type(ZeroValueType)      :: ZeroValue ! finite values to represent small or negative values
    end type settingType

    type(settingType), target :: setting

contains
!%
!%==========================================================================
!% PUBLIC
!%==========================================================================
!%
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
            write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"

        call json%initialize()
        call json%load(filename = trim(setting%File%setting_file))

    !% --- load the settings from the json file

    !% AC method
        !% --- ACmethod anomaly settings
        call json%get('ACmethod.dtau', real_value, found)
        setting%ACmethod%dtau = real_value
        if (.not. found) then
            write(*,"(A)") "Error - json file - setting " // 'ACmethod.dtau not found'
            write(*,"(A)") "this is first item in json file, which may indicate formatting problem in file"
            stop 970984
        end if

        call json%get('ACmethod.Anomaly.DensityLowCutoff', real_value, found)
        setting%ACmethod%Anomaly%DensityLowCutoff = real_value
        if (.not. found) stop "Error - json file - setting " // 'ACmethod.Anomaly.DensityLowCutoff not found'
        call json%get('ACmethod.Anomaly.FullPipeFactor', real_value, found)
        setting%ACmethod%Anomaly%FullPipeFactor = real_value
        if (.not. found) stop "Error - json file - setting " // 'ACmethod.Anomaly.FullPipeFactor not found'
        call json%get('ACmethod.Anomaly.OpenPipeFactor', real_value, found)
        setting%ACmethod%Anomaly%OpenPipeFactor = real_value
        if (.not. found) stop "Error - json file - setting " // 'ACmethod.Anomaly.OpenPipeFactor not found'
        call json%get('ACmethod.Anomaly.UseDensityCorrection', logical_value, found)
        setting%ACmethod%Anomaly%UseDensityCorrection = logical_value
        if (.not. found) stop "Error - json file - setting " // 'ACmethod.Anomaly.UseDensityCorrection not found'
        call json%get('ACmethod.Anomaly.DensityHighCutoff', real_value, found)
        setting%ACmethod%Anomaly%DensityHighCutoff = real_value
        if (.not. found) stop "Error - json file - setting " // 'ACmethod.Anomaly.DensityHighCutoff not found'

        !% --- ACmethod implicit stencil coefficients
        call json%get('ACmethod.ImplicitCoef.a1', real_value, found)
        setting%ACmethod%ImplicitCoef%a1 = real_value
        if (.not. found) stop "Error - json file - setting " // 'ACmethod.ImplicitCoef.a1 not found'
        call json%get('ACmethod.ImplicitCoef.a2', real_value, found)
        setting%ACmethod%ImplicitCoef%a2 = real_value
        if (.not. found) stop "Error - json file - setting " // 'ACmethod.ImplicitCoef.a2 not found'
        call json%get('ACmethod.ImplicitCoef.a3', real_value, found)
        setting%ACmethod%ImplicitCoef%a3 = real_value
        if (.not. found) stop "Error - json file - setting " // 'ACmethod.ImplicitCoef.a3 not found'

        !% --- AC method CFL settings
        call json%get('ACmethod.CFL.CFLmax', real_value, found)
        setting%ACmethod%CFL%CFLmax = real_value
        if (.not. found) stop "Error - json file - setting " // 'ACmethod.CFL.CFLmax not found'
        call json%get('ACmethod.CFL.CFLmax', real_value, found)
        setting%ACmethod%CFL%CFLmax = real_value
        if (.not. found) stop "Error - json file - setting " // 'ACmethod.CFL.CFLmax not found'

        !% --- AC method celerity settings
        call json%get('ACmethod.Celerity.RC', real_value, found)
        setting%ACmethod%Celerity%RC = real_value
        if (.not. found) stop "Error - json file - setting " // 'ACmethod.Celerity.RC not found'

        !% --- AC method convergence settings
        call json%get('ACmethod.Convergence.Habsolute', real_value, found)
        setting%ACmethod%Convergence%Habsolute = real_value
        if (.not. found) stop "Error - json file - setting " // 'ACmethod.Convergence.Habsolute not found'
        call json%get('ACmethod.Convergence.Hrelative', real_value, found)
        setting%ACmethod%Convergence%Hrelative = real_value
        if (.not. found) stop "Error - json file - setting " // 'ACmethod.Convergence.Hrelative not found'
        call json%get('ACmethod.Convergence.Qabsolute', real_value, found)
        setting%ACmethod%Convergence%Qabsolute = real_value
        if (.not. found) stop "Error - json file - setting " // 'ACmethod.Convergence.Qabsolute not found'
        call json%get('ACmethod.Convergence.Qrelative', real_value, found)
        setting%ACmethod%Convergence%Qrelative = real_value
        if (.not. found) stop "Error - json file - setting " // 'ACmethod.Convergence.Qrelative not found'

        !% --- AC method iter settings
        call json%get('ACmethod.Iter.Firststep', integer_value, found)
        setting%ACmethod%Iter%Firststep = integer_value
        if (.not. found) stop "Error - json file - setting " // 'ACmethod.Iter.Firststep not found'
        call json%get('ACmethod.Iter.Max', integer_value, found)
        setting%ACmethod%Iter%Max = integer_value
        if (.not. found) stop "Error - json file - setting " // 'ACmethod.Iter.Max not found'
        call json%get('ACmethod.Iter.Min', integer_value, found)
        setting%ACmethod%Iter%Min = integer_value
        if (.not. found) stop "Error - json file - setting " // 'ACmethod.Iter.Min not found'

        !% --- AC method switch settings
        call json%get('ACmethod.Switch.Area', real_value, found)
        setting%ACmethod%Switch%Area = real_value
        if (.not. found) stop "Error - json file - setting " // 'ACmethod.Switch.Area not found'
        call json%get('ACmethod.Switch.Buffer', real_value, found)
        setting%ACmethod%Switch%Buffer = real_value
        if (.not. found) stop "Error - json file - setting " // 'ACmethod.Switch.Buffer not found'
        call json%get('ACmethod.Switch.Depth', real_value, found)
        setting%ACmethod%Switch%Depth = real_value
        if (.not. found) stop "Error - json file - setting " // 'ACmethod.Switch.Depth not found'

    !% Adjustments
        !% --- Adjust flowrate settings
        call json%get('Adjust.Flowrate.Apply', logical_value, found)
        setting%Adjust%Flowrate%Apply = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Adjust.Flowrate.Apply not found'
        call json%get('Adjust.Flowrate.approach', c, found)
        call util_lower_case(c)
        if (c == 'vshape') then
            setting%Adjust%Flowrate%approach = vshape
        else
            stop "Error, Adjust.Flowrate.approach not compatible"
        end if
        call json%get('Adjust.Flowrate.Coef', real_value, found)
        setting%Adjust%Flowrate%Coef = real_value
        if (.not. found) stop "Error - json file - setting " // 'Adjust.Flowrate.Coef not found'

        !% --- Adjust head settings
        call json%get('Adjust.Head.Apply', logical_value, found)
        setting%Adjust%Head%Apply = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Adjust.Head.Apply not found'
        call json%get('Adjust.Head.approach', c, found)
        call util_lower_case(c)
        if (c == 'vshape_surcharge_only') then
            setting%Adjust%Head%approach = vshape_surcharge_only
        else
            stop "Error, Adjust.Head.approach not compatible"
        end if
        call json%get('Adjust.Head.Coef', real_value, found)
        setting%Adjust%Head%Coef = real_value
        if (.not. found) stop "Error - json file - setting " // 'Adjust.Head.Coef not found'

    !% BC
        !% --- BC settings
        call json%get('BC.slots', real_value, found)
        setting%BC%slots = real_value
        if (.not. found) stop "Error - json file - setting " // 'BC.slots not found'
        call json%get('BC.disableInterpolation', logical_value, found)
        setting%BC%disableInterpolation = logical_value
        if (.not. found) stop "Error - json file - setting " // 'BC.disableInterpolation not found'

    !% Case
        !%--- Case name settings
        call json%get('CaseName.Long', c, found)
        if (.not. found) stop "Error - json file - setting " // 'CaseName.Long not found'
        setting%CaseName%Long = trim(c)
        call json%get('CaseName.Short', c, found)
        if (.not. found) stop "Error - json file - setting " // 'CaseName.Short not found'
        setting%CaseName%Short = trim(c)

    !% Constant
        !% --- Constant settings
        call json%get('Constant.gravity', real_value, found)
        setting%Constant%gravity = real_value
        if (.not. found) stop "Error - json file - setting " // 'Constant.gravity not found'

    !% Discretization
        !% -- Nominal element length adjustment
        call json%get('Discretization.NominalElemLength', real_value, found)
        setting%Discretization%NominalElemLength = real_value
        !% -- Minimum CC element length
        call json%get('Discretization.MinElemLengthFactor', real_value, found)
        setting%Discretization%MinElemLengthFactor = real_value
        if (.not. found) stop "Error - json file - setting " // 'Discretization.MinElemLengthFactor not found'
        !% -- Minimum CC element length algorithm
        call json%get('Discretization.MinElemLengthMethod', c, found)
        call util_lower_case(c)
        if (c == 'elemlengthadjust') then
            setting%Discretization%MinElemLengthMethod = ElemLengthAdjust
        elseif (c == 'rawelemlength') then
            setting%Discretization%MinElemLengthMethod = RawElemLength
        else
            stop "Error, Discretization.MinElemLengthMethod not compatible. See data_keys.f90"
        end if

        call json%get('Discretization.LinkShortingFactor', real_value, found)
        setting%Discretization%LinkShortingFactor = real_value
        if (.not. found) stop "Error - json file - setting " // 'Discretization.LinkShortingFactor not found'

    !% epsilon
        !% --- Eps Settings
        call json%get('Eps.FroudeJump', real_value, found)
        setting%Eps%FroudeJump = real_value
        if (.not. found) stop "Error - json file - setting " // 'Eps.FroudeJump not found'
        call json%get('Eps.InflowDepthIncreaseFroudeLimit', real_value, found)
        setting%Eps%InflowDepthIncreaseFroudeLimit = real_value
        if (.not. found) stop "Error - json file - setting " // 'Eps.InflowDepthIncreaseFroudeLimit not found'

    !% Face interpolation
        !% --- FaceInterp settings
        call json%get('FaceInterp.DownJBFaceInterp', c, found)
        call util_lower_case(c)
        if (c == 'static') then
            setting%FaceInterp%DownJBFaceInterp = static
        else if (c == 'dynamic') then
            setting%FaceInterp%DownJBFaceInterp = dynamic
        else
            print *, "Error, the setting '" // trim(c) // "' is not supported for DownJBFaceInterp"
            stop
        end if

    !% Files
        !% --- (filenames and paths should not be read)
        !% -- links and nodes input files exist
        call json%get('File.links_input_file_exist', logical_value, found)
        setting%File%links_input_file_exist = logical_value
        if (.not. found) stop "Error - json file - setting " // 'File.links_input_file_exist not found'
        call json%get('File.nodes_input_file_exist', logical_value, found)
        setting%File%nodes_input_file_exist = logical_value
        if (.not. found) stop "Error - json file - setting " // 'File.nodes_input_file_exist not found'

        call json%get('File.force_folder_creation', logical_value, found)
        setting%File%force_folder_creation = logical_value
        if (.not. found) stop "Error - json file - setting " // 'File.force_folder_creation not found'

        call json%get('File.library_folder', c, found)
        setting%File%library_folder = trim(c)
        if (.not. found) stop "Error - json file - setting "// 'File.library_folder not found'

    !% Junctions
        !%--- Junction Settings
        call json%get('Junction.kFactor', real_value, found)
        setting%Junction%kFactor = real_value
        if (.not. found) stop "Error - json file - setting " // 'Junction.kFactor not found'

        call json%get('Junction.HeadCoef', real_value, found)
        setting%Junction%HeadCoef = real_value
        if (.not. found) stop "Error - json file - setting " // 'Junction.HeeadCoef not found'

        call json%get('Junction.CFLlimit', real_value, found)
        setting%Junction%CFLlimit = real_value
        if (.not. found) stop "Error - json file - setting " // 'Junction.CFLlimit not found'

        call json%get('Junction.FunStorageN', integer_value, found)
        setting%Junction%FunStorageN = integer_value
        if (.not. found) stop "Error - json file - setting " // 'Junction.CFLlimit not found'

        call json%get('Junction.isDynamic', logical_value, found)
        setting%Junction%isDynamic = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Limiter.Junction.isDynamic not found'

    !% Limiters
        !% --- Limiter settings
        call json%get('Limiter.BC.approach', c, found)
        call util_lower_case(c)
        if (c == 'froudenumber') then
            setting%Limiter%BC%approach = FroudeNumber
        else
            stop "Error, Limiter.BC.approach not compatible"
        end if
        if (.not. found) stop "Error - json file - setting " // 'Limiter.BC.approach not found'
        call json%get('Limiter.BC.FroudeInflowMaximum', real_value, found)
        setting%Limiter%BC%FroudeInflowMaximum = real_value
        if (.not. found) stop "Error - json file - setting " // 'Limiter.BC.FroudeInflowMaximum not found'
        call json%get('Limiter.BC.UseInflowLimiter', logical_value, found)
        setting%Limiter%BC%UseInflowLimiter = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Limiter.BC.UseInflowLimiter not found'

        call json%get('Limiter.Channel.LargeDepthFactor', real_value, found)
        setting%Limiter%Channel%LargeDepthFactor = real_value
        if (.not. found) stop "Error - json file - setting " // 'Limiter.Channel.LargeDepthFactor not found'

        call json%get('Limiter.Flowrate.FaceVolumeTransport', real_value, found)
        setting%Limiter%Flowrate%FaceVolumeTransport = real_value
        if (.not. found) stop "Error - json file - setting " // 'Limiter.Flowrate.FaceVolumeTransport not found'
        call json%get('Limiter.Flowrate.UseFaceVolumeTransport', logical_value, found)
        setting%Limiter%Flowrate%UseFaceVolumeTransport = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Limiter.Flowrate.UseFaceVolumeTransport not found'

        call json%get('Limiter.InterpWeight.Maximum', real_value, found)
        setting%Limiter%InterpWeight%Maximum = real_value
        if (.not. found) stop "Error - json file - setting " // 'Limiter.InterpWeight.Maximum not found'
        call json%get('Limiter.InterpWeight.Minimum', real_value, found)
        setting%Limiter%InterpWeight%Minimum = real_value
        if (.not. found) stop "Error - json file - setting " // 'Limiter.InterpWeight.Minimum not found'

        call json%get('Limiter.Velocity.Maximum', real_value, found)
        setting%Limiter%Velocity%Maximum = real_value
        if (.not. found) stop "Error - json file - setting " // 'Limiter.Velocity.Maximum not found'
        call json%get('Limiter.Velocity.UseLimitMax', logical_value, found)
        setting%Limiter%Velocity%UseLimitMax = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Limiter.Velocity.UseLimitMax not found'

        call json%get('Limiter.Dt.Minimum', real_value, found)
        setting%Limiter%Dt%Minimum = real_value
        if (.not. found) stop "Error - json file - setting " // 'Limiter.Dt.Minimum not found'
        call json%get('Limiter.Dt.UseLimitMin', logical_value, found)
        setting%Limiter%Dt%UseLimitMin = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Limiter.Dt.UseLimitMin not found'

    !% Links
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
            stop 
        end if
        if (.not. found) stop "Error - json file - setting " // 'Link.DefaultInitDepthType not found'

        call json%get('Link.PropertiesFile', c, found)
        setting%Link%PropertiesFile = c
        if (.not. found) stop "Error - json file - setting " // 'Link.PropertiesFile not found'

    !% Output
        !% --- Report settings
        call json%get('Output.report', logical_value, found)
        setting%Output%report = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Output.report not found'

        call json%get('Output.reportStartTime', real_value, found)
        setting%Output%reportStartTime = real_value
        if (.not. found) stop "Error - json file - setting " // 'Output.reportStartTime not found'

        call json%get('Output.reportDt', real_value, found)
        setting%Output%reportDt = real_value
        if (.not. found) stop "Error - json file - setting " // 'Output.reportDt not found'

        call json%get('Output.reportStep', integer_value, found)
        setting%Output%reportStep = integer_value
        if (.not. found) stop "Error - json file - setting " // 'Output.reportStep not found'

       ! Load BIPQuick settings
        call json%get('Output.reportTimeUnits', c, found)
        call util_lower_case(c)
        if (c == 'seconds') then
            setting%Output%reportTimeUnits = InSeconds
        else if (c == 'minutes') then
            setting%Output%reportTimeUnits = InMinutes
        else if (c == 'hours') then
            setting%Output%reportTimeUnits = InHours
        else if (c == 'days') then
            setting%Output%reportTimeUnits = InDays
        else
            print *, "Error, the setting '" // trim(c) // "' is not supported for Output.reportTimeUnits"
            stop 4201
        end if
        if (.not. found) stop "Error - json file - setting " // 'Output.reportTimeUnits'


        call json%get('Output.StoredLevels', integer_value, found)
        setting%Output%StoredLevels = integer_value
        if (.not. found) stop "Error - json file - setting " // 'Output.StoredLevels not found'

        call json%get('Output.StoredFileNames', integer_value, found)
        setting%Output%StoredFileNames = integer_value
        if (.not. found) stop "Error - json file - setting " // 'Output.StoredFileNames not found'

        call json%get('Output.UseFileNameFile', logical_value, found)
        setting%Output%UseFileNameFile = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Output.UseFileNameFile not found'

        !% --- HACK NOT USED?
        !call json%get('Output.slots', integer_value, found)
        !setting%Output%slots = integer_value
        !if (.not. found) stop "Error - json file - setting " // 'Output.slots not found'

        !% --- suppression of output
        call json%get('Output.suppress_MultiLevel_Output', logical_value, found)
        setting%Output%suppress_MultiLevel_Output = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Output.suppress_MultiLevel_Output not found'


        !% --- link and node text files
        call json%get('Output.print_links_csv', logical_value, found)
        setting%Output%print_links_csv = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Output.print_links_csv not found'

        call json%get('Output.print_nodes_csv', logical_value, found)
        setting%Output%print_nodes_csv = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Output.print_nodes_csv not found'

        !% --- command line
        call json%get('Output.CommandLine.quiet', logical_value, found)
        setting%Output%CommandLine%quiet = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Output.CommandLine.quiet not found'
        call json%get('Output.CommandLine.interval', integer_value, found)
        setting%Output%CommandLine%interval = integer_value
        if (.not. found) stop "Error - json file - setting " // 'Output.CommandLine.interval not found'

        !% -- verbose or non-verbose run
        call json%get('Output.Verbose', logical_value, found)
        setting%Output%Verbose = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Output.Verbose not found'

        !% --- warnings
        call json%get('Output.Warning', logical_value, found)
        setting%Output%Warning = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Output.Warning not found'

        !% -------------------------
        !% --- data out ----
        !%
        !% --- Area
        call json%get('Output.DataOut.isAreaOut', logical_value, found)
        setting%Output%DataOut%isAreaOut = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Output.DataOut.isAreaOut not found'
        !% --- Depth
        call json%get('Output.DataOut.isDepthOut', logical_value, found)
        setting%Output%DataOut%isDepthOut = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Output.DataOut.isDepthOut not found'
        !% --- Flowrate
        call json%get('Output.DataOut.isFlowrateOut', logical_value, found)
        setting%Output%DataOut%isFlowrateOut = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Output.DataOut.isFlowrateOut not found'
        !% --- Froude Number
        call json%get('Output.DataOut.isFroudeNumberOut', logical_value, found)
        setting%Output%DataOut%isFroudeNumberOut = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Output.DataOut.isFroudeNumberOut not found'
        !% --- Head
        call json%get('Output.DataOut.isHeadOut', logical_value, found)
        setting%Output%DataOut%isHeadOut = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Output.DataOut.isHeadOut not found'
        !% --- Hydraulic Radius
        call json%get('Output.DataOut.isHydRadiusOut', logical_value, found)
        setting%Output%DataOut%isHydRadiusOut = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Output.DataOut.isHydRadiusOut not found'
        !% --- Perimeter
        call json%get('Output.DataOut.isPerimeterOut', logical_value, found)
        setting%Output%DataOut%isPerimeterOut = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Output.DataOut.isPerimeterOut not found'
        !% --- Slot width
        call json%get('Output.DataOut.isSlotWidthOut', logical_value, found)
        setting%Output%DataOut%isSlotWidthOut = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Output.DataOut.isSlotWidthOut not found'
        !% --- Slot depth
        call json%get('Output.DataOut.isSlotDepthOut', logical_value, found)
        setting%Output%DataOut%isSlotDepthOut = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Output.DataOut.isSlotDepthOut not found'
        !% --- Top Width
        call json%get('Output.DataOut.isTopWidthOut', logical_value, found)
        setting%Output%DataOut%isTopWidthOut = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Output.DataOut.isTopWidthOut not found'
        !% --- Velocity
        call json%get('Output.DataOut.isVelocityOut', logical_value, found)
        setting%Output%DataOut%isVelocityOut = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Output.DataOut.isVelocityOut not found'
        !% --- Volume
        call json%get('Output.DataOut.isVolumeOut', logical_value, found)
        setting%Output%DataOut%isVolumeOut = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Output.DataOut.isVolumeOut not found'

        !-------------------------
        ! --- Orifice settings
        call json%get('Orifice.SharpCrestedWeirCoefficient', real_value, found)
        setting%Orifice%SharpCrestedWeirCoefficient = real_value
        if (.not. found) stop "Error - setting " // 'Orifice.SharpCrestedWeirCoefficient not found'
        call json%get('Orifice.TransverseWeirExponent', real_value, found)
        setting%Orifice%TransverseWeirExponent = real_value
        if (.not. found) stop "Error - setting " // 'Orifice.TransverseWeirExponent not found'
        call json%get('Orifice.VillemonteCorrectionExponent', real_value, found)
        setting%Orifice%VillemonteCorrectionExponent = real_value
        if (.not. found) stop "Error - setting " // 'Orifice.VillemonteCorrectionExponent not found'

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
        if (.not. found) stop "Error - json file - setting " // 'Partitioning.PartitioningMethod not found'

    !% Preissman slot
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
        if (.not. found) stop "Error - json file - setting " // 'PreissmannSlot.PreissmannSlotMethod not found'
        call json%get('PreissmannSlot.CelerityFactor', real_value, found)
        setting%PreissmannSlot%CelerityFactor = real_value
        if (.not. found) stop "Error - json file - setting " // 'PreissmannSlot.CelerityFactor not found'

    !% Simulation controls
        !% --- Hydrology and hydraulics simulation controls
        call json%get('Simulation.useHydrology', logical_value, found)
        setting%Simulation%useHydrology = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Simulation.useHydrology not found'
        call json%get('Simulation.useHydraulics', logical_value, found)
        setting%Simulation%useHydraulics = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Simulation.useHydraulics not found'
        call json%get('Simulation.useSWMMC', logical_value, found)
        if (logical_value) then
            setting%Simulation%useSWMMC = 1
        else
            setting%Simulation%useSWMMC = 0
        end if
        if (.not. found) stop "Error - json file - setting " // 'Simulation.useSWMMC not found'

    !% Small volume
        !% --- SmallVolume settings
        call json%get('SmallVolume.DepthCutoff', real_value, found)
        setting%SmallVolume%DepthCutoff = real_value
        if (.not. found) stop "Error - json file - setting " // 'SmallVolume.DepthCutoff not found'
        call json%get('SmallVolume.ManningsN', real_value, found)
        setting%SmallVolume%ManningsN = real_value
        if (.not. found) stop "Error - json file - setting " // 'SmallVolume.ManningsN not found'
        call json%get('SmallVolume.MinimumArea', real_value, found)
        setting%SmallVolume%MinimumArea = real_value
        if (.not. found) stop "Error - json file - setting " // 'SmallVolume.MinimumArea not found'
        call json%get('SmallVolume.MinimumHydRadius', real_value, found)
        setting%SmallVolume%MinimumHydRadius = real_value
        if (.not. found) stop "Error - json file - setting " // 'SmallVolume.MinimumHydRadius not found'
        call json%get('SmallVolume.MinimumPerimeter', real_value, found)
        setting%SmallVolume%MinimumPerimeter = real_value
        if (.not. found) stop "Error - json file - setting " // 'SmallVolume.MinimumPerimeter not found'
        call json%get('SmallVolume.MinimumTopwidth', real_value, found)
        setting%SmallVolume%MinimumTopwidth = real_value
        if (.not. found) stop "Error - json file - setting " // 'SmallVolume.MinimumTopwidth not found'
        call json%get('SmallVolume.UseSmallVolumes', logical_value, found)
        setting%SmallVolume%UseSmallVolumes = logical_value
        if (.not. found) stop "Error - json file - setting " // 'SmallVolume.UseSmallVolumes not found'

    !% Solver
        !% ---Solver Settings
        ! call json%get('Solver.crk2', real_value, found)
        ! setting%Solver%crk2 = [real_value, real_value]
        ! if (.not. found) stop "Error - json file - setting " // 'Solver.crk2 not found'

        !% --- Momentum source
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
        if (.not. found) stop "Error - json file - setting " // 'Solver.MomentumSourceMethod not found'

        !% --- Preissmann slot
        call json%get('Solver.PreissmanSlot', logical_value, found)
        setting%Solver%PreissmanSlot = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Solver.PreissmanSlot not found'

        !% --- solver selection
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
        if (.not. found) stop "Error - json file - setting " // 'Solver.SolverSelect not found'

        !% --- solver switch fraction
        call json%get('Solver.SwitchFractionDn', real_value, found)
        setting%Solver%SwitchFractionDn = real_value
        if (.not. found) stop "Error - json file - setting " // 'Solver.SwitchFractionDn not found'
        call json%get('Solver.SwitchFractionUp', real_value, found)
        setting%Solver%SwitchFractionUp = real_value
        if (.not. found) stop "Error - json file - setting " // 'Solver.SwitchFractionUp not found'

    !% Time
        !% --- Time settings
        call json%get('Time.matchHydrologyStep', logical_value, found)
        setting%Time%matchHydrologyStep = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Time.matchHydrologyStep not found'
        call json%get('Time.Start', real_value, found)
        setting%Time%Start = real_value
        if (.not. found) stop "Error - json file - setting " // 'Time.Start not found'
        call json%get('Time.Now', real_value, found)
        setting%Time%Now = real_value
        if (.not. found) stop "Error - json file - setting " // 'Time.Now not found'
        call json%get('Time.End', real_value, found)
        setting%Time%End = real_value
        if (.not. found) stop "Error - json file - setting " // 'Time.End not found'
        call json%get('Time.Dt', real_value, found)
        setting%Time%Dt = real_value
        if (.not. found) stop "Error - json file - setting " // 'Time.Dt not found'
        call json%get('Time.Step', integer_value, found)
        setting%Time%Step = integer_value
        if (.not. found) stop "Error - json file - setting " // 'Time.Step not found'
        call json%get('Time.Hydraulics.Dt', real_value, found)
        setting%Time%Hydraulics%Dt = real_value
        if (.not. found) stop "Error - json file - setting " // 'Time.Hydraulics.Dt not found'
        call json%get('Time.Hydraulics.Step', integer_value, found)
        setting%Time%Hydraulics%Step = integer_value
        if (.not. found) stop "Error - json file - setting " // 'Time.Hydraulics.Step not found'
        call json%get('Time.Hydrology.Dt', real_value, found)
        setting%Time%Hydrology%Dt = real_value
        if (.not. found) stop "Error - json file - setting " // 'Time.Hydrology.Dt not found'
        call json%get('Time.Hydrology.Step', integer_value, found)
        setting%Time%Hydrology%Step = integer_value
        if (.not. found) stop "Error - json file - setting " // 'Time.Hydrology.Step not found'
        call json%get('Time.DtTol', real_value, found)
        setting%Time%DtTol = real_value
        if (.not. found) stop "Error - json file - setting " // 'Time.Hydrology.DtTol not found'
        call json%get('Time.DateTimeStamp', c, found)
        setting%Time%DateTimeStamp = c
        if (.not. found) stop "Error - json file - setting " // 'Time.DateTimeStamp not found'

        ! NOTE: these are NOT initialized because we set the times before the json file is read
        !call json%get('Time.Real.EpochStartSeconds', integer_value, found)
        !setting%Time%Real%EpochStartSeconds = integer_value
        !if (.not. found) stop "Error - json file - setting " // 'Time.Real.EpochStartSeconds not found'
        !call json%get('Time.Real.EpochNowSeconds', integer_value, found)
        !setting%Time%Real%EpochNowSeconds = integer_value
        !if (.not. found) stop "Error - json file - setting " // 'Time.Real.EpochNowSeconds not found'
        !call json%get('Time.CPU.EpochStartSeconds', integer_value, found)

        ! NOTE: these are NOT initialized because we set the times before the json file is read
        !setting%Time%CPU%EpochStartSeconds = real_value
        !if (.not. found) stop "Error - json file - setting " // 'Time.CPU.EpochStartSeconds not found'
        !call json%get('Time.CPU.EpochNowSeconds', real_value, found)
        !setting%Time%CPU%EpochNowSeconds = real_value
        !if (.not. found) stop "Error - json file - setting " // 'Time.CPU.EpochNowSeconds not found'
        !call json%get('Time.CPU.EpochFinishSeconds', real_value, found)
        !setting%Time%CPU%EpochFinishSeconds = real_value
        !if (.not. found) stop "Error - json file - setting " // 'Time.CPU.EpochFinishSeconds not found'

    !% Weir
        !% --- Transverse weir settings
        call json%get('Weir.Transverse.WeirExponent', real_value, found)
        setting%Weir%Transverse%WeirExponent = real_value
        if (.not. found) stop "Error - json file - setting " // 'Weir.Transverse.WeirExponent not found'
        call json%get('Weir.Transverse.WeirContractionFactor', real_value, found)
        setting%Weir%Transverse%WeirContractionFactor = real_value
        if (.not. found) stop "Error - json file - setting " // 'Weir.Transverse.WeirContractionFactor not found'
        call json%get('Weir.Transverse.SideFlowWeirCrestExponent', real_value, found)
        setting%Weir%Transverse%SideFlowWeirCrestExponent = real_value
        if (.not. found) stop "Error - json file - setting " // 'Weir.Transverse.SideFlowWeirCrestExponent not found'
        call json%get('Weir.Transverse.VillemonteCorrectionExponent', real_value, found)
        setting%Weir%Transverse%VillemonteCorrectionExponent = real_value
        if (.not. found) stop "Error - json file - setting " // 'Weir.Transverse.VillemonteCorrectionExponent not found'

        !% --- Sideflow weir settings
        call json%get('Weir.SideFlow.WeirExponent', real_value, found)
        setting%Weir%SideFlow%WeirExponent = real_value
        if (.not. found) stop "Error - json file - setting " // 'Weir.SideFlow.WeirExponent not found'
        call json%get('Weir.SideFlow.WeirContractionFactor', real_value, found)
        setting%Weir%SideFlow%WeirContractionFactor = real_value
        if (.not. found) stop "Error - json file - setting " // 'Weir.SideFlow.WeirContractionFactor not found'
        call json%get('Weir.SideFlow.SideFlowWeirCrestExponent', real_value, found)
        setting%Weir%SideFlow%SideFlowWeirCrestExponent = real_value
        if (.not. found) stop "Error - json file - setting " // 'Weir.SideFlow.SideFlowWeirCrestExponent not found'
        call json%get('Weir.SideFlow.VillemonteCorrectionExponent', real_value, found)
        setting%Weir%SideFlow%VillemonteCorrectionExponent = real_value
        if (.not. found) stop "Error - json file - setting " // 'Weir.SideFlow.VillemonteCorrectionExponent not found'

        !% --- V-notch weir settings
        call json%get('Weir.VNotch.WeirExponent', real_value, found)
        setting%Weir%VNotch%WeirExponent = real_value
        if (.not. found) stop "Error - json file - setting " // 'Weir.VNotch.WeirExponent not found'
        call json%get('Weir.VNotch.WeirContractionFactor', real_value, found)
        setting%Weir%VNotch%WeirContractionFactor = real_value
        if (.not. found) stop "Error - json file - setting " // 'Weir.VNotch.WeirContractionFactor not found'
        call json%get('Weir.VNotch.SideFlowWeirCrestExponent', real_value, found)
        setting%Weir%VNotch%SideFlowWeirCrestExponent = real_value
        if (.not. found) stop "Error - json file - setting " // 'Weir.VNotch.SideFlowWeirCrestExponent not found'
        call json%get('Weir.VNotch.VillemonteCorrectionExponent', real_value, found)
        setting%Weir%VNotch%VillemonteCorrectionExponent = real_value
        if (.not. found) stop "Error - json file - setting " // 'Weir.VNotch.VillemonteCorrectionExponent not found'

        !% --- Trapezoidal weir settings
        call json%get('Weir.Trapezoidal.WeirExponent', real_value, found)
        setting%Weir%Trapezoidal%WeirExponent = real_value
        if (.not. found) stop "Error - json file - setting " // 'Weir.Trapezoidal.WeirExponent not found'
        call json%get('Weir.Trapezoidal.WeirContractionFactor', real_value, found)
        setting%Weir%Trapezoidal%WeirContractionFactor = real_value
        if (.not. found) stop "Error - json file - setting " // 'Weir.Trapezoidal.WeirContractionFactor not found'
        call json%get('Weir.Trapezoidal.SideFlowWeirCrestExponent', real_value, found)
        setting%Weir%Trapezoidal%SideFlowWeirCrestExponent = real_value
        if (.not. found) stop "Error - json file - setting " // 'Weir.Trapezoidal.SideFlowWeirCrestExponent not found'
        call json%get('Weir.Trapezoidal.VillemonteCorrectionExponent', real_value, found)
        setting%Weir%Trapezoidal%VillemonteCorrectionExponent = real_value
        if (.not. found) stop "Error - json file - setting " // 'Weir.Trapezoidal.VillemonteCorrectionExponent not found'

    !% Variable time step
        !% ---variable time step settings
        call json%get('VariableDT.Apply', logical_value, found)
        setting%VariableDT%Apply = logical_value
        if (.not. found) stop "Error - json file - setting " // 'VariableDT.Apply not found'
        call json%get('VariableDT.CFL_hi_max', real_value, found)
        setting%VariableDT%CFL_hi_max = real_value
        if (.not. found) stop "Error - json file - setting " // 'VariableDT.CFL_hi_max not found'
        call json%get('VariableDT.CFL_target', real_value, found)
        setting%VariableDT%CFL_target = real_value
        if (.not. found) stop "Error - json file - setting " // 'VariableDT.CFL_target not found'
        call json%get('VariableDT.CFL_lo_max', real_value, found)
        setting%VariableDT%CFL_lo_max = real_value
        if (.not. found) stop "Error - json file - setting " // 'VariableDT.CFL_lo_max not found'
        call json%get('VariableDT.decreaseFactor', real_value, found)
        setting%VariableDT%decreaseFactor = real_value
        if (.not. found) stop "Error - json file - setting " // 'VariableDT.decreaseFactor not found'
        call json%get('VariableDT.increaseFactor', real_value, found)
        setting%VariableDT%increaseFactor = real_value
        if (.not. found) stop "Error - json file - setting " // 'VariableDT.increaseFactor not found'
        call json%get('VariableDT.NstepsForCheck', integer_value, found)
        setting%VariableDT%NstepsForCheck = integer_value
        if (.not. found) stop "Error - json file - setting " // 'VariableDT.NstepsForCheck not found'

    !% Zero Value
        !% --- ZeroValue settings
        call json%get('ZeroValue.Area', real_value, found)
        setting%ZeroValue%Area = real_value
        if (.not. found) stop "Error - json file - setting " // 'ZeroValue.Area not found'
        call json%get('ZeroValue.Depth', real_value, found)
        setting%ZeroValue%Depth = real_value
        if (.not. found) stop "Error - json file - setting " // 'ZeroValue.Depth not found'
        call json%get('ZeroValue.Topwidth', real_value, found)
        setting%ZeroValue%Topwidth = real_value
        if (.not. found) stop "Error - json file - setting " // 'ZeroValue.Topwidth not found'
        call json%get('ZeroValue.UseZeroValues', logical_value, found)
        setting%ZeroValue%UseZeroValues = logical_value
        if (.not. found) stop "Error - json file - setting " // 'ZeroValue.UseZeroValues not found'
        call json%get('ZeroValue.Volume', real_value, found)
        setting%ZeroValue%Volume = real_value
        if (.not. found) stop "Error - json file - setting " // 'ZeroValue.Volume not found'



    !% Debug
        !% --- Debug Settings
        call json%get('Debug.Tests', logical_value, found)
        setting%Debug%Tests = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.Tests not found'
        call json%get('Debug.File.adjust', logical_value, found)
        setting%Debug%File%adjust = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.adjust not found'
        call json%get('Debug.File.boundary_conditions', logical_value, found)
        setting%Debug%File%boundary_conditions = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.boundary_conditions not found'
        call json%get('Debug.File.c_library', logical_value, found)
        setting%Debug%File%c_library = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.c_library not found'
        call json%get('Debug.File.define_globals', logical_value, found)
        setting%Debug%File%define_globals = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.define_globals not found'
        call json%get('Debug.File.define_indexes', logical_value, found)
        setting%Debug%File%define_indexes = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.define_indexes not found'
        call json%get('Debug.File.define_keys', logical_value, found)
        setting%Debug%File%define_keys = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.define_keys not found'
        call json%get('Debug.File.define_settings', logical_value, found)
        setting%Debug%File%define_settings = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.define_settings not found'
        call json%get('Debug.File.define_types', logical_value, found)
        setting%Debug%File%define_types = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.define_types not found'
        call json%get('Debug.File.diagnostic_elements', logical_value, found)
        setting%Debug%File%diagnostic_elements = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.diagnostic_elements not found'
        call json%get('Debug.File.discretization', logical_value, found)
        setting%Debug%File%discretization = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.discretization not found'
        call json%get('Debug.File.face', logical_value, found)
        setting%Debug%File%face = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.face not found'
        call json%get('Debug.File.geometry', logical_value, found)
        setting%Debug%File%geometry = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.geometry not found'
        call json%get('Debug.File.interface', logical_value, found)
        setting%Debug%File%interface = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.interface not found'
        call json%get('Debug.File.initial_condition', logical_value, found)
        setting%Debug%File%initial_condition = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.initial_condition not found'
        call json%get('Debug.File.initialization', logical_value, found)
        setting%Debug%File%initialization = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.initialization not found'
        call json%get('Debug.File.jump', logical_value, found)
        setting%Debug%File%jump = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.jump not found'
        call json%get('Debug.File.lowlevel_rk2', logical_value, found)
        setting%Debug%File%lowlevel_rk2 = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.lowlevel_rk2 not found'
        call json%get('Debug.File.network_define', logical_value, found)
        setting%Debug%File%network_define = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.network_define not found'
        call json%get('Debug.File.orifice_elements', logical_value, found)
        setting%Debug%File%orifice_elements = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.orifice_elements not found'
        call json%get('Debug.File.pack_mask_arrays', logical_value, found)
        setting%Debug%File%pack_mask_arrays = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.pack_mask_arrays not found'
        call json%get('Debug.File.partitioning', logical_value, found)
        setting%Debug%File%partitioning = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.partitioning not found'
        call json%get('Debug.File.pump_elements', logical_value, found)
        setting%Debug%File%pump_elements = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.pump_elements not found'
        call json%get('Debug.File.rectangular_channel', logical_value, found)
        setting%Debug%File%rectangular_channel = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.rectangular_channel not found'
        call json%get('Debug.File.trapezoidal_channel', logical_value, found)
        setting%Debug%File%trapezoidal_channel = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.trapezoidal_channel not found'
        call json%get('Debug.File.runge_kutta2', logical_value, found)
        setting%Debug%File%runge_kutta2 = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.runge_kutta2 not found'
        call json%get('Debug.File.timeloop', logical_value, found)
        setting%Debug%File%timeloop = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.timeloop not found'
        call json%get('Debug.File.update', logical_value, found)
        setting%Debug%File%update = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.update not found'
        call json%get('Debug.File.utility_allocate', logical_value, found)
        setting%Debug%File%utility_allocate = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.utility_allocate not found'
        call json%get('Debug.File.utility_deallocate', logical_value, found)
        setting%Debug%File%utility_deallocate = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.utility_deallocate not found'
        call json%get('Debug.File.utility_array', logical_value, found)
        setting%Debug%File%utility_array = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.utility_array not found'
        call json%get('Debug.File.utility_datetime', logical_value, found)
        setting%Debug%File%utility_datetime = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.utility_datetime not found'
        call json%get('Debug.File.utility_interpolate', logical_value, found)
        setting%Debug%File%utility_interpolate = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.utility_interpolate not found'
        call json%get('Debug.File.utility_output', logical_value, found)
        setting%Debug%File%utility_output = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.utility_output not found'
        call json%get('Debug.File.utility_string', logical_value, found)
        setting%Debug%File%utility_string = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.utility_string not found'
        call json%get('Debug.File.utility', logical_value, found)
        setting%Debug%File%utility = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.utility not found'
        call json%get('Debug.File.weir_elements', logical_value, found)
        setting%Debug%File%weir_elements = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.weir_elements not found'
        call json%get('Debug.File.output', logical_value, found)
        setting%Debug%File%output = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.File.output not found'
        call json%get('Debug.FileGroup.all', logical_value, found)
        setting%Debug%FileGroup%all = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.FileGroup.all not found'
        call json%get('Debug.FileGroup.definitions', logical_value, found)
        setting%Debug%FileGroup%definitions = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.FileGroup.definitions not found'
        call json%get('Debug.FileGroup.finalization', logical_value, found)
        setting%Debug%FileGroup%finalization = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.FileGroup.finalization not found'
        call json%get('Debug.FileGroup.geometry', logical_value, found)
        setting%Debug%FileGroup%geometry = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.FileGroup.geometry not found'
        call json%get('Debug.FileGroup.initialization', logical_value, found)
        setting%Debug%FileGroup%initialization = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.FileGroup.initialization not found'
        call json%get('Debug.FileGroup.interface', logical_value, found)
        setting%Debug%FileGroup%interface = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.FileGroup.interface not found'
        call json%get('Debug.FileGroup.output', logical_value, found)
        setting%Debug%FileGroup%output = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.FileGroup.output not found'
        call json%get('Debug.FileGroup.timeloop', logical_value, found)
        setting%Debug%FileGroup%timeloop = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.FileGroup.timeloop not found'
        call json%get('Debug.FileGroup.utility', logical_value, found)
        setting%Debug%FileGroup%utility = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.FileGroup.utility not found'
        call json%get('Debug.Setup', logical_value, found)
        setting%Debug%Setup = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.Setup not found'
        call json%get('Debug.Output', logical_value, found)
        setting%Debug%Output = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Debug.Output not found'

        call def_update_debug_options()

    !% Profiler
        !% --- Profile settings
        call json%get('Profile.YN', logical_value, found)
        setting%Profile%YN = logical_value
        if (.not. found) stop "Error - json file - setting " // 'Profile.YN not found'

     !% finished
        call json%destroy()
        if (json%failed()) stop "JSON failed to destroy"

        if (setting%Debug%File%define_settings) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine def_load_settings
!%
!%==========================================================================
!%==========================================================================
!%
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
!%
!%==========================================================================
!%==========================================================================
!%
end module define_settings
