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
    !%===================================================================
    !% Third Level Types
    !%===================================================================

    !% setting%ACmethod%Anomaly
    type ACmethodAnomalyType
        ! density anomaly correction to handle residual (or non-convergence) of the AC
        logical :: UseDensityCorrection = .false.
        ! Note: setting the lowcutoff to zero can cause small truncation error velocities
        ! to cause waves that buildup over time into a sloshing flow. The lowcutoff avoid this
        real(8) :: DensityHighCutoff = 0.1
        real(8) :: DensityLowCutoff = 1e-10
        ! fraction of residual that is corrected - generally should be 1.0
        real(8) :: FullPipeFactor = 1.0
        ! fraction of residual that is corrected - generally should be 1.0
        real(8) :: OpenPipeFactor = 1.0
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
        ! AC convergence change in absolute L2 norm when cutting off AC
        real(8) :: Habsolute = 1e-5
        ! AC convergence relative change in L2 norm when cutting off AC
        real(8) :: Hrelative = 1e-2
        ! AC convergence change in absolute L2 norm when cutting off AC
        real(8) :: Qabsolute = 1e-5
        ! AC convergence relative change in L2 norm when cutting off AC
        real(8) :: Qrelative = 1e-2
    end type ACmethodConvergenceType

    !% setting%ACmethod%Iter
    type ACmethodIterType
        ! allows more iterations in first time step
        integer :: Firststep = 100
        ! cutoff for AC iterations without convergence
        integer :: Max = 100
        ! reset by code so that itermin * dtau >  dt
        integer :: Min = 3    
    end type ACmethodIterType

    !% setting%ACmethod%Switch
    type ACmethodSwitchType
        ! switch to AC solver if area/areaMax >  0.9
        real(8) :: Area = 0.9
        ! 5% buffer for the switch
        real(8) :: Buffer = 0.05
        ! switch to AC solver if depth/depthMax >  0.9
        real(8) :: Depth = 0.9
    end type ACmethodSwitchType

    !% setting%Adjust%Flowrate
    type AdjustFlowrateType
        logical :: ApplyYN = .true.
        integer :: Approach = vshape
        real(8) :: Coef = 1.0
    end type AdjustFlowrateType

    !% setting%Adjust%Head
    type AdjustHeadType
        logical :: ApplyYN = .true.
        integer :: Approach = vshape_surcharge_only
        real(8) :: Coef = 1.0
    end type AdjustHeadType

    ! !% setting%Adjust%WidthDepth
    ! type AdjustWidthDepthType
    !     logical :: ApplyYN = .true.
    ! end type AdjustWidthDepthType

    !% setting%Output%CommandLine
    type CommandLineType
        logical :: quietYN = .false.
        integer :: interval = 1 !% steps between output
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

    !% setting%Limiter%ArraySize
    type LimiterArraySizeType
        integer :: TemporalInflows = 10
        integer :: TotalInflows = 50
    end type LimiterArraySizeType

    !% setting%Limiter%BC
    type LimiterBCType
        logical :: UseInflowLimiterYN = .true.
        integer :: Approach = FroudeNumber
        ! max value of Fr at inflow
        real(8) :: FroudeInflowMaximum = 1.5
    end type LimiterBCType

    !% setting%Limiter%Channel
    type LimiterChannelType
        real(8) :: LargeDepthFactor = 10.0
    end type LimiterChannelType

    !% setting%Limiter%Dt
    type LimiterDtType
        logical :: UseLimitMinYN = .true.
        real(8) :: Minimum     = 1e-6
    end type LimiterDtType

    !% setting%Limiter%Flowrate
    type LimiterFlowrateType
        logical :: UseFaceVolumeTransportYN = .true.
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
        logical :: UseLimitMaxYN = .true.
        real(8) :: Maximum = 10.0 ! m/s
    end type LimiterVelocityType

    !% setting%Output%Report
    type ReportType
        logical :: useSWMMinpYN = .true.
        logical :: provideYN = .true.
        logical :: suppress_MultiLevel_Output = .false.
        real(8) :: StartTime = 0.0
        real(8) :: TimeInterval = 300.0
        integer :: ThisStep
        integer :: TimeUnits = InHours
    end type ReportType

    !% setting%Time% ...Hydraulics, Hydrology, Dry
    !% these is initialized in define_settings_defaults
    type TimeStepType
       real(8) :: Dt
       real(8) :: LastTime
       real(8) :: NextTime
       integer :: Step
    end type TimeStepType

    !% setting%File%UnitNumber
    type UnitNumberType
        integer :: inp_file
        integer :: out_file
        integer :: rpt_file
        integer :: setting_file
        !integer :: links_input_file
        !integer :: nodes_input_file
        ! integer :: debug_setup_linkR_file
        ! integer :: debug_setup_linkI_file
        ! integer :: debug_setup_nodeR_file
        ! integer :: debug_setup_nodeI_file
        ! integer :: debug_setup_nodeYN_file
        ! !integer :: outputML_combined_file
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

    !% setting%Time%WallClock
    type WallClockType
        integer (kind=8):: Start = 0
        integer (kind=8):: Now  = 0
        integer (kind=8):: End = 0
        integer (kind=8):: CountRate = 0
        integer (kind=8):: HydraulicsStart= 0
        integer (kind=8):: HydraulicsStop = 0
        integer (kind=8):: HydraulicsCumulative = 0
        integer (kind=8):: LoopOutputStart = 0
        integer (kind=8):: LoopOutputStop = 0
        integer (kind=8):: LoopOutputCumulative = 0
        integer (kind=8):: HydrologyStart = 0
        integer (kind=8):: HydrologyStop  = 0
        integer (kind=8):: HydrologyCumulative = 0
        integer (kind=8):: TimeMarchStart = 0
        integer (kind=8):: TimeMarchEnd = 0
        integer (kind=8):: InitializationEnd = 0
        integer (kind=8):: FinalOutputStart = 0
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
        logical :: BIPquick            = .false.
        logical :: boundary_conditions = .false.
        logical :: c_library           = .false.
        logical :: define_globals      = .false.
        logical :: define_indexes      = .false.
        logical :: define_keys         = .false.
        logical :: define_settings     = .false.
        logical :: define_types        = .false.
        logical :: diagnostic_elements = .false.
        logical :: discretization      = .false.
        logical :: face                = .false.
        logical :: finalization        = .false.
        logical :: geometry            = .false.
        logical :: initial_condition   = .false.
        logical :: initialization      = .false.
        logical :: interface           = .false.
        logical :: jump                = .false.
        logical :: lowlevel_rk2        = .false.
        logical :: network_define      = .false.
        logical :: orifice_elements    = .false.
        logical :: output              = .false.
        logical :: pack_mask_arrays    = .false.
        logical :: partitioning        = .false.
        logical :: pump_elements       = .false.
        logical :: rectangular_channel = .false.
        logical :: trapezoidal_channel = .false.
        logical :: runge_kutta2        = .false.
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

    !%===================================================================
    ! Second Level Types
    !%===================================================================

    ! setting%ACmethodType
    type ACmethodType
        real(8) :: dtau = 1.0
        type(ACmethodAnomalyType)      :: Anomaly
        type(ACmethodImplicitCoefType) :: ImplicitCoef
        type(ACmethodCFLType)          :: CFL
        type(ACmethodCelerityType)     :: Celerity
        type(ACmethodConvergenceType)  :: Convergence
        type(ACmethodIterType)         :: Iter
        type(ACmethodSwitchType)       :: Switch
    end type ACmethodType

    ! setting%Adjust
    type AdjustType
        type(AdjustFlowrateType)   :: Flowrate
        type(AdjustHeadType)       :: Head
        !type(AdjustWidthDepthType) :: WidthDepth
    end type AdjustType

    ! setting%BC
    type BCPropertiesType
        integer :: TimeSlotsStored        = 10
        logical :: disableInterpolationYN = .false.
    end type BCPropertiesType

    ! setting%CaseName
    type CaseNameType
        character(256) ::    Long = 'default_CaseName'
        character(16)  ::    Short = 'default'
        character(31)  ::    withTimeStamp
    end type CaseNameType

    ! setting%Constant
    type ConstantType
        real(8) :: gravity = 9.81 ! m^2/s
    end type ConstantType

    !% setting%Debug
    type DebugType
        type(DebugFileYNType) :: File
        type(DebugFileGroupYNType) :: FileGroup
        logical :: SetupYN
        logical :: OutputYN
    end type DebugType

    !% setting%Discretization
    type DiscretizationType
        real(8) :: LinkShortingFactor  = 0.33
        real(8) :: MinElemLengthFactor = 0.50
        integer :: MinElemLengthMethod = ElemLengthAdjust
        real(8) :: NominalElemLength   = 10.0
    end type DiscretizationType

    ! setting%Eps
    type EpsilonType
        !% +- small non-dimensional range for hyd jump discrimination
        real(8) :: FroudeJump = 0.1
        !% Fractional increase in depth under froude limitation
        real(8) :: InflowDepthIncreaseFroudeLimit = 0.1
    end type EpsilonType

    !% setting%FaceInterp
    type FaceInterpType
        integer :: DownJBFaceInterp = dynamic
    end type FaceInterpType

    !% setting%File
    type FileType
        logical              :: UseCommandLineFoldersYN  = .true.
        logical              :: force_folder_creationYN = .true.
        !% standard files and folders
        character(len=256)   :: base_folder = ""
        character(len=256)   :: library_folder = ""
        character(len=256)   :: output_folder= "" !
        character(len=256)   :: output_timestamp_subfolder = ""
        character(len=256)   :: output_temp_subfolder = "temp"
        character(len=256)   :: project_folder = "" ! project path
        character(len=256)   :: setting_file = "" ! path to settings JSON file
        character(len=256)   :: input_kernel = "" ! main part of input file name
        character(len=256)   :: output_kernel= "" ! main part ouf output file name
        character(len=256)   :: inp_file = "" ! path to SWMM input (.inp) file
        character(len=256)   :: rpt_file = "" ! path to SWMM report (.rpt) file
        character(len=256)   :: out_file = "" ! path to SWMM output (.out) file
        
        !% for multi-level output
        character(len=256)   :: outputML_Link_kernel = "link"
        character(len=256)   :: outputML_Node_kernel = "node"
        character(len=256)   :: outputML_combinedfile_kernel = "combined_output" !% filenames that combine across images
        character(len=256)   :: outputML_filename_file = 'output_filenames.txt'  !% list of interim filenames used in combined output
        character(len=256)   :: outputML_control_file = 'control.unf'  !% global parameters that control the processing of combined output
        integer              :: outputML_Ncombined_file_written = 0
        integer              :: outputML_total_timelevels_written = 0

        !% file handling
        integer              :: last_unit = 1000 !% starting point for assigning unit numbers
        type(UnitNumberType) :: UnitNumber

        !% for csv dump output
        ! character(len=256) :: debug_setup_link_folder = ""
        ! character(len=256) :: debug_setup_node_folder = ""
        ! character(len=256) :: debug_output_link_folder = ""
        ! character(len=256) :: debug_output_node_folder = ""
        ! character(len=256) :: debug_output_elemR_folder = ""
        ! character(len=256) :: debug_output_faceR_folder = ""
        ! character(len=256) :: debug_output_summary_folder = ""
        ! character(len=256) :: swmm5_output_link_folder = ""
        ! character(len=256) :: swmm5_output_node_folder = ""
        ! character(len=256) :: links_input_file = "links_input.csv"
        ! character(len=256) :: nodes_input_file = "nodes_input.csv"
        ! character(len=256) :: debug_setup_linkR_file = "linkR.csv"
        ! character(len=256) :: debug_setup_linkI_file = "linkI.csv"
        ! character(len=256) :: debug_setup_nodeR_file = "nodeR.csv"
        ! character(len=256) :: debug_setup_nodeI_file = "nodeI.csv"
        ! character(len=256) :: debug_setup_nodeYN_file = "nodeYN.csv"
        
        !logical :: links_input_file_exist = .false.
        !logical :: nodes_input_file_exist = .false.
    end type FileType

    ! setting%Junction
    type JunctionType
        logical :: isDynamicYN    = .false.
        real(8) :: CFLlimit     = 0.5   !% limiter on CFL to control dynamic junction
        integer :: FunStorageN  = 10    !% number of curve entries for functional storage   
        real(8) :: HeadCoef     = 1.0   !% junction branch head coef for diagnostic junction (must be > 0)
        real(8) :: kFactor      = 0.0   !% default entrance/exit losses at junction branch (use 0.0 as needs debugging)
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
    end type LinkType

    !% setting%Orifice
    type OrificeType
        real(8) :: SharpCrestedWeirCoefficient
        real(8) :: TransverseWeirExponent
        real(8) :: VillemonteCorrectionExponent
    end type OrificeType

    !% setting%Output
    type OutputType
        logical :: UseFileNameFile = .false.
        logical :: OutputElementsExist = .false.
        logical :: OutputFacesExist = .false.
        logical :: Verbose = .true.
        logical :: Warning = .true.
        integer :: LastLevel = 0
        integer :: MaxExpectedLevels = 0
        integer :: NumberOfWriteSteps = 0
        integer :: NumberOfTimeLevelSaved = 0
        integer :: StoredLevels = 100        
        integer :: StoredFileNames = 2
        type(CommandLineType) :: CommandLine
        type(DataOutType) :: DataOut
        type(ReportType) :: Report
    end type OutputType

    !% setting%Partitioning
    type PartitioningType
        integer :: PartitioningMethod = BQuick
    endtype PartitioningType

    !% setting%PreissmannSlot
    type PreissmannSlotType
        integer :: PreissmannSlotMethod = VariableSlot
        real(8) :: CelerityFactor = 1.0
        real(8) :: DesiredTimeStep = 0.01
    end type PreissmannSlotType

    !% setting%Profile
    type ProfileType
        logical :: useYN = .false.
        !logical :: Tests = .false.
        !type(ProfileFileYNType) :: File
        !type(ProfileFileGroupYNType) :: FileGroup
        !logical :: Input
        !logical :: Output
    end type ProfileType

    !% setting%Simulation
    type SimulationType
        logical :: stopAfterInitializionYN = .false.
        logical :: useHydrology = .true.
        logical :: useHydraulics = .true.
    end type SimulationType

    ! setting%SmallVolume
    type SmallVolumeType
        ! Dont using small volumes for weir case.
        ! Needed to be changed later SmallVolumeType
        logical :: UseSmallVolumesYN = .true.
        real(8) :: DepthCutoff = 0.01 ! m
        real(8) :: ManningsN = 0.01
        real(8) :: MinimumArea = 0.005 ! m^2
        real(8) :: MinimumHydRadius = 0.009 ! m
        real(8) :: MinimumPerimeter = 0.52 ! m
        real(8) :: MinimumTopwidth = 0.5 ! m
    end type SmallVolumeType

    ! setting%Solver
    type SolverType
        logical :: PreissmannSlot = .true.
        logical :: QinterpWithLocalHeadGradient = .true.
        logical :: SubtractReferenceHead = .true.
        integer :: MomentumSourceMethod = T00
        integer :: SolverSelect = ETM
        real(8) :: SwitchFractionDn = 0.8
        real(8) :: SwitchFractionUp = 0.9
        real(8) :: QHgradFactor = 0.5
        real(8) :: ReferenceHead = zeroR
        real(8) :: AverageZbottom = zeroR
        real(8) :: MaxZbottom = zeroR
        real(8), dimension(2) :: crk2 = [0.5, 1.0]
    end type SolverType

    type TestCaseType
        logical       :: UseTestCaseYN = .false.
        character(64) :: TestName
    end type TestCaseType

    type SWMMinputType
        real(8) :: ReportStartTimeEpoch = nullvalueR
        real(8) :: ReportTimeInterval = nullvalueR
        real(8) :: StartEpoch = nullvalueR
        real(8) :: EndEpoch = nullvalueR
        real(8) :: WetStep = nullvalueR
        real(8) :: DryStep = nullvalueR
        real(8) :: RouteStep = nullvalueR
        real(8) :: TotalDuration = nullvalueR
    end type SWMMinputType    

    ! setting%Time
    type TimeType
        logical            :: useSWMMinpYN = .true.
        logical            :: matchHydrologyStep = .true.
        real(8)            :: DtTol = 1e-1   !% tolerance when comparing Now time and accumulation of time steps.
        character(14)      :: DateTimeStamp
        integer            :: Step  !% count of the number of steps (either hydrology, hydraulics, or combined)
        real(8)            :: Start
        real(8)            :: Now
        real(8)            :: End
        real(8)            :: StartEpoch = nullvalueR
        real(8)            :: EndEpoch   = nullvalueR
        type(TimeStepType) :: Hydraulics
        type(TimeStepType) :: Hydrology
        type(WallClockType):: WallClock
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
        logical :: ApplyYN = .true.
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

    !%===================================================================
    !% First Level Type (setting)
    !%===================================================================

    type settingType
        logical                  :: JSON_FoundFileYN = .false.
        logical                  :: JSON_CheckAllInputYN = .true.
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
        type(SWMMinputType)      :: SWMMinput ! storage of SWMM *.inp data
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
    subroutine define_settings_default()
        !% -----------------------------------------------------------------
        !% Description:
        !% For settings that where subtypes cannot be initialized in the
        !% define_settings structure (which occurs when a type is used for
        !% more than one setting, e.g. WeirConstantType) the default values 
        !% are provided below
        !% -----------------------------------------------------------------

    setting%Time%Hydraulics%Dt = 10.0
    setting%Time%Hydrology%Dt = 600.0

    setting%Weir%Transverse%WeirExponent = 1.5
    setting%Weir%Transverse%WeirContractionFactor = 0.1
    setting%Weir%Transverse%SideFlowWeirCrestExponent = 1.0
    setting%Weir%Transverse%VillemonteCorrectionExponent = 0.385

    setting%Weir%SideFlow%WeirExponent = 167
    setting%Weir%SideFlow%WeirContractionFactor = 0.1
    setting%Weir%SideFlow%SideFlowWeirCrestExponent = 0.83
    setting%Weir%SideFlow%VillemonteCorrectionExponent = 0.385

    setting%Weir%Vnotch%WeirExponent = 2.5
    setting%Weir%Vnotch%WeirContractionFactor = 1.0
    setting%Weir%Vnotch%SideFlowWeirCrestExponent = 1.0
    setting%Weir%Vnotch%VillemonteCorrectionExponent = 0.385

    setting%Weir%Trapezoidal%WeirExponent = 1.5
    setting%Weir%Trapezoidal%WeirContractionFactor = 1.0
    setting%Weir%Trapezoidal%SideFlowWeirCrestExponent = 1.0
    setting%Weir%Trapezoidal%VillemonteCorrectionExponent = 0.385

    end subroutine define_settings_default
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine define_settings_load()
        !%------------------------------------------------------------------
        !% Description:
        !%    Loads setting values from external JSON file.
        !%------------------------------------------------------------------
            character(kind=json_CK, len=:), allocatable :: c
            real(8) :: real_value
            integer :: integer_value
            logical :: logical_value
            logical :: found
            logical, pointer :: jsoncheck
            type(json_file) :: json
            character(64) :: subroutine_name = 'def_load_settings'
        !%------------------------------------------------------------------
        !% Preliminaries
            if (setting%Debug%File%define_settings) &
                write(*,"(A,i5,A)") '*** enter ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
        !%------------------------------------------------------------------
        !% load the hard-coded default values
        call define_settings_default()

        !% exit if there is no setting.json file
        if (.not. setting%JSON_FoundFileYN) return

        !% set up and read the json file       
        call json%initialize()
        call json%load(filename = trim(setting%File%setting_file))

        !% get the default setting for checking the json file
        jsoncheck => setting%JSON_CheckAllInputYN

        !% look to see if the JSON file was found
        call json%get('JSON_FoundFileYN', logical_value, found)
        if (found) then 
            !% the file is found even if the variable is false!
            setting%JSON_FoundFileYN = .true.
        else
            write(*,"(A)") "Error - json file - setting " // 'JSON_FoundFileYN not found'
            write(*,"(A)") "...This should be first item in json file."
            write(*,"(A)") "...There may be a formatting problem in the json file,"
            write(*,"(A)") "...or this line may be missing. The most common cause "
            write(*,"(A)") "...is the input JSON filename had a typographical error"
            write(*,"(A)") "...in its name or path."
            stop 970984
        end if

        !% get the level of input checking
        call json%get('JSON_CheckAllInputYN', logical_value, found)
        if (found) setting%JSON_CheckAllInputYN= logical_value

        if (jsoncheck) then 
            write(*,"(A)") "All required input from JSON file will be checked"
            write(*,"(A)") "to see if it is present. If you would like to "
            write(*,"(A)") "accept the defaults for missing items, please"
            write(*,"(A)") "change the JSON_CheckAllInputYN to false."
        end if

    !% ACmethod. =====================================================================
        !% --- ACmethod
        !%                       ACmethod.dtau
        call json%get('ACmethod.dtau', real_value, found)
        if (found) setting%ACmethod%dtau = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting "// 'ACmethod.dtau not found'

        !% --- ACmethod.Anomaly
        !%                       Anomaly.UseDensityCorrectionn
        call json%get('ACmethod.Anomaly.UseDensityCorrection', logical_value, found)
        if (found) setting%ACmethod%Anomaly%UseDensityCorrection = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'ACmethod.Anomaly.UseDensityCorrection not found'
  
        !%                       Anomaly.DensityLowCutoff
        call json%get('ACmethod.Anomaly.DensityLowCutoff', real_value, found)
        if (found) setting%ACmethod%Anomaly%DensityLowCutoff = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'ACmethod.Anomaly.DensityLowCutoff not found'

        !%                       Anomaly.DensityHighCutoff              
        call json%get('ACmethod.Anomaly.DensityHighCutoff', real_value, found)
        if (found) setting%ACmethod%Anomaly%DensityHighCutoff = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'ACmethod.Anomaly.DensityHighCutoff not found'
       
        !%                       Anomaly.FullPipeFactor
        call json%get('ACmethod.Anomaly.FullPipeFactor', real_value, found)
        if (found) setting%ACmethod%Anomaly%FullPipeFactor = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'ACmethod.Anomaly.FullPipeFactor not found'
        
        !%                       Anomaly.OpenPipeFactor
        call json%get('ACmethod.Anomaly.OpenPipeFactor', real_value, found)
        if (found) setting%ACmethod%Anomaly%OpenPipeFactor = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'ACmethod.Anomaly.OpenPipeFactor not found'
     
        !% --- ACmethod.ImplicitCoef 
        !%                       ImplicitCoef.a1
        call json%get('ACmethod.ImplicitCoef.a1', real_value, found)
        if (found) setting%ACmethod%ImplicitCoef%a1 = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'ACmethod.ImplicitCoef.a1 not found'
        
        !%                       ImplicitCoef.a2
        call json%get('ACmethod.ImplicitCoef.a2', real_value, found)
        if (found) setting%ACmethod%ImplicitCoef%a2 = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'ACmethod.ImplicitCoef.a2 not found'
        
        !%                       ImplicitCoef.a3
        call json%get('ACmethod.ImplicitCoef.a3', real_value, found)
        if (found) setting%ACmethod%ImplicitCoef%a3 = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'ACmethod.ImplicitCoef.a3 not found'

        !% --- ACmethod.CFL
        !%                       CFL.CFLmax
        call json%get('ACmethod.CFL.CFLmax', real_value, found)
        if (found) setting%ACmethod%CFL%CFLmax = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'ACmethod.CFL.CFLmax not found'
        
        !%                       CFL.CFLsmall
        call json%get('ACmethod.CFL.CFLsmall', real_value, found)
        if (found) setting%ACmethod%CFL%CFLsmall = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'ACmethod.CFL.CFLsmall not found'

        !% --- ACmethod.Celerity
        !%                       Celerity.RC
        call json%get('ACmethod.Celerity.RC', real_value, found)
        if (found) setting%ACmethod%Celerity%RC = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'ACmethod.Celerity.RC not found'

        !% --- ACmethod.Convergence
        !%                       Convergence.Habsolute
        call json%get('ACmethod.Convergence.Habsolute', real_value, found)
        if (found) setting%ACmethod%Convergence%Habsolute = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'ACmethod.Convergence.Habsolute not found'
        
        !%                       Convergence.Hrelative
        call json%get('ACmethod.Convergence.Hrelative', real_value, found)
        if (found) setting%ACmethod%Convergence%Hrelative = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'ACmethod.Convergence.Hrelative not found'
        
        !%                       Convergence.Qabsolute
        call json%get('ACmethod.Convergence.Qabsolute', real_value, found)
        if (found) setting%ACmethod%Convergence%Qabsolute = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'ACmethod.Convergence.Qabsolute not found'
        
        !%                       Convergence.Qrelative
        call json%get('ACmethod.Convergence.Qrelative', real_value, found)
        if (found) setting%ACmethod%Convergence%Qrelative = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'ACmethod.Convergence.Qrelative not found'

        !% --- ACmethod.Iter
        !%                       Iter.Firststep
        call json%get('ACmethod.Iter.Firststep', integer_value, found)
        if (found) setting%ACmethod%Iter%Firststep = integer_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'ACmethod.Iter.Firststep not found'
        
        !%                       Iter.Max
        call json%get('ACmethod.Iter.Max', integer_value, found)
        if (found) setting%ACmethod%Iter%Max = integer_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'ACmethod.Iter.Max not found'
        
        !%                       Iter.Min
        call json%get('ACmethod.Iter.Min', integer_value, found)
        if (found) setting%ACmethod%Iter%Min = integer_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'ACmethod.Iter.Min not found'

        !% --- ACmethod.Switch 
        !%                       Switch.Area
        call json%get('ACmethod.Switch.Area', real_value, found)
        if (found) setting%ACmethod%Switch%Area = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'ACmethod.Switch.Area not found'
        
        !%                       Switch.Buffer
        call json%get('ACmethod.Switch.Buffer', real_value, found)
        if (found) setting%ACmethod%Switch%Buffer = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'ACmethod.Switch.Buffer not found'
        
        !%                       Switch.Depth
        call json%get('ACmethod.Switch.Depth', real_value, found)
        if (found) setting%ACmethod%Switch%Depth = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'ACmethod.Switch.Depth not found'

    !% Adjust. =====================================================================
        !% --- Adjust.Flowrate
        !%                       Flowrate.ApplyYN
        call json%get('Adjust.Flowrate.ApplyYN', logical_value, found)
        if (found) setting%Adjust%Flowrate%ApplyYN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Adjust.Flowrate.Apply not found'
        
        !%                       Flowrate.Approach
        call json%get('Adjust.Flowrate.Approach', c, found)
        if (found) then 
            call util_lower_case(c)
            if (c == 'vshape') then
                setting%Adjust%Flowrate%Approach = vshape
            else
                write(*,"(A)") 'Error - json file - setting.Adjust.Flowrate.Approach of ',trim(c)
                write(*,"(A)") '..is not in allowed options of:'
                write(*,"(A)") '...vshape '
                stop 98705
            end if
        end if
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Adjust.Flowrate.Approach not found'
        
        !%                       Flowrate.Coef
        call json%get('Adjust.Flowrate.Coef', real_value, found)
        if (found) setting%Adjust%Flowrate%Coef = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Adjust.Flowrate.Coef not found'

        !% --- Adjust head settings
        !%                       Head.ApplyYN
        call json%get('Adjust.Head.ApplyYN', logical_value, found)
        if (found) setting%Adjust%Head%ApplyYN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Adjust.Head.Apply not found'
        
        !%                       Head.Approach
        call json%get('Adjust.Head.Approach', c, found)
        if (found) then 
            call util_lower_case(c)
            if (c == 'vshape_surcharge_only') then
                setting%Adjust%Head%approach = vshape_surcharge_only    
            else
                write(*,"(A)") 'Error - json file - settingAdjust.Head.Approach of ',trim(c)
                write(*,"(A)") '..is not in allowed options of:'
                write(*,"(A)") '...vshape_surcharge_only '
                stop 566339
            end if
        end if
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " //  'Adjust.Head.Approach not found'
        
        !%                       Head.Coef
        call json%get('Adjust.Head.Coef', real_value, found)
        if (found) setting%Adjust%Head%Coef = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Adjust.Head.Coef not found'

    !% BC. =====================================================================
        !%                       BC.TimeSlotsStored
        call json%get('BC.TimeSlotsStored', real_value, found)
        if (found) setting%BC%TimeSlotsStored = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'BC.TimeSlotsStored not found'
        
        !%                       BC.disableInterpolationYN
        call json%get('BC.disableInterpolationYN', logical_value, found)
        if (found) setting%BC%disableInterpolationYN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'BC.disableInterpolationYN not found'

    !% Case. =====================================================================
        !%                       CaseName.Long
        call json%get('CaseName.Long', c, found)
        if (found) setting%CaseName%Long = trim(c)
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'CaseName.Long not found' 
        
        !%                       CaseName.Short
        call json%get('CaseName.Short', c, found)
        if (found) setting%CaseName%Short = trim(c)
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'CaseName.Short not found'
        
        !% do not read           CaseName.withTimeStamp

    !% Constant. =====================================================================
        !%                       Constant.gravity
        call json%get('Constant.gravity', real_value, found)
        if (found) setting%Constant%gravity = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Constant.gravity not found'

    !% Discretization. =====================================================================
        !% -- Nominal element length adjustment
        !%                      Discretization.LinkShortingFactor
        call json%get('Discretization.LinkShortingFactor', real_value, found)
        if (found) setting%Discretization%LinkShortingFactor = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Discretization.LinkShortingFactor not found'
        
        !%                       Discretization.MinElemLengthFactor
        call json%get('Discretization.MinElemLengthFactor', real_value, found)
        if (found) setting%Discretization%MinElemLengthFactor = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Discretization.MinElemLengthFactor not found'
        
        !%                       Discretization.MinElemLengthMethod
        call json%get('Discretization.MinElemLengthMethod', c, found)
        if (found) then 
            call util_lower_case(c)
            if (c == 'elemlengthadjust') then
                setting%Discretization%MinElemLengthMethod = ElemLengthAdjust
            elseif (c == 'rawelemlength') then
                setting%Discretization%MinElemLengthMethod = RawElemLength
            else
                write(*,"(A)") 'Error - json file - setting.Discretization.MinElemLengthMethod of ',trim(c)
                write(*,"(A)") '..is not in allowed options of:'
                write(*,"(A)") '...elemlengthadjust,  rawelemlength'
                stop 33986
            end if
        end if
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Discretization.MinElemLengthMethod not found'
  
        !%                       Discretization.NominalElemLength
        call json%get('Discretization.NominalElemLength', real_value, found)
        if (found) setting%Discretization%NominalElemLength = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " //'Discretization.NominalElemLength not found'

    !% Eps. =====================================================================
        !% --- Eps Settings
        !%                       Eps.FroudeJump
        call json%get('Eps.FroudeJump', real_value, found)
        if (found) setting%Eps%FroudeJump = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Eps.FroudeJump not found'
        
        !%                       Eps.InflowDepthIncreaseFroudeLimit
        call json%get('Eps.InflowDepthIncreaseFroudeLimit', real_value, found)
        if (found) setting%Eps%InflowDepthIncreaseFroudeLimit = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Eps.InflowDepthIncreaseFroudeLimit not found'

    !% FaceInterp. =====================================================================
        !%                       FaceInterp.DownJBFaceInterp
        call json%get('FaceInterp.DownJBFaceInterp', c, found)
        if (found) then 
            call util_lower_case(c)
            if (c == 'static') then
                setting%FaceInterp%DownJBFaceInterp = static
            else if (c == 'dynamic') then
                setting%FaceInterp%DownJBFaceInterp = dynamic
            else
                write(*,"(A)") 'Error - json file - setting.FaceInterp.DownJBFaceInterp of ',trim(c)
                write(*,"(A)") '..is not in allowed options of:'
                write(*,"(A)") '...static, dynamic'
                stop 993895
            end if
        end if
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'FaceInterp.DownJBFaceInterp not found'

    !% File. =====================================================================
        !% --- (filenames and paths should not be read)
        !%                       File.UseCommandLineFoldersYN
        call json%get('File.UseCommandLineFoldersYN', logical_value, found)
        if (found) setting%File%UseCommandLineFoldersYN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'File.UseCommandLineFoldersYN not found'

        !%                       File.force_folder_creationYN
        call json%get('File.force_folder_creationYN', logical_value, found)
        if (found) setting%File%force_folder_creationYN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'File.force_folder_creationYN not found'

        if (.not. setting%File%UseCommandLineFoldersYN) then
            !%                       File.base_folder
            call json%get('File.base_folder', c, found)
            if (found) setting%File%base_folder = trim(c)
            if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting "// 'File.base_folder not found'
        
            !%                       File.library_folder
            call json%get('File.library_folder', c, found)
            if (found) setting%File%library_folder = trim(c)
            if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting "// 'File.library_folder not found'
    
            !%                       File.output_folder
            call json%get('File.output_folder', c, found)
            if (found) setting%File%output_folder = trim(c)
            if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting "// 'File.output_folder not found'

            !%                       File.project_folder
            call json%get('File.project_folder', c, found)
            if (found) setting%File%project_folder = trim(c)
            if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting "// 'File.project_folder not found'

            !%                       File.inp_file
            call json%get('File.inp_file', c, found)
            if (found) setting%File%inp_file = trim(c)
            if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting "// 'File.inp_file not found'
        end if
        !%                       File.output_temp_subfolder
        call json%get('File.output_temp_subfolder', c, found)
        if (found) setting%File%output_temp_subfolder = trim(c)
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting "// 'File.output_temp_subfolder not found'

        !% do not read           File.output_temp_subfolder
        !% do not read           File.setting_file
        !% do not read           File.rpt_file
        !% do not read           File.out_file
        !% do not read           File.outputML... (any of these)
        !% do not read           File.last_unit
        !% do not read           File.UnitNumber...

    !% Junctions. =====================================================================
        !%                       Junction.isDynamicYN
        call json%get('Junction.isDynamicYN', logical_value, found)
        if (found) setting%Junction%isDynamicYN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Limiter.Junction.isDynamicYN not found'

        !%                       Junction.CFLlimit
        call json%get('Junction.CFLlimit', real_value, found)
        if (found) setting%Junction%CFLlimit = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Junction.CFLlimit not found'

        !%                       Junction.FunStorageN
        call json%get('Junction.FunStorageN', integer_value, found)
        if (found) setting%Junction%FunStorageN = integer_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Junction.CFLlimit not found'
        
        !%                       Junction.HeadCoef
        call json%get('Junction.HeadCoef', real_value, found)
        if (found) setting%Junction%HeadCoef = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Junction.HeeadCoef not found'

        !%                       Junction.kFactor
        call json%get('Junction.kFactor', real_value, found)
        if (found) setting%Junction%kFactor = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Junction.kFactor not found'

   
    !% Limiter. =====================================================================
        !%                       Limiter.ArraySize.TemporalInflows
        call json%get('Limiter.ArraySize.TemporalInflows', integer_value, found)
        if (found) setting%Limiter%ArraySize%TemporalInflows = integer_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Limiter.ArraySize.TemporalInflows not found'

        !%                       Limiter.ArraySize.TotalInflows
        call json%get('Limiter.ArraySize.TotalInflows', integer_value, found)
        if (found) setting%Limiter%ArraySize%TotalInflows = integer_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Limiter.ArraySize.TotalInflows not found'

        !%                       Limiter.BC.UseInflowLimiterYN
        call json%get('Limiter.BC.UseInflowLimiterYN', logical_value, found)
        if (found) setting%Limiter%BC%UseInflowLimiterYN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Limiter.BC.UseInflowLimiterYN not found'

        !%                       Limiter.BC.Approach
        call json%get('Limiter.BC.Approach', c, found)
        if (found) then 
            call util_lower_case(c)
            if (c == 'froudenumber') then
                setting%Limiter%BC%Approach = FroudeNumber
            else
                write(*,"(A)") 'Error - json file - setting.Limiter.BC.Approach of ',trim(c)
                write(*,"(A)") '..is not in allowed options of:'
                write(*,"(A)") '...froudenumber'
                stop 837735
            end if
        end if
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Limiter.BC.Approach not found'
        
        !%                       Limiter.BC.FroudeInflowMaximum
        call json%get('Limiter.BC.FroudeInflowMaximum', real_value, found)
        if (found) setting%Limiter%BC%FroudeInflowMaximum = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Limiter.BC.FroudeInflowMaximum not found'
        
        !%                       Limiter.Channel.LargeDepthFactor
        call json%get('Limiter.Channel.LargeDepthFactor', real_value, found)
        if (found) setting%Limiter%Channel%LargeDepthFactor = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Limiter.Channel.LargeDepthFactor not found'
        
        !%                       Limiter.Dt.Minimum
        call json%get('Limiter.Dt.Minimum', real_value, found)
        if (found) setting%Limiter%Dt%Minimum = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Limiter.Dt.Minimum not found'
        
        !%                       Limiter.Dt.UseLimitMinYN
        call json%get('Limiter.Dt.UseLimitMinYN', logical_value, found)
        if (found) setting%Limiter%Dt%UseLimitMinYN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Limiter.Dt.UseLimitMinYN not found'

        !%                       Limiter.Flowrate.FaceVolumeTransport
        call json%get('Limiter.Flowrate.FaceVolumeTransport', real_value, found)
        if (found) setting%Limiter%Flowrate%FaceVolumeTransport = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Limiter.Flowrate.FaceVolumeTransport not found'
        
        !%                       Limiter.Flowrate.UseFaceVolumeTransportYN
        call json%get('Limiter.Flowrate.UseFaceVolumeTransportYN', logical_value, found)
        if (found) setting%Limiter%Flowrate%UseFaceVolumeTransportYN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Limiter.Flowrate.UseFaceVolumeTransportYN not found'

        !%                       Limiter.InterpWeight.Maximum
        call json%get('Limiter.InterpWeight.Maximum', real_value, found)
        if (found) setting%Limiter%InterpWeight%Maximum = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Limiter.InterpWeight.Maximum not found'
        
        !%                       Limiter.InterpWeight.Minimum
        call json%get('Limiter.InterpWeight.Minimum', real_value, found)
        if (found) setting%Limiter%InterpWeight%Minimum = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Limiter.InterpWeight.Minimum not found'

        !%                       Limiter.Velocity.Maximum
        call json%get('Limiter.Velocity.Maximum', real_value, found)
        if (found) setting%Limiter%Velocity%Maximum = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Limiter.Velocity.Maximum not found'
        
        !%                       Limiter.Velocity.UseLimitMaxYN
        call json%get('Limiter.Velocity.UseLimitMaxYN', logical_value, found)
        if (found) setting%Limiter%Velocity%UseLimitMaxYN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Limiter.Velocity.UseLimitMaxYN not found'


    !% Link. =====================================================================
        !%                      Link.DefaultInitDepthType       
        call json%get('Link.DefaultInitDepthType', c, found)
        if (found) then            
            call util_lower_case(c)
            if (c == 'linear') then
                setting%Link%DefaultInitDepthType = LinearlyVarying
            else if (c == 'uniform') then
                setting%Link%DefaultInitDepthType = Uniform
            else if (c == 'exponential') then
                setting%Link%DefaultInitDepthType = ExponentialDecay
            else
                write(*,"(A)") 'Error - json file - setting.Link.DefaultInitDepthType of ',trim(c)
                write(*,"(A)") '..is not in allowed options of:'
                write(*,"(A)") '... linear, uniform, exponential'
                stop 93775
            end if
        end if
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Link.DefaultInitDepthType not found'

    !% Orifice. =====================================================================
        !%                       SharpCrestedWeirCoefficeint
        call json%get('Orifice.SharpCrestedWeirCoefficient', real_value, found)
        if (found) setting%Orifice%SharpCrestedWeirCoefficient = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - setting " // 'Orifice.SharpCrestedWeirCoefficient not found'
        
        !%                       TransverseWeirExponent
        call json%get('Orifice.TransverseWeirExponent', real_value, found)
        if (found) setting%Orifice%TransverseWeirExponent = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - setting " // 'Orifice.TransverseWeirExponent not found'
        
        !%                       VillemonteCorrectionExponent
        call json%get('Orifice.VillemonteCorrectionExponent', real_value, found)
        if (found) setting%Orifice%VillemonteCorrectionExponent = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - setting " // 'Orifice.VillemonteCorrectionExponent not found'

    !% Output. =====================================================================
        !%                       Output.UseFileNameFile
        call json%get('Output.UseFileNameFile', logical_value, found)
        if (found)setting%Output%UseFileNameFile = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.UseFileNameFile not found'
       
        !% do not read           OutputElementsExist
        !% do not read           OutputFacesExist

        !%                       Output.Verbose
        call json%get('Output.Verbose', logical_value, found)
        if (found) setting%Output%Verbose = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.Verbose not found'

        !%                       Output.Warning
        call json%get('Output.Warning', logical_value, found)
        if (found) setting%Output%Warning = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.Warning not found'
        
        !% do not read           Output.LastLeve
        !% do not read           Output.MaxExpectedLevels

        !%                       Output.StoredLevels
        call json%get('Output.StoredLevels', integer_value, found)
        if (found) setting%Output%StoredLevels = integer_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.StoredLevels not found'

        !%                       Output.StoredFileNames
        call json%get('Output.StoredFileNames', integer_value, found)
        if (found) setting%Output%StoredFileNames = integer_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.StoredFileNames not found'

        !% --- Output.CommandLine
        !%                       CommandLine.quietYN
        call json%get('Output.CommandLine.quietYN', logical_value, found)
        if (found) setting%Output%CommandLine%quietYN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.CommandLine.quietYN not found'
        
        !%                       CommandLine.interval
        call json%get('Output.CommandLine.interval', integer_value, found)
        if (found) setting%Output%CommandLine%interval = integer_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.CommandLine.interval not found'

        !% --- Output.DataOut
        !%                       Dataout.isAreaOut                   
        call json%get('Output.DataOut.isAreaOut', logical_value, found)
        if (found) setting%Output%DataOut%isAreaOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isAreaOut not found'
        
        !%                       Dataout.isDepthOut
        call json%get('Output.DataOut.isDepthOut', logical_value, found)
        if (found) setting%Output%DataOut%isDepthOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isDepthOut not found'
        
        !%                       Dataout.isFlowrateOut
        call json%get('Output.DataOut.isFlowrateOut', logical_value, found)
        if (found) setting%Output%DataOut%isFlowrateOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isFlowrateOut not found'
        
        !%                       Dataout.isFroudeNumberOut
        call json%get('Output.DataOut.isFroudeNumberOut', logical_value, found)
        if (found) setting%Output%DataOut%isFroudeNumberOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isFroudeNumberOut not found'
        
        !%                       Dataout.isHeadOut
        call json%get('Output.DataOut.isHeadOut', logical_value, found)
        if (found) setting%Output%DataOut%isHeadOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isHeadOut not found'
        
        !%                       Dataout.isHydRadiusOut
        call json%get('Output.DataOut.isHydRadiusOut', logical_value, found)
        if (found) setting%Output%DataOut%isHydRadiusOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isHydRadiusOut not found'
        
        !%                       Dataout.isPerimeterOut
        call json%get('Output.DataOut.isPerimeterOut', logical_value, found)
        if (found) setting%Output%DataOut%isPerimeterOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isPerimeterOut not found'
        
        !%                       Dataout.isSlotWidthOut
        call json%get('Output.DataOut.isSlotWidthOut', logical_value, found)
        if (found) setting%Output%DataOut%isSlotWidthOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isSlotWidthOut not found'
        
        !%                       Dataout.isSlotDepthOut
        call json%get('Output.DataOut.isSlotDepthOut', logical_value, found)
        if (found) setting%Output%DataOut%isSlotDepthOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isSlotDepthOut not found'
        
        !%                       Dataout.isTopWidthOut
        call json%get('Output.DataOut.isTopWidthOut', logical_value, found)
        if (found) setting%Output%DataOut%isTopWidthOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isTopWidthOut not found'
        
        !%                       Dataout.isVelocityOut
        call json%get('Output.DataOut.isVelocityOut', logical_value, found)
        if (found) setting%Output%DataOut%isVelocityOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isVelocityOut not found'
        
        !%                       Dataout.isVolume
        call json%get('Output.DataOut.isVolumeOut', logical_value, found)
        if (found) setting%Output%DataOut%isVolumeOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isVolumeOut not found'

        !%                       Dataout.isWaveSpeedOut
        call json%get('Output.DataOut.isWaveSpeedOut', logical_value, found)
        if (found) setting%Output%DataOut%isWaveSpeedOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isWaveSpeedOut not found'
             
        !% --- Report settings
        !%                       Report.useSWMMinpYN
        call json%get('Output.Report.useSWMMinpYN', logical_value, found)
        if (found) setting%Output%Report%useSWMMinpYN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.Report.useSWMMinpYN not found'

        !%                       Report.provideYN
        call json%get('Output.Report.provideYN', logical_value, found)
        if (found) setting%Output%Report%provideYN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.Report.provideYN not found'

        !%                       Report.suppress_MultiLevel_Output
        call json%get('Output.Report.suppress_MultiLevel_Output', logical_value, found)
        if (found) setting%Output%Report%suppress_MultiLevel_Output = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.Report.suppress_MultiLevel_Output not found'

        !%                       Report.StartTime
        call json%get('Output.Report.StartTime', real_value, found)
        if (found) setting%Output%Report.StartTime = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.Report.StartTime not found'

        !%                       Report.TimeInterval
        call json%get('Output.Report.TimeInterval', real_value, found)
        if (found) setting%Output%Report.TimeInterval = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.Report.TimeInterval not found'

        !% do not read           Report.ThisStep
    
        !%                       Report.TimeUnits
        call json%get('Output.Report.TimeUnits', c, found)
        if (found) then
            call util_lower_case(c)
            if (c == 'seconds') then
                setting%Output%Report%TimeUnits = InSeconds
            else if (c == 'minutes') then
                setting%Output%Report%TimeUnits = InMinutes
            else if (c == 'hours') then
                setting%Output%Report%TimeUnits = InHours
            else if (c == 'days') then
                setting%Output%Report%TimeUnits = InDays
            else
                write(*,"(A)") 'Error - json file - setting.Output.Report.TimeUnits of ',trim(c)
                write(*,"(A)") '..is not in allowed options of:'
                write(*,"(A)") '... seconds, minutes, hours, days'
                stop 386668
            end if
        end if
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.Report.TimeUnits'

    !% Partitioning.  =====================================================================
        !%                       Partitioning.PartitioningMethod
        call json%get('Partitioning.PartitioningMethod', c, found)
        if (found) then 
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
                write(*,"(A)") 'Error - json file - setting.Partitioning.PartitioningMethod of ',trim(c)
                write(*,"(A)") '..is not in allowed options of:'
                write(*,"(A)") '... seconds, minutes, hours, days'
                stop 73785
            end if
        end if
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Partitioning.PartitioningMethod not found'
  
    !% PreissmannSlot. =====================================================================
        !%                       PreissmannSlotMethod
        call json%get('PreissmannSlot.PreissmannSlotMethod', c, found)
        if (found) then 
            call util_lower_case(c)
            if (c == 'staticslot') then
                setting%PreissmannSlot%PreissmannSlotMethod = StaticSlot
            else if (c == 'variableslot') then
                setting%PreissmannSlot%PreissmannSlotMethod = VariableSlot
            else
                write(*,"(A)") 'Error - json file - setting.PreissmannSlot%PreissmannSlotMethod of ',trim(c)
                write(*,"(A)") '..is not in allowed options of:'
                write(*,"(A)") '... staticslot, variableslot'
                stop 22496
            end if
        end if
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'PreissmannSlot.PreissmannSlotMethod not found'
        
        !%                      CelerityFactor          
        call json%get('PreissmannSlot.CelerityFactor', real_value, found)
        if (found) setting%PreissmannSlot%CelerityFactor = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'PreissmannSlot.CelerityFactor not found'

        !%                      DesiredTimeStep
        call json%get('PreissmannSlot.DesiredTimeStep', real_value, found)
        if (found) setting%PreissmannSlot%DesiredTimeStep = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'PreissmannSlot.DesiredTimeStep not found'

    !% Profile. =====================================================================
        !%                       Profile.useYN
        call json%get('Profile.useYN', logical_value, found)
        if (found) setting%Profile%useYN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Profile.useYN not found'

    !% Simulation. =====================================================================
        !%                       useHydrology
        call json%get('Simulation.useHydrology', logical_value, found)
        if (found) setting%Simulation%useHydrology = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Simulation.useHydrology not found'
        
        !%                       useHydraulics
        call json%get('Simulation.useHydraulics', logical_value, found)
        if (found) setting%Simulation%useHydraulics = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Simulation.useHydraulics not found'
  
    !% SmallVolume. =====================================================================
        !%                       UseSmallVolumesYN
        call json%get('SmallVolume.UseSmallVolumesYN', logical_value, found)
        if (found) setting%SmallVolume%UseSmallVolumesYN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'SmallVolume.UseSmallVolumesYN not found'
        
        !%                       DepthCutoff
        call json%get('SmallVolume.DepthCutoff', real_value, found)
        if (found) setting%SmallVolume%DepthCutoff = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'SmallVolume.DepthCutoff not found'
        
        !%                       ManningsN
        call json%get('SmallVolume.ManningsN', real_value, found)
        if (found) setting%SmallVolume%ManningsN = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'SmallVolume.ManningsN not found'
        
        !%                       MinimumArea
        call json%get('SmallVolume.MinimumArea', real_value, found)
        if (found) setting%SmallVolume%MinimumArea = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'SmallVolume.MinimumArea not found'
       
        !%                       MinimumHydRadius
        call json%get('SmallVolume.MinimumHydRadius', real_value, found)
        if (found) setting%SmallVolume%MinimumHydRadius = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'SmallVolume.MinimumHydRadius not found'
        
        !%                       MinimumPerimeter
        call json%get('SmallVolume.MinimumPerimeter', real_value, found)
        if (found) setting%SmallVolume%MinimumPerimeter = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'SmallVolume.MinimumPerimeter not found'
        
        !%                       MinimumTopwidth
        call json%get('SmallVolume.MinimumTopwidth', real_value, found)
        if (found)  setting%SmallVolume%MinimumTopwidth = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'SmallVolume.MinimumTopwidth not found'
        

    !% Solver. =====================================================================
        !%                       Solver.PreissmannSlot
        call json%get('Solver.PreissmannSlot', logical_value, found)
        if (found) setting%Solver%PreissmannSlot = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Solver.PreissmannSlot not found'
     
        !%                       Solver.QinterpWithLocalHeadGradient
        call json%get('Solver.QinterpWithLocalHeadGradient', logical_value, found)
        if (found) setting%Solver%QinterpWithLocalHeadGradient = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Solver.QinterpWithLocalHeadGradient not found'
        
        !%                       Solver.QHgradFactor
        call json%get('Solver.QHgradFactor', real_value, found)
        if (found)  setting%Solver%QHgradFactor = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Solver.QHgradFactor not found'

        !%                       Solver.SubtractReferenceHead
        call json%get('Solver.SubtractReferenceHead', logical_value, found)
        if (found) setting%Solver%SubtractReferenceHead = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Solver.SubtractReferenceHead not found'
        
        !%                       Solver.MomentumSourceMethod                  
        call json%get('Solver.MomentumSourceMethod', c, found)
        if (found) then 
            call util_lower_case(c)
            if (c == 't00') then
                setting%Solver%MomentumSourceMethod = T00
            else if (c == 't10') then
                setting%Solver%MomentumSourceMethod = T10
            else if (c == 't20') then
                setting%Solver%MomentumSourceMethod = T20
            else
                write(*,"(A)") 'Error - json file - setting.Solver.MomentumSourceMethod of ',trim(c)
                write(*,"(A)") '..is not in allowed options of:'
                write(*,"(A)") '... T00, T10, T20'
                stop 110985
            end if
        end if
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Solver.MomentumSourceMethod not found'

        !%                       Solver.SolverSelect
        call json%get('Solver.SolverSelect', c, found)
        call util_lower_case(c)
        if (found) then 
            if (c == 'etm') then
                setting%Solver%SolverSelect = ETM
            else if (c == 'etm_ac') then
                setting%Solver%SolverSelect = ETM_AC
            else if (c == 'ac') then
                setting%Solver%SolverSelect = AC
            else
                write(*,"(A)") 'Error - json file - setting.Solver.SolverSelect of ',trim(c)
                write(*,"(A)") '..is not in allowed options of:'
                write(*,"(A)") '... ETM, ETM_AC, AC'
                stop 937546
            end if
        end if
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Solver.SolverSelect not found'

        !%                       Solver.SwitchFractionDn
        call json%get('Solver.SwitchFractionDn', real_value, found)
        if (found)  setting%Solver%SwitchFractionDn = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Solver.SwitchFractionDn not found'
        
        !%                       Solver.SwitchFractionU
        call json%get('Solver.SwitchFractionUp', real_value, found)
        if (found)   setting%Solver%SwitchFractionUp = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Solver.SwitchFractionUp not found'

        !% do not read           Solver.crk2

    !% TestCase.  =====================================================================
        !%                       TestCase.UseTestCaseYN
        call json%get('TestCase.UseTestCaseYN', logical_value, found)
        if (found) setting%TestCase%UseTestCaseYN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'TestCase.UseTestCaseYN not found'
        
        !%                       TestCase.TestName
        if (setting%TestCase%UseTestCaseYN) then 
            call json%get('TestCase.TestName', c, found)
            if (found) setting%TestCase%TestName = trim(c)
            if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'TestCase.TestName not found'           
        end if

    !% SWMMinput.  =====================================================================
        !% do not read           ReportStartTimeEpoch
        !% do not read           ReportTimeInterval
        !% do not read           StartEpoch
        !% do not read           EndEpoch
        !% do not read           WetStep
        !% do not read           DryStep
        !% do not read           RouteStep
        !% do not read           TotalDuratoin
        
        
    !% Time.  =====================================================================
        !%                       useSWMMinpYN
        call json%get('Time.useSWMMinpYN', logical_value, found)
        if (found) setting%Time%useSWMMinpYN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Time.useSWMMinpYN not found'

        !%                       matchHydrologyStep
        call json%get('Time.matchHydrologyStep', logical_value, found)
        if (found) setting%Time%matchHydrologyStep = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Time.matchHydrologyStep not found'

        !%                      Time.DtTol
        call json%get('Time.DtTol', real_value, found)
        if (found) setting%Time%DtTol = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Time.DtTol not found'

        !% do not read          Time.DateTimeStamp
        !% do not read          Time.Start
        !% do not read          Time.Now
        !% do not read          Time.End

        !%                      Time.StartEpoch
        call json%get('Time.StartEpoch', real_value, found)
        if (found) setting%Time%StartEpoch = real_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Time.StartEpoch not found'
    
        !%                      Time.EndEpoch
        call json%get('Time.EndEpoch', real_value, found)
        if (found) setting%Time%EndEpoch = real_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Time.EndEpoch not found'

        !%                      Time.Hydraulics.Dt
        call json%get('Time.Hydraulics.Dt', real_value, found)
        if (found) setting%Time%Hydraulics%Dt = real_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Time.Hydraulics.Dt not found'
    
        !% do not read          Time.Hydraulics.LastTime
        !% do not read          Time.Hydraulics.NextTime
        !% do not read          Time.Hydraulics.Step

        !%                      Time.Hydrology.Dt
        if (found) call json%get('Time.Hydrology.Dt', real_value, found)
        setting%Time%Hydrology%Dt = real_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Time.Hydrology.Dt not found'

        !% do not read          Time.Hydrology.LastTime
        !% do not read          Time.Hydrology.NextTime
        !% do not read          Time.Hydrology.Step

        !% do not read          Time.WallClock.ClockStart
        !% do not read          Time.WallClock.ClockLoopStart
        !% do not read          Time.WallClock.ClockNow
        !% do not read          Time.WallClock.ClockCountRate

        !% do not read          Time.CPU.EpochStartSeconds
        !% do not read          Time.CPU.EpochNowSeconds
        !% do not read          TIme.CPU.EpochFinishSeconds

    !% Weir.=====================================================================
        !%                       Transverse.WeirExponent
        call json%get('Weir.Transverse.WeirExponent', real_value, found)
        if (found) setting%Weir%Transverse%WeirExponent = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Weir.Transverse.WeirExponent not found'
        
        !%                       Transverse.WeirContractionFactor
        call json%get('Weir.Transverse.WeirContractionFactor', real_value, found)
        if (found) setting%Weir%Transverse%WeirContractionFactor = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Weir.Transverse.WeirContractionFactor not found'
        
        !%                       Transverse.SideFlowWeirCrestExponent
        call json%get('Weir.Transverse.SideFlowWeirCrestExponent', real_value, found)
        if (found) setting%Weir%Transverse%SideFlowWeirCrestExponent = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Weir.Transverse.SideFlowWeirCrestExponent not found'
        
        !%                       Transverse.VillemonteCorrectionExponent
        call json%get('Weir.Transverse.VillemonteCorrectionExponent', real_value, found)
        if (found) setting%Weir%Transverse%VillemonteCorrectionExponent = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Weir.Transverse.VillemonteCorrectionExponent not found'

        !%                       SideFlow.WeirExponent
        call json%get('Weir.SideFlow.WeirExponent', real_value, found)
        if (found) setting%Weir%SideFlow%WeirExponent = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Weir.SideFlow.WeirExponent not found'
        
        !%                       SideFlow.WeirContractionFactor
        call json%get('Weir.SideFlow.WeirContractionFactor', real_value, found)
        if (found) setting%Weir%SideFlow%WeirContractionFactor = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Weir.SideFlow.WeirContractionFactor not found'
        
        !%                       SideFlow.SideFlowWeirCrestExponent
        call json%get('Weir.SideFlow.SideFlowWeirCrestExponent', real_value, found)
        if (found) setting%Weir%SideFlow%SideFlowWeirCrestExponent = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Weir.SideFlow.SideFlowWeirCrestExponent not found'
        
        !%                       SideFlow.VillemonteCorrectionExponent
        call json%get('Weir.SideFlow.VillemonteCorrectionExponent', real_value, found)
        if (found) setting%Weir%SideFlow%VillemonteCorrectionExponent = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Weir.SideFlow.VillemonteCorrectionExponent not found'

        !%                       VNotch.WeirExponent
        call json%get('Weir.VNotch.WeirExponent', real_value, found)
        if (found) setting%Weir%VNotch%WeirExponent = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Weir.VNotch.WeirExponent not found'
        
        !%                       VNotch.WeirContractionFactor
        call json%get('Weir.VNotch.WeirContractionFactor', real_value, found)
        if (found) setting%Weir%VNotch%WeirContractionFactor = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Weir.VNotch.WeirContractionFactor not found'
        
        !%                       VNotch.SideFlowWeirCrestExponent
        call json%get('Weir.VNotch.SideFlowWeirCrestExponent', real_value, found)
        if (found) setting%Weir%VNotch%SideFlowWeirCrestExponent = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Weir.VNotch.SideFlowWeirCrestExponent not found'
        
        !%                       VNotch.VillemonteCorrectionExponent
        call json%get('Weir.VNotch.VillemonteCorrectionExponent', real_value, found)
        setting%Weir%VNotch%VillemonteCorrectionExponent = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Weir.VNotch.VillemonteCorrectionExponent not found'

        !%                       Trapezoidal.WeirExponent
        call json%get('Weir.Trapezoidal.WeirExponent', real_value, found)
        if (found) setting%Weir%Trapezoidal%WeirExponent = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Weir.Trapezoidal.WeirExponent not found'
        
        !%                       Trapezoidal.WeirContractionFactor
        call json%get('Weir.Trapezoidal.WeirContractionFactor', real_value, found)
        if (found) setting%Weir%Trapezoidal%WeirContractionFactor = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Weir.Trapezoidal.WeirContractionFactor not found'
        
        !%                       Trapezoidal.SideFlowWeirCrestExponent
        call json%get('Weir.Trapezoidal.SideFlowWeirCrestExponent', real_value, found)
        if (found) setting%Weir%Trapezoidal%SideFlowWeirCrestExponent = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Weir.Trapezoidal.SideFlowWeirCrestExponent not found'
        
        !%                       Trapezoidal.VillemonteCorrectionExponent
        call json%get('Weir.Trapezoidal.VillemonteCorrectionExponent', real_value, found)
        if (found) setting%Weir%Trapezoidal%VillemonteCorrectionExponent = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Weir.Trapezoidal.VillemonteCorrectionExponent not found'

    !% VariableDt. =====================================================================
        !%                       ApplyYN
        call json%get('VariableDT.ApplyYN', logical_value, found)
        if (found) setting%VariableDT%ApplyYN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'VariableDT.Apply not found'
        
        !%                       CFL_hi_max
        call json%get('VariableDT.CFL_hi_max', real_value, found)
        if (found) setting%VariableDT%CFL_hi_max = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'VariableDT.CFL_hi_max not found'
        
        !%                       CFL_target
        call json%get('VariableDT.CFL_target', real_value, found)
        if (found) setting%VariableDT%CFL_target = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'VariableDT.CFL_target not found'
        
        !%                       CFL_lo_max
        call json%get('VariableDT.CFL_lo_max', real_value, found)
        if (found) setting%VariableDT%CFL_lo_max = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'VariableDT.CFL_lo_max not found'
        
        !%                       decreaseFactor
        call json%get('VariableDT.decreaseFactor', real_value, found)
        if (found) setting%VariableDT%decreaseFactor = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'VariableDT.decreaseFactor not found'
        
        !%                       increaseFactor
        call json%get('VariableDT.increaseFactor', real_value, found)
        if (found) setting%VariableDT%increaseFactor = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'VariableDT.increaseFactor not found'
        
        !%                       NstepsForCheck
        call json%get('VariableDT.NstepsForCheck', integer_value, found)
        if (found) setting%VariableDT%NstepsForCheck = integer_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'VariableDT.NstepsForCheck not found'

        !% do not read           LastCheckStep

    !% ZeroValue.=====================================================================
        !%                       UseZeroValues
        call json%get('ZeroValue.UseZeroValues', logical_value, found)
        if (found) setting%ZeroValue%UseZeroValues = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'ZeroValue.UseZeroValues not found'
        
        !%                       Area
        call json%get('ZeroValue.Area', real_value, found)
        if (found) setting%ZeroValue%Area = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'ZeroValue.Area not found'
        
        !%                       Depth
        call json%get('ZeroValue.Depth', real_value, found)
        if (found) setting%ZeroValue%Depth = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'ZeroValue.Depth not found'
        
        !%                       Topwidth
        call json%get('ZeroValue.Topwidth', real_value, found)
        if (found) setting%ZeroValue%Topwidth = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'ZeroValue.Topwidth not found'

        !%                       Volume
        call json%get('ZeroValue.Volume', real_value, found)
        if (found) setting%ZeroValue%Volume = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'ZeroValue.Volume not found'

    !% Debug.File =====================================================================
        !%                       
        call json%get('Debug.File.adjust', logical_value, found)
        if (found) setting%Debug%File%adjust = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.adjust not found'
        
        !%                       
        call json%get('Debug.File.BIPquick', logical_value, found)
        if (found) setting%Debug%File%BIPquick = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.BIPquick not found'
        
        !%                       
        call json%get('Debug.File.boundary_conditions', logical_value, found)
        if (found) setting%Debug%File%boundary_conditions = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.boundary_conditions not found'
        
        !%                       
        call json%get('Debug.File.c_library', logical_value, found)
        if (found) setting%Debug%File%c_library = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.c_library not found'
        
        !%                       
        call json%get('Debug.File.define_globals', logical_value, found)
        if (found) setting%Debug%File%define_globals = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.define_globals not found'
        
        !%                       
        call json%get('Debug.File.define_indexes', logical_value, found)
        if (found) setting%Debug%File%define_indexes = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.define_indexes not found'
        
        !%                       
        call json%get('Debug.File.define_keys', logical_value, found)
        if (found) setting%Debug%File%define_keys = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.define_keys not found'
        
        !%                       
        call json%get('Debug.File.define_settings', logical_value, found)
        if (found) setting%Debug%File%define_settings = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.define_settings not found'
        
        !%                       
        call json%get('Debug.File.define_types', logical_value, found)
        if (found) setting%Debug%File%define_types = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.define_types not found'
        
        !%                       
        call json%get('Debug.File.diagnostic_elements', logical_value, found)
        if (found) setting%Debug%File%diagnostic_elements = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.diagnostic_elements not found'
        
        !%                       
        call json%get('Debug.File.discretization', logical_value, found)
        if (found) setting%Debug%File%discretization = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.discretization not found'
        
        !%                       
        call json%get('Debug.File.face', logical_value, found)
        if (found) setting%Debug%File%face = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.face not found' 
                
        !%                       
        call json%get('Debug.File.finalization', logical_value, found)
        if (found) setting%Debug%File%finalization = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.finalization not found'
        
        !%                       
        call json%get('Debug.File.geometry', logical_value, found)
        if (found) setting%Debug%File%geometry = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.geometry not found'
        
        !%                       
        call json%get('Debug.File.initial_condition', logical_value, found)
        if (found) setting%Debug%File%initial_condition = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.initial_condition not found'
        
        !%                       
        call json%get('Debug.File.initialization', logical_value, found)
        if (found) setting%Debug%File%initialization = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.initialization not found'
        
        !%                       
        call json%get('Debug.File.interface', logical_value, found)
        if (found) setting%Debug%File%interface = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.interface not found'
        
        !%                       
        call json%get('Debug.File.jump', logical_value, found)
        if (found) setting%Debug%File%jump = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.jump not found'
        
        !%                       
        call json%get('Debug.File.lowlevel_rk2', logical_value, found)
        if (found) setting%Debug%File%lowlevel_rk2 = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.lowlevel_rk2 not found'
        
        !%                       
        call json%get('Debug.File.network_define', logical_value, found)
        if (found) setting%Debug%File%network_define = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.network_define not found'
        
        !%                       
        call json%get('Debug.File.orifice_elements', logical_value, found)
        if (found) setting%Debug%File%orifice_elements = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.orifice_elements not found'
        
        !%                       
        call json%get('Debug.File.output', logical_value, found)
        if (found) setting%Debug%File%output = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.output not found'
               
        !%                       
        call json%get('Debug.File.pack_mask_arrays', logical_value, found)
        if (found) setting%Debug%File%pack_mask_arrays = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.pack_mask_arrays not found'
        
        !%                       
        call json%get('Debug.File.partitioning', logical_value, found)
        if (found) setting%Debug%File%partitioning = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.partitioning not found'
        
        !%                       
        call json%get('Debug.File.pump_elements', logical_value, found)
        if (found) setting%Debug%File%pump_elements = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.pump_elements not found'
        
        !%                       
        call json%get('Debug.File.rectangular_channel', logical_value, found)
        if (found) setting%Debug%File%rectangular_channel = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.rectangular_channel not found'
        
        !%                       
        call json%get('Debug.File.trapezoidal_channel', logical_value, found)
        if (found) setting%Debug%File%trapezoidal_channel = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.trapezoidal_channel not found'
        
        !%                       
        call json%get('Debug.File.runge_kutta2', logical_value, found)
        if (found) setting%Debug%File%runge_kutta2 = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.runge_kutta2 not found'
        
        !%                       
        call json%get('Debug.File.timeloop', logical_value, found)
        if (found) setting%Debug%File%timeloop = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.timeloop not found'
        
        !%                       
        call json%get('Debug.File.update', logical_value, found)
        if (found) setting%Debug%File%update = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.update not found'
 
        !%                       
        call json%get('Debug.File.utility', logical_value, found)
        if (found) setting%Debug%File%utility = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.utility not found'
                
        !%                       
        call json%get('Debug.File.utility_allocate', logical_value, found)
        if (found) setting%Debug%File%utility_allocate = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.utility_allocate not found'
        
        !%                       
        call json%get('Debug.File.utility_deallocate', logical_value, found)
        if (found) setting%Debug%File%utility_deallocate = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.utility_deallocate not found'
        
        !%                       
        call json%get('Debug.File.utility_array', logical_value, found)
        if (found) setting%Debug%File%utility_array = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.utility_array not found'
        
        !%                       
        call json%get('Debug.File.utility_datetime', logical_value, found)
        if (found) setting%Debug%File%utility_datetime = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.utility_datetime not found'
        
        !%                       
        call json%get('Debug.File.utility_interpolate', logical_value, found)
        if (found) setting%Debug%File%utility_interpolate = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.utility_interpolate not found'
        
        !%                       
        call json%get('Debug.File.utility_output', logical_value, found)
        if (found) setting%Debug%File%utility_output = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.utility_output not found'
        
        !%                       
        call json%get('Debug.File.utility_string', logical_value, found)
        if (found) setting%Debug%File%utility_string = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.utility_string not found'
        
        !%                       
        call json%get('Debug.File.weir_elements', logical_value, found)
        if (found) setting%Debug%File%weir_elements = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.File.weir_elements not found'
        

    !% Debug.FileGroup =====================================================================
        !%                       
        call json%get('Debug.FileGroup.all', logical_value, found)
        if (found) setting%Debug%FileGroup%all = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.FileGroup.all not found'
        
        !%                       
        call json%get('Debug.FileGroup.definitions', logical_value, found)
        if (found) setting%Debug%FileGroup%definitions = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.FileGroup.definitions not found'
        
        !%                       
        call json%get('Debug.FileGroup.finalization', logical_value, found)
        if (found) setting%Debug%FileGroup%finalization = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.FileGroup.finalization not found'
        
        !%                       
        call json%get('Debug.FileGroup.geometry', logical_value, found)
        if (found) setting%Debug%FileGroup%geometry = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.FileGroup.geometry not found'
        
        !%                       
        call json%get('Debug.FileGroup.initialization', logical_value, found)
        if (found) setting%Debug%FileGroup%initialization = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.FileGroup.initialization not found'
        
        !%                       
        call json%get('Debug.FileGroup.interface', logical_value, found)
        if (found) setting%Debug%FileGroup%interface = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.FileGroup.interface not found'
        
        !%                       
        call json%get('Debug.FileGroup.output', logical_value, found)
        if (found) setting%Debug%FileGroup%output = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.FileGroup.output not found'
        
        !%                       
        call json%get('Debug.FileGroup.timeloop', logical_value, found)
        if (found) setting%Debug%FileGroup%timeloop = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.FileGroup.timeloop not found'
        
        !%                       
        call json%get('Debug.FileGroup.utility', logical_value, found)
        if (found) setting%Debug%FileGroup%utility = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.FileGroup.utility not found'
        
    !% Debug. =====================================================================
        !%                       'Debug.SetupYN
        call json%get('Debug.SetupYN', logical_value, found)
        if (found) setting%Debug%SetupYN = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.Setup not found'
        
        !%                       Debug.OutputYN
        call json%get('Debug.OutputYN', logical_value, found)
        if (found) setting%Debug%OutputYN = logical_value
        !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.Output not found'

        call define_settings_update_debug()
 
     !% finished
        call json%destroy()
        if (json%failed()) stop "JSON failed to destroy"

        if (setting%Debug%File%define_settings) &
            write(*,"(A,i5,A)") '*** leave ' // trim(subroutine_name) // " [Processor ", this_image(), "]"
    end subroutine define_settings_load
!%
!%==========================================================================
!%==========================================================================
!%
    subroutine define_settings_update_debug()
        !% -----------------------------------------------------------------
        !% Description:
        !% Uses the Debug%FileGroup... = true option to turn on the debug
        !% for all the files in that group
        !% -----------------------------------------------------------------
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
    end subroutine define_settings_update_debug
!%
!%==========================================================================
!%==========================================================================
!%
end module define_settings
