module define_settings
    !%==========================================================================
    !% SWMM5+ release, version 1.0.0
    !% 20230608
    !% Hydraulics engine that links with EPA SWMM-C
    !% June 8, 2023
    !%
    !% Description:
    !% Defines types used in the setting%...%... structure
    !%
    !% Methods
    !% Settings are stored in the setting object which is public.
    !% Settings can be defined in three ways:
    !%    1) via flags at command line
    !%    2) via JSON file
    !%    3) via default values
    !% If the settings are not defined via JSON file, default values
    !% should be defined in the declaration of the derived type
    !% associated to the setting. Settings flags (e.g., verbose -v)
    !% overwrite any value defined via default or JSON values.
    !%
    !% To understand the setting structure, programmers should first
    !% familiarize themselves with the First Level Types at the bottom 
    !% of this file
    !%==========================================================================

    use json_module
    use define_keys
    use define_globals
    use utility_string, only: util_lower_case

    implicit none
    public

!%==========================================================================
!% PUBLIC TYPES
!%==========================================================================
!%
    !%===================================================================
    !% Third Level Types
    !%===================================================================

    !% setting%Adjust%Flowrate
    type AdjustFlowrateType
        logical :: ApplyYN = .true.
        integer :: Approach = vshape !% vshape, doesnotexist
        real(8) :: Coef = 1.0d0
        real(8) :: SmallDepthMultiplier = 3.0d0
    end type AdjustFlowrateType

    !% setting%Adjust%Head
    type AdjustHeadType
        logical :: ApplyYN = .false.  !% BUG IN CODE FOR ADJUSTMENT 
        !% --- only use vshape_surcharge_CC.  
        integer :: Approach = vshape_surcharge !% vshape_surcharge, doesnotexist
        real(8) :: Coef = 1.0d0
    end type AdjustHeadType

    !% setting%Output%CommandLine
    type CommandLineType
        logical :: quietYN = .false.
        integer :: interval = 1000 !% steps between commandline output
    end type CommandLineType

    !% setting%Time%CPU
    !% NOT USER SETTINGS
    type CPUTimeType
        real (8) :: EpochStartSeconds = 0.0d0
        real (8) :: EpochNowSeconds  = 0.0d0
        real (8) :: EpochFinishSeconds  = 0.0d0
    end type CPUTimeType

    !% setting%Output%DataOut
    type DataOutType
        logical :: isAreaOut                = .true.
        logical :: isDepthOut               = .true.
        logical :: isFlowrateOut            = .true.
        logical :: isFlowrateAvgOut         = .false. 
        logical :: isFluxConsOut            = .false.
        logical :: isFroudeNumberOut        = .false.
        logical :: isHeadOut                = .true.
        logical :: isHydRadiusOut           = .false.
        logical :: isPerimeterOut           = .false.
        logical :: isManningsNout           = .false.
        logical :: isSlotWidthOut           = .false.
        logical :: isSlotDepthOut           = .false.
        logical :: isTopWidthOut            = .false.
        logical :: isVelocityOut            = .true.
        logical :: isVolumeOut              = .true.
        logical :: isVolumeConsOut          = .false.
        logical :: isVolumeOverflowOut      = .true.
        logical :: isVolumePondedOut        = .true.
        logical :: isWaveSpeedOut           = .false.
        logical :: isPreissmannCelerityOut  = .false.
        logical :: isPreissmannNumberOut    = .false.
        Logical :: isElemLengthOut          = .true.  
        Logical :: isElemBottomSlopeOut     = .false.      
        Logical :: isElemBreathMaxOut       = .false.      
        Logical :: isElemFullAreaOut        = .false.      
        Logical :: isElemFullDepthOut       = .false.  
        Logical :: isElemManningsOut        = .false.      
        Logical :: isElemZBottomOut         = .true.  
        Logical :: isElemZCrownOut          = .true.  
        Logical :: isLinkLengthOut          = .false.  
        Logical :: isLinkAdjustedLengthOut  = .false.          
        Logical :: isLinkInletOffsetOut     = .false.      
        Logical :: isLinkOutletOffsetOut    = .false.          
        Logical :: isLinkSlopeOut           = .false.  
        Logical :: isLinkZBottomUpOut       = .false.      
        Logical :: isLinkZBottomDownOut     = .false.      
        Logical :: isNodeZBottomOut         = .false.  
        Logical :: isNodeFullDepthOut       = .false.  
    end type DataOutType

    !% setting%Solver%ForceMain
    type ForceMainType
        logical :: AllowForceMainTF      = .true.  !% when false, FM is forced to standard pipe with default ManningsN
        logical :: UseSWMMinputMethodTF  = .true.
        logical :: FMallClosedConduitsTF = .false. 
        logical :: errorCheck_RoughnessTF     = .true. !% check scale of FM coef is consistent with HW or DW
        logical :: HazenWilliams_equivalent_isSetTF = .false.  !% NOT A USER SETTING
        integer :: Default_method = HazenWilliams
        real(8) :: Default_HazenWilliams_coef = 120          !% 120  is concrete
        real(8) :: Default_DarcyWeisbach_roughness_mm = 0.36 !% 0.36 is concrete, good joints 
        real(8) :: Default_ManningsN = 0.03  !% used on FM elements when AllowForceMainTF = false
        real(8) :: minimum_slope  = 1.0d-3  !% minimum slope in HW computation
    end type ForceMainType

    !% setting%Limiter%Dt
    type LimiterDtType
        logical :: FailOnMinYN = .true.
        logical :: FailOnMaxYN = .false.
        logical :: UseLimitMinYN = .true.
        real(8) :: Minimum     = 1.0d-3
        logical :: UseLimitMaxYN = .true.
        real(8) :: Maximum     = 86400.0d0
    end type LimiterDtType

    !% setting%Limiter%InterpWeight
    type LimiterInterpWeightType
        real(8) :: Maximum = 1.0d8
        real(8) :: Minimum = 1.0d-8
    end type LimiterInterpWeightType

    !% setting%Limiter%Velocity
    type LimiterVelocityType
        logical :: UseLimitMaxYN = .true.
        real(8) :: Maximum = 10.0d0 ! m/s
        logical :: ZeroMinimumVelocitiesYN = .true. 
        real(8) :: Minimum = 1.0d-8 ! m/s below this set to zero
    end type LimiterVelocityType

    !% setting%Limiter%Flowrate
    type LimiterFlowrateType
        logical :: UseSWMMlinkValueYN = .false. !% --- SWMM input file hard limit (do not use)
        logical :: UseLocalVolumeYN   = .true.  !% --- limit outflow volumes by local volume factor
        real(8) :: LocalVolumeFactor  = 0.49d0   !% --- fraction of local volume plus lateral inflow that can exit
    end type LimiterFlowrateType

    !% setting%Solver%ManningsN
    type ManningsNtype
        real(8) :: alpha = 1.0d0
        real(8) :: beta  = 1.0d0
        real(8) :: FlowReversalFactor = 100.d0 !% multiplier for Manning's n on flow reversal or small velocity
        real(8) :: SmallVelocity      = 1.0d-3 !% small velocity below which Manning's n in increased.
    end type ManningsNtype

    !% setting%Junction%Overflow
    type OverflowType
        !% --- global values are initialized in define_settings_default()
        !%     changes made here have no effect!
        real(8) :: OrificeLength !% length of overflow orifice (m)
        real(8) :: OrificeHeight!% height of overflow orifice (m)
        real(8) :: WeirLengthFactor !% multiplier of the default overflow weir length that is derived from junction plan area
        real(8) :: CbroadCrestedWeir
        real(8) :: coef1 
        real(8) :: coef2
        real(8) :: coef3
        real(8) :: coef4
    end type

    !% setting%Junction%PlanArea
    type PlanAreaType
        logical :: UseLargeBranchStorageTF = .true.  !% estimated plan area for junctions without defined plan area (false causes blowup for large diameter conduits)
        real(8) :: AreaMinimum               !% NOT A USER SETTING !% the surfaceArea_minimum from EPA-SWMM, minimum for ImpliedStorage plan area
        real(8) :: AreaFactorMaximum = 100.d0 !% maximum multiplier by which branch-implied plan are can be greater than surfaceAreaMinimum
        real(8) :: LargeBranchDepth  = 1.0d0 !% Depth (m) used for breadth in plan area computation fo large branch width.
                                             !% Note that if the LargeBranchDepth is too large, then trapezodial and irregular channels will have
                                             !% large topwidths that will cause large plan areas at their junctions.
    end type

    !% setting%Solver%PreissmannSlot
    type PreissmannSlotType
        logical :: useSlotTF = .true.
        !% Allowable values: DynamicSlot, StaticSlot, SplitDynamicSlot
        integer :: Method = DynamicSlot
        real(8) :: TargetCelerity = 15.0d0
        real(8) :: Alpha = 3.0d0
        real(8) :: DecayRate = 5.0
        real(8) :: initPNminimum = 3.0d0
    end type PreissmannSlotType

    !% setting%Output%Report
    type ReportType
        logical :: useSWMMinpYN = .true.
        logical :: provideYN = .true.
        logical :: useHDF5   = .true.
        logical :: useCSV    = .true.
        logical :: suppress_MultiLevel_Output = .false.
        real(8) :: StartTime = 0.0d0
        real(8) :: TimeInterval = 300.0d0
        integer :: ThisStep                  !% NOT A USER SETTING
        integer :: TimeUnits = InHours
    end type ReportType

    !% setting%Time% ...Hydraulics, Hydrology, Dry
    !% these is initialized in define_settings_defaults
    type TimeStepType
       real(8) :: Dt
       real(8) :: LastTime  !% NOT A USER SETTING
       real(8) :: NextTime  !% NOT A USER SETTING
       integer(kind=8) :: Step  !% NOT A USER SETTING
    end type TimeStepType

    !% setting%File%UnitNumber
    !% NOT USER SETTINGS
    type UnitNumberType
        integer :: inp_file
        integer :: out_file
        integer :: rpt_file
        integer :: setting_file
        integer :: outputML_filename_file
        integer :: outputML_control_file
    end type

    !% setting%Time%WallClock
    !% NOT USER SETTINGS
    type WallClockType
        integer         :: StoreInterval = 1000
        integer         :: LastStepStored
        integer (kind=8):: Start = 0
        integer (kind=8):: Now  = 0
        integer (kind=8):: End = 0
        integer (kind=8):: LastTimeStored
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
        integer (kind=8):: SharedStart = 0
        integer (kind=8):: SharedStop = 0
        integer (kind=8):: SharedCumulative = 0
        integer (kind=8):: TimeMarchStart = 0
        integer (kind=8):: TimeMarchEnd = 0
        integer (kind=8):: InitializationEnd = 0
        integer (kind=8):: PartitionEnd = 0
        integer (kind=8):: FinalOutputStart = 0
        integer (kind=8):: SharedStart_A = 0
        integer (kind=8):: SharedStop_A = 0
        integer (kind=8):: SharedCumulative_A = 0
        integer (kind=8):: SharedStart_B = 0
        integer (kind=8):: SharedStop_B = 0
        integer (kind=8):: SharedCumulative_B = 0
        integer (kind=8):: SharedStart_C = 0
        integer (kind=8):: SharedStop_C = 0
        integer (kind=8):: SharedCumulative_C = 0
    end type

    !% setting%Weir% ...Transverse, SideFlow, VNotch, Trapezoidal
    !%      default values provided in define_settings_default
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

    !%===================================================================
    ! Second Level Types
    !%===================================================================

    ! setting%Adjust
    type AdjustType
        type(AdjustFlowrateType)   :: Flowrate
        type(AdjustHeadType)       :: Head
        real(8), dimension(12)     :: Temperature   !% from SWMM input [ADJUSTMENTS] 
        real(8), dimension(12)     :: Evaporation   !% from SWMM input [ADJUSTMENTS]
        real(8), dimension(12)     :: Rainfall      !% from SWMM input [ADJUSTMENTS] 
        real(8), dimension(12)     :: Conductivity  !% from SWMM input [ADJUSTMENTS] 
    end type AdjustType

    ! setting%AirTracking
    type AirTrackingType
        logical :: UseAirTrackingYN = .false.
    end type AirTrackingType

    ! setting%BC
    type BCPropertiesType
        integer :: TimeSlotsStored        = 1000
        logical :: disableInterpolationYN = .false.
        real(8) :: smallestTimeInterval   = 86400.d0  !% NOT A USER SETTING
    end type BCPropertiesType

    ! setting%CaseName
    type CaseNameType
        character(256) ::    Long = 'default_CaseName'
        character(16)  ::    Short = 'default'
        character(31)  ::    withTimeStamp   !% NOT A USER SETTING
    end type CaseNameType

    !% setting%Climate
    type ClimateType
        logical :: useHydraulicsEvaporationTF = .true. !% user may turn off hydrauilcs evaporation
        real(8) :: HydraulicsOnlyIntervalHours = 1.d0 !% hours between updating climate if no runoff
        real(8) :: LastTimeUpdate = 0.d0 !% time (seconds) of last update NOT A USER SETTING
        real(8) :: NextTimeUpdate = 0.d0 !% time (seconds) of next update NOT A USER SETTING
        real(8) :: EvapRate = 0.d0 !% current evaporation rate in m/s
    end type ClimateType

    ! setting%Constant
    type ConstantType
        real(8) :: gravity = 9.81d0 ! m^2/s
        real(8) :: pi = 4.d0*datan(1.d0)   !% NOT A USER SETTING
        real(8) :: water_temperature = 20.d0
        real(8) :: water_kinematic_viscosity = 1.0d-6  !% NOT A USER SETTING
    end type ConstantType

    ! setting%Control
    !% structure to setup a simple time-setting based controls for links
    !% NOT USER SETTINGS
    type ControlType !% NOT A USER SETTING
        integer              :: NumControl = 0          !% number of links controlled  NOT A USER SETTING
        integer, allocatable :: LinkIdx(:)              !% index of the links controlled  NOT A USER SETTING            
        integer, allocatable :: ElemIdx(:)              !% index of the elements controlled  NOT A USER SETTING
        real(8), allocatable :: TimeArray(:,:)          !% arrays of time to change setting. Each column is for a each of the links NOT A USER SETTING
        real(8), allocatable :: SettingsArray(:,:)      !% arrays of settings. Each column is for a each of the links NOT A USER SETTING
        character(len=:), allocatable :: Links(:)       !% link names from SWMM5 input file to be controlled NOT A USER SETTING
    end type ControlType

    !% setting%Crash !% NOT USER SETTING(see util_crash_initialize)
    type CrashType
        real(8) :: DepthMax  !% maximum depth exceedence in 1 cell that causes crash
        real(8) :: HeadMax   !% maximum head exceedence in 1 cell that causes crash
        real(8) :: FlowrateMax   !% maximum flow exceedence in 1 cell that causes crash
        real(8) :: PercentVelocityAtLimit  !% percentage of cells at velocity limit that causes crash
    end type CrashType

    type CulvertType
        logical :: UseCulvertsTF = .true.   !% false allows culverts to be treated as conduits
    end type CulvertType

    !% setting%Debug
    type DebugType
        logical :: isGlobalVolumeBalance = .true. !% conducts global volume balance at each time step
        real(8) :: GlobalVolumeBalance = 0.d0 !% NOT A USER SETTING
        logical :: checkIsNanTF = .false. !% only need true if util_utext_checkIsNan is called (not standard)
        !% THESE debugFile WILL BE OBSOLETE
        type(DebugFileYNType) :: File
        !type(DebugFileGroupYNType) :: FileGroup
        !logical :: SetupYN = .false.
        !logical :: OutputYN = .false.
    end type DebugType

    !% setting%Discretization
    type DiscretizationType
        !logical :: StopOnLengthAdjustTF = .false.  !% Can be used to force code to stop if link lengths are adjusted.
        !% NOTE Channel overflow not tested and is disabled as of 20230508
        logical :: AllowChannelOverflowTF = .false. !% if true, then open channels (CC) can overflow (lose water) NOT IN EPA SWMM
        !logical :: AdjustLinkLengthForJunctionBranchYN = .false.          !% OBSOLETE DO NOT USE TRUE -- if true then JB (junction branch) length is subtracted from link length
        !real(8) :: JunctionBranchLengthFactor  = 1.d0    !% MUST USE 1.0   !% fraction of NominalElemLength used for JB
        ! real(8) :: MinElemLengthFactor = 0.5d0           !% define the minimum allowable fraction of an element size to help with the cfl
        !integer :: MinElemLengthMethod = ElemLengthAdjust
        logical :: UseNominalElemLength = .true.
        logical :: DistributeOpenChannelInflowsTF = .false.
        real(8) :: NominalElemLength   = 10.0d0
        integer :: MinElementPerLink   = 3               !% force a minimum number of elements per link
        logical :: UseEquivalentOrifice = .false.        !% replace small conduits with equivalent orifice
        real(8) :: EquivalentOrificeDischargeCoeff = 0.6 !% discharge coefficient of the equivalent orifice
        real(8) :: MinLinkLength        = 10.0d0         !% below the min link length, replace with equivalent orifice
        real(8) :: FullConduitTopwidthDepthFraction = 0.95d0  !% fraction of full depth used for full topwidth
    end type DiscretizationType

    ! setting%Eps
    type EpsilonType
        !% epsilon relative to machine precision
        real(8) :: Machine = 0.d0  !% computed in code
        !% +- small non-dimensional range for hyd jump discrimination
        real(8) :: FroudeJump = 0.1d0
        !% Fractional increase in depth under froude limitation
        !rm 20220207brh real(8) :: InflowDepthIncreaseFroudeLimit = 0.1
        !% small velocity
        real(8) :: Velocity = 1.0d-6
        !% small head difference (tolerance for no flow)
        real(8) :: Head = 1.0d-6
        !% small time (prevent division by zero)
        real(8) :: TimeStep = 1.0d-6
    end type EpsilonType

    !% setting%File
    type FileType
        logical              :: UseCommandLineFoldersYN  = .true.
        logical              :: force_folder_creationYN = .true.
        logical              :: duplicate_input_file = .true.  !% NOT A USER SETTING, should always be true
        !% standard files and folders
        character(len=256)   :: base_folder = ""
        !character(len=256)   :: library_folder = "/usr/local/lib"  !% Use this for standard release
        character(len=256)   :: library_folder = "build"            !% Use this for developer
        character(len=256)   :: output_folder= "" !
        character(len=256)   :: output_timestamp_subfolder = ""
        character(len=256)   :: output_temp_subfolder = "temp"
        character(len=256)   :: project_folder = "" ! project path
        character(len=256)   :: setting_file = "" ! path to settings JSON file
        character(len=256)   :: input_kernel = "" ! main part of input file name   NOT A USER SETTING
        character(len=256)   :: output_kernel= "" ! main part ouf output file name  NOT A USER SETTING
        character(len=256)   :: inp_file = "" ! path to SWMM input (.inp) file
        character(len=256)   :: rpt_file = "" ! path to SWMM report (.rpt) file  NOT A USER SETTING
        character(len=256)   :: out_file = "" ! path to SWMM output (.out) file  NOT A USER SETTING
        
        !% for multi-level output NOT USER SETTINGS
        character(len=256)   :: outputML_Link_kernel = "link"
        character(len=256)   :: outputML_Node_kernel = "node"
        character(len=256)   :: outputML_combinedfile_kernel = "combined_output" !% filenames that combine across images
        character(len=256)   :: outputML_filename_file = 'output_filenames.txt'  !% list of interim filenames used in combined output
        character(len=256)   :: outputML_control_file = 'control.unf'  !% global parameters that control the processing of combined output
        integer              :: outputML_Ncombined_file_written = 0
        integer              :: outputML_total_timelevels_written = 0

        !% file handling !% NOT USER SETTINGS
        integer              :: last_unit = 1000 !% starting point for assigning unit numbers
        type(UnitNumberType) :: UnitNumber  
    end type FileType

    ! setting%Junction
    type JunctionType
        !% 20230405 brh  Default is FALSE to force JM, and TRUE to Force Storage 
        logical :: ForceNodesJM = .false.  !% forces CC nodes between two conduits to be nJM rather than nJ2 faces
        !%                                 !% note CC nodes will only be nJ2 faces if SurchargeDepth = InfinitExtraDepthValue
                                           !% ONLY SET TRUE FOR ALGORITHM TESTING
        !% NOTE ForceStorage must be true as of 20230507. Future extension may include junction solution that does
        !% not require the minimum surface area of the ImpliedStorage type
        logical :: ForceStorage = .true.        !% forces nJM junctions without explicit storage to have implied storage
        integer :: FunStorageN  = 10            !% number of curve entries for functional storage   
        real(8) :: kFactor      = 0.5d0         !% default entrance/exit losses at junction branch (use 0.5 )
        real(8) :: OutflowDampingFactor = 0.99d0 !% 0 to <1, with 0 being strictly energy-based outflow, 1 being strictly extrapolated from neighbor (fails); 0.99 recommended
        real(8) :: ZeroHeadDiffValue = 1.0d-8   !% Head difference (m) that results in zero outflow in a branch. 
        real(8) :: ZeroOutflowValue  = 1.0d-8   !% Outflow (m^3/s) that results in zero outflow in a branch
        real(8) :: InfiniteExtraDepthValue = 999.d0  !% Surcharge Depth if this value or higher is treated as impossible to overflow

        !% Ponding ScaleFactor is multiplier of junction/storage length scale (sqrt of area) to get minimum length scale of ponding
        real(8) :: PondingScaleFactor = 10.d0 
        type(OverflowType) :: Overflow
        type(PlanAreaType) :: PlanArea
    end type JunctionType

    ! setting%Limiter
    type LimiterType
        real(8) :: NormalDepthInfinite   = 1000.d0     ! value used when normal depth would be infinite.
        type(LimiterInterpWeightType) :: InterpWeight
        type(LimiterVelocityType)     :: Velocity
        type(LimiterFlowrateType)     :: Flowrate
        type(LimiterDtType)           :: Dt
        !% --- limiter%VolumeFractionInTimeStep: limits momentum/velocity in RK solution to a rate consistent with the
        !%     VolumeFractionInTimeStep * volume / dt. This helps prevent oscillations in very small depths.
        !%     Should generally be less than one.  1/4 seems to work well.
        real(8)                       :: VolumeFractionInTimeStep = 0.5d0  
    end type LimiterType

    type LinkType
        ! Allowable: UniformDepth, LinearlyVaryingDepth, IncreasingDepth,  FixedHead
        integer :: DefaultInitDepthType = LinearlyVaryingDepth 
        !% --- if users include an unnecessarily large full depth, this can affect the precision in 
        !%     look up tables. The following allows global replacement of excessive maximum depths.
        logical :: OpenChannelLimitDepthYN = .false.  !% Y/N overwrite for the max depth from SWMM.inp file
        real(8) :: OpenChannelFullDepth = 100.d0      !% overwrite the max depth of open channels in SWMM.inp file
    end type

    !% setting%Orifice
    type OrificeType
        logical :: ForceOrificeNodesToJM = .false.  !% true only used for algorithm design
        real(8) :: SharpCrestedWeirCoefficient = 0.414
        real(8) :: TransverseWeirExponent = 1.5
        real(8) :: VillemonteCorrectionExponent = 0.385
    end type OrificeType

    !% setting%Output
    type OutputType
        logical :: UseFileNameFile = .false.
        logical, allocatable :: ElementsExist_byImage(:)[:]  !% NOT A USER SETTING
        logical, allocatable :: FacesExist_byImage(:)[:]     !% NOT A USER SETTING
        logical :: ElementsExist_global                      !% NOT A USER SETTING 
        logical :: FacesExist_global                         !% NOT A USER SETTING  
        logical :: BarrelsExist = .false.                    !% NOT A USER SETTING
        logical :: Verbose = .true.
        logical :: Warning = .true.
        integer :: LastLevel = 0                             !% NOT A USER SETTING
        integer :: MaxExpectedLevels = 0                     !% NOT A USER SETTING
        integer :: NumberOfWriteSteps = 0                    !% NOT A USER SETTING
        integer :: NumberOfTimeLevelSaved = 0                !% NOT A USER SETTING
        integer :: StoredLevels = 100        
        integer :: StoredFileNames = 100
        integer :: ElemHeadIndex = 0                         !% NOT A USER SETTING
        integer :: FaceUpHeadIndex = 0                       !% NOT A USER SETTING
        integer :: faceDnHeadIndex = 0                       !% NOT A USER SETTING
        integer :: ElemFlowAvgIndex= 0                       !% NOT A USER SETTING 
        integer (kind=8) :: MemoryStorageMax = 29000000  
        type(CommandLineType) :: CommandLine
        type(DataOutType) :: DataOut
        type(ReportType) :: Report
    end type OutputType

    !% setting%Partitioning
    type PartitioningType
        !% Allowable values: BQquick
        integer :: PartitioningMethod = BQuick
        logical :: PhantomLinkAdjust  = .true.
    endtype PartitioningType

    !% setting%Profile
    type ProfileType
        logical :: useYN = .false.
        !logical :: Tests = .false.
        !type(ProfileFileYNType) :: File
        !type(ProfileFileGroupYNType) :: FileGroup
        !logical :: Input
        !logical :: Output
    end type ProfileType

    !% setting%Pump 
    type PumpSettingType 
        real(8) :: RampupTime = 5.d0  !% seconds for pump ramp up. If too small, instabilities can result (0 causes seg fault)
        real(8) :: MinShutoffTime = 60.d0  !% seconds for pump to be idle before starting again 
    end type PumpSettingType

    !% setting%Simulation
    type SimulationType
        logical :: AllowReverseGradientInitialConditionsTF = .false.
        logical :: stopAfterInitializationYN = .false.
        logical :: useHydrology    = .true.
        logical :: useHydraulics   = .true.
        logical :: useSpinUp       = .false.
        logical :: stopAfterSpinUp = .false.
        real(8) :: SpinUpDays      = 10.d0
    end type SimulationType

    ! setting%SmallDepth
    type SmallDepthType
        logical :: useMomentumCutoffYN = .false. !% causes reversed flow with JB 20230410brh
        real(8) :: MomentumDepthCutoff = 0.01d0 ! m  !% prior to 20230327 using 0.03
        real(8) :: ManningsN = 0.1d0
        real(8) :: PumpVolumeFactor = 0.5d0  !% 
        real(8) :: LateralInflowSmallDepth = 0.05d0
        real(8) :: BCInflowSmallDepth = 0.05d0
        real(8) :: InflowSmallDepthMultiplier = 2.d0
    end type SmallDepthType

    ! setting%Solver
    type SolverType
        !logical :: PreissmannSlot = .true.
        logical :: SubtractReferenceHead = .false.
        integer :: MomentumSourceMethod = T10 
        integer :: SolverSelect = ETM
        real(8) :: SwitchFractionDn = 0.8d0
        real(8) :: SwitchFractionUp = 0.9d0
        real(8) :: ReferenceHead = zeroR
        real(8) :: AverageZbottom = zeroR               !% NOT A USER SETTING
        real(8) :: MaxZbottom = zeroR                   !% NOT A USER SETTING
        real(8) :: MinZbottom = zeroR                   !% NOT A USER SETTING
        real(8), dimension(2)    :: crk2 = [0.5d0, 1.0d0]  !% NOT A USER SETTING
        type(ForceMainType)      :: ForceMain
        type(ManningsNtype)      :: ManningsN
        type(PreissmannSlotType) :: PreissmannSlot
    end type SolverType


    !% storage for data read from SWMM input file
    !% NOT USER SETTINGS
    !% setting%SWMMinput...
    type SWMMinputType
        !% --- following SWMM *.inp file Rossman user manual
        integer :: FlowUnitsType
        integer :: InfiltrationType
        integer :: FlowRoutingType
        integer :: LinkOffsetsType
        integer :: ForceMainEquationType
        logical :: IgnoreRainfall
        logical :: IgnoreSnowmelt
        logical :: IgnoreGroundwater
        logical :: IgnoreRDII
        logical :: IgnoreRouting
        logical :: IgnoreQuality
        logical :: AllowPonding
        logical :: SteadyState_Skip
        real(8) :: SteadyState_System_FlowrateTolerance
        real(8) :: SteadyState_Lateral_FlowrateTolerance
        real(8) :: StartEpoch           = nullvalueR
        real(8) :: EndEpoch             = nullvalueR
        real(8) :: ReportStartTimeEpoch = nullvalueR
        integer :: Sweep_Start_Day = 1
        integer :: Sweep_End_Day   = 365
        integer :: DryDaysBeforeStart = 0
        real(8) :: ReportTimeInterval      = nullvalueR
        real(8) :: Hydrology_WetStep       = nullvalueR
        real(8) :: Hydrology_DryStep       = nullvalueR
        real(8) :: Hydraulics_RouteStep    = nullvalueR
        real(8) :: RoutingStep_LengtheningTime
        real(8) :: RoutingStep_CourantFactor
        real(8) :: RoutingStep_Minimum
        integer :: InertialDampingType
        integer :: NormalFlowLimiterType
        real(8) :: SurfaceArea_Minimum
        real(8) :: ConduitSlope_Minimum
        integer :: NumberOfTrials_Maximum
        real(8) :: Head_ConvergenceTolerance
        integer :: NumberParallelThreads
        logical :: TempDirectory_Provided
        !% --- other SWMM configuration values
        integer :: ControlRuleStep = nullvalueI 
        integer :: SurchargeMethod
        real(8) :: TotalDuration = nullvalueR
        !% --- object counts
        integer :: N_gage = zeroI
        integer :: N_subcatch = zeroI
        integer :: N_node = zeroI
        integer :: N_link = zeroI
        integer :: N_pollutant = zeroI
        integer :: N_landuse = zeroI
        integer :: N_timepattern = zeroI
        integer :: N_curve = zeroI
        integer :: N_tseries = zeroI
        integer :: N_control = zeroI
        integer :: N_transect = zeroI
        integer :: N_aquifer = zeroI
        integer :: N_unithyd = zeroI
        integer :: N_snowmelt = zeroI
        integer :: N_shape = zeroI
        integer :: N_LID = zeroI
        integer :: N_junction = zeroI
        integer :: N_outfall  = zeroI
        integer :: N_storage  = zeroI
        integer :: N_divider = zeroI
        !% --- other counts
        integer :: N_groundwater = zeroI
        !integer :: N_link_transect = zeroI
        integer :: N_transect_depth_items = zeroI
    end type SWMMinputType 
    
    !% setting%TestCase
    type TestCaseType  !% NOT WORKING YET
        logical       :: UseTestCaseYN = .false.
        character(64) :: TestName
    end type TestCaseType

    ! setting%Time
    type TimeType
        logical            :: useSWMMinpYN = .true.
        logical            :: matchHydrologyStep = .true.
        real(8)            :: DtTol = 1.0d-1   !% tolerance when comparing Now time and accumulation of time steps.
        character(14)      :: DateTimeStamp  !% NOT A USER SETTING
        integer(kind=8)    :: Step  !% NOT A USER SETTING count of the number of steps (either hydrology, hydraulics, or combined)
        real(8)            :: Start !% NOT A USER SETTING
        real(8)            :: Now  ! % NOT A USER SETTING
        real(8)            :: End   !% NOT A USER SETTING
        real(8)            :: StartEpoch = nullvalueR  !% NOT A USER SETTING
        real(8)            :: EndEpoch   = nullvalueR  !% NOT A USER SETTING
        type(TimeStepType) :: Hydraulics
        type(TimeStepType) :: Hydrology
        type(TimeStepType) :: ControlRule  !% NOT A USER SETTING
        type(WallClockType):: WallClock    !% NOT A USER SETTING
        type(CPUTimeType)  :: CPU          !% NOT A USER SETTING
    end type TimeType

    !% setting%Weir  !% values set in define_settings_default
    type WeirType
        logical :: ForceWeirNodesToJM = .false.  !% true only used for algorithm design
        type(WeirConstantType) :: Transverse
        type(WeirConstantType) :: SideFlow
        type(WeirConstantType) :: VNotch
        type(WeirConstantType) :: Trapezoidal
    end type WeirType

    !% setting%VariableDT
    type VariableDTType
        logical :: ApplyYN = .true.
        logical :: limitByBC_YN = .true.  !% limits time step by BC step for inflows/head
        real(8) :: CFL_hi = 0.6d0
        real(8) :: CFL_target = 0.5d0
        real(8) :: CFL_lo = 0.4d0
        real(8) :: CFL_inflow_max = 0.4d0
        real(8) :: decreaseFactor = 0.8d0  
        real(8) :: increaseFactor = 1.2d0 
        real(8) :: InitialDt = 1.d0  !% recommend 1 second to reduce chance of initial instability (should be less than pump rampup time)
        integer :: NstepsForCheck = 10
        integer(kind=8) :: LastCheckStep = 0  !% NOT A USER SETTING
    end type VariableDTType

    !% setting%ZeroValue
    !% Note that Depth is the setting users should change
    type ZerovalueType
        logical :: UseZeroValues = .true.
        real(8) :: Area = 1.d-12 ! m^2 -- NOT A USER SETTING
        real(8) :: Depth = 1.d-4 ! m
        real(8) :: Slope = 1.e-6 ! prevents zero values (may be + or -)
        real(8) :: Topwidth = 1.d-12 ! m -- NOT A USER SETTING
        real(8) :: Volume = 1.d-12 ! m^3 -- NOT A USER SETTING
        real(8) :: VolumeResetLevel !m^3 -- NOT A USER SETTING
        real(8) :: Velocity = 1.d-3
    end type ZerovalueType

    !%===================================================================
    !% First Level Type (setting)
    !%===================================================================

    type settingType
        logical                  :: JSON_FoundFileYN = .false.
        logical                  :: JSON_CheckAllInputYN = .true.
        type(AdjustType)         :: Adjust
        type(AirTrackingType)    :: AirTracking
        type(BCPropertiesType)   :: BC
        type(CaseNameType)       :: CaseName ! name of case
        type(ClimateType)        :: Climate  ! climate controls
        type(ConstantType)       :: Constant ! Constants
        type(ControlType)        :: Control  ! Control data structure
        type(CrashType)          :: Crash    !% conditions where code is considered crashin
        type(CulvertType)        :: Culvert  
        type(DebugType)          :: Debug
        type(DiscretizationType) :: Discretization
        type(EpsilonType)        :: Eps ! epsilons used to provide bandwidth for comparisons
        type(FileType)           :: File
        type(JunctionType)       :: Junction
        type(LimiterType)        :: Limiter ! maximum and minimum limiters
        type(LinkType)           :: Link
        type(OrificeType)        :: Orifice
        type(OutputType)         :: Output
        type(PartitioningType)   :: Partitioning
        !type(PreissmannSlotType) :: PreissmannSlot
        type(ProfileType)        :: Profile
        type(PumpSettingType)    :: Pump
        type(SimulationType)     :: Simulation
        type(SmallDepthType)     :: SmallDepth! controls for small (non-zero) depths
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

    setting%Time%Hydrology%Dt = 600.0d0

    setting%Weir%Transverse%WeirExponent = 1.5d0
    setting%Weir%Transverse%WeirContractionFactor = 0.1d0
    setting%Weir%Transverse%SideFlowWeirCrestExponent = 1.0d0
    setting%Weir%Transverse%VillemonteCorrectionExponent = 0.385d0

    setting%Weir%SideFlow%WeirExponent = 1.67d0
    setting%Weir%SideFlow%WeirContractionFactor = 0.1d0
    setting%Weir%SideFlow%SideFlowWeirCrestExponent = 0.83d0
    setting%Weir%SideFlow%VillemonteCorrectionExponent = 0.385d0

    setting%Weir%Vnotch%WeirExponent = 2.5d0
    setting%Weir%Vnotch%WeirContractionFactor = 1.0d0
    setting%Weir%Vnotch%SideFlowWeirCrestExponent = 1.0d0
    setting%Weir%Vnotch%VillemonteCorrectionExponent = 0.385d0

    setting%Weir%Trapezoidal%WeirExponent = 1.5d0
    setting%Weir%Trapezoidal%WeirContractionFactor = 1.0d0
    setting%Weir%Trapezoidal%SideFlowWeirCrestExponent = 1.0d0
    setting%Weir%Trapezoidal%VillemonteCorrectionExponent = 0.385d0

    Setting%Orifice%SharpCrestedWeirCoefficient = 0.414d0
    Setting%Orifice%TransverseWeirExponent = 1.5d0
    Setting%Orifice%VillemonteCorrectionExponent = 0.385d0

    !% --- overflow orifice for junctions based on broad-crested weir using
    !%     non-dimensional formulation of Brater and King Eq. 5.10 and Table 5.1
    setting%Junction%Overflow%OrificeLength     = 2.0d0 !% length of overflow orifice (m)
    setting%Junction%Overflow%OrificeHeight     = 0.2d0 !% height of overflow orifice (m)
    setting%Junction%Overflow%WeirLengthFactor  = 1.d0  !% multiplier of the default overflow weir length that is derived from junction plan area
    setting%Junction%Overflow%CbroadCrestedWeir = 1.5d0 !% Brater and King, Table 5.1   (m^{1/2}/s)   
    !% --- coef1 used in Q = CLH^{3/2} for broad crested weir where L = 2(pi A)^{1/2}
    setting%Junction%Overflow%coef1 = twoR * setting%Junction%Overflow%CbroadCrestedWeir * sqrt(setting%Constant%pi)
    !% --- coef2 is for dQ/dH = (3/2) C L H^{1/2} = (3/2)(2 C * sqrt(pi)) * sqrt(A) * sqrt(H)
    setting%Junction%Overflow%coef2 = threehalfR * setting%Junction%Overflow%coef1 * setting%Junction%Overflow%WeirLengthFactor
    setting%Junction%Overflow%coef3 = twothirdR  * sqrt(twoR * setting%Constant%gravity)
    setting%Junction%Overflow%coef4 = sqrt(twoR * setting%Constant%gravity)

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
            character(kind=json_CK, len=:), allocatable :: c, cvec(:)
            integer              :: ii, integer_value, n_controls, n_cols, n_rows,var_type
            integer              :: len_max, n_cols1, n_rows1
            integer, allocatable :: ilen(:)
            real(8)              :: real_value
            real(8), allocatable :: rvec(:)
            logical :: logical_value
            logical :: found
            logical, pointer :: jsoncheck
            type(json_file)  :: json
            type(json_core)  :: core
            type(json_value),pointer :: PointerMatrix,row

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
            write(*,"(A)") "USER CONFIGURATION ERROR - json file - setting " // 'JSON_FoundFileYN not found'
            write(*,"(A)") "...This should be first item in json file."
            write(*,"(A)") "...There may be a formatting problem inside the json file,"
            write(*,"(A)") "...itself (e.g., mismatched braces or missing comma) or"
            write(*,"(A)") "...this line may be missing from the file. The most common"
            write(*,"(A)") "...cause is a typographical error in the name of, or path."
            write(*,"(A)") "...to, the JSON file."
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
                stop 987051
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
            if (c == 'vshape_surcharge') then
                    setting%Adjust%Head%approach = vshape_surcharge
            else
                write(*,"(A)") 'Error - json file - setting.Adjust.Head.Approach of ',trim(c)
                write(*,"(A)") '..is not in allowed options of:'
                write(*,"(A)") '...vshape_surcharge '
                stop 566339
            end if
        end if
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " //  'Adjust.Head.Approach not found'
        
        !%                       Head.Coef
        call json%get('Adjust.Head.Coef', real_value, found)
        if (found) setting%Adjust%Head%Coef = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Adjust.Head.Coef not found'

    !% AirTracking. ============================================================
    !% --- AirTracking.UseAirTrackingYN
        !%                       UseAirTrackingYN
        call json%get('AirTracking.UseAirTrackingYN', logical_value, found)
        if (found) setting%AirTracking%UseAirTrackingYN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'AirTracking.UseAirTrackingYN not found'

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

    !% Climate. =====================================================================
        !%                        Climate.useHydraulicsEvaporationTF
        call json%get('Climate.useHydraulicsEvaporationTF',logical_value, found)
        if (found) setting%Climate%useHydraulicsEvaporationTF = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Climate.useHydraulicsEvaporationTF not found'

        !%                          Climate.HydraulicsOnlyIntervalHours
        call json%get('Climate.HydraulicsOnlyIntervalHours', real_value, found)
        if (found) setting%Climate%HydraulicsOnlyIntervalHours = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Climate.HydraulicsOnlyIntervalHours not found'
   
        !% do not read Climate.LastTimeUpdate or Climate.NextTimeUpdate

    !% Constant. =====================================================================
        !%                       Constant.gravity
        call json%get('Constant.gravity', real_value, found)
        if (found) setting%Constant%gravity = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Constant.gravity not found'
    
        !%                       Constant.water_temperature
        call json%get('Constant.water_temperature', real_value, found)
        if (found) setting%Constant%water_temperature = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Constant.water_temperature not found'
        
        !%                       Constant.water_kinematic_viscosity
        call json%get('Constant.water_kinematic_viscosity', real_value, found)
        if (found) setting%Constant%water_kinematic_viscosity = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Constant.water_kinematic_viscosity not found'


    !% Control. =====================================================================
    !%                       Control.NumControl
        call json%info('Control.Links',found,var_type,n_controls)
        if (found) setting%Control%NumControl = n_controls
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Control.Links not found'

        
        if (n_controls > zeroI) then
            !%                       Control.Links
            call json%get('Control.Links', cvec, ilen=ilen)
            len_max = maxval(ilen)
            allocate(character(len_max) :: setting%Control%Links(n_controls))
            setting%Control%Links = cvec

            !%                       Control.TimeArray
            !% assuming data stored by column, each column has the same number of elements, and is the same data type
            call json%info('Control.TimeArray',found,var_type,n_cols)
            if (.not. found) stop 'error: Control.TimeArray not found'
            if(n_cols /= n_controls) stop 'error: Control.TimeArray has more/less columns than number of controlled links'

            call json%info('Control.TimeArray(1)',found,var_type,n_rows)
            if (.not. found) stop 'error: Control.TimeArray(1) setting%Control%ElemIdxnot found'

            !% get a pointer to the wind matrix:
            call json%get('Control.TimeArray',PointerMatrix)
            if (.not. associated(PointerMatrix)) stop "Error - json file - setting " // 'Control.TimeArray not found'

            allocate(setting%Control%TimeArray(n_rows,n_cols))
            !% grab each column of the PointerMatrix:
            do ii=1,n_cols
                call core%get_child(PointerMatrix,ii,row)
                if (.not. associated(row)) stop 'error: TimeArray column not found'
                call core%get(row,rvec)
                if (.not. allocated(rvec)) stop 'error: could not get TimeArray real column'
                if (size(rvec)/=n_rows) stop 'error: Control.TimeArray column is wrong size'
                setting%Control%TimeArray(:,ii) = rvec
                deallocate(rvec)
                nullify(row)
            end do
            nullify(PointerMatrix)

            !%                       Control.SettingsArray
            !% assuming data stored by column, each column has the same number of elements, and is the same data type
            call json%info('Control.SettingsArray',found,var_type,n_cols1)
            if (.not. found) stop 'error: Control.SettingsArray not found'
            if(n_cols1 /= n_controls) stop 'error: Control.SettingsArray has more/less columns than number of controlled links'

            call json%info('Control.SettingsArray(1)',found,var_type,n_rows1)
            if(n_rows1 /= n_rows) stop 'error: Control.SettingsArray has more/less rows than the Control.TimeArray'
            if (.not. found) stop 'error: Control.SettingsArray(1) not found'

            !% get a pointer to the PointerMatrix:
            call json%get('Control.SettingsArray',PointerMatrix)
            if (.not. associated(PointerMatrix)) stop "Error - json file - setting " // 'Control.SettingsArray not found'

            allocate(setting%Control%SettingsArray(n_rows1,n_cols1))
            !% grab each column of the PointerMatrix:
            do ii=1,n_cols1
                call core%get_child(PointerMatrix,ii,row)
                if (.not. associated(row)) stop "Error - json file - setting " // 'Control.SettingsArray column not found'
                call core%get(row,rvec)
                if (.not. allocated(rvec)) stop "Error - json file - setting " // 'Could not get SettingsArray real column'
                if (size(rvec)/=n_rows1) stop "Error - json file - setting " // 'Control.SettingsArray is wrong size'
                setting%Control%SettingsArray(:,ii) = rvec
                deallocate(rvec)
                nullify(row)
            end do
            nullify(PointerMatrix)

        end if

    !% Crash =====================================================================
        !%                       Crash. are set by code    


    !% Culvert =====================================================================
    !%                          Culvert.UseCulvertsTF
        call json%get('Culvert.UseCulvertsTF',logical_value, found)
        if (found) setting%Culvert%UseCulvertsTF = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Culvert.UseCulvertsTF not found'

        
    !% Discretization. =====================================================================

        !% -- Channel overflow
        !%                      Discretization.AllowChannelOverflowTF
        call json%get('Discretization.AllowChannelOverflowTF', logical_value, found)
        if (found) setting%Discretization%AllowChannelOverflowTF = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Discretization.AllowChannelOverflowTF not found'
      
        !%                       Discretization.UseNominalElemLength
        call json%get('Discretization.UseNominalElemLength', logical_value, found)
        if (found) setting%Discretization%UseNominalElemLength = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " //'Discretization.UseNominalElemLength not found'

        !%                       Discretization.UseEquivalentOrifice
        call json%get('Discretization.UseEquivalentOrifice', logical_value, found)
        if (found) setting%Discretization%UseEquivalentOrifice = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " //'Discretization.UseEquivalentOrifice not found'

        !%                       Discretization.DistributeOpenChannelInflowsTF
        call json%get('Discretization.DistributeOpenChannelInflowsTF', logical_value, found)
        if (found) setting%Discretization%DistributeOpenChannelInflowsTF = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " //'Discretization.DistributeOpenChannelInflowsTF not found'

        !%                       Discretization.EquivalentOrificeDischargeCoeff
        call json%get('Discretization.EquivalentOrificeDischargeCoeff', real_value, found)
        if (found) setting%Discretization%EquivalentOrificeDischargeCoeff = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " //'Discretization.EquivalentOrificeDischargeCoeff not found'

        !%                       Discretization.NominalElemLength
        call json%get('Discretization.NominalElemLength', real_value, found)
        if (found) setting%Discretization%NominalElemLength = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " //'Discretization.NominalElemLength not found'

        !%                       Discretization.MinLinkLength
        call json%get('Discretization.MinLinkLength', real_value, found)
        if (found) setting%Discretization%MinLinkLength = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " //'Discretization.MinLinkLength not found'

        !%                       Discretization.MinElementPerLink
        call json%get('Discretization.MinElementPerLink', integer_value, found)
        if (found) setting%Discretization%MinElementPerLink = integer_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " //'Discretization.MinElementPerLink not found'

    !% Eps. =====================================================================
        !% --- Eps Settings
        !%                       Eps.FroudeJump
        call json%get('Eps.FroudeJump', real_value, found)
        if (found) setting%Eps%FroudeJump = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Eps.FroudeJump not found'
        

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
        
        !%                       Junction.ForceNodesJM
        call json%get('Junction.ForceNodesJM', logical_value, found)
        if (found) setting%Junction%ForceNodesJM = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Junction.ForceNodesJM not found'

        !%                       Junction.ForceStorage
        call json%get('Junction.ForceStorage', logical_value, found)
        if (found) setting%Junction%ForceStorage = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Junction.ForceStorage not found'

        !%                       Junction.FunStorageN
        call json%get('Junction.FunStorageN', integer_value, found)
        if (found) setting%Junction%FunStorageN = integer_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Junction.CFLlimit not found'
        
        !%                       Junction.kFactor
        call json%get('Junction.kFactor', real_value, found)
        if (found) setting%Junction%kFactor = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Junction.kFactor not found'

        !%                       Junction.OutflowDampingFactor
        call json%get('Junction.OutflowDampingFactor', real_value, found)
        if (found) setting%Junction%OutflowDampingFactor = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Junction.OutflowDampingFactor not found'

        !%                       Junction.InfiniteExtraDepthValue
        call json%get('Junction.InfiniteExtraDepthValue', real_value, found)
        if (found) setting%Junction%InfiniteExtraDepthValue = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Junction.InfiniteExtraDepthValue not found'

         !%                       Junction.PondingScaleFactor
        call json%get('Junction.PondingScaleFactor', real_value, found)
        if (found) setting%Junction%PondingScaleFactor = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Junction.PondingScaleFactor not found'



        !%                       Junction.Overflow.OrificeLength
        call json%get('Junction.Overflow.OrificeLength', real_value, found)
        if (found) setting%Junction%Overflow%OrificeLength = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Junction.Overflow.OrificeLength not found'

        !%                       Junction.Overflow.OrificeHeight
        call json%get('Junction.Overflow.OrificeHeight', real_value, found)
        if (found) setting%Junction%Overflow%OrificeHeight = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Junction.Overflow.OrificeHeight not found'

        !%                       Junction.Overflow.WeirLengthFactor
        call json%get('Junction.Overflow.WeirLengthFactor', real_value, found)
        if (found) setting%Junction%Overflow%WeirLengthFactor= real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Junction.Overflow.WeirLengthFactor not found'


        !%                       Junction.Overflow.CbroadCrestedWeir
        call json%get('Junction.Overflow.CbroadCrestedWeir', real_value, found)
        if (found) setting%Junction%Overflow%CbroadCrestedWeir = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Junction.Overflow.CbroadCrestedWeir not found'

        !%                       Junction.Overflow.coef1
        call json%get('Junction.Overflow.coef1', real_value, found)
        if (found) setting%Junction%Overflow%coef1 = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Junction.Overflow.coef1 not found'

        !%                       Junction.Overflow.coef2
        call json%get('Junction.Overflow.x', real_value, found)
        if (found) setting%Junction%Overflow%coef2 = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Junction.Overflow.coef2 not found'

        !%                       Junction.Overflow.coef3
        call json%get('Junction.Overflow.coef3', real_value, found)
        if (found) setting%Junction%Overflow%coef3 = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Junction.Overflow.coef3 not found'

        !%                       Junction.Overflow.coef4
        call json%get('Junction.Overflow.coef4', real_value, found)
        if (found) setting%Junction%Overflow%coef4 = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Junction.Overflow.coef4 not found'


        !%                       Junction.UseLargeBranchStorageTF
        call json%get('Junction.PlanArea.UseLargeBranchStorageTF', logical_value, found)
        if (found) setting%Junction%PlanArea%UseLargeBranchStorageTF = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Junction.PlanArea.UseLargeBranchStorageTF not found'

        !%                       Junction.PlanArea_Minimum
        call json%get('Junction.PlanArea.AreaMinimum', real_value, found)
        if (found) setting%Junction%PlanArea%AreaMinimum = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Junction.PlanArea5.AreaMinimum not found'

        !%                       Junction.BreadthFactor
        call json%get('Junction.PlanArea.AreaFactorMaximum', real_value, found)
        if (found) setting%Junction%PlanArea%AreaFactorMaximum = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Junction.PlanArea.AreaFactorMaximum not found'

        !%                       Junction.LargeBranchDepth
        call json%get('Junction.PlanArea.LargeBranchDepth', real_value, found)
        if (found) setting%Junction%PlanArea%LargeBranchDepth = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Junction.PlanArea.LargeBranchDepth not found'


    !% Limiter. =====================================================================
    
        !%                       Limiter.Dt.Minimum
        call json%get('Limiter.Dt.Minimum', real_value, found)
        if (found) setting%Limiter%Dt%Minimum = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Limiter.Dt.Minimum not found'
        
        !%                       Limiter.Dt.UseLimitMinYN
        call json%get('Limiter.Dt.UseLimitMinYN', logical_value, found)
        if (found) setting%Limiter%Dt%UseLimitMinYN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Limiter.Dt.UseLimitMinYN not found'

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

         !%                       Limiter.VolumeFractionInTimeStep
        call json%get('Limiter.VolumeFractionInTimeStep', real_value, found)
        if (found) setting%Limiter%VolumeFractionInTimeStep = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Limiter.VolumeFractionInTimeStep not found'

        !%                       Limiter.Flowrate.UseSWMMlinkValueYN
        call json%get('Limiter.Flowrate.UseSWMMlinkValueYN', logical_value, found)
        if (found) setting%Limiter%Flowrate%UseSWMMlinkValueYN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Limiter.Flowrate.UseSWMMlinkValueYN not found'

        !%                       Limiter.Flowrate.UseLocalVolumeYN
        call json%get('Limiter.Flowrate.UseLocalVolumeYN', logical_value, found)
        if (found) setting%Limiter%Flowrate%UseLocalVolumeYN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Limiter.Flowrate.UseLocalVolumeYN not found'

        !%                       Limiter.Flowrate.LocalVolumeFactor
        call json%get('Limiter.Flowrate.LocalVolumeFactor', real_value, found)
        if (found) setting%Limiter%Flowrate%LocalVolumeFactor = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Limiter.Flowrate.LocalVolumeFactornot found'
        

    !% Link. =====================================================================
        !%                      Link.DefaultInitDepthType       
        call json%get('Link.DefaultInitDepthType', c, found)
        if (found) then            
            call util_lower_case(c)
            if (c == 'linear') then
                setting%Link%DefaultInitDepthType = LinearlyVaryingDepth
            else if (c == 'uniform') then
                setting%Link%DefaultInitDepthType = UniformDepth
            else if (c == 'exponential') then
                setting%Link%DefaultInitDepthType = ExponentialDepth
            else if (c == 'fixedhead') then
                setting%Link%DefaultInitDepthTYpe = FixedHead
            else
                write(*,"(A)") 'Error - json file - setting.Link.DefaultInitDepthType of ',trim(c)
                write(*,"(A)") '..is not in allowed options of:'
                write(*,"(A)") '... linear, uniform, exponential, fixedhead'
                stop 93775
            end if
        end if
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Link.DefaultInitDepthType not found'

        !%                      Link.OpenChannelLimitDepthYN
        call json%get('Link.OpenChannelLimitDepthYN', logical_value, found)
        if (found)setting%Link%OpenChannelLimitDepthYN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Link.OpenChannelLimitDepthYN not found'

        !%                      Link.OpenChannelFullDepth
        call json%get('Link.OpenChannelFullDepth', real_value, found)
        if (found) setting%Link%OpenChannelFullDepth = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - setting " // 'Link.OpenChannelFullDepth not found'

    !% Orifice. =====================================================================
        !%                      ForceOrificeNodesToJM
        call json%get('Orifice.ForceOrificeNodesToJM', logical_value, found)
        if (found) setting%Orifice%ForceOrificeNodesToJM = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Orifice.ForceOrificeNodesToJM not found'

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
       
        !% do not read           ElementsExist_byImage
        !% do not read           FacesExist_byImage

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

        !%                       Dataout.isFlowrateAvgOut
        call json%get('Output.DataOut.isFlowrateAvgOut', logical_value, found)
        if (found) setting%Output%DataOut%isFlowrateAVgOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isFlowrateAvgOut not found'

        !%                       Dataout.isFluxConsOut
        call json%get('Output.DataOut.isFluxConsOut', logical_value, found)
        if (found) setting%Output%DataOut%isFluxConsOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isFluxConsOut not found'
        
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

        !%                       Dataout.isManningsNout
        call json%get('Output.DataOut.isManningsNout', logical_value, found)
        if (found) setting%Output%DataOut%isManningsNout = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isManningsNout not found'
        
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

        !%                       Dataout.isVolumeConsOut
        call json%get('Output.DataOut.isVolumeConsOut', logical_value, found)
        if (found) setting%Output%DataOut%isVolumeConsOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isVolumeConsOut not found'

        !%                       Dataout.isVolumeOverflowOut
        call json%get('Output.DataOut.isVolumeOverflowOut', logical_value, found)
        if (found) setting%Output%DataOut%isVolumeOverflowOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isVolumeOverflowOut not found'

        !%                       Dataout.isVolumePondedOut
        call json%get('Output.DataOut.isVolumePondedOut', logical_value, found)
        if (found) setting%Output%DataOut%isVolumePondedOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isVolumePondedOut not found'

        !%                       Dataout.isWaveSpeedOut
        call json%get('Output.DataOut.isWaveSpeedOut', logical_value, found)
        if (found) setting%Output%DataOut%isWaveSpeedOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isWaveSpeedOut not found'

        !%                       Dataout.isPreissmannCelerityOut
        call json%get('Output.DataOut.isPreissmannCelerityOut', logical_value, found)
        if (found) setting%Output%DataOut%isPreissmannCelerityOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isPreissmannCelerityOut not found'

        !%                       Dataout.isPreissmannNumberOut
        call json%get('Output.DataOut.isPreissmannNumberOut', logical_value, found)
        if (found) setting%Output%DataOut%isPreissmannNumberOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isPreissmannNumberOut not found'
             
        !%                       Dataout.isElemLengthOut
        call json%get('Output.DataOut.isElemLengthOut', logical_value, found)
        if (found) setting%Output%DataOut%isElemLengthOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isElemLengthOut not found'
         
         !%                       Dataout.isElemBottomSlopeOut
        call json%get('Output.DataOut.isElemBottomSlopeOut', logical_value, found)
        if (found) setting%Output%DataOut%isElemBottomSlopeOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isElemBottomSlopeOut not found'
         
         !%                       Dataout.isElemBreathMaxOut
        call json%get('Output.DataOut.isElemBreathMaxOut', logical_value, found)
        if (found) setting%Output%DataOut%isElemBreathMaxOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isElemBreathMaxOut not found'

        !%                       Dataout.isElemFullAreaOut
        call json%get('Output.DataOut.isElemFullAreaOut', logical_value, found)
        if (found) setting%Output%DataOut%isElemFullAreaOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isElemFullAreaOut not found'

        !%                       Dataout.isElemFullDepthOut
        call json%get('Output.DataOut.isElemFullDepthOut', logical_value, found)
        if (found) setting%Output%DataOut%isElemFullDepthOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isElemFullDepthOut not found'

        !%                       Dataout.isElemManningsOut
        call json%get('Output.DataOut.isElemManningsOut', logical_value, found)
        if (found) setting%Output%DataOut%isElemManningsOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isElemManningsOut not found'

        !%                       Dataout.isElemZBottomOut
        call json%get('Output.DataOut.isElemZBottomOut', logical_value, found)
        if (found) setting%Output%DataOut%isElemZBottomOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isElemZBottomOut not found'

        !%                       Dataout.isElemZCrownOut
        call json%get('Output.DataOut.isElemZCrownOut', logical_value, found)
        if (found) setting%Output%DataOut%isElemZCrownOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isElemZCrownOut not found'

        !%                       Dataout.isLinkLengthOut
        call json%get('Output.DataOut.isLinkLengthOut', logical_value, found)
        if (found) setting%Output%DataOut%isLinkLengthOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isLinkLengthOut not found'

        !%                       Dataout.isLinkAdjustedLengthOut
        call json%get('Output.DataOut.isLinkAdjustedLengthOut', logical_value, found)
        if (found) setting%Output%DataOut%isLinkAdjustedLengthOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isLinkAdjustedLengthOut not found'

        !%                       Dataout.isLinkInletOffsetOut
        call json%get('Output.DataOut.isLinkInletOffsetOut', logical_value, found)
        if (found) setting%Output%DataOut%isLinkInletOffsetOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isLinkInletOffsetOut not found'

        !%                       Dataout.isLinkOutletOffsetOut
        call json%get('Output.DataOut.isLinkOutletOffsetOut', logical_value, found)
        if (found) setting%Output%DataOut%isLinkOutletOffsetOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isLinkOutletOffsetOut not found'

        !%                       Dataout.isLinkSlopeOut
        call json%get('Output.DataOut.isLinkSlopeOut', logical_value, found)
        if (found) setting%Output%DataOut%isLinkSlopeOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isLinkSlopeOut not found'

        !%                       Dataout.isLinkZBottomUpOut
        call json%get('Output.DataOut.isLinkZBottomUpOut', logical_value, found)
        if (found) setting%Output%DataOut%isLinkZBottomUpOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isLinkZBottomUpOut not found'
                 
        !%                       Dataout.isLinkZBottomDownOut
        call json%get('Output.DataOut.isLinkZBottomDownOut', logical_value, found)
        if (found) setting%Output%DataOut%isLinkZBottomDownOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isLinkZBottomDownOut not found'

        !%                       Dataout.isNodeZBottomOut
        call json%get('Output.DataOut.isNodeZBottomOut', logical_value, found)
        if (found) setting%Output%DataOut%isNodeZBottomOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isNodeZBottomOut not found'

        !%                       Dataout.isNodeFullDepthOut
        call json%get('Output.DataOut.isNodeFullDepthOut', logical_value, found)
        if (found) setting%Output%DataOut%isNodeFullDepthOut = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.DataOut.isNodeFullDepthOut not found'
    


        !% --- Report settings
        !%                       Report.useSWMMinpYN
        call json%get('Output.Report.useSWMMinpYN', logical_value, found)
        if (found) setting%Output%Report%useSWMMinpYN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.Report.useSWMMinpYN not found'

        !%                       Report.provideYN
        call json%get('Output.Report.provideYN', logical_value, found)
        if (found) setting%Output%Report%provideYN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.Report.provideYN not found'

        !%                       Report.useHDF5
        call json%get('Output.Report.useHDF5', logical_value, found)
        if (found) setting%Output%Report%useHDF5 = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.Report.useHDF5 not found'

        !%                       Report.useCSV
        call json%get('Output.Report.useCSV', logical_value, found)
        if (found) setting%Output%Report%useCSV = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Output.Report.useCSV not found'


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

        !%                       PhantomLinkAdjust
        call json%get('Partitioning.PhantomLinkAdjust', logical_value, found)
        if (found) setting%Partitioning%PhantomLinkAdjust = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Partitioning.PhantomLinkAdjust not found'
  
    !% Profile. =====================================================================
        !%                       Profile.useYN
        call json%get('Profile.useYN', logical_value, found)
        if (found) setting%Profile%useYN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Profile.useYN not found'

    !% Simulation. =====================================================================
        !%                       AllowReverseGradientInitialConditionsTF
        call json%get('Simulation.AllowReverseGradientInitialConditionsTF', logical_value, found)
        if (found) setting%Simulation%AllowReverseGradientInitialConditionsTF = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Simulation.AllowReverseGradientInitialConditionsTF not found'
          !%                       useHydrology
        call json%get('Simulation.useHydrology', logical_value, found)
        if (found) setting%Simulation%useHydrology = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Simulation.useHydrology not found'
        
        !%                       useHydraulics
        call json%get('Simulation.useHydraulics', logical_value, found)
        if (found) setting%Simulation%useHydraulics = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Simulation.useHydraulics not found'

        !%                       stopAfterInitializationYN
        call json%get('Simulation.stopAfterInitializationYN', logical_value, found)
        if (found) setting%Simulation%stopAfterInitializationYN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Simulation.stopAfterInitializationYN not found'

        !%                       useSpinUp
        call json%get('Simulation.useSpinUp', logical_value, found)
        if (found) setting%Simulation%useSpinUp = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Simulation.useSpinUp not found'

        !%                       stopAfterSpinUp
        call json%get('Simulation.stopAfterSpinUp', logical_value, found)
        if (found) setting%Simulation%stopAfterSpinUp = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Simulation.stopAfterSpinUp not found'

        !%                       SpinUpDays
        call json%get('Simulation.SpinUpDays', real_value, found)
        if (found) setting%Simulation%SpinUpDays = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Simulation.SpinUpDays not found'
  
    !% SmallDepth. =====================================================================    

                !%                       useMomentumCutoffYN
        call json%get('SmallDepth.useMomentumCutoffYN', logical_value, found)
        if (found) setting%SmallDepth%useMomentumCutoffYN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'SmallDepth.useMomentumCutoffYN not found'
        
        !%                       MomentumDepthCutoff
        call json%get('SmallDepth.MomentumDepthCutoff', real_value, found)
        if (found) setting%SmallDepth%MomentumDepthCutoff = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'SmallDepth.MomentumDepthCutoff not found'
        
        !%                       ManningsN
        call json%get('SmallDepth.ManningsN', real_value, found)
        if (found) setting%SmallDepth%ManningsN = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'SmallDepth.ManningsN not found'
        
        !%                       PumpVolumeFactor
        call json%get('SmallDepth.PumpVolumeFactor', real_value, found)
        if (found) setting%SmallDepth%PumpVolumeFactor = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'SmallDepth.PumpVolumeFactor not found'

        !%                       LateralInflowSmallDepth
        call json%get('SmallDepth.LateralInflowSmallDepth', real_value, found)
        if (found) setting%SmallDepth%LateralInflowSmallDepth = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'SmallDepth.LateralInflowSmallDepth not found'

    !% Solver. =====================================================================
       
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
            else if (c == 't2l0') then
                setting%Solver%MomentumSourceMethod = T2L0
            ! else if (c == 't10s2') then
            !     setting%Solver%MomentumSourceMethod = T10s2
            ! else if (c == 'ta1') then
            !     setting%Solver%MomentumSourceMethod = TA1
            ! else if (c == 'ta2') then
            !     setting%Solver%MomentumSourceMethod = TA2
            else
                write(*,"(A)") 'Error - json file - setting.Solver.MomentumSourceMethod of ',trim(c)
                write(*,"(A)") '..is not in allowed options of:'
                write(*,"(A)") '... T00, T10, T20, T2L0'
                stop 110985
            end if
        end if
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Solver.MomentumSourceMethod not found'     

        !%                       Solver.SwitchFractionDn
        call json%get('Solver.SwitchFractionDn', real_value, found)
        if (found)  setting%Solver%SwitchFractionDn = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Solver.SwitchFractionDn not found'
        
        !%                       Solver.SwitchFractionU
        call json%get('Solver.SwitchFractionUp', real_value, found)
        if (found)   setting%Solver%SwitchFractionUp = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Solver.SwitchFractionUp not found'

        !% do not read           Solver.crk2

        !% do not read          Solver.ReferenceHead

      
    !% Solver.ForceMain =====================================================================
        !%                       Solver.ForceMain.AllowForceMainTF
        call json%get('Solver.ForceMain.AllowForceMainTF', logical_value, found)
        if (found) setting%Solver%ForceMain%AllowForceMainTF = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Solver.ForceMain.AllowForceMainTF not found'

        !%                       Solver.ForceMain.UseSWMMinputMethodTF
        call json%get('Solver.ForceMain.UseSWMMinputMethodTF', logical_value, found)
        if (found) setting%Solver%ForceMain%UseSWMMinputMethodTF = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Solver.ForceMain.x\UseSWMMinputMethodTF not found'

        !%                       Solver.ForceMain.FMallClosedConduitsTF
        call json%get('Solver.ForceMain.FMallClosedConduitsTF', logical_value, found)
        if (found) setting%Solver%ForceMain%FMallClosedConduitsTF = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Solver.ForceMain.FMallClosedConduitsTF not found'

         !%                       Solver.ForceMain.FMallClosedConduitsTF
        call json%get('Solver.ForceMain.errorCheck_RoughnessTF', logical_value, found)
        if (found) setting%Solver%ForceMain%errorCheck_RoughnessTF = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Solver.ForceMain.errorCheck_RoughnessTF not found'
  
        !%                       Solver.ForceMain.Default_method
        call json%get('Solver.ForceMain.Default_method', c, found)
        if (found) then
            call util_lower_case(c)
            select case (trim(c))
            case ('hazenwilliams')
                setting%Solver%ForceMain%Default_method = HazenWilliams
            case ('darcyweisbach')
                setting%Solver%ForceMain%Default_method = DarcyWeisbach
            case default 
                write(*,"(A)") 'Error - json file - setting.Solver.ForceMain.Default_method of ',trim(c)
                write(*,"(A)") '..is not in allowed options of:'
                write(*,"(A)") '... HazenWilliams, DarcyWeisbach'
                stop 9375467
            end select
        end if
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Solver.ForceMain.Default_method  not found'     
      
        !%                       Solver.ForceMain.Default_HazenWilliams_coef
        call json%get('Solver.ForceMain.Default_HazenWilliams_coef', real_value, found)
        if (found)  setting%Solver%ForceMain%Default_HazenWilliams_coef = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Solver.ForceMain.Default_HazenWilliams_coef not found'

        !%                       Solver.ForceMain.Default_DarcyWeisbach_roughness_mm
        call json%get('Solver.ForceMain.Default_DarcyWeisbach_roughness_mm', real_value, found)
        if (found)  setting%Solver%ForceMain%Default_DarcyWeisbach_roughness_mm = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Solver.ForceMain.Default_DarcyWeisbach_roughness_mm not found'

        !%                       Solver.ForceMain.Default_ManningsN
        call json%get('Solver.ForceMain.Default_ManningsN', real_value, found)
        if (found)  setting%Solver%ForceMain%Default_ManningsN = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Solver.ForceMain.Default_ManningsN not found'

        !%                       Solver.ForceMain.minimum_slope
        call json%get('Solver.ForceMain.minimum_slope', real_value, found)
        if (found)  setting%Solver%ForceMain%minimum_slope = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Solver.ForceMain.minimum_slope not found'

    !% Solver.PreissmannSlot. =====================================================================
         !%                       Solver.PreissmannSlot.useSlotTF
        call json%get('Solver.PreissmannSlot.useSlotTF', logical_value, found)
        if (found) setting%Solver%PreissmannSlot%useSlotTF = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Solver.PreissmannSlot.useSlotTF not found'

        !%                       Solver.PreissmannSlot.SlotMethod
        call json%get('Solver.PreissmannSlot.Method', c, found)
        if (found) then 
            call util_lower_case(c)
            if (c == 'staticslot') then
                setting%Solver%PreissmannSlot%Method = StaticSlot
            else if (c == 'dynamicslot') then
                setting%Solver%PreissmannSlot%Method = DynamicSlot
            else if (c == 'splitdynamicslot') then
                setting%Solver%PreissmannSlot%Method = SplitDynamicSlot
            else
                write(*,"(A)") 'Error - json file - setting.Solver.PreissmannSlot%Method of ',trim(c)
                write(*,"(A)") '..is not in allowed options of:'
                write(*,"(A)") '... StaticSlot, DynamicSlot, SplitDynamicSlot'
                stop 22496
            end if
        end if
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Solver.PreissmannSlot.SlotMethod not found'

        !%                      Solver.PreissmannSlot.TargetCelerity
        call json%get('Solver.PreissmannSlot.TargetCelerity', real_value, found)
        if (found) setting%Solver%PreissmannSlot%TargetCelerity = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Solver.PreissmannSlot.TargetCelerity not found'

        !%                      Alpha
        call json%get('Solver.PreissmannSlot.Alpha', real_value, found)
        if (found) setting%Solver%PreissmannSlot%Alpha = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Solver.PreissmannSlot.Alpha not found'

        !%                      DecayRate
        call json%get('Solver.PreissmannSlot.DecayRate', real_value, found)
        if (found) setting%Solver%PreissmannSlot%DecayRate = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Solver.PreissmannSlot.DecayRate not found'

      !%                      Minimum initial preissmann number
        call json%get('Solver.PreissmannSlot.initPNminimum', real_value, found)
        if (found) setting%Solver%PreissmannSlot%initPNminimum = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Solver.PreissmannSlot.initPNminimum not found'


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

        !%                      ForceWeirNodesToJM
        call json%get('Weir.ForceWeirNodesToJM', logical_value, found)
        if (found) setting%Weir%ForceWeirNodesToJM = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Weir.ForceWeirNodesToJM not found'

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
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'VariableDT.ApplyYN not found'
        
         !%                       limitByBC_YN
        call json%get('VariableDT.limitByBC_YN', logical_value, found)
        if (found) setting%VariableDT%limitByBC_YN = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'VariableDT.limitByBC_YN not found'

        !%                       CFL_hi
        call json%get('VariableDT.CFL_hi', real_value, found)
        if (found) setting%VariableDT%CFL_hi = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'VariableDT.CFL_hi not found'
        
        !%                       CFL_target
        call json%get('VariableDT.CFL_target', real_value, found)
        if (found) setting%VariableDT%CFL_target = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'VariableDT.CFL_target not found'
        
        !%                       CFL_lo
        call json%get('VariableDT.CFL_lo', real_value, found)
        if (found) setting%VariableDT%CFL_lo = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'VariableDT.CFL_lo not found'
        
         !%                       CFL_inflow_max
        call json%get('VariableDT.CFL_inflow_max', real_value, found)
        if (found) setting%VariableDT%CFL_lo = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'VariableDT.CFL_inflow_max not found'
        
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

          !%                       InitialDt
        call json%get('VariableDT.InitialDt', real_value, found)
        if (found) setting%VariableDT%InitialDt = real_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'VariableDT.InitialDt not found'

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

        !% do not read VolumeResetLevel

    !% Debug. =====================================================================
        !%                       isGlobalVolumeBalance 
        call json%get('Debug.isGlobalVolumeBalance', logical_value, found)
        if (found) setting%Debug%isGlobalVolumeBalance = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.isGlobalVolumeBalance not found'
        
        !%                       checkIsNanTF
        call json%get('Debug.checkIsNanTF', logical_value, found)
        if (found) setting%Debug%checkIsNanTF = logical_value
        if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.checkIsNanTF not found'
        
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
        ! call json%get('Debug.FileGroup.all', logical_value, found)
        ! if (found) setting%Debug%FileGroup%all = logical_value
        ! !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.FileGroup.all not found'
        
        ! !%                       
        ! call json%get('Debug.FileGroup.definitions', logical_value, found)
        ! if (found) setting%Debug%FileGroup%definitions = logical_value
        ! !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.FileGroup.definitions not found'
        
        ! !%                       
        ! call json%get('Debug.FileGroup.finalization', logical_value, found)
        ! if (found) setting%Debug%FileGroup%finalization = logical_value
        ! !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.FileGroup.finalization not found'
        
        ! !%                       
        ! call json%get('Debug.FileGroup.geometry', logical_value, found)
        ! if (found) setting%Debug%FileGroup%geometry = logical_value
        ! !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.FileGroup.geometry not found'
        
        ! !%                       
        ! call json%get('Debug.FileGroup.initialization', logical_value, found)
        ! if (found) setting%Debug%FileGroup%initialization = logical_value
        ! !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.FileGroup.initialization not found'
        
        ! !%                       
        ! call json%get('Debug.FileGroup.interface', logical_value, found)
        ! if (found) setting%Debug%FileGroup%interface = logical_value
        ! !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.FileGroup.interface not found'
        
        ! !%                       
        ! call json%get('Debug.FileGroup.output', logical_value, found)
        ! if (found) setting%Debug%FileGroup%output = logical_value
        ! !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.FileGroup.output not found'
        
        ! !%                       
        ! call json%get('Debug.FileGroup.timeloop', logical_value, found)
        ! if (found) setting%Debug%FileGroup%timeloop = logical_value
        ! !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.FileGroup.timeloop not found'
        
        ! !%                       
        ! call json%get('Debug.FileGroup.utility', logical_value, found)
        ! if (found) setting%Debug%FileGroup%utility = logical_value
        ! !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.FileGroup.utility not found'
        
    !% Debug. =====================================================================
        !%                       'Debug.SetupYN
        ! call json%get('Debug.SetupYN', logical_value, found)
        ! if (found) setting%Debug%SetupYN = logical_value
        ! !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.Setup not found'
        
        ! !%                       Debug.OutputYN
        ! call json%get('Debug.OutputYN', logical_value, found)
        ! if (found) setting%Debug%OutputYN = logical_value
        ! !if ((.not. found) .and. (jsoncheck)) stop "Error - json file - setting " // 'Debug.Output not found'

        ! call define_settings_update_debug()
 
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
    ! subroutine define_settings_update_debug()
    !     !% -----------------------------------------------------------------
    !     !% Description:
    !     !% Uses the Debug%FileGroup... = true option to turn on the debug
    !     !% for all the files in that group
    !     !% -----------------------------------------------------------------
    !     if (setting%Debug%FileGroup%all) then
    !         setting%Debug%FileGroup%definitions = .true.
    !         setting%Debug%FileGroup%finalization = .true.
    !         setting%Debug%FileGroup%geometry = .true.
    !         setting%Debug%FileGroup%initialization = .true.
    !         setting%Debug%FileGroup%interface = .true.
    !         setting%Debug%FileGroup%output  = .true.
    !         setting%Debug%FileGroup%timeloop  = .true.
    !         setting%Debug%FileGroup%utility = .true.
    !     end if
    !     if (setting%Debug%FileGroup%definitions) then
    !         setting%Debug%File%define_globals = .true.
    !         setting%Debug%File%define_indexes = .true.
    !         setting%Debug%File%define_keys = .true.
    !         setting%Debug%File%define_settings = .true.
    !         setting%Debug%File%define_types = .true.
    !     end if
    !     if (setting%Debug%FileGroup%finalization) then
    !         setting%Debug%File%finalization = .true.
    !     end if
    !     if (setting%Debug%FileGroup%geometry) then
    !         setting%Debug%File%geometry = .true.
    !         setting%Debug%File%rectangular_channel = .true.
    !         setting%Debug%File%trapezoidal_channel = .true.
    !     end if
    !     if (setting%Debug%FileGroup%initialization) then
    !         setting%Debug%File%discretization = .true.
    !         setting%Debug%File%initial_condition = .true.
    !         setting%Debug%File%initialization = .true.
    !         setting%Debug%File%network_define = .true.
    !         setting%Debug%File%partitioning = .true.
    !         setting%Debug%File%BIPquick = .true.
    !         setting%Debug%File%pack_mask_arrays = .true.
    !     end if
    !     if (setting%Debug%FileGroup%interface) then
    !         setting%Debug%File%c_library = .true.
    !         setting%Debug%File%interface = .true.
    !     end if
    !     if (setting%Debug%FileGroup%output) then
    !         setting%Debug%File%output = .true.
    !     end if
    !     if (setting%Debug%FileGroup%timeloop) then
    !         setting%Debug%File%adjust = .true.
    !         setting%Debug%File%boundary_conditions = .true.
    !         setting%Debug%File%diagnostic_elements = .true.
    !         setting%Debug%File%face = .true.
    !         setting%Debug%File%jump = .true.
    !         setting%Debug%File%lowlevel_rk2 = .true.
    !         setting%Debug%File%orifice_elements = .true.
    !         setting%Debug%File%pump_elements = .true.
    !         setting%Debug%File%runge_kutta2 = .true.
    !         setting%Debug%File%timeloop = .true.
    !         setting%Debug%File%update = .true.
    !         setting%Debug%File%weir_elements = .true.
    !     end if
    !     if (setting%Debug%FileGroup%utility) then
    !         setting%Debug%File%utility_allocate = .true.
    !         setting%Debug%File%utility_deallocate = .true.
    !         setting%Debug%File%utility_array = .true.
    !         setting%Debug%File%utility_datetime = .true.
    !         setting%Debug%File%utility_interpolate = .true.
    !         setting%Debug%File%utility_output = .true.
    !         setting%Debug%File%utility_string = .true.
    !         setting%Debug%File%utility = .true.
    !     end if
    ! end subroutine define_settings_update_debug
!%==========================================================================
!% End of module
!%==========================================================================
end module define_settings
