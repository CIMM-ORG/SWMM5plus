module define_api_keys

    implicit none

    !% --------------------------------------------------------
    !% Interface keys
    !%
    !%   The following UPPERCASE variables must have the same values of
    !%   the variables in the EPA-SWMM source code. The especific
    !%   C file that is associated to each enumerator is
    !%   written at the beginning of each 'enum' declaration.
    !%
    !%   The lowercase variables must have the same enum values
    !%   as those in api.h
    !%
    !%   The variables have been renamed with respect to the
    !%   EPA-SWMM source code to avoid name conflicts. Concretely,
    !%   the 'API_' substring was added to the name at the
    !%   the beginning of the API key
    !%
    !%   Variables written in caps are extracted from native
    !%   EPA-SWMM, whereas lower case vars are used in the
    !%   SWMM5+ fortran code
    !% --------------------------------------------------------

    character(len=32), allocatable :: reverseKey_api(:)

    ! !% SWMM Objects ($API_DIR/src/enums.h)
    ! enum, bind(c)
    !     enumerator :: API_NODE = 2
    !     enumerator :: API_LINK
    ! end enum

    !% SWMM Objects ($API_DIR/src/enums.h)
    enum, bind(c)
        enumerator :: API_GAGE = 0      ! not used yet
        enumerator :: API_SUBCATCH      ! 1
        enumerator :: API_NODE          ! 2
        enumerator :: API_LINK          ! 3
        enumerator :: API_POLLUT        ! 4
        enumerator :: API_LANDUSE       ! 5
        enumerator :: API_TIMEPATTERN   ! 6
        enumerator :: API_CURVE         ! 7
        enumerator :: API_TSERIES       ! 8
        enumerator :: API_CONTROL       ! 9
        enumerator :: API_TRANSECT      ! 10
        enumerator :: API_AQUIFER       ! 11
        enumerator :: API_UNITHYD       ! 12
        enumerator :: API_SNOWMELT      ! 13
        enumerator :: API_SHAPE         ! 14
        enumerator :: API_LID           ! 15
        enumerator :: API_MAX_OBJ_TYPES ! 16
    end enum

    !% SWMM Node types ($API_DIR/src/enums.h) -> NodeType
    enum, bind(c)
        enumerator :: API_JUNCTION = 0
        enumerator :: API_OUTFALL       ! 1
        enumerator :: API_STORAGE       ! 2
        enumerator :: API_DIVIDER       ! 3
    end enum

    !% SWMM Link types ($API_DIR/src/enums.h) -> LinkType
    enum, bind(c)
        enumerator :: API_CONDUIT = 0
        enumerator :: API_PUMP          ! 1
        enumerator :: API_ORIFICE       ! 2
        enumerator :: API_WEIR          ! 3
        enumerator :: API_OUTLET        ! 4
    end enum

    !% SWMM File types ($API_DIR/src/enums.h)
    enum, bind(c)
        enumerator :: API_RAINFALL_FILE = 0 ! not used yet
        enumerator :: API_RUNOFF_FILE   ! 1 not used yet
        enumerator :: API_HOTSTART_FILE ! 2 not used yet
        enumerator :: API_RDII_FILE     ! 3 not used yet
        enumerator :: API_INFLOWS_FILE  ! 4 not used yet
        enumerator :: API_OUTFLOWS_FILE ! 5 not used yet
    end enum

    !% SWMM File usage types
    enum, bind(c)
        enumerator :: API_NO_FILE = 0   ! not used yet
        enumerator :: API_SCRATCH_FILE  ! 1 not used yet
        enumerator :: API_USE_FILE      ! 2 not used yet
        enumerator :: API_SAVE_FILE     ! 3 not used yet
    end enum

    !% SWMM Rain gage data types
    enum, bind(c)
        enumerator :: API_RAIN_TSERIES = 0  ! not used yet
        enumerator :: API_RAIN_FILE         ! 1 not used yet
    end enum

    !% SWMM XSECT_TYPES ($API_DIR/src/enums.h -> XsectType)
    enum, bind(c)
        enumerator :: API_DUMMY = 0
        enumerator :: API_CIRCULAR          ! 1
        enumerator :: API_FILLED_CIRCULAR   ! 2
        enumerator :: API_RECT_CLOSED       ! 3
        enumerator :: API_RECT_OPEN         ! 4
        enumerator :: API_TRAPEZOIDAL       ! 5
        enumerator :: API_TRIANGULAR        ! 6
        enumerator :: API_PARABOLIC         ! 7
        enumerator :: API_POWERFUNC         ! 8
        enumerator :: API_RECT_TRIANG       ! 9
        enumerator :: API_RECT_ROUND        ! 10
        enumerator :: API_MOD_BASKET        ! 11
        enumerator :: API_HORIZ_ELLIPSE     ! 12
        enumerator :: API_VERT_ELLIPSE      ! 13
        enumerator :: API_ARCH              ! 14
        enumerator :: API_EGGSHAPED         ! 15
        enumerator :: API_HORSESHOE         ! 16
        enumerator :: API_GOTHIC            ! 17
        enumerator :: API_CATENARY          ! 18
        enumerator :: API_SEMIELLIPTICAL    ! 19
        enumerator :: API_BASKETHANDLE      ! 20
        enumerator :: API_SEMICIRCULAR      ! 21
        enumerator :: API_IRREGULAR         ! 22
        enumerator :: API_CUSTOM            ! 23
        enumerator :: API_FORCE_MAIN        ! 24
    end enum

    !% SWMM measurement unit types
    enum, bind(c)
        enumerator :: API_US = 0
        enumerator :: API_SI        ! 1
    end enum

    !% SWMM Flow Units type
    enum, bind(c)
        enumerator :: API_CFS = 0
        enumerator :: API_GPM       ! 1
        enumerator :: API_MGD       ! 2
        enumerator :: API_CMS       ! 3
        enumerator :: API_LPS       ! 4
        enumerator :: API_MLD       ! 5
    end enum

    !% SWMM Concentration units type
    enum, bind(c)
        enumerator :: API_MG = 0
        enumerator :: API_UG        ! 1
        enumerator :: API_COUNT     ! 2
    end enum

    !% SWMM quantities requiring unit conversions
    enum, bind(c)
        enumerator :: API_RAINFALL = 0
        enumerator :: API_RAINDEPTH     ! 1
        enumerator :: API_EVAPRATE      ! 2
        enumerator :: API_LENGTH        ! 3
        enumerator :: API_LANDAREA      ! 4
        enumerator :: API_VOLUME        ! 5
        enumerator :: API_WINDSPEED     ! 6
        enumerator :: API_TEMPERATURE   ! 7
        enumerator :: API_MASS          ! 8
        enumerator :: API_GWFLOW        ! 9
        enumerator :: API_FLOW          ! 10
    end enum

    !% SWMM computed subcatchment quantities
    enum, bind(c)
        enumerator :: API_SUBCATCH_RAINFALL = 0
        enumerator :: API_SUBCATCH_SNOWDEPTH        ! 1
        enumerator :: API_SUBCATCH_EVAP             ! 2
        enumerator :: API_SUBCATCH_INFIL            ! 3
        enumerator :: API_SUBCATCH_RUNOFF           ! 4
        enumerator :: API_SUBCATCH_GW_FLOW          ! 5
        enumerator :: API_SUBCATCH_GW_ELEV          ! 6
        enumerator :: API_SUBCATCH_SOIL_MOISTURE    ! 7
        enumerator :: API_SUBCATCH_WASHOFF          ! 8
    end enum

    !% SWMM node quantites
    enum, bind(c)
        enumerator :: API_NODE_DEPTH = 0
        enumerator :: API_NODE_HEAD         ! 1
        enumerator :: API_NODE_VOLUME       ! 2
        enumerator :: API_NODE_LATFLOW      ! 3
        enumerator :: API_NODE_INFLOW       ! 4
        enumerator :: API_NODE_OVERFLOW     ! 5
        enumerator :: API_NODE_QUAL         ! 6
    end enum

    !% SWMM Link result type
    enum, bind(c)
        enumerator :: API_LINK_FLOW = 0
        enumerator :: API_LINK_DEPTH        ! 1
        enumerator :: API_LINK_VELOCITY     ! 2
        enumerator :: API_LINK_VOLUME       ! 3
        enumerator :: API_LINK_CAPACITY     ! 4
        enumerator :: API_LINK_QUAL         ! 5
    end enum

    !% SWMM system-wide flow quantities
    enum, bind(c)
        enumerator :: API_SYS_TEMPERATURE = 0
        enumerator :: API_SYS_RAINFALL      ! 1
        enumerator :: API_SYS_SNOWDEPTH     ! 2
        enumerator :: API_SYS_INFIL         ! 3
        enumerator :: API_SYS_RUNOFF        ! 4
        enumerator :: API_SYS_DWFLOW        ! 5
        enumerator :: API_SYS_GWFLOW        ! 6
        enumerator :: API_SYS_IIFLOW        ! 7
        enumerator :: API_SYS_EXFLOW        ! 8
        enumerator :: API_SYS_INFLOW        ! 9
        enumerator :: API_SYS_FLOODING      ! 10
        enumerator :: API_SYS_OUTFLOW       ! 11
        enumerator :: API_SYS_STORAGE       ! 12
        enumerator :: API_SYS_EVAP          ! 13
        enumerator :: API_SYS_PET           ! 14
    end enum

    !% SWMM conduit flow classifications
    enum, bind(c)
        enumerator :: API_DRY = 0
        enumerator :: API_UP_DRY            ! 1
        enumerator :: API_DN_DRY            ! 2
        enumerator :: API_SUBCRITICAL       ! 3   
        enumerator :: API_SUPCRITICAL       ! 4
        enumerator :: API_UP_CRITICAL       ! 5
        enumerator :: API_DN_CRITICAL       ! 6
        enumerator :: API_MAX_FLOW_CLASSES  ! 7
        enumerator :: API_UP_FULL           ! 8
        enumerator :: API_DN_FULL           ! 9
        enumerator :: API_ALL_FULL          ! 10
    end enum

    !% SWMM  runoff flow categories
    enum, bind(c)
        enumerator :: API_RUNOFF_RAINFALL = 0
        enumerator :: API_RUNOFF_EVAP       ! 1
        enumerator :: API_RUNOFF_INFIL      ! 2
        enumerator :: API_RUNOFF_RUNOFF     ! 3
        enumerator :: API_RUNOFF_DRAINS     ! 4
        enumerator :: API_RUNOFF_RUNON      ! 5
    end enum

    !% SWMM pollutant loading categories
    enum, bind(c)
        enumerator :: API_BUILDUP_LOAD = 0
        enumerator :: API_DEPOSITION_LOAD   ! 1
        enumerator :: API_SWEEPING_LOAD     ! 2
        enumerator :: API_BMP_REMOVAL_LOAD  ! 3
        enumerator :: API_INFIL_LOAD        ! 4
        enumerator :: API_RUNOFF_LOAD       ! 5
        enumerator :: API_FINAL_LOAD        ! 6
    end enum

    !% SWMM input data options - rainfall
    enum, bind(c)
        enumerator :: API_RAINFALL_INTENSITY = 0
        enumerator :: API_RAINFALL_VOLUME       ! 1
        enumerator :: API_CUMULATIVE_RAINFALL   ! 2
    end enum

    !% SWMM input data options -- temperature
    enum, bind(c)
        enumerator :: API_NO_TEMP = 0
        enumerator :: API_TSERIES_TEMP  ! 1
        enumerator :: API_FILE_TEMP     ! 2
    end enum

    !% SWMM input data options -- wind
    enum, bind(c)
        enumerator :: API_MONTHLY_WIND = 0
        enumerator :: API_FILE_WIND     ! 1

    end enum

    !% SWMM input data options -- evaporation 
    enum, bind(c)
        enumerator :: API_CONSTANT_EVAP = 0
        enumerator :: API_MONTHLY_EVAP      ! 1
        enumerator :: API_TIMESERIES_EVAP   ! 2
        enumerator :: API_TEMPERATURE_EVAP  ! 3
        enumerator :: API_FILE_EVAP         ! 4
        enumerator :: API_RECOVERY          ! 5
        enumerator :: API_DRYONLY           ! 6
    end enum

    !% SWMM input data options -- normalizing
    enum, bind(c)
        enumerator :: API_PER_AREA = 0
        enumerator :: API_PER_CURB      ! 1
    end enum

    !% SWMM input data options -- buildup type
    enum, bind(c)
        enumerator :: API_NO_BUILDUP = 0
        enumerator :: API_POWER_BUILDUP     ! 1
        enumerator :: API_EXPON_BUILDUP     ! 2
        enumerator :: API_SATUR_BUILDUP     ! 3
        enumerator :: API_EXTERNAL_BUILDUP  ! 4
    end enum

    !% SWMM input data options -- wash off type
    enum, bind(c)
        enumerator :: API_NO_WASHOFF = 0
        enumerator :: API_EXPON_WASHOFF     ! 1
        enumerator :: API_RATING_WASHOFF    ! 2
        enumerator :: API_EMC_WASHOFF       ! 3
    end enum

    !% SWMM input data options -- SubArea type
    enum, bind(c)
        enumerator :: API_IMPERV0 = 0
        enumerator :: API_IMPERV1       ! 1
        enumerator :: API_PERV          ! 2
    end enum

    !% SWMM input data options -- Runoff Routing type
    enum, bind(c)
        enumerator :: API_TO_OUTLET = 0
        enumerator :: API_TO_IMPERV     ! 1
        enumerator :: API_TO_PERV       ! 2
    end enum

    !% SWMM input data options -- Routine model type
    enum, bind(c)
        enumerator :: API_NO_ROUTING = 0
        enumerator :: API_SF            ! 1
        enumerator :: API_KW            ! 2
        enumerator :: API_EKW           ! 3
        enumerator :: API_DW            ! 4
    end enum

    !% SWMM input data options --  Force main type
    enum, bind(c)
        enumerator :: API_H_W = 0
        enumerator :: API_D_W       ! 1
    end enum

    !% SWMM input data options -- Offset Type
    enum, bind(c)
        enumerator :: API_DEPTH_OFFSET = 0
        enumerator :: API_ELEV_OFFSET       ! 1
    end enum

    !% SWMM input data options -- Kinematic wave method type
    enum, bind(c)
        enumerator :: API_NORMAL = 0
        enumerator :: API_MODIFIED      ! 1
    end enum

    !% SWMM input data options -- Compatibility type
    enum, bind(c)
        enumerator :: API_SWMM5 = 0
        enumerator :: API_SWMM3     ! 1
        enumerator :: API_SWMM4     ! 2
    end enum

    !% SWMM input data options -- Normal flow type
    enum, bind(c)
        enumerator :: API_SLOPE = 0
        enumerator :: API_FROUDE        ! 1
        enumerator :: API_BOTH          ! 2
    end enum

    !% SWMM input data options -- Inertial damping type
    enum, bind(c)
        enumerator :: API_NO_DAMPING = 0
        enumerator :: API_PARTIAL_DAMPING   ! 1
        enumerator :: API_FULL_DAMPING      ! 2
    end enum
    
    !% SWMM input data options -- Surcharge Method
    enum, bind(c)
        enumerator :: API_EXTRAN = 0
        enumerator :: API_SLOT      ! 1
    end enum

    !% SWMM input data options -- Inflow type
    enum, bind(c)
        enumerator :: API_EXTERNAL_INFLOW = 0
        enumerator :: API_DRY_WEATHER_INFLOW    ! 1
        enumerator :: API_WET_WEATHER_INFLOW    ! 2
        enumerator :: API_GROUNDWATER_INFLOW    ! 3
        enumerator :: API_RDII_INFLOW           ! 4
        enumerator :: API_FLOW_INFLOW           ! 5
        enumerator :: API_CONCEN_INFLOW         ! 6
        enumerator :: API_MASS_INFLOW           ! 7
    end enum

    !% SWMM PATTERN TYPES ($API_DIR/src/enums.h -> PatternType)
    enum, bind(c)
        enumerator :: API_MONTHLY_PATTERN = 0
        enumerator :: API_DAILY_PATTERN         ! 1
        enumerator :: API_HOURLY_PATTERN        ! 2
        enumerator :: API_WEEKEND_PATTERN       ! 3
    end enum

    !% SWMM outfall types
    enum, bind(c)
        enumerator :: API_FREE_OUTFALL = 0
        enumerator :: API_NORMAL_OUTFALL        ! 1
        enumerator :: API_FIXED_OUTFALL         ! 2
        enumerator :: API_TIDAL_OUTFALL         ! 3
        enumerator :: API_TIMESERIES_OUTFALL    ! 4
    end enum

    !% SWMM input data options -- Storage type
    enum, bind(c)
        enumerator :: API_TABULAR = 0
        enumerator :: API_FUNCTIONAL        ! 1
    end enum

    !% SWMM input data options -- Reactor Type
    enum, bind(c)
        enumerator :: API_CSTR = 0
        enumerator :: API_PLUG      ! 1
    end enum
        
    !% SWMM input data options -- Treatment Type
    enum, bind(c)
        enumerator :: API_REMOVAL = 0
        enumerator :: API_CONCEN        ! 1
    end enum
        
    !% SWMM input data options -- Divider Type
    enum, bind(c)
        enumerator :: API_CUTOFF_DIVIDER = 0
        enumerator :: API_TABULAR_DIVIDER       ! 1
        enumerator :: API_WEIR_DIVIDER          ! 2
        enumerator :: API_OVERFLOW_DIVIDER      ! 3
    end enum

    !% SWMM Pump Curve types ($API_DIR/src/enums.h) -> PumpCurveType
    enum, bind(c)
        enumerator :: API_TYPE1_PUMP = 0
        enumerator :: API_TYPE2_PUMP        ! 1
        enumerator :: API_TYPE3_PUMP        ! 2
        enumerator :: API_TYPE4_PUMP        ! 3
        enumerator :: API_IDEAL_PUMP        ! 4
    end enum

    !% SWMM Orifice types ($API_DIR/src/enums.h) -> OrificeType
    enum, bind(c)
        enumerator :: API_SIDE_ORIFICE = 0
        enumerator :: API_BOTTOM_ORIFICE    ! 1
    end enum

    !% SWMM Weir types ($API_DIR/src/enums.h) -> WeirType
    enum, bind(c)
        enumerator :: API_TRANSVERSE_WEIR = 0
        enumerator :: API_SIDEFLOW_WEIR     ! 1
        enumerator :: API_VNOTCH_WEIR       ! 2
        enumerator :: API_TRAPEZOIDAL_WEIR  ! 3
        enumerator :: API_ROADWAY_WEIR      ! 4
    end enum
 
    !% SWMM Curve types ($API_DIR/src/enums.h -> ObjectType)
    enum, bind(c)
        enumerator :: API_STORAGE_CURVE = 0         !% 0 surf. area v. depth for storage node
        enumerator :: API_DIVERSION_CURVE           !% 1 diverted flow v. inflow for divider node
        enumerator :: API_TIDAL_CURVE               !% 2 water elev. v. hour of day for outfall
        enumerator :: API_RATING_CURVE              !% 3 flow rate v. head for outlet link
        enumerator :: API_CONTROL_CURVE             !% 4 control setting v. controller variable
        enumerator :: API_SHAPE_CURVE               !% 5 width v. depth for custom x-section
        enumerator :: API_WEIR_CURVE                !% 6 discharge coeff. v. head for weir
        enumerator :: API_PUMP1_CURVE               !% 7 flow v. wet well volume for pump
        enumerator :: API_PUMP2_CURVE               !% 8 flow v. depth for pump (discrete)
        enumerator :: API_PUMP3_CURVE               !% 9 flow v. head for pump (continuous)
        enumerator :: API_PUMP4_CURVE               !% 10 flow v. depth for pump (continuous)
    end enum       

    !% brh 20211208 commented this and replaced with completel list
    ! !% SWMM Computed node quantities ($API_DIR/src/enums.h) -> MAX_NODE_RESULTS
    ! !% HACK: These keys are used for SWMM5 outlet type as well
    ! enum, bind(c)
    !     enumerator :: API_NODE_DEPTH = 0
    !     enumerator :: API_NODE_HEAD
    ! end enum

    ! !% SWMM Table types ($API_DIR/src/enums.h -> ObjectType)
    ! enum, bind(c)
    !     enumerator :: API_TIMEPATTERN = 6
    !     enumerator :: API_CURVE
    !     enumerator :: API_TSERIES
    ! end enum

    !% brh20211208s -- this looks like obsolete code
    !% API VARS
    !enum, bind(c)
    !    enumerator :: API_NODES_WITH_EXTINFLOW = 1000
    !    enumerator :: API_NODES_WITH_DWFINFLOW
    !end enum
    !% brh20211208e

    !% API Node Attributes
    !% the "nodef" is used to prevent clashes with uppercase API_ constants
    !% Note that api.h must be changed if any enum changed!
    enum, bind(c)
        enumerator :: api_nodef_start = 1
        enumerator :: api_nodef_ID    !2
        enumerator :: api_nodef_type   ! 3
        enumerator :: api_nodef_outfall_type  !4
        enumerator :: api_nodef_invertElev    !5
        enumerator :: api_nodef_initDepth     !6
        enumerator :: api_nodef_StorageConstant   !7
        enumerator :: api_nodef_StorageCoeff      !8
        enumerator :: api_nodef_StorageExponent   !9
        enumerator :: api_nodef_StorageCurveID    !10
        enumerator :: api_nodef_extInflow_tSeries !11
        enumerator :: api_nodef_extInflow_tSeries_x1  !12
        enumerator :: api_nodef_extInflow_tSeries_x2  !13
        enumerator :: api_nodef_extInflow_basePat_idx     !14
        !% brh20211707s 
        !enumerator :: api_node_extInflow_baseline    !moved
        enumerator :: api_nodef_extInflow_basePat_type   !15
        enumerator :: api_nodef_extInflow_baseline    !16
        !% brh202211207e
        enumerator :: api_nodef_extInflow_sFactor        !17
        enumerator :: api_nodef_has_extInflow            !18
        enumerator :: api_nodef_dwfInflow_monthly_pattern  !19
        enumerator :: api_nodef_dwfInflow_daily_pattern    !20
        enumerator :: api_nodef_dwfInflow_hourly_pattern   !21
        enumerator :: api_nodef_dwfInflow_weekend_pattern  !22
        enumerator :: api_nodef_dwfInflow_avgvalue         !23
        enumerator :: api_nodef_has_dwfInflow   !24
        !% brh20211207s 
        enumerator :: api_nodef_newDepth        !25
        !% brh20211207e 
        enumerator :: api_nodef_fullDepth       !26  x!24
        enumerator :: api_nodef_inflow          !27  x!25
        enumerator :: api_nodef_volume          !28  x!26
        enumerator :: api_nodef_overflow        !29  x!27
        enumerator :: api_nodef_rptFlag     !30
        !% brh 20220521
        enumerator :: api_nodef_hasFlapGate  !31
        enumerator :: api_nodef_end         !32

        enumerator :: api_linkf_start     ! 33
        enumerator :: api_linkf_ID        ! 34
        enumerator :: api_linkf_subIndex  ! 35
        enumerator :: api_linkf_direction ! 36
        enumerator :: api_linkf_node1     ! 37
        enumerator :: api_linkf_node2     ! 38
        enumerator :: api_linkf_offset1   ! 39
        enumerator :: api_linkf_offset2   ! 40
        enumerator :: api_linkf_q0        ! 41
        enumerator :: api_linkf_flow      ! 42
        enumerator :: api_linkf_depth     ! 43
        enumerator :: api_linkf_volume    ! 44
        enumerator :: api_linkf_froude    ! 45
        enumerator :: api_linkf_setting   ! 46
        enumerator :: api_linkf_left_slope       ! 47
        enumerator :: api_linkf_right_slope      ! 48  
        enumerator :: api_linkf_weir_end_contractions ! 49
        enumerator :: api_linkf_weir_side_slope   ! 50
        enumerator :: api_linkf_curveid      ! 51
        enumerator :: api_linkf_discharge_coeff1  ! 52
        enumerator :: api_linkf_discharge_coeff2  ! 53
        enumerator :: api_linkf_initSetting       ! 54
        enumerator :: api_linkf_yOn               ! 55
        enumerator :: api_linkf_yOff              ! 56
        enumerator :: api_linkf_conduit_roughness ! 57
        enumerator :: api_linkf_conduit_length    ! 58
        enumerator :: api_linkf_rptFlag           ! 59
        !% brh20211207s
        enumerator :: api_linkf_commonBreak      ! 60  ! must be end of common types
        !% brh20211207e
        ! --- special elements attributes
        enumerator :: api_linkf_type         ! 61
        enumerator :: api_linkf_sub_type     ! 62
        enumerator :: api_linkf_typeBreak    ! 63  ! must be end of special types
        ! --- xsect attributes 
        enumerator :: api_linkf_xsect_type   ! 64
        enumerator :: api_linkf_geometry     ! 65
        enumerator :: api_linkf_xsect_wMax   ! 66
        enumerator :: api_linkf_xsect_yBot   ! 67
        enumerator :: api_linkf_xsect_yFull  ! 68
        enumerator :: api_linkf_transectidx  ! 69
        enumerator :: api_linkf_end     ! 70  ! must be end of the xsect attributes

        enumerator :: api_transectf_start   ! 71
        enumerator :: api_transectf_ID      ! 72
        enumerator :: api_transectf_yFull   ! 73
        enumerator :: api_transectf_aFull   ! 74
        enumerator :: api_transectf_rFull   ! 75
        enumerator :: api_transectf_wMax    ! 76
        enumerator :: api_transectf_ywMax   ! 77
        enumerator :: api_transectf_sMax    ! 78
        enumerator :: api_transectf_aMax    ! 79
        enumerator :: api_transectf_lengthFactor    ! 80
        enumerator :: api_transectf_roughness       ! 81
        enumerator :: api_transectf_end             ! 82

        enumerator :: api_keyslastplusone  ! 83
    end enum

    ! !% API link attributes
    ! enum, bind(c)
    !     enumerator :: api_linkf_ID = 1
    !     enumerator :: api_linkf_subIndex  ! 2
    !     enumerator :: api_linkf_direction ! 3
    !     enumerator :: api_linkf_node1     ! 4
    !     enumerator :: api_linkf_node2     ! 5
    !     enumerator :: api_linkf_offset1   ! 6
    !     enumerator :: api_linkf_offset2   ! 7
    !     enumerator :: api_linkf_q0        ! 8
    !     enumerator :: api_linkf_flow      ! 9
    !     enumerator :: api_linkf_depth     ! 10
    !     enumerator :: api_linkf_volume    ! 11
    !     enumerator :: api_linkf_froude    ! 12
    !     enumerator :: api_linkf_setting   ! 13
    !     enumerator :: api_linkf_left_slope       ! 14
    !     enumerator :: api_linkf_right_slope      ! 15  
    !     enumerator :: api_linkf_weir_end_contractions ! 16
    !     enumerator :: api_linkf_weir_side_slope   ! 17
    !     enumerator :: api_linkf_curveid      ! 18
    !     enumerator :: api_linkf_discharge_coeff1  ! 19
    !     enumerator :: api_linkf_discharge_coeff2  ! 20
    !     enumerator :: api_linkf_initSetting       ! 21
    !     enumerator :: api_linkf_yOn               ! 22
    !     enumerator :: api_linkf_yOff              ! 23
    !     enumerator :: api_linkf_conduit_roughness ! 24
    !     enumerator :: api_linkf_conduit_length    ! 25
    !     !% brh20211207s
    !     enumerator :: api_linkf_commonBreak      ! 26  ! must be end of common types
    !     !% brh20211207e
    !     ! --- special elements attributes
    !     enumerator :: api_linkf_type         ! 27
    !     enumerator :: api_linkf_sub_type     ! 28
    !     enumerator :: api_linkf_typeBreak    ! 29  ! must be end of special types
    !     ! --- xsect attributes 
    !     enumerator :: api_linkf_xsect_type   ! 30
    !     enumerator :: api_linkf_geometry     ! 31
    !     enumerator :: api_linkf_xsect_wMax   ! 32
    !     enumerator :: api_linkf_xsect_yBot   ! 33
    !     enumerator :: api_linkf_xsect_yFull  ! 34
    !     enumerator :: api_linkf_transectidx  ! 35
    !     enumerator :: api_linkf_totalEnd     ! 36  ! must be end of the enum

    ! end enum


    !% API table attributes
    enum, bind(c)
        enumerator :: api_table_ID = 1
        enumerator :: api_table_type
        enumerator :: api_table_refers_to
        enumerator :: api_table_end
    end enum

    !% API Link Output attributes
    enum, bind(c)
        enumerator :: api_output_link_depth = 0
        enumerator :: api_output_link_flow
        enumerator :: api_output_link_volume
        enumerator :: api_output_link_direction
    end enum

    !% API Node Output attributes
    enum, bind(c)
        enumerator :: api_output_node_depth = 0
        enumerator :: api_output_node_volume
        enumerator :: api_output_node_latflow
        enumerator :: api_output_node_inflow
    end enum

    ! Datetime resolution types
    enum, bind(c)
        enumerator :: api_monthly = 1
        enumerator :: api_weekend
        enumerator :: api_daily
        enumerator :: api_hourly
    end enum

    contains   
    !%==========================================================================
    !%
        subroutine define_apikeys_reverse ()
            !%------------------------------------------------------------------
            !% Description:
            !% creates the reverseKey_api global that provide the string 
            !% name for the keys defined in define_keys
            !%------------------------------------------------------------------
            !% Declarations
            !%------------------------------------------------------------------
            !% Preliminaries
            !%------------------------------------------------------------------
            !% allocate the global space for the reverse keys
            allocate(reverseKey_api(api_keyslastplusone))
    
            !% define the reverse keys
            reverseKey_api(api_nodef_start) = 'api_nodef_start'
            reverseKey_api(api_nodef_ID) = 'api_nodef_ID'
            reverseKey_api(api_nodef_type) = 'api_nodef_type'
            reverseKey_api(api_nodef_outfall_type) = 'api_nodef_outfall_type'
            reverseKey_api(api_nodef_invertElev) = 'api_nodef_invertElev'
            reverseKey_api(api_nodef_initDepth) = 'api_nodef_initDepth'
            reverseKey_api(api_nodef_StorageConstant) = 'api_nodef_StorageConstant'
            reverseKey_api(api_nodef_StorageCoeff) = 'api_nodef_StorageCoeff'
            reverseKey_api(api_nodef_StorageExponent) = 'api_nodef_StorageExponent'
            reverseKey_api( api_nodef_StorageCurveID) = ' api_nodef_StorageCurveID'
            reverseKey_api(api_nodef_extInflow_tSeries) = 'api_nodef_extInflow_tSeries'
            reverseKey_api(api_nodef_extInflow_tSeries_x1) = 'api_nodef_extInflow_tSeries_x1'
            reverseKey_api(api_nodef_extInflow_tSeries_x2) = 'api_nodef_extInflow_tSeries_x2'
            reverseKey_api(api_nodef_extInflow_basePat_idx) = 'api_nodef_extInflow_basePat_idx'
            reverseKey_api(api_nodef_extInflow_basePat_type) = 'api_nodef_extInflow_basePat_type'
            reverseKey_api(api_nodef_extInflow_baseline) = 'api_nodef_extInflow_baseline'
            reverseKey_api(api_nodef_extInflow_sFactor) = 'api_nodef_extInflow_sFactor'
            reverseKey_api(api_nodef_has_extInflow) = 'api_nodef_has_extInflow'
            reverseKey_api(api_nodef_dwfInflow_monthly_pattern) = 'api_nodef_dwfInflow_monthly_pattern'
            reverseKey_api(api_nodef_dwfInflow_daily_pattern) = 'api_nodef_dwfInflow_daily_pattern'
            reverseKey_api(api_nodef_dwfInflow_hourly_pattern) = 'api_nodef_dwfInflow_hourly_pattern'
            reverseKey_api(api_nodef_dwfInflow_weekend_pattern) = 'api_nodef_dwfInflow_weekend_pattern'
            reverseKey_api(api_nodef_dwfInflow_avgvalue) = 'api_nodef_dwfInflow_avgvalue'
            reverseKey_api(api_nodef_has_dwfInflow) = 'api_nodef_has_dwfInflow'
            reverseKey_api(api_nodef_newDepth) = 'api_nodef_newDepth'
            reverseKey_api(api_nodef_fullDepth) = 'api_nodef_fullDepth'
            reverseKey_api(api_nodef_inflow) = 'api_nodef_inflow'
            reverseKey_api(api_nodef_volume) = 'api_nodef_volume'
            reverseKey_api(api_nodef_overflow) = 'api_nodef_overflow'
            reverseKey_api(api_nodef_rptFlag) = 'api_nodef_rptFlag'
            reverseKey_api(api_nodef_end) = 'api_nodef_end'
            reverseKey_api(api_linkf_start) = 'api_linkf_start'
            reverseKey_api(api_linkf_ID) = 'api_linkf_ID'
            reverseKey_api(api_linkf_subIndex) = 'api_linkf_subIndex'
            reverseKey_api(api_linkf_direction) = 'api_linkf_direction'
            reverseKey_api(api_linkf_node1) = 'api_linkf_node1'
            reverseKey_api(api_linkf_node2) = 'api_linkf_node2'
            reverseKey_api(api_linkf_offset1) = 'api_linkf_offset1'
            reverseKey_api(api_linkf_offset2) = 'api_linkf_offset2'
            reverseKey_api(api_linkf_q0) = 'api_linkf_q0'
            reverseKey_api(api_linkf_flow) = 'api_linkf_flow'
            reverseKey_api(api_linkf_depth) = 'api_linkf_depth'
            reverseKey_api(api_linkf_volume) = 'api_linkf_volume'
            reverseKey_api(api_linkf_froude) = 'api_linkf_froude'
            reverseKey_api(api_linkf_setting) = 'api_linkf_setting'
            reverseKey_api(api_linkf_left_slope) = 'api_linkf_left_slope'
            reverseKey_api(api_linkf_right_slope) = 'api_linkf_right_slope'
            reverseKey_api(api_linkf_weir_end_contractions) = 'api_linkf_weir_end_contractions'
            reverseKey_api(api_linkf_weir_side_slope) = 'api_linkf_weir_side_slope'
            reverseKey_api(api_linkf_curveid) = 'api_linkf_curveid'
            reverseKey_api(api_linkf_discharge_coeff1) = 'api_linkf_discharge_coeff1'
            reverseKey_api(api_linkf_discharge_coeff2) = 'api_linkf_discharge_coeff2'
            reverseKey_api(api_linkf_initSetting) = 'api_linkf_initSetting'
            reverseKey_api(api_linkf_yOn) = 'api_linkf_yOn'
            reverseKey_api(api_linkf_yOff) = 'api_linkf_yOff'
            reverseKey_api(api_linkf_conduit_roughness) = 'api_linkf_conduit_roughness'
            reverseKey_api(api_linkf_conduit_length) = 'api_linkf_conduit_length'
            reverseKey_api(api_linkf_rptFlag) = 'api_linkf_rptFlag'
            reverseKey_api(api_linkf_commonBreak) = 'api_linkf_commonBreak'
            reverseKey_api(api_linkf_type) = 'api_linkf_type'
            reverseKey_api(api_linkf_sub_type ) = 'api_linkf_sub_type'
            reverseKey_api(api_linkf_typeBreak ) = 'api_linkf_typeBreak'
            reverseKey_api(api_linkf_xsect_type ) = 'api_linkf_xsect_type' 
            reverseKey_api(api_linkf_geometry) = 'api_linkf_geometry'
            reverseKey_api(api_linkf_xsect_wMax) = 'api_linkf_xsect_wMax'
            reverseKey_api(api_linkf_xsect_yBot) = 'api_linkf_xsect_yBot'
            reverseKey_api(api_linkf_xsect_yFull) = 'api_linkf_xsect_yFull'
            reverseKey_api(api_linkf_transectidx) = 'api_linkf_transectidx'
            reverseKey_api(api_linkf_end) = 'api_linkf_end'

            reverseKey_api(api_transectf_ID) = "api_transectf_ID"
            reverseKey_api(api_transectf_yFull) = "api_transectf_yFull"
            reverseKey_api(api_transectf_aFull) = "api_transectf_aFull"
            reverseKey_api(api_transectf_rFull) = "api_transectf_rFull"
            reverseKey_api(api_transectf_wMax) = "api_transectf_wMax"
            reverseKey_api(api_transectf_ywMax) = "api_transectf_ywMax"
            reverseKey_api(api_transectf_sMax) = "api_transectf_sMax"
            reverseKey_api(api_transectf_aMax) = "api_transectf_aMax"
            reverseKey_api(api_transectf_lengthFactor) = "api_transectf_lengthFactor"
            reverseKey_api(api_transectf_roughness) = "api_transectf_roughness"
            reverseKey_api(api_transectf_end) = "api_transectf_end"

            reverseKey_api(api_keyslastplusone) = 'api_keyslastplusone'
        !%------------------------------------------------------------------
        !% Closing
    end subroutine define_apikeys_reverse
!%    
!%==========================================================================
!%==========================================================================
!%           
end module define_api_keys