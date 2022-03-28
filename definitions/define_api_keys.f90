module define_api_keys

    implicit none

    !% --------------------------------------------------------
    !% Interface keys
    !%
    !%   The following variables must have the same values of
    !%   the variables in the EPA-SWMM source code. The especific
    !%   C file that is associated to each enumerator is
    !%   written at the beginning of each 'enum' declaration.
    !%
    !%   The variables have been renamed with respect to the
    !%   EPA-SWMM source code to avoid name conflicts. Concretely,
    !%   the 'API_' substring was added to the name at the
    !%   the beginning of the API key
    !%
    !%   Variables written in caps are extracted from native
    !%   EPA-SWMM, whereas lower case vars are added to the
    !%   EPA-SWMM
    !% --------------------------------------------------------

    ! !% SWMM Objects ($API_DIR/src/enums.h)
    ! enum, bind(c)
    !     enumerator :: API_NODE = 2
    !     enumerator :: API_LINK
    ! end enum

    !% SWMM Objects ($API_DIR/src/enums.h)
    enum, bind(c)
        enumerator :: API_GAGE = 0      ! not used yet
        enumerator :: API_SUBCATCH
        enumerator :: API_NODE
        enumerator :: API_LINK
        enumerator :: API_POLLUT        ! not used yet
        enumerator :: API_LANDUSE       ! not used yet
        enumerator :: API_TIMEPATTERN   !
        enumerator :: API_CURVE
        enumerator :: API_TSERIES
        enumerator :: API_CONTROL
        enumerator :: API_AQUIFER
        enumerator :: API_UNITHYD
        enumerator :: API_SNOWMELT
        enumerator :: API_SHAPE
        enumerator :: API_LID
        enumerator :: API_MAX_OBJ_TYPES
    end enum

    !% SWMM Node types ($API_DIR/src/enums.h) -> NodeType
    enum, bind(c)
        enumerator :: API_JUNCTION = 0
        enumerator :: API_OUTFALL
        enumerator :: API_STORAGE
        enumerator :: API_DIVIDER
    end enum

    !% SWMM Link types ($API_DIR/src/enums.h) -> LinkType
    enum, bind(c)
        enumerator :: API_CONDUIT = 0
        enumerator :: API_PUMP
        enumerator :: API_ORIFICE
        enumerator :: API_WEIR
        enumerator :: API_OUTLET
    end enum

    !% SWMM File types ($API_DIR/src/enums.h)
    enum, bind(c)
        enumerator :: API_RAINFALL_FILE = 0 ! not used yet
        enumerator :: API_RUNOFF_FILE   ! not used yet
        enumerator :: API_HOTSTART_FILE ! not used yet
        enumerator :: API_RDII_FILE     ! not used yet
        enumerator :: API_INFLOWS_FILE  ! not used yet
        enumerator :: API_OUTFLOWS_FILE ! not used yet
    end enum

    !% SWMM File usage types
    enum, bind(c)
        enumerator :: API_NO_FILE = 0       ! not used yet
        enumerator :: API_SCRATCH_FILE  ! not used yet
        enumerator :: API_USE_FILE      ! not used yet
        enumerator :: API_SAVE_FILE     ! not used yet
    end enum

    !% SWMM Rain gage data types
    enum, bind(c)
        enumerator :: API_RAIN_TSERIES = 0  ! not used yet
        enumerator :: API_RAIN_FILE     ! not used yet
    end enum

    !% SWMM XSECT_TYPES ($API_DIR/src/enums.h -> XsectType)
    enum, bind(c)
        enumerator :: API_DUMMY = 0
        enumerator :: API_CIRCULAR
        enumerator :: API_FILLED_CIRCULAR
        enumerator :: API_RECT_CLOSED
        enumerator :: API_RECT_OPEN
        enumerator :: API_TRAPEZOIDAL
        enumerator :: API_TRIANGULAR
        enumerator :: API_PARABOLIC
        enumerator :: API_POWERFUNC
        enumerator :: API_RECT_TRIANG
        enumerator :: API_RECT_ROUND
        enumerator :: API_MOD_BASKET
        enumerator :: API_HORIZ_ELLIPSE
        enumerator :: API_VERT_ELLIPSE
        enumerator :: API_ARCH
        enumerator :: API_EGGSHAPED
        enumerator :: API_HORSESHOE
        enumerator :: API_GOTHIC
        enumerator :: API_CATENARY
        enumerator :: API_SEMIELLIPTICAL
        enumerator :: API_BASKETHANDLE
        enumerator :: API_SEMICIRCULAR
        enumerator :: API_IRREGULAR
        enumerator :: API_CUSTOM
        enumerator :: API_FORCE_MAIN
    end enum

    !% SWMM measurement unit types
    enum, bind(c)
        enumerator :: API_US = 0
        enumerator :: API_SI
    end enum

    !% SWMM Flow Units type
    enum, bind(c)
        enumerator :: API_CFS = 0
        enumerator :: API_GPM
        enumerator :: API_MGD
        enumerator :: API_CMS
        enumerator :: API_LPS
        enumerator :: API_MLD
    end enum

    !% SWMM Concentration units type
    enum, bind(c)
        enumerator :: API_MG = 0
        enumerator :: API_UG
        enumerator :: API_COUNT
    end enum

    !% SWMM quantities requiring unit conversions
    enum, bind(c)
        enumerator :: API_RAINFALL = 0
        enumerator :: API_RAINDEPTH
        enumerator :: API_EVAPRATE
        enumerator :: API_LENGTH
        enumerator :: API_LANDAREA
        enumerator :: API_VOLUME
        enumerator :: API_WINDSPEED
        enumerator :: API_TEMPERATURE
        enumerator :: API_MASS
        enumerator :: API_GWFLOW
        enumerator :: API_FLOW
    end enum

    !% SWMM computed subcatchment quantities
    enum, bind(c)
        enumerator :: API_SUBCATCH_RAINFALL = 0
        enumerator :: API_SUBCATCH_SNOWDEPTH
        enumerator :: API_SUBCATCH_EVAP
        enumerator :: API_SUBCATCH_INFIL
        enumerator :: API_SUBCATCH_RUNOFF
        enumerator :: API_SUBCATCH_GW_FLOW
        enumerator :: API_SUBCATCH_GW_ELEV
        enumerator :: API_SUBCATCH_SOIL_MOISTURE
        enumerator :: API_SUBCATCH_WASHOFF
    end enum

    !% SWMM node quantites
    enum, bind(c)
        enumerator :: API_NODE_DEPTH = 0
        enumerator :: API_NODE_HEAD
        enumerator :: API_NODE_VOLUME
        enumerator :: API_NODE_LATFLOW
        enumerator :: API_NODE_INFLOW
        enumerator :: API_NODE_OVERFLOW
        enumerator :: API_NODE_QUAL
        enumerator :: API_
    end enum

    !% SWMM Link result type
    enum, bind(c)
        enumerator :: API_LINK_FLOW = 0
        enumerator :: API_LINK_DEPTH
        enumerator :: API_LINK_VELOCITY
        enumerator :: API_LINK_VOLUME
        enumerator :: API_LINK_CAPACITY
        enumerator :: API_LINK_QUAL
    end enum

    !% SWMM system-wide flow quantities
    enum, bind(c)
        enumerator :: API_SYS_TEMPERATURE = 0
        enumerator :: API_SYS_RAINFALL
        enumerator :: API_SYS_SNOWDEPTH
        enumerator :: API_SYS_INFIL
        enumerator :: API_SYS_RUNOFF
        enumerator :: API_SYS_DWFLOW
        enumerator :: API_SYS_GWFLOW
        enumerator :: API_SYS_IIFLOW
        enumerator :: API_SYS_EXFLOW
        enumerator :: API_SYS_INFLOW
        enumerator :: API_SYS_FLOODING
        enumerator :: API_SYS_OUTFLOW
        enumerator :: API_SYS_STORAGE
        enumerator :: API_SYS_EVAP
        enumerator :: API_SYS_PET
    end enum

    !% SWMM conduit flow classifications
    enum, bind(c)
        enumerator :: API_DRY = 0
        enumerator :: API_UP_DRY
        enumerator :: API_DN_DRY
        enumerator :: API_SUBCRITICAL
        enumerator :: API_SUPCRITICAL
        enumerator :: API_UP_CRITICAL
        enumerator :: API_DN_CRITICAL
        enumerator :: API_MAX_FLOW_CLASSES
        enumerator :: API_UP_FULL
        enumerator :: API_DN_FULL
        enumerator :: API_ALL_FULL
    end enum

    !% SWMM  runoff flow categories
    enum, bind(c)
        enumerator :: API_RUNOFF_RAINFALL = 0
        enumerator :: API_RUNOFF_EVAP
        enumerator :: API_RUNOFF_INFIL
        enumerator :: API_RUNOFF_RUNOFF
        enumerator :: API_RUNOFF_DRAINS
        enumerator :: API_RUNOFF_RUNON
    end enum

    !% SWMM pollutant loading categories
    enum, bind(c)
        enumerator :: API_BUILDUP_LOAD = 0
        enumerator :: API_DEPOSITION_LOAD
        enumerator :: API_SWEEPING_LOAD
        enumerator :: API_BMP_REMOVAL_LOAD
        enumerator :: API_INFIL_LOAD
        enumerator :: API_RUNOFF_LOAD
        enumerator :: API_FINAL_LOAD
    end enum

    !% SWMM input data options - rainfall
    enum, bind(c)
        enumerator :: API_RAINFALL_INTENSITY = 0
        enumerator :: API_RAINFALL_VOLUME
        enumerator :: API_CUMULATIVE_RAINFALL
    end enum

    !% SWMM input data options -- temperature
    enum, bind(c)
        enumerator :: API_NO_TEMP = 0
        enumerator :: API_TSERIES_TEMP
        enumerator :: API_FILE_TEMP
    end enum

    !% SWMM input data options -- wind
    enum, bind(c)
        enumerator :: API_MONTHLY_WIND = 0
        enumerator :: API_FILE_WIND

    end enum

    !% SWMM input data options -- evaporation 
    enum, bind(c)
        enumerator :: API_CONSTANT_EVAP = 0
        enumerator :: API_MONTHLY_EVAP
        enumerator :: API_TIMESERIES_EVAP
        enumerator :: API_TEMPERATURE_EVAP
        enumerator :: API_FILE_EVAP
        enumerator :: API_RECOVERY
        enumerator :: API_DRYONLY
    end enum

    !% SWMM input data options -- normalizing
    enum, bind(c)
        enumerator :: API_PER_AREA = 0
        enumerator :: API_PER_CURB
    end enum

    !% SWMM input data options -- buildup type
    enum, bind(c)
        enumerator :: API_NO_BUILDUP = 0
        enumerator :: API_POWER_BUILDUP
        enumerator :: API_EXPON_BUILDUP
        enumerator :: API_SATUR_BUILDUP
        enumerator :: API_EXTERNAL_BUILDUP
    end enum

    !% SWMM input data options -- wash off type
    enum, bind(c)
        enumerator :: API_NO_WASHOFF = 0
        enumerator :: API_EXPON_WASHOFF
        enumerator :: API_RATING_WASHOFF
        enumerator :: API_EMC_WASHOFF
    end enum

    !% SWMM input data options -- SubArea type
    enum, bind(c)
        enumerator :: API_IMPERV0 = 0
        enumerator :: API_IMPERV1
        enumerator :: API_PERV
    end enum

    !% SWMM input data options -- Runoff Routing type
    enum, bind(c)
        enumerator :: API_TO_OUTLET = 0
        enumerator :: API_TO_IMPERV
        enumerator :: API_TO_PERV
    end enum

    !% SWMM input data options -- Routine model type
    enum, bind(c)
        enumerator :: API_NO_ROUTING = 0
        enumerator :: API_SF
        enumerator :: API_KW
        enumerator :: API_EKW
        enumerator :: API_DW
    end enum

    !% SWMM input data options --  Force main type
    enum, bind(c)
        enumerator :: API_H_W = 0
        enumerator :: API_D_W
    end enum

    !% SWMM input data options -- Offset Type
    enum, bind(c)
        enumerator :: API_DEPTH_OFFSET = 0
        enumerator :: API_ELEV_OFFSET
    end enum

    !% SWMM input data options -- Kinematic wave method type
    enum, bind(c)
        enumerator :: API_NORMAL = 0
        enumerator :: API_MODIFIED
    end enum

    !% SWMM input data options -- Compatibility type
    enum, bind(c)
        enumerator :: API_SWMM5 = 0
        enumerator :: API_SWMM3
        enumerator :: API_SWMM4
    end enum

    !% SWMM input data options -- Normal flow type
    enum, bind(c)
        enumerator :: API_SLOPE = 0
        enumerator :: API_FROUDE
        enumerator :: API_BOTH
    end enum

    !% SWMM input data options -- Inertial damping type
    enum, bind(c)
        enumerator :: API_NO_DAMPING = 0
        enumerator :: API_PARTIAL_DAMPING
        enumerator :: API_FULL_DAMPING
    end enum
    
    !% SWMM input data options -- Surcharge Method
    enum, bind(c)
        enumerator :: API_EXTRAN = 0
        enumerator :: API_SLOT
    end enum

    !% SWMM input data options -- Inflow type
    enum, bind(c)
        enumerator :: API_EXTERNAL_INFLOW = 0
        enumerator :: API_DRY_WEATHER_INFLOW
        enumerator :: API_WET_WEATHER_INFLOW
        enumerator :: API_GROUNDWATER_INFLOW
        enumerator :: API_RDII_INFLOW
        enumerator :: API_FLOW_INFLOW
        enumerator :: API_CONCEN_INFLOW
        enumerator :: API_MASS_INFLOW
    end enum

    !% SWMM PATTERN TYPES ($API_DIR/src/enums.h -> PatternType)
    enum, bind(c)
        enumerator :: API_MONTHLY_PATTERN = 0
        enumerator :: API_DAILY_PATTERN
        enumerator :: API_HOURLY_PATTERN
        enumerator :: API_WEEKEND_PATTERN
    end enum

    !% SWMM outfall types
    enum, bind(c)
        enumerator :: API_FREE_OUTFALL = 0
        enumerator :: API_NORMAL_OUTFALL
        enumerator :: API_FIXED_OUTFALL
        enumerator :: API_TIDAL_OUTFALL
        enumerator :: API_TIMESERIES_OUTFALL
    end enum

    !% SWMM input data options -- Storage type
    enum, bind(c)
        enumerator :: API_TABULAR = 0
        enumerator :: API_FUNCTIONAL
    end enum

    !% SWMM input data options -- Reactor Type
    enum, bind(c)
        enumerator :: API_CSTR = 0
        enumerator :: API_PLUG
    end enum
        
    !% SWMM input data options -- Treatment Type
    enum, bind(c)
        enumerator :: API_REMOVAL = 0
        enumerator :: API_CONCEN
    end enum
        
    !% SWMM input data options -- Divider Type
    enum, bind(c)
        enumerator :: API_CUTOFF_DIVIDER = 0
        enumerator :: API_TABULAR_DIVIDER
        enumerator :: API_WEIR_DIVIDER
        enumerator :: API_OVERFLOW_DIVIDER
    end enum

    !% SWMM Pump Curve types ($API_DIR/src/enums.h) -> PumpCurveType
    enum, bind(c)
        enumerator :: API_TYPE1_PUMP = 0
        enumerator :: API_TYPE2_PUMP
        enumerator :: API_TYPE3_PUMP
        enumerator :: API_TYPE4_PUMP
        enumerator :: API_IDEAL_PUMP
    end enum

    !% SWMM Orifice types ($API_DIR/src/enums.h) -> OrificeType
    enum, bind(c)
        enumerator :: API_SIDE_ORIFICE = 0
        enumerator :: API_BOTTOM_ORIFICE
    end enum

    !% SWMM Weir types ($API_DIR/src/enums.h) -> WeirType
    enum, bind(c)
        enumerator :: API_TRANSVERSE_WEIR = 0
        enumerator :: API_SIDEFLOW_WEIR
        enumerator :: API_VNOTCH_WEIR
        enumerator :: API_TRAPEZOIDAL_WEIR
        enumerator :: API_ROADWAY_WEIR
    end enum
 
    !% SWMM Curve types ($API_DIR/src/enums.h -> ObjectType)
    enum, bind(c)
        enumerator :: API_STORAGE_CURVE = 0         !% surf. area v. depth for storage node
        enumerator :: API_DIVERSION_CURVE           !% diverted flow v. inflow for divider node
        enumerator :: API_TIDAL_CURVE               !% water elev. v. hour of day for outfall
        enumerator :: API_RATING_CURVE              !% flow rate v. head for outlet link
        enumerator :: API_CONTROL_CURVE             !% control setting v. controller variable
        enumerator :: API_SHAPE_CURVE               !% width v. depth for custom x-section
        enumerator :: API_WEIR_CURVE                !% discharge coeff. v. head for weir
        enumerator :: API_PUMP1_CURVE               !% flow v. wet well volume for pump
        enumerator :: API_PUMP2_CURVE               !% flow v. depth for pump (discrete)
        enumerator :: API_PUMP3_CURVE               !% flow v. head for pump (continuous)
        enumerator :: API_PUMP4_CURVE               !% flow v. depth for pump (continuous)
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
    enum, bind(c)
        enumerator :: api_nodef_ID = 1
        enumerator :: api_nodef_type   ! 2
        enumerator :: api_nodef_outfall_type  !3
        enumerator :: api_nodef_invertElev    !4
        enumerator :: api_nodef_initDepth     !5
        enumerator :: api_nodef_StorageConstant   !6
        enumerator :: api_nodef_StorageCoeff      !7
        enumerator :: api_nodef_StorageExponent   !8
        enumerator :: api_nodef_StorageCurveID    !9
        enumerator :: api_nodef_extInflow_tSeries !10
        enumerator :: api_nodef_extInflow_tSeries_x1  !11
        enumerator :: api_nodef_extInflow_tSeries_x2  !12
        enumerator :: api_nodef_extInflow_basePat_idx     !13
        !% brh20211707s 
        !enumerator :: api_node_extInflow_baseline    !moved
        enumerator :: api_nodef_extInflow_basePat_type   !14
        enumerator :: api_nodef_extInflow_baseline    !15
        !% brh202211207e
        enumerator :: api_nodef_extInflow_sFactor        !16
        enumerator :: api_nodef_has_extInflow            !17
        enumerator :: api_nodef_dwfInflow_monthly_pattern  !18
        enumerator :: api_nodef_dwfInflow_daily_pattern    !19
        enumerator :: api_nodef_dwfInflow_hourly_pattern   !20
        enumerator :: api_nodef_dwfInflow_weekend_pattern  !21
        enumerator :: api_nodef_dwfInflow_avgvalue         !22
        enumerator :: api_nodef_has_dwfInflow   !23
        !% brh20211207s 
        enumerator :: api_nodef_newDepth        !24
        !% brh20211207e 
        enumerator :: api_nodef_fullDepth       !25  x!24
        enumerator :: api_nodef_inflow          !26  x!25
        enumerator :: api_nodef_volume          !27  x!26
        enumerator :: api_nodef_overflow        !28  x!27
        !% brh20211207s
        enumerator :: api_nodef_rptFlag         !29
        !% brh20211207e
    end enum

    !% API link attributes
    enum, bind(c)
        enumerator :: api_linkf_ID = 1
        enumerator :: api_linkf_subIndex  ! 2
        enumerator :: api_linkf_direction ! 3
        enumerator :: api_linkf_node1     ! 4
        enumerator :: api_linkf_node2     ! 5
        enumerator :: api_linkf_offset1   ! 6
        enumerator :: api_linkf_offset2   ! 7
        enumerator :: api_linkf_q0        ! 8
        enumerator :: api_linkf_flow      ! 9
        enumerator :: api_linkf_depth     ! 10
        enumerator :: api_linkf_volume    ! 11
        enumerator :: api_linkf_froude    ! 12
        enumerator :: api_linkf_setting   ! 13
        enumerator :: api_linkf_left_slope       ! 14
        enumerator :: api_linkf_right_slope      ! 15  
        enumerator :: api_linkf_weir_end_contractions ! 16
        enumerator :: api_linkf_weir_side_slope   ! 17
        enumerator :: api_linkf_curveid      ! 18
        enumerator :: api_linkf_discharge_coeff1  ! 19
        enumerator :: api_linkf_discharge_coeff2  ! 20
        enumerator :: api_linkf_conduit_roughness ! 21
        enumerator :: api_linkf_conduit_length    ! 22
        !% brh20211207s
        enumerator :: api_linkf_rptFlag      ! 23
        !% brh20211207e
        ! --- special elements attributes
        enumerator :: api_linkf_type         ! 24
        enumerator :: api_linkf_sub_type     ! 25
        ! --- xsect attributes 
        enumerator :: api_linkf_xsect_type   ! 26
        enumerator :: api_linkf_geometry     ! 27
        enumerator :: api_linkf_xsect_wMax   ! 28
        enumerator :: api_linkf_xsect_yBot   ! 29
        enumerator :: api_linkf_xsect_yFull  ! 30

    end enum

    !% API table attributes
    enum, bind(c)
        enumerator :: api_table_ID = 1
        enumerator :: api_table_type
        enumerator :: api_table_refers_to
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

end module define_api_keys