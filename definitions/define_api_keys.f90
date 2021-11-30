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

    !% SWMM Objects ($API_DIR/src/objects.c)
    enum, bind(c)
        enumerator :: API_NODE = 2
        enumerator :: API_LINK
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

    !% SWMM Computed node quantities ($API_DIR/src/enums.h) -> MAX_NODE_RESULTS
    !% HACK: These keys are used for SWMM5 outlet type as well
    enum, bind(c)
        enumerator :: API_NODE_DEPTH = 0
        enumerator :: API_NODE_HEAD
    end enum

    !% SWMM Table types ($API_DIR/src/enums.h -> ObjectType)
    enum, bind(c)
        enumerator :: API_TIMEPATTERN = 6
        enumerator :: API_CURVE
        enumerator :: API_TSERIES
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

    !% SWMM XSECT_TYPES ($API_DIR/src/enums.h -> XsectType)
    enum, bind(c)
        enumerator :: API_CIRCULAR = 1
        enumerator :: API_FILLED_CIRCULAR
        enumerator :: API_RECT_CLOSED
        enumerator :: API_RECT_OPEN
        enumerator :: API_TRAPEZOIDAL
        enumerator :: API_TRIANGULAR
        enumerator :: API_PARABOLIC

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

    !% API VARS
    enum, bind(c)
        enumerator :: API_NODES_WITH_EXTINFLOW = 1000
        enumerator :: API_NODES_WITH_DWFINFLOW
    end enum

    !% API Node Attributes
    enum, bind(c)
        enumerator :: api_node_ID = 1
        enumerator :: api_node_type
        enumerator :: api_node_outfall_type
        enumerator :: api_node_invertElev
        enumerator :: api_node_initDepth
        enumerator :: api_node_StorageConstant
        enumerator :: api_node_StorageCoeff
        enumerator :: api_node_StorageExponent
        enumerator :: api_node_StorageCurveID
        enumerator :: api_node_extInflow_tSeries
        enumerator :: api_node_extInflow_tSeries_x1
        enumerator :: api_node_extInflow_tSeries_x2
        enumerator :: api_node_extInflow_basePat
        enumerator :: api_node_extInflow_baseline
        enumerator :: api_node_extInflow_basePat_type
        enumerator :: api_node_extInflow_sFactor
        enumerator :: api_node_has_extInflow
        enumerator :: api_node_dwfInflow_monthly_pattern
        enumerator :: api_node_dwfInflow_daily_pattern
        enumerator :: api_node_dwfInflow_hourly_pattern
        enumerator :: api_node_dwfInflow_weekend_pattern
        enumerator :: api_node_dwfInflow_avgvalue
        enumerator :: api_node_has_dwfInflow
        enumerator :: api_node_fullDepth
        enumerator :: api_node_inflow
        enumerator :: api_node_volume
        enumerator :: api_node_overflow
    end enum

    !% API link attributes
    enum, bind(c)
        enumerator :: api_link_ID = 1
        enumerator :: api_link_subIndex
        enumerator :: api_link_node1
        enumerator :: api_link_node2
        enumerator :: api_link_offset1
        enumerator :: api_link_offset2
        enumerator :: api_link_q0
        enumerator :: api_link_flow
        enumerator :: api_link_depth
        enumerator :: api_link_volume
        enumerator :: api_link_froude
        enumerator :: api_link_setting
        enumerator :: api_link_left_slope
        enumerator :: api_link_right_slope
        enumerator :: api_weir_end_contractions
        enumerator :: api_weir_side_slope
        enumerator :: api_link_curveid
        enumerator :: api_discharge_coeff1
        enumerator :: api_discharge_coeff2
        enumerator :: api_conduit_roughness
        enumerator :: api_conduit_length
        ! --- special elements attributes
        enumerator :: api_link_type
        enumerator :: api_weir_type
        enumerator :: api_orifice_type
        enumerator :: api_outlet_type
        enumerator :: api_pump_type
        ! --- xsect attributes
        enumerator :: api_link_xsect_type
        enumerator :: api_link_geometry
        enumerator :: api_link_xsect_wMax
        enumerator :: api_link_xsect_yBot
        enumerator :: api_link_xsect_yFull
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