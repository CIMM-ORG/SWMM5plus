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
    !% --------------------------------------------------------

    !% SWMM Objects ($API_DIR/src/objects.c)
    enum, bind(c)
        enumerator :: API_NODE = 2
        enumerator :: API_LINK
    end enum

    !% SWMM Table types ($API_DIR/src/enums.h -> ObjectType)
    enum, bind(c)
        enumerator :: API_TIMEPATTERN = 6
        enumerator :: API_CURVES
        enumerator :: API_TSERIES
    end enum

    !% SWMM XSECT_TYPES ($API_DIR/src/enums.h -> XsectType)
    enum, bind(c)
        enumerator :: API_RECT_CLOSED = 3
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
        enumerator :: api_node_invertElev
        enumerator :: api_node_initDepth
        enumerator :: api_node_extInflow_tSeries
        enumerator :: api_node_extInflow_basePat
        enumerator :: api_node_extInflow_baseline
        enumerator :: api_node_extInflow_sFactor
        enumerator :: api_node_has_extInflow
        enumerator :: api_node_dwfInflow_monthly_pattern
        enumerator :: api_node_dwfInflow_daily_pattern
        enumerator :: api_node_dwfInflow_hourly_pattern
        enumerator :: api_node_dwfInflow_weekend_pattern
        enumerator :: api_node_dwfInflow_avgvalue
        enumerator :: api_node_has_dwfInflow
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
        enumerator :: api_link_q0
        enumerator :: api_link_flow
        enumerator :: api_link_depth
        enumerator :: api_link_volume
        enumerator :: api_link_froude
        enumerator :: api_link_setting
        enumerator :: api_link_left_slope
        enumerator :: api_link_right_slope
        enumerator :: api_conduit_roughness
        enumerator :: api_conduit_length
        ! --- xsect attributes
        enumerator :: api_link_type
        enumerator :: api_link_xsect_type
        enumerator :: api_link_geometry
        enumerator :: api_link_xsect_wMax
        enumerator :: api_link_xsect_yBot
    end enum

    ! Datetime resolution types
    enum, bind(c)
        enumerator :: api_monthly = 1
        enumerator :: api_daily
        enumerator :: api_hourly
        enumerator :: api_weekend
    end enum

end module define_api_keys