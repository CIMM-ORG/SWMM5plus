! module define_keys
!
! Provides relationship between integers and keys used for different
! data types.
!
! For example, the elem2I(:,:) array has column ei_elem_type that provides
! the index to the element type for each element. The possible values are
! provided below as eChannel, ePipe, etc.
!
!==========================================================================
module define_keys

    implicit none

    !% BRH 20210608 New data keyes needed for identifiying time loops
    enum, bind(c)
        enumerator :: hydrology = 1     !% indicates hydrology loop
        enumerator :: hydraulics        !% indicates hydraulics loop
        enumerator :: ALLtm             !% indictes all time marching types
        enumerator :: ETM               !% Explicit time march
        enumerator :: ETM_AC            !% Explicit time march and AC
        enumerator :: AC                !% AC
        !% different link and their geometry types (HACK: probably could be consolidated to element types)
        enumerator :: lchannel          !% channel link  
        enumerator :: lpipe             !% pipe link
        enumerator :: lweir             !% weir link
        enumerator :: lTransverseWeir   !% transverse weir link
        enumerator :: lSideFlowWeir     !% sideflow weir link
        enumerator :: lRoadWayWeir      !% roadway weir link
        enumerator :: lVnotchWeir       !% vnotch weir link
        enumerator :: lTrapezoidalWeir  !% trapezoidal weir link
        enumerator :: lOrifice          !% orifice link
        enumerator :: lBottomOrifice    !% bottom orifice link
        enumerator :: lSideOrifice      !% side orifice link
        enumerator :: lPump             !% pump link
        enumerator :: lType1Pump        !% type 1 pump link
        enumerator :: lType2Pump        !% type 2 pump link
        enumerator :: lType3Pump        !% type 3 pump link
        enumerator :: lType4Pump        !% type 4 pump link
        enumerator :: lRectangular      !% rectangular link geometry
        enumerator :: lParabolic        !% parabolic link geometry
        enumerator :: lTrapezoidal      !% trapezoidal link geometry
        enumerator :: lTriangular       !% triangular link geometry
        enumerator :: lWidthDepth       !% widthdepth link geometry
        enumerator :: lCircular         !% circular link geometry
        !% different link roughness types
        enumerator :: lManningsN        !% ManningsN roughness type
        enumerator :: lCD               !% drag coefficient roughness type
        !% data types for link lengths adjustments
        enumerator :: NoAdjust          !% no link length adjustment has done
        enumerator :: OneSideAdjust     !% one sided link length adjustment has done
        enumerator :: BothSideAdjust    !% both sided link length adjustment has done
        !% different node types (HACK: probably could be consolidated to element types)
        enumerator :: nJ2               !% junction node with 2 links
        enumerator :: nJm               !% junction node with multiple links
        enumerator :: nStorage          !% stroage node
        enumerator :: nBCdn             !% downstream BC node
        enumerator :: nBCup             !% upstream BC node
        !% SWMM5+ elements types
        enumerator :: CC                !% conduit or channel element
        enumerator :: weir              !% weir element
        enumerator :: orifice           !% orifice element
        enumerator :: pump              !% pump element
        enumerator :: JM                !% junction main element
        enumerator :: JB                !% junction branch element
        enumerator :: storage           !% storage element
        enumerator :: manhole           !% manhole elemen (HACK: not sure if we need this)
        !% SWMM5+ CC geometry types
        enumerator :: rectangular       !% rectangular channel or conduit
        enumerator :: circular          !% circular channel or conduit
        enumerator :: triangular        !% triangular channel or conduit
        enumerator :: trapezoidal       !% trapezoidal channel or conduit
        enumerator :: parabolic         !% parabolic channel or conduit
        !% SWMM5+ CC roughness type
        enumerator :: ManningsN         ! ID for mannings n for roughness_type
        enumerator :: CD                ! ID for using drag coefficient for roughness_type
        !% SWMM5+ element types based on time marching
        enumerator :: diagnostic        !% diagnostic element
        enumerator :: time_march        !% indicates a time marched 
        !% SWMM5+ special element types
        enumerator :: transverse_weir   !% transverse weir type 
        enumerator :: side_flow         !% sideflow weir type
        enumerator :: roadway_weir      !% roadway weir type  
        enumerator :: vnotch_weir       !% vnotch weir type
        enumerator :: trapezoidal_weir  !% trapezoidal weir type
        enumerator :: bottom_orifice    !% bottom orifice type
        enumerator :: side_orifice      !% side orifice type
        enumerator :: type1_Pump        !% type 1 pump type
        enumerator :: type2_Pump        !% type 2 pump type
        enumerator :: type3_Pump        !% type 3 pump type
        enumerator :: type4_Pump        !% type 4 pump type
        !% BC face types
        enumerator :: BCup              !% Up BC face
        enumerator :: BCdn              !% Dn BC face
        enumerator :: critical          !% critical depth BC designator
        enumerator :: depth             !% depth (subcritical) BC designator
        enumerator :: flowrate          !% flowrate BC designato
        !% face jump types
        enumerator :: jump_none             !% type of hydraulic jump
        enumerator :: jump_from_upstream    !% type of hydraulic jump
        enumerator :: jump_from_downstream  !% type of hydraulic jump
        !% type of momentum source
        enumerator :: T00               !% type of momentum source
        enumerator :: T10               !% type of momentum source
        enumerator :: T20               !% type of momentum source
        !% MISC keys
        enumerator :: doesnotexist      !% used in various places to indicate something doesn't exist
        enumerator :: vshape            !% type of adjustment
        enumerator :: vshape_surcharge_only     !% type of adjustment
        enumerator :: FroudeNumber      !% data types for limiter BC approach
    end enum


    ! data types for bcdata
    integer, parameter :: bc_updn_downstream = 1
    integer, parameter :: bc_updn_upstream = 0
    integer, parameter :: bc_category_elevation = 0
    integer, parameter :: bc_category_inflowrate = 1

    
    ! data types for Partitioing Algorithm type (setting%Partitioning%PartitioningMethod)
    enum, bind(c)
        enumerator :: Default = 1
        enumerator :: BQuick
        enumerator :: Random
        enumerator :: BLink
    end enum

    !% rm 20210610 brh because of issues with ALLtm key
    ! ! data types for solver (setting%Solver%SolverSelect)
    ! enum, bind(c)
    !     enumerator :: ETM = 1
    !     enumerator :: ETM_AC
    !     enumerator :: AC
    ! end enum

    ! SWMM objects
    enum, bind(c)
        enumerator :: SWMM_NODE = 2
        enumerator :: SWMM_LINK
    end enum

    ! SWMM Table types
    enum, bind(c)
        enumerator :: SWMM_TIMEPATTERN = 6
        enumerator :: SWMM_CURVES
        enumerator :: SWMM_TSERIES
    end enum

    ! API VARS
    enum, bind(c)
        enumerator :: API_NODES_WITH_EXTINFLOW = 1000
        enumerator :: API_NODES_WITH_DWFINFLOW
    end enum

    ! SWMM XSECT_TYPES
    enum, bind(c)
        enumerator :: SWMM_RECT_CLOSED = 3
        enumerator :: SWMM_RECT_OPEN
        enumerator :: SWMM_TRAPEZOIDAL
        enumerator :: SWMM_TRIANGULAR
        enumerator :: SWMM_PARABOLIC
    end enum

    ! SWMM PATTERN TYPES
    enum, bind(c)
        enumerator :: SWMM_MONTHLY_PATTERN = 0
        enumerator :: SWMM_DAILY_PATTERN
        enumerator :: SWMM_HOURLY_PATTERN
        enumerator :: SWMM_WEEKEND_PATTERN
    end enum

    ! API Node Attributes
    enum, bind(c)
        enumerator :: node_ID = 1
        enumerator :: node_type
        enumerator :: node_invertElev
        enumerator :: node_initDepth
        enumerator :: node_extInflow_tSeries
        enumerator :: node_extInflow_basePat
        enumerator :: node_extInflow_baseline
        enumerator :: node_extInflow_sFactor
        enumerator :: node_has_extInflow
        enumerator :: node_dwfInflow_monthly_pattern
        enumerator :: node_dwfInflow_daily_pattern
        enumerator :: node_dwfInflow_hourly_pattern
        enumerator :: node_dwfInflow_weekend_pattern
        enumerator :: node_dwfInflow_avgvalue
        enumerator :: node_has_dwfInflow
        enumerator :: node_inflow
        enumerator :: node_volume
        enumerator :: node_overflow
    end enum
    integer, parameter :: num_node_attributes = node_overflow

    ! nodeInflow attributes
    enum, bind(c)
        enumerator :: nf_extInflow_tSeries
        enumerator :: nf_extInflow_basePat
        enumerator :: nf_extInflow_baseline
        enumerator :: nf_extInflow_sFactor
        enumerator :: nf_dwfInflow_avgvalue
        enumerator :: nf_dwfInflow_monthly_pattern
        enumerator :: nf_dwfInflow_daily_pattern
        enumerator :: nf_dwfInflow_hourly_pattern
        enumerator :: nf_dwfInflow_weekend_pattern
    end enum
    integer, parameter :: num_nodeInflow_attributes = nf_dwfInflow_weekend_pattern

    ! API link attributes
    enum, bind(c)
        enumerator :: link_ID = 1
        enumerator :: link_subIndex
        enumerator :: link_node1
        enumerator :: link_node2
        enumerator :: link_q0
        enumerator :: link_flow
        enumerator :: link_depth
        enumerator :: link_volume
        enumerator :: link_froude
        enumerator :: link_setting
        enumerator :: link_left_slope
        enumerator :: link_right_slope
        enumerator :: conduit_roughness
        enumerator :: conduit_length
        ! --- xsect attributes
        enumerator :: link_type
        enumerator :: link_xsect_type
        enumerator :: link_geometry
        enumerator :: link_xsect_wMax
        enumerator :: link_xsect_yBot
    end enum
    integer, parameter :: num_link_attributes = conduit_length
    integer, parameter :: num_link_xsect_attributes = link_xsect_yBot - num_link_attributes
    integer, parameter :: num_total_link_attributes = num_link_attributes + num_link_xsect_attributes

    ! Table types
    enum, bind(c)
        enumerator :: tseries_table = 1
        enumerator :: curve_table
        enumerator :: tinflow
    end enum
    
    ! Datetime resolution types
    enum, bind(c)
        enumerator :: monthly = 1
        enumerator :: daily
        enumerator :: hourly
        enumerator :: weekend
    end enum

end module define_keys
