! module data_keys
!
! Provides relationship between integers and keys used for different
! data types.
!
! For example, the elem2I(:,:) array has column ei_elem_type that provides
! the index to the element type for each element. The possible values are
! provided below as eChannel, ePipe, etc.
!
!==========================================================================
module data_keys

    implicit none

    ! data types for elemI(:,ei_meta_elem_type).
    enum, bind(c)
        enumerator :: eHQ2 = 1 ! ID for 2 face HQ meta-element
        enumerator :: eHQm ! ID for multi-face HQ meta-element
        enumerator :: eQonly ! ID for a Q only meta-element
        enumerator :: eHonly ! ID for a H only meta-element
        enumerator :: eNonHQ ! ID for a non HQ meta-element
    end enum

    ! data types for elemI(:,ei_elem_type). faceI(:,e_type_u), faceI(:,e_type_d)
    enum, bind(c)
        enumerator :: eChannel = 1 ! ID for an open channel element
        enumerator :: ePipe ! ID for an pipe
        enumerator :: eJunctionChannel ! ID for a junction
        enumerator :: eJunctionPipe ! ID for a junction
        enumerator :: eCulvert ! ID for a culvert in an open channel
        enumerator :: ePump ! ID for a pump
        enumerator :: eValve ! ID for a valve
        enumerator :: eOrifice ! ID for an orifice
        enumerator :: eWeir ! ID for a weir
        enumerator :: eStorage ! ID for a storage
        enumerator :: eBCup ! ID for face upstream BC
        enumerator :: eBCdn ! ID for face downstream BC
    end enum

    ! date types for elemI(:,ei_weir_elem_type)
    enum, bind(c)
        enumerator :: eTransverseWeir = 1 ! ID for rectangular transverse weir
        enumerator :: eSideFlowWeir ! ID for rectangular sideflow weir
        enumerator :: eRoadWayWeir ! ID for rectangular roadway weir
        enumerator :: eVnotchWeir ! ID for triangular v-notch weir
        enumerator :: eTrapezoidalWeir ! ID for trapezoidal weir
    end enum

    ! date types for elemI(:,ei_orif_elem_type)
    enum, bind(c)
        enumerator :: eBottomOrifice = 1 ! ID for bottom orifice
        enumerator :: eSideOrifice ! ID for side orifice
    end enum

    ! date types for elemI(:,ei_pump_elem_type)
    enum, bind(c)
        enumerator :: eType1Pump = 1 ! ID for Type 1 pump
        enumerator :: eType2Pump ! ID for Type 2 pump
        enumerator :: eType3Pump ! ID for Type 3 pump
        enumerator :: eType4Pump ! ID for Type 4 pump
    end enum

    ! data types for faceI(:,fi_type)
    integer, parameter :: fChannel = eChannel ! ID for open channel on both sides
    integer, parameter :: fPipe = ePipe ! ID for pipe on both sides
    integer, parameter :: fWeir = eWeir ! ID for pipe on both sides
    integer, parameter :: fOrifice = eOrifice
    integer, parameter :: fPump = ePump
    integer, parameter :: fMultiple = eJunctionChannel ! ID for moderation by separate up/dn element types
    integer, parameter :: fBCup = eBCup ! ID for face upstream BC
    integer, parameter :: fBCdn = eBCdn ! ID for face downstream BC

    ! date types for elemI(:,ei_geometry)
    enum, bind(c)
        enumerator :: eRectangular = 1 ! ID for rectangular chanel, weir
        enumerator :: eParabolic ! ID for parabolic channel
        enumerator :: eTrapezoidal ! ID for trapezoidal channel, weir
        enumerator :: eTriangular ! ID for triangular channel, weir
        enumerator :: eWidthDepth ! ID for general geometry by data pairs
        enumerator :: eCircular ! ID for circular pipe, orifice
    end enum

    ! data types for elemI(:,ei_roughness_type)
    enum, bind(c)
        enumerator :: eManningsN = 1 ! ID for mannings n for roughness_type
        enumerator :: eCD ! ID for using drag coefficient for roughness_type
    end enum

    ! data types for faceI(:,jump_type)
    enum, bind(c)
        enumerator :: jump_none = 0 ! ID for no jump
        enumerator :: jump_downstream ! ID for jump in nominal downstream direction
        enumerator :: jump_upstream ! ID for jump in nominal upstream direction
    end enum

    ! data types for nodeI(:,ni_node_type)
    enum, bind(c)
        enumerator :: nJ2 = 1 ! ID for junction with 2 links
        enumerator :: nJm ! ID for junction with multiple links
        enumerator :: nStorage ! ID for stroage units
        enumerator :: nBCdn ! ID for downstream BC
        enumerator :: nBCup ! iD for upstream BC
    end enum

    ! data types for nodeI(:,ni_assigned) for assignment to faces and links
    integer, parameter :: nUnassigned = 998877
    integer, parameter :: nAssigned = 1
    integer, parameter :: nDeferred = -1

    ! data types for linkI(:,li_link_type)
    ! note that these must correspond to element types
    integer, parameter :: lchannel = eChannel ! ID for link that is open channel
    integer, parameter :: lpipe = ePipe ! ID for link that is pipe
    integer, parameter :: lweir = eWeir ! ID for link that is weir
    integer, parameter :: lOrifice = eOrifice ! ID for link that is orifice
    integer, parameter :: lPump = ePump ! ID for link that is pump

    ! date types for linkI(:,li_weir_type)
    integer, parameter :: lTransverseWeir = eTransverseWeir ! ID for rectangular transverse weir
    integer, parameter :: lSideFlowWeir = eSideFlowWeir ! ID for rectangular sideflow weir
    integer, parameter :: lRoadWayWeir = eRoadWayWeir ! ID for rectangular roadway weir
    integer, parameter :: lVnotchWeir = eVnotchWeir ! ID for triangular v-notch weir
    integer, parameter :: lTrapezoidalWeir = eTrapezoidalWeir ! ID for trapezoidal weir

    ! date types for linkI(:,li_orif_type)
    integer, parameter :: lBottomOrifice = eBottomOrifice ! ID for bottom orifice
    integer, parameter :: lSideOrifice = eSideOrifice ! ID for side orifice

    ! date types for elemI(:,li_pump_type)
    integer, parameter :: lType1Pump = eType1Pump ! ID for Type 1 pump
    integer, parameter :: lType2Pump = eType2Pump ! ID for Type 2 pump
    integer, parameter :: lType3Pump = eType3Pump ! ID for Type 3 pump
    integer, parameter :: lType4Pump = eType4Pump ! ID for Type 4 pump

    ! data types for linkI(:,li_geometry) (must corresponde with ei_geometry)
    integer, parameter :: lRectangular = eRectangular ! ID for link that rectangular channel, weir
    integer, parameter :: lParabolic = eParabolic ! ID for parabolic channel
    integer, parameter :: lTrapezoidal = eTrapezoidal ! ID for trapezoidal channel, weir
    integer, parameter :: lTriangular = eTriangular ! ID for triangle channel, weir
    integer, parameter :: lWidthDepth = eWidthDepth ! ID for general geometry by data pairs
    integer, parameter :: lCircular = eCircular ! ID for circular pipe, orifice

    ! data types for linkII(:,li_roughness_type)
    integer, parameter :: lManningsN = eManningsN ! ID for mannings n for roughness_type
    integer, parameter :: lCD = eCD ! ID for using drag coefficient for roughness_type

    ! data types for linkI(:,li_assigned) for assignment to faces and links
    integer, parameter :: lUnassigned = 998877
    integer, parameter :: lAssigned = 1
    integer, parameter :: lDeferred = -1

    ! data types for bcdata
    integer, parameter :: bc_updn_downstream = 1
    integer, parameter :: bc_updn_upstream = 0
    integer, parameter :: bc_category_elevation = 0
    integer, parameter :: bc_category_inflowrate = 1

    ! default number of elements for different node types
    integer, parameter :: N_elem_nJ2 = 0 ! 2-link nodes are assigned to a single face
    integer, parameter :: N_elem_nJm = 7 ! M-link nodes are assigned a maximum of 7 elements
    integer, parameter :: N_elem_nStorage = 1 ! Storage nodes are assigned to 1 element
    integer, parameter :: N_elem_nBCdn = 1 ! Downstream BC nodes are assigned to 1 element
    integer, parameter :: N_elem_nBCup = 1 ! Upstream BC nodes are assigned to 1 element

    ! data types for Partitioing Algorithm type (setting%Partitioning%PartitioningMethod)
    enum, bind(c)
        enumerator :: Default = 1
        enumerator :: BIPquick
    end enum

    ! data types for Momentum Source type (setting%Solver%MomentumSourceM)
    enum, bind(c)
        enumerator :: T00 = 1
        enumerator :: T10
        enumerator :: T20
    end enum

    ! data types for solver (setting%Solver%SolverSelect)
    enum, bind(c)
        enumerator :: SVE = 1
        enumerator :: SVE_AC
        enumerator :: AC
    end enum

    ! data types for adjust approach (setting%Adjust)
    enum, bind(c)
        enumerator :: vshape = 1
        enumerator :: smoothall
    end enum

    ! data types for limiter BC approach (setting%Limiter%BC%approach)
    enum, bind(c)
        enumerator :: FroudeNumber = 1
    end enum

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
end module data_keys