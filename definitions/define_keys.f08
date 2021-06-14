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

    !% BRH 20210608 New data keyes needed for identifiying time loops
    enum, bind(c)
        enumerator :: hydrology = 1 !% indicates hydrology loop
        enumerator :: hydraulics    !% indicates hydraulics loop
        enumerator :: ALLtm         !% indictes all time marching types
        enumerator :: ETM
        enumerator :: ETM_AC
        enumerator :: AC
    end enum

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
        enumerator :: eJunctionMain ! ID for a junction
        enumerator :: eJunctionBranch ! ID for a junction
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
    integer, parameter :: fMultiple = eJunctionBranch ! ID for moderation by separate up/dn element types
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
    integer, parameter :: lUnassigned = 998877  ! Can't use the nullValueI because it causes a circular condition
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

    ! default for edge and non-edge node
    integer, parameter :: EdgeNode    = 1 ! Edge node of a partition
    integer, parameter :: nonEdgeNode = 0 ! Upstream BC nodes are assigned to 1 element

    ! data types from initial depth type in links
    enum, bind(c)
        enumerator :: Uniform = 1
        enumerator :: LinearlyVarying
        enumerator :: ExponentialDecay
    end enum
    
    ! data types for Partitioing Algorithm type (setting%Partitioning%PartitioningMethod)
    enum, bind(c)
        enumerator :: Default = 1
        enumerator :: BQuick
        enumerator :: Random
        enumerator :: BLink
    end enum

    

    ! data types for link lengths adjustments
    enum, bind(c)
        enumerator :: NoAdjust = 1
        enumerator :: OneSideAdjust
        enumerator :: BothSideAdjust
    end enum

    ! data types for Momentum Source type (setting%Solver%MomentumSourceM)
    enum, bind(c)
        enumerator :: T00 = 1
        enumerator :: T10
        enumerator :: T20
    end enum

    !% rm 20210610 brh because of issues with ALLtm key
    ! ! data types for solver (setting%Solver%SolverSelect)
    ! enum, bind(c)
    !     enumerator :: ETM = 1
    !     enumerator :: ETM_AC
    !     enumerator :: AC
    ! end enum

    ! data types for adjust approach (setting%Adjust)
    enum, bind(c)
        enumerator :: vshape = 1
        enumerator :: vshape_surcharge_only
    end enum

    ! data types for limiter BC approach (setting%Limiter%BC%approach)
    enum, bind(c)
        enumerator :: FroudeNumber = 1
    end enum

end module define_keys
