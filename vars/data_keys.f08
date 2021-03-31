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

    use globals

    implicit none

    ! data types for elemI(:,ei_meta_elem_type).
    integer, parameter :: eHQ2               = 1   ! ID for 2 face HQ meta-element
    integer, parameter :: eHQm               = 2   ! ID for multi-face HQ meta-element
    integer, parameter :: eQonly             = 3   ! ID for a Q only meta-element
    integer, parameter :: eHonly             = 4   ! ID for a H only meta-element
    integer, parameter :: eNonHQ             = 5   ! ID for a non HQ meta-element

    ! data types for elemI(:,ei_elem_type). faceI(:,e_type_u), faceI(:,e_type_d)
    integer, parameter :: eChannel          = 1   ! ID for an open channel element
    integer, parameter :: ePipe             = 2   ! ID for an pipe
    integer, parameter :: eJunctionChannel  = 3   ! ID for a junction
    integer, parameter :: eJunctionPipe     = 4   ! ID for a junction
    integer, parameter :: eCulvert          = 5   ! ID for a culvert in an open channel
    integer, parameter :: ePump             = 6   ! ID for a pump
    integer, parameter :: eValve            = 7   ! ID for a valve
    integer, parameter :: eOrifice          = 8   ! ID for an orifice
    integer, parameter :: eWeir             = 9   ! ID for a weir
    integer, parameter :: eStorage          = 10  ! ID for a storage
    integer, parameter :: eBCup             = 11  ! ID for face upstream BC
    integer, parameter :: eBCdn             = 12  ! ID for face downstream BC

    ! date types for elemI(:,ei_weir_elem_type)
    integer, parameter :: eTransverseWeir   = 1  ! ID for rectangular transverse weir
    integer, parameter :: eSideFlowWeir     = 2  ! ID for rectangular sideflow weir
    integer, parameter :: eRoadWayWeir      = 3  ! ID for rectangular roadway weir
    integer, parameter :: eVnotchWeir       = 4  ! ID for triangular v-notch weir
    integer, parameter :: eTrapezoidalWeir  = 5  ! ID for trapezoidal weir

    ! date types for elemI(:,ei_orif_elem_type)
    integer, parameter :: eBottomOrifice    = 1  ! ID for bottom orifice
    integer, parameter :: eSideOrifice      = 2  ! ID for side orifice

    ! date types for elemI(:,ei_pump_elem_type)
    integer, parameter :: eType1Pump        = 1  ! ID for Type 1 pump
    integer, parameter :: eType2Pump        = 2  ! ID for Type 2 pump
    integer, parameter :: eType3Pump        = 3  ! ID for Type 3 pump
    integer, parameter :: eType4Pump        = 4  ! ID for Type 4 pump

    ! data types for faceI(:,fi_type)
    integer, parameter :: fChannel          = eChannel  ! ID for open channel on both sides
    integer, parameter :: fPipe             = ePipe  ! ID for pipe on both sides
    integer, parameter :: fWeir             = eWeir  ! ID for pipe on both sides
    integer, parameter :: fOrifice          = eOrifice
    integer, parameter :: fPump             = ePump
    integer, parameter :: fMultiple         = eJunctionChannel  ! ID for moderation by separate up/dn element types
    integer, parameter :: fBCup             = eBCup   ! ID for face upstream BC
    integer, parameter :: fBCdn             = eBCdn   ! ID for face downstream BC

    ! date types for elemI(:,ei_geometry)
    integer, parameter :: eRectangular      = 1  ! ID for rectangular chanel, weir
    integer, parameter :: eParabolic        = 2  ! ID for parabolic channel
    integer, parameter :: eTrapezoidal      = 3  ! ID for trapezoidal channel, weir
    integer, parameter :: eTriangular       = 4  ! ID for triangular channel, weir
    integer, parameter :: eWidthDepth       = 5  ! ID for general geometry by data pairs
    integer, parameter :: eCircular         = 6  ! ID for circular pipe, orifice

    ! data types for elemI(:,ei_roughness_type)
    integer, parameter :: eManningsN        = 1   ! ID for mannings n for roughness_type
    integer, parameter :: eCD               = 2   ! ID for using drag coefficient for roughness_type

    ! data types for faceI(:,jump_type)
    integer, parameter :: jump_none         = 0   ! ID for no jump
    integer, parameter :: jump_downstream   = 1   ! ID for jump in nominal downstream direction
    integer, parameter :: jump_upstream     = 2   ! ID for jump in nominal upstream direction

    ! data types for nodeI(:,ni_node_type)
    integer, parameter :: nJ2               = 1     ! ID for junction with 2 links
    integer, parameter :: nJm               = 2     ! ID for junction with multiple links
    integer, parameter :: nStorage          = 3     ! ID for stroage units
    integer, parameter :: nBCdn             = 4     ! ID for downstream BC
    integer, parameter :: nBCup             = 5     ! iD for upstream BC

    ! data types for nodeI(:,ni_assigned) for assignment to faces and links
    integer, parameter :: nUnassigned       = nullvalueI
    integer, parameter :: nAssigned         = 1
    integer, parameter :: nDeferred         = -1

    ! data types for linkI(:,li_link_type)
    ! note that these must correspond to element types
    integer, parameter :: lchannel         = eChannel     ! ID for link that is open channel
    integer, parameter :: lpipe            = ePipe        ! ID for link that is pipe
    integer, parameter :: lweir            = eWeir        ! ID for link that is weir
    integer, parameter :: lOrifice         = eOrifice     ! ID for link that is orifice
    integer, parameter :: lPump            = ePump        ! ID for link that is pump

    ! date types for linkI(:,li_weir_type)
    integer, parameter :: lTransverseWeir  = eTransverseWeir  ! ID for rectangular transverse weir
    integer, parameter :: lSideFlowWeir    = eSideFlowWeir    ! ID for rectangular sideflow weir
    integer, parameter :: lRoadWayWeir     = eRoadWayWeir     ! ID for rectangular roadway weir
    integer, parameter :: lVnotchWeir      = eVnotchWeir      ! ID for triangular v-notch weir
    integer, parameter :: lTrapezoidalWeir = eTrapezoidalWeir ! ID for trapezoidal weir

    ! date types for linkI(:,li_orif_type)
    integer, parameter :: lBottomOrifice   = eBottomOrifice     ! ID for bottom orifice
    integer, parameter :: lSideOrifice     = eSideOrifice       ! ID for side orifice

    ! date types for elemI(:,li_pump_type)
    integer, parameter :: lType1Pump       = eType1Pump  ! ID for Type 1 pump
    integer, parameter :: lType2Pump       = eType2Pump  ! ID for Type 2 pump
    integer, parameter :: lType3Pump       = eType3Pump  ! ID for Type 3 pump
    integer, parameter :: lType4Pump       = eType4Pump  ! ID for Type 4 pump

    ! data types for linkI(:,li_geometry) (must corresponde with ei_geometry)
    integer, parameter :: lRectangular     = eRectangular     ! ID for link that rectangular channel, weir
    integer, parameter :: lParabolic       = eParabolic       ! ID for parabolic channel
    integer, parameter :: lTrapezoidal     = eTrapezoidal     ! ID for trapezoidal channel, weir
    integer, parameter :: lTriangular      = eTriangular      ! ID for triangle channel, weir
    integer, parameter :: lWidthDepth      = eWidthDepth      ! ID for general geometry by data pairs
    integer, parameter :: lCircular        = eCircular        ! ID for circular pipe, orifice

    ! data types for linkII(:,li_roughness_type)
    integer, parameter :: lManningsN       = eManningsN   ! ID for mannings n for roughness_type
    integer, parameter :: lCD              = eCD   ! ID for using drag coefficient for roughness_type

    ! data types for linkI(:,li_assigned) for assignment to faces and links
    integer, parameter :: lUnassigned      = nullvalueI
    integer, parameter :: lAssigned        = 1
    integer, parameter :: lDeferred        = -1

    ! data types for bcdata
    integer, parameter :: bc_updn_downstream        = 1
    integer, parameter :: bc_updn_upstream          = 0
    integer, parameter :: bc_category_elevation     = 0
    integer, parameter :: bc_category_inflowrate    = 1

    ! data types for Momentum Source type (setting%Solver%MomentumSourceM)
    integer, parameter :: T00 = 1
    integer, parameter :: T10 = 2
    integer, parameter :: T20 = 3

    ! data types for solver (setting%Solver%SolverSelect)
    integer, parameter :: SVE = 1
    integer, parameter :: SVE_AC = 2
    integer, parameter :: AC = 3

    ! data types for adjust approach (setting%Adjust)
    integer, parameter :: vshape = 1
    integer, parameter :: smoothall = 2

    ! data types for limiter BC approach (setting%Limiter%BC%approach)
    integer, parameter :: FroudeNumber = 1

    !==========================================================================
    ! END OF MODULE data_keys
    !==========================================================================
end module data_keys