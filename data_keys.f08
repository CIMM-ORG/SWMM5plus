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
    
    ! data types for elemI(:,ei_elem_type). faceI(:,e_type_u), faceI(:,e_type_d)
    integer, parameter :: eChannel          = 1   ! ID for an open channel element
    integer, parameter :: ePipe             = 2   ! ID for an pipe
    integer, parameter :: eJunctionChannel  = 3   ! ID for a junction 
    integer, parameter :: eJunctionPipe     = 4   ! ID for a junction 
    integer, parameter :: eCulvert          = 5   ! ID for a culvert in an open channel
    integer, parameter :: ePumpChannel      = 6   ! ID for a pump
    integer, parameter :: ePumpPipe         = 7   ! ID for a pump
    integer, parameter :: eValveChannel     = 8   ! ID for a valve
    integer, parameter :: eValvePipe        = 9   ! ID for a valve
    integer, parameter :: eOrificeChannel   = 10  ! ID for an orifice
    integer, parameter :: eOrificePipe      = 11  ! ID for an orifice
    integer, parameter :: eWeirChannel      = 12  ! ID for a weir
    integer, parameter :: eWeirPipe         = 13  ! ID for a weir
    integer, parameter :: eBCup             = 14  ! ID for face upstream BC
    integer, parameter :: eBCdn             = 15  ! ID for face downstream BC
    
    ! data types for faceI(:,fi_type)
    integer, parameter :: fChannel          = eChannel  ! ID for open channel on both sides
    integer, parameter :: fPipe             = ePipe  ! ID for pipe on both sides
    integer, parameter :: fMultiple         = eJunctionChannel  ! ID for moderation by separate up/dn element types
    integer, parameter :: fBCup             = eBCup   ! ID for face upstream BC
    integer, parameter :: fBCdn             = eBCdn   ! ID for face downstream BC    
    
    ! date types for elemI(:,ei_geometry)
    integer, parameter :: eRectangularChannel = 1  ! ID for rectangular channel
    integer, parameter :: eParabolicChannel   = 2  ! ID for parabolic channel
    integer, parameter :: eTrapezoidalChannel = 3  ! ID for trapezoidal channel
    integer, parameter :: eWidthDepthChannel  = 4  ! ID for general geometry by data pairs
    integer, parameter :: eCircularOrifice    = 5  ! ID for orifice
    integer, parameter :: eTransverseWeir     = 6  ! ID for transverse weir
    integer, parameter :: eSideFlowWeir       = 7  ! ID for sideflow weir
    integer, parameter :: eRoadwayWeir        = 8  ! ID for roadway weir
    integer, parameter :: eVnotchWeir         = 9  ! ID for V_notch weir
    integer, parameter :: eTrapezoidalWeir    = 10 ! ID for trapezoidal weir
    

    ! data types for elemI(:,ei_roughness_type)
    integer, parameter :: eManningsN    = 1   ! ID for mannings n for roughness_type
    integer, parameter :: eCD           = 2   ! ID for using drag coefficient for roughness_type
 
    ! data types for faceI(:,jump_type)
    integer, parameter :: jump_none       = 0   ! ID for no jump
    integer, parameter :: jump_downstream = 1   ! ID for jump in nominal downstream direction
    integer, parameter :: jump_upstream   = 2   ! ID for jump in nominal upstream direction
 
    ! data types for nodeI(:,ni_node_type)
    integer, parameter :: nJ2            = 1     ! ID for junction with 2 links
    integer, parameter :: nJm            = 2     ! ID for junction with multiple links
    integer, parameter :: nBCdn          = 3     ! ID for downstream BC 
    integer, parameter :: nBCup          = 4     ! iD for upstream BC 
    
    ! data types for nodeI(:,ni_assigned) for assignment to faces and links
    integer, parameter :: nUnassigned   = nullvalueI
    integer, parameter :: nAssigned     = 1
    integer, parameter :: nDeferred     = -1
    
    ! data types for linkI(:,li_link_type)
    ! note that these must correspond to element types
    integer, parameter :: lchannel      = eChannel     ! ID for link that is open channel
    integer, parameter :: lpipe         = ePipe     ! ID for link that is pipe
    
   ! data types for linkI(:,li_geometry) (must corresponde with ei_geometry)
    integer, parameter :: lRectangularChannel = eRectangularChannel ! ID for link that rectangular channel
    integer, parameter :: lParabolicChannel   = eParabolicChannel   ! ID for parabolic channel
    integer, parameter :: lTrapezoidalChannel = eTrapezoidalChannel ! ID for trapezoidal channel
    integer, parameter :: lWidthDepthChannel  = eWidthDepthChannel  ! ID for general geometry by data pairs
    integer, parameter :: lCircularOrifice    = eCircularOrifice    ! ID for orifice
    integer, parameter :: lTransverseWeir     = eTransverseWeir     ! ID for transverse weir
    integer, parameter :: lSideFlowWeir       = eSideFlowWeir       ! ID for sideflow weir
    integer, parameter :: lRoadwayWeir        = eRoadwayWeir        ! ID for roadway weir
    integer, parameter :: lVnotchWeir         = eVnotchWeir         ! ID for V_notch weir
    integer, parameter :: lTrapezoidalWeir    = eTrapezoidalWeir    ! ID for trapezoidal weir

    ! data types for linkII(:,li_roughness_type)
    integer, parameter :: lManningsN    = eManningsN   ! ID for mannings n for roughness_type
    integer, parameter :: lCD           = eCD   ! ID for using drag coefficient for roughness_type
 
    ! data types for linkI(:,li_assigned) for assignment to faces and links
    integer, parameter :: lUnassigned   = nullvalueI
    integer, parameter :: lAssigned     = 1
    integer, parameter :: lDeferred     = -1

    ! data types for bcdata
    integer, parameter :: bc_updn_downstream = 1
    integer, parameter :: bc_updn_upstream   = 0   
    integer, parameter :: bc_category_elevation = 0    
    integer, parameter :: bc_category_inflowrate = 1    

!==========================================================================
! END OF MODULE data_keys
!==========================================================================    
 end module data_keys
