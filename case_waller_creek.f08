! module case_waller_creek
!
! Test case for Waller_Creek case.
!
!==========================================================================
!
 module case_waller_creek
! 
    use allocate_storage
    use array_index
    use bc
    use data_keys
    use globals
    use setting_definition
    use read_width_depth

    implicit none
    
    private
    
    public :: case_waller_creek_initialize

    integer :: debuglevel = 0
    
 contains
!
!========================================================================== 
!==========================================================================
!
 subroutine case_waller_creek_initialize &
    (channel_length, channel_breadth, subdivide_length, lowerZ, upperZ, &
    initial_flowrate, depth_upstream, depth_dnstream,   &
    ManningsN, roughness_type, idepth_type,                             &
    linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName,     &
    bcdataDn, bcdataUp)
!
! initialize the link-node system and boundary conditions for a simple channel
! 
 character(64) :: subroutine_name = 'case_waller_creek_initialize'
 
 real,  intent(in)  :: channel_length, channel_breadth, subdivide_length
 real,  intent(in)  :: lowerZ, upperZ, ManningsN, initial_flowrate
 real,  intent(in)  :: depth_upstream, depth_dnstream
 
 integer, intent(in):: roughness_type, idepth_type
 
 integer,   dimension(:,:), allocatable, target, intent(out)    :: linkI 
 integer,   dimension(:,:), allocatable, target, intent(out)    :: nodeI
 
 real,      dimension(:,:), allocatable, target, intent(out)    :: linkR 
 real,      dimension(:,:), allocatable, target, intent(out)    :: nodeR 
 
 logical,   dimension(:,:), allocatable, target, intent(out)    :: linkYN
 logical,   dimension(:,:), allocatable, target, intent(out)    :: nodeYN
 
 type(string), dimension(:), allocatable, target, intent(out)   :: linkName 
 type(string), dimension(:), allocatable, target, intent(out)   :: nodeName
 
 type(bcType), dimension(:), allocatable, intent(out) :: bcdataUp, bcdataDn

 integer    :: ntimepoint, ndnstreamBC, nupstreamBC
 
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
! Boundary conditions
 ntimepoint = 2
 nupstreamBC = 1
 ndnstreamBC = 1
 
! check if  
 
 call bc_allocate &
    (bcdataDn, bcdataUp, ndnstreamBC, nupstreamBC, ntimepoint) 

! assign values
! downstream is default to elevation
 bcdataDn(1)%NodeID = 2
 bcdataDn(1)%TimeArray(1)     = setting%Time%StartTime 
 bcdataDn(1)%TimeArray(2)     = setting%Time%EndTime + 100.0 !s
 bcdataDn(1)%ValueArray(1)    = lowerZ +  depth_dnstream  ! m
 bcdataDn(1)%ValueArray(2)    = lowerZ +  depth_dnstream ! m

! upstream is default to flowrate
 bcdataUp(1)%NodeID = 1
 bcdataUp(1)%TimeArray(1)  = setting%Time%StartTime
 bcdataUp(1)%TimeArray(2)  = setting%Time%EndTime + 100.0 !s
 bcdataUp(1)%ValueArray(1) = initial_flowrate  ! m^3/s
 bcdataUp(1)%ValueArray(2) = initial_flowrate  ! m^3/2
    
 call case_waller_creek_links_and_nodes &
    (channel_length, channel_breadth, subdivide_length, lowerZ, upperZ, &
    initial_flowrate, depth_upstream, depth_dnstream, ManningsN, &
    roughness_type,  idepth_type, &
    linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName)

 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine case_waller_creek_initialize
!
!========================================================================== 
!
! PRIVATE BELOW HERE
!
!==========================================================================
!
subroutine define_geometry (geometry_downstream_minimum_length, &
    Waller_Creek_cellsize_target, n_rows_in_file_node)
 
 character(64) :: subroutine_name = 'define_geometry'
 
 real, dimension(:), allocatable, target, intent(inout) :: Length
 real, dimension(:), allocatable, target, intent(inout) :: xDistance
 
 real, intent(in) :: geometry_downstream_minimum_length
 real, intent(in) :: n_rows_in_file_node
 
 !private
 real :: oldX, oldL
 integer :: NXold, ncell
 
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 !stretching out the last cell
 if (geometry_downstream_minimum_length > 0.0) then
    oldX = xDistance(size(xDistance))
    oldL = Length(size(Length))
    if (geometry_downstream_minimum_length > oldL) then
        Length(size(Length)) = geometry_downstream_minimum_length
        xDistance(size(xDistance)) = oldX &
            + 0.5*(geometry_downstream_minimum_length-oldL)
    endif
 endif
 
 !splitting domain into smaller cells
 ncell = 0
 NXold = n_rows_in_file_node
 if (Waller_Creek_cellsize_target > 0.0) then
    !temporary creation of zbottom and xvalue on face for establishing 
    !interpolated zbottom throughout the smaller cells
    
 endif
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine define_geometry
!
!========================================================================== 
!==========================================================================
!
subroutine widthdepth_pair_consistency (NX, widthDepthData, cellType)
 
 character(64) :: subroutine_name = 'widthdepth_pair_consistency'
 
 integer, dimension(:), allocatable :: ID
 integer, dimension(:), allocatable :: numberPairs
 
 real, dimension(:), allocatable :: ManningsN
 real, dimension(:), allocatable :: Length
 real, dimension(:), allocatable :: zBottom
 real, dimension(:), allocatable :: xDistance
 real, dimension(:), allocatable :: Breadth
 
 real, dimension(:,:,:), allocatable :: widthDepthData
 
 real, dimension(:,:,:), allocatable :: dWidth
 real, dimension(:,:,:), allocatable :: dDepth
 
 character(len=:), allocatable :: cellType(:)
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 

 allocate(dWidth(size(widthDepthData,1),size(widthDepthData,2)))
 dWidth(:,:) = 0.0
 allocate(dDepth(size(widthDepthData,1),size(widthDepthData,2)))
 dDepth(:,:) = 0.0
 
 do i=1, size(numberPairs)
    do j=1, numberPairs(i)
        dWidth(i,1:j-1) = widthDepthData(i,2:j,1) - widthDepthData(i,1:j-1,1)
        dDepth(i,1:j-1) = widthDepthData(i,2:j,2) - widthDepthData(i,1:j-1,2)
    enddo
 enddo
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine widthdepth_pair_consistency
!
!========================================================================== 
!==========================================================================
!
 subroutine case_waller_creek_links_and_nodes &
    (channel_length, channel_breadth, subdivide_length, lowerZ, upperZ, &
     initial_flowrate, depth_upstream, depth_dnstream, ManningsN,       &
     roughness_type, idepth_type, &
     linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName)
!
! creates a simple rectangular channel with 1 link and 2 nodes
! 
 character(64) :: subroutine_name = 'case_waller_creek_links_and_nodes'
 
 real,  intent(in)  :: channel_length, channel_breadth, subdivide_length
 real,  intent(in)  :: lowerZ, upperZ, ManningsN, initial_flowrate
 real,  intent(in)  :: depth_upstream, depth_dnstream
 
 integer, intent(in):: roughness_type, idepth_type
 
 integer,   dimension(:,:), allocatable, target, intent(out)    :: linkI 
 integer,   dimension(:,:), allocatable, target, intent(out)    :: nodeI
 
 real,      dimension(:,:), allocatable, target, intent(out)    :: linkR 
 real,      dimension(:,:), allocatable, target, intent(out)    :: nodeR 
 
 logical,   dimension(:,:), allocatable, target, intent(out)    :: linkYN
 logical,   dimension(:,:), allocatable, target, intent(out)    :: nodeYN
 
 type(string), dimension(:), allocatable, target, intent(out)   :: linkName 
 type(string), dimension(:), allocatable, target, intent(out)   :: nodeName
 
 integer :: ii
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 N_link = 1
 N_node = 2
 
 call allocate_linknode_storage &
    (linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName)
 
! assign the indexes
 linkI(:,li_idx) = (/ (ii, ii=1,N_link) /)
 nodeI(:,ni_idx) = (/ (ii, ii=1,N_node) /)
 
! assign no names for links
 do ii=1,N_link
    linkName(ii)%str = "Channel" 
 end do
    
! assign zeros for accumulators
 nodeI(:,ni_N_link_d) = 0    
 nodeI(:,ni_N_link_u) = 0   
 
! assign uniform physical data 
 linkI(:,li_roughness_type)  = roughness_type
 linkR(:,lr_Roughness)       = ManningsN
 
! designate the upstream nodes
 nodeI(1,ni_node_type) = nBCup

 nodeR(1,nr_Zbottom) = upperZ
 
 nodeName(1)%str = 'UpstreamBC'
    
! designate the downstream node
 nodeI(2,ni_node_type) = nBCdn

 nodeR(2,nr_Zbottom) = lowerZ
 
 nodeName(2)%str = 'DownstreamBC'

! assign the link types
 linkI(:,li_link_type) = lChannel

! assign all as rectangular channels
 linkI(:,li_geometry) = lWidthDepth

! assign the link position and mappings

 linkI(1,li_Mnode_u) = 1 ! map to upstream node
 linkI(1,li_Mnode_d) = 2 ! map to downstream node

 linkR(1,lr_Length)          = channel_length
 linkR(1,lr_BreadthScale)    = channel_breadth 
 linkR(1,lr_ElementLength)   = subdivide_length
 linkR(1,lr_InitialFlowrate) = initial_flowrate
 linkR(1,lr_InitialUpstreamDepth)    = depth_upstream
 linkR(1,lr_InitialDnstreamDepth)    = depth_dnstream
 linkI(1,li_InitialDepthType)        = idepth_type
 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) then
    print *
    print *, subroutine_name,'-----------------------------------'
    print *, 'link info'
    print *, linkI(:,li_idx), ' idx'
    print *, linkI(:,li_link_type), ' type'
    print *, linkI(:,li_Mnode_u) , ' upstream node'
    print *, linkI(:,li_Mnode_d) , ' downstream node'
    print *,
    print *, 'node info'
    print *, nodeI(:,ni_idx), ' idx'
    print *, nodeI(:,ni_node_type), ' type'
    print *, nodeI(:,ni_N_link_d), 'number of downstream links'
    print *, nodeI(:,ni_N_link_u), 'number of upstream links'
 endif
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine case_waller_creek_links_and_nodes
!
!========================================================================== 
! END OF MODULE case_waller_creek
!==========================================================================
 end module case_waller_creek
