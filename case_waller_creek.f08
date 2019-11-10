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
subroutine widthdepth_pair_auxiliary (NX, widthDepthData, cellType)
 
 character(64) :: subroutine_name = 'widthdepth_pair_auxiliary'
 
 integer, dimension(:), allocatable :: ID
 integer, dimension(:), allocatable :: numberPairs
 
 real, dimension(:), allocatable :: ManningsN
 real, dimension(:), allocatable :: Length
 real, dimension(:), allocatable :: zBottom
 real, dimension(:), allocatable :: xDistance
 real, dimension(:), allocatable :: Breadth
 
 real, dimension(:,:,:), allocatable :: widthDepthData
 
 real, dimension(:,:), allocatable :: dWidth
 real, dimension(:,:), allocatable :: dDepth
 
 character(len=:), allocatable :: cellType(:)
 
 real, pointer :: width(:,:), depth(:,:), area(:,:), areaTBL(:,:)
 real, pointer :: dWidth(:,:), dDepth(:,:), angle(:,:), perimeterBL(:,:)
 
 integer :: ii,jj
 integer :: eIn1
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 eIn1 = size(widthDepthData,2)
 
 width       => widthDepthData (:,:, wd_widthThisLayer)
 depth       => widthDepthData (:,:, wd_depthAtLayerTop)
 area        => widthDepthData (:,:, wd_areaThisLayer)
 areaTBL     => widthDepthData (:,:, wd_areaTotalBelowThisLayer)
 dWidth      => widthDepthData (:,:, wd_Dwidth)
 dDepth      => widthDepthData (:,:, wd_Ddepth)
 angle       => widthDepthData (:,:, angle)
 perimeterBL => widthDepthData (:,:, perimeterBelowThisLayer)
 
 ! lowest layer is triangular
 area(:,1) = onehalfR * width(:,1) * depth(:,1)
        
 ! the area in this layer
 area(:,2:eIn1) = onehalfR * (width(:,2:eIn1) + width(:,1:eIn1-1) &
                           * (depth(:,2:eIn1) + depth(:,1:eIn1-1))

        
 ! set areas to zero above the uppermost pair
 where (cellType(:) == 'widthdepth_pair')
        area(:, numberPairs(:):eIn1) = 0.0
 endwhere
 
 ! store width and depth differences 
 dWidth(:,1) = width(:,1)
 dDepth(:,1) = depth(:,1)
 
 ! delta width between top and bottom of this layer
 dWidth(:,2:eIn1-1) = width(:,2:eIn1-1) - width:,1:eIn1-2
 ! delta depth between top and bottom of this layer
 dDepth(:,2:eIn1-1) = depth(:,2:eIn1-1) - depth:,1:eIn1-2
 
 where(dWidth > setting%Method%AdjustWidthDepth%SmallWidth)
    ! pairs that are not 90 degree angles
    angle = atan(2.0 * dDepth / dWidth)
 elsewhere
    ! a near-90 degree angle will have an infinite tangent.
    ! we handle this case by setting all these angles to pi/2 - small value 
    angle = pi/onehalfR - setting%Method%AdjustWidthDepth%angleMinimum
 endwhere
 
 ! accumulated area of the trapezoids to the ii width-depth level
 areaTBL(:,1) = 0.0
 areaTBL(:,2:eIn1) = areaTBL(:,1:eIn1-1) + area(:,1:eIn1-1)
 
 ! check that the setting maximum area value is greater than any accumulated 
 ! area at the uppermost level.
 if (setting%Method%AdjustWidthDepth%areaMaximum < maxval(areaTBL(:, eIn1))) then
    setting%Method%AdjustWidthDepth%areaMaximum = 2.0 * maxval(areaTBL(:, eIn1))
 endif
 
 ! perimeter below this layer
 perimeterBL(:,1) = 0.0
 perimeterBL(:,2:eIn1) = perimeterBL(:,1:eIn1-1) &
        + 2.0 * sqrt(dDepth(:,1:eIn1-1)**2.0 + 0.5 * dWidth(:,1:eIn1-1)**2.0)
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine widthdepth_pair_auxiliary
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
 
 real, dimension(:,:), allocatable :: dWidth
 real, dimension(:,:), allocatable :: dDepth
 
 character(len=:), allocatable :: cellType(:)
 
 integer :: ii,jj, nfix
 integer :: width = wd_widthAtLayerTop
 integer :: depth = wd_depthAtLayerTop
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 

 allocate(dWidth(size(widthDepthData,1),size(widthDepthData,2)))
 dWidth(:,:) = 0.0
 allocate(dDepth(size(widthDepthData,1),size(widthDepthData,2)))
 dDepth(:,:) = 0.0
 
 nfix = 0
 do ii=1, size(numberPairs)
    do jj=1, numberPairs(i)
        !compute the difference in width across each level
        dWidth(ii,1:jj-1) = widthDepthData(ii,2:jj,width) &
                            - widthDepthData(ii,1:jj-1,width)
        !compute difference in depth across eacg level
        dDepth(ii,1:jj-1) = widthDepthData(ii,2:jj,depth) &
                            - widthDepthData(ii,1:jj-1,depth)
    enddo
 enddo
 
 !negative values indicate non-monotonic behavior that can be fixed.
 nfix = nfix + count(dWidth < 0.0 .or. dDepth < 0.0)
 
 if ((setting%Method%AdjustWidthDepth == .true.) .and. (nfix > 0)) then
    call widthdepth_pair_fix(widthDepthData)
 endif
 
 !check that the width-depth pairs cover enough depth and fix with vertical walls
 !width-Depth matrix has an extra cell when allocated
 do ii=1, size(numberPairs)
    if (maxval(widthDepthData(ii,:,depth) &
        < setting%Method%AdjustWidthDepth%DepthMaxExpected) then
        
        widthDepthData(ii,numberPairs(ii)+1,depth) 
                        = 2.0*setting%Method%AdjustWidthDepth%DepthMaxExpected
                        
        widthDepthData(ii,numberPairs(ii)+1,width) 
                        = widthDepthData(ii,numberPairs(ii),1) 
                        
        !an additional pair has been added at this element
        numberPairs(ii) = numberPairs(ii) + 1
    endif
 enddo
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine widthdepth_pair_consistency
!
!========================================================================== 
!==========================================================================
!
subroutine widthdepth_pair_fix (widthDepthData)
 
 character(64) :: subroutine_name = 'widthdepth_pair_fix'
 
 real, pointer :: up2W(:,:), up1W(:,:), lowW(:,:)
 real, pointer :: up2D(:,:), up1D(:,:), lowD(:,:)
 
 integer :: width = wd_widthAtLayerTop
 integer :: depth = wd_depthAtLayerTop
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 !width 2 layers above
 up2W => widthDepthData(:, 3:size(widthDepthData,2)-1, width)
 !width 1 layer above
 up1W => widthDepthData(:, 2:size(widthDepthData,2)-2, width)
 !width this layer (below)
 lowW => widthDepthData(:, 1:size(widthDepthData,2)-3, width)

 !depth 2 layers above
 up2D => widthDepthData(:, 3:size(widthDepthData,2)-1, depth)
 !depth 1 layer above
 up1D => widthDepthData(:, 2:size(widthDepthData,2)-2, depth)
 !depth this layer (below)
 lowD => widthDepthData(:, 1:size(widthDepthData,2)-3, depth)
 
 !fix inconsistent depth
 where (up1D <= lowD)
    where (up2D <= lowD)
        ! multi-level inconsistency - simple expansion
        up1D = lowD*(1.0 + setting%Method%AdjustWidthDepth%AdjustFraction)
    elsewhere
        !single-level inconsistency
        where ((up1W > lowW) .and. (up2W > up1W))
            !width is consistent - use linear interpolation
            up1D = lowD + (up2D - lowD)*(up1W - lowW)/(up2W - lowW)
        elsewhere
            !both depth and width are inconsistent - simple expansion
            up1D = lowD*(1.0 + setting%Method%AdjustWidthDepth%AdjustFraction)
        endwhere
    endwhere
 endwhere
 
 !fix inconsistent widht
 where (up1W < lowW)
    where (up2W < lowW)
        ! multi-level inconsistency - simple expansion
        up1W = lowW*(1.0 + setting%Method%AdjustWidthDepth%AdjustFraction)
    elsewhere
        !single-level inconsistency
        where ((up1D > lowD) .and. (up2D > up1D))
            !width is consistent - use linear interpolation
            up1W = lowW + (up2W - lowW)*(up1D - lowD)/(up2D - lowD)
        elsewhere
            !both widht and width are inconsistent - simple expansion
            up1W = lowW*(1.0 + setting%Method%AdjustWidthDepth%AdjustFraction)
        endwhere
    endwhere
 endwhere
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine widthdepth_pair_fix
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
