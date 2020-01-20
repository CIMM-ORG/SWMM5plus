! module case_waller_creek
!
! Test case for Waller_Creek case.
!
! @EhsanMadadi
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
    use utility

    implicit none
    
    private
    
    public :: case_waller_creek_initialize
    public :: nonmonotonic_subdivide
    public :: widthdepth_pair_auxiliary

    integer :: debuglevel = 0
    
 contains
!
!========================================================================== 
!==========================================================================
!
 subroutine case_waller_creek_initialize &
    ()
!
! initialize the link-node system and boundary conditions for a simple channel
! 
 character(64) :: subroutine_name = 'case_waller_creek_initialize'
 
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 

 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine case_waller_creek_initialize
!
!========================================================================== 
!
! PRIVATE BELOW HERE
!
!==========================================================================
!
! subroutine define_geometry (geometry_downstream_minimum_length, &
!     n_rows_in_file_node)
!  
!  character(64) :: subroutine_name = 'define_geometry'
!  
!  real, dimension(:), allocatable, target, intent(inout) :: Length
!  real, dimension(:), allocatable, target, intent(inout) :: xDistance
!  integer, dimension(:), allocatable, target, intent(inout) :: nadd
!  real, dimension(:), allocatable, target, intent(inout) :: xface
!  real, dimension(:), allocatable, target, intent(inout) :: dx
!  
!  real, dimension(:,:,:), allocatable :: outputWidthDepthData
!  
!  real, intent(in) :: geometry_downstream_minimum_length
!  real, intent(in) :: n_rows_in_file_node
!  
!  real, pointer :: nWidth(:,:), nDepth(:,:), nArea(:,:), nAreaTBL(:,:)
!  real, pointer :: ndWidth(:,:), ndDepth(:,:), nAngle(:,:), nPerimeterBL(:,:)
!  
!  !private
!  real :: oldX, oldL
!  integer :: NXold, ncell, NX
!  
!   
! !-------------------------------------------------------------------------- 
!  if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
!  
!  !stretching out the last cell
!  if (geometry_downstream_minimum_length > 0.0) then
!     oldX = xDistance(size(xDistance))
!     oldL = Length(size(Length))
!     if (geometry_downstream_minimum_length > oldL) then
!         Length(size(Length)) = geometry_downstream_minimum_length
!         xDistance(size(xDistance)) = oldX &
!             + 0.5*(geometry_downstream_minimum_length-oldL)
!     endif
!  endif
!  
!  !splitting domain into smaller cells
!  ncell = 0
!  NX = n_rows_in_file_node
!  NXold = NX
!  if (setting%Method%AdjustWidthDepth%cellSizeTarget > 0.0) then
!     ! check widthdepth pair geometry for consistency
!     call widthdepth_pair_consistency (NX, widthDepthData, cellType)
!     ! compute additional geometry data for widthdepth pairs
!     call widthdepth_pair_auxiliary (NX, widthDepthData, cellType)
!     
!     ! cycle through to find the number of cells to add at each cross-section
!     allocate(nadd(size(Length,1))
!     nadd(:) = 0
!     where(Length > (1.5*setting%Method%AdjustWidthDepth%cellSizeTarget))
!         nadd = int(Length/setting%Method%AdjustWidthDepth%cellSizeTarget)
!     elsewhere
!         nadd = 1
!     endwhere
!     
!     ncell = ncell + sum(nadd)
!     
!     NX = ncell
!     allocate(outputWidthDepthData(NX, size(widthDepthData,2), wd_idx_max), stat=allocation_status, errmsg=emsg)
!     outputWidthDepthData(:,:,:) = 0.0
!     
!     allocate(xface(size(Length,1))
!     allocate(dx(size(Length,1))
!     
!     xface = xDistance - 0.5 * Length
!     xface(NXold) = xDistance(NXold-1) + Length(NXold-1)
!     
!     dx = Length/nadd
!     
!     nWidth       => outputWidthDepthData (:,:, wd_widthThisLayer)
!     nDepth       => outputWidthDepthData (:,:, wd_depthAtLayerTop)
!     nArea        => outputWidthDepthData (:,:, wd_areaThisLayer)
!     nAreaTBL     => outputWidthDepthData (:,:, wd_areaTotalBelowThisLayer)
!     ndWidth      => outputWidthDepthData (:,:, wd_Dwidth)
!     ndDepth      => outputWidthDepthData (:,:, wd_Ddepth)
!     nAngle       => outputWidthDepthData (:,:, angle)
!     nPerimeterBL => outputWidthDepthData (:,:, perimeterBelowThisLayer)
!     
!     outputLength = dx
!     outputXDistance = xface + 0.5*dx
!     outputManningsN = ManningsN
!     
!     
!  endif
!  
!  
!  
!  
!  
!  
!  
!  
!  
!  
!  
!  
!  
!  
!  if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
!  end subroutine define_geometry
!
!========================================================================== 
!==========================================================================
!
 subroutine nonmonotonic_subdivide &
     (ID, numberPairs, ManningsN, Length, zBottom, xDistance,                  &
      Breadth, widthDepthData, cellType, NX, faceZBottom, max_number_of_pairs, &
      newID, newNumberPairs, newManningsN, newLength, newZBottom,              &
      newXDistance, newBreadth, newWidthDepthData, newCellType,                &
      subdivide_length_check)
!
! initialize the link-node system and boundary conditions for a simple channel
! 
 character(64) :: subroutine_name = 'nonmonotonic_subdivide'
 
 
 integer, intent(inout) :: NX
 integer, intent(inout) :: max_number_of_pairs
 integer :: newNX = 0
 integer :: ii, jj
 
 integer :: allocation_status
 character(len=99) :: emsg
 
 integer, intent(in) :: ID(:)
 integer, intent(in) :: numberPairs(:)
 real, intent(in)    :: ManningsN(:)
 real, intent(in)    :: Length(:)
 real, intent(in)    :: zBottom(:)
 real, intent(in)    :: xDistance(:)
 real, intent(in)    :: Breadth(:)
 real, intent(in)    :: widthDepthData(:,:,:)
 character, intent(in)    :: cellType(:)
 
 real, intent(in)    :: subdivide_length_check
 
 real,    dimension(:),     allocatable :: faceZBottom
 real,    dimension(:),     allocatable :: temp1
 real,    dimension(:),     allocatable :: temp2
 integer, dimension(:),     allocatable :: isnonmonotonic
 
 integer, dimension(:),     allocatable :: newID
 integer, dimension(:),     allocatable :: newNumberPairs
 real,    dimension(:),     allocatable :: newManningsN
 real,    dimension(:),     allocatable :: newLength
 real,    dimension(:),     allocatable :: newZBottom
 real,    dimension(:),     allocatable :: newXDistance
 real,    dimension(:),     allocatable :: newBreadth
 real,    dimension(:,:,:), allocatable :: newWidthDepthData
 character(len=:),          allocatable :: newCellType(:)
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 allocate(isnonmonotonic(NX))
 isnonmonotonic= 0
 
 allocate(faceZBottom(NX+1))
 faceZbottom(:) = 0.0
 
! find the z at the faces
 call face_zbottom(faceZbottom, zbottom, Length, NX)
 
 print *, "checking for nonmonotonic elements"
 
! determine number of non-monotonic zbottom product of signs is negative
 allocate(temp1(NX))
 allocate(temp2(NX))
 temp1(1:NX) = faceZbottom (1:NX) - zbottom(1:NX)
 temp2(1:NX) = zbottom (1:NX) - faceZbottom(2:NX+1)

 ! set array for integer counting of non-monotonic cells
 where (temp1*temp2 < zeroR)
    isnonmonotonic = 1
 endwhere
 
! counting the new elements
 newNX = NX + twoI * count(isnonmonotonic/= 0)
 
 allocate(newID(newNX), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)
 newID(:) = 0
 
 allocate(newNumberPairs(newNX), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)
 newNumberPairs(:) = 0
 
 allocate(newManningsN(newNX), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)
 newManningsN(:) = 0.0
 
 allocate(newLength(newNX), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)
 newLength(:) = 0.0
 
 allocate(newZBottom(newNX), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)
 newZBottom(:) = 0.0
 
 allocate(newXDistance(newNX), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)
 newXDistance(:) = 0.0
 
 allocate(newBreadth(newNX), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)
 newBreadth(:) = 0.0
 
!+1 is a ghost cell for the chack that the width-depth pairs cover enough
!depth and fixing of vertical walls
 allocate(newWidthDepthData(newNX, max_number_of_pairs+1, wd_idx_max), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)
 newWidthDepthData(:,:,:) = 0.0
 
 allocate(character(100):: newCellType(newNX), stat=allocation_status, errmsg=emsg)
 call utility_check_allocation (allocation_status, emsg)
 
! print('Checking for non-monotonic zbottom, 1=found')
 jj = 1
 do ii = 1, NX
    if (isnonmonotonic(ii) /= 0) then
        print*, "nonmonotonic element ==", ii, faceZbottom(ii), zbottom(ii), faceZbottom(ii+1)
        newID          (jj)   = jj
        newID          (jj+1) = jj + 1
        newID          (jj+2) = jj + 2
        newNumberPairs (jj)   = numberPairs (ii)
        newNumberPairs (jj+1) = numberPairs (ii)
        newNumberPairs (jj+2) = numberPairs (ii)
        newManningsN   (jj)   = ManningsN   (ii)
        newManningsN   (jj+1) = ManningsN   (ii)
        newManningsN   (jj+2) = ManningsN   (ii)
        newLength      (jj)   = onehalfR*Length (ii)                           &
            - onehalfR * subdivide_length_check
        newLength      (jj+1) = subdivide_length_check
        newLength      (jj+2) = onehalfR*Length (ii)                           &
            - onehalfR * subdivide_length_check
        newZBottom     (jj)   = faceZbottom (ii)                               &
            - onehalfR * (faceZbottom (ii) - zBottom (ii))
        newZBottom     (jj+1) = zBottom     (ii)
        newZBottom     (jj+2) = zBottom (ii)                                   &
            - onehalfR * (zBottom (ii) - faceZbottom (ii))
        newXDistance   (jj)   = onehalfR*xDistance (ii)                        &
            - onehalfR*Length (ii) - onehalfR * subdivide_length_check
        newXDistance   (jj+1) = newXDistance   (jj) + subdivide_length_check
        newXDistance   (jj+2) = newXDistance   (jj+1)                          &
            + onehalfR*Length (ii) - onehalfR * subdivide_length_check
        newBreadth     (jj)   = Breadth     (ii)
        newBreadth     (jj+1) = Breadth     (ii)
        newWidthDepthData(jj,:,:)   = widthDepthData (ii,:,:)
        newWidthDepthData(jj+1,:,:) = widthDepthData (ii,:,:)
        newCellType    (jj)   = cellType    (ii)
        newCellType    (jj+1) = cellType    (ii)
        jj = jj + 3
    else
        newID          (jj) = jj
        newNumberPairs (jj) = numberPairs (ii)
        newManningsN   (jj) = ManningsN   (ii)
        newLength      (jj) = Length      (ii)
        newZBottom     (jj) = zBottom     (ii)
        newXDistance   (jj) = xDistance   (ii)
        newBreadth     (jj) = Breadth     (ii)
        newWidthDepthData(jj,:,:) = widthDepthData (ii,:,:)
        newCellType    (jj) = cellType    (ii)
        jj = jj + 1
    endif
 enddo
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine nonmonotonic_subdivide
!
!========================================================================== 
!==========================================================================
!
 subroutine face_zbottom (faceZbottom, zBottom, Length, NX)
!
! initialize the link-node system and boundary conditions for a simple channel
! 
 character(64) :: subroutine_name = 'face_zbottom'
 
 real, intent(inout) :: faceZbottom(:)
 real, intent(in)    :: zBottom(:)
 real, intent(in)    :: Length(:)
 integer, intent(in) :: NX
 
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 faceZbottom (2:NX) = &
         (zBottom(1:NX-1)*Length(2:NX) + zBottom(2:NX)*Length(1:NX-1)) &
        /(Length(2:NX) + Length(1:NX-1))

 faceZbottom (1)  = zBottom(1)
 faceZbottom (NX+1) = zBottom(NX)
 

 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine face_zbottom
!
!========================================================================== 
!==========================================================================
!
subroutine widthdepth_pair_auxiliary (widthDepthData, cellType, numberPairs)
 
 character(64) :: subroutine_name = 'widthdepth_pair_auxiliary'
 
 integer, intent(inout) :: numberPairs(:)
 
 real, target, intent(inout) :: widthDepthData(:,:,:)
 
 character, intent(in) :: cellType(:)
 
 real, pointer :: width(:,:), depth(:,:), area(:,:), areaTBL(:,:)
 real, pointer :: dWidth(:,:), dDepth(:,:), angle(:,:), perimeterBL(:,:)
 
 integer :: ii,jj
 integer :: eIn1
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 eIn1 = size(widthDepthData,2)
 
 width       => widthDepthData (:,:, wd_widthAtLayerTop)
 depth       => widthDepthData (:,:, wd_depthAtLayerTop)
 area        => widthDepthData (:,:, wd_areaThisLayer)
 areaTBL     => widthDepthData (:,:, wd_areaTotalBelowThisLayer)
 dWidth      => widthDepthData (:,:, wd_Dwidth)
 dDepth      => widthDepthData (:,:, wd_Ddepth)
 angle       => widthDepthData (:,:, wd_angle)
 perimeterBL => widthDepthData (:,:, wd_perimeterBelowThisLayer)
 
 ! lowest layer is triangular
 area(:,1) = onehalfR * width(:,1) * depth(:,1)
        
 ! the area in this layer
 area(:,2:eIn1) = onehalfR * (width(:,2:eIn1) + width(:,1:eIn1-1)) &
                           * (depth(:,2:eIn1) - depth(:,1:eIn1-1))

        
 ! set areas to zero above the uppermost pair
 do ii = 1, eIn1
    if(cellType(ii) == 'widthdepth_pair') then
        area(ii, numberPairs(ii):eIn1) = zeroR
    endif
 enddo
 
 ! store width and depth differences 
 dWidth(:,1) = width(:,1)
 dDepth(:,1) = depth(:,1)
 
 ! delta width between top and bottom of this layer
 dWidth(:,2:eIn1-1) = width(:,2:eIn1-1) - width(:,1:eIn1-2)
 ! delta depth between top and bottom of this layer
 dDepth(:,2:eIn1-1) = depth(:,2:eIn1-1) - depth(:,1:eIn1-2)
 
 where(dWidth > setting%Method%AdjustWidthDepth%SmallWidth)
    ! pairs that are not 90 degree angles
    angle = atan(twoR * dDepth / dWidth)
 elsewhere
    ! a near-90 degree angle will have an infinite tangent.
    ! we handle this case by setting all these angles to pi/2 - small value 
    angle = pi/twoR - setting%Method%AdjustWidthDepth%angleMinimum
 endwhere
 
 ! accumulated area of the trapezoids to the ii width-depth level
 areaTBL(:,1) = zeroR
 areaTBL(:,2:eIn1) = areaTBL(:,1:eIn1-1) + area(:,1:eIn1-1)
 
 ! check that the setting maximum area value is greater than any accumulated 
 ! area at the uppermost level.
 if (setting%Method%AdjustWidthDepth%areaMaximum < maxval(areaTBL(:, eIn1))) then
    setting%Method%AdjustWidthDepth%areaMaximum = twoR * maxval(areaTBL(:, eIn1))
 endif
 
 ! perimeter below this layer
 perimeterBL(:,1) = zeroR
 perimeterBL(:,2:eIn1) = perimeterBL(:,1:eIn1-1) &
        + twoR * sqrt(dDepth(:,1:eIn1-1)**twoR + onehalfR * dWidth(:,1:eIn1-1)**twoR)
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine widthdepth_pair_auxiliary
!
!========================================================================== 
!==========================================================================
!
subroutine widthdepth_pair_consistency (widthDepthData, numberPairs)
 
 character(64) :: subroutine_name = 'widthdepth_pair_consistency'
 
 integer, intent(inout) :: numberPairs(:)
 
 real, target, intent(inout) :: widthDepthData(:,:,:)
 
 real, dimension(:,:), allocatable :: dWidth
 real, dimension(:,:), allocatable :: dDepth
 
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
    do jj=1, numberPairs(ii)
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
 
 if ((setting%Method%AdjustWidthDepth%Apply .eqv. .true.) .and. (nfix > 0)) then
    call widthdepth_pair_fix(widthDepthData)
 endif
 
 !check that the width-depth pairs cover enough depth and fix with vertical walls
 !width-Depth matrix has an extra cell when allocated
 do ii=1, size(numberPairs)
    if (maxval(widthDepthData(ii,:,depth)) &
        < setting%Method%AdjustWidthDepth%DepthMaxExpected) then
        
        widthDepthData(ii,numberPairs(ii)+1,depth) &
                        = 2.0*setting%Method%AdjustWidthDepth%DepthMaxExpected
                        
        widthDepthData(ii,numberPairs(ii)+1,width) &
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
 
 real, target, intent(inout) :: widthDepthData(:,:,:)
 
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
    ()
!
! creates a simple rectangular channel with 1 link and 2 nodes
! 
 character(64) :: subroutine_name = 'case_waller_creek_links_and_nodes'
 
  
!-------------------------------------------------------------------------- 
 if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name 
 
 
 if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
 end subroutine case_waller_creek_links_and_nodes
!
!========================================================================== 
! END OF MODULE case_waller_creek
!==========================================================================
 end module case_waller_creek
