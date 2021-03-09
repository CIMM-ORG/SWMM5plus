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

    integer :: debuglevel = 1

contains
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine case_waller_creek_initialize &
        (channel_length, channel_breadth, channel_topwidth, subdivide_length, &
        lowerZ, initial_flowrate, init_depth, depth_upstream, &
        depth_dnstream, ManningsN, roughness_type, idepth_type,             &
        linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName,     &
        bcdataDn, bcdataUp, &
        wdID, widthDepthData)
        !
        ! initialize the link-node system and boundary conditions for a simple channel
        !
        character(64) :: subroutine_name = 'case_waller_creek_initialize'

        real(8),  intent(in)  :: channel_length(:), channel_breadth(:)
        real(8),  intent(in)  :: channel_topwidth(:), subdivide_length(:)
        real(8),  intent(in)  :: lowerZ(:), initial_flowrate(:)
        real(8),  intent(in)  :: depth_upstream(:), depth_dnstream(:), init_depth(:)
        real(8),  intent(in)  :: ManningsN(:)

        integer, intent(in):: roughness_type, idepth_type(:)

        integer, target, intent(in out)    :: wdID(:)
        real(8),    target, intent(in out)    :: widthDepthData(:,:,:)
        !integer, target, intent(in out)    :: wdnumberPairs(:)
        !real(8),    target, intent(in out)    :: wdxDistance(:)
        !type(string), target, intent(in out)   :: wdcellType(:)

        integer,   dimension(:,:), allocatable, target, intent(out)    :: linkI
        integer,   dimension(:,:), allocatable, target, intent(out)    :: nodeI

        real(8),      dimension(:,:), allocatable, target, intent(out)    :: linkR
        real(8),      dimension(:,:), allocatable, target, intent(out)    :: nodeR

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
        ! upstream is default to flowrate
        bcdataUp(1)%NodeID = 1
        bcdataUp(1)%TimeArray(1)  = setting%Time%StartTime
        bcdataUp(1)%TimeArray(2)  = setting%Time%EndTime + 1000.0 !s
        bcdataUp(1)%ValueArray(1) = initial_flowrate(1)  ! m^3/s
        bcdataUp(1)%ValueArray(2) = initial_flowrate(1)  ! m^3/s

        ! downstream is default to elevation
        bcdataDn(1)%NodeID = 567
        bcdataDn(1)%TimeArray(1)     = setting%Time%StartTime
        bcdataDn(1)%TimeArray(2)     = setting%Time%EndTime + 1000.0 !s
        bcdataDn(1)%ValueArray(1)    = lowerZ(size(lowerZ)) +  depth_dnstream(size(depth_dnstream)) ! m
        bcdataDn(1)%ValueArray(2)    = lowerZ(size(lowerZ)) +  depth_dnstream(size(depth_dnstream)) ! m

        print *, "lowerZ = ", lowerZ(size(lowerZ))
        !call test_case_initial_condition_setup (linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName, &
        !    bcdataUp, bcdataDn, channel_length, channel_breadth)

        call case_waller_creek_links_and_nodes &
            (channel_length, channel_breadth, channel_topwidth, subdivide_length, &
            lowerZ, initial_flowrate, depth_upstream, &
            depth_dnstream, ManningsN, roughness_type,  idepth_type, &
            linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName, &
            wdID, widthDepthData, bcdataUp, bcdataDn)

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine case_waller_creek_initialize
    !
    !==========================================================================
    !
    ! PRIVATE BELOW HERE
    !
    !==========================================================================
    !
    subroutine case_waller_creek_links_and_nodes &
        (channel_length, channel_breadth, channel_topwidth, subdivide_length, &
        lowerZ, initial_flowrate, depth_upstream, &
        depth_dnstream, ManningsN, roughness_type, idepth_type, &
        linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName, &
        wdID, widthDepthData, bcdataUp, bcdataDn)

        character(64) :: subroutine_name = 'case_waller_creek_links_and_nodes'

        real(8),  intent(in)  :: channel_length(:), channel_breadth(:)
        real(8),  intent(in)  :: channel_topwidth(:), subdivide_length(:)
        real(8),  intent(in)  :: lowerZ(:), ManningsN(:), initial_flowrate(:)
        real(8),  intent(in)  :: depth_upstream(:), depth_dnstream(:)!, init_depth(:)

        integer, intent(in):: roughness_type, idepth_type(:)

        integer, target, intent(in out)    :: wdID(:)
        real(8),    target, intent(in out)    :: widthDepthData(:,:,:)
        !integer, target, intent(in out)    :: wdnumberPairs(:)
        !real(8),    target, intent(in out)    :: wdxDistance(:)
        !type(string), target, intent(in out)   :: wdcellType(:)

        integer,   dimension(:,:), allocatable, target, intent(out)    :: linkI
        integer,   dimension(:,:), allocatable, target, intent(out)    :: nodeI

        real(8),      dimension(:,:), allocatable, target, intent(out)    :: linkR
        real(8),      dimension(:,:), allocatable, target, intent(out)    :: nodeR

        logical,   dimension(:,:), allocatable, target, intent(out)    :: linkYN
        logical,   dimension(:,:), allocatable, target, intent(out)    :: nodeYN

        type(string), dimension(:), allocatable, target, intent(out)   :: linkName
        type(string), dimension(:), allocatable, target, intent(out)   :: nodeName
        
        type(bcType), intent(in out) :: bcdataUp(:), bcdataDn(:) !for setting up the depth at each link
        integer :: mm, ii


        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        N_link = wdID(size(wdID)) !wdID is the newID acquired from the text reader
        N_node = N_link + 1

        call allocate_linknode_storage &
            (linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName)

        ! assign the indexes
        linkI(:,li_idx) = (/ (ii, ii=1,N_link) /)
        nodeI(:,ni_idx) = (/ (ii, ii=1,N_node) /)
        
        ! assign names for links
        do ii=1,N_link
            linkName(ii)%str = 'widthdepth_pair'
        end do

        ! Justin: Define the connectivity here. Need to modify it and move it to somewhere else
        ! assign zeros for accumulators
        nodeI(:,ni_N_link_d) = 0
        nodeI(:,ni_N_link_u) = 0

        ! assign uniform physical data
        linkI(:,li_roughness_type)  = roughness_type
        linkR(:,lr_Roughness)       = ManningsN
        
        
        nodeI(1,ni_node_type) = nBCup

        nodeR(1,nr_Zbottom) = lowerZ(1) + (lowerZ(1) - lowerZ(2))

        nodeName(1)%str = 'UpstreamBC'

        do ii=2,N_node
            nodeI(ii,ni_node_type) = nJ2

            nodeR(ii,nr_Zbottom) = lowerZ(ii-1)

            nodeName(ii)%str = 'Junction'
        end do

        ! designate the downstream node
        nodeI(N_node,ni_node_type) = nBCdn

        nodeName(N_node)%str = 'DownstreamBC'

        ! assign the link types
        linkI(:,li_link_type) = lChannel

        ! assign all as rectangular channels
        linkI(:,li_geometry) = lWidthDepth

        ! assign the link position and mappings

        do ii=1,N_link
            linkI(ii,li_Mnode_u) = ii 
            linkI(ii,li_Mnode_d) = ii +1
        end do

        do mm=1,N_link
            linkR(mm,lr_Length)          = channel_length(mm)
            linkR(mm,lr_BreadthScale)    = channel_breadth(mm)
            linkR(mm,lr_TopWidth)        = channel_topwidth(mm)
            linkR(mm,lr_ElementLength)   = subdivide_length(mm)
            linkR(mm,lr_InitialFlowrate) = initial_flowrate(mm)
            linkI(mm,li_InitialDepthType)= idepth_type(mm)
        enddo
        
        call test_case_initial_depth_setup(wdID, linkR, nodeR, linkI, nodeI, linkName, nodeName, &
            bcdataUp, bcdataDn, ManningsN, widthDepthData)

        !linkR(:  ,lr_InitialDepth)         = init_depth(:)
        !linkR(:  ,lr_InitialDnstreamDepth) = depth_dnstream(:) !Justin: Need to change
        !linkR(:  ,lr_InitialUpstreamDepth) = depth_upstream(:) !Justin: Need to change


        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine case_waller_creek_links_and_nodes
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine test_case_initial_depth_setup &
        (wdID, linkR, nodeR, linkI, nodeI, linkName, nodeName, &
        bcdataUp, bcdataDn, ManningsN, widthDepthData)
        character(64) :: subroutine_name = 'test_case_initial_depth_setup'
        
        integer, intent(in)    :: wdID(:)   
        integer,   target, intent(in out)    :: linkI(:,:)
        integer,   target, intent(in out)    :: nodeI(:,:)

        real(8),      target, intent(in out)    :: linkR(:,:)
        real(8),      target, intent(in out)    :: nodeR(:,:)

        real(8),  intent(in)  :: ManningsN(:)

        type(string), target, intent(out)   :: linkName(:)
        type(string), target, intent(out)   :: nodeName(:)

        real(8),    target, intent(in out)    :: widthDepthData(:,:,:)
        type(bcType), intent(in out) :: bcdataUp(:), bcdataDn(:)

        real(8), pointer :: areaTotalBelowThisLayer(:,:)
        real(8), pointer :: perimeterBelowThisLayer(:,:)
        real(8), pointer :: depthTBL(:,:)
        real(8), dimension(:), allocatable    :: temp !temp array to store A**(5/3)/P**(2/3)

        integer :: mm, idx 
        real(8)    :: flowrate, temp_manning, temp_slope
        real(8)    :: ME_RHS !Manning's eq RHS
        
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name
        
        !Note: This is just for single channel, we need better a subroutine to describe the cumulative flowrate from upstream to downstream 
        N_link = size(wdID,1)

        if (size(bcdataUp%NodeID,1) == 1) then
            flowrate = bcdataUp(1)%ValueArray(1) !Use the first timestep flowrate to solve the initial condition
        else
            print *, "For the time being, the initial condition computation for link only supports single channel test case."
        endif

        ! Use that flowrate at upstream BC to estimate the initial depth
        ! Roughly estimation by using Manning's equation

        areaTotalBelowThisLayer => widthDepthData (:,:, wd_areaTotalBelowThisLayer)
        depthTBL                => widthDepthData (:,:, wd_depthTotalBelowThisLayer)
        perimeterBelowThisLayer => widthDepthData (:,:, wd_perimeterBelowThisLayer)
        !rule: widthDepthData(link_ID, Layers, Variable)

        ! In this loop, find the layer (A**(5/3))/(P**(2/3)) closest to (n*Q)/(S**0.5), and use the corresponding depth at the layer as the initial depth

        allocate(temp(size(areaTotalBelowThisLayer,2)), source=nullvalueR) !initialize the temp array with nullvalueR (-998877)

        do mm=1,N_link-1

            temp_manning    = linkR(mm,lr_Roughness) !Local roughness at the link
            temp_slope      = (nodeR(mm,nr_Zbottom) - nodeR(mm+1,nr_Zbottom)) / linkR(mm,lr_Length)

            if ( temp_slope <=0 )  print *, "Link=", mm, " Local slope is adverse, ", temp_slope
            
            ME_RHS          = temp_manning * flowrate / temp_slope**0.5
            temp            = (areaTotalBelowThisLayer(mm,:)**(5.0/3.0)) / (perimeterBelowThisLayer(mm,:)**(2.0/3.0))
            idx             = minloc(abs(temp - ME_RHS), DIM=1)
            ! find the index of the closest layer 

            linkR(mm  ,lr_InitialDepth)         = depthTBL(mm, idx) 
            linkR(mm  ,lr_InitialDnstreamDepth) = depthTBL(mm, idx)
            linkR(mm  ,lr_InitialUpstreamDepth) = depthTBL(mm, idx)
            !print *, "Link = ", mm, "Initial Depth = ", linkR(mm, lr_InitialDepth)
        enddo

        !linkR(mm  ,lr_InitialDepth)         = bcdataDn(1)%ValueArray(1)
        !linkR(mm  ,lr_InitialDnstreamDepth) = bcdataDn(1)%ValueArray(1)
        !linkR(mm  ,lr_InitialUpstreamDepth) = bcdataDn(1)%ValueArray(1)

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name

    end subroutine test_case_initial_depth_setup


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
        
        character(64) :: subroutine_name = 'nonmonotonic_subdivide'


        integer, intent(inout) :: NX
        integer, intent(inout) :: max_number_of_pairs
        integer :: newNX = 0
        integer :: ii, jj, mm

        integer :: allocation_status
        character(len=99) :: emsg

        integer, intent(in) :: ID(:)
        integer, intent(in) :: numberPairs(:)
        real(8), intent(in)    :: ManningsN(:)
        real(8), intent(in)    :: Length(:)
        real(8), intent(in)    :: zBottom(:)
        real(8), intent(in)    :: xDistance(:)
        real(8), intent(in)    :: Breadth(:)
        real(8), intent(in)    :: widthDepthData(:,:,:)
        type(string), intent(in out)   :: cellType(:)

        real(8), intent(in)    :: subdivide_length_check

        real(8),    dimension(:),      allocatable :: faceZBottom
        real(8),    dimension(:),      allocatable :: temp1
        real(8),    dimension(:),      allocatable :: temp2
        integer, dimension(:),      allocatable :: isnonmonotonic

        integer, dimension(:),      allocatable :: newID
        integer, dimension(:),      allocatable :: newNumberPairs
        real(8),    dimension(:),      allocatable :: newManningsN
        real(8),    dimension(:),      allocatable :: newLength
        real(8),    dimension(:),      allocatable :: newZBottom
        real(8),    dimension(:),      allocatable :: newXDistance
        real(8),    dimension(:),      allocatable :: newBreadth
        real(8),    dimension(:,:,:),  allocatable :: newWidthDepthData
        type(string), dimension(:), allocatable :: newCellType(:)

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        allocate(isnonmonotonic(NX))
        isnonmonotonic= 0

        allocate(faceZBottom(NX+1))
        faceZbottom(:) = 0.0

        ! find the z at the faces
        call face_zbottom(faceZbottom, zbottom, Length, NX)

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

        allocate(newCellType(newNX), stat=allocation_status, errmsg=emsg)
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
                if (Length(ii) <= threeR * subdivide_length_check) then ! modify here, the original way to subdivide can cause negative length value
                    newLength   (jj)     =  Length(ii) / threeR
                    newLength   (jj+1)   =  Length(ii) / threeR
                    newLength   (jj+2)   =  Length(ii) / threeR

                    newXDistance(jj)     = xDistance(ii) - twoR * Length(ii) / threeR
                    newXDistance(jj+1)   = xDistance(ii) - oneR * Length(ii) / threeR
                    newXDistance(jj+2)   = xDistance(ii)
                    
                else    ! length is longer than 3 sub
                    newLength      (jj)   = subdivide_length_check
                    newLength      (jj+1) = Length(ii) - twoR * subdivide_length_check
                    newLength      (jj+2) = subdivide_length_check
                
                    newXDistance   (jj)   = xDistance(ii) - Length(ii) + oneR * subdivide_length_check
                    newXDistance   (jj+1) = xDistance(ii) - oneR * subdivide_length_check
                    newXDistance   (jj+2) = xDistance(ii)
                endif
                newZBottom     (jj)   = faceZbottom (ii)    - onehalfR * (faceZbottom (ii) - zBottom (ii))
                newZBottom     (jj+1) = zBottom     (ii)
                newZBottom     (jj+2) = zBottom     (ii)    - onehalfR * (zBottom (ii) - faceZbottom (ii))
                
                !newXDistance   (jj)   = onehalfR * xDistance (ii) - onehalfR * Length (ii)      - onehalfR * subdivide_length_check
                !newXDistance   (jj+1) = newXDistance   (jj)     + subdivide_length_check
                !newXDistance   (jj+2) = newXDistance   (jj+1)   + onehalfR * Length (ii)      - onehalfR * subdivide_length_check
                
                newBreadth     (jj)   = Breadth     (ii)
                newBreadth     (jj+1) = Breadth     (ii)
                newBreadth     (jj+2) = Breadth     (ii)
                newWidthDepthData(jj,:,:)   = widthDepthData (ii,:,:)
                newWidthDepthData(jj+1,:,:) = widthDepthData (ii,:,:)
                newWidthDepthData(jj+2,:,:) = widthDepthData (ii,:,:)
                newCellType      (jj)%str   = cellType    (ii)%str
                newCellType      (jj+1)%str = cellType    (ii)%str
                newCellType      (jj+2)%str = cellType    (ii)%str
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
                newCellType      (jj)%str = cellType       (ii)%str
                jj = jj + 1
            endif
        enddo

        if ((debuglevel > 0) .or. (debuglevelall > 0)) then
            print *, "NEW)  ID,    LENGTH,    Z,    XDistance"
            do mm = 1, NX
                print *,  newID(mm), newLength(mm), newZBottom(mm), newXDistance(mm)
            enddo
        endif

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

        real(8), intent(inout) :: faceZbottom(:)
        real(8), intent(in)    :: zBottom(:)
        real(8), intent(in)    :: Length(:)
        integer, intent(in) :: NX
        integer :: ii

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        faceZbottom (2:NX) = &
            (zBottom(1:NX-1)*Length(2:NX) + zBottom(2:NX)*Length(1:NX-1)) &
            /(Length(2:NX) + Length(1:NX-1))

        faceZbottom (1)  = zBottom(1)
        faceZbottom (NX+1) = zBottom(NX)

        print *, "zBottom", "      faceZbottom"
        do ii=1,NX
            print *, zBottom(ii), faceZbottom(ii)
        enddo

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine face_zbottom
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine widthdepth_pair_auxiliary (widthDepthData, cellType, numberPairs)

        character(64) :: subroutine_name = 'widthdepth_pair_auxiliary'

        integer, intent(inout) :: numberPairs(:)

        real(8), target, intent(inout) :: widthDepthData(:,:,:)

        type(string), intent(in out)   :: cellType(:)

        real(8), pointer :: width(:,:), depth(:,:), area(:,:), areaTBL(:,:)
        real(8), pointer :: dWidth(:,:), dDepth(:,:), angle(:,:), perimeterBL(:,:)
        real(8), pointer :: rH(:,:), gammaBTL(:,:), depthTBL(:,:)
        real(8), pointer :: area_difference(:,:), local_difference(:,:)

        integer :: ii,jj, kk
        integer :: eIn1

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        eIn1 = size(widthDepthData,2) ! JY: The number of read-in pairs

        width       => widthDepthData (:,:, wd_widthAtLayerTop)
        depth       => widthDepthData (:,:, wd_depthAtLayerTop)
        area        => widthDepthData (:,:, wd_areaThisLayer)
        areaTBL     => widthDepthData (:,:, wd_areaTotalBelowThisLayer)
        depthTBL    => widthDepthData (:,:, wd_depthTotalBelowThisLayer)
        dWidth      => widthDepthData (:,:, wd_Dwidth)
        dDepth      => widthDepthData (:,:, wd_Ddepth)
        angle       => widthDepthData (:,:, wd_angle)
        perimeterBL => widthDepthData (:,:, wd_perimeterBelowThisLayer)
        rH          => widthDepthData (:,:, wd_rHBelowThisLayer)
        gammaBTL    => widthDepthData (:,:, wd_gammaBelowThisLayer)
        area_difference  => widthDepthData (:,:, wd_area_difference)
        local_difference => widthDepthData (:,:, wd_local_difference)

        ! lowest layer is triangular
        area(:,1) = onehalfR * width(:,1) * depth(:,1) !Justin: I'm not sure if this is ok, could cause some discrepancies
        !print *, "width", width(1,:)
        ! the area in this layer
        area(:,2:eIn1) = onehalfR * (width(:,2:eIn1) + width(:,1:eIn1-1)) &
            * (depth(:,2:eIn1) - depth(:,1:eIn1-1)) !Justin: Assuming each section is trapezoidal shape

        ! set areas to zero above the uppermost pair
        do ii = 1, size(widthDepthData,1)
            if(cellType(ii)%str == 'widthdepth_pair') then
                area(ii, numberPairs(ii):eIn1) = zeroR
            endif
        enddo

        ! store width and depth differences
        !dWidth(:,1) = width(:,1) !Justin: This is not right
        dDepth(:,1) = depth(:,1)
        
        ! delta width between top and bottom of this layer
        dWidth(:,1:eIn1-2) = width(:,2:eIn1-1) - width(:,1:eIn1-2)
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
        !areaTBL(:,2) = area(:,2)
        !areaTBL(:,3:eIn1) = areaTBL(:,2:eIn1-1) + area(:,2:eIn1-1)
        areaTBL(:,2:eIn1) = areaTBL(:,1:eIn1-1) + area(:,1:eIn1-1)
        !print *, "area size 1 = ", size(area,1)
        !print *, "area size 2 = ", size(area,2)
        do ii=1,size(area,1)
            do kk= 2, eIn1
            widthDepthData (ii,kk, wd_areaTotalBelowThisLayer) = sum(widthDepthData (ii,1:kk, wd_areaThisLayer))
            enddo
        enddo
        !print *, widthDepthData (1,:, wd_areaTotalBelowThisLayer)
        

        !print *, "eTn1 = ", eIn1
        !do kk=2,eIn1
        !    areaTBL(:,kk) =  sum(area(:,1:kk),2)
        !enddo
        !print *, "area", areaTBL(1,:)
        !print *,"area   = ", size(widthDepthData (1,:, wd_areaThisLayer))
        !print *,"areaTBL= ", size(widthDepthData (1,:, wd_areaTotalBelowThisLayer))

        


        depthTBL(:,1) = depth(:,1)
        depthTBL(:,2:eIn1) = depthTBL(:,1:eIn1-1) + depth(:,1:eIn1-1)
        ! check that the setting maximum area value is greater than any accumulated
        ! area at the uppermost level.
        if (setting%Method%AdjustWidthDepth%areaMaximum < maxval(areaTBL(:, eIn1))) then
            setting%Method%AdjustWidthDepth%areaMaximum = twoR * maxval(areaTBL(:, eIn1))
        endif


        !print *, "dDepth",dDepth(1,:)
        !print *, "dWidth", dWidth(1,:)
        ! perimeter below this layer
        !perimeterBL(:,1) = width(:,1)
        !perimeterBL(:,2:eIn1) = perimeterBL(:,1:eIn1-1) &
        !    + twoR * sqrt(dDepth(:,1:eIn1-1)**twoR + onehalfR * dWidth(:,1:eIn1-1)**twoR)
        widthDepthData (:,1, wd_perimeterBelowThisLayer) = width(:,1)
        ! Justin: I didn't use pointers in the for loop due to overflow issue
        do ii=1,size(width,1)
            do kk=2,eIn1
                widthDepthData (ii,kk, wd_perimeterBelowThisLayer) = widthDepthData (ii,kk-1, wd_perimeterBelowThisLayer) &
                    + twoR * sqrt(widthDepthData (ii,kk, wd_Ddepth)**twoR + (onehalfR * widthDepthData (ii,kk, wd_Dwidth))**twoR)
            enddo
        enddo

        !print *, "areaTBL", areaTBL(1,:)
        rH(:,1) = zeroR
        rH(:,2:eIn1) = areaTBL(:,2:eIn1)/perimeterBL(:,2:eIn1)
        !print *, "rH", rH(1,:)
        gammaBTL = areaTBL * rH**twothirdR
        !print *, "gammaBTL", gammaBTL(1,:)
        area_difference = area - areaTBL
        local_difference = area_difference - area

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine widthdepth_pair_auxiliary
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine widthdepth_pair_consistency (widthDepthData, numberPairs)

        character(64) :: subroutine_name = 'widthdepth_pair_consistency'

        integer, intent(inout) :: numberPairs(:)

        real(8), target, intent(inout) :: widthDepthData(:,:,:)

        real(8), dimension(:,:), allocatable :: dWidth
        real(8), dimension(:,:), allocatable :: dDepth

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

        real(8), target, intent(inout) :: widthDepthData(:,:,:)

        real(8), pointer :: up2W(:,:), up1W(:,:), lowW(:,:)
        real(8), pointer :: up2D(:,:), up1D(:,:), lowD(:,:)

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
    ! END OF MODULE case_waller_creek
    !==========================================================================
end module case_waller_creek
