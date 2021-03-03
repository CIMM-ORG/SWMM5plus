! module test_cases
!
! Calling routines for custom test cases, e.g. calls the case_simple_channel
! functions to setup a single channel reach.
!
!==========================================================================
!
module trajkovic_cases
    !
    use allocate_storage
    use array_index
    use bc
    use control
    use data_keys
    use globals
    use setting_definition
    use utility
    use xsect_tables

    implicit none

    private

    public :: trajkovic_cases_setup
    public :: trajkovic_cases_initialize

    integer :: debuglevel = 1

contains
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine trajkovic_cases_setup &
        (init_depth, depth_dnstream, depth_upstream, lowerZ, upperZ, link_length, &
        link_breadth, subdivide_length, flowrate, area, velocity, Froude,         &
        ManningsN, idepth_type, link_geometry, inlet_offset, outlet_offset,       &
        full_depth, cDis1)

        character(64) :: subroutine_name = 'trajkovic_cases_setup'

        real, intent(inout) :: init_depth(:), depth_dnstream(:), depth_upstream(:)
        real, intent(inout) :: lowerZ(:), upperZ(:), link_length(:), link_breadth(:)
        real, intent(inout) :: subdivide_length(:), flowrate(:), area(:), velocity(:)
        real, intent(inout) :: Froude(:), ManningsN(:), full_depth(:), inlet_offset(:)
        real, intent(inout) :: cDis1(:), outlet_offset(:)

        integer, intent(inout) :: idepth_type(:)
        integer, intent(inout) :: link_geometry(:)

        real,    dimension(:), allocatable :: link_slope
        integer, dimension(:), allocatable :: link_type, subdivide_elements

        real :: CFL, ManningsNBuffer,total_length
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        allocate(link_slope(N_link))
        allocate(link_type(N_link))
        allocate(subdivide_elements(N_link))

        !% pipe u/s of gate 1          = 1.5m
        !% length of gate1             = 0.02m
        !% pipe from gate1 to gate2    = 8.26m
        !% length of gate2             = 0.02m
        !% pipe extension d/s of gate2 = 1.4m
        !% buffet pipe for dissipation = 5.0m
        total_length = 16.2

        ! set up flow and time step for differen subcases
        ! tests that ran:  Fr = 0.25, 0.5
        Froude       = 0.5   ! determines flowrate and slope to get Froude
        CFL          = 0.5   ! determines dt from subdivide_length

        idepth_type     = 1  
        ManningsN       = 0.008
        ManningsNBuffer = 0.016  ! ManningsN at the buffer link

        link_geometry       = lCircular
        link_slope          = 0.027
        depth_upstream      = 0.02
        depth_dnstream      = 0.02
        init_depth          = 0.02
        flowrate            = 0.0013
        inlet_offset        = 0.0
        outlet_offset       = 0.0
        full_depth          = 0.1
        cDis1               = 0.0

        !% case specific channel and subdivide lengths
        !% pipe upstream for gate 1
        link_length(1)        = 1.5
        link_type(1)          = lPipe
        subdivide_elements(1) = 6
        subdivide_length(1)   = link_length(1) / subdivide_elements(1)

        !% gate 1
        link_length(2)        = 0.02
        link_type(2)          = lOrifice
        subdivide_elements(2) = 1
        subdivide_length(2)   = link_length(2) / subdivide_elements(2)
        cDis1(2)              = 0.86

        !% pipe in between gate 1 and 2
        link_length(3)        = 8.26
        link_type(3)          = lPipe
        subdivide_elements(3) = 40
        subdivide_length(3)   = link_length(3) / subdivide_elements(3)

        !% gate 2
        link_length(4)        = 0.02
        link_type(4)          = lOrifice
        subdivide_elements(4) = 1
        subdivide_length(4)   = link_length(4) / subdivide_elements(4)
        cDis1(4)              = 0.86

        !% pipe downstream of gate 2
        link_length(5)        = 1.4
        link_type(5)          = lPipe
        subdivide_elements(5) = 6
        subdivide_length(5)   = link_length(5) / subdivide_elements(5)
        

        !% buffer pipe
        link_length(6)        = 5.0
        link_type(6)          = lPipe
        subdivide_elements(6) = 20
        subdivide_length(6)   = link_length(6) / subdivide_elements(6)
        link_slope(6)         = 0.0
        ManningsN(6)          = ManningsNBuffer
        full_depth(6)         = 0.5
        inlet_offset(6)       = -0.4
        outlet_offset(6)      = 0.0
 
        !% setup zbottom for links from downstream to upstream
        !% buffer pipe
        lowerZ(6) = -0.4
        upperZ(6) = -0.4

        !% pipe downstream of gate 2
        lowerZ(5) = 0.0
        upperZ(5) = lowerZ(5) + link_slope(5) * link_length(5)

        !% gate 2
        lowerZ(4) = upperZ(5)
        upperZ(4) = lowerZ(4) + link_slope(4) * link_length(4)

        !% pipe in between gate 1 and 2
        lowerZ(3) = upperZ(4)
        upperZ(3) = lowerZ(3) + link_slope(3) * link_length(3)

        !% gate 1
        lowerZ(2) = upperZ(3)
        upperZ(2) = lowerZ(2) + link_slope(2) * link_length(2)

        !% pipe upstream for gate 1
        lowerZ(1) = upperZ(2)
        upperZ(1) = lowerZ(1) + link_slope(1) * link_length(1)


        !% assign small values for trajkovic cases
        setting%SmallVolume%DepthCutoff      = 1E-03
        !% set minimum topwidth as 5% of pipe radius
        setting%SmallVolume%MinimumTopwidth  = 0.05 * onehalfR * full_depth(1)
        setting%SmallVolume%MinimumArea      = onehalfR * setting%SmallVolume%DepthCutoff * &
                                                          setting%SmallVolume%MinimumTopwidth
        setting%SmallVolume%MinimumPerimeter = setting%SmallVolume%MinimumTopwidth
        setting%SmallVolume%MinimumHydRadius = setting%SmallVolume%MinimumArea / setting%SmallVolume%MinimumPerimeter

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name
    end subroutine trajkovic_cases_setup
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine trajkovic_cases_initialize &
        (link_length, link_breadth, subdivide_length, lowerZ, upperZ,      &
        initial_flowrate, depth_upstream, depth_dnstream, initial_depth,   &
        inlet_offset, outlet_offset, discharge_coefficient1, full_depth,   &
        ManningsN, idepth_type, linkR, nodeR, linkI, nodeI,linkYN, nodeYN, &
        linkName, nodeName, bcdataDn, bcdataUp, gateSetting)
        !
        ! initialize the link-node system and boundary conditions for trajkovic_cases
        !
        character(64) :: subroutine_name = 'trajkovic_cases_initialize'

        real,  intent(in)  :: link_length(:), link_breadth(:), subdivide_length(:)
        real,  intent(in)  :: lowerZ(:), upperZ(:),  initial_flowrate(:)
        real,  intent(in)  :: depth_upstream(:), depth_dnstream(:), initial_depth(:)
        real,  intent(in)  :: inlet_offset(:), discharge_coefficient1(:)
        real,  intent(in)  :: full_depth(:), ManningsN(:), outlet_offset(:)

        integer, intent(in):: idepth_type(:)

        integer,   dimension(:,:), allocatable, target, intent(out)    :: linkI
        integer,   dimension(:,:), allocatable, target, intent(out)    :: nodeI

        real,      dimension(:,:), allocatable, target, intent(out)    :: linkR
        real,      dimension(:,:), allocatable, target, intent(out)    :: nodeR

        logical,   dimension(:,:), allocatable, target, intent(out)    :: linkYN
        logical,   dimension(:,:), allocatable, target, intent(out)    :: nodeYN

        type(string), dimension(:), allocatable, target, intent(out)   :: linkName
        type(string), dimension(:), allocatable, target, intent(out)   :: nodeName

        type(bcType), dimension(:), allocatable, intent(out) :: bcdataUp, bcdataDn
        type(controlType), dimension(:), allocatable, intent(out) :: gateSetting

        integer    :: ntimepoint, ndnstreamBC, nupstreamBC, nChanges
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        !% Boundary conditions
        ntimepoint = 2
        nupstreamBC = 1
        ndnstreamBC = 1

        call bc_allocate &
            (bcdataDn, bcdataUp, ndnstreamBC, nupstreamBC, ntimepoint)

        !% assign values
        !% downstream is default to elevation
        bcdataDn(1)%NodeID = 7
        bcdataDn(1)%TimeArray(1)     = setting%Time%StartTime
        bcdataDn(1)%TimeArray(2)     = setting%Time%EndTime + 100.0 !s
        bcdataDn(1)%ValueArray(1)    = 0.03   ! m !% from pipeAC. needs further investigation
        bcdataDn(1)%ValueArray(2)    = 0.03   ! m

        ! upstream is default to flowrate
        bcdataUp(1)%NodeID = 1
        bcdataUp(1)%TimeArray(1)  = setting%Time%StartTime
        bcdataUp(1)%TimeArray(2)  = setting%Time%EndTime + 100.0 !s
        bcdataUp(1)%ValueArray(1) = initial_flowrate(1)  ! m^3/s
        bcdataUp(1)%ValueArray(2) = initial_flowrate(1)  ! m^3/2

        !% assign control values
        N_Gates  = 2
        nChanges = 2

        call control_allocate (gateSetting, N_Gates, nChanges)

        !% assign control values to gate 1
        gateSetting(1)%LinkId              = 2
        gateSetting(1)%HeightStart         = 0.014
        gateSetting(1)%HeightNow           = 0.0
        gateSetting(1)%AreaNow             = 0.0
        gateSetting(1)%AreaPrior           = 0.0
        gateSetting(1)%GateTimeChange1     = nullvalueR
        gateSetting(1)%GateTimeChange2     = nullvalueR
        gateSetting(1)%GateHeightChange1   = 0.014
        gateSetting(1)%GateHeightChange2   = 0.014
        gateSetting(1)%HeightMinimum       = 1E-6
        gateSetting(1)%GateSpeed           = 0.01
        gateSetting(1)%CanMove             = .false.
        gateSetting(1)%MovedThisStep       = .false.

        !% assign control values to gate 2
        gateSetting(2)%LinkId              = 4
        gateSetting(2)%HeightStart         = 0.1
        gateSetting(2)%HeightNow           = 0.0
        gateSetting(2)%AreaNow             = 0.0
        gateSetting(2)%AreaPrior           = 0.0
        gateSetting(2)%GateTimeChange1     = 120.0
        gateSetting(2)%GateTimeChange2     = 150.0
        gateSetting(2)%GateHeightChange1   = 1E-6
        gateSetting(2)%GateHeightChange2   = 0.028
        gateSetting(2)%HeightMinimum       = 1E-6
        gateSetting(2)%GateSpeed           = 0.01
        gateSetting(2)%CanMove             = .true.
        gateSetting(2)%MovedThisStep       = .false.

        call trajkovic_cases_link_node &
            (link_length, link_breadth, subdivide_length, lowerZ, upperZ,      &
            initial_flowrate, depth_upstream, depth_dnstream, initial_depth,   &
            inlet_offset, outlet_offset, discharge_coefficient1, full_depth,   &
            ManningsN, idepth_type, linkR, nodeR, linkI, nodeI,linkYN, nodeYN, &
            linkName, nodeName)

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine trajkovic_cases_initialize
    !
    !==========================================================================
    !
    ! PRIVATE BELOW HERE
    !
    !==========================================================================
    !
    subroutine trajkovic_cases_link_node &
        (link_length, link_breadth, subdivide_length, lowerZ, upperZ,      &
        initial_flowrate, depth_upstream, depth_dnstream, initial_depth,   &
        inlet_offset, outlet_offset, discharge_coefficient1, full_depth,   &
        ManningsN, idepth_type, linkR, nodeR, linkI, nodeI,linkYN, nodeYN, &
        linkName, nodeName)
        !
        ! creates the link-node system for trajkovic_cases
        !
        character(64) :: subroutine_name = 'trajkovic_cases_link_node'

        real,  intent(in)  :: link_length(:), link_breadth(:), subdivide_length(:)
        real,  intent(in)  :: lowerZ(:), upperZ(:),  initial_flowrate(:)
        real,  intent(in)  :: depth_upstream(:), depth_dnstream(:), initial_depth(:)
        real,  intent(in)  :: inlet_offset(:), discharge_coefficient1(:)
        real,  intent(in)  :: full_depth(:), ManningsN(:), outlet_offset(:)

        integer, intent(in):: idepth_type(:)

        integer,   dimension(:,:), allocatable, target, intent(out)    :: linkI
        integer,   dimension(:,:), allocatable, target, intent(out)    :: nodeI

        real,      dimension(:,:), allocatable, target, intent(out)    :: linkR
        real,      dimension(:,:), allocatable, target, intent(out)    :: nodeR

        logical,   dimension(:,:), allocatable, target, intent(out)    :: linkYN
        logical,   dimension(:,:), allocatable, target, intent(out)    :: nodeYN

        type(string), dimension(:), allocatable, target, intent(out)   :: linkName
        type(string), dimension(:), allocatable, target, intent(out)   :: nodeName

        integer :: ii, mm
        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        call allocate_linknode_storage &
            (linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName)

        ! assign the indexes
        linkI(:,li_idx) = (/ (ii, ii=1,N_link) /)
        nodeI(:,ni_idx) = (/ (ii, ii=1,N_node) /)

        ! assign no names for links
        linkName(1)%str = 'Pipe up'
        linkName(2)%str = 'Gated Orifice 1'
        linkName(3)%str = 'Pipe middle'
        linkName(4)%str = 'Gated Orifice 2'
        linkName(5)%str = 'Pipe dn'
        linkName(6)%str = 'Pipe buffer'

        ! assign zeros for accumulators
        nodeI(:,ni_N_link_d) = 0
        nodeI(:,ni_N_link_u) = 0

        ! assign orifice type
        linkI(2,li_orif_type) = lSideOrifice
        linkI(4,li_orif_type) = lSideOrifice

        ! assign uniform physical data
        linkI(:,li_roughness_type)  = lManningsN
        linkR(:,lr_Roughness)       = ManningsN

        ! designate the upstream nodes
        nodeI(1,ni_node_type) = nBCup
        nodeR(1,nr_Zbottom)   = upperZ(1)
        nodeName(1)%str = 'UpStreamBC'

        do ii = 2,N_node-1
            nodeI(ii,ni_node_type) = nJ2
            nodeR(ii,nr_Zbottom)   = lowerZ(ii-1) 
            ! print*, ii, nodeR(ii,nr_Zbottom) , sum(channel_length(1:N_link)) - (ii-1)* channel_length(ii)
        end do

        ! designate the downstream node
        nodeI(N_node,ni_node_type) = nBCdn
        nodeR(N_node,nr_Zbottom)   = lowerZ(6)
        nodeName(N_node)%str = 'DownstreamBC'

        ! assign all as rectangular channels
        linkI(:,li_geometry)  = lCircular
        linkI(:,li_link_type) = lPipe
        linkI(2,li_link_type) = lOrifice
        linkI(4,li_link_type) = lOrifice

        ! assign the link position and mappings
        do mm = 1, N_link
            linkI(mm,li_Mnode_u) = mm     ! map to upstream node
            linkI(mm,li_Mnode_d) = mm + 1 ! map to downstream node
            linkR(mm,lr_Length)          = link_length(mm)
            linkR(mm,lr_BreadthScale)    = full_depth(mm)
            linkR(mm,lr_ElementLength)   = subdivide_length(mm)
            linkR(mm,lr_InitialFlowrate) = initial_flowrate(mm)
            linkR(mm,lr_InitialDepth)    = initial_depth(mm)
            linkR(mm,lr_InitialUpstreamDepth)    = depth_upstream(mm)
            linkR(mm,lr_InitialDnstreamDepth)    = depth_dnstream(mm)
            linkR(mm,lr_InletOffset)             = inlet_offset(mm)
            linkR(mm,lr_OutletOffset)            = outlet_offset(mm)
            linkR(mm,lr_DischargeCoeff1)         = discharge_coefficient1(mm)
            linkR(mm,lr_FullDepth)               = full_depth(mm)
            linkI(mm,li_InitialDepthType)        = idepth_type(mm)
        enddo

        if ((debuglevel > 0) .or. (debuglevelall > 0)) then
            print *
            print *, subroutine_name,'-----------------------------------'
            print *, 'link info'
            print *, linkI(:,li_idx), ' idx'
            print *, linkI(:,li_link_type), ' type'
            print *, linkI(:,li_Mnode_u) , ' upstream node'
            print *, linkI(:,li_Mnode_d) , ' downstream node'
            print *, linkR(:,lr_InitialDnstreamDepth), 'downstream depth'
            print *, linkR(:,lr_InitialUpstreamDepth), 'upstream depth'
            print *, linkR(:,lr_InletOffset), ' inletOffset'
            print *, linkR(:,lr_OutletOffset), ' outletOffset'
            print *, linkR(:,lr_Length), 'length' 
            print *, 'node info'
            print *, nodeI(:,ni_idx), ' idx'
            print *, nodeI(:,ni_node_type), ' type'
            print *, nodeR(:,nr_Zbottom),   ' zbottom'
            ! print *, nodeI(:,ni_N_link_d), 'number of downstream links'
            ! print *, nodeI(:,ni_Mlink_d1), 'downstream1 link'
            ! print *, nodeI(:,ni_N_link_u), 'number of upstream links'
            ! print *, nodeI(:,ni_Mlink_u1), 'upstream1 link'
        endif
        
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine trajkovic_cases_link_node
    !
    !==========================================================================
    ! END OF MODULE trajkovic_cases
    !==========================================================================
end module trajkovic_cases
