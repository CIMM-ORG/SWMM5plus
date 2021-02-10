!==========================================================================
!
program main
    use allocate_storage
    use array_index
    use bc
    use checking
    use data_keys
    use debug
    use diagnostic
    use globals
    use initialization
    use initial_condition
    use junction
    use inflow
    use interface
    use network_graph
    use network_define
    use output
    use setting_definition
    use type_definitions
    ! use test_cases
    use time_loop
    use utility

    implicit none

    !%  elem2# are the values for elements that have only 2 faces
    real,       dimension(:,:), allocatable, target    :: elem2R       ! real data for elements with 2 faces
    integer,    dimension(:,:), allocatable, target    :: elem2I       ! integer data for elements with 2 faces
    logical,    dimension(:,:), allocatable, target    :: elem2YN      ! logical data for elements with 2 faces

    type(string), dimension(:), allocatable, target    :: elem2Name    ! array of character strings

    !%  elemM# are the values for elements that have more than 2 faces
    real,       dimension(:,:), allocatable, target    :: elemMR       ! real data for elements with multi faces
    integer,    dimension(:,:), allocatable, target    :: elemMI       ! integer data for elements with multi faces
    logical,    dimension(:,:), allocatable, target    :: elemMYN      ! logical data for elements with multi faces

    type(string), dimension(:), allocatable, target    :: elemMName    ! array of character strings

    !%  face# are the values for faces (always bounded by 2 elements)
    real,       dimension(:,:), allocatable, target    :: faceR       ! real data for faces
    integer,    dimension(:,:), allocatable, target    :: faceI       ! integer data for faces
    logical,    dimension(:,:), allocatable, target    :: faceYN      ! logical data for face

    type(string), dimension(:), allocatable, target    :: faceName    ! array of character strings

    !%  links are the building blocks from SWMM link-node formulation
    real,       dimension(:,:), allocatable, target    :: linkR       ! real data for links
    integer,    dimension(:,:), allocatable, target    :: linkI       ! integer data for links
    logical,    dimension(:,:), allocatable, target    :: linkYN      ! logical data for links

    type(string), dimension(:), allocatable, target    :: linkName    ! array of character strings

    !%  nodes are the building blocks from teh SWMM link-node formulation
    real,       dimension(:,:), allocatable, target    :: nodeR       ! real data for nodes
    integer,    dimension(:,:), allocatable, target    :: nodeI       ! integer data for nodes
    logical,    dimension(:,:), allocatable, target    :: nodeYN      ! logical data for nodes

    type(string), dimension(:), allocatable, target    :: nodeName   ! array of character strings

    !%  bcdata are structures containing boundary condition data
    type(bcType), dimension(:), allocatable :: bcdataUp, bcdataDn

    !%  debug output file information
    type(debugfileType),  dimension(:),   allocatable :: debugfile

    !%  diagnostic information
    type(diagnosticType), allocatable, dimension(:)   :: diagnostic

    !% threaded output files
    type(threadedfileType), allocatable, dimension(:) :: threadedfile

    integer, dimension(:),      allocatable :: wdID
    integer, dimension(:),      allocatable :: wdNumberPairs
    real,    dimension(:),      allocatable :: wdManningsN
    real,    dimension(:),      allocatable :: wdLength
    real,    dimension(:),      allocatable :: wdZBottom
    real,    dimension(:),      allocatable :: wdXDistance
    real,    dimension(:),      allocatable :: wdBreadth
    real,    dimension(:,:,:),  allocatable :: wdWidthDepthData
    type(string), dimension(:), allocatable :: wdCellType(:)

    type(graph) :: swmm_graph
    type(tseries) :: ts_ups
    integer :: ii, jj

    !--------------------------------------------------------------------------
    print *, ''
    print *, '====================================================='
    print *, 'starting main program'
    print *, ''

    !%  simulation controls
    call setting_default

    !===========================================================
    !%  hard-code setting for test cases

    setting%TestCase%UseTestCase = .false.
    ! setting%TestCase%TestName = 'simple_channel_001'
    ! setting%TestCase%TestName = 'y_channel_002'
    ! setting%TestCase%TestName = 'simple_weir_003'
    ! setting%TestCase%TestName = 'simple_orifice_004'
    ! setting%TestCase%TestName = 'y_storage_channel_005'
    ! setting%TestCase%TestName = 'waller_creek'

    !%  hard-code for debug output
    setting%Debugout%SuppressAllFiles  = .true. ! use this to easily suppress debug files

    setting%Debugout%SuppressTimeStep  = .true. ! use the next 3 to suppress headers
    setting%Debugout%SuppressTimeValue = .true. ! which can make debug files easier
    setting%Debugout%SuppressNdat      = .true. ! to read (but less useful)

    setting%Debugout%elem2R = .true.   ! select arrays to have debug output
    setting%Debugout%elemMR = .true.   ! select arrays to have debug output
    setting%Debugout%faceR  = .true.   ! note that not all are implemented

    setting%OutputThreadedLink%SuppressAllFiles = .true.

    setting%OutputThreadedLink%UseThisOutput = .true.
    setting%OutputThreadedLink%area = .true.
    setting%OutputThreadedLink%flowrate = .true.
    setting%OutputThreadedLink%velocity = .true.
    setting%OutputThreadedLink%eta = .true.
    setting%OutputThreadedLink%depth = .true.

    !!===========================================================

    !% bookkeeping routines
    call utility_get_datetime_stamp (setting%Time%DateTimeStamp)

    call debug_initialize (debugfile)

    call checking_consistency

    call initialize_arrayindex

    !% custom setup for hard-code test cases
    if (setting%TestCase%UseTestCase) then
        stop
        ! call test_case_initiation &
        !     (linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName, &
        !     bcdataDn, bcdataUp, &
        !     wdID, wdNumberPairs, wdManningsN, wdLength, wdZBottom, wdXDistance, &
        !     wdBreadth, wdWidthDepthData, wdCellType)
    else
        ! --------------------
        ! --- Initialize C API
        ! --------------------

        call initialize_api()

        ! Retrieve system properties from SWMM C
        call initialize_linknode_arrays &
            (linkI, nodeI, linkR, nodeR, linkYN, nodeYN, linkName, nodeName)

        ! Load inflows from SWMM C
        call inflow_load_inflows(nodeI, nodeR)

        ! Create system graph
        swmm_graph = get_network_graph()

        ! call finalize_api()

        ! Allocate boundary conditions
        nodeI(1:N_BCdnstream, ni_temp1) = pack(nodeI(:,ni_idx),nodeI(:,ni_node_type) == nBCdn)
        N_BCdnstream = count(nodeI(:,ni_node_type) == nBCdn)
        N_BCupstream = count(nodeI(:,ni_node_type) == nBCup)

        call bc_allocate(bcdataDn, bcdataUp)

        print *, "Setting up BC upstream"
        do ii = 1, N_BCupstream
            print *, "BC upstream", ii, '/', N_BCupstream
            jj = nodes_with_extinflow%array(ii)
            ts_ups = all_tseries(ext_inflows(ii)%t_series)
            bcdataUp(ii)%NodeID = jj
            bcdataUp(ii)%TimeArray = ts_ups%table%data(1)%array(1:ts_ups%table%tsize(1))
            bcdataUp(ii)%ValueArray = ts_ups%table%data(2)%array(1:ts_ups%table%tsize(2))
        enddo

        print *, "Setting up BC downstream"
        do ii = 1, N_BCdnstream
            print *, "BC dnstream", ii, '/', N_BCdnstream
            bcdataDn(ii)%NodeID = nodeI(ii, ni_temp1)
            allocate(bcdataDn(ii)%TimeArray(2))
            allocate(bcdataDn(ii)%ValueArray(2))
            bcdataDn(ii)%TimeArray = (/0.0, real((setting%time%endtime - setting%time%starttime)*dble(secsperday))/)
            bcdataDn(ii)%ValueArray = nr_Zbottom
        enddo
        ! --------------------
        ! --- Finalize C API
        ! --------------------
        ! print *, 'error - code only designed for use with test cases'
        ! stop
    end if

    !% create the network of elements from link and node data
    call network_initiation &
        (linkI, linkR, linkYN, linkName, &
        nodeI, nodeR, nodeYN, nodeName, &
        elem2R, elem2I, elem2YN, elem2Name, &
        elemMR, elemMI, elemMYN, elemMName, &
        faceR,  faceI,  faceYN,  faceName, swmm_graph)
    print *, 'in main'


    print *, "Start Time", setting%time%starttime
    print *, "End Time", setting%time%endtime
    print *, "Time step", setting%time%dt
    print *, "Number of elements", sum(linkI(:, li_N_element))

    !% check the boundary condition data arrays are correctly defined
    call bc_checks(bcdataUp, bcdataDn, elem2I, faceI, nodeI )

    !% set the initial conditions throughout
    call initial_condition_setup &
        (elem2R, elem2I, elem2YN, elemMR, elemMI, elemMYN, faceR, faceI, faceYN, &
        linkR, linkI, nodeR, nodeI, bcdataDn, bcdataUp, real(setting%Time%StartTime), &
        wdID, wdNumberPairs, wdManningsN, wdLength, wdZBottom, wdXDistance, &
        wdBreadth, wdWidthDepthData, wdCellType)

    !% check consistency of the smallvolume setup
    call checking_smallvolume_consistency (elem2R, elemMR)

    ! initialize the diagnostics
    call diagnostic_initialize &
        (diagnostic, elem2R, elem2I, elemMR, elemMI, faceR, &
        bcdataUp, bcdataDn)

    !% setting a zero starting condition is useful for robustness tests
    print *, 'in main setting flowrate and velocity to 0'
    elem2R(:,e2r_Velocity) = 0.0
    elem2R(:,e2r_Flowrate) = 0.0
    elemMR(:,eMr_Velocity) = 0.0
    elemMR(:,eMr_Flowrate) = 0.0
    elemMR(:,eMr_FlowrateUp(:)) = 0.0
    elemMR(:,eMr_FlowrateDn(:)) = 0.0
    elemMR(:,eMr_VelocityDn(:)) = 0.0
    elemMR(:,eMr_VelocityUp(:)) = 0.0
    faceR(1:size(faceR,1)-1,fr_Velocity_d) = 0.0
    faceR(1:size(faceR,1)-1,fr_Velocity_u) = 0.0
    faceR(1:size(faceR,1)-1,fr_Flowrate) = 0.0

    ! initialize output by threaded link
    call output_threaded_by_link_initialize (threadedfile)

    !%  time marching of continuity and momentum
    call time_marching &
        (elem2R, elemMR, faceR, elem2I, elemMI, faceI, elem2YN, elemMYN, faceYN, &
        bcdataDn, bcdataUp, linkI, debugfile, diagnostic, threadedfile, &
        wdID, wdNumberPairs, wdManningsN, wdLength, wdZBottom, wdXDistance, &
        wdBreadth, wdWidthDepthData, wdCellType)

    !% uncomment this if you want a final debug output
    ! call debug_output &
    !    (debugfile, &
    !     elem2R, elem2I, elem2YN, elemMR, elemMI, elemMYN, faceR, faceI, faceYN, &
    !     bcdataUp, bcdataDn)

    !
    !=========================================================
    ! FINAL CHECKING
    !
    !%  check that index arrays were not altered during execution
    call initialize_arrayindex_status

    !%  close out the debug files
    call debug_finalize(debugfile)

    print *
    print *, 'finished main program'
    print *, '====================='
    print *, char(7)  ! sound the system beep

end program main
!==========================================================================
