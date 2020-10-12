! module test_cases
!
! Calling routines for custom test cases, e.g. calls the case_simple_channel
! functions to setup a single channel reach.
!
!==========================================================================
!
module test_cases
    !
    use array_index
    use bc
    use case_simple_channel
    use case_simple_orifice
    use case_simple_weir
    use case_y_channel
    use case_y_storage_channel
    use case_waller_creek
    use case_simple_pipe
    use data_keys
    use globals
    use read_width_depth
    use setting_definition
    use utility
    use xsect_tables

    implicit none

    private

    public :: test_case_initiation

    integer :: debuglevel = 0

contains
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine test_case_initiation &
        (linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName,     &
        bcdataDn, bcdataUp, newID, newNumberPairs, newManningsN, newLength, &
        newZBottom, newXDistance, newBreadth, newWidthDepthData, newCellType)

        character(64) :: subroutine_name = 'test_case_initiation'

        integer,      dimension(:,:), allocatable, intent(out) :: linkI
        integer,      dimension(:,:), allocatable, intent(out) :: nodeI
        real,         dimension(:,:), allocatable, intent(out) :: linkR
        real,         dimension(:,:), allocatable, intent(out) :: nodeR
        logical,      dimension(:,:), allocatable, intent(out) :: linkYN
        logical,      dimension(:,:), allocatable, intent(out) :: nodeYN
        type(string), dimension(:),   allocatable, intent(out) :: linkName
        type(string), dimension(:),   allocatable, intent(out) :: nodeName
        type(bcType), dimension(:),   allocatable, intent(out) :: bcdataUp, bcdataDn

        real, dimension(:), allocatable :: depth_dnstream, depth_upstream, init_depth
        real, dimension(:), allocatable :: subdivide_length, channel_length
        real, dimension(:), allocatable :: channel_breadth, channel_topwidth
        real, dimension(:), allocatable :: lowerZ, upperZ, flowrate
        real, dimension(:), allocatable :: area, velocity,  Froude, ManningsN
        real, dimension(:), allocatable :: fullDepth, inletOffset, sideSlope
        real, dimension(:), allocatable :: parabolaValue, leftSlope, rightSlope
        real, dimension(:), allocatable :: dischargeCoefficient1, dischargeCoefficient2
        real, dimension(:), allocatable :: endContractions

        integer, dimension(:), allocatable :: idepth_type
        integer, dimension(:), allocatable :: channel_geometry

        real :: CFL

        integer :: first_step, last_step, display_interval, mm

        real :: climit, cvel, uz, lz

        !Waller Creek
        integer, dimension(:), allocatable :: init_ID
        integer, dimension(:), allocatable :: init_numberPairs
        real,    dimension(:), allocatable :: init_ManningsN
        real,    dimension(:), allocatable :: init_Length
        real,    dimension(:), allocatable :: init_zBottom
        real,    dimension(:), allocatable :: init_xDistance
        real,    dimension(:), allocatable :: init_Breadth
        real,    dimension(:,:,:),  allocatable :: init_widthDepthData
        type(string), dimension(:), allocatable :: init_cellType
        real, dimension(:),         allocatable :: faceZBottom

        integer, dimension(:),      allocatable :: newID
        integer, dimension(:),      allocatable :: newNumberPairs
        real,    dimension(:),      allocatable :: newManningsN
        real,    dimension(:),      allocatable :: newLength
        real,    dimension(:),      allocatable :: newZBottom
        real,    dimension(:),      allocatable :: newXDistance
        real,    dimension(:),      allocatable :: newBreadth
        real,    dimension(:,:,:),  allocatable, target :: newWidthDepthData
        type(string), dimension(:), allocatable :: newCellType(:)


        real :: inflowBC, heightBC, Waller_Creek_initial_depth
        real :: geometry_downstream_minimum_length
        real :: Waller_Creek_cellsize_target

        logical :: geometry_add_downstream_buffer

        integer :: unit = 11
        integer :: n_rows_in_file_node = 0
        integer :: max_number_of_pairs = 0

        real :: subdivide_length_check = 10.0

        integer :: ii

        real, pointer :: widthAtLayerTop(:,:), depthAtLayerTop(:,:), areaThisLayer(:,:)
        real, pointer :: areaTotalBelowThisLayer(:,:), dWidth(:,:)
        real, pointer :: dDepth(:,:), angle(:,:), perimeterBelowThisLayer(:,:)
        real, pointer :: area_difference(:,:), local_difference(:,:)

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        select case (setting%TestCase%TestName)

            !% Write a new case statement for each unique test case
          case ('simple_channel_001')

            N_link = 1
            N_node = 1
            N_BCupstream = 1
            N_BCdnstream = 1

            !% create the local variables that must be populated to set up the test case
            call control_variable_allocation &
                (init_depth, depth_dnstream, depth_upstream, lowerZ, upperZ, channel_length, &
                channel_breadth, channel_topwidth, subdivide_length, flowrate, area,        &
                velocity,  Froude, ManningsN, idepth_type, channel_geometry, parabolaValue, &
                leftSlope, rightSlope, sideslope, inletOffset, dischargeCoefficient1,       &
                dischargeCoefficient2, fullDepth, endContractions)

            ! step controls
            display_interval = 1000
            first_step = 1
            last_step  =  10000 ! note 1000 is good enough to show blow up or not, 10000 is smooth

            ! set up flow and time step for differen subcases
            ! tests that ran:  Fr = 0.25, 0.5
            Froude       = 0.25   ! determines flowrate and slope to get Froude
            CFL          = 0.25  ! determines dt from subdivide_length

            ! keep these physics fixed
            channel_breadth = 3.0

            ! assign geometry type for links
            do ii=1,N_link
                channel_geometry(ii) = lRectangular
            end do

            depth_upstream  = 1.0
            depth_dnstream  = 1.0
            init_depth      = 1.0
            idepth_type     = 1  !1 = uniform, 2=linear, 3=exponential decay
            ManningsN       = 0.03
            channel_length    = 10000.0
            lowerZ          = 1.0
            subdivide_length = 5000.0

            ! rectangular geometry
            parabolaValue = zeroR
            leftSlope     = zeroR
            rightSlope    = zeroR

            call froude_driven_setup &
                (upperZ(1), area(1), flowrate(1), velocity(1),  &
                Froude(1),  channel_breadth(1), channel_topwidth(1), ManningsN(1), &
                channel_length(1), lowerZ(1),  init_depth(1), &
                channel_geometry(1), parabolaValue(1), leftSlope(1), rightSlope(1))

            call this_setting_for_time_and_steps &
                (CFL, velocity, init_depth, subdivide_length, &
                first_step, last_step, display_interval,2)

            call case_simple_channel_initialize &
                (channel_length(1), channel_breadth(1), channel_topwidth(1), subdivide_length(1), &
                lowerZ(1), upperZ(1), flowrate(1), init_depth(1), depth_upstream(1), depth_dnstream(1), &
                ManningsN(1), lManningsN, idepth_type(1),                                   &
                linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName,    &
                bcdataDn, bcdataUp)

            if (.not. setting%Debugout%SuppressAllFiles) then
                call write_testcase_setup_file &
                    (Froude, CFL, flowrate, velocity, init_depth, depth_upstream,   &
                    depth_dnstream, channel_breadth, channel_topwidth, area, &
                    channel_length, subdivide_length, &
                    lowerZ, upperZ, ManningsN)
            endif

            !print *, flowrate, depth_dnstream
            !stop

            !% Write a new case statement for each unique test case
          case ('waller_creek')

            !open the Waller Creek depth list
            open(newunit=unit, file='WLR_WidthDepthList.txt', status='OLD')

            ! get the number of links and number of pairs per each link from the file
            n_rows_in_file_node = read_number_of_cells(unit)
            max_number_of_pairs = read_max_number_of_pairs(unit)

            ! read the entire data from the width-depth list
            call read_widthdepth_pairs &
                (unit, init_ID, init_numberPairs, init_ManningsN, init_Length,    &
                init_zBottom, init_xDistance, init_Breadth, init_widthDepthData, &
                init_cellType)

            ! subdivide length for checking the length of nonmonotonic elements
            subdivide_length_check = 10.0

            ! dividing nun monotonic elements
            call nonmonotonic_subdivide &
                (init_ID, init_numberPairs, init_ManningsN, init_Length,          &
                init_zBottom, init_xDistance, init_Breadth, init_widthDepthData, &
                init_cellType, n_rows_in_file_node, faceZBottom,                 &
                max_number_of_pairs, newID, newNumberPairs, newManningsN,        &
                newLength, newZBottom, newXDistance, newBreadth,                 &
                newWidthDepthData, newCellType, subdivide_length_check)

            N_link = newID(size(newID))
            N_node = N_link + 1
            N_BCupstream = 1
            N_BCdnstream = 1

            !% create the local variables that must be populated to set up the test case
            call control_variable_allocation &
                (init_depth, depth_dnstream, depth_upstream, lowerZ, upperZ, channel_length, &
                channel_breadth, channel_topwidth, subdivide_length, flowrate, area,        &
                velocity,  Froude, ManningsN, idepth_type, channel_geometry, parabolaValue, &
                leftSlope, rightSlope, sideslope, inletOffset, dischargeCoefficient1,       &
                dischargeCoefficient2, fullDepth, endContractions)

            ! step controls
            display_interval = 1000
            first_step = 1
            last_step  =  1000 ! note 1000 is good enough to show blow up or not, 10000 is smooth

            ! set up flow and time step for differen subcases
            Froude(:)    = 0.25   ! determines flowrate and slope to get Froude
            CFL          = 0.6  ! determines dt from subdivide_length

            ! keep these physics fixed
            channel_breadth     = newBreadth

            ! assign geometry type for links
            do ii=1,N_link
                channel_geometry(ii) = lWidthDepth
            end do

            depth_upstream(:)   = 0.3
            depth_dnstream(:)   = 0.3
            init_depth(:)       = 0.3
            idepth_type         = 1  !1 = uniform, 2=linear, 3=exponential decay
            ManningsN           = newManningsN
            channel_length      = newLength
            lowerZ              = newZBottom
            subdivide_length(:) = 5000.0

            ! rectangular geometry
            parabolaValue = zeroR
            leftSlope     = zeroR
            rightSlope    = zeroR

            !calculate the geometry related information from widthDepth information
            !and store it at the same matrix
            call widthdepth_pair_auxiliary (newWidthDepthData, newCellType, newNumberPairs)

            call Initial_condition_for_width_depth_system &
                (upperZ, area, flowrate, velocity, Froude,  channel_breadth, &
                channel_topwidth, ManningsN, channel_length, lowerZ, &
                init_depth, newWidthDepthData)

            call this_setting_for_time_and_steps &
                (CFL, velocity, init_depth, subdivide_length, first_step, last_step, &
                display_interval, 2)

            call case_waller_creek_initialize &
                (channel_length, channel_breadth, channel_topwidth, subdivide_length, &
                lowerZ, upperZ, flowrate, init_depth, depth_upstream, depth_dnstream, &
                ManningsN, lManningsN, idepth_type,                                   &
                linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName,    &
                bcdataDn, bcdataUp, newID, newNumberPairs, &
                newXDistance, newWidthDepthData, newCellType)

            if (.not. setting%Debugout%SuppressAllFiles) then
                call write_testcase_setup_file &
                    (Froude, CFL, flowrate, velocity, init_depth, depth_upstream,   &
                    depth_dnstream, channel_breadth, channel_topwidth, area, &
                    channel_length, subdivide_length, &
                    lowerZ, upperZ, ManningsN)
            endif

          case ('y_channel_002')

            N_link = 3
            N_node = 4
            N_BCupstream = 2
            N_BCdnstream = 1

            call control_variable_allocation &
                (init_depth, depth_dnstream, depth_upstream, lowerZ, upperZ, channel_length, &
                channel_breadth, channel_topwidth, subdivide_length, flowrate, area,        &
                velocity,  Froude, ManningsN, idepth_type, channel_geometry, parabolaValue, &
                leftSlope, rightSlope, sideslope, inletOffset, dischargeCoefficient1,       &
                dischargeCoefficient2, fullDepth, endContractions)

            ! step controls
            display_interval = 1000
            first_step = 1
            last_step  = 10000! note 1000 is good enough to show blow up or not, 10000 is smooth

            ! set up flow and time step for differen subcases
            ! tests that ran:  Fr = 0.25, 0.5
            Froude(1)       = 0.8   ! determines flowrate and slope to get Froude

            Froude(2)       = 0.8   ! determines flowrate and slope to get Froude

            Froude(3)       = 0.8   ! determines flowrate and slope to get Froude

            CFL          = 0.6  ! determines dt from subdivide_length

            init_depth(1:3)    = 1.0
            depth_dnstream(1)  = 1.0
            depth_upstream(1)  = 1.0 ! junction

            depth_dnstream(2:3) = depth_upstream(1) ! junction should be consistent

            depth_upstream(2)  = 1.0 ! upstream bc right
            depth_upstream(3)  = 1.0 ! upstream bc left
            idepth_type     = 1  !1 = uniform, 2=linear, 3=exponential decay
            ManningsN       = 0.03

            channel_breadth(1)   = 3.0
            channel_breadth(2)   = 3.0
            channel_breadth(3)   = 3.0

            ! assign geometry type for links
            do ii=1,N_link
                channel_geometry(ii) = lRectangular
            end do

            channel_length(1)    = 1000.0
            channel_length(2)    = 1000.0
            channel_length(3)    = 1000.0

            lowerZ(1)           = 1.0
            subdivide_length(1) = 100.0
            subdivide_length(2) = 100.0
            subdivide_length(3) = 100.0

            ! rectangular geometry
            parabolaValue = zeroR
            leftSlope     = zeroR
            rightSlope    = zeroR


            ! get consistent bottom Z values for the desired Froude number in each link
            do mm=1,N_link
                if (mm==1) then
                    ! start with the Z for the inflow link
                    lz = lowerZ(1)
                end if
                call froude_driven_setup &
                    (uz, area(mm), flowrate(mm), velocity(mm), Froude(mm),         &
                    channel_breadth(mm), channel_topwidth(mm), ManningsN(mm),     &
                    channel_length(mm), lz, init_depth(mm), channel_geometry(mm), &
                    parabolaValue(mm), leftSlope(mm), rightSlope(mm) )
                select case (mm)
                  case (1)
                    ! the upstream z of the downstream link becomes the lower z of the upstream links
                    lz = uz
                    upperZ(1) = uz
                  case (2,3)
                    lowerZ(mm) = upperZ(1)
                    upperZ(mm) = uz
                    lz = upperZ(1)
                end select
            end do

            call this_setting_for_time_and_steps &
                (CFL, velocity, init_depth, subdivide_length, first_step, last_step, &
                display_interval, 2)

            call case_y_channel_initialize &
                (channel_length, channel_breadth, channel_topwidth, subdivide_length, &
                lowerZ, upperZ, flowrate, init_depth, depth_upstream, depth_dnstream,       &
                ManningsN, lManningsN, idepth_type,                             &
                linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName, &
                bcdataDn, bcdataUp)

            if (.not. setting%Debugout%SuppressAllFiles) then
                call write_testcase_setup_file &
                    (Froude, CFL, flowrate, velocity, init_depth, depth_upstream,   &
                    depth_dnstream, channel_breadth, channel_topwidth, area, &
                    channel_length, subdivide_length, &
                    lowerZ, upperZ, ManningsN)
            endif


            !print *, flowrate
            !print *, linkR(:,lr_InitialFlowrate)
            !print *, trim(subroutine_name)
            !stop
          case ('simple_weir_003')

            N_link = 3
            N_node = 4
            N_BCupstream = 1
            N_BCdnstream = 1

            !% create the local variables that must be populated to set up the test case
            call control_variable_allocation &
                (init_depth, depth_dnstream, depth_upstream, lowerZ, upperZ, channel_length, &
                channel_breadth, channel_topwidth, subdivide_length, flowrate, area,        &
                velocity,  Froude, ManningsN, idepth_type, channel_geometry, parabolaValue, &
                leftSlope, rightSlope, sideslope, inletOffset, dischargeCoefficient1,       &
                dischargeCoefficient2, fullDepth, endContractions)

            ! step controls
            display_interval = 100
            first_step = 1
            last_step  =  30000

            ! set up flow and time step for differen subcases
            ! tests that ran:  Fr = 0.25, 0.5

            ! This is from case_y_channel
            Froude       = 0.25   ! determines flowrate and slope to get Froude
            CFL          = 0.25  ! determines dt from subdivide_length

            ! keep these physics fixed
            idepth_type        = 1  !1 = uniform, 2=linear, 3=exponential decay
            ManningsN          = 0.03

            channel_geometry(1) = lRectangular
            lowerZ(1)          = 1.0
            depth_dnstream(1)  = 1.0e-2
            depth_upstream(1)  = 1.0e-2
            init_depth(1)      = 1.0e-2

            channel_geometry(2) = lTriangular
            depth_dnstream(2)  = 1.0e-2         !This is the depth in weir
            depth_upstream(2)  = 1.0e-2         !This is the depth in weir
            init_depth(2)      = 1.0e-2

            channel_geometry(3) = lRectangular
            depth_dnstream(3)  = 0.25
            depth_upstream(3)  = 0.25
            init_depth(3)      = 0.25

            channel_breadth      = 3.0
            channel_length       = 1000.0
            channel_length(2)    = 1        !This is Weir Length

            subdivide_length     = 500.0
            subdivide_length(2)  = 1        !We are not subdividing weir element. So this value is same as weir length

            leftSlope             = 1.0
            rightSlope            = 1.0
            sideSlope             = nullValueR
            inletOffset           = nullValueR
            dischargeCoefficient1 = nullValueR
            dischargeCoefficient2 = nullValueR
            fullDepth             = nullValueR
            endContractions       = nullValueR

            ! wier settings
            inletOffset(2)             = 1.0
            dischargeCoefficient1(2)   = 1.40
            fullDepth(2)               = 1.5
            sideSlope(2)               = 1.0

            ! get consistent bottom Z values for the desired Froude number in each link
            do mm=1,N_link
                if (mm==1) then
                    ! start with the Z for the inflow link
                    lz = lowerZ(1)
                end if
                select case (mm)
                  case (1)
                    call froude_driven_setup &
                        (uz, area(mm), flowrate(mm), velocity(mm),                               &
                        Froude(mm), channel_breadth(mm), channel_topwidth(mm), &
                        ManningsN(mm), channel_length(mm), &
                        lz, init_depth(mm), channel_geometry(mm), &
                        parabolaValue(mm), leftSlope(mm), rightSlope(mm) )
                    ! the upstream z of the downstream link becomes the lower z of the upstream links
                    lz = uz
                    upperZ(1) = uz
                  case (2)
                    call weir_setup &
                        (uz, area(mm), flowrate(mm), velocity(mm), sideslope(mm), Froude(mm), &
                        channel_breadth(mm), ManningsN(mm), channel_length(mm), lz,           &
                        depth_upstream(mm) )
                    lowerZ(mm) = upperZ(1)
                    upperZ(mm) = uz
                    lz = uz
                  case (3)
                    call froude_driven_setup &
                        (uz, area(mm), flowrate(mm), velocity(mm),                               &
                        Froude(mm), channel_breadth(mm), channel_topwidth(mm), &
                        ManningsN(mm), channel_length(mm), &
                        lz, init_depth(mm), channel_geometry(mm), &
                        parabolaValue(mm), leftSlope(mm), rightSlope(mm) )
                    lowerZ(mm) = upperZ(2)
                    upperZ(mm) = uz
                end select
            end do

            call this_setting_for_time_and_steps &
                (CFL, velocity, depth_upstream, subdivide_length, &
                first_step, last_step, display_interval,2)

            call case_simple_weir_initialize &
                (channel_length, channel_breadth, subdivide_length, lowerZ, upperZ,    &
                flowrate, depth_upstream, depth_dnstream, init_depth,                 &
                sideslope, inletOffset, dischargeCoefficient1, dischargeCoefficient2, &
                fullDepth, endContractions, ManningsN, lManningsN, idepth_type,       &
                linkR, nodeR, linkI, nodeI,linkYN, nodeYN, linkName, nodeName,        &
                bcdataDn, bcdataUp)

            if (.not. setting%Debugout%SuppressAllFiles) then
                call write_testcase_setup_file &
                    (Froude, CFL, flowrate, velocity, init_depth, depth_upstream,   &
                    depth_dnstream, channel_breadth, channel_topwidth, area, &
                    channel_length, subdivide_length, &
                    lowerZ, upperZ, ManningsN)
            endif

          case ('simple_orifice_004')

            N_link = 3
            N_node = 4
            N_BCupstream = 1
            N_BCdnstream = 1

            !% create the local variables that must be populated to set up the test case
            call control_variable_allocation &
                (init_depth, depth_dnstream, depth_upstream, lowerZ, upperZ, channel_length, &
                channel_breadth, channel_topwidth, subdivide_length, flowrate, area,        &
                velocity,  Froude, ManningsN, idepth_type, channel_geometry, parabolaValue, &
                leftSlope, rightSlope, sideslope, inletOffset, dischargeCoefficient1,       &
                dischargeCoefficient2, fullDepth, endContractions)

            ! step controls
            display_interval = 100
            first_step = 1
            ! last_step  =  6500
            last_step  =  30000

            ! set up flow and time step for differen subcases
            ! tests that ran:  Fr = 0.25, 0.5

            ! This is from case_y_channel
            Froude       = 0.25   ! determines flowrate and slope to get Froude
            CFL          = 0.25  ! determines dt from subdivide_length

            ! keep these physics fixed
            idepth_type        = 1  !1 = uniform, 2=linear, 3=exponential decay
            ManningsN          = 0.03

            channel_geometry(1) = lRectangular
            lowerZ(1)          = 1.0
            depth_dnstream(1)  = 1.0e-2
            depth_upstream(1)  = 1.0e-2
            init_depth(1)      = 1.0e-2

            channel_geometry(2) = lCircular
            depth_dnstream(2)  = 1.0e-2         !This is the depth in weir
            depth_upstream(2)  = 1.0e-2         !This is the depth in weir
            init_depth(2)      = 1.0e-2

            channel_geometry(3) = lRectangular
            depth_dnstream(3)  = 0.25
            depth_upstream(3)  = 0.25
            init_depth(3)      = 0.25

            channel_breadth      = 3.0
            channel_breadth(2)   = 1.5
            channel_length       = 1000.0
            channel_length(2)    = 1        !This is Weir Length

            subdivide_length     = 500.0
            subdivide_length(2)  = 1        !We are not subdividing weir element. So this value is same as weir length

            leftSlope             = nullValueR
            rightSlope            = nullValueR
            sideSlope             = nullValueR
            inletOffset           = nullValueR
            dischargeCoefficient1 = nullValueR
            dischargeCoefficient2 = nullValueR
            fullDepth             = nullValueR
            endContractions       = nullValueR

            ! wier settings
            inletOffset(2)             = 1.0
            dischargeCoefficient1(2)   = 0.6
            fullDepth(2)               = 1.5

            ! get consistent bottom Z values for the desired Froude number in each link
            do mm=1,N_link
                if (mm==1) then
                    ! start with the Z for the inflow link
                    lz = lowerZ(1)
                end if
                select case (mm)
                  case (1)
                    call froude_driven_setup &
                        (uz, area(mm), flowrate(mm), velocity(mm),                               &
                        Froude(mm), channel_breadth(mm), channel_topwidth(mm), &
                        ManningsN(mm), channel_length(mm), &
                        lz, init_depth(mm), channel_geometry(mm), &
                        parabolaValue(mm), leftSlope(mm), rightSlope(mm) )
                    ! the upstream z of the downstream link becomes the lower z of the upstream links
                    lz = uz
                    upperZ(1) = uz
                  case (2)
                    call orifice_setup &
                        (uz, area(mm), flowrate(mm), velocity(mm), sideslope(mm), Froude(mm), &
                        channel_breadth(mm), ManningsN(mm), channel_length(mm), lz,           &
                        depth_upstream(mm) )
                    lowerZ(mm) = upperZ(1)
                    upperZ(mm) = uz
                    lz = uz
                  case (3)
                    call froude_driven_setup &
                        (uz, area(mm), flowrate(mm), velocity(mm),                               &
                        Froude(mm), channel_breadth(mm), channel_topwidth(mm), &
                        ManningsN(mm), channel_length(mm), &
                        lz, init_depth(mm), channel_geometry(mm), &
                        parabolaValue(mm), leftSlope(mm), rightSlope(mm) )
                    lowerZ(mm) = upperZ(2)
                    upperZ(mm) = uz
                end select
            end do

            call this_setting_for_time_and_steps &
                (CFL, velocity, depth_upstream, subdivide_length, &
                first_step, last_step, display_interval,2)

            call case_simple_orifice_initialize &
                (channel_length, channel_breadth, subdivide_length, lowerZ, upperZ,    &
                flowrate, depth_upstream, depth_dnstream, init_depth,                 &
                sideslope, inletOffset, dischargeCoefficient1, dischargeCoefficient2, &
                fullDepth, endContractions, ManningsN, lManningsN, idepth_type,       &
                linkR, nodeR, linkI, nodeI,linkYN, nodeYN, linkName, nodeName,        &
                bcdataDn, bcdataUp)

            if (.not. setting%Debugout%SuppressAllFiles) then
                call write_testcase_setup_file &
                    (Froude, CFL, flowrate, velocity, init_depth, depth_upstream,   &
                    depth_dnstream, channel_breadth, channel_topwidth, area, &
                    channel_length, subdivide_length, &
                    lowerZ, upperZ, ManningsN)
            endif

          case ('y_storage_channel_005')

            N_link = 3
            N_node = 4
            N_BCupstream = 2
            N_BCdnstream = 1

            call control_variable_allocation &
                (init_depth, depth_dnstream, depth_upstream, lowerZ, upperZ, channel_length, &
                channel_breadth, channel_topwidth, subdivide_length, flowrate, area,        &
                velocity,  Froude, ManningsN, idepth_type, channel_geometry, parabolaValue, &
                leftSlope, rightSlope, sideslope, inletOffset, dischargeCoefficient1,       &
                dischargeCoefficient2, fullDepth, endContractions)

            ! step controls
            display_interval = 100
            first_step = 1
            last_step  = 1000! note 1000 is good enough to show blow up or not, 10000 is smooth

            ! set up flow and time step for differen subcases
            ! tests that ran:  Fr = 0.25, 0.5
            Froude(1)       = 0.8   ! determines flowrate and slope to get Froude
            Froude(2)       = 0.8  ! determines flowrate and slope to get Froude
            Froude(3)       = 0.8   ! determines flowrate and slope to get Froude

            CFL          = 0.6  ! determines dt from subdivide_length

            init_depth(1:3)    = 1.0
            depth_dnstream(1)  = 1.0
            depth_upstream(1)  = 1.0 ! junction

            depth_dnstream(2:3) = depth_upstream(1) ! junction should be consistent

            depth_upstream(2)  = 1.0 ! upstream bc right
            depth_upstream(3)  = 1.0 ! upstream bc left
            idepth_type     = 1  !1 = uniform, 2=linear, 3=exponential decay
            ManningsN       = 0.03

            channel_breadth(1)   = 3.0
            channel_breadth(2)   = 3.0
            channel_breadth(3)   = 3.0

            ! assign geometry type for links
            do ii=1,N_link
                channel_geometry(ii) = lRectangular
            end do

            channel_length(1)    = 1000.0
            channel_length(2)    = 1000.0
            channel_length(3)    = 1000.0

            lowerZ(1)           = 1.0
            subdivide_length(1) = 100.0
            subdivide_length(2) = 100.0
            subdivide_length(3) = 100.0

            ! rectangular geometry
            parabolaValue = zeroR
            leftSlope     = zeroR
            rightSlope    = zeroR


            ! get consistent bottom Z values for the desired Froude number in each link
            do mm=1,N_link
                if (mm==1) then
                    ! start with the Z for the inflow link
                    lz = lowerZ(1)
                end if
                call froude_driven_setup &
                    (uz, area(mm), flowrate(mm), velocity(mm), Froude(mm),         &
                    channel_breadth(mm), channel_topwidth(mm), ManningsN(mm),     &
                    channel_length(mm), lz, init_depth(mm), channel_geometry(mm), &
                    parabolaValue(mm), leftSlope(mm), rightSlope(mm) )
                select case (mm)
                  case (1)
                    ! the upstream z of the downstream link becomes the lower z of the upstream links
                    lz = uz
                    upperZ(1) = uz
                  case (2,3)
                    lowerZ(mm) = upperZ(1)
                    upperZ(mm) = uz
                    lz = upperZ(1)
                end select
            end do

            call this_setting_for_time_and_steps &
                (CFL, velocity, init_depth, subdivide_length, first_step, last_step, &
                display_interval, 2)

            call case_y_storage_channel_initialize &
                (channel_length, channel_breadth, channel_topwidth, subdivide_length, &
                lowerZ, upperZ, flowrate, init_depth, depth_upstream, depth_dnstream,       &
                ManningsN, lManningsN, idepth_type,                             &
                linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName, &
                bcdataDn, bcdataUp)

            if (.not. setting%Debugout%SuppressAllFiles) then
                call write_testcase_setup_file &
                    (Froude, CFL, flowrate, velocity, init_depth, depth_upstream,   &
                    depth_dnstream, channel_breadth, channel_topwidth, area, &
                    channel_length, subdivide_length, &
                    lowerZ, upperZ, ManningsN)
            endif

        !% Write a new case statement for each unique test case
          case ('simple_pipe_006')

            N_link = 1
            N_node = 1
            N_BCupstream = 1
            N_BCdnstream = 1

            !% create the local variables that must be populated to set up the test case
            call control_variable_allocation &
                (init_depth, depth_dnstream, depth_upstream, lowerZ, upperZ, channel_length, &
                channel_breadth, channel_topwidth, subdivide_length, flowrate, area,        &
                velocity,  Froude, ManningsN, idepth_type, channel_geometry, parabolaValue, &
                leftSlope, rightSlope, sideslope, inletOffset, dischargeCoefficient1,       &
                dischargeCoefficient2, fullDepth, endContractions)

            ! step controls
            display_interval = 10
            first_step = 1
            last_step  = 100 ! note 1000 is good enough to show blow up or not, 10000 is smooth

            ! set up flow and time step for differen subcases
            ! tests that ran:  Fr = 0.25, 0.5
            Froude       = 0.25  ! determines flowrate and slope to get Froude
            CFL          = 0.25  ! determines dt from subdivide_length

            ! keep these physics fixed
            channel_breadth = 2.0

            ! assign geometry type for links
            do ii=1,N_link
                channel_geometry(ii) = lCircular
            end do

            depth_upstream  = 1.75
            depth_dnstream  = 1.75
            init_depth      = 1.75
            idepth_type     = 1  !1 = uniform, 2=linear, 3=exponential decay
            ManningsN       = 0.016
            channel_length  = 10000.0
            lowerZ          = 1.0
            subdivide_length = 5000.0
            fullDepth       = 2.0    

            ! rectangular geometry
            parabolaValue = zeroR
            leftSlope     = zeroR
            rightSlope    = zeroR

            call froude_driven_pipe_setup &
                (upperZ(1), area(1), fullDepth(1), flowrate(1), velocity(1),  &
                Froude(1),  channel_breadth(1), channel_topwidth(1), ManningsN(1), &
                channel_length(1), lowerZ(1),  init_depth(1))

            call this_setting_for_time_and_steps &
                (CFL, velocity, init_depth, subdivide_length, &
                first_step, last_step, display_interval,2)

            call case_simple_pipe_initialize &
                (channel_length(1), channel_breadth(1), channel_topwidth(1), subdivide_length(1), &
                lowerZ(1), upperZ(1), flowrate(1), init_depth(1), depth_upstream(1), depth_dnstream(1), &
                ManningsN(1), lManningsN, idepth_type(1), fullDepth(1),                                   &
                linkR, nodeR, linkI, nodeI, linkYN, nodeYN, linkName, nodeName,    &
                bcdataDn, bcdataUp)

            if (.not. setting%Debugout%SuppressAllFiles) then
                call write_testcase_setup_file &
                    (Froude, CFL, flowrate, velocity, init_depth, depth_upstream,   &
                    depth_dnstream, channel_breadth, channel_topwidth, area, &
                    channel_length, subdivide_length, &
                    lowerZ, upperZ, ManningsN)
            endif

            !print *, flowrate, depth_dnstream
            !stop
            
          case default
            print *, setting%TestCase%TestName
            print *, 'error: no valid test case of ',&
                trim(setting%TestCase%TestName),' in ',subroutine_name
            stop
        end select

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine test_case_initiation
    !
    !==========================================================================
    !
    ! PRIVATE BELOW HERE
    !
    !==========================================================================
    !
    subroutine control_variable_allocation &
        (init_depth, depth_dnstream, depth_upstream, lowerZ, upperZ, channel_length, &
        channel_breadth, channel_topwidth, subdivide_length, flowrate, area,        &
        velocity,  Froude, ManningsN, idepth_type, channel_geometry, parabolaValue, &
        leftSlope, rightSlope, sideSlope, inletOffset, dischargeCoefficient1,       &
        dischargeCoefficient2, fullDepth, endContractions)

        character(64) :: subroutine_name = 'control_variable_allocation'

        real, dimension(:), allocatable, intent(out) :: depth_dnstream, depth_upstream, init_depth
        real, dimension(:), allocatable, intent(out) :: subdivide_length, channel_length
        real, dimension(:), allocatable, intent(out) :: channel_breadth, channel_topwidth
        real, dimension(:), allocatable, intent(out) :: lowerZ, upperZ, flowrate
        real, dimension(:), allocatable, intent(out) :: area, velocity, Froude, ManningsN
        real, dimension(:), allocatable, intent(out) :: parabolaValue, leftSlope, rightSlope
        real, dimension(:), allocatable, intent(out) :: inletOffset, sideSlope
        real, dimension(:), allocatable, intent(out) :: dischargeCoefficient1, dischargeCoefficient2
        real, dimension(:), allocatable, intent(out) :: fullDepth, endContractions

        integer, dimension(:), allocatable, intent(out) :: idepth_type
        integer, dimension(:), allocatable, intent(out) :: channel_geometry

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        allocate(init_depth(N_link))
        allocate(depth_dnstream(N_link))
        allocate(depth_upstream(N_link))
        allocate(lowerZ(N_link))
        allocate(upperZ(N_link))
        allocate(channel_length(N_link))
        allocate(channel_breadth(N_link))
        allocate(channel_topwidth(N_link))
        allocate(subdivide_length(N_link))
        !allocate(initial_flowrate(N_link))
        allocate(area(N_link))
        allocate(velocity(N_link))
        allocate(flowrate(N_link))
        allocate(Froude(N_link))
        allocate(ManningsN(N_link))
        allocate(channel_geometry(N_link))
        allocate(parabolaValue(N_link))
        allocate(leftSlope(N_link))
        allocate(rightSlope(N_link))
        allocate(sideSlope(N_link))
        allocate(inletOffset(N_link))
        allocate(dischargeCoefficient1(N_link))
        allocate(dischargeCoefficient2(N_link))
        allocate(fullDepth(N_link))
        allocate(idepth_type(N_link))
        allocate(endContractions(N_link))

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine control_variable_allocation
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine this_setting_for_time_and_steps &
        (CFL, velocity, depth, subdivide_length, first_step, last_step, &
        display_interval, dt_significant_digits)

        character(64) :: subroutine_name = 'this_setting_for_time_and_steps'

        real,  intent(in) :: CFL, velocity(:), depth(:), subdivide_length(:)

        integer, intent(in) :: first_step, last_step, display_interval, dt_significant_digits

        real,  dimension(size(velocity)) :: dtSet, CFLset

        real       :: dtmin
        integer    :: dtscale

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        ! use the same CFL in every link
        CFLset = CFL
        ! get the set of time step (dt) base on every branch
        dtSet = get_dt_from_CFL (CFL, velocity, depth, subdivide_length)
        ! get the minimum dt value
        dtmin  = minval(dtSet)
        ! get the largest n for 10^n relative to the dtmin
        dtscale = utility_scale_of_number(dtmin)


        setting%Time%dt = utility_round_to_significant_digits(dtmin,dt_significant_digits)

        setting%Step%Current = 1

        setting%Step%First = first_step
        setting%Step%Final = last_step

        setting%Debugout%DisplayInterval = display_interval
        setting%OutputThreadedLink%DisplayInterval = display_interval

        setting%Time%StartTime = 0.0
        setting%Time%EndTime = setting%Time%StartTime  &
            + setting%Time%dt * (setting%Step%Final - setting%Step%First + 1)

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine this_setting_for_time_and_steps
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine froude_driven_setup &
        (upperZ, area, flowrate, velocity,  &
        Froude,  breadth, topwidth, ManningsN, total_length, &
        lowerZ, depth, channel_geometry, parabolaValue, leftSlope, rightSlope)

        character(64) :: subroutine_name = 'froude_driven_setup'

        real,  intent(out)    :: area, flowrate, velocity, upperZ, topwidth
        real,  intent(in)     :: Froude,  breadth, ManningsN, lowerZ, total_length
        real,  intent(in)     :: depth
        real,  intent(in)     :: parabolaValue, leftSlope, rightSlope

        real :: perimeter, rh, slope
        integer, intent(in) :: channel_geometry

        real :: xParabola , hDepth



        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        select case (channel_geometry)
          case(lRectangular)
            hDepth = depth
            topwidth = breadth
            area = hDepth * breadth
            perimeter = 2.0 * hDepth + breadth

          case(lParabolic)
            hDepth = (twoR / threeR) * depth
            topwidth = twoR * sqrt(abs(depth/parabolaValue))
            area = twothirdR * topwidth * depth
            xParabola = fourR * depth / topwidth
            perimeter = onehalfR * topwidth &
                *( sqrt( oneR + xParabola**twoR ) + oneR / xParabola &
                * log ( xParabola + sqrt( oneR + xParabola**twoR )) )

          case(lTrapezoidal)
            area = (breadth + onehalfR * (leftSlope + rightSlope) * depth ) * depth
            topwidth = breadth + (leftSlope + rightSlope) * depth
            hDepth = area / topwidth
            perimeter = breadth + depth * (sqrt(oneR + leftSlope**twoR ) &
                + sqrt(oneR + rightSlope**twoR))

          case(lTriangular)
            area = onehalfR * (leftSlope + rightSlope) * depth ** twoR
            topwidth = (leftSlope + rightSlope) * depth
            hDepth = onehalfR * depth
            perimeter = depth * (sqrt(oneR + leftSlope**twoR) + sqrt(oneR + rightSlope**twoR))

        end select

        rh = area / perimeter
        velocity = Froude * sqrt(grav * hDepth)
        flowrate = area * velocity
        slope = (velocity * ManningsN / (rh**(2.0/3.0)) )**2
        upperZ = lowerZ + slope * total_length


        ! print *,'-----------------'
        ! print *, area, 'area'
        ! print *, hDepth, 'hDepth'
        ! print *, perimeter, 'perimeter'
        ! print *, rh, 'rh'
        ! print *, velocity, 'velocity'
        ! print *, flowrate, 'flowrate'
        ! print *, slope, 'slope'
        ! print *, upperZ, 'upperZ', lowerZ, 'lowerZ'
        ! print *, total_length, 'total_length'
        ! print *, slope*total_length, 'slope*total_length'
        ! print *,'-----------------'

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine froude_driven_setup
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine weir_setup &
        (upperZ, area, flowrate, velocity, sideslope, Froude,  breadth, &
        ManningsN, total_length, lowerZ, depth)

        character(64) :: subroutine_name = 'weir_setup'

        real,  intent(out)    :: area, flowrate, velocity, upperZ
        real,  intent(in)     :: Froude,  breadth, ManningsN, lowerZ, total_length
        real,  intent(in)     :: depth, sideslope

        real :: perimeter, rh, slope

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name


        ! These needed to be changed when the weir is surcharged
        area        = sideslope* depth ** twoR
        perimeter   = twoR * depth * sqrt(1 + sideslope ** 2)
        !rh          = area / perimeter
        !velocity    = (setting%Weir%WeirDischargeCoeff * setting%Weir%WeirSideSlope * depth ** 2.5) / area
        !flowrate    = area * velocity
        !Hard coading this to run the weir without any water
        rh          = 0.0
        velocity    = 0.0
        flowrate    = 0.0
        slope       = zeroR
        upperZ      = lowerZ

        ! print *,'-----------------'
        ! print *, area, 'area'
        ! print *, perimeter, 'perimeter'
        ! print *, rh, 'rh'
        ! print *, velocity, 'velocity'
        ! print *, flowrate, 'flowrate'
        ! print *, slope, 'slope'
        ! print *, upperZ, 'upperZ', lowerZ, 'lowerZ'
        ! print *, total_length, 'total_length'
        ! print *, slope*total_length, 'slope*total_length'
        ! print *,'-----------------'

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine weir_setup
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine orifice_setup &
        (upperZ, area, flowrate, velocity, sideslope, Froude,  breadth, &
        ManningsN, total_length, lowerZ, depth)

        character(64) :: subroutine_name = 'orifice_setup'

        real,  intent(out)    :: area, flowrate, velocity, upperZ
        real,  intent(in)     :: Froude,  breadth, ManningsN, lowerZ, total_length
        real,  intent(in)     :: depth, sideslope

        real :: perimeter, rh, slope

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name


        ! These needed to be changed when the weir is surcharged
        area        = zeroR
        perimeter   = zeroR
        !for now testing everything with zero
        rh          = zeroR
        velocity    = zeroR
        flowrate    = zeroR
        slope       = zeroR
        upperZ      = lowerZ

        ! print *,'-----------------'
        ! print *, area, 'area'
        ! print *, perimeter, 'perimeter'
        ! print *, rh, 'rh'
        ! print *, velocity, 'velocity'
        ! print *, flowrate, 'flowrate'
        ! print *, slope, 'slope'
        ! print *, upperZ, 'upperZ', lowerZ, 'lowerZ'
        ! print *, total_length, 'total_length'
        ! print *, slope*total_length, 'slope*total_length'
        ! print *,'-----------------'

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine orifice_setup
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine froude_driven_pipe_setup &
        (upperZ, area, fulldepth, flowrate, velocity, Froude,  breadth, &
        topwidth, ManningsN, total_length, lowerZ, depth)

        character(64) :: subroutine_name = 'froude_driven_pipe_setup'

        real,  intent(out)    :: area, topwidth, flowrate, velocity, upperZ
        real,  intent(in)     :: Froude, breadth, ManningsN, lowerZ, fulldepth
        real,  intent(in)     :: depth, total_length

        real :: perimeter, rh, slope, AoverAfull, YoverYfull

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name


        YoverYfull  = depth / fulldepth

        if (YoverYfull .GT. 1) then
            print*, 'Warning: Initial water depth is greater than pipe full depth'
        endif

        area        = (onefourthR * pi * fulldepth ** twoR) * table_lookup(YoverYfull, ACirc, NACirc)
        topwidth    = fulldepth * table_lookup(YoverYfull, WCirc, NWCirc)
        rh          = onefourthR * fulldepth * table_lookup(YoverYfull, RCirc, NRCirc)
        perimeter   = area / rh
        velocity    = Froude * sqrt(grav * depth)
        flowrate    = area * velocity
        slope       = (velocity * ManningsN / (rh**(2.0/3.0)) )**2
        upperZ      = lowerZ + slope * total_length

        print *,'-----------------'
        print *, area, 'area'
        print *, topwidth, 'topwidth'
        print *, perimeter, 'perimeter'
        print *, rh, 'rh'
        print *, velocity, 'velocity'
        print *, flowrate, 'flowrate'
        print *, slope, 'slope'
        print *, upperZ, 'upperZ', lowerZ, 'lowerZ'
        print *, total_length, 'total_length'
        print *, slope*total_length, 'slope*total_length'
        print *,'-----------------'

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine froude_driven_pipe_setup
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine Initial_condition_for_width_depth_system &
        (upperZ, area, flowrate, velocity,  &
        Froude,  breadth, topwidth, ManningsN, total_length, &
        lowerZ, depth, widthDepthData)

        character(64) :: subroutine_name = 'Initial_condition_for_width_depth_system'

        real,  intent(out)    :: area(:), flowrate(:), velocity(:), upperZ(:), topwidth(:)
        real,  intent(in)     :: Froude(:),  breadth(:), ManningsN(:), lowerZ(:), total_length(:)
        real,  intent(in)     :: depth(:)
        !
        real,    target, intent(in out)    :: widthDepthData(:,:,:)

        real :: perimeter(size(depth,1)), rh(size(depth,1)), slope(size(depth,1))
        real :: AA, BB, CC, DD
        integer :: ind, ii
        real :: temp1(size(depth,1)), hDepth(size(depth,1))

        real, pointer :: widthAtLayerTop(:,:), depthAtLayerTop(:,:), areaThisLayer(:,:)
        real, pointer :: areaTotalBelowThisLayer(:,:), dWidth(:,:)
        real, pointer :: dDepth(:,:), angle(:,:), perimeterBelowThisLayer(:,:)
        real, pointer :: area_difference(:,:), local_difference(:,:), depthTBL(:,:)

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        widthAtLayerTop         => widthDepthData (:,:, wd_widthAtLayerTop)
        depthAtLayerTop         => widthDepthData (:,:, wd_depthAtLayerTop)
        areaThisLayer           => widthDepthData (:,:, wd_areaThisLayer)
        areaTotalBelowThisLayer => widthDepthData (:,:, wd_areaTotalBelowThisLayer)
        depthTBL                => widthDepthData (:,:, wd_depthTotalBelowThisLayer)
        dWidth                  => widthDepthData (:,:, wd_Dwidth)
        dDepth                  => widthDepthData (:,:, wd_Ddepth)
        angle                   => widthDepthData (:,:, wd_angle)
        perimeterBelowThisLayer => widthDepthData (:,:, wd_perimeterBelowThisLayer)
        area_difference         => widthDepthData (:,:, wd_area_difference)
        local_difference        => widthDepthData (:,:, wd_local_difference)

        do ii= 1, N_link
            temp1(:) = depthTBL(ii,:)-depth(ii)
            ind = minloc(abs(temp1(:)), DIM=1)

            area (ii) =  areaTotalBelowThisLayer(ii, ind)
            AA = oneR/tan(angle(ii,ind))
            BB = widthAtLayerTop(ii,ind) - dWidth(ii,ind)
            CC = - area_difference(ii,ind)
            DD = (-BB + sqrt(BB**twoR - fourR*AA*CC))/(twoR*AA)

            hdepth(ii) = DD + depthAtLayerTop(ii,ind) - dDepth(ii,ind)
            topwidth(ii) = widthAtLayerTop(ii,ind) - (dDepth(ii,ind)-DD) &
                *dWidth(ii,ind)/dDepth(ii,ind)
            perimeter(ii) = perimeterBelowThisLayer(ii,ind) + twoR * DD/sin(angle(ii,ind))
        enddo
        rh = area / perimeter
        velocity = Froude * sqrt(grav * hDepth)
        flowrate = area * velocity
        slope = (velocity * ManningsN / (rh**(2.0/3.0)) )**2
        upperZ = lowerZ + slope * total_length

        !  do ii= 1, N_link
        !     print *, "elelment No. =", ii
        !     print *, area(ii)
        !     print *, perimeter(ii)
        !     print *, rh(ii)
        !     print *, velocity(ii)
        !     print *, flowrate(ii)
        !     print *, slope(ii)
        !     print *, upperZ(ii), lowerZ(ii)
        !     print *, total_length(ii)
        !     print *, slope(ii)*total_length(ii)
        !  enddo

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine Initial_condition_for_width_depth_system
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine write_testcase_setup_file &
        (Froude, CFL, flowrate, velocity, init_depth, depth_upstream, depth_dnstream, breadth,  &
        topwidth, area, total_length, subdivide_length, lowerZ, upperZ, ManningsN)

        character(64) :: subroutine_name = ' write_testcase_setup_file'

        real,  intent(in)  :: CFL
        real,  intent(in)  :: Froude(:),  flowrate(:), velocity(:),  breadth(:), topwidth(:)
        real,  intent(in)  :: area(:), total_length(:), subdivide_length(:), lowerZ(:), upperZ(:)
        real,  intent(in)  :: ManningsN(:), depth_upstream(:), depth_dnstream(:), init_depth(:)

        integer        :: UnitNumber

        character(64)  :: thisFilePath, thisFileStatus, thisFileName
        character(256) :: thisFileWriteName

        logical        :: thisFileisOpen = .false.

        integer                :: open_status
        character(len=512)     :: emsg

        !--------------------------------------------------------------------------
        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ',subroutine_name

        open_status = 0

        UnitNumber = outputfile_next_unitnumber
        outputfile_next_unitnumber = outputfile_next_unitnumber+1

        thisFileName   = trim(setting%TestCase%TestName)
        thisFilePath   = trim(setting%DebugOut%FolderPath) &
            // trim(setting%Debugout%FolderName) // '/'
        thisFileStatus = 'new'
        thisFileIsOpen     = .true.

        thisFileWriteName  = trim(thisFilePath) // &
            trim(thisFileName) // &
            trim(setting%Time%DateTimeStamp) //&
            '.txt'

        ! print *, trim(setting%TestCase%TestName)
        ! print *, trim(setting%DebugOut%FolderPath)
        ! print *, trim(setting%Debugout%FolderName)
        !
        ! print *, trim(thisFileName)
        ! print *, trim(thisFilePath)
        ! print *, trim(thisFileWriteName)
        ! stop
        !
        open(unit=UnitNumber, &
            file=trim(thisFileWriteName), &
            status = 'new', &
            access = 'sequential', &
            form   = 'formatted', &
            action = 'write', &
            iostat = open_status)

        emsg = 'file exists: file open failed in '//trim(subroutine_name) &
            // '; filename = '//trim(thisFileWriteName)
        call utility_check_fileopen (open_status, emsg)

        write(UnitNumber,*) trim(setting%TestCase%TestName)
        write(UnitNumber,*) trim(setting%Time%DateTimeStamp)
        write(UnitNumber,*)
        write(UnitNumber,*) Froude  ,'=Froude'
        write(UnitNumber,*) CFL     ,'=CFL combined'
        write(UnitNumber,*) velocity * setting%Time%Dt / subdivide_length,'=CFL advective'
        write(UnitNumber,*) sqrt(grav * init_depth) * setting%Time%DT / subdivide_length,'=CFL barotropic'
        write(UnitNumber,*)
        write(UnitNumber,*) flowrate, '=flowrate'
        write(UnitNumber,*) velocity, '=velocity'
        write(UnitNumber,*) setting%Time%Dt,' = dt'
        write(UnitNumber,*)
        write(UnitNumber,*) init_depth       ,'=init depth'
        write(UnitNumber,*) depth_upstream   ,'=depth upstream'
        write(UnitNumber,*) depth_dnstream   ,'=depth downstream'
        write(UnitNumber,*) breadth ,'=breadth'
        write(UnitNumber,*) topwidth ,'=topwidth'
        write(UnitNumber,*) area    ,'=area'
        write(UnitNumber,*) total_length ,'=total_length'
        write(UnitNumber,*) subdivide_length ,'=subdivide_length'
        write(UnitNumber,*) area * subdivide_length,'=element_volume'
        write(UnitNumber,*) lowerZ,'=lowerZ'
        write(UnitNumber,*) upperZ,'=upperZ'
        write(UnitNumber,*) (upperZ - lowerZ )/ total_length,'=slope'
        write(UnitNumber,*)
        write(UnitNumber,*) ManningsN,'=ManningsN'
        write(UnitNumber,*)
        write(UnitNumber,*) setting%Step%First,'=first step'
        write(UnitNumber,*) setting%Step%Final,'=last step'
        write(UnitNumber,*) setting%Time%StartTime,'=start time'
        write(UnitNumber,*) setting%Time%EndTime,'=end time'
        write(UnitNumber,*)

        close(UnitNumber)
        outputfile_next_unitnumber = outputfile_next_unitnumber-1

        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ',subroutine_name
    end subroutine  write_testcase_setup_file
    !
    !==========================================================================
    !==========================================================================
    !
    elemental function get_dt_from_CFL &
        (CFL, velocity, depth, element_length) &
        result (dt)

        ! character(64) :: subroutine_name = 'get_dt_from_CFL'

        real,  intent(in) :: CFL, velocity, depth, element_length
        real :: dt

        !--------------------------------------------------------------------------

        dt = CFL * onehalfR * element_length / (velocity + sqrt(grav * depth))

    end function get_dt_from_CFL
    !
    !==========================================================================
    ! END OF MODULE test_cases
    !==========================================================================
end module test_cases
