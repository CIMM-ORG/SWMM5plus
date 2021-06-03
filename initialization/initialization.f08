module initialization
    use interface
    use partitioning
    use discretization
    use utility_allocate
    use utility_array
    use define_indexes
    use define_keys
    use define_globals
    use network_define
    use utility, only: utility_export_linknode_csv
    use define_settings, only: setting

    implicit none

!-----------------------------------------------------------------------------
!
! Description:
!   General initialization of data structures (not including network)
!
! Method:
!    Creates the arrays index structures that are used for accessing data.
!    Arguably, this could be done more simply but we want the fundamental
!    column indexes in array_index to be parameters rather than variables. By
!    using parameters we reduce the possibility of accidentally changing a
!    column definition.
!
!-----------------------------------------------------------------------------

    private

    public :: initialize_all

contains

    subroutine initialize_all()

        call load_settings(setting%Paths%setting)
        if (this_image() == 1) then
            call execute_command_line("if [ -d debug ]; then rm -r debug; fi && mkdir debug")
        end if

        call read_arguments()

        if (setting%Verbose) print *, "Simulation Starts"

        call initialize_api()
        call initialize_linknode_arrays()
        call initialize_partition_coarray()

        sync all
        call network_initiation()
        sync all

    end subroutine initialize_all

    subroutine initialize_linknode_arrays()
    !-----------------------------------------------------------------------------
    !
    ! Description:
    !   Retrieves data from SWMM C interface and populates link and node tables
    !
    !-----------------------------------------------------------------------------

        integer       :: i, total_n_links
        logical       :: l1, l2
        character(64) :: subroutine_name = 'initialize_linknode_arrays'

    !-----------------------------------------------------------------------------

        if (setting%Debug%File%initialization) print *, '*** enter ', subroutine_name

        if (.not. api_is_initialized) then
            print *, "ERROR: API is not initialized"
            stop
        end if

        ! Allocate storage for link & node tables
        call allocate_linknode_storage()

        nodeI(:,ni_N_link_u) = 0
        nodeI(:,ni_N_link_d) = 0

        do i = 1, N_link
            linkI(i,li_idx) = i
            linkI(i,li_link_type) = get_link_attribute(i, link_type)
            linkI(i,li_geometry) = get_link_attribute(i, link_geometry)
            linkI(i,li_Mnode_u) = get_link_attribute(i, link_node1) + 1 ! node1 in C starts from 0
            linkI(i,li_Mnode_d) = get_link_attribute(i, link_node2) + 1 ! node2 in C starts from 0

            ! HACK This is a temporary hardcode until Gerardo can populate this column from the CFL condition
            linkI(i, li_N_element) = 10

            nodeI(linkI(i,li_Mnode_d), ni_N_link_u) = nodeI(linkI(i,li_Mnode_d), ni_N_link_u) + 1
            nodeI(linkI(i,li_Mnode_d), ni_idx_base1 + nodeI(linkI(i,li_Mnode_d), ni_N_link_u)) = i
            nodeI(linkI(i,li_Mnode_u), ni_N_link_d) = nodeI(linkI(i,li_Mnode_u), ni_N_link_d) + 1
            nodeI(linkI(i,li_Mnode_u), ni_idx_base2 + nodeI(linkI(i,li_Mnode_u), ni_N_link_d)) = i

            linkI(i,li_InitialDepthType) = 1 ! TODO - get from params file
            linkR(i,lr_Length) = get_link_attribute(i, conduit_length)

            ! linkR(i,lr_TopWidth): defined in network_define.f08
            linkR(i,lr_BreadthScale) = get_link_attribute(i, link_xsect_wMax)
            ! linkR(i,lr_Slope): defined in network_define.f08
            linkR(i,lr_LeftSlope) = get_link_attribute(i, link_left_slope)
            linkR(i,lr_RightSlope) = get_link_attribute(i, link_right_slope)
            linkR(i,lr_Roughness) = get_link_attribute(i, conduit_roughness)
            linkR(i,lr_InitialFlowrate) = get_link_attribute(i, link_q0)
            linkR(i,lr_InitialUpstreamDepth) = get_node_attribute(linkI(i,li_Mnode_u), node_initDepth)
            linkR(i,lr_InitialDnstreamDepth) = get_node_attribute(linkI(i,li_Mnode_d), node_initDepth)
            linkR(i,lr_InitialDepth) = (linkR(i,lr_InitialDnstreamDepth) + linkR(i,lr_InitialUpstreamDepth)) / 2.0
        end do

        do i = 1, N_node
            total_n_links = nodeI(i,ni_N_link_u) + nodeI(i,ni_N_link_d)
            nodeI(i, ni_idx) = i
            if (get_node_attribute(i, node_type) == oneI) then ! OUTFALL
                nodeI(i, ni_node_type) = nBCdn
            else if (total_n_links == twoI) then
                nodeI(i, ni_node_type) = nJ2
            else if (total_n_links > twoI) then
                nodeI(i, ni_node_type) = nJm
            end if
            ! Nodes with nBCup are defined in inflow.f08 -> (inflow_load_inflows)
            l1 = get_node_attribute(i, node_has_extInflow) == 1
            l2 = get_node_attribute(i, node_has_dwfInflow) == 1
            if (l1 .or. l2) then
                !nodeYN(i, nYN_has_inflow) = .true.
                if (total_n_links == 1) then
                    nodeI(i, ni_node_type) = nBCup
                end if
            end if


            nodeR(i,nr_InitialDepth) = get_node_attribute(i, node_initDepth)
            nodeR(i,nr_Zbottom) = get_node_attribute(i, node_invertElev)
        end do

        if (setting%Debug%File%initialization) then
            if (this_image() == 1) then
                call utility_export_linknode_csv()
            end if
        end if

        if (setting%Debug%File%initialization)  print *, '*** leave ', subroutine_name
    end subroutine initialize_linknode_arrays
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine initialize_partition_coarray()
    !-----------------------------------------------------------------------------
    !
    ! Description:
    !   This subroutine calls the public subroutine from the utility module,
    !   partitioning.f08. It also calls a public subroutine from the temporary
    !   coarray_partition.f08 utility module that defines how big the coarrays
    !   must be.
    !
    !-----------------------------------------------------------------------------
        character(64) :: subroutine_name = 'initialize_partition'
    !-----------------------------------------------------------------------------

        if (setting%Debug%File%initialization) print *, '*** enter ', subroutine_name
        ! adjust the length and calculate the number/length of elements in each link
        call link_length_adjust()
        call nominal_discretization()
        !% In order to keep the main() clean, move the following two subroutines here, BIPquick can be removed
        call execute_partitioning()
        call coarray_length_calculation()
        call allocate_elemX_faceX()
        call allocate_columns()

        if (setting%Debug%File%initialization)  print *, '*** leave ', subroutine_name

    end subroutine initialize_partition_coarray

    subroutine read_arguments()
        integer :: ii
        logical :: arg_param = .false.
        character(len=8) :: param
        character(len=256) :: arg

        do ii = 1, iargc()
            call getarg(ii, arg)
            if (.not. arg_param) then
                param = arg
                if (ii == 1) then
                    if (arg(:1) == '-') then
                        print *, "ERROR: it is necessary to define the path to the .inp file"
                        stop
                    end if
                    setting%Paths%inp = arg
                elseif ((trim(arg) == "-s") .or. & ! user provides settings file
                        ((trim(arg) == "--test") .or. (trim(arg) == "-t"))) then  ! hard coded test case
                    arg_param = .true.
                elseif ((trim(arg) == '--verbose') .or. (trim(arg) == "-v")) then
                    setting%Verbose = .true.
                elseif ((trim(arg) == '--warnings-off') .or. (trim(arg) == "-woff")) then
                    setting%Warning = .false.
                else
                    write(*, *) 'The argument ' // trim(arg) // ' is unsupported'
                    stop
                end if
            else
                arg_param = .false.
                if (trim(param) == '-s') then
                    setting%Paths%setting = arg
                elseif ((trim(arg) == "--test") .or. (trim(arg) == "-t")) then
                    setting%TestCase%UseTestCase = .true.
                    setting%TestCase%TestName = trim(arg)
                    if (trim(arg) == 'simple_channel') then
                    else if (trim(arg) == 'simple_orifice') then
                    else
                        write(*, *) 'The test case ' // trim(arg) // ' is unsupported. Please use one of the following:'
                        print *, new_line('')
                        print *, "simple_channel, simple_orifice, simple_pipe"
                        print *, "simple_weir, swashes, waller_creek"
                        print *, "y_channel, y_storage_channel"
                        stop
                    end if
                elseif (trim(param) == '--run-tests') then
                end if
            end if
        end do
    end subroutine read_arguments

end module initialization
