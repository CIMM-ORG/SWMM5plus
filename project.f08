module project

    use type_definitions
    use objects
    use interface
    use selectors
    use globals
    use data_keys
    use array_index
    use network_graph
    use initialization
    use bc

    implicit none

contains

    subroutine project_open(linkI, nodeI, linkR, nodeR, linkYN, nodeYN, linkName, nodeName, bcdataUp, bcdataDn)
        integer,   dimension(:,:), allocatable, target, intent(out)    :: linkI
        integer,   dimension(:,:), allocatable, target, intent(out)    :: nodeI
        real(8),      dimension(:,:), allocatable, target, intent(out)    :: linkR
        real(8),      dimension(:,:), allocatable, target, intent(out)    :: nodeR
        logical,   dimension(:,:), allocatable, target, intent(out)    :: linkYN
        logical,   dimension(:,:), allocatable, target, intent(out)    :: nodeYN
        type(string), dimension(:), allocatable, target, intent(out)   :: linkName
        type(string), dimension(:), allocatable, target, intent(out)   :: nodeName
        type(bcType), dimension(:), allocatable, target, intent(out)   :: bcdataUp, bcdataDn

        integer :: ii, jj, tssize

        ! --------------------
        ! --- Initialize C API
        ! --------------------

        call initialize_api()

        ! Retrieve system properties from SWMM C
        call initialize_linknode_arrays &
            (linkI, nodeI, linkR, nodeR, linkYN, nodeYN, linkName, nodeName)

        ! Load inflows from SWMM C
        call inflow_load_inflows(nodeI, nodeR, bcdataDn, bcdataUp)

        ! Create system graph
        swmm_graph = get_network_graph()

        call finalize_api()
        ! --------------------
        ! --- Finalize C API
        ! --------------------
    end subroutine project_open

    subroutine project_close(bcdataDn, bcdataUp)
        type(bcType), dimension(:), allocatable, target, intent(inout)   :: bcdataUp, bcdataDn
        call free_graph(swmm_graph)
        call free_interface()
        call free_bc(bcdataDn, bcdataUp)
    end subroutine project_close

end module project