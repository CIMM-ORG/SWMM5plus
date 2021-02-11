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
        real,      dimension(:,:), allocatable, target, intent(out)    :: linkR
        real,      dimension(:,:), allocatable, target, intent(out)    :: nodeR
        logical,   dimension(:,:), allocatable, target, intent(out)    :: linkYN
        logical,   dimension(:,:), allocatable, target, intent(out)    :: nodeYN
        type(string), dimension(:), allocatable, target, intent(out)   :: linkName
        type(string), dimension(:), allocatable, target, intent(out)   :: nodeName
        type(bcType), dimension(:), allocatable, target, intent(out)   :: bcdataUp, bcdataDn

        type(tseries) :: ts_ups
        integer :: ii, jj

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
            allocate(bcdataUp(ii)%TimeArray(2))
            allocate(bcdataUp(ii)%ValueArray(2))
            bcdataUp(ii)%TimeArray = (/0.0, real(setting%time%endtime)/)
            bcdataUp(ii)%ValueArray = 1
        enddo

        print *, "Setting up BC downstream"
        do ii = 1, N_BCdnstream
            print *, "BC dnstream", ii, '/', N_BCdnstream
            bcdataDn(ii)%NodeID = nodeI(ii, ni_temp1)
            allocate(bcdataDn(ii)%TimeArray(2))
            allocate(bcdataDn(ii)%ValueArray(2))
            bcdataDn(ii)%TimeArray = (/0.0, real(setting%time%endtime)/)
            bcdataDn(ii)%ValueArray = nr_Zbottom
        enddo
        ! --------------------
        ! --- Finalize C API
        ! --------------------
    end subroutine project_open

    subroutine project_close(bcdataDn, bcdataUp)
        type(bcType), dimension(:), allocatable, target, intent(out)   :: bcdataUp, bcdataDn
        call free_graph(swmm_graph)
        call free_interface()
        call free_bc(bcdataDn, bcdataUp)
    end subroutine project_close



end module project