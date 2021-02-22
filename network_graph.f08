module network_graph

    use errors
    use interface
    use dynamic_array
    use globals
    use array_index
    use data_keys
    use setting_definition
    use inflow
    use type_definitions

    implicit none

    integer, private, parameter :: debuglevel = 0
contains

    function new_graph(num_vertices)
        integer, intent(in) :: num_vertices
        type(graph) :: new_graph
        character(64) :: subroutine_name  = 'new_graph'
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** enter ', subroutine_name
        new_graph%num_vertices = num_vertices
        allocate(new_graph%g(num_vertices))
        allocate(new_graph%in_degree(num_vertices))
        new_graph%in_degree(:) = 0
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ', subroutine_name
    end function

    subroutine add_graph_link(g, source, destination, link_id)
        type(graph), intent(inout) :: g
        integer, intent(in) :: source, destination, link_id
        character(64) :: subroutine_name  = 'add_graph_link'
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** enter ', subroutine_name
        call dyna_integer_append(g%g(source)%neighbors, destination)
        call dyna_integer_append(g%g(source)%link_id, link_id)
        call dyna_real_append(g%g(source)%neighbor_flows, dble(0.0))
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ', subroutine_name
    end subroutine add_graph_link

    subroutine free_graph(g)
        type(graph), intent(inout) :: g
        character(64) :: subroutine_name  = 'free_graph'
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** enter ', subroutine_name
        deallocate(g%g)
        deallocate(g%in_degree)
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ', subroutine_name
    end subroutine free_graph

    function get_network_graph()
        type(graph) :: get_network_graph
        integer :: i, src, dest
        character(64) :: subroutine_name  = 'get_network_graph'
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** enter ', subroutine_name

        if (.not. api_is_initialized) then
            print *, MSG_API_NOT_INITIALIZED
            stop
        end if

        get_network_graph = new_graph(num_nodes)

        do i = 1, num_links
            src = int(get_link_attribute(i, link_node1)) + 1
            dest = int(get_link_attribute(i, link_node2)) + 1
            get_network_graph%in_degree(dest) = get_network_graph%in_degree(dest) + 1
            call add_graph_link(get_network_graph, src, dest, i)
        end do
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ', subroutine_name
    end function get_network_graph

    subroutine traverse_graph_flow(g, i, flow)
        type(graph), intent(inout) :: g ! graph
        integer, intent(in) :: i ! root node
        real(8), intent(in) :: flow ! flow value

        type(integer_array) :: n_nodes ! stack of nodes that haven't been traversed
        integer :: k
        character(64) :: subroutine_name  = 'traverse_graph_flow'
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** enter ', subroutine_name

        call dyna_integer_append(n_nodes, i)
        do while (n_nodes%len > 0)
            k = dyna_integer_pop(n_nodes)
            if (g%g(k)%neighbor_flows%len == 0) cycle
            g%g(k)%neighbor_flows%array = g%g(k)%neighbor_flows%array + flow
            call dyna_integer_extend(n_nodes, g%g(k)%neighbors)
        end do
        if ((debuglevel > 0) .or. (debuglevelall > 0))  print *, '*** leave ', subroutine_name
    end subroutine traverse_graph_flow

    subroutine traverse_cfl_condition(g, linkR, nodeR, linkI, nodeI)
        type(graph), intent(inout) :: g
        integer, allocatable, target, intent(inout) :: linkI(:,:)
        integer, allocatable, target, intent(inout) :: nodeI(:,:)
        real(8), allocatable, target, intent(inout) :: linkR(:,:)
        real(8), allocatable, target, intent(inout) :: nodeR(:,:)
        integer :: i, j, link_id
        real(8) :: Q, N_R, SLP, Y, ML, MR, BT, TOL, A, P, F, DFDY, K0, K1, K2
        real(8), allocatable :: velocities(:)
        character(64) :: subroutine_name  = 'traverse_cfl_condition'
        if ((debuglevel > -1) .or. (debuglevelall > 0))  print *, '*** enter ', subroutine_name

        allocate(velocities(num_links))
        velocities(:) = 0

        do i = 1, g%num_vertices
            do j = 1, g%g(i)%neighbors%len
                link_id = g%g(i)%link_id%array(j)
                if ((linkI(link_id, li_link_type) .ne. lchannel) .and. (linkI(link_id, li_link_type) .ne. lpipe)) then
                    linkI(link_id, li_N_element) = 1
                    cycle
                end if
                Q = g%g(i)%neighbor_flows%array(j)
                print *, "max_inflow", Q
                N_R = linkR(link_id, lr_Roughness)
                SLP = linkR(link_id, lr_Slope)
                Y = 1.0
                TOL = 1000
                if (linkI(link_id, li_geometry) == lTrapezoidal) then
                    BT = linkR(link_id, lr_BreadthScale)
                    ML = linkR(link_id, lr_LeftSlope)
                    MR = linkR(link_id, lr_RightSlope)
                    K0 = SLP**0.5/N_R
                    K1 = ML + MR
                    K2 = (1.0+ML**2.0)**0.5 + (1.0+MR**2.0)**0.5
                    do while (ABS(TOL) > 1E-6)
                        F = K0*(Y*BT+Y**2.0*K1/2.0)**(5.0/3.0)/(BT+Y*K2)**(2.0/3.0) - Q
                        DFDY = 5.0*K0*(BT+K1*Y)*(BT*Y+K1*Y**2.0/2.0)**(2.0/3.0)/(3.0*(BT+K2*Y)**(2.0/3.0)) - &
                            2.0*K0*K2*(BT*Y+K1*Y**2.0/2.0)**(5.0/3.0)/(3.0*(BT+K2*Y)**(5.0/3.0))
                        ! print *, "F", F, "DDF", DFDY, "TOL", ABS(TOL), "Y", Y
                        TOL = - F/DFDY
                        Y = Y + TOL
                    end do
                    ! solve V with Manning
                    A = Y*BT + Y**2.0*K1/2.0
                    P = BT + Y*K2
                    velocities(link_id) = K0*(A/P)**(2.0/3.0)
                    ! print *, "Slope", SLP, "BREADTH", BT
                    ! print *, "K0", K0, "K1", K1, "K2", K2
                    ! print*, "ROughness", N_R, "L slope", ML, "R slope", MR
                    ! print *, "Wet Area", A, "Wet Perimeter", P, "Depth", Y, "VEL", velocities(link_id)
                else
                    print *, MSG_FEATURE_NOT_COMPATIBLE
                    stop
                end if
            end do
        end do

        do i = 1, num_links
            if (linkR(i, lr_Length) == 0) then
                cycle
            end if
            velocities(i) = linkR(i, lr_Length) / velocities(i) ! dt*num_elem
        end do

        setting%time%dt = minval(velocities) ! minimum dt
        setting%step%final = int(setting%time%endtime / setting%time%dt)

        linkI(:, li_N_element) = ceiling(velocities/setting%time%dt)
        linkR(:, lr_ElementLength) = linkR(:, lr_Length) / linkI(:, li_N_element)
        deallocate(velocities)

        print *, "Start Time", setting%time%starttime
        print *, "End Time", setting%time%endtime
        print *, "Time step", setting%time%dt
        print *, "Number of elements", sum(linkI(:, li_N_element))
        print *, "Minmax number of elements", minval(linkI(:, li_N_element)), maxval(linkI(:, li_N_element))
        print *, "Minmax inflows", minval(nodeR(:, nr_maxinflow)), maxval(nodeR(:, nr_maxinflow))
        
        if ((debuglevel > -1) .or. (debuglevelall > 0))  print *, '*** leave ', subroutine_name
    end subroutine traverse_cfl_condition
end module network_graph