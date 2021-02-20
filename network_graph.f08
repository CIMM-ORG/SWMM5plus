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
        call dyna_real_append(g%g(source)%neighbor_flows, real(0.0))
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
        real(4), intent(in) :: flow ! flow value

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
        real(4), intent(inout) :: linkR(:,:)
        real(4), intent(in) :: nodeR(:,:)
        integer, intent(inout) :: linkI(:,:)
        integer, intent(in) :: nodeI(:,:)
        integer :: i, j, link_id
        real(4) :: Q, N_R, SLP, Y, ML, MR, BT, TOL, A, P, F, DDF
        real(4), allocatable :: velocities(:)
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
                N_R = linkR(link_id, lr_Roughness)
                SLP = linkR(link_id, lr_Slope)
                BT = linkR(link_id, lr_BreadthScale)
                Y = 0
                TOL = 1000
                if (linkI(link_id, li_geometry) == lTrapezoidal) then
                    ML = linkI(link_id, lr_LeftSlope)
                    MR = linkI(link_id, lr_RightSlope)
                    do while (ABS(TOL) > 10**(-6))
                        F = (SLP**(1/2)*(Y*BT + Y**2*ML/2 + Y**2*MR/2)**(5/3)) / &
                            (N_R*(BT+(Y**2+(Y*ML)**2)**0.5+(Y**2+(MR*Y)**2)**0.5)**(2/3)) - Q
                        DDF = 5*SLP**0.5*(BT + ML*Y + MR*Y)*(BT*Y + ML*Y**2/2 + MR*Y**2/2)**(2/3) / &
                            (3*N_R*(BT + ((ML*Y)**2 + Y**2)**0.5 + ((MR*Y)**2 + Y**2)**0.5)**(2/3)) - &
                            2*SLP**-.5*((2*ML**2*Y+2*Y)/(2*((ML*Y)**2+Y**2)**0.5) + (2*MR**2*Y+2*Y)/(2*((MR*Y)**2+Y**2)**0.5)) &
                            *(BT*Y+ML*Y**2/2+MR*Y**2/2)**(5/3) / (3*N_R*(BT+((ML*Y)**2 + Y**2)**0.5+((MR*Y)**2 + Y**2)**0.5)**(5/3))
                        TOL = - F/DDF
                        Y = Y + TOL
                    end do
                    ! solve V with Manning
                    A = Y*BT + Y**2*ML/2 + Y**2*MR/2
                    P = BT + (Y**2 + (Y*ML)**2)**0.5 + (Y**2 + (Y**MR)**2)**0.5
                    velocities(link_id) = 1/N_R * (A/P)**(2/3) * SLP**0.5
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