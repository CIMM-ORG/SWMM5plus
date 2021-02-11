module network_graph

    use objects

    implicit none

contains
    function new_graph(num_vertices)
        integer, intent(in) :: num_vertices
        type(graph) :: new_graph
        allocate(new_graph%g(num_vertices))
    end function

    subroutine add_edge(g, source, destination)
        type(graph), intent(inout) :: g
        integer, intent(in) :: source, destination

        type(graph_node) :: new_node
        new_node%node_id = destination
        new_node%next_node = g%g(source)
        g%g(source) = new_node
    end subroutine

    subroutine free_graph(g)
        type(graph), intent(inout) :: g
        deallocate(g%g)
    end subroutine



end module network_graph