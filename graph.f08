module network_graph

    implicit none

    type graph_node
        integer :: node_id
        type(graph_node), pointer :: next_node
    end type graph_node

    type graph
        type(graph_node), allocatable, dimension(:) :: g ! graph linked lists
        integer :: num_vertices
    end type graph

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
        integer :: i
        deallocate(g%g)
    end subroutine



end module network_graph