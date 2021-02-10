module selectors

    use dynamic_array

    implicit none

    ! Temporal dynamics arrays
    type(integer_array) :: nodes_with_extinflow
    type(integer_array) :: nodes_with_dwfinflow
    type(integer_array) :: nodes_with_inflow

end module selectors