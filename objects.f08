module objects

    use linked_list
    use tables

    implicit none
    public

    integer, dimension(2) :: nobjects ! vector with the count of objects
    type(table), dimension(:), allocatable :: Curve
    type(table), dimension(:), allocatable :: TSeries

end module objects