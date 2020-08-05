module linked_list

    use globals

    implicit none

    type node
        type(node), pointer :: previous => null()
        type(node), pointer :: next => null()
        double precision :: x
        double precision :: y
    contains
        final :: node_deallocator
    end type node

    type list
        integer :: len = 0
        type(node), pointer :: head => null()
        type(node), pointer :: tail => null()
    contains
        procedure :: append
    end type list

    integer, private :: debuglevel = 0

contains

    elemental subroutine node_deallocator(this)
        type(node), intent(inout) :: this
        if (associated(this%next)) nullify(this%next)
    end subroutine node_deallocator
    subroutine append(this, x, y)
        class(list) :: this
        double precision :: x, y
        type(node), pointer :: new_node_ref
        character(64) :: subroutine_name =  "append (linked_list.f08)"

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ', subroutine_name

        allocate(new_node_ref)
        new_node_ref%x = x
        new_node_ref%y = y

        if (this%len == 0) then
            this%head => new_node_ref
            this%tail => new_node_ref
        else if (this%len == 1) then
            this%head%next => new_node_ref
            this%tail => new_node_ref
            new_node_ref%previous => this%head
        else
            new_node_ref%previous => this%tail
            this%tail => new_node_ref
        end if

        this%len = this%len + 1

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ', subroutine_name

    end subroutine append
end module linked_list
