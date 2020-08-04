module linked_list

    use globals

    implicit none

    type node
        type(node), pointer :: previous => null()
        type(node), pointer :: next => null()
        real :: item
    contains
        final :: node_deallocator
    end type node

    type list
        integer :: len = 0
        integer :: i = 0
        type(node), pointer :: root => null()
        type(node), pointer :: current => null()
    end type list

    integer, private :: debuglevel = 0

contains

    elemental subroutine node_deallocator(this)
        type(node), intent(inout) :: this
        if (associated(this%next)) nullify(this%next)
    end subroutine node_deallocator
    subroutine append(this, x)
        type(list), intent(inout) :: this
        real, intent(in) :: x
        type(node), target :: new_node
        type(node), pointer :: new_node_ref
        character(64) :: subroutine_name =  "append (linked_list)"

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ', subroutine_name

        allocate(new_node_ref)
        if (.not. associated(this%root)) then
            allocate(this%root)
        endif
        if (.not. associated(this%current)) then
            allocate(this%current)
        endif

        new_node_ref => new_node
        new_node%item = x

        if (this%len == 0) then
            this%root = new_node_ref
            this%current = this%root
            this%i = this%i + 1
        else
            if (.not. associated(this%current%next)) then
                allocate(this%current%next)
            endif
            this%current%next = new_node_ref
        end if

        this%len = this%len + 1

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ', subroutine_name

    end subroutine append

    subroutine next(this)
        type(list), intent(inout) :: this
        character(64) :: subroutine_name =  "next (linked_list)"

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** enter ', subroutine_name
        if (this%i == this%len) then
            write(*,*) "Error (i > LEN)"
            stop
        end if
        this%current = this%current%next
        this%i = this%i + 1

        if ((debuglevel > 0) .or. (debuglevelall > 0)) print *, '*** leave ', subroutine_name
    end subroutine next
end module linked_list
