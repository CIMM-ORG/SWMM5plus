module dynamic_array

    implicit none

    integer :: MAX_DARRAY_SIZE = 10

    type real_array
        integer :: max_size = 0
        integer :: len = 0
        real, allocatable :: array(:)
    end type real_array

    type integer_array
        integer :: max_size = 0
        integer :: len = 0
        integer, allocatable :: array(:)
    end type integer_array

contains

    subroutine dyna_real_append(this, x)
        type(real_array), intent(inout) :: this
        real, intent(in) :: x
        real, allocatable :: resized_arr(:)

        if (this%max_size == 0) then
            allocate(this%array(MAX_DARRAY_SIZE))
            this%max_size = MAX_DARRAY_SIZE
        else if (this%len == this%max_size) then
            allocate(resized_arr(this%max_size * 2))
            resized_arr(1:this%max_size) = this%array(1:this%max_size)
            this%max_size = this%max_size * 2
            call free_real_array(this)
            this%array = resized_arr
        end if

        this%len = this%len + 1
        this%array(this%len) = x
    end subroutine dyna_real_append

    subroutine dyna_real_extend(this, a)
        type(real_array), intent(inout):: this
        type(real_array), intent(in) :: a
        integer :: i
        do i = 1, a%len
            call dyna_real_append(this, a%array(i))
        end do
    end subroutine dyna_real_extend

    subroutine dyna_integer_append(this, x)
        type(integer_array), intent(inout) :: this
        integer, intent(in) :: x
        integer, allocatable :: resized_arr(:)

        if (this%max_size == 0) then
            allocate(this%array(MAX_DARRAY_SIZE))
            this%max_size = MAX_DARRAY_SIZE
        else if (this%len == this%max_size) then
            allocate(resized_arr(this%max_size * 2))
            resized_arr(1:this%max_size) = this%array(1:this%max_size)
            this%max_size = this%max_size * 2
            call free_integer_array(this)
            this%array = resized_arr
        end if

        this%len = this%len + 1
        this%array(this%len) = x
    end subroutine dyna_integer_append

    subroutine dyna_integer_extend(this, a)
        type(integer_array), intent(inout):: this
        type(integer_array), intent(in) :: a
        integer :: i
        do i = 1, a%len
            call dyna_integer_append(this, a%array(i))
        end do
    end subroutine dyna_integer_extend

    function dyna_real_pop(this)
        type(real_array), intent(inout) :: this
        real :: dyna_real_pop
        dyna_real_pop = this%array(this%len)
        this%len = this%len - 1
    end function dyna_real_pop

    function dyna_integer_pop(this)
        type(integer_array), intent(inout) :: this
        integer :: dyna_integer_pop
        dyna_integer_pop = this%array(this%len)
        this%len = this%len - 1
    end function dyna_integer_pop

    subroutine free_real_array(this)
        type(real_array), intent(inout) :: this
        deallocate(this%array)
    end subroutine

    subroutine free_integer_array(this)
        type(integer_array), intent(inout) :: this
        deallocate(this%array)
    end subroutine
end module
