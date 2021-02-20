module dynamic_array

    use type_definitions

    implicit none

    integer :: MAX_DARRAY_SIZE = 10

    integer, private :: debuglevel = 0

contains

    subroutine dyna_real_append(this, x)
        type(real_array), intent(inout) :: this
        real(4), intent(in) :: x
        real(4), allocatable :: resized_arr(:)
        character(64) :: subroutine_name = 'dyna_real_append'

        if (debuglevel > 0) print *, '*** enter ',subroutine_name

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
        if (debuglevel > 0)  print *, '*** leave ',subroutine_name
    end subroutine dyna_real_append

    subroutine dyna_real_extend(this, a)
        type(real_array), intent(inout):: this
        type(real_array), intent(in) :: a
        integer :: i
        character(64) :: subroutine_name = 'dyna_real_extend'

        if (debuglevel > 0) print *, '*** enter ',subroutine_name
        do i = 1, a%len
            call dyna_real_append(this, a%array(i))
        end do
        if (debuglevel > 0)  print *, '*** leave ',subroutine_name
    end subroutine dyna_real_extend

    subroutine dyna_integer_append(this, x)
        type(integer_array), intent(inout) :: this
        integer, intent(in) :: x
        integer, allocatable :: resized_arr(:)
        character(64) :: subroutine_name = 'dyna_integer_append'

        if (debuglevel > 0) print *, '*** enter ',subroutine_name

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
        if (debuglevel > 0)  print *, '*** leave ',subroutine_name
    end subroutine dyna_integer_append

    subroutine dyna_integer_extend(this, a)
        type(integer_array), intent(inout):: this
        type(integer_array), intent(in) :: a
        integer :: i
        character(64) :: subroutine_name = 'dyna_integer_extend'

        if (debuglevel > 0) print *, '*** enter ',subroutine_name
        do i = 1, a%len
            call dyna_integer_append(this, a%array(i))
        end do
        if (debuglevel > 0)  print *, '*** leave ',subroutine_name
    end subroutine dyna_integer_extend

    function dyna_real_pop(this)
        type(real_array), intent(inout) :: this
        real(4) :: dyna_real_pop
        character(64) :: subroutine_name = 'dyna_real_pop'

        if (debuglevel > 0) print *, '*** enter ',subroutine_name
        dyna_real_pop = this%array(this%len)
        this%len = this%len - 1
        if (debuglevel > 0)  print *, '*** leave ',subroutine_name
    end function dyna_real_pop

    function dyna_integer_pop(this)
        type(integer_array), intent(inout) :: this
        integer :: dyna_integer_pop
        character(64) :: subroutine_name = 'dyna_integer_pop'

        if (debuglevel > 0) print *, '*** enter ',subroutine_name
        dyna_integer_pop = this%array(this%len)
        this%len = this%len - 1
        if (debuglevel > 0)  print *, '*** leave ',subroutine_name
    end function dyna_integer_pop

    subroutine free_real_array(this)
        type(real_array), intent(inout) :: this
        character(64) :: subroutine_name = 'free_real_array'

        if (debuglevel > 0) print *, '*** enter ',subroutine_name
        deallocate(this%array)
        if (debuglevel > 0)  print *, '*** leave ',subroutine_name
    end subroutine

    subroutine free_integer_array(this)
        type(integer_array), intent(inout) :: this
        character(64) :: subroutine_name = 'free_integer_array'

        if (debuglevel > 0) print *, '*** enter ',subroutine_name
        deallocate(this%array)
        if (debuglevel > 0)  print *, '*** leave ',subroutine_name
    end subroutine
end module
