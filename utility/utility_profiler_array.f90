module utility_prof_array

    !use define_types
    !use define_globals
    !use utility_allocate
    !use utility_deallocate
    type f_array
        integer :: max_size = 0
        integer :: len = 0
        real, allocatable :: arr(:)
    end type f_array
contains
    subroutine util_free_arr(this)
        !move to utility_deallocate
        type(f_array), intent(inout) :: this
        deallocate(this%arr)
    end subroutine

    subroutine util_prof_append(this, x)
        type(f_array), intent(inout) :: this
        real, intent(in) :: x
        real, allocatable :: resized_arr(:)

        !call util_allocate_profiler(this)
        if (this%max_size == 0) then
            allocate(this%arr(100))
            this%max_size = 100
        else if (this%len == this%max_size) then
            allocate(resized_arr(this%max_size * 2))
            resized_arr(1:this%max_size) = this%arr(1:this%max_size)
            this%max_size = this%max_size * 2
            call util_free_arr(this)
            this%arr = resized_arr
        end if

        this%len = this%len + 1
        this%arr(this%len) = x
    end subroutine util_prof_append
end module utility_prof_array

