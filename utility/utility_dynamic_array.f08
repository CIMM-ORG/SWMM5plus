module util_dynamic_array
    type f_array
        integer :: max_size = 0
        integer :: len = 0
        real, allocatable :: arr(:)
    end type f_array
contains
    subroutine util_free_arr(this)
        type(f_array), intent(inout) :: this
        deallocate(this%arr)
    end subroutine

    subroutine append(this, x)
        type(f_array), intent(inout) :: this
        real, intent(in) :: x
        real, allocatable :: resized_arr(:)

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
    end subroutine append
end module

! program main
!     use Array
!     implicit none
!     type(f_array) :: x
!     integer i

!     do i = 1 , 25
!         call append(x, 2.0*i)
!         print *, x%arr(x%len), x%len, x%max_size
!     end do

!     call free_arr(x)
! end program main