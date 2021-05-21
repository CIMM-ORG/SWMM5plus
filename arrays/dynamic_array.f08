module dynamic_array

    use type_definitions
    use setting_definition, only: setting

    implicit none

    integer :: MAX_DARRAY_SIZE = 10

contains

    subroutine dyna_real_append(this, x)
    !-----------------------------------------------------------------------------
    !
    ! Description:
    !
    !
    ! Method:
    !    
    !
    !-----------------------------------------------------------------------------
        type(real_array), intent(inout) :: this
        real(8), intent(in) :: x
        real(8), allocatable :: resized_arr(:)
        character(64) :: subroutine_name = 'dyna_real_append'

        if (setting%Debug%File%dynamic_array) print *, '*** enter ',subroutine_name

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
        if (setting%Debug%File%dynamic_array)  print *, '*** leave ',subroutine_name
    end subroutine dyna_real_append

    subroutine dyna_real_extend(this, a)
    !-----------------------------------------------------------------------------
    !
    ! Description:
    !
    !
    ! Method:
    !    
    !
    !-----------------------------------------------------------------------------
        type(real_array), intent(inout):: this
        type(real_array), intent(in) :: a
        integer :: i
        character(64) :: subroutine_name = 'dyna_real_extend'

        if (setting%Debug%File%dynamic_array) print *, '*** enter ',subroutine_name
        do i = 1, a%len
            call dyna_real_append(this, a%array(i))
        end do
        if (setting%Debug%File%dynamic_array)  print *, '*** leave ',subroutine_name
    end subroutine dyna_real_extend

    subroutine dyna_integer_append(this, x)
    !-----------------------------------------------------------------------------
    !
    ! Description:
    !
    !
    ! Method:
    !    
    !
    !-----------------------------------------------------------------------------
        type(integer_array), intent(inout) :: this
        integer, intent(in) :: x
        integer, allocatable :: resized_arr(:)
        character(64) :: subroutine_name = 'dyna_integer_append'

        if (setting%Debug%File%dynamic_array) print *, '*** enter ',subroutine_name

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
        if (setting%Debug%File%dynamic_array)  print *, '*** leave ',subroutine_name
    end subroutine dyna_integer_append

    subroutine dyna_integer_extend(this, a)
    !-----------------------------------------------------------------------------
    !
    ! Description:
    !
    !
    ! Method:
    !    
    !
    !-----------------------------------------------------------------------------
        type(integer_array), intent(inout):: this
        type(integer_array), intent(in) :: a
        integer :: i
        character(64) :: subroutine_name = 'dyna_integer_extend'

        if (setting%Debug%File%dynamic_array) print *, '*** enter ',subroutine_name
        do i = 1, a%len
            call dyna_integer_append(this, a%array(i))
        end do
        if (setting%Debug%File%dynamic_array)  print *, '*** leave ',subroutine_name
    end subroutine dyna_integer_extend

    function dyna_real_pop(this)
        type(real_array), intent(inout) :: this
        real(8) :: dyna_real_pop
        character(64) :: subroutine_name = 'dyna_real_pop'

        if (setting%Debug%File%dynamic_array) print *, '*** enter ',subroutine_name
        dyna_real_pop = this%array(this%len)
        this%len = this%len - 1
        if (setting%Debug%File%dynamic_array)  print *, '*** leave ',subroutine_name
    end function dyna_real_pop

    function dyna_integer_pop(this)
        type(integer_array), intent(inout) :: this
        integer :: dyna_integer_pop
        character(64) :: subroutine_name = 'dyna_integer_pop'

        if (setting%Debug%File%dynamic_array) print *, '*** enter ',subroutine_name
        dyna_integer_pop = this%array(this%len)
        this%len = this%len - 1
        if (setting%Debug%File%dynamic_array)  print *, '*** leave ',subroutine_name
    end function dyna_integer_pop

    subroutine free_real_array(this)
    !-----------------------------------------------------------------------------
    !
    ! Description:
    !
    !
    ! Method:
    !    
    !
    !-----------------------------------------------------------------------------
        type(real_array), intent(inout) :: this
        character(64) :: subroutine_name = 'free_real_array'

        if (setting%Debug%File%dynamic_array) print *, '*** enter ',subroutine_name
        deallocate(this%array)
        if (setting%Debug%File%dynamic_array)  print *, '*** leave ',subroutine_name
    end subroutine

    subroutine free_integer_array(this)
    !-----------------------------------------------------------------------------
    !
    ! Description:
    !
    !
    ! Method:
    !    
    !
    !-----------------------------------------------------------------------------
        type(integer_array), intent(inout) :: this
        character(64) :: subroutine_name = 'free_integer_array'

        if (setting%Debug%File%dynamic_array) print *, '*** enter ',subroutine_name
        deallocate(this%array)
        if (setting%Debug%File%dynamic_array)  print *, '*** leave ',subroutine_name
    end subroutine
end module
