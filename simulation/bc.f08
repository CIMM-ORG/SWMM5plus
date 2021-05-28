! module bc
!
! Setup and apply boundary conditions. All the BC functions and subroutines
! should be located in this module.
!
! Note that BC should only be applied on faces connected to an elem2 element.
! That is, a BC cannot be associated with a junction branch.
!
!==========================================================================
!
module bc
    !
    ! boundary condition definitions and enforcement
    !
    use array_index
    use data_keys
    use globals
    use setting_definition, only: setting
    use type_definitions
    use utility

    implicit none

    private

    public :: bc_allocate
    public :: free_bc

    integer :: debuglevel = 0

contains

    
    subroutine free_bc(bcdataDn, bcdataUp)
        type(bcType), dimension(:), allocatable, intent(inout)   :: bcdataUp, bcdataDn
        integer :: ii
        do ii = 1, N_BCdnstream
            deallocate(bcdataDn(ii)%TimeArray)
            deallocate(bcdataDn(ii)%ValueArray)
        end do
        do ii = 1, N_BCupstream
            deallocate(bcdataUp(ii)%TimeArray)
            deallocate(bcdataUp(ii)%ValueArray)
        end do
        deallocate(bcdataUp)
        deallocate(bcdataDn)
    end subroutine free_bc

end module bc