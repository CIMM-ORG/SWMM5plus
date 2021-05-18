!
! module network_define
!
! Handles relationship between coarse link-node network and high-resolution
! element-face network.
!
! This module defines all the indexes and mappings
!
!==========================================================================
!
module network_define
    !
    use allocate_storage
    use assign_index
    use initialization
    use data_keys
    use globals
    use interface

    implicit none

    private

    public :: network_initiation

contains
    !
    !==========================================================================
    !==========================================================================
    !
    subroutine network_initiation ()
    !-----------------------------------------------------------------------------
    !
    ! Initializes a element-face network from a link-node network.
    !   Requires network links and nodes before execution 
    !
    !-----------------------------------------------------------------------------

        character(64) :: subroutine_name = 'network_initiation'

    !-----------------------------------------------------------------------------

        if (setting%Debug%File%network_define) print *, '*** enter ',subroutine_name


        ! !%   get the slope of each link given the node Z values
        ! call network_get_link_slope

        if (setting%Debug%File%network_define) print *, '*** leave ',subroutine_name

    end subroutine network_initiation

end module network_define