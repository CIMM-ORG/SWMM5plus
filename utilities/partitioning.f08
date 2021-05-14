module partitioning

    use array_index
    use data_keys
    use globals
    use setting_definition, only: setting
    use BIPquick

    implicit none

!***********************************************************************************
!
! Description:
!   This module controls the partitioning algorithm used.  Currently the options are
!       - BIPquick
!       - Default
!
!***********************************************************************************

    private

    public :: default_partitioning, partitioning_algorithm_check


contains

subroutine partitioning_algorithm_check()
    print*, setting%Partitioning%UseBIPquick, setting%Partitioning%UseDefault
    if ( (setting%Partitioning%UseBIPquick .eqv. .true.) .and. (setting%Partitioning%UseDefault .eqv. .true.) )  then
        print*, "There are two partitioning algorithms being used"
        stop
    else if ( (setting%Partitioning%UseBIPquick .eqv. .false.) .and. (setting%Partitioning%UseDefault .eqv. .false.) ) then
        print*, "No partitioning algorithms have been specified, default partitioning will be used"
        setting%Partitioning%UseDefault = .true.
    else
        if ( setting%Partitioning%UseBIPquick .eqv. .true. ) then
            print*, "Using BIPquick Partitioning"
        else if ( setting%Partitioning%UseDefault .eqv. .true. ) then
            print*, "Using Default Partitioning"
        end if
    end if
end subroutine partitioning_algorithm_check

subroutine default_partitioning()
    integer :: ii, num_nJm_nodes
    num_nJm_nodes = count_nJm_nodes()
    
    ! do ii = 1, size(P_nodeI,1)
    !     print*, P_nodeI(ii, :)
    ! end do

end subroutine default_partitioning

function count_nJm_nodes() result(num_nJm_nodes)
    integer :: num_nJm_nodes
    ! This subroutine iterates through the nodeI array and counts the number of instances of ni_node_type == nJm
    num_nJm_nodes = 1
    
end function count_nJm_nodes


end module partitioning