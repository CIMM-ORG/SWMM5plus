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
    integer :: ii 

    do ii = 1, size(P_nodeI,1)
        print*, P_nodeI(ii, :)
    end do

end subroutine default_partitioning


end module partitioning