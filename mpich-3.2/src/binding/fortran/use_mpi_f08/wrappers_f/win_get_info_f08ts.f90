!   -*- Mode: Fortran; -*-
!
!   (C) 2014 by Argonne National Laboratory.
!   See COPYRIGHT in top-level directory.
!
subroutine MPI_Win_get_info_f08(win, info_used, ierror)
    use, intrinsic :: iso_c_binding, only : c_int
    use :: mpi_f08, only : MPI_Win, MPI_Info
    use :: mpi_c_interface_types, only : c_Win, c_Info
    use :: mpi_c_interface, only : MPIR_Win_get_info_c

    implicit none

    type(MPI_Win), intent(in) :: win
    type(MPI_Info), intent(out) :: info_used
    integer, optional, intent(out) :: ierror

    integer(c_Win) :: win_c
    integer(c_Info) :: info_used_c
    integer(c_int) :: ierror_c

    if (c_int == kind(0)) then
        ierror_c = MPIR_Win_get_info_c(win%MPI_VAL, info_used%MPI_VAL)
    else
        win_c = win%MPI_VAL
        ierror_c = MPIR_Win_get_info_c(win_c, info_used_c)
        info_used%MPI_VAL = info_used_c
    end if

    if (present(ierror)) ierror = ierror_c

end subroutine MPI_Win_get_info_f08
