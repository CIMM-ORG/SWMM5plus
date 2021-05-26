!   -*- Mode: Fortran; -*-
!
!   (C) 2014 by Argonne National Laboratory.
!   See COPYRIGHT in top-level directory.
!
subroutine MPI_Dist_graph_neighbors_f08(comm, maxindegree, sources, sourceweights, &
    maxoutdegree, destinations, destweights, ierror)
    use, intrinsic :: iso_c_binding, only : c_int
    use :: mpi_f08, only : MPI_Comm
    use :: mpi_c_interface, only : c_Comm
    use :: mpi_c_interface, only : MPIR_Dist_graph_neighbors_c

    implicit none

    type(MPI_Comm), intent(in) :: comm
    integer, intent(in) :: maxindegree
    integer, intent(in) :: maxoutdegree
    integer, intent(out) :: sources(maxindegree)
    integer, intent(out) :: destinations(maxoutdegree)
    integer, intent(out) :: sourceweights(maxindegree)
    integer, intent(out) :: destweights(maxoutdegree)
    integer, optional, intent(out) :: ierror

    integer(c_Comm) :: comm_c
    integer(c_int) :: maxindegree_c
    integer(c_int) :: maxoutdegree_c
    integer(c_int) :: sources_c(maxindegree)
    integer(c_int) :: destinations_c(maxoutdegree)
    integer(c_int) :: sourceweights_c(maxindegree)
    integer(c_int) :: destweights_c(maxoutdegree)
    integer(c_int) :: ierror_c

    if (c_int == kind(0)) then
        ierror_c = MPIR_Dist_graph_neighbors_c(comm%MPI_VAL, maxindegree, sources, sourceweights, maxoutdegree, &
            destinations, destweights)
    else
        comm_c = comm%MPI_VAL
        maxindegree_c = maxindegree
        maxoutdegree_c = maxoutdegree
        ierror_c = MPIR_Dist_graph_neighbors_c(comm_c, maxindegree_c, sources_c, sourceweights_c, maxoutdegree_c, &
            destinations_c, destweights_c)
        sources = sources_c
        sourceweights = sourceweights_c
        destinations = destinations_c
        destweights = destweights_c
    end if

    if (present(ierror)) ierror = ierror_c

end subroutine MPI_Dist_graph_neighbors_f08
