/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2006 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

int main(int argc, char **argv)
{
    int *ptr = NULL;
    int rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        (*ptr)++;       /* will segfault */
        /* never get below here, just present to prevent dead-code elimination */
        printf("*ptr=%d\n", (*ptr));
    }
    MPI_Finalize();
    return 0;
}
