/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2007 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "mpiimpl.h"

#define PMI2U_Malloc MPIU_Malloc
#define PMI2U_Free MPIU_Free
#define PMI2U_Strdup MPIU_Strdup
#define PMI2U_Snprintf MPL_snprintf
#define PMI2U_Strncpy MPIU_Strncpy
#define PMI2U_Strnapp MPIU_Strnapp
#define PMI2U_Assert MPIU_Assert
#define PMI2U_Exit MPL_exit
#define PMI2U_Info MPID_Info
#define PMI2U_Memcpy MPIU_Memcpy

