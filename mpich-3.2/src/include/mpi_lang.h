/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 *
 * Portions of this code were written by Microsoft. Those portions are
 * Copyright (c) 2007 Microsoft Corporation. Microsoft grants
 * permission to use, reproduce, prepare derivative works, and to
 * redistribute to others. The code is licensed "as is." The User
 * bears the risk of using it. Microsoft gives no express warranties,
 * guarantees or conditions. To the extent permitted by law, Microsoft
 * excludes the implied warranties of merchantability, fitness for a
 * particular purpose and non-infringement.
 */

#ifndef MPI_LANG_H_INCLUDED
#define MPI_LANG_H_INCLUDED

#include "mpi_attr.h"

#if defined(__cplusplus)
extern "C" {
#endif

/*E
  Language bindings for MPI

  A few operations in MPI need to know how to marshal the callback into the calling
  lanuage calling convention. The marshaling code is provided by a thunk layer which
  implements the correct behavior.  Examples of these callback functions are the
  keyval attribute copy and delete functions.

  Module:
  Attribute-DS
  E*/

/*
 * Support bindings for Attribute copy/del callbacks
 * Consolidate Comm/Type/Win attribute functions together as the handle type is the same
 * use MPI_Comm for the prototypes
 */
typedef
int
(MPID_Attr_copy_proxy)(
    MPI_Comm_copy_attr_function* user_function,
    int handle,
    int keyval,
    void* extra_state,
    MPIR_AttrType attrib_type,
    void* attrib,
    void** attrib_copy,
    int* flag
    );

typedef
int
(MPID_Attr_delete_proxy)(
    MPI_Comm_delete_attr_function* user_function,
    int handle,
    int keyval,
    MPIR_AttrType attrib_type,
    void* attrib,
    void* extra_state
    );

void
MPIR_Keyval_set_proxy(
    int keyval,
    MPID_Attr_copy_proxy copy_proxy,
    MPID_Attr_delete_proxy delete_proxy
    );

#if defined(__cplusplus)
}
#endif

#endif /* MPI_LANG_H_INCLUDED */
