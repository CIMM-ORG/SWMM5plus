/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef MPLTRMEM_H_INCLUDED
#define MPLTRMEM_H_INCLUDED

/* FIXME: Consider an option of specifying __attribute__((malloc)) for
   gcc - this lets gcc-style compilers know that the returned pointer
   does not alias any pointer prior to the call.
 */
void MPL_trinit(int);
void *MPL_trmalloc(size_t, int, const char[]);
void MPL_trfree(void *, int, const char[]);
int MPL_trvalid(const char[]);
int MPL_trvalid2(const char[],int,const char[]);
void MPL_trspace(size_t *, size_t *);
void MPL_trid(int);
void MPL_trlevel(int);
void MPL_trDebugLevel(int);
void *MPL_trcalloc(size_t, size_t, int, const char[]);
void *MPL_trrealloc(void *, size_t, int, const char[]);
void *MPL_trstrdup(const char *, int, const char[]);
void MPL_TrSetMaxMem(size_t);

/* Make sure that FILE is defined */
#include <stdio.h>
void MPL_trdump(FILE *, int);
void MPL_trSummary(FILE *, int);
void MPL_trdumpGrouped(FILE *, int);

#endif /* !defined(MPLTRMEM_H_INCLUDED) */
