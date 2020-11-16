#include <stdarg.h>

#include "arrays.h"

// --- Array ---
// --- Table ---

void init_api_table(API_Table *a, int dimension)
{
    int i;
    a->table = (double **) malloc(dimension*sizeof(double*));
    for (int i = 0; i < dimension; i++)
    {
        a->table[i] = (double *) malloc(INITIAL_ARRAY_SIZE*sizeof(double));
    }
    a->used = 0;
    a->size = INITIAL_ARRAY_SIZE;
    a->dim = dimension;
}

void append_api_table(API_Table *a, size_t size, ...)
{
    int i;
    double x;
    va_list valist;

    va_start(valist, size); // size == a->dim

    if (a->used == a->size) {
        a->size *= 2;
        for (i = 0; i < a->dim; i++)
        {
            a->table[i] = (double *) realloc(a->table, a->size * sizeof(double));
        }
    }
    for (i = 0; i < a->dim; i++)
    {
        x = va_arg(valist, double);
        a->table[i][a->used++] = x;
    }

   va_end(valist);
}

void trim_api_table(API_Table *a)
{
    int i;
    a->size = a->used;
    for (i = 0; i < a->dim; i++)
    {
        a->table[i] = (double *) realloc(a->table, a->size * sizeof(double));
    }
}

void free_api_table(API_Table *a)
{
    int i;
    for (i = 0; i < a->dim; i++)
    {
        free(a->table[i]);
    }
    free(a->table);
    a->table = NULL;
    a->used = a->size = 0;
}