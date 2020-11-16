#include <stdlib.h>
#include <stdio.h>

#define INITIAL_ARRAY_SIZE 50

typedef struct {
  double ** table;
  int used;
  int size;
  int dim;
} API_Table;

void init_api_table(API_Table *a, int dimension);
void append_api_table(API_Table *a, size_t size, ...);
void trim_api_table(API_Table *a);
void free_api_table(API_Table *a);