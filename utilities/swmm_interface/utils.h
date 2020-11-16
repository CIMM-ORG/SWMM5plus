#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "headers.h"

//-----------------------------------------------------------------------------
//  Shared variables
//-----------------------------------------------------------------------------
static char *Tok[MAXTOKS];             // String tokens from line of input
static int  Ntokens;                   // Number of tokens in line of input

struct UNode
{
    // Any data type can be stored in this node
    void  *data;
    struct UNode *next;
};

void get_num_hydrographs(void* f_api);
int getTokens(char *s);
void printInt(void *n);
void printFloat(void *f);
void printList(struct UNode *node, void (*fptr)(void *));
void push(struct UNode** head_ref, void *new_data, size_t data_size);