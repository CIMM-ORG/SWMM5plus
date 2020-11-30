#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "headers.h"

//-----------------------------------------------------------------------------
//  Shared variables
//-----------------------------------------------------------------------------
static char *Tok[MAXTOKS];             // String tokens from line of input
static int  Ntokens;                   // Number of tokens in line of input

struct Snode {
    int val;
    struct Snode* next;
};

typedef struct Stacks {
    int len;
    struct Snode* last;
} Stack;

typedef struct g_node {
  int vertex;
  G_node * next;
} G_node;

typedef struct graph {
  int num_vertices;
  int* visited;
  G_node** adjLists;
} Graph;


void get_num_hydrographs(void* f_api);
int getTokens(char *s);
void printInt(void *n);
void printFloat(void *f);
void printList(struct UNode *node, void (*fptr)(void *));
void push(struct UNode** head_ref, void *new_data, size_t data_size);

// Graph
G_node* createNode(int v);
Graph* createGraph(int vertices);
void addEdge(struct Graph* graph, int src, int dest);