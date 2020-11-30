#include "utils.h"

// ----------------------------
// Data Structures
// ----------------------------

// ** Stacks

void s_append(Stack* s, int val)
{
  struct Snode* n = malloc(sizeof(struct Snode));
  n->val = val;
  n->next = s->last;
  s->last = n;
  s->len ++;
}

int s_pop(Stack* s)
{
  int val;
  struct Snode* tmp;

  if (s->len > 0)
  {
    val = s->last->val;
    tmp = s->last;
    s->last = s->last->next;
    s->len --;
    free(tmp);
    return val;
  }
  else
  {
    printf("Error, the stack is empty\n");
    return -1;
  }
}

Stack* create_stack()
{
  Stack* s = malloc(sizeof(Stack));
  s->len = 0;
  s->last = NULL;
  return s;
}

// ** Graphs

G_node* createNode(int v) {
  G_node* newNode = malloc(sizeof(G_node));
  newNode->vertex = v;
  newNode->next = NULL;
  return newNode;
}

// Create graph
Graph* createGraph(int vertices) {
  int i;
  Graph* graph = malloc(sizeof(Graph));
  graph->numVertices = vertices;
  graph->adjLists = malloc(vertices * sizeof(G_node*));
  graph->visited = malloc(vertices * sizeof(int));

  for (i = 0; i < vertices; i++) {
    graph->adjLists[i] = NULL;
    graph->visited[i] = 0;
  }
  return graph;
}

void addEdge(struct Graph* graph, int src, int dest) {
  struct node* newNode = createNode(dest);
  newNode->next = graph->adjLists[src];
  graph->adjLists[src] = newNode;
}