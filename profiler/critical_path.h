#ifndef _CRITICAL_PATH_H
#define _CRITICAL_PATH_H

typedef struct node node;
typedef struct edge edge;
typedef struct Graph Graph;
typedef struct longest longest;

struct node {
  char name[32];
};

struct edge {
  node from;
  node to;
  int weight; // in ms
  int removed; // indicate whether it has been removed
  int ignore;
  int msgsize;
};

/* G = (V, E) */
struct Graph {
  node *nodes;
  edge *edges;
  int edge_count;
  int node_count;
};

struct longest {
  node nd;
  int len;
};

void find_critical_path( vertex **graphs, int *sizes, int numranks );

#endif
