#ifndef _CRITICAL_PATH_H
#define _CRITICAL_PATH_H

typedef struct Edge Edge;
typedef struct Graph Graph;
typedef struct AdjList AdjList;
typedef struct Index2Name Index2Name;

struct Index2Name {
  int index;
  char name[ 32 ];
};

struct AdjList {
  char src[ 32 ];
  Edge *edges;
  int NumEdges;
};

struct Edge {
  char dest[ 32 ];
  int time; // in ms
  int msgsize;
};

struct Graph {
  int NumNodes;
  AdjList **adjlist;
};

void find_critical_path( vertex **graphs, int *sizes, int numranks );

#endif
