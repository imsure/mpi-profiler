
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "vector.h"
#include "helpers.h"
#include "critical_path.h"

#define DEBUG 1

Graph G;
Index2Name *i2n;

static int GetNodeIndex( char *name );
static void toposort();

static void print_graph()
{
  int i, j;
  printf( "Printing graph......\n" );
  printf( "Number of nodes: %d\n", G.NumNodes );

  for ( i = 0; i < G.NumNodes; ++i ) {
    AdjList *adj = G.adjlist[ i ];
    printf( "%s ==> ", adj->src );
    for ( j = 0; j < adj->NumEdges; ++j ) {
      Edge e = adj->edges[ j ];
      printf( "(%s,%d,%d) ", e.dest, e.time, e.msgsize );
    }
    putchar( '\n' );
  }
}

static void count_nodes( vertex **graphs, int *sizes,int num, int *node_count )
{
  int i, j, collective_count = 0;

  /* count number of nodes in the graph */
  *node_count += 2; // 2 counts for MPI_Init and MPI_Finalize
  for ( i = 1; i < sizes[0]-1; ++i ) {
    if ( !is_sendrecv_oper( graphs[0][i].name ) ) {
      ++collective_count; // count number of collectives
    }
    ++(*node_count);
  }

  for ( i = 1; i < num; ++i ) {
    // Ignore MPI_Init and MPI_Finalize and all collectives
    // since they are all shared by ranks.
    *node_count += sizes[i] - collective_count - 2;
  }
}

static void MatchSendRecv( vertex **graphs, int *vec_sizes, int numranks )
{
  int i, j, msg_size;
  vertex *g;
  char *name;
  int *indexes = (int *) malloc( numranks * sizeof(int) );
  AdjList *adj;
  
  for ( i = 0; i < numranks; ++i ) {
    indexes[ i ] = 0; // always traverse from the beginning.
  }

  for ( i = 0; i < numranks; ++i ) {
    g = graphs[ i ];
    for ( j = 0; j < vec_sizes[i]-1; ++j ) {
      name = g[ j ].name;
      if ( is_send_oper(name) ) {
	vertex *g2;
	int k, rank2, latency;

	msg_size = g[ j ].msg_size;
	latency = sendrecv_latency( msg_size );
	rank2 = g[ j ].receiver_rank;
	g2 = graphs[ rank2 ]; // get the graph of the receiver rank.

	// Go through the graph of receiver rank, find the match receving
	// operations, either MPI_Recv or MPI_Wait.
	for ( k = indexes[rank2]; k < vec_sizes[rank2]; ++k ) {
	  if ( is_block_recv_oper( g2[k].name ) ) {
	    if ( g[j].sender_rank == g2[k].sender_rank &&
		 g[j].receiver_rank == g2[k].receiver_rank &&
		 g[j].tag == g2[k].tag ) { // find a match
	      
	      adj = (AdjList *) malloc( sizeof(AdjList) );
	      strcpy( adj->src, name ); // sender name
	      adj->NumEdges = 2; // there is always two edges for a sender
	      adj->edges = (Edge *) malloc( adj->NumEdges * sizeof(Edge) );

	      // edge in the same rank: computation edge
	      strcpy( adj->edges[0].dest, graphs[i][j+1].name );
	      adj->edges[0].time = (int ) round( (graphs[i][j+1].start_time -
						  graphs[i][j].end_time) * 1000 );

	      // edge in the same rank: latency edge
	      strcpy( adj->edges[1].dest, g2[k].name );
	      adj->edges[1].time = latency;
	      adj->edges[1].msgsize = msg_size;

	      G.adjlist[ GetNodeIndex(name) ] = adj;

	      indexes[ rank2 ]++;
	      break;
	    }
	  }
	  indexes[ rank2 ]++;
	}
      }
    }
    // Reset indexes for the next graph.
    for ( j = 0; j < numranks; ++j ) {
      indexes[ j ] = 0; 
    }
  }
}

static void BuildAdjList( vertex **graphs, int *sizes, int numranks )
{
  int i, j;
  AdjList *adj = (AdjList *) malloc( sizeof(AdjList) );

  // build adj list for MPI_Init
  if ( is_collective_oper(graphs[0][1].name) ||
       strcmp(graphs[0][1].name, "MPI_Finalize") == 0 ) {
    if ( G.adjlist[ GetNodeIndex("MPI_Init") ] == NULL ) {
      adj = (AdjList *) malloc( sizeof(AdjList) );
      adj->NumEdges = 1;
      adj->edges = (Edge *) malloc( 1 * sizeof(Edge) );
      strcpy( adj->src, "MPI_Init" );
      strcpy( adj->edges[0].dest, graphs[0][1].name );
      G.adjlist[ GetNodeIndex("MPI_Init") ] = adj;
      G.adjlist[GetNodeIndex("MPI_Init")]->edges[0].time = 0;
    }
    for ( i = 0; i < numranks; ++i ) {
      int time = (int) round( (graphs[i][1].start_time -
			       graphs[i][0].end_time) * 1000 );
      if ( time > G.adjlist[GetNodeIndex("MPI_Init")]->edges[0].time )
	G.adjlist[GetNodeIndex("MPI_Init")]->edges[0].time = time;
    }
  } else {
    strcpy( adj->src, "MPI_Init" );
    adj->edges = (Edge *) malloc( numranks * sizeof(Edge) );
    for ( i = 0; i < numranks; ++i ) {
      strcpy( adj->edges[i].dest, graphs[i][1].name );
      adj->edges[i].time = (int) round( (graphs[i][1].start_time -
					 graphs[i][0].end_time) * 1000 );
    }
    adj->NumEdges = numranks;
    G.adjlist[ GetNodeIndex("MPI_Init") ] = adj;
  }

  MatchSendRecv( graphs, sizes, numranks );
  for ( i = 0; i < numranks; ++i ) {
    // skip MPI_Init and MPI_Finalize
    for ( j = 1; j < sizes[i]-1; ++j ) {
      if ( is_recv_oper(graphs[i][j].name) ) { // recvs
	adj = (AdjList *) malloc( sizeof(AdjList) );
	strcpy( adj->src, graphs[i][j].name );
	adj->NumEdges = 1; // always one outgoing edge for recvs
	adj->edges = (Edge *) malloc( adj->NumEdges * sizeof(Edge) );
	strcpy( adj->edges[0].dest, graphs[i][j+1].name );
	adj->edges[0].time = (int) round( (graphs[i][j+1].start_time -
					    graphs[i][j].end_time) * 1000 );
	G.adjlist[ GetNodeIndex(graphs[i][j].name) ] = adj;
      } else if ( is_collective_oper(graphs[i][j].name) ) { // collectives
	
	if ( G.adjlist[ GetNodeIndex(graphs[i][j].name) ] == NULL ) {
	  adj = (AdjList *) malloc( sizeof(AdjList) );
	  strcpy( adj->src, graphs[i][j].name );
	  if ( is_collective_oper(graphs[i][j+1].name) ||
	       strcmp(graphs[i][j+1].name, "MPI_Finalize") == 0 ) {
	    adj->NumEdges = 1;
	  } else {
	    adj->NumEdges = numranks; // 'numranks' fanout for collectives
	  }
	  adj->edges = (Edge *) malloc( adj->NumEdges * sizeof(Edge) );
	  strcpy( adj->edges[0].dest, graphs[i][j+1].name );
	  G.adjlist[ GetNodeIndex(graphs[i][j].name) ] = adj;
	  G.adjlist[GetNodeIndex(graphs[i][j].name)]->edges[0].time = 0;
	}
	int time = (int) round( (graphs[i][j+1].start_time -
				 graphs[i][j].end_time) * 1000 );

	if ( is_collective_oper(graphs[i][j+1].name) ||
	     strcmp(graphs[i][j+1].name, "MPI_Finalize") == 0 ) {
	  if ( time > G.adjlist[GetNodeIndex(graphs[i][j].name)]->edges[0].time )
	    G.adjlist[GetNodeIndex(graphs[i][j].name)]->edges[0].time = time;
	} else {
	  strcpy( G.adjlist[GetNodeIndex(graphs[i][j].name)]->edges[i].dest, graphs[i][j+1].name );
	  G.adjlist[GetNodeIndex(graphs[i][j].name)]->edges[i].time = time;
	}
      }
    }
  }

  adj = (AdjList *) malloc( sizeof(AdjList) );
  strcpy( adj->src, "MPI_Finalize" );
  adj->NumEdges = 0;
  G.adjlist[ GetNodeIndex("MPI_Finalize") ] = adj;

  print_graph();
}

static void SetupIndex2Name( vertex **graphs, int *sizes, int numranks )
{
  int i, j, collective_count = 0, count = 0;

  i2n[ count ].index = count;
  strcpy( i2n[ count ].name, "MPI_Init" );
  count++;
  i2n[ count ].index = count;
  strcpy( i2n[ count ].name, "MPI_Finalize" );
  count++;
  
  for ( i = 1; i < sizes[0]-1; ++i ) {
    if ( !is_sendrecv_oper( graphs[0][i].name ) ) {
      i2n[ count ].index = count;
      strcpy( i2n[ count ].name, graphs[0][i].name );
      count++;
    }
  }

  for ( i = 0; i < numranks; ++i ) {
    // Ignore MPI_Init and MPI_Finalize and all collectives
    // since they are all shared by ranks.
    for ( j = 1; j < sizes[i]-1; ++j ) {
      if ( is_sendrecv_oper( graphs[i][j].name ) ) {
	i2n[ count ].index = count;
	strcpy( i2n[ count ].name, graphs[i][j].name );
	count++;
      }
    }
  }

#if DEBUG
  printf( "Number of nodes: %d\n", count );
  printf( "Node name : index\n" );
  for ( i = 0; i < count; ++i ) {
    printf( "%s : %d\n", i2n[i].name, i2n[i].index );
  }
#endif
}

static int GetNodeIndex( char *name )
{
  int i;
  for ( i = 0; i < G.NumNodes; ++i ) {
    if ( strcmp(i2n[i].name, name) == 0)
      return i2n[i].index;
  }
}

/**
 * graphs: an array of graph defined in profile.c
 * sizes: number of nodes in each graph
 * num: size of the array 'graphs' and 'sizes'
 */
static void build_graph( vertex **graphs, int *sizes, int numranks )
{
  int node_count = 0, i, j, id = 0;
  int *indexes;
  char *name;
  vertex *g;

  count_nodes( graphs, sizes, numranks, &node_count );
  G.NumNodes = node_count;
  i2n = (Index2Name *) malloc( node_count * sizeof(Index2Name) );
  SetupIndex2Name( graphs, sizes, numranks );

  G.adjlist = (AdjList **) malloc( G.NumNodes * sizeof(AdjList*) );
  for ( i = 0; i < G.NumNodes; ++i )
    G.adjlist[ i ] = NULL;
  BuildAdjList( graphs, sizes, numranks );
}

void find_critical_path( vertex **graphs, int *sizes, int numranks )
{
  int i, j, prev, len;
  FILE *cpath;
  
  build_graph( graphs, sizes, numranks );
}
