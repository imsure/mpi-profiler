
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "vector.h"
#include "helpers.h"
#include "critical_path.h"

#define DEBUG 0

longest *lon;
node *path;
int path_len;
Graph G;

static void toposort();

static void print_graph()
{
  int i,j;
  for ( i = 0; i < G.node_count; ++i ) {
    printf( "Node %s\n", G.nodes[i].name );
  }
  for ( i = 0; i < G.edge_count; ++i ) {
    printf( "Edge %s -> %s (%d removed=%d ignore=%d)\n", G.edges[i].from.name,
	    G.edges[i].to.name, G.edges[i].weight, G.edges[i].removed, G.edges[i].ignore );
  }
}

static void count_nodes_edges( vertex **graphs, int *sizes,int num,
			       int *node_count, int *edge_count )
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

  /* count number of edges in the graph */
  for ( i = 0; i < num; ++i ) {
    // count the number of edges within each local graph
    *edge_count += sizes[i] - 1;
    for ( j = 1; j < sizes[i] - 1; ++j ) {
      if ( is_send_oper( graphs[i][j].name ) ) {
	++(*edge_count); // count 1 if it is MPI_Send or MPI_Isend
      }
    }
  }
}

static void remove_redundent_edges( int numranks )
{
  int i, j, max, maxindex;
  node from, to;

  for ( i = 1; i < G.edge_count; ++i ) {
    from = G.edges[i].from;
    to = G.edges[i].to;

    if ( is_collective_oper( from.name ) &&
	 is_collective_oper( to.name ) ) {
      max = G.edges[i].weight;
      maxindex = i;

      for ( j = i+1; j < G.edge_count; ++j ) {
	if ( strcmp(from.name, G.edges[j].from.name) == 0 &&
	     strcmp(to.name, G.edges[j].to.name) == 0 ) {
	  if ( G.edges[j].weight > max ) {
	    max = G.edges[j].weight;
	    maxindex = j;
	  }
	}
      }

      for ( j = i; j < G.edge_count; ++j ) {
	if ( strcmp(from.name, G.edges[j].from.name) == 0 &&
	     strcmp(to.name, G.edges[j].to.name) == 0 ) {
	  if ( j != maxindex ) {
	    G.edges[j].ignore = 1;
	    G.edges[j].removed = 1;
	  }
	}
      }      
    }
  }

#if DEBUG
  print_graph();
#endif
}

/**
 * graphs: an array of graph defined in profile.c
 * sizes: number of nodes in each graph
 * num: size of the array 'graphs' and 'sizes'
 */
static void build_graph( vertex **graphs, int *sizes, int num )
{
  int node_count = 0, edge_count = 0, i, j, id = 0;
  int *indexes;
  char *name;
  vertex *g;

  count_nodes_edges( graphs, sizes, num, &node_count, &edge_count );
#if DEBUG
  printf( "Number of nodes: %d\n", node_count );
  printf( "Number of edges: %d\n", edge_count );
#endif

  G.nodes = (node *) malloc( node_count * sizeof(node) );
  G.edges = (edge *) malloc( edge_count * sizeof(edge) );
  node_count = 0; // reset
  edge_count = 0; // reset
  
  /* build nodes */
  
  for ( i = 0; i < sizes[0]; ++i ) {
    strcpy( G.nodes[ node_count++ ].name, graphs[0][i].name );
  }
      
  for ( i = 1; i < num; ++i ) {
    for ( j = 1; j < sizes[i] - 1; ++j ) {
      if ( is_sendrecv_oper( graphs[i][j].name ) ) {
	strcpy( G.nodes[ node_count++ ].name, graphs[i][j].name );
      }
    }
  }

  /* build edges */
  for ( i = 0; i < num; ++i ) {
    for ( j = 0; j < sizes[i]-1; ++j ) {
      strcpy( G.edges[ edge_count ].from.name, graphs[i][j].name );
      strcpy( G.edges[ edge_count ].to.name, graphs[i][j+1].name );
      G.edges[ edge_count ].weight = (int) round( graphs[i][j+1].start_time -
						  graphs[i][j].end_time ) * 1000; // ms
      G.edges[ edge_count ].removed = 0;
      G.edges[ edge_count ].ignore = 0;
      ++edge_count;
    }
  }

  /* Third pass: collect edges between local graphs. */

  indexes = (int *) malloc( num * sizeof(int) );
  for ( i = 0; i < num; ++i ) {
    indexes[ i ] = 0; // always traverse from the beginning.
  }

  for ( i = 0; i < num; ++i ) {
    g = graphs[ i ];
    for ( j = 0; j < sizes[i]; ++j ) {
      name = g[ j ].name;
      if ( is_send_oper(name) ) {
	vertex *g2;
	int k, rank2, latency;

	latency = sendrecv_latency( g[j].msg_size );
	rank2 = g[ j ].receiver_rank;
	g2 = graphs[ rank2 ]; // get the graph of the receiver rank.

	// Go through the graph of receiver rank, find the match receving
	// operations, either MPI_Recv or MPI_Wait.
	for ( k = indexes[rank2]; k < sizes[rank2]; ++k ) {
	  if ( is_block_recv_oper( g2[k].name ) ) {
	    if ( g[j].sender_rank == g2[k].sender_rank &&
		 g[j].receiver_rank == g2[k].receiver_rank &&
		 g[j].tag == g2[k].tag ) { // find a match

	      strcpy( G.edges[ edge_count ].from.name, g[j].name );
	      strcpy( G.edges[ edge_count ].to.name, g2[k].name );
	      G.edges[ edge_count ].weight = latency;
	      G.edges[ edge_count ].removed = 0;
	      G.edges[ edge_count ].ignore = 0;
	      G.edges[ edge_count ].msgsize = g[j].msg_size;

	      ++edge_count;

	      indexes[ rank2 ]++;
	      break;
	    }
	  }
	  indexes[ rank2 ]++;
	}
      }
    }
    // Reset indexes for the next graph.
    for ( j = 0; j < num; ++j ) {
      indexes[ j ] = 0; 
    }
  }

  G.node_count = node_count;
  G.edge_count = edge_count;
  
#if DEBUG
  //print_graph();
#endif
}

static node * get_neighbors( node *nd, int *size )
{
  int i;
  node *ns = (node *) malloc( 25 * sizeof(node) );
  *size = 0;

  for ( i = 0; i < G.edge_count; ++i ) {
    if ( strcmp(nd->name, G.edges[i].from.name) == 0 ) {
      if ( G.edges[i].ignore == 0 ) { // skip ignored edges.
	ns[ (*size)++ ] = G.edges[ i ].to;
	G.edges[ i ].removed = 1; // remove the edge
      }
    }
  }

  return ns;
}

static int has_incoming_edges( node *nd )
{
  int i;
  for ( i = 0; i < G.edge_count; ++i ) {
    if ( strcmp(nd->name, G.edges[i].to.name) == 0  &&
	 G.edges[i].removed == 0 ) {
      return 1;
    }
  }
  return 0;
}

/**
 * Topologically sorting the graph G.
 */
static void toposort()
{
  int i, j, tmp_count, topo_count = 0, size;
  node *neighbors;
  
  // list of nodes in topological order
  node *topo_nodes = (node *) malloc( G.node_count * sizeof(node) );
  // list of nodes with no incoming edges
  node *tmp_nodes = (node *) malloc( G.node_count * sizeof(node) );

  /* Initially, tmp_nodes has only MPI_Init. */
  tmp_nodes[0] = G.nodes[0];
  tmp_count = 1;

  while ( tmp_count > 0 ) {
    topo_nodes[ topo_count ] = tmp_nodes[ --tmp_count ];
    neighbors = get_neighbors( &topo_nodes[topo_count], &size );
    for ( i = 0; i < size; ++i ) {
      //printf( "neighbor %d: %s\n", i, neighbors[i].name );
      if ( !has_incoming_edges( &neighbors[i] ) ) {
	tmp_nodes[ tmp_count++ ] = neighbors[i];
      }
    }
    topo_count++;
  }

#if DEBUG
  printf( "topo_count = %d\n", topo_count );
  for ( i = 0; i < topo_count; ++i ) {
    printf( "%s -> ", topo_nodes[i].name );
  }
  putchar( '\n' );
#endif

  G.nodes = topo_nodes; // replace G.nodes with topologically ording nodes.
  /* restore removed flag for edges */
  for ( i = 0; i < G.edge_count; ++i ) {
    if ( G.edges[i].ignore == 0 ) {
      G.edges[i].removed = 0;
    }
  }
}

static int get_record( node *nd )
{
  int i;
  for ( i = 0; i < G.node_count; ++i) {
    if ( strcmp(nd->name, lon[i].nd.name) == 0 ) {
      return lon[i].len;
    }
  }
}

static int get_incomings( node *nd )
{
  int i, max = 0, len, prev;
  node *ns = (node *) malloc( 25 * sizeof(node) );
  
  for ( i = 0; i < G.edge_count; ++i ) {
    if ( strcmp(nd->name, G.edges[i].to.name) == 0  &&
	 G.edges[i].ignore == 0 ) {
      prev = get_record( &G.edges[i].from );
      len = prev + G.edges[i].weight;
      if ( len > max )
	max = len;
    }
  }

  return max;
}

void find_critical_path( vertex **graphs, int *sizes, int numranks )
{
  int i, j, prev, len;
  FILE *cpath;
  
  build_graph( graphs, sizes, numranks );
  remove_redundent_edges( numranks );  
  toposort();

  lon = (longest *) malloc( G.node_count * sizeof(longest) );
  path = (node *) malloc( G.node_count * sizeof(node ) );
  path_len = 0;
  
  //print_graph();
  lon[ 0 ].nd = G.nodes[ 0 ];
  lon[ 0 ].len = 0;
  for ( i = 1; i < G.node_count; ++i ) {
    lon[ i ].nd = G.nodes[ i ];
    lon[ i ].len = get_incomings( &G.nodes[i] );
  }

  for ( i = 0; i < G.node_count; ++i ) {
    //printf( "node %s longest %d\n", lon[i].nd.name, lon[i].len );
  }

  path[path_len] = lon[G.node_count-1].nd;

  while ( strcmp(path[path_len].name, "MPI_Init") != 0 ) {
    for ( j = 0; j < G.edge_count; ++j ) {
      if ( strcmp(path[path_len].name, G.edges[j].to.name) == 0  &&
	   G.edges[j].ignore == 0 ) {
	prev = get_record( &G.edges[j].from );
	len = prev + G.edges[j].weight;
	if ( len ==  get_record( &(path[path_len]) )) {
	  path[++path_len] = G.edges[j].from;
	}
      }
    }
  }

  for ( i = 0; i <= path_len; ++i ) {
    //printf( "path: %s\n", path[i].name );
  }

  cpath = fopen( "critPath.out", "w" );
  for ( i = path_len; i >= 1; --i ) {
    if (strcmp(path[i].name, "MPI_Init")) {
      char *name = path[i].name;
      char mpi[30];
      char rank[5];
      int cnt = 0;
      for (j = 0; j < strlen(name); ++j) {
	if (name[j] != '_' && cnt < 2) {
	  mpi[j] = name[j];
	} else if (name[j] == '_' ) {
	  cnt++;
	  if (cnt < 2)
	    mpi[j] = name[j];
	  else {
	    mpi[j] = 0;
	    break;
	  }
	}
      }
      if (is_sendrecv_oper(path[i].name)) {
	int index = 0;
	for ( j = j+1; j < strlen(name); ++j) {
	  if (name[j] != '_') {
	    rank[index++] = name[j];
	  } else break;
	}
	rank[ index ] = 0;
	fprintf( cpath, "%s %s\n", mpi, rank );
      } else {
	fprintf( cpath, "%s -1\n", mpi );
      }
    }
    else {
      fprintf( cpath, "%s -1\n", path[i].name );
    }
    for ( j = 0; j < G.edge_count; ++j ) {
      if ( strcmp(path[i].name, G.edges[j].from.name) == 0 &&
	   strcmp(path[i-1].name, G.edges[j].to.name) == 0 &&
	   G.edges[j].ignore == 0 ) {
	if (is_send_oper(path[i].name)) {
	  fprintf( cpath, "%d\n", G.edges[j].msgsize );
	}
	else {
	  fprintf( cpath, "%d\n", G.edges[j].weight/1000 );
	}
	break;
      }
    }
  }
  fprintf( cpath, "%s -1\n", path[0].name );
}
