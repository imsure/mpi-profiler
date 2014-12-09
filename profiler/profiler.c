/**
 * A simple MPI Profiler.
 *
 * Author: Shuo Yang
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "vector.h"
#include "helpers.h"
#include "critical_path.h"

#ifndef _EXTERN_C_
#ifdef __cplusplus
#define _EXTERN_C_ extern "C"
#else /* __cplusplus */
#define _EXTERN_C_
#endif /* __cplusplus */
#endif /* _EXTERN_C_ */

#ifdef MPICH_HAS_C2F
_EXTERN_C_ void *MPIR_ToPointer(int);
#endif // MPICH_HAS_C2F

#ifdef PIC
/* For shared libraries, declare these weak and figure out which one was linked
   based on which init wrapper was called.  See mpi_init wrappers.  */
#pragma weak pmpi_init
#pragma weak PMPI_INIT
#pragma weak pmpi_init_
#pragma weak pmpi_init__
#endif /* PIC */

#define DEBUG 0

#define MPI_Init_Index 0
#define MPI_Send_Index 1

#define NO_RANK -1 // rank for collective operations

_EXTERN_C_ void pmpi_init(MPI_Fint *ierr);
_EXTERN_C_ void PMPI_INIT(MPI_Fint *ierr);
_EXTERN_C_ void pmpi_init_(MPI_Fint *ierr);
_EXTERN_C_ void pmpi_init__(MPI_Fint *ierr);

static int myrank, numranks;
static int sendrecv_id = 0; // ID for send/recv operations
static int collective_id = 0; // ID for collective operations
static vector graph; // vector holding the MPI task graph for every rank
static int *vec_sizes; // array to hold sizes of vectors of each rank
static vertex **graphs; // array of vertex *, each of which point to a local graph
extern node *path;
extern int path_len;

double init_times[200]; 
int init_count = 0;

double send_times[200]; 
int send_count = 0;

double recv_times[200]; 
int recv_count = 0;

double wait_times[200]; 
int wait_count = 0;

double barr_times[200]; 
int barr_count = 0;

double Isend_times[200]; 
int Isend_count = 0;

double Irecv_times[200]; 
int Irecv_count = 0;

MPI_Datatype vertex_type;

int get_mpi_datatype_size( MPI_Datatype datatype )
{
  return 0;
}

/* ================== C Wrappers for MPI_Init ================== */
_EXTERN_C_ int PMPI_Init(int *argc, char ***argv);
_EXTERN_C_ int MPI_Init(int *argc, char ***argv) { 
  int _wrap_py_return_val = 0;
  vertex v;
  MPI_Datatype types[7] = { MPI_CHAR, MPI_DOUBLE, MPI_DOUBLE,
			    MPI_INT, MPI_INT, MPI_INT, MPI_INT };
  int blocklens[7] = { 32, 1, 1, 1, 1, 1, 1 };
  MPI_Aint disp[7];
  
  v.start_time = MPI_Wtime();
  _wrap_py_return_val = PMPI_Init(argc, argv);
  v.end_time = MPI_Wtime();

  init_times[init_count++] = (v.end_time - v.start_time) * 1000;

  /*
   * Set up the MPI datatype for vertex structure.
   */
  MPI_Address( &v.name, &disp[0] );
  MPI_Address( &v.start_time, &disp[1] );
  MPI_Address( &v.end_time, &disp[2] );
  MPI_Address( &v.sender_rank, &disp[3] );
  MPI_Address( &v.receiver_rank, &disp[4] );
  MPI_Address( &v.tag, &disp[5] );
  MPI_Address( &v.msg_size, &disp[6] );
  /* make relative */
  disp[ 1 ] -= disp[ 0 ];
  disp[ 2 ] -= disp[ 0 ];
  disp[ 3 ] -= disp[ 0 ];
  disp[ 4 ] -= disp[ 0 ];
  disp[ 5 ] -= disp[ 0 ];
  disp[ 6 ] -= disp[ 0 ];
  disp[ 0 ] = 0;

  MPI_Type_create_struct( 7, blocklens, disp, types, &vertex_type );
  MPI_Type_commit( &vertex_type );

  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
  MPI_Comm_size( MPI_COMM_WORLD, &numranks );

  vec_sizes = (int *) malloc( numranks * sizeof(int) );
  graphs = (vertex **) malloc( numranks * sizeof(vertex *) );

  init_vector( &graph, INIT_CAPACITY );
  sprintf( v.name, "%s", __FUNCTION__ );
  append_vector( &graph, &v );
    
  return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Send ================== */
_EXTERN_C_ int PMPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
_EXTERN_C_ int MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) { 
  int _wrap_py_return_val = 0, size;
  vertex v;

  v.start_time = MPI_Wtime();
  _wrap_py_return_val = PMPI_Send(buf, count, datatype, dest, tag, comm);
  v.end_time = MPI_Wtime();

  send_times[send_count++] = (v.end_time - v.start_time) * 1000;

  v.sender_rank = myrank;
  v.receiver_rank = dest;
  v.tag = tag;
  MPI_Type_size( datatype, &size ); // get size of message in bytes
  v.msg_size = count * size;
  sprintf( v.name, "%s_%d_%d", __FUNCTION__, myrank, sendrecv_id++ );
  
  append_vector( &graph, &v );
  
  return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Recv ================== */
_EXTERN_C_ int PMPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
_EXTERN_C_ int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status) { 
  int _wrap_py_return_val = 0;
  vertex v;

  v.start_time = MPI_Wtime();
  _wrap_py_return_val = PMPI_Recv(buf, count, datatype, source, tag, comm, status);
  v.end_time = MPI_Wtime();

  recv_times[recv_count++] = (v.end_time - v.start_time) * 1000;

  sprintf( v.name, "%s_%d_%d", __FUNCTION__, myrank, sendrecv_id++ );
  v.sender_rank = source;
  v.receiver_rank = myrank;
  v.tag = tag;
  
  append_vector( &graph, &v );

  return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Barrier ================== */
_EXTERN_C_ int PMPI_Barrier(MPI_Comm comm);
_EXTERN_C_ int MPI_Barrier(MPI_Comm comm) { 
  int _wrap_py_return_val = 0;
  vertex v;

  v.start_time = MPI_Wtime();
  _wrap_py_return_val = PMPI_Barrier(comm);
  v.end_time = MPI_Wtime();

  barr_times[barr_count++] = (v.end_time - v.start_time) * 1000;

  sprintf( v.name, "%s_%d", __FUNCTION__, collective_id++ );
  append_vector( &graph, &v );
  
  return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Alltoall ================== */
_EXTERN_C_ int PMPI_Alltoall(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm);
_EXTERN_C_ int MPI_Alltoall(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) { 
  int _wrap_py_return_val = 0;
  vertex v;
  
  v.start_time = MPI_Wtime();
  _wrap_py_return_val = PMPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
  v.end_time = MPI_Wtime();

  sprintf( v.name, "%s_%d", __FUNCTION__, collective_id++ );
  append_vector( &graph, &v );

  return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Scatter ================== */
_EXTERN_C_ int PMPI_Scatter(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
_EXTERN_C_ int MPI_Scatter(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) { 
  int _wrap_py_return_val = 0;
  vertex v;

  v.start_time = MPI_Wtime();
  _wrap_py_return_val = PMPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
  v.end_time = MPI_Wtime();

  sprintf( v.name, "%s_%d", __FUNCTION__, collective_id++ );
  append_vector( &graph, &v );

  return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Gather ================== */
_EXTERN_C_ int PMPI_Gather(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
_EXTERN_C_ int MPI_Gather(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) { 
  int _wrap_py_return_val = 0;
  vertex v;

  v.start_time = MPI_Wtime();
  _wrap_py_return_val = PMPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
  v.end_time = MPI_Wtime();

  sprintf( v.name, "%s_%d", __FUNCTION__, collective_id++ );
  append_vector( &graph, &v );

  return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Reduce ================== */
_EXTERN_C_ int PMPI_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
_EXTERN_C_ int MPI_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) { 
  int _wrap_py_return_val = 0;
  vertex v;
  
  v.start_time = MPI_Wtime();
  _wrap_py_return_val = PMPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
  v.end_time = MPI_Wtime();

  sprintf( v.name, "%s_%d", __FUNCTION__, collective_id++ );
  append_vector( &graph, &v );

  return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Allreduce ================== */
_EXTERN_C_ int PMPI_Allreduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
_EXTERN_C_ int MPI_Allreduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) { 
  int _wrap_py_return_val = 0;
  vertex v;
  
  v.start_time = MPI_Wtime();
  _wrap_py_return_val = PMPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
  v.end_time = MPI_Wtime();

  sprintf( v.name, "%s_%d", __FUNCTION__, collective_id++ );
  append_vector( &graph, &v );

  return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Isend ================== */
_EXTERN_C_ int PMPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
_EXTERN_C_ int MPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request) { 
  int _wrap_py_return_val = 0, size;
  vertex v;

  v.start_time = MPI_Wtime();
  _wrap_py_return_val = PMPI_Isend(buf, count, datatype, dest, tag, comm, request);
  v.end_time = MPI_Wtime();

  Isend_times[Isend_count++] = (v.end_time - v.start_time) * 1000;

  v.sender_rank = myrank;
  v.receiver_rank = dest;
  v.tag = tag;
  MPI_Type_size( datatype, &size ); // get size of message in bytes
  v.msg_size = count * size;
  sprintf( v.name, "%s_%d_%d", __FUNCTION__, myrank, sendrecv_id++ );

  append_vector( &graph, &v );

  return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Irecv ================== */
_EXTERN_C_ int PMPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request);
_EXTERN_C_ int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request) { 
  int _wrap_py_return_val = 0;
  vertex v;
  
  v.start_time = MPI_Wtime();
  _wrap_py_return_val = PMPI_Irecv(buf, count, datatype, source, tag, comm, request);
  v.end_time = MPI_Wtime();

  Irecv_times[Irecv_count++] = (v.end_time - v.start_time) * 1000;

  sprintf( v.name, "%s_%d_%d", __FUNCTION__, myrank, sendrecv_id++ );
  append_vector( &graph, &v );

  return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Wait ================== */
_EXTERN_C_ int PMPI_Wait(MPI_Request *request, MPI_Status *status);
_EXTERN_C_ int MPI_Wait(MPI_Request *request, MPI_Status *status) { 
  int _wrap_py_return_val = 0;
  vertex v;

  v.start_time = MPI_Wtime();
  _wrap_py_return_val = PMPI_Wait(request, status);
  v.end_time = MPI_Wtime();

  wait_times[wait_count++] = (v.end_time - v.start_time) * 1000;

  sprintf( v.name, "%s_%d_%d", __FUNCTION__, myrank, sendrecv_id++ );
  v.sender_rank = status->MPI_SOURCE;
  v.receiver_rank = myrank;
  v.tag = status->MPI_TAG;

  append_vector( &graph, &v );

  return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Waitall ================== */
_EXTERN_C_ int PMPI_Waitall(int count, MPI_Request *array_of_requests, MPI_Status *array_of_statuses);
_EXTERN_C_ int MPI_Waitall(int count, MPI_Request *array_of_requests, MPI_Status *array_of_statuses) { 
  int _wrap_py_return_val = 0, i;

  /* convert waitall into waits. */
  for (i = 0; i < count; ++i) {
    MPI_Wait( &array_of_requests[i], &array_of_statuses[i] );
  }

  return _wrap_py_return_val;
}

static void collect_graphs()
{
  int i,j;

  /* First, send size of vector to rank 0 from other ranks. */
  if (myrank != 0) {
    PMPI_Send( &graph.size, 1, MPI_INT, 0, 99, MPI_COMM_WORLD );
  } else { // for rank 0
    vec_sizes[ 0 ] = graph.size;
    graphs[ 0 ] = graph.vs; // assign local graph for rank 0
    //printf( "size of graph of rank %d: %d\n", myrank, vec_sizes[0] );
    for (i = 1; i < numranks; ++i) { // receive size of vector from other ranks.
      PMPI_Recv( &vec_sizes[i], 1, MPI_INT, i, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      //printf( "size of graph of rank %d: %d\n", i, vec_sizes[i] );
      /* Allocate space for graphs of other ranks. */
      graphs[ i ] = (vertex *) malloc( sizeof(vertex) * vec_sizes[i] );
    }
  }

  if (myrank == 0) {
    for (i = 1; i < numranks; ++i) { // receive graph from other ranks.
      // grapns[i] holds local graph for the rank i.
      PMPI_Recv( graphs[i], vec_sizes[i], vertex_type, i, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
    }
#if DEBUG
    for (i = 0; i < numranks; ++i) { 
      for (j = 0; j < vec_sizes[i]; ++j) {
	printf( "%s %d ", graphs[i][j].name, graphs[i][j].msg_size );
      }
      putchar('\n');
    }
#endif
  } else { // other ranks send graph to rank 0
    PMPI_Send( graph.vs, graph.size, vertex_type, 0, 99, MPI_COMM_WORLD );
  }
}

static void cleanup()
{
  int i;
  /* Clean up ... */
  MPI_Type_free( &vertex_type ); // clean up the constructed MPI type.
  free_vector( &graph ); // free vector
  if (myrank == 0) {
    for (i = 1; i < numranks; ++i) {
      free( graphs[i] );
    }
    free( vec_sizes );
    free( graphs );
  }
}

static void merge_graphs( FILE *dot )
{
  int i, j, msg_size, time_spent;
  vertex *g;
  char *name;
  int *indexes; // indexes holding the start index from which to traverse.

  fprintf( dot, "digraph {\n" );  
  /* First pass: collect vertices and attach labels. */

  // Go through graph for rank 0 and output every vertex to file.
  g = graphs[ 0 ];
  for ( j = 0; j < vec_sizes[0]; ++j ) {
    name = g[ j ].name;
    fprintf( dot, "\t%s [label=\"%s\"];\n", name, get_label(name) );
  }

  // Go through graphs of other ranks, only collect send/recv operations.
  for ( i = 1; i < numranks; ++i ) {
    g = graphs[ i ];
    for ( j = 0; j < vec_sizes[i]; ++j ) {
      name = g[ j ].name;
      if ( is_sendrecv_oper( name ) ) {
	fprintf( dot, "\t%s [label=\"%s\"];\n", name, get_label(name) );
      }
    }
  }

  fprintf( dot, "\n" );
  /* Second pass: collect edges local to individual graphs.
     Ignore MPI_Init and Finalize and all collectives in this pass. */

  for ( i = 0; i < numranks; ++i ) {
    fprintf( dot, "\tsubgraph cluster_%d {\n", i );
    fprintf( dot, "\t\tlabel = \"rank %d\"\n", i );
    fprintf( dot, "\t\tcolor = white\n" );
    g = graphs[ i ];
    for ( j = 1; j < vec_sizes[i] - 2; ++j ) {
      if ( is_sendrecv_oper(g[j].name) &&
	   is_sendrecv_oper(g[j+1].name) ) {
	// time spent between two vertex u -> v is from the end of
	// u to the start of v.
	time_spent = (int) round( g[j+1].start_time - g[j].end_time );

        int t, is_cri = 0;
	for ( t = path_len-1; t > 0; --t ) {
	  if (strcmp(g[j].name, path[t].name) == 0 &&
	      strcmp(g[j+1].name, path[t-1].name) == 0 ) {
	    is_cri = 1;
	    fprintf( dot, "\t\t%s -> %s [label=%d,color=red];\n", g[j].name, g[j+1].name, time_spent);
	  }
	}
	if (!is_cri)
	  fprintf( dot, "\t\t%s -> %s [label=%d];\n", g[j].name, g[j+1].name, time_spent);
      }
    }
    fprintf( dot, "\t}\n" );
  }

  fprintf( dot, "\n" );
  /* Third pass: collect edges between local graphs. */

  indexes = (int *) malloc( numranks * sizeof(int) );
  for ( i = 0; i < numranks; ++i ) {
    indexes[ i ] = 0; // always traverse from the beginning.
  }

  for ( i = 0; i < numranks; ++i ) {
    g = graphs[ i ];
    for ( j = 0; j < vec_sizes[i]; ++j ) {
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

	      int t, is_cri = 0;
	      for ( t = path_len-1; t > 0; --t ) {
		if (strcmp(g[j].name, path[t].name) == 0 &&
		    strcmp(g2[k].name, path[t-1].name) == 0 ) {
		  is_cri = 1;
		  fprintf( dot, "\t%s -> %s [label=\"%d(%d)\",color=red];\n",
			   g[j].name, g2[k].name, latency, msg_size);
		}
	      }
	      if (!is_cri)
		fprintf( dot, "\t%s -> %s [label=\"%d(%d)\"];\n",
			 g[j].name, g2[k].name, latency, msg_size);

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

  fprintf( dot, "\n" );
  /* Final pass: connect vertices in the cluster to collectives. */
  for ( i = 0; i < numranks; ++i ) {
    g = graphs[ i ];
    for ( j = 0; j < vec_sizes[i] - 1; ++j ) {
      if ( !is_sendrecv_oper(g[j].name) ||
	   !is_sendrecv_oper(g[j+1].name) ) {
	// time spent between two vertex u -> v is from the end of
	// u to the start of v.
	time_spent = (int) round( g[j+1].start_time - g[j].end_time );

	int t, is_cri = 0;
	for ( t = path_len; t > 0; --t ) {
	  if (strcmp(g[j].name, path[t].name) == 0 &&
	      strcmp(g[j+1].name, path[t-1].name) == 0 ) {
	    is_cri = 1;
	    fprintf( dot, "\t%s -> %s [label=%d,color=red];\n", g[j].name, g[j+1].name, time_spent);
	  }
	}
	if (!is_cri)
	  fprintf( dot, "\t%s -> %s [label=%d];\n", g[j].name, g[j+1].name, time_spent);
      }
    }
  }  

  fprintf( dot, "}\n" );
  fclose(dot);
  system( "dot -Tpng graph.dot > graph.png" );
  system( "open graph.png" );
}

/* ================== C Wrappers for MPI_Finalize ================== */
_EXTERN_C_ int PMPI_Finalize();
_EXTERN_C_ int MPI_Finalize() { 
  int _wrap_py_return_val = 0, i, j;
  vertex v;
  FILE *dotfile, *stats;

  v.start_time = MPI_Wtime(); // For MPI_Finalize, no need to record end time.
  sprintf( v.name, "%s", __FUNCTION__ );
  append_vector( &graph, &v );

  collect_graphs(); // Let rank 0 collects graphs from other ranks.

  if ( myrank == 0 ) {
    find_critical_path( graphs, vec_sizes, numranks );
    dotfile = fopen( "graph.dot", "w" );
    // Let rank 0 merge local graphs into one graph and output it to dot file.
    merge_graphs( dotfile );

  }

  if (myrank == 0) {
    stats = fopen( "stats.dat", "w" );
    fprintf( stats, "Function\tInvocations\tMean\tMin\tMedian\tMax\n" );
  }

  int init_total;
  double init_sum = 0, init_min = 0, init_max = 0;
  double init_avg, init_min_r, init_max_r;

  for (i = 0; i < init_count; ++i) {
    init_sum += init_times[i];
    if (init_times[i] > init_max)
      init_max = init_times[i];
    if (init_times[i] < init_min)
      init_min = init_times[i];	
  }

  PMPI_Reduce(&init_count, &init_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  PMPI_Reduce(&init_sum, &init_avg, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (init_total > 0)
    init_avg = init_avg / init_total;
  else
    init_avg = 0;
  PMPI_Reduce(&init_min, &init_min_r, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
  PMPI_Reduce(&init_max, &init_max_r, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  if (myrank == 0)
    fprintf( stats, "MPI_Init\t%d\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n",
	     init_total, init_avg, init_min_r, 0.0, init_max_r );

  int send_total;
  double send_sum = 0, send_min = 0, send_max = 0;
  double send_avg, send_min_r, send_max_r;

  for (i = 0; i < send_count; ++i) {
    send_sum += send_times[i];
    if (send_times[i] > send_max)
      send_max = send_times[i];
    if (send_times[i] < send_min)
      send_min = send_times[i];	
  }

  PMPI_Reduce(&send_count, &send_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  PMPI_Reduce(&send_sum, &send_avg, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (send_total > 0)
    send_avg = send_avg / send_total;
  else
    send_avg = 0;
  PMPI_Reduce(&send_min, &send_min_r, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
  PMPI_Reduce(&send_max, &send_max_r, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  if (myrank == 0)
    fprintf( stats, "MPI_Send\t%d\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n",
	     send_total, send_avg, send_min_r, 0.0, send_max_r );

  int recv_total;
  double recv_sum = 0, recv_min = 0, recv_max = 0;
  double recv_avg, recv_min_r, recv_max_r;

  for (i = 0; i < recv_count; ++i) {
    recv_sum += recv_times[i];
    if (recv_times[i] > recv_max)
      recv_max = recv_times[i];
    if (recv_times[i] < recv_min)
      recv_min = recv_times[i];	
  }
  PMPI_Reduce(&recv_count, &recv_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  PMPI_Reduce(&recv_sum, &recv_avg, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (recv_total > 0)
    recv_avg = recv_avg / recv_total;
  else
    recv_avg = 0;
  PMPI_Reduce(&recv_min, &recv_min_r, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
  PMPI_Reduce(&recv_max, &recv_max_r, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  if (myrank == 0)
    fprintf( stats, "MPI_Recv\t%d\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n",
	     recv_total, recv_avg, recv_min_r, 0.0, recv_max_r );

  int wait_total;
  double wait_sum = 0, wait_min = 0, wait_max = 0;
  double wait_avg, wait_min_r, wait_max_r;

  for (i = 0; i < wait_count; ++i) {
    wait_sum += wait_times[i];
    if (wait_times[i] > wait_max)
      wait_max = wait_times[i];
    if (wait_times[i] < wait_min)
      wait_min = wait_times[i];	
  }
  PMPI_Reduce(&wait_count, &wait_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  PMPI_Reduce(&wait_sum, &wait_avg, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (wait_total > 0)
    wait_avg = wait_avg / wait_total;
  else
    wait_avg = 0;
  PMPI_Reduce(&wait_min, &wait_min_r, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
  PMPI_Reduce(&wait_max, &wait_max_r, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  if (myrank == 0)
    fprintf( stats, "MPI_Wait\t%d\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n",
	     wait_total, wait_avg, wait_min_r, 0.0, wait_max_r );

  int barr_total;
  double barr_sum = 0, barr_min = 0, barr_max = 0;
  double barr_avg, barr_min_r, barr_max_r;

  for (i = 0; i < barr_count; ++i) {
    barr_sum += barr_times[i];
    if (barr_times[i] > barr_max)
      barr_max = barr_times[i];
    if (barr_times[i] < barr_min)
      barr_min = barr_times[i];	
  }
  PMPI_Reduce(&barr_count, &barr_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  PMPI_Reduce(&barr_sum, &barr_avg, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (barr_total > 0)
    barr_avg = barr_avg / barr_total;
  else
    barr_avg = 0;
  PMPI_Reduce(&barr_min, &barr_min_r, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
  PMPI_Reduce(&barr_max, &barr_max_r, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  if (myrank == 0)
    fprintf( stats, "MPI_Barrier\t%d\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n",
	     barr_total, barr_avg, barr_min_r, 0.0, barr_max_r );


  int Isend_total;
  double Isend_sum = 0, Isend_min = 0, Isend_max = 0;
  double Isend_avg, Isend_min_r, Isend_max_r;

  for (i = 0; i < Isend_count; ++i) {
    Isend_sum += Isend_times[i];
    if (Isend_times[i] > Isend_max)
      Isend_max = Isend_times[i];
    if (Isend_times[i] < Isend_min)
      Isend_min = Isend_times[i];	
  }
  PMPI_Reduce(&Isend_count, &Isend_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  PMPI_Reduce(&Isend_sum, &Isend_avg, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (Isend_total > 0)
    Isend_avg = Isend_avg / Isend_total;
  else
    Isend_avg = 0;
  PMPI_Reduce(&Isend_min, &Isend_min_r, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
  PMPI_Reduce(&Isend_max, &Isend_max_r, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  if (myrank == 0)
    fprintf( stats, "MPI_Isend\t%d\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n",
	     Isend_total, Isend_avg, Isend_min_r, 0.0, Isend_max_r );

  int Irecv_total;
  double Irecv_sum = 0, Irecv_min = 0, Irecv_max = 0;
  double Irecv_avg, Irecv_min_r, Irecv_max_r;

  for (i = 0; i < Irecv_count; ++i) {
    Irecv_sum += Irecv_times[i];
    if (Irecv_times[i] > Irecv_max)
      Irecv_max = Irecv_times[i];
    if (Irecv_times[i] < Irecv_min)
      Irecv_min = Irecv_times[i];	
  }
  PMPI_Reduce(&Irecv_count, &Irecv_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  PMPI_Reduce(&Irecv_sum, &Irecv_avg, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (Irecv_total > 0)
    Irecv_avg = Irecv_avg / Irecv_total;
  else
    Irecv_avg = 0;
  PMPI_Reduce(&Irecv_min, &Irecv_min_r, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
  PMPI_Reduce(&Irecv_max, &Irecv_max_r, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  if (myrank == 0)
    fprintf( stats, "MPI_Irecv\t%d\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n",
	     Irecv_total, Irecv_avg, Irecv_min_r, 0.0, Irecv_max_r );


  if (myrank == 0)
    fclose( stats );
  cleanup();
  _wrap_py_return_val = PMPI_Finalize();

  return _wrap_py_return_val;
}

