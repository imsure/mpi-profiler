
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vector.h"

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
  sprintf( v.name, "%s:-1", __FUNCTION__ );
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

  sprintf( v.name, "%s:%d:%d", __FUNCTION__, myrank, sendrecv_id++ );
  v.sender_rank = myrank;
  v.receiver_rank = dest;
  v.tag = tag;
  MPI_Type_size( datatype, &size ); // get size of message in bytes
  v.msg_size = count * size;
  
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

  sprintf( v.name, "%s:%d:%d", __FUNCTION__, myrank, sendrecv_id++ );
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

  sprintf( v.name, "%s:%d:%d", __FUNCTION__, NO_RANK, collective_id++ );
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

  sprintf( v.name, "%s:%d:%d", __FUNCTION__, NO_RANK, collective_id++ );
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

  sprintf( v.name, "%s:%d:%d", __FUNCTION__, NO_RANK, collective_id++ );
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

  sprintf( v.name, "%s:%d:%d", __FUNCTION__, NO_RANK, collective_id++ );
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

  sprintf( v.name, "%s:%d:%d", __FUNCTION__, NO_RANK, collective_id++ );
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

  sprintf( v.name, "%s:%d:%d", __FUNCTION__, NO_RANK, collective_id++ );
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

  sprintf( v.name, "%s:%d:%d", __FUNCTION__, myrank, sendrecv_id++ );
  v.sender_rank = myrank;
  v.receiver_rank = dest;
  v.tag = tag;
  MPI_Type_size( datatype, &size ); // get size of message in bytes
  v.msg_size = count * size;

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

  sprintf( v.name, "%s:%d:%d", __FUNCTION__, myrank, sendrecv_id++ );
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

  sprintf( v.name, "%s:%d:%d", __FUNCTION__, myrank, sendrecv_id++ );
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

/* ================== C Wrappers for MPI_Finalize ================== */
_EXTERN_C_ int PMPI_Finalize();
_EXTERN_C_ int MPI_Finalize() { 
  int _wrap_py_return_val = 0, i, j;
  vertex v;

  v.start_time = MPI_Wtime(); // For MPI_Finalize, no need to record end time.
  sprintf( v.name, "%s:%d", __FUNCTION__, NO_RANK );
  append_vector( &graph, &v );

  /* printf( "size of vector: %d\n", graph.size ); */
  /* for (i = 0; i < graph.size; ++i) { */
  /*   printf( "%s %d ", graph.vs[i].name, graph.vs[i].msg_size ); */
  /* } */
  /* putchar( '\n' ); */

  /* First, send size of vector to rank 0 from other ranks. */
  if (myrank != 0) {
    PMPI_Send( &graph.size, 1, MPI_INT, 0, 99, MPI_COMM_WORLD );
  } else { // for rank 0
    vec_sizes[ 0 ] = graph.size;
    graphs[ 0 ] = graph.vs; // assign local graph for rank 0
    printf( "size of graph of rank %d: %d\n", myrank, vec_sizes[0] );
    for (i = 1; i < numranks; ++i) { // receive size of vector from other ranks.
      PMPI_Recv( &vec_sizes[i], 1, MPI_INT, i, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      printf( "size of graph of rank %d: %d\n", i, vec_sizes[i] );
      /* Allocate space for graphs of other ranks. */
      graphs[ i ] = (vertex *) malloc( sizeof(vertex) * vec_sizes[i] );
    }
  }

  if (myrank == 0) {
    for (i = 1; i < numranks; ++i) { // receive graph from other ranks.
      PMPI_Recv( graphs[i], vec_sizes[i], vertex_type, i, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      for (j = 0; j < vec_sizes[i]; ++j) {
	printf( "%s %d ", graphs[i][j].name, graphs[i][j].msg_size );
      }
      putchar('\n');
    }
  } else { // other ranks send graph to rank 0
    PMPI_Send( graph.vs, graph.size, vertex_type, 0, 99, MPI_COMM_WORLD );
  }

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
  _wrap_py_return_val = PMPI_Finalize();

  return _wrap_py_return_val;
}

