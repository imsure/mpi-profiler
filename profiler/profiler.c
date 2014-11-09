
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

_EXTERN_C_ void pmpi_init(MPI_Fint *ierr);
_EXTERN_C_ void PMPI_INIT(MPI_Fint *ierr);
_EXTERN_C_ void pmpi_init_(MPI_Fint *ierr);
_EXTERN_C_ void pmpi_init__(MPI_Fint *ierr);

typedef struct vertex vertex;
typedef struct sendrecv sendrecv;

/* Metadata related to the send/recv operations. */
struct sendrecv {
  int sender_rank; // rank id of the sender
  int receiver_rank; // rank id of the receiver
  int tag; // message tag
  int msg_size; // size of message
};

/* Represents a vertex in the MPI task graph. */
struct vertex {
  char name[ 50 ]; // name of the operation
  double start_time; // time of the start of the operation
  double end_time; // time of the end of the operation
  sendrecv * sr;
  vertex * next; // next vertex
};


static vertex head;
static vertex * tail = &head;
static int myrank;

/* file that contains local task graph for the current rank. */
static FILE * local_graph = NULL;
static char fname[20];

static void init_tail( const char * func_name )
{
  tail->next = (vertex *) malloc( sizeof(vertex) );
  tail = tail->next;
  tail->sr = NULL;
  tail->next = NULL;
  strcpy( tail->name, func_name );
}

int get_mpi_datatype_size( MPI_Datatype datatype )
{
  return 0;
}

/* ================== C Wrappers for MPI_Init ================== */
_EXTERN_C_ int PMPI_Init(int *argc, char ***argv);
_EXTERN_C_ int MPI_Init(int *argc, char ***argv) { 
  int _wrap_py_return_val = 0;
  char rank[5];

  strcpy( head.name, __FUNCTION__ );

  head.start_time = MPI_Wtime();
  _wrap_py_return_val = PMPI_Init(argc, argv);
  head.end_time = MPI_Wtime();
  
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

  sprintf( fname, "rank%d.out", myrank ); // construct file name for output
  local_graph = fopen( fname, "w" );
  fprintf( local_graph, "myrank=%d\n", myrank );
  
  return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Send ================== */
_EXTERN_C_ int PMPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
_EXTERN_C_ int MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) { 
  int _wrap_py_return_val = 0, size;

  init_tail( __FUNCTION__ );

  /* Recording send/recv operation. */
  tail->sr = (sendrecv *) malloc( sizeof(sendrecv) );
  tail->sr->sender_rank = myrank;
  tail->sr->receiver_rank = dest;
  tail->sr->tag = tag;
  MPI_Type_size( datatype, &size );
  tail->sr->msg_size = count * size;

  tail->start_time = MPI_Wtime();
  _wrap_py_return_val = PMPI_Send(buf, count, datatype, dest, tag, comm);
  tail->end_time = MPI_Wtime();
  
  return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Recv ================== */
_EXTERN_C_ int PMPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
_EXTERN_C_ int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status) { 
  int _wrap_py_return_val = 0;

  init_tail( __FUNCTION__ );

  /* Recording send/recv operation. */
  tail->sr = (sendrecv *) malloc( sizeof(sendrecv) );
  tail->sr->sender_rank = source;
  tail->sr->receiver_rank = myrank;
  tail->sr->tag = tag;
  // ignore message size on the receiver side since it is already stored on the sender side.
  tail->sr->msg_size = -1; 

  tail->start_time = MPI_Wtime();
  _wrap_py_return_val = PMPI_Recv(buf, count, datatype, source, tag, comm, status);
  tail->end_time = MPI_Wtime();
  
  return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Barrier ================== */
_EXTERN_C_ int PMPI_Barrier(MPI_Comm comm);
_EXTERN_C_ int MPI_Barrier(MPI_Comm comm) { 
  int _wrap_py_return_val = 0;

  init_tail( __FUNCTION__ );
  
  tail->start_time = MPI_Wtime();
  _wrap_py_return_val = PMPI_Barrier(comm);
  tail->end_time = MPI_Wtime();
  
  return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Alltoall ================== */
_EXTERN_C_ int PMPI_Alltoall(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm);
_EXTERN_C_ int MPI_Alltoall(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) { 
  int _wrap_py_return_val = 0;

  init_tail( __FUNCTION__ );
  
  tail->start_time = MPI_Wtime();
  _wrap_py_return_val = PMPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
  tail->end_time = MPI_Wtime();
  
  return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Scatter ================== */
_EXTERN_C_ int PMPI_Scatter(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
_EXTERN_C_ int MPI_Scatter(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) { 
  int _wrap_py_return_val = 0;

  init_tail( __FUNCTION__ );

  tail->start_time = MPI_Wtime();
  _wrap_py_return_val = PMPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
  tail->end_time = MPI_Wtime();
  
  return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Gather ================== */
_EXTERN_C_ int PMPI_Gather(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
_EXTERN_C_ int MPI_Gather(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) { 
  int _wrap_py_return_val = 0;

  init_tail( __FUNCTION__ );

  tail->start_time = MPI_Wtime();
  _wrap_py_return_val = PMPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
  tail->end_time = MPI_Wtime();

  return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Reduce ================== */
_EXTERN_C_ int PMPI_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
_EXTERN_C_ int MPI_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) { 
  int _wrap_py_return_val = 0;

  init_tail( __FUNCTION__ );
  
  tail->start_time = MPI_Wtime();
  _wrap_py_return_val = PMPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
  tail->end_time = MPI_Wtime();

  return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Allreduce ================== */
_EXTERN_C_ int PMPI_Allreduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
_EXTERN_C_ int MPI_Allreduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) { 
  int _wrap_py_return_val = 0;

  init_tail( __FUNCTION__ );

  tail->start_time = MPI_Wtime();
  _wrap_py_return_val = PMPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
  tail->end_time = MPI_Wtime();

  return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Isend ================== */
_EXTERN_C_ int PMPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
_EXTERN_C_ int MPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request) { 
  int _wrap_py_return_val = 0;

  init_tail( __FUNCTION__ );

  /* Recording send/recv operation. */
  tail->sr = (sendrecv *) malloc( sizeof(sendrecv) );
  tail->sr->sender_rank = myrank;
  tail->sr->receiver_rank = dest;
  tail->sr->tag = tag;
  tail->sr->msg_size = count * sizeof(datatype);

  tail->start_time = MPI_Wtime();
  _wrap_py_return_val = PMPI_Isend(buf, count, datatype, dest, tag, comm, request);
  tail->end_time = MPI_Wtime();

  return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Irecv ================== */
_EXTERN_C_ int PMPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request);
_EXTERN_C_ int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request) { 
  int _wrap_py_return_val = 0;

  init_tail( __FUNCTION__ );
  
  tail->start_time = MPI_Wtime();
  _wrap_py_return_val = PMPI_Irecv(buf, count, datatype, source, tag, comm, request);
  tail->end_time = MPI_Wtime();

  return _wrap_py_return_val;
}

/* ================== C Wrappers for MPI_Wait ================== */
_EXTERN_C_ int PMPI_Wait(MPI_Request *request, MPI_Status *status);
_EXTERN_C_ int MPI_Wait(MPI_Request *request, MPI_Status *status) { 
  int _wrap_py_return_val = 0;

  init_tail( __FUNCTION__ );

  /* Recording send/recv operation. */
  tail->sr = (sendrecv *) malloc( sizeof(sendrecv) );
  tail->sr->sender_rank = status->MPI_SOURCE;
  tail->sr->receiver_rank = myrank;
  tail->sr->tag = status->MPI_TAG;
  // ignore message size on the receiver side since it is already stored on the sender side.
  tail->sr->msg_size = -1; 

  tail->start_time = MPI_Wtime();
  _wrap_py_return_val = PMPI_Wait(request, status);
  tail->end_time = MPI_Wtime();
  
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
  int _wrap_py_return_val = 0;

  init_tail( __FUNCTION__ );

  tail->start_time = MPI_Wtime();
  _wrap_py_return_val = PMPI_Finalize();
  tail->end_time = MPI_Wtime();

  /* Traverse the local task graph and write each vertex to
     the local output file, one vertex per line. */
  vertex * h = &head;
  while (h != NULL) {
    if (h->sr == NULL) {
      fprintf( local_graph, "vertex=%s start_time=%lf end_time=%lf\n",
	       h->name, h->start_time, h->end_time );
    } else {
      fprintf( local_graph, "vertex=%s start_time=%lf end_time=%lf \
sender_rank=%d receiver_rank=%d tag=%d msg_size=%d\n",
	       h->name, h->start_time, h->end_time, h->sr->sender_rank,
	       h->sr->receiver_rank, h->sr->tag, h->sr->msg_size );
    }
    h = h->next;
  }
  fclose( local_graph );

  return _wrap_py_return_val;
}

