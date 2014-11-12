/**
 * A simple vector interface for holding an array of vertex which
 * represents a MPI task graph.
 */

#ifndef _VECTOR_H_
#define _VECTOR_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define INIT_CAPACITY 100

typedef struct vertex vertex;
typedef struct vector vector;

/**
 * A vector for holding an array of 'vertex'.
 */
struct vector {
  int size; // the current size of the vector (number of elements)
  int capacity; // capacity of the vector
  vertex * vs; // array of vertices.
};

/**
 * Represents a vertex in the MPI task graph.
 * A vertex is a MPI operation. 
 */
struct vertex {
  /* The name of the vertex. It must uniquely identify the vertex
     in the (merged) task graph. We encode name as follows:
     1. For MPI_Init & MPI_Finalize, the encoding is:
     
        MPI_Init
	MPI_Finalize
	
	The final merged graph will have only one MPI_Init &
	MPI_Finalize.
	
     2. For MPI collectives, the encoding is:
     
        collective_name_counter
	
	Ex. MPI_Barrier_1 means MPI_Barrier was the first collective
	operation called by all ranks.
	MPI_Scatter_2 means MPI_Scatter was the second collective
	operation called by all ranks.
	We will have a counter for all collectives.
	
     3. For other supported MPI operations, the encoding is:
     
        operation_name_rank_number_id

	Ex. MPI_Send_0_1 means MPI_Send was called by rank 0 with id 1.
	id is used to identify operations inside a rank. It increments by 1
	each time an MPI operation other than MPI_Init, MPI_Finalized or
	MPI Collective is called.
  */
  char name[ 32 ];
  double start_time; // start time of the operation
  double end_time; // end time of the operation

  /*
   * The following fields are for sending/recving operations only.
   */
  int sender_rank; // rank id of the sender
  int receiver_rank; // rank id of the receiver
  int tag; // message tag
  int msg_size; // size of message
};

/* Initialize vector 'vec' with capacity of 'capacity'. */
void init_vector( vector * vec, int capacity );
/* Append vertex 'vt' to the vector 'vec'. */
void append_vector( vector * vec, vertex * vt );
/* Double the capacity of the vector 'vec' if it is full. */
void resize_vector( vector * vec );
/* Cleaup up vector 'vec' */
void free_vector( vector * vec );

#endif // _VECTOR_H_
