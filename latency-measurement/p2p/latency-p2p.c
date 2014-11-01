/**
 * Measuring message latency of Peer to peer communication in MPI.
 *
 * Author: Shuo Yang
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MIN_SIZE 32 // 32 bytes
#define MAX_SIZE 1024*1024*32 // 32 MB
#define ROUND_TRIP 100
#define NUM_CASES 4
#define TAG 99

/**
 * Measure latency between MPI_Send and MPI_Recv.
 */
void latency_send_recv( unsigned char * msg, int msg_size,
		       int myrank, int roundtrip, FILE * output )
{
  int i;
  double start_time, end_time;
  MPI_Status status;

  if (myrank == 0) {
    start_time = MPI_Wtime();
    for (i = 0; i < roundtrip; ++i) {
      // send message to rank 1. assume MPI_send not blocking
      MPI_Send( msg, msg_size, MPI_UNSIGNED_CHAR, 1, TAG, MPI_COMM_WORLD );
      // waiting to receive the same message from rank 1
      MPI_Recv( msg, msg_size, MPI_UNSIGNED_CHAR, 1, TAG, MPI_COMM_WORLD, &status );
    }
    end_time = MPI_Wtime();

    fprintf( output, "%d,%.3lf\n", msg_size, (end_time - start_time)*1000/roundtrip/2 );
  } else if (myrank == 1) {
    for (i = 0; i < roundtrip; ++i) {
      // waiting to receive the message from rank 0
      MPI_Recv( msg, msg_size, MPI_UNSIGNED_CHAR, 0, TAG, MPI_COMM_WORLD, &status );
      // send the same message to rank 0
      MPI_Send( msg, msg_size, MPI_UNSIGNED_CHAR, 0, TAG, MPI_COMM_WORLD );
    }
  }
  free( msg );
}

/**
 * Measure latency between MPI_Send and MPI_Irecv.
 */
void latency_send_Irecv( unsigned char * msg, int msg_size,
			 int myrank, int roundtrip, FILE * output )
{
  int i;
  double start_time, end_time;
  MPI_Status status;
  MPI_Request request;

  if (myrank == 0) {
    start_time = MPI_Wtime();
    for (i = 0; i < roundtrip; ++i) {
      // send message to rank 1. assume MPI_send not blocking
      MPI_Send( msg, msg_size, MPI_UNSIGNED_CHAR, 1, TAG, MPI_COMM_WORLD );
      // waiting to receive the same message from rank 1
      MPI_Irecv( msg, msg_size, MPI_UNSIGNED_CHAR, 1, TAG, MPI_COMM_WORLD, &request );
      MPI_Wait( &request, &status );
    }
    end_time = MPI_Wtime();

    fprintf( output, "%d,%.3lf\n", msg_size, (end_time - start_time)*1000/roundtrip/2 );
  } else if (myrank == 1) {
    for (i = 0; i < roundtrip; ++i) {
      // waiting to receive the message from rank 0
      MPI_Irecv( msg, msg_size, MPI_UNSIGNED_CHAR, 0, TAG, MPI_COMM_WORLD, &request );
      MPI_Wait( &request, &status );
      // send the same message to rank 0
      MPI_Send( msg, msg_size, MPI_UNSIGNED_CHAR, 0, TAG, MPI_COMM_WORLD );
    }
  }
  free( msg );
}

/**
 * Measure latency between MPI_Isend and MPI_Recv.
 */
void latency_Isend_recv( unsigned char * msg, int msg_size,
			 int myrank, int roundtrip, FILE * output )
{
  int i;
  double start_time, end_time;
  MPI_Status status;
  MPI_Request request;

  if (myrank == 0) {
    start_time = MPI_Wtime();
    for (i = 0; i < roundtrip; ++i) {
      // send message to rank 1.
      MPI_Isend( msg, msg_size, MPI_UNSIGNED_CHAR, 1, TAG, MPI_COMM_WORLD, &request );
      MPI_Wait( &request, &status );
      // waiting to receive the same message from rank 1
      MPI_Recv( msg, msg_size, MPI_UNSIGNED_CHAR, 1, TAG, MPI_COMM_WORLD, &status );
    }
    end_time = MPI_Wtime();

    fprintf( output, "%d,%.3lf\n", msg_size, (end_time - start_time)*1000/roundtrip/2 );
  } else if (myrank == 1) {
    for (i = 0; i < roundtrip; ++i) {
      // waiting to receive the message from rank 0
      MPI_Recv( msg, msg_size, MPI_UNSIGNED_CHAR, 0, TAG, MPI_COMM_WORLD, &status );
      // send the same message to rank 0
      MPI_Isend( msg, msg_size, MPI_UNSIGNED_CHAR, 0, TAG, MPI_COMM_WORLD, &request );
      MPI_Wait( &request, &status );
    }
  }
  free( msg );
}

/**
 * Measure latency between MPI_Isend and MPI_Irecv.
 */
void latency_Isend_Irecv( unsigned char * msg, int msg_size,
			 int myrank, int roundtrip, FILE * output )
{
  int i;
  double start_time, end_time;
  MPI_Status status;
  MPI_Request request;

  if (myrank == 0) {
    start_time = MPI_Wtime();
    for (i = 0; i < roundtrip; ++i) {
      // send message to rank 1.
      MPI_Isend( msg, msg_size, MPI_UNSIGNED_CHAR, 1, TAG, MPI_COMM_WORLD, &request );
      MPI_Wait( &request, &status );
      // waiting to receive the same message from rank 1
      MPI_Irecv( msg, msg_size, MPI_UNSIGNED_CHAR, 1, TAG, MPI_COMM_WORLD, &request );
      MPI_Wait( &request, &status );
    }
    end_time = MPI_Wtime();

    fprintf( output, "%d,%.3lf\n", msg_size, (end_time - start_time)*1000/roundtrip/2 );
  } else if (myrank == 1) {
    for (i = 0; i < roundtrip; ++i) {
      // waiting to receive the message from rank 0
      MPI_Irecv( msg, msg_size, MPI_UNSIGNED_CHAR, 0, TAG, MPI_COMM_WORLD, &request );
      MPI_Wait( &request, &status );
      // send the same message to rank 0
      MPI_Isend( msg, msg_size, MPI_UNSIGNED_CHAR, 0, TAG, MPI_COMM_WORLD, &request );
      MPI_Wait( &request, &status );
    }
  }
  free( msg );
}


int main( int argc, char *argv[] )
{
  unsigned char * msg;
  int len, i, j, myrank, msg_size, roundtrip = ROUND_TRIP;
  char hostname[ MPI_MAX_PROCESSOR_NAME ];
  FILE * outputs[ NUM_CASES ];
  char * fnames[ NUM_CASES ];  
 
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
  MPI_Get_processor_name( hostname, &len );

  if (myrank == 0) {
    fnames[ 0 ] = "latency-sendrecv.csv";
    fnames[ 1 ] = "latency-sendIrecv.csv";
    fnames[ 2 ] = "latency-Isendrecv.csv";
    fnames[ 3 ] = "latency-IsendIrecv.csv";

    for (i = 0; i < NUM_CASES; ++i)
      outputs[ i ] = fopen( fnames[i], "w" );
  }

  //printf( "Rank %d is running on %s\n", myrank, hostname );

  for (j = 0; j < NUM_CASES; ++j) {
    for (msg_size = MIN_SIZE; msg_size <= MAX_SIZE; msg_size *= 2) {
      msg = (unsigned char *) malloc( msg_size*sizeof(unsigned char) );
      for (i = 0; i < msg_size; ++i) {
	msg[ i ] = 0xff;
      }

      switch ( j ) {
      case 0:
	latency_send_recv( msg, msg_size, myrank, roundtrip, outputs[j] );
	break;
      case 1:
	latency_send_Irecv( msg, msg_size, myrank, roundtrip, outputs[j] );
	break;
      case 2:
	latency_Isend_recv( msg, msg_size, myrank, roundtrip, outputs[j] );
	break;
      case 3:
	latency_Isend_Irecv( msg, msg_size, myrank, roundtrip, outputs[j] );
	break;
      }
    }
  }

  if (myrank == 0) {
    for (i = 0; i < NUM_CASES; ++i)
      fclose( outputs[i] );
  }

  MPI_Finalize();
  return 0;
}
