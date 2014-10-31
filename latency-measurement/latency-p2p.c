/**
 * Peer to peer communication in MPI.
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MIN_SIZE 256
#define MAX_SIZE 1024*1024*64

int main( int argc, char *argv[] )
{
  unsigned char * msg;
  int len, i, myrank, msg_size, repetition = 100;
  MPI_Status status;
  double start_time, end_time;
  char hostname[ MPI_MAX_PROCESSOR_NAME ];
  FILE * output;

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
  MPI_Get_processor_name( hostname, &len );

  if (myrank == 0) {
    output = fopen( "latency-sendrecv.csv", "w" );
  }

  //printf( "Rank %d is running on %s\n", myrank, hostname );
  for (msg_size = MIN_SIZE; msg_size <= MAX_SIZE; msg_size *= 2) {
    msg = (unsigned char *) malloc( msg_size*sizeof(unsigned char) );
    for (i = 0; i < msg_size; ++i) {
      msg[ i ] = 0xff;
    }
  
    if (myrank == 0) {
      start_time = MPI_Wtime();
      for (i = 0; i < repetition; ++i) {
	// send message to rank 1. assume MPI_send not blocking
	MPI_Send( msg, msg_size, MPI_UNSIGNED_CHAR, 1, 99, MPI_COMM_WORLD );
	// waiting to receive the same message from rank 1
	MPI_Recv( msg, msg_size, MPI_UNSIGNED_CHAR, 1, 99, MPI_COMM_WORLD, &status );
      }
      end_time = MPI_Wtime();

      fprintf( output, "%d,%.3lf\n", msg_size, (end_time - start_time)/repetition/2 );
    } else if (myrank == 1) {
      for (i = 0; i < repetition; ++i) {
	// waiting to receive the message from rank 0
	MPI_Recv( msg, msg_size, MPI_UNSIGNED_CHAR, 0, 99, MPI_COMM_WORLD, &status );
	// send the same message to rank 0
	MPI_Send( msg, msg_size, MPI_UNSIGNED_CHAR, 0, 99, MPI_COMM_WORLD );
      }
    }
    free( msg );
  }
  
  MPI_Finalize();
  return 0;
}
