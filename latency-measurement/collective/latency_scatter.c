/**
 * Measuring latency of MPI_Alltoall.
 *
 * Author: Shuo Yang
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SKIP_SIZE 10
#define ITER_SIZE 100
#define MSG_SIZE_MIN 32
#define MSG_SIZE_MAX 1024*32

#define DEBUG 0

int main( int argc, char *argv[] )
{
  int i, myrank, skip_num, iter_num, num_nodes, msg_size;
  double latency, t_start, t_end, sum, avg_latency;
  char * send_buf;
  char * recv_buf;

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
  MPI_Comm_size( MPI_COMM_WORLD, &num_nodes );

  skip_num = SKIP_SIZE;
  iter_num = ITER_SIZE;

  for (msg_size = MSG_SIZE_MIN; msg_size <= MSG_SIZE_MAX; msg_size *= 2) {
    if (myrank == 0) {
      int buf_size = msg_size * num_nodes * sizeof(char) ;
      send_buf = (char *) malloc( buf_size );
      memset( send_buf, 0x7f, buf_size );
    }
    recv_buf = (char *) malloc( msg_size * sizeof(char) );
  
    MPI_Barrier( MPI_COMM_WORLD );

    /* A few iterations are run without timing to ignore
       any start-up overheads. */
    for (i = 0; i < skip_num; ++i) {
      MPI_Scatter( send_buf, msg_size, MPI_CHAR, recv_buf, msg_size, MPI_CHAR, 0, MPI_COMM_WORLD );
#if DEBUG
      printf( "rank %d: the first byte of receive buffer is: %d\n", myrank, recv_buf[0] );
#endif
    }

    t_start = MPI_Wtime();
    for (i = 0; i < iter_num; ++i) {
      MPI_Scatter( send_buf, msg_size, MPI_CHAR, recv_buf, msg_size, MPI_CHAR, 0, MPI_COMM_WORLD );
    }
    t_end = MPI_Wtime();

    latency = (t_end - t_start) * 1e3 / iter_num;
#if DEBUG
    printf( "rank: %d\tscatter latency: %.3lf (ms)\n", myrank, latency );
#endif
    MPI_Reduce( &latency, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

    if (myrank == 0) {
      avg_latency = sum / num_nodes;
      printf( "number of nodes: %d\tmessage size: %d\tscatter latency: %.3lf (ms)\n",
	      num_nodes, msg_size, avg_latency );
      free( send_buf );
    }
    free( recv_buf );
    MPI_Barrier( MPI_COMM_WORLD );
  }
  
  MPI_Finalize();
}
