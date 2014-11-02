/**
 * Measuring latency of MPI_Barrier.
 *
 * Author: Shuo Yang
 */

#include "mpi.h"
#include <stdio.h>

#define SKIP_SIZE 10
#define ITER_SIZE 1000

int main( int argc, char *argv[] )
{
  int i, myrank, skip_num, iter_num, num_nodes;
  double latency, t_start, t_end, sum, avg_latency;

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
  MPI_Comm_size( MPI_COMM_WORLD, &num_nodes );

  skip_num = SKIP_SIZE;
  iter_num = ITER_SIZE;

  /* A few iterations are run without timing to ignore
     any start-up overheads. */
  for (i = 0; i < skip_num; ++i) {
    MPI_Barrier( MPI_COMM_WORLD );
  }

  t_start = MPI_Wtime();
  for (i = 0; i < iter_num; ++i) {
    MPI_Barrier( MPI_COMM_WORLD );
  }
  t_end = MPI_Wtime();

  latency = (t_end - t_start) * 1e3 / iter_num;
  MPI_Reduce( &latency, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

  if (myrank == 0) {
    avg_latency = sum / num_nodes;
    printf( "number of nodes: %d\tbarrier latency: %.3lf (ms)\n", num_nodes, avg_latency );
  }

  MPI_Finalize();
}
