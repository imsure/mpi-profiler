/**
 * Test program for MPI Profiler.
 */

#include "mpi.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#define MYTAG 99
#define DEBUG 0

int main( int argc, char *argv[] )
{
  double *buf_0_1, *buf_2_1;
  char *buf_0_2;
  int *scatter_buf, *recvbuf;
  int myrank, numranks, bufsize_0_1, bufsize_2_1;
  int bufsize_0_2, bufsize_scatter, i;
  MPI_Status status;
  MPI_Request request;

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
  MPI_Comm_size( MPI_COMM_WORLD, &numranks );

  if ( myrank == 0 ) {
    if ( numranks != 3 ) { // make sure there are exactly 3 ranks
      fprintf( stderr, "There should be exactly three ranks!\n", numranks );
      MPI_Abort( MPI_COMM_WORLD, -1 );
    }
  }

  bufsize_0_1 = 1000;
  bufsize_2_1 = 5000;
  bufsize_0_2 = 3000;
  // Allocate buffer for message passing between rank 0 and 1.
  buf_0_1 = (double *) malloc( bufsize_0_1 * sizeof(double) );
  // Allocate buffer for message passing between rank 1 and 2.
  buf_2_1 = (double *) malloc( bufsize_2_1 * sizeof(double) );
  // Allocate buffer for message passing between rank 0 and 2.
  buf_0_2 = (char *) malloc( bufsize_0_2 * sizeof(char) );

  if ( myrank == 0 ) {
    buf_0_1[ 0 ] = 0.1;
    // rank 0 send buf_0_1 to rank 1.
    MPI_Send( buf_0_1, bufsize_0_1, MPI_DOUBLE, 1, MYTAG, MPI_COMM_WORLD );
#if DEBUG
    printf( "Rank %d sends data to rank 1 with first element %lf.\n", myrank, buf_0_1[0] );
#endif

    // rank 0 recevie data back from rank 1.
    MPI_Irecv( buf_0_1, bufsize_0_1, MPI_DOUBLE, 1, MYTAG, MPI_COMM_WORLD, &request );
    MPI_Wait( &request, &status );
#if DEBUG
    printf( "Rank %d receives data from rank 1 with first element %lf.\n", myrank, buf_0_1[0] );
#endif

    // rank 0 send data to rank 2 repeately.
    MPI_Send( buf_0_2, bufsize_0_2, MPI_CHAR, 2, MYTAG, MPI_COMM_WORLD );
    MPI_Send( buf_0_2, bufsize_0_2, MPI_CHAR, 2, MYTAG, MPI_COMM_WORLD );
    MPI_Send( buf_0_2, bufsize_0_2, MPI_CHAR, 2, MYTAG, MPI_COMM_WORLD );
  }

  if ( myrank == 2 ) {
    buf_2_1[ 0 ] = 2.1;
    // rank 2 sends buf_2_1 to rank 1.
    MPI_Send( buf_2_1, bufsize_2_1, MPI_DOUBLE, 1, MYTAG, MPI_COMM_WORLD );
#if DEBUG
    printf( "Rank %d sends data to rank 1 with first element %lf.\n", myrank, buf_2_1[0] );
#endif

    // rank 2 receives data back from rank 1.
    MPI_Recv( buf_2_1, bufsize_2_1, MPI_DOUBLE, 1, MYTAG, MPI_COMM_WORLD, &status );
#if DEBUG
    printf( "Rank %d receives data from rank 1 with first element %lf.\n", myrank, buf_2_1[0] );
#endif

    // rank 2 receives data from rank 0.
    MPI_Recv( buf_0_2, bufsize_0_2, MPI_CHAR, 0, MYTAG, MPI_COMM_WORLD, &status );
    MPI_Recv( buf_0_2, bufsize_0_2, MPI_CHAR, 0, MYTAG, MPI_COMM_WORLD, &status );
    MPI_Recv( buf_0_2, bufsize_0_2, MPI_CHAR, 0, MYTAG, MPI_COMM_WORLD, &status );
  }

  if ( myrank == 1 ) { 
    // rank 1 receves data from rank 0 
    MPI_Recv( buf_0_1, bufsize_0_1, MPI_DOUBLE, 0, MYTAG, MPI_COMM_WORLD, &status );
#if DEBUG
    printf( "Rank %d receives data from rank 0 with first element %lf.\n", myrank, buf_0_1[0] );
#endif
    // rank 1 receves data from rank 2
    MPI_Irecv( buf_2_1, bufsize_2_1, MPI_DOUBLE, 2, MYTAG, MPI_COMM_WORLD, &request );
    MPI_Wait( &request, &status );
#if DEBUG
    printf( "Rank %d receives data from rank 2 with first element %lf.\n", myrank, buf_2_1[0] );
#endif    

    buf_0_1[ 0 ] = 1.0;
    // rank 1 sends data back to rank 0.
    MPI_Send( buf_0_1, bufsize_0_1, MPI_DOUBLE, 0, MYTAG, MPI_COMM_WORLD );

    buf_2_1[ 0 ] = 1.2;
    // rank 1 sends data back to rank 2 and wait.
    MPI_Isend( buf_2_1, bufsize_2_1, MPI_DOUBLE, 2, MYTAG, MPI_COMM_WORLD, &request );
    MPI_Wait( &request, &status );
  }

  // Meet together and continue
  MPI_Barrier( MPI_COMM_WORLD );

  /* Start scattering ... */

  bufsize_scatter = 1000;
  recvbuf = (int *) malloc( bufsize_scatter * sizeof(int) );
  if ( myrank == 0 ) {
    scatter_buf = (int *) malloc( numranks * bufsize_scatter * sizeof(int) );
    for (i = 0; i < numranks; ++i) {
      scatter_buf[ bufsize_scatter*i ] = i + 1000;
    }
  }

  MPI_Scatter( scatter_buf, bufsize_scatter, MPI_INT, recvbuf,
	       bufsize_scatter, MPI_INT, 0, MPI_COMM_WORLD );
#if DEBUG
  printf( "Rank %d receives scattered data from rank 0 with first element being %d.\n"
	  , myrank, recvbuf[ 0 ] );
#endif

  /* End scattering ... */

  MPI_Finalize();
  return 0;
}
