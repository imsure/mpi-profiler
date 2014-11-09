/**
 * Test program for MPI Profiler.
 */

#include "mpi.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#define MYTAG 99
#define BUF_SIZE 1000

int main( int argc, char *argv[] )
{
  char * char_buf;
  double * double_buf;
  int myrank, numranks, bufsize, i;
  MPI_Status status;

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
  MPI_Comm_size( MPI_COMM_WORLD, &numranks );

  bufsize = BUF_SIZE;
  char_buf = (char *)malloc( bufsize * sizeof(char) );
  double_buf = (double *)malloc( bufsize * sizeof(double) );

  if (myrank == 0) {
    for (i = 1; i < numranks; ++i) {
      MPI_Send( char_buf, bufsize, MPI_CHAR, i, MYTAG, MPI_COMM_WORLD );
      sleep( 1 );
    }
  } else {
    MPI_Recv( char_buf, bufsize, MPI_CHAR, 0, MYTAG, MPI_COMM_WORLD, &status );
    sleep( 2 );
  }

  MPI_Barrier( MPI_COMM_WORLD );

  MPI_Finalize();
  return 0;
}
