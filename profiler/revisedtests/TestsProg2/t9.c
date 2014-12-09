#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
 
/* A simple test of Reduce with all choices of root process */
int main( int argc, char *argv[] )
{
    int errs = 0;
    int rank, size, root;
    int *sendbuf, *recvbuf, i;
    int minsize = 2, count;
    MPI_Comm comm;
 
    MPI_Init( &argc, &argv );
 
    comm = MPI_COMM_WORLD;
    /* Determine the sender and receiver */
    MPI_Comm_rank( comm, &rank );
    MPI_Comm_size( comm, &size );

    /* Check that we run on exactly two processors */
    if (size != 4) {
      printf("this test uses 4 processors\n");
      MPI_Finalize();
      exit(0);
    }
 
    sendbuf = (int *)malloc(sizeof(int) );
    recvbuf = (int *)malloc(sizeof(int) );
    *sendbuf = rank;
    sleep(rank);
    MPI_Reduce( sendbuf, recvbuf, 1, MPI_INT, MPI_SUM, 0, comm );
 
    MPI_Finalize();
    return errs;
}
