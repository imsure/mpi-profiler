#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1
#endif

int main( int argc, char *argv[] )
{
    int rank, size;
    int chunk = 2;
    int i;
    int *sb;
    int *rb;
    int status, gstatus;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    if (size != 4) {
      printf("this test uses 4 processors\n");
      MPI_Finalize();	       /* Quit if there is only one processor */
      exit(0);
    }

    sb = (int *)malloc(size*chunk*sizeof(int));
    if ( !sb ) {
        perror( "can't allocate send buffer" );fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }
    rb = (int *)malloc(size*chunk*sizeof(int));
    if ( !rb ) {
        perror( "can't allocate recv buffer");fflush(stderr);
        free(sb);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    for ( i=0 ; i < size*chunk ; ++i ) {
        sb[i] = i;
        rb[i] = 0;
    }
    sleep(4-rank);
    status = MPI_Alltoall(sb, chunk, MPI_INT, rb, chunk, MPI_INT, MPI_COMM_WORLD);
    free(sb);
    free(rb);
    MPI_Finalize();
    return(EXIT_SUCCESS);
}

