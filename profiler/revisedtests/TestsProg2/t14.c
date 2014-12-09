#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

/* Gather data from a vector to contiguous */

int main( int argc, char **argv )
{
    MPI_Datatype vec;
    MPI_Comm comm;
    double *vecin, *vecout;
    int minsize = 2, count;
    int root, i, n, stride, errs = 0,chk = 0;
    int rank, size;
 
    MPI_Init( &argc, &argv );
    comm = MPI_COMM_WORLD;
    /* Determine the sender and receiver */
    MPI_Comm_rank( comm, &rank );
    MPI_Comm_size( comm, &size );
 
    if (size != 4) {
      printf("this test uses 4 processors\n");
      MPI_Finalize();            /* Quit if there is only one processor */
      exit(0);
    }

    for (root=0; root<size; root++) {
        for (count = 1; count < 30; count = count * 2) {
            n = 12;
            stride = 10;
            vecin = (double *)malloc( n * stride * size * sizeof(double) );
            vecout = (double *)malloc( size * n * sizeof(double) );
 
            MPI_Type_vector( n, 1, stride, MPI_DOUBLE, &vec );
            MPI_Type_commit( &vec );

            for (i=0; i<n*stride; i++) vecin[i] =-2;
            for (i=0; i<n; i++) vecin[i*stride] = rank * n + i;

	    if (rank == 0)
	      sleep(2);

            MPI_Gather( vecin, 1, vec, vecout, n, MPI_DOUBLE, root, comm );
			chk++;

	    if (rank == size-1)
	      sleep(2);

            if (rank == root) {
                for (i=0; i<n*size; i++) {
                    if (vecout[i] != i) {
                        errs++;
                        if (errs < 2) {
                            fprintf( stderr, "vecout[%d]=%d\n", i, (int)vecout[i] );fflush(stderr);
                        }
                    }
                }
            }
            MPI_Type_free( &vec );
            free( vecin );
            free( vecout );
        }
    }
 
    MPI_Finalize();
    return 0;
}

