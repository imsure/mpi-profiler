#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"
#include <unistd.h>

int main(int argc, char* argv[]) {
  int x, y, np, me, i;
  int tag = 42;
  MPI_Status  status;
  MPI_Request request, request2;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  /* Check that we run on exactly two processors */
  if (np != 2) {
    printf("this test uses 2 processors\n");
    MPI_Finalize();	       /* Quit if there is only one processor */
    exit(0);
  }

  if (me == 0) {
    MPI_Isend(&x, 1, MPI_INT, 1, tag, MPI_COMM_WORLD, &request);
    sleep(8);
    MPI_Wait(&request, &status);
  } else {
    MPI_Recv (&y, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
    sleep(5);
  }
  
  MPI_Finalize();
  exit(0);
}
