#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"
#include <unistd.h>

int main(int argc, char* argv[]) {
  int x, y, np, me, i, tag = 2;
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  /* Check that we run on exactly four processors */
  if (np != 4) {
    printf("this test uses 4 processors\n");
    MPI_Finalize();	       /* Quit if there is only one processor */
    exit(0);
  }

  x = me;
  if (me == 0) {
    sleep(1);
    MPI_Send(&x, 1, MPI_INT, 1, tag, MPI_COMM_WORLD);
  }
  else if (me != 3) {
    MPI_Recv (&y, 1, MPI_INT, me-1, tag, MPI_COMM_WORLD, &status);
    sleep(1);
    MPI_Send(&x, 1, MPI_INT, me+1, tag, MPI_COMM_WORLD);
  }
  else {
    MPI_Recv (&y, 1, MPI_INT, 2, tag, MPI_COMM_WORLD, &status);
    sleep(1);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  exit(0);
}
