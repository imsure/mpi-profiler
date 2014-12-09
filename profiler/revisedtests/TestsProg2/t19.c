#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"
#include <unistd.h>

int main(int argc, char* argv[]) {
  int x, y, np, me, i;
  int tag = 42;
  MPI_Status  status;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);



  /* Check that we run on exactly three processors */
  if (np != 3) {
    printf("this test uses 3 processors\n");
    MPI_Finalize();	       /* Quit if there is only one processor */
    exit(0);
  }

  //rank 0: mpi_init->mpi_send(to rank 1)->mpi_barrier->mpi_finalize
  //rank1: mpi_init->mpi_recv(from rank 0)->mpi_send(to rank 2)->mpi_recv(from rank 2)->mpi_barrier->mpi_finalize
  //rank2: mpi_init->mpi_recv(from rank 1)->mpi_send(to rank 1)->mpi_barrier->mpi_finalize

  if (me == 0) {
    sleep(1);
    MPI_Send(&x, 1, MPI_INT, 1, tag, MPI_COMM_WORLD);
    sleep(1);
  }
  else if (me == 1) {
    MPI_Recv (&y, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
    sleep(1);
    MPI_Send(&x, 1, MPI_INT, 2, tag, MPI_COMM_WORLD);
    MPI_Recv (&y, 1, MPI_INT, 2, tag, MPI_COMM_WORLD, &status);
    sleep(1);
  } else {   // me == 2
    MPI_Recv (&y, 1, MPI_INT, 1, tag, MPI_COMM_WORLD, &status);
    sleep(1);
    MPI_Send(&x, 1, MPI_INT, 1, tag, MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  exit(0);
}
