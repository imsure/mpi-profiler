/**
 * Test result: Failed.
 * Error: *** glibc detected *** ./app: malloc():
 * smallbin double linked list corrupted: 0x0000000000e95db0 ***
 */

#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"
#include <unistd.h>

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  sleep(3);
  MPI_Finalize();
  exit(0);
}
