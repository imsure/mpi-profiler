/**
 * Test result: Failed.
 * Error: *** glibc detected *** ./app: malloc():
 * smallbin double linked list corrupted: 0x0000000000e95db0 ***
 */

#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Finalize();
  exit(0);
}
