/**
 * Helper functions for MPI Profiler.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

const char * get_label( const char *name )
{
  int i, len = strlen( name ), counter = 0;
  char *label = (char *) malloc( 25 * sizeof(char) );

  for ( i = 0; i <= len; ++i ) {
    if ( name[ i ] == '_' ) {
      counter++;
    }
    
    if ( counter < 2 ) {
      label[ i ] = name[ i ];
    } else {
      label[ i ] = 0;
      break;
    }
  }
  return label;
}

int is_sendrecv_oper( const char *name )
{
  if ( strncmp(name, "MPI_Send", strlen("MPI_Send")) == 0 ) {
    return 1;
  }
  if ( strncmp(name, "MPI_Isend", strlen("MPI_Isend")) == 0 ) {
    return 1;
  }
  if ( strncmp(name, "MPI_Recv", strlen("MPI_Recv")) == 0 ) {
    return 1;
  }
  if ( strncmp(name, "MPI_Irecv", strlen("MPI_Irecv")) == 0 ) {
    return 1;
  }
  if ( strncmp(name, "MPI_Wait", strlen("MPI_Wait")) == 0 ) {
    return 1;
  }

  return 0;
}

int is_send_oper( const char *name )
{
  if ( strncmp(name, "MPI_Send", strlen("MPI_Send")) == 0 ) {
    return 1;
  }
  if ( strncmp(name, "MPI_Isend", strlen("MPI_Isend")) == 0 ) {
    return 1;
  }

  return 0;
}

int is_block_recv_oper( const char *name )
{
  if ( strncmp(name, "MPI_Recv", strlen("MPI_Recv")) == 0 ) {
    return 1;
  }
  if ( strncmp(name, "MPI_Wait", strlen("MPI_Wait")) == 0 ) {
    return 1;
  }

  return 0;
}

/* Return message latency in ms given input byte size. */
int sendrecv_latency( int size )
{
  int latency;
  latency = (int) round( 8.595e-06 * size + 3.262e-01 );
  return latency;
}
