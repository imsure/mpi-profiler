/**
 * A simple vector interface. Only support init, append and resize.
 */

#include "vector.h"

/* Initialize vector 'vec' with capacity of 'capacity'. */
void init_vector( vector * vec, int capacity )
{
  if (vec != NULL) {
    vec->size = 0;
    vec->capacity = capacity;
    vec->vs = (vertex *) malloc( capacity * sizeof(vertex) );
  }
}

/* Append vertex 'vt' to the vector 'vec'. */
void append_vector( vector * vec, vertex * vt )
{
  if (vec->size == vec->capacity) {
    resize_vector( vec );
  }

  strcpy( vec->vs[ vec->size ].name, vt->name ) ;
  vec->vs[ vec->size ].start_time = vt->start_time;
  vec->vs[ vec->size ].end_time = vt->end_time;
  vec->vs[ vec->size ].sender_rank = vt->sender_rank;
  vec->vs[ vec->size ].receiver_rank = vt->receiver_rank;
  vec->vs[ vec->size ].tag = vt->tag;
  vec->vs[ vec->size ].msg_size = vt->msg_size;
  ++(vec->size);
}
  
/* Double the capacity of the vector 'vec' if it is full. */
void resize_vector( vector * vec )
{
  vec->capacity *= 2;
  vec->vs = (vertex *) realloc( vec->vs,
				vec->capacity * sizeof(vertex) );
}

/* Clean up the vector 'vec' */
void free_vector( vector * vec )
{
  free( vec->vs );
}
