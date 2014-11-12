#ifndef _HELPERS_H_
#define _HELPERS_H_

/* Get label of MPI operations by name. */
const char * get_label( const char *name );

/* Decide if the 'name' is a send/recv related operation. */
int is_sendrecv_oper( const char *name );

/* Decide if the 'name' is a sending operation. */
int is_send_oper( const char *name );

/* Decide if the 'name' is a blocking receving operation. */
int is_block_recv_oper( const char *name );
  
#endif
