
/**************************************************************************
 *
 *  An introductory MPI program using the default communicator
 *  MPI_COMM_WORLD to print a message to the screen on each rank.
 *
 *  Access to MPI types and constants is via use of the module "mpi"
 *  While we must provide an error flag (here "ifail"), in practice
 *  we never have to check its value on return.
 *
 *  Note that some programs may use the older idiom of
 *  include 'mpif.h'
 *  instead of "use mpi". Some compilers may still require the older
 *  form.
 *
 *************************************************************************/

#include <mpi.h>
#include <stdio.h>

int main (int argc, char* argv[]){

  int rank, size, len;
  char name[MPI_MAX_PROCESSOR_NAME];

  MPI_Init(&argc,&argv);

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Get_processor_name(name, &len);
  printf("Hello from rank %2d of %d on %s\n",rank, size,name);

  MPI_Finalize();
  
  return 0;

} 
