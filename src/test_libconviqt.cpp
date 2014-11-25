// Include the autotools-provided configuration macros
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// This is set in config.h
#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include <iostream>

#include "conviqt.hpp"

int main( int argc, char **argv ) {

  beam b;
  sky s;
  pointing pnt;
  convolver cnv(s,b);

  int ntasks=0, rank=0, ierr;
  
#ifdef HAVE_MPI
  ierr = MPI_Init( &argc, &argv );
  ierr = MPI_Comm_size( MPI_COMM_WORLD, &ntasks );
  ierr = MPI_Comm_rank( MPI_COMM_WORLD, &rank );
#endif

  cnv.convolve( pnt );

  std::cout << rank << " : Successfully tested libconviqt!" << std::endl;

#ifdef HAVE_MPI
  ierr = MPI_Finalize();
#endif

  return 0;
}

