// Include the autotools-provided configuration macros
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "mpi.h"

#include <iostream>
#include <iomanip>
#include <exception>
#include <sstream>

#include "conviqt.hpp"

using namespace conviqt;

int main( int argc, char **argv ) {

  MPI_Comm comm_world = MPI_COMM_WORLD;
  MPI_Comm comm_self = MPI_COMM_SELF;

  std::cout << "Checking that MPI is not yet initialized." << std::endl;

  int flag;
  MPI_Initialized( &flag );
  if ( flag ) throw std::runtime_error( "ERROR: MPI was already initialized" );
  
  std::cout << "Initializing MPI" << std::endl;

  if ( MPI_Init( &argc, &argv ) ) throw std::runtime_error( "ERROR: Failed to initialize MPI" );

  int ntasks=0, rank=0;
  if ( MPI_Comm_size( comm_world, &ntasks ) ) throw std::runtime_error( "ERROR: Failed get MPI communicator size" );
  if ( MPI_Comm_rank( comm_world, &rank ) ) throw std::runtime_error( "ERROR: Failed to get MPI rank" );

  std::cout << std::setprecision( 16 );

  beam b;
  sky s;
  detector d( "LFITEST" );
  pointing pnt_world, pnt_self;

  d.set_epsilon( 1.32495160e-04 );

  long lmax=32; // default 5000
  long beamlmax=lmax;
  long beammmax=4; // default 14;
  bool pol=true; // default false
  double fwhm=4.0;
  std::string beamfile( "../data/mb_lfi_30_27_x_rescaled.alm" );
  std::string skyfile( "../data/slm.fits" );

  //b.read( beamlmax, beammmax, pol, beamfile, comm_world );
  //s.read( lmax, pol, skyfile, fwhm, comm_world );

  // Populate the pointing array

  long ntheta = 3;
  long nphi = 3;
  long npsi = 3;
  long nsamp_self = ntheta*nphi*npsi;
  long nsamp_world = nsamp_self / ntasks;
  long first_sample = rank*nsamp_world;
  if (rank == ntasks-1) nsamp_world = nsamp_self - (ntasks-1)*nsamp_world;

  pnt_self.alloc( 5 * nsamp_self );
  //pnt_world.alloc( 5 * nsamp_world );

  std::vector<double> theta(nsamp_world);
  std::vector<double> phi(nsamp_world);
  std::vector<double> psi(nsamp_world);

  long isamp=0;
  for ( long itheta=0; itheta < ntheta; ++itheta ) {
    double dtheta = itheta * (pi / (double)ntheta);
    for ( long iphi=0; iphi < nphi; ++iphi ) {
      double dphi = iphi * (twopi / (double)nphi);
      for ( long ipsi=0; ipsi < npsi; ++ipsi ) {
	double dpsi = ipsi * (pi / (double)npsi);
	theta[isamp] = dtheta;
	phi[isamp] = dphi;
	psi[isamp] = dpsi;
	isamp++;
      }
    }
  }
  
  for (long row=0; row<nsamp_self; ++row) {
    pnt_self[row*5+0] = phi[row]; // longitude
    pnt_self[row*5+1] = theta[row]; // latitude
    pnt_self[row*5+2] = psi[row]; // position angle
    pnt_self[row*5+3] = 0; // TOD
    pnt_self[row*5+4] = row; // time
  }

  /*
  for (long row=0; row<nsamp_world; ++row) {
    pnt_world[row*5+0] = phi[row+first_sample]; // longitude
    pnt_world[row*5+1] = theta[row+first_sample]; // latitude
    pnt_world[row*5+2] = psi[row+first_sample]; // position angle
    pnt_world[row*5+3] = 0; // TOD
    pnt_world[row*5+4] = row+first_sample; // time
  }
  */
    
  long Nbetafac=10; // 2400
  long MCSamples=0 ;
  long lmaxOut=32; // 3000
  long order=3; // 5
  
  //std::cout << rank << " : Instantiating convolver." << std::endl;

  std::cout << rank << " : Convolving self" << std::endl;
  //convolver cnv_self( &s, &b, &d, pol, lmax, beammmax, Nbetafac, MCSamples, lmaxOut, order, comm_self );
  //cnv_self.convolve( pnt_self );

  /*
    std::cout << rank << " : Convolving world" << std::endl;
    convolver cnv_world( &s, &b, &d, pol, lmax, beammmax, Nbetafac, MCSamples, lmaxOut, order, comm_world );
    cnv_world.convolve( pnt_world );
  }
  */

  /*
  int counts[ntasks];
  int my_count = nsamp_world*5;
  int ierr = MPI_Gather( &my_count, 1, MPI_INT, counts, 1, MPI_INT, 0, comm_world );
  if (ierr != 0) throw std::runtime_error( "Failed to gather counts" );
  int displs[ntasks];
  displs[0] = 0;
  for (int irank=1; irank<ntasks; ++irank) displs[irank] = displs[irank-1] + counts[irank-1];
  */

  /*
  pointing pnt_world_tot;
  pnt_world_tot.alloc( 5*nsamp_self );
  ierr = MPI_Gatherv( &(pnt_world[0]), nsamp_world*5, MPI_DOUBLE, &(pnt_world_tot[0]), counts, displs, MPI_DOUBLE, 0, comm_world );
  if (ierr != 0) throw std::runtime_error( "Failed to gather convolved data" );
  */

  /*
  if ( rank == 0 ) {
    
    std::cout << "Convolved TOD:" << std::endl;
    for ( long row=0; row < nsamp_self; ++row ) {
      std::cout << pnt_world_tot[row*5+4] << " " << pnt_world_tot[row*5+3] << " " << pnt_self[row*5+4] << " " << pnt_self[row*5+3] << std::endl;
    }
    
    if ( fabs( pnt_world_tot[ 0*5+3] -  0.8546349819096275 ) > 1e-6 ) {
      std::ostringstream o;
      o << std::setprecision( 16 );
      o << "Row 0 should be 0.8546349819096275, not " << pnt_world_tot[0*5+3] << " or " << pnt_self[0*5+3];
      throw std::runtime_error( o.str() );
    }
    if ( fabs( pnt_world_tot[10*5+3] + 25.53734467183137   ) > 1e-4 ) {
	std::ostringstream o;
	o << std::setprecision( 16 );
	o << "Row 10 should be -25.53734467183137, not " << pnt_world_tot[10*5+3] << " or " << pnt_self[10*5+3];
	throw std::runtime_error( o.str() );
    }
    if ( fabs( pnt_world_tot[15*5+3] + 76.04945574990082   ) > 1e-4 ) {
      std::ostringstream o;
      o << std::setprecision( 16 );
      o << "Row 15 should be -76.04945574990082, not " << pnt_world_tot[15*5+3] << " or " << pnt_self[15*5+3];
      throw std::runtime_error( o.str() );
    }
  }
  */
    
  if ( MPI_Finalize() ) throw std::runtime_error( "ERROR: Failed to finalize MPI" );

  return 0;
}

