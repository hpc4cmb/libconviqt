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

  long lmax=24; // default 5000
  long beamlmax=lmax;
  long beammmax=5; // default 14;
  long lmaxOut=24; // 3000
  long order=13; // 5
  bool pol=true; // default false
  double fwhm=4.0;
  std::string beamfile( "../data/mb_lfi_30_27_x_rescaled.alm" );
  std::string skyfile( "../data/slm.fits" );

  b.read( beamlmax, beammmax, pol, beamfile, comm_world );
  s.read( lmax, pol, skyfile, fwhm, comm_world );

  // Populate the pointing array

  long ntheta = 100;
  long nphi = 10;
  long npsi = 10;
  long nsamp_self = ntheta*nphi*npsi;
  long nsamp_world = nsamp_self / ntasks;
  long first_sample = rank*nsamp_world;
  if (rank == ntasks-1) nsamp_world = nsamp_self - (ntasks-1)*nsamp_world;

  pnt_self.alloc( 5 * nsamp_self );
  pnt_world.alloc( 5 * nsamp_world );

  std::vector<double> theta(nsamp_self);
  std::vector<double> phi(nsamp_self);
  std::vector<double> psi(nsamp_self);

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

  for (long row=0; row<nsamp_world; ++row) {
    pnt_world[row*5+0] = phi[row+first_sample]; // longitude
    pnt_world[row*5+1] = theta[row+first_sample]; // latitude
    pnt_world[row*5+2] = psi[row+first_sample]; // position angle
    pnt_world[row*5+3] = 0; // TOD
    pnt_world[row*5+4] = row+first_sample; // time
  }
    
  std::cout << rank << " : Convolving self" << std::endl;
  convolver cnv_self( &s, &b, &d, pol, lmax, beammmax, lmaxOut, order, comm_self );
  if (rank == 0) cnv_self.convolve( pnt_self );

  MPI_Barrier( comm_world );

  std::cout << rank << " : Convolving world" << std::endl;
  convolver cnv_world( &s, &b, &d, pol, lmax, beammmax, lmaxOut, order, comm_world );
  cnv_world.convolve( pnt_world );

  std::vector<int> counts(ntasks);
  int my_count = nsamp_world*5;
  int ierr = MPI_Gather( &my_count, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm_world );
  if (ierr != 0) throw std::runtime_error( "Failed to gather counts" );
  std::vector<int> displs(ntasks);
  displs[0] = 0;
  for (int irank=1; irank<ntasks; ++irank) displs[irank] = displs[irank-1] + counts[irank-1];

  pointing pnt_world_tot;
  pnt_world_tot.alloc( 5*nsamp_self );
  ierr = MPI_Gatherv( &(pnt_world[0]), nsamp_world*5, MPI_DOUBLE, &(pnt_world_tot[0]), counts.data(), displs.data(), MPI_DOUBLE, 0, comm_world );
  if (ierr != 0) throw std::runtime_error( "Failed to gather convolved data" );

  if ( rank == 0 ) {

    /*
    std::cout << "Convolved TOD:" << std::endl;
    for ( long row=0; row < nsamp_self; ++row ) {
      std::cout << pnt_world_tot[row*5+4] << " " << pnt_world_tot[row*5+3] << " " << pnt_self[row*5+4] << " " << pnt_self[row*5+3] << std::endl;
    }
    */

    for ( long row=0; row < nsamp_self; ++row ) {
      if ( fabs( pnt_world_tot[row*5+3] - pnt_self[row*5+3] ) > 1e-4 ) {
	std::ostringstream o;
	o << std::setprecision( 16 );
	o << "Row " << row << " ( " << theta[row] << ", " << phi[row] << ", " << psi[row] << ") differs:  " << pnt_self[row*5+3] << " != " << pnt_world_tot[row*5+3];
	//throw std::runtime_error( o.str() );
	std::cout << o.str() << std::endl;
      }
    }

  }
    
  if ( MPI_Finalize() ) throw std::runtime_error( "ERROR: Failed to finalize MPI" );

  return 0;
}

