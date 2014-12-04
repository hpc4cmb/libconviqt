// Include the autotools-provided configuration macros
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// This is set in config.h
//#ifdef HAVE_MPI
#include "mpi.h"
//#endif

#include <iostream>
#include <exception>

#include "conviqt.hpp"

int main( int argc, char **argv ) {

  MPI_Comm comm = MPI_COMM_WORLD;

  int flag;
  MPI_Initialized( &flag );
  if ( flag ) throw std::runtime_error( "ERROR: MPI was already initialized" );
  
  if ( MPI_Init( &argc, &argv ) ) throw std::runtime_error( "ERROR: Failed to initialize MPI" );

  int ntasks=0, rank=0;
  if ( MPI_Comm_size( comm, &ntasks ) ) throw std::runtime_error( "ERROR: Failed get MPI communicator size" );
  if ( MPI_Comm_rank( comm, &rank ) ) throw std::runtime_error( "ERROR: Failed to get MPI rank" );

  beam b;
  sky s;
  detector d( "LFITEST" );
  pointing pnt;

  d.set_epsilon( 1.32495160e-04 );

  for ( int ipol=0; ipol<2; ++ipol ) {
    
    long lmax=32; // default 5000
    long beamlmax=lmax;
    long beammmax=4; // default 14;
    bool pol=(ipol==0); // default false
    double fwhm=4.0;
    std::string beamfile( "../data/mb_lfi_30_27_x_rescaled.alm" );
    std::string skyfile( "../data/slm.fits" );
    
    b.read( beamlmax, beammmax, pol, beamfile, comm );
    s.read( lmax, pol, skyfile, fwhm, comm );

    // Populate the pointing array

    long ntheta = 3;
    long nphi = 3;
    long npsi = 3;
    long nsamp = ntheta*nphi*npsi;

    pnt.alloc( 5 * nsamp );

    long row=0;
    for ( long itheta=0; itheta < ntheta; ++itheta ) {
      double theta = itheta * (pi / (double)ntheta);
      for ( long iphi=0; iphi < nphi; ++iphi ) {
	double phi = iphi * (twopi / (double)nphi);
	for ( long ipsi=0; ipsi < npsi; ++ipsi ) {
	  double psi = ipsi * (pi / (double)npsi);
	  pnt[row*5+0] = phi; // longitude
	  pnt[row*5+1] = theta; // latitude
	  pnt[row*5+2] = psi; // position angle
	  pnt[row*5+3] = 0; // TOD
	  pnt[row*5+4] = row; // time
	  ++row;
	}
      }
    }
    
    long Nbetafac=10; // 2400
    long MCSamples=0 ;
    long lmaxOut=32; // 3000
    long order=3; // 5
    
    convolver cnv( &s, &b, &d, pol, lmax, beammmax, Nbetafac, MCSamples, lmaxOut, order, comm );

    cnv.convolve( pnt );

    if ( rank == 0 ) {
    
      std::cout << "Convolved TOD:" << std::endl;
      for ( long row=0; row < nsamp; ++row ) {
	std::cout << pnt[row*5+4] << " " << pnt[row*5+3] << std::endl;
      }

      if ( pol ) {
	if ( fabs( pnt[ 0*5+3] -  0.854635 ) > 1e-6 ) throw std::runtime_error( "Row 0 should be 0.854635, not "  + std::to_string(pnt[ 0*5+3]) );
	if ( fabs( pnt[10*5+3] + 25.5373   ) > 1e-4 ) throw std::runtime_error( "Row 10 should be -25.5373, not " + std::to_string(pnt[10*5+3]) );
	if ( fabs( pnt[15*5+3] + 76.0495   ) > 1e-4 ) throw std::runtime_error( "Row 15 should be -76.0495, not " + std::to_string(pnt[15*5+3]) );
      } else {
	if ( fabs( pnt[ 0*5+3] -  0.854585 ) > 1e-6 ) throw std::runtime_error( "Row 0 should be 0.854585, not "  + std::to_string(pnt[ 0*5+3]) );
	if ( fabs( pnt[10*5+3] + 25.2015   ) > 1e-4 ) throw std::runtime_error( "Row 10 should be -25.2015, not " + std::to_string(pnt[10*5+3]) );
	if ( fabs( pnt[15*5+3] + 76.1472   ) > 1e-4 ) throw std::runtime_error( "Row 15 should be -76.1472, not " + std::to_string(pnt[15*5+3]) );
      }
    
    }
    
  } // ipol loop

  if ( MPI_Finalize() ) throw std::runtime_error( "ERROR: Failed finalize MPI" );

  return 0;
}

