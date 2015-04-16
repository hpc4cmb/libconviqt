// Include the autotools-provided configuration macros
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// This is set in config.h
//#ifdef HAVE_MPI
#include "mpi.h"
//#endif

//#include <iostream>
//#include <exception>
#include <math.h>

#include "cconviqt.h"

int main( int argc, char **argv ) {

  MPI_Comm comm = MPI_COMM_WORLD;
  int flag;
  int ntasks=0, rank=0;

  char *det_id = "LFITEST";
  void *beam = conviqt_beam_new();
  void *sky = conviqt_sky_new();
  void *det = conviqt_detector_new_with_id( det_id );
  void *pnt = conviqt_pointing_new();
  double *ppnt;
  void *cnv;

  long lmax=32; // default 5000
  long beamlmax=lmax;
  long beammmax=4; // default 14;
  char pol;
  long ipol;
  double fwhm=4.0;
  char *beamfile = "../data/mb_lfi_30_27_x_rescaled.alm";
  char *skyfile = "../data/slm.fits";
  double epsilon = 1.32495160e-04;

  long ntheta = 3;
  long nphi = 3;
  long npsi = 3;
  long nsamp = ntheta*nphi*npsi;
  long row;
  long itheta, iphi, ipsi;
  double theta, phi, psi;
  
  long Nbetafac=10; // 2400
  long MCSamples=0 ;
  long lmaxOut=32; // 3000
  long order=3; // 5

  char calibrate=1;

  int err;

  MPI_Initialized( &flag );
  if ( flag ) {
    printf( "ERROR: MPI was already initialized\n" );
    return -1;
  }
  
  if ( MPI_Init( &argc, &argv ) ) {
    printf( "ERROR: Failed to initialize MPI\n" );
    return -1;
  }

  if ( MPI_Comm_size( comm, &ntasks ) ) {
    printf( "ERROR: Failed get MPI communicator size\n" );
    return -1;
  }

  if ( MPI_Comm_rank( comm, &rank ) ) {
    printf( "ERROR: Failed to get MPI rank\n" );
    return -1;
  }

  err = conviqt_detector_set_epsilon( det, epsilon );
  if ( err != 0 ) {
    printf( "Failed to set epsilon\n" );
    return -1;
  }

  err = conviqt_pointing_alloc( pnt, 5 * nsamp );
  if ( err != 0 ) {
    printf( "Failed to allocate pointing\n" );
    return -1;
  }
  ppnt = conviqt_pointing_data( pnt );

  for ( ipol=1; ipol<2; ++ipol ) {

    pol=(ipol==0); // default false

    printf( "pol = %i\n", pol );
    
    conviqt_beam_read( beam, beamlmax, beammmax, pol, beamfile, &comm );
    conviqt_sky_read( sky, lmax, pol, skyfile, fwhm, &comm );

    // Populate the pointing array

    row = 0;
    for ( itheta=0; itheta < ntheta; ++itheta ) {
      theta = itheta * (M_PI / (double)ntheta);
      for ( iphi=0; iphi < nphi; ++iphi ) {
	//phi = iphi * (M_2_PI / (double)nphi);
	phi = (double)iphi * (2. * M_PI / (double)nphi);
	for ( ipsi=0; ipsi < npsi; ++ipsi ) {
	  psi = ipsi * (M_PI / (double)npsi);
	  ppnt[row*5+0] = phi; // longitude
	  ppnt[row*5+1] = theta; // latitude
	  ppnt[row*5+2] = psi; // position angle
	  ppnt[row*5+3] = 0; // TOD
	  ppnt[row*5+4] = row; // time
	  ++row;
	}
      }
    }

    printf( "Contents of pnt BEFORE convolving:\n" );
    for ( row=0; row<nsamp; ++row ) {
      printf( "%i : %f %f %f %f %f\n", row, ppnt[row*5+0], ppnt[row*5+1], ppnt[row*5+2], ppnt[row*5+3], ppnt[row*5+4] );
    }
    
    cnv = conviqt_convolver_new( sky, beam, det, pol, lmax, beammmax, Nbetafac, MCSamples, lmaxOut, order, &comm );

    conviqt_convolver_convolve( cnv, pnt, calibrate );

    ppnt = conviqt_pointing_data( pnt );

    printf( "Contents of pnt AFTER convolving:\n" );
    for ( row=0; row<nsamp; ++row ) {
      printf( "%i : %f %f %f %f %f\n", row, ppnt[row*5+0], ppnt[row*5+1], ppnt[row*5+2], ppnt[row*5+3], ppnt[row*5+4] );
    }
    
    if ( rank == 0 ) {
    
      printf( "Convolved TOD:\n" );
      for ( row=0; row < nsamp; ++row ) {
	printf( "%f %f\n",  ppnt[row*5+4], ppnt[row*5+3] );
      }

      if ( pol ) {
	if ( fabs( ppnt[ 0*5+3] -  0.8546349819096275 ) > 1e-6 ) {
          printf( "Row 0 should be 0.8546349819096275, not %f\n", ppnt[ 0*5+3] );
          return -1;
        }
	if ( fabs( ppnt[10*5+3] + 25.53734467183137   ) > 1e-4 ) {
          printf( "Row 10 should be -25.53734467183137, not %f\n", ppnt[ 10*5+3] );
          return -1;
        }
	if ( fabs( ppnt[15*5+3] + 76.04945574990082   ) > 1e-4 ) {
          printf( "Row 15 should be -76.04945574990082, not %f\n", ppnt[ 15*5+3] );
          return -1;
        }
      } else {
	if ( fabs( ppnt[ 0*5+3] -  0.8545846415739397 ) > 1e-6 ) {
          printf( "Row 0 should be 0.8545846415739397, not %f\n", ppnt[ 0*5+3] );
          return -1;
        }
	if ( fabs( ppnt[10*5+3] + 25.20150061107036   ) > 1e-4 ) {
          printf( "Row 10 should be -25.20150061107036, not %f\n", ppnt[ 10*5+3] );
          return -1;
        }
	if ( fabs( ppnt[15*5+3] + 76.14723911261254   ) > 1e-4 ) {
          printf( "Row 15 should be -76.14723911261254, not %f\n", ppnt[ 15*5+3] );
          return -1;
        }
      }
    
    }
    
  } // ipol loop

  if ( MPI_Finalize() ) {
    printf( "ERROR: Failed finalize MPI\n" );
    return -1;
  }

  return 0;
}

