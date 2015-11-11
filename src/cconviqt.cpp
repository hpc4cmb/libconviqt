#include "mpi.h"

#include "cconviqt.h"
#include "conviqt.hpp"

// The C interface must not throw exceptions at calling programs that may not be able to handle them

extern "C" {
  
  void *conviqt_beam_new() { return new(std::nothrow) conviqt::beam; }
  
  int conviqt_beam_del( void *ptr ) {
    try {
      conviqt::beam *ref = reinterpret_cast< conviqt::beam * >( ptr );
      delete ref;
    } catch (...) {
      return -1;
    }
    return 0;
  }
  
  int conviqt_beam_read( void *ptr, long beamlmax, long beammmax, char pol, char *infile_beam, MPI_Comm *comm ) {
    try {
      conviqt::beam *ref = reinterpret_cast< conviqt::beam * >( ptr );
      //return ref->read( beamlmax, beammmax, pol, infile_beam, *comm );
      ref->read( beamlmax, beammmax, pol, infile_beam, MPI_COMM_WORLD );
    } catch (...) {
      return -1;
    }
    return 0;
  }

  void *conviqt_sky_new() { return new(std::nothrow) conviqt::sky; }

  int conviqt_sky_del( void *ptr ) {
    try {
      conviqt::sky *ref = reinterpret_cast< conviqt::sky * >( ptr );
      delete ref;
    } catch (...) {
      return -1;
    }
    return 0;
  }
  
  int conviqt_sky_read( void *ptr, long skylmax, char pol, char *infile_sky, double fwhm_deconv_sky, MPI_Comm *comm ) {
    try {
      conviqt::sky *ref = reinterpret_cast< conviqt::sky * >( ptr );
      //return ref->read( skylmax, pol, infile_sky, fwhm_deconv_sky, *comm );
      ref->read( skylmax, pol, infile_sky, fwhm_deconv_sky, MPI_COMM_WORLD );
    } catch (...) {
      return -1;
    }
    return 0;
  }
  
  void *conviqt_detector_new() { return new(std::nothrow) conviqt::detector; }

  void *conviqt_detector_new_with_id( char *det_id ) { return new(std::nothrow) conviqt::detector(det_id); }

  int conviqt_detector_del( void *ptr ) {
    try {
      conviqt::detector *ref = reinterpret_cast< conviqt::detector * >( ptr );
      delete ref;
    } catch (...) {
      return -1;
    }
    return 0;
  }
  
  int conviqt_detector_set_epsilon( void *ptr, double epsilon ) {
    try {
      conviqt::detector *ref = reinterpret_cast< conviqt::detector * >( ptr );
      ref->set_epsilon( epsilon );
    } catch (...) {
      return -1;
    }
    return 0;
  }

  int conviqt_detector_get_epsilon( void *ptr, double *epsilon ) {
    try {
      conviqt::detector *ref = reinterpret_cast< conviqt::detector * >( ptr );
      *epsilon = ref->get_epsilon();
    } catch (...) {
      return -1;
    }
    return 0;
  }

  int conviqt_detector_set_id( void *ptr, char *det_id ) {
    try {
      conviqt::detector *ref = reinterpret_cast< conviqt::detector * >( ptr );
      ref->set_id( det_id );
    } catch (...) {
      return -1;
    }
    return 0;
  }

  int conviqt_detector_get_id( void *ptr, char *det_id ) {
    try {
      conviqt::detector *ref = reinterpret_cast< conviqt::detector * >( ptr );
      strcpy( det_id, ref->get_id().c_str() );
    } catch (...) {
      return -1;
    }
    return 0;
  }

  void *conviqt_pointing_new() { return new(std::nothrow) conviqt::pointing; }

  int conviqt_pointing_del( void *ptr ) {
    try {
      conviqt::pointing *ref = reinterpret_cast< conviqt::pointing * >( ptr );
      delete ref;
    } catch (...) {
      return -1;
    }
    return 0;
  }
  
  int conviqt_pointing_alloc( void *ptr, long n ) {
    std::cout << "Allocating " << n << " doubles for pointing" << std::endl;
    try {
      conviqt::pointing *ref = reinterpret_cast< conviqt::pointing * >( ptr );
      ref->alloc( n );
    } catch (...) {
      return -1;
    }
    return 0;
  }

  double *conviqt_pointing_data( void *ptr ) {
    try {
      conviqt::pointing *ref = reinterpret_cast< conviqt::pointing * >( ptr );
      return ref->begin();
    } catch (...) {
      return NULL;
    }
  }

  void *conviqt_convolver_new( void *skyptr, void *beamptr, void *detptr, char pol, long lmax, long beammmax, long nbetafac, long mcsamples, long lmaxout, long order, MPI_Comm *comm ) {

    std::cout << "conviqt_convolver_new called with pol = " << int(pol) << ", lmax = " << lmax << ", beammmax = " << beammmax << ", nbetafac = " << nbetafac << ", mcsamples = " << mcsamples << ", lmaxout = " << lmaxout << ", order = " << order << std::endl;

    if ( lmax > LMAXMAX || beammmax > LMAXMAX || nbetafac > LMAXMAX || mcsamples > LMAXMAX || lmaxout> LMAXMAX ) {
      std::cerr << "Suspiciously large convolver parameters: lmax = " << lmax << ", beammax = " << beammmax << ", nbetafac = " << nbetafac << ", mcsamples = " << mcsamples << ", lmaxout = " << lmaxout << std::endl;
      return NULL;
    }

    conviqt::sky *skyref = reinterpret_cast< conviqt::sky * >( skyptr );
    conviqt::beam *beamref = reinterpret_cast< conviqt::beam * >( beamptr );
    conviqt::detector *detref = reinterpret_cast< conviqt::detector * >( detptr );

    return new(std::nothrow) conviqt::convolver( skyref, beamref, detref, pol, lmax, beammmax, nbetafac, mcsamples, lmaxout, order, MPI_COMM_WORLD );
  }

  int conviqt_convolver_del( void *ptr ) {
    try {
      conviqt::convolver *ref = reinterpret_cast< conviqt::convolver * >( ptr );
      delete ref;
    } catch (...) {
      return -1;
    }
    return 0;
  }
  
  int conviqt_convolver_set_sky( void *cnvptr, void *skyptr ) {
    try {
      conviqt::convolver *cnvref = reinterpret_cast< conviqt::convolver * >( cnvptr );
      conviqt::sky *skyref = reinterpret_cast< conviqt::sky * >( skyptr );
      cnvref->set_sky( skyref );
    } catch (...) {
      return -1;
    }
    return 0;
  }

  int conviqt_convolver_set_beam( void *cnvptr, void *beamptr ) {
    try {
      conviqt::convolver *cnvref = reinterpret_cast< conviqt::convolver * >( cnvptr );
      conviqt::beam *beamref = reinterpret_cast< conviqt::beam * >( beamptr );
      cnvref->set_beam( beamref );
    } catch (...) {
      return -1;
    }
    return 0;
  }

  int conviqt_convolver_set_detector( void *cnvptr, void *detptr ) {
    try {
      conviqt::convolver *cnvref = reinterpret_cast< conviqt::convolver * >( cnvptr );
      conviqt::detector *detref = reinterpret_cast< conviqt::detector * >( detptr );
      cnvref->set_detector( detref );
    } catch (...) {
      return -1;
    }
    return 0;
  }

  int conviqt_convolver_convolve( void *cnvptr, void *pntptr, char calibrate ) {
    try {
      conviqt::convolver *cnvref = reinterpret_cast< conviqt::convolver * >( cnvptr );
      conviqt::pointing *pntref = reinterpret_cast< conviqt::pointing * >( pntptr );

      std::cout << "conviqt_convolver_convolve called. pnt.size() = " << pntref->size() << std::endl;
      for ( int i = 0; i<10; ++i ) {
        std::cout << "pnt[" << i << "] = " << (*pntref)[i] << std::endl;
      }

      cnvref->convolve( *pntref, calibrate );
    } catch (...) {
      std::cout << "convolve threw an exception:" << std::endl;
      throw;
      return -1;
    }
    return 0;
  }


}
