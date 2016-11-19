#include "conviqt.hpp"

// This file will contain the I/O facilities: read beam, read sky, read detector

// An end user will have the choice between using these facilities or providing inputs directly

using namespace levels;

namespace conviqt {

int beam::read( long beamlmax, long beammmax, bool beampol, std::string infile_beam, MPI_Comm comm ) {

  MPI_Manager mpiMgr( comm );
  int rank;
  MPI_Comm_rank( comm, &rank );
  
  if ( CMULT_VERBOSITY > 1 ) std::cout << "Reading " << infile_beam << " on task = " << rank << std::endl;

  pol = beampol;
  if ( rank == 0 ) {
    int lmax_file = -1;
    int mmax_file = -1;
    if ( beamlmax < 0 || beammmax < 0 ) {
      if ( pol )
	get_almsize_pol( infile_beam, lmax_file, mmax_file );
      else
	get_almsize( infile_beam, lmax_file, mmax_file );
    }
    lmax = beamlmax;
    mmax = beammmax;
    if ( lmax < 0 ) lmax = lmax_file;
    if ( mmax < 0 ) mmax = mmax_file;
  }
  mpiMgr.bcast( lmax );
  mpiMgr.bcast( mmax );
  fname = infile_beam;

  blmT_.Set( lmax, mmax );
  if ( pol ) {
    blmG_.Set( lmax, mmax );
    blmC_.Set( lmax, mmax );
  }

  long blmsize = (mmax+1)*(mmax+2) + 2*(mmax+1)*(lmax-mmax);

  if ( rank == 0 ) {
    pol ? read_Alm_from_fits( infile_beam, blmT_, blmG_, blmC_, lmax, mmax )
      : read_Alm_from_fits( infile_beam, blmT_, lmax, mmax );
    if ( CMULT_VERBOSITY > 1 ) std::cout << "Done reading beam alms on task = " << rank << std::endl;
  }

  mpiMgr.bcastRaw( &blmT_(0,0).re, blmsize, 0 );
  if ( pol ) {
    mpiMgr.bcastRaw( &blmG_(0,0).re, blmsize, 0 );
    mpiMgr.bcastRaw( &blmC_(0,0).re, blmsize, 0 );
  }

  return 0;
}

Alm< xcomplex<float> > & beam::blmT( void ) {
  return blmT_;
}

Alm< xcomplex<float> > & beam::blmG( void ) {
  if ( !pol ) throw std::runtime_error( "Requested polarized components of an unpolarized beam" );
  
  return blmG_;
}

Alm< xcomplex<float> > & beam::blmC( void ) {
  if ( !pol ) throw std::runtime_error( "Requested polarized components of an unpolarized beam" );

  return blmC_;
}

int sky::read( long skylmax, bool skypol, std::string infile_sky, double fwhm_deconv_sky, MPI_Comm comm ) {

  MPI_Manager mpiMgr( comm );
  int rank;
  MPI_Comm_rank( comm, &rank );
  
  if ( CMULT_VERBOSITY > 1 ) std::cout << "Reading " << infile_sky << " on task = " << rank << std::endl;

  pol = skypol;
  if ( rank == 0 ) {
    if ( skylmax < 0 ) {
      int lmax_file = -1;
      int mmax_file = -1;
      if ( pol )
	get_almsize_pol( infile_sky, lmax_file, mmax_file );
      else
	get_almsize( infile_sky, lmax_file, mmax_file );
      lmax = lmax_file;
    } else {
      lmax = skylmax;
    }
  }
  mpiMgr.bcast( lmax );
  fname = infile_sky;
  fwhm_deconv = fwhm_deconv_sky;
  
  slmT_.Set( lmax, lmax );
  if ( pol ) {
    slmG_.Set( lmax, lmax );
    slmC_.Set( lmax, lmax );
  }

  long slmsize = (lmax+1)*(lmax+2);

  if ( rank == 0 ) { 
    pol ? read_Alm_from_fits( infile_sky, slmT_, slmG_, slmC_, lmax, lmax)
      : read_Alm_from_fits (infile_sky, slmT_, lmax, lmax);
    if ( CMULT_VERBOSITY > 1 ) std::cout << "Done reading sky alms on task = " << rank << std::endl;
  }

  mpiMgr.bcastRaw( &slmT_(0,0).re, slmsize, 0 );
  if ( pol ) {
    mpiMgr.bcastRaw( &slmG_(0,0).re, slmsize, 0 );
    mpiMgr.bcastRaw( &slmC_(0,0).re, slmsize, 0 );
  }
  
  double fwhm = arcmin2rad * fwhm_deconv;
  pol ? smoothWithGauss (slmT_, slmG_, slmC_, -fwhm) : smoothWithGauss (slmT_, -fwhm);

  return 0;
}

Alm< xcomplex<float> > & sky::slmT( void ) {
  return slmT_;
}

Alm< xcomplex<float> > & sky::slmG( void ) {
  if ( !pol ) throw std::runtime_error( "Requested polarized components of an unpolarized sky" );
  
  return slmG_;
}

Alm< xcomplex<float> > & sky::slmC( void ) {
  if ( !pol ) throw std::runtime_error( "Requested polarized components of an unpolarized sky" );

  return slmC_;
}

} // namespace conviqt
