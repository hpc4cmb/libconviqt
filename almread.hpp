#ifndef ALMREAD_HPP

#define ALMREAD_HPP

#include <mpi.h>

#include "alm.h"
#include "alm_dmcio.h"
#include "alm_powspec_tools.h"
#include "xcomplex.h"

#include <iostream>
#include <string>
#include <stdexcept>
#include <fstream>
#include <algorithm>

using namespace std;

template<typename T> int write_alm_bin(string fname, Alm<xcomplex<T> > almT, Alm<xcomplex<T> > almG, Alm<xcomplex<T> > almC, long lmax, long mmax ) {

  long nalm = Alm<xcomplex<T> >::Num_Alms(lmax, mmax);
  size_t nwrite = nalm * sizeof( xcomplex<T> ) / sizeof(char);

  ofstream outfile( fname, ios::out | ios::binary );
  outfile.write( (char*) &( almT(0,0).re ), nwrite );
  outfile.write( (char*) &( almG(0,0).re ), nwrite );
  outfile.write( (char*) &( almC(0,0).re ), nwrite );
  outfile.close();

  return 0;
}

template<typename T> int read_alm_bin(string fname, Alm<xcomplex<T> > &almT, Alm<xcomplex<T> > &almG, Alm<xcomplex<T> > &almC, long lmax, long mmax ) {

  long nalm = Alm<xcomplex<T> >::Num_Alms(lmax, mmax);
  size_t nread = nalm * sizeof( xcomplex<T> ) / sizeof(char);

  almT.Set(lmax, mmax);
  almG.Set(lmax, mmax);
  almC.Set(lmax, mmax);

  ifstream infile( fname, ios::in | ios::binary );
  infile.read( (char*) &( almT(0,0).re ), nread );
  infile.read( (char*) &( almG(0,0).re ), nread );
  infile.read( (char*) &( almC(0,0).re ), nread );
  infile.close();

  return 0;
}

template<typename T> double compare_alm( Alm<xcomplex<T> > &alm1, Alm<xcomplex<T> > &alm2 ) {
  int lmax = min( alm1.Lmax(),  alm2.Lmax() );
  int mmax = min( alm1.Mmax(),  alm2.Mmax() );

  double rms = 0;
  long nnz = 0;
  for ( int ell=0; ell<lmax+1; ++ell ) {
    for ( int m=0; m<mmax+1; ++m ) {
      if ( m > ell ) break;
      double norm = alm1(ell,m).norm();
      if ( norm != 0 ) {
        double diff_norm = (alm1(ell,m) - alm2(ell,m)).norm();
        rms += diff_norm / norm;
        ++nnz;
      }
    }
  }

  return sqrt( rms / nnz );
}

#endif
