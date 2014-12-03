#ifndef __CONVIQT_HPP__

// Include the autotools-provided configuration macros
#ifdef HAVE_CONFIG_H
#include"config.h"
#endif

#include <new> // std::bad_alloc
#include <utility> // std::swap

// Healpix headers
#include <alm.h>

// LevelS cxxsupport and c_utils headers
#include "mpi_support.h"
#include "alm_fitsio.h"
#include "arr.h"
#include "lsconstants.h"
#include "alm_powspec_tools.h"
#include "wigner.h"
#include "fftpack_support.h"

// convolver accuracy parameter
#define conv_acc 1.e-30

// convolver verbosity
#define CMULT_VERBOSITY 1

// The declarations go here

#define __CONVIQT_HPP__

class beam {
public :
  int read( long beamlmax, long beammmax, bool pol, std::string infile_beam, MPI_Comm comm=MPI_COMM_WORLD );
  Alm< xcomplex<float> > & blmT( void );
  Alm< xcomplex<float> > & blmG( void );
  Alm< xcomplex<float> > & blmC( void );
private :
  Alm< xcomplex<float> > blmT_, blmG_, blmC_;
  long lmax, mmax;
  bool pol;
  std::string fname;
};

class sky {
public :
  int read( long skylmax, bool pol, std::string infile_sky, double fwhm_deconv_sky=0, MPI_Comm comm=MPI_COMM_WORLD );
  Alm< xcomplex<float> > & slmT( void );
  Alm< xcomplex<float> > & slmG( void );
  Alm< xcomplex<float> > & slmC( void );
private :
  Alm<xcomplex<float>> slmT_, slmG_, slmC_;
  long lmax;
  bool pol;
  double fwhm_deconv;
  std::string fname;
};

class pointing : public arr<double> {
  // The storage is 1D but there are 5 virtual columns:
  // row*5 + 0 = longitude (theta)
  // row*5 + 1 = latitude (phi)
  // row*5 + 2 = position angle (psi)
  // row*5 + 3 = signal
  // row*5 * 4 = time
public :
  ;
private :
  ;
};


class detector {
public :
  detector( ) {};
  detector( std::string det_id ) : det_id_(det_id) {};
  void set_epsilon(double epsilon) { epsilon_ = epsilon; };
  double get_epsilon() {return epsilon_;};
  void set_id( std::string det_id ) { det_id_ = det_id; };
  std::string get_id() { return det_id_; };
private :
  std::string det_id_;
  double epsilon_;
};

class convolver {  
public :
  convolver( sky *s, beam *b, detector *d, bool pol=false, long lmax=5000, long beammax=14, long Nbetafac=2400, long MCSamples=0, long lmaxOut=3000, long order=5, MPI_Comm comm=MPI_COMM_WORLD );
  int convolve( pointing & pnt, bool calibrate=true );
  int set_sky( sky *s );
  int set_beam( beam *b );
  int set_detector( detector *d );

  void thetaDeltaThetacm( double thetaini, int cores, int corenum, double &theta, double &deltatheta );
  void deltaTheta2( long iival, double thetaini, arr<double> &dbeta );
  void ratiobetacalcgreatercm( long &ratiobeta, double theta, int cores, arr<double> &corethetaarr );  
  void ratiobetacalcsmallercm( long &ratiobeta, double theta, int cores, arr<double> &corethetaarr );

private :
  MPI_Manager mpiMgr;
  sky *s;
  beam *b;
  detector *d;
  bool pol;
  long lmax, beammmax, Nbetafac, MCSamples, lmaxOut, order;

  void weight_ncm( double x, arr<double> &wgt );
  
  void conviqt_hemiscm_v4( arr3<xcomplex<double> > &tod1, arr3<xcomplex<double> > &tod2, long NThetaIndex1, arr<double> &rthetas );

  void conviqt_hemiscm_pol_v4( arr3<xcomplex<double> > &tod1, arr3<xcomplex<double> > &tod2, long NThetaIndex1, arr<double> &rthetas );
  
  void conviqt_hemiscm_single( arr3<xcomplex<double> > &tod1, long NThetaIndex1, arr<double> &rthetas );

  void conviqt_hemiscm_pol_single( arr3<xcomplex<double> > &tod1, long NThetaIndex1, arr<double> &rthetas );

  void conviqt_hemiscm_pol_fill( arr2<xcomplex<double> > &tod1, arr2<xcomplex<double> > &tod2, double beta, arr<long> &effM, long &countdlm, int corenum, arr<double> &dlmelements);

  void todAnnulus(arr2<xcomplex<double> > &tod1, arr2<xcomplex<double> > &Cmm, arr<double> &cs, arr<double> &sn, arr<double> &cs0, arr<double> &sn0);
  
  void todAnnulus_v3(arr3<xcomplex<double> > &tod1, arr3<xcomplex<double> > &Cmm, arr<double> &cs, arr<double> &sn, arr<double> &cs0, arr<double> &sn0, long NThetaIndex1);

  void conviqt_tod_loop_fill(arr<long> &lowerIndex, arr<long> &upperIndex, arr<double> &outpntarr, arr2<xcomplex<double> > &TODAsym, long thetaIndex, arr<long> &itheta0, long max_order, double inv_delta_theta, double theta0, double inv_delta_phi, double phioffset, long nphi, long ioffset, long npoints, long ntod, arr<double> &TODValue, arr<int> &iphi0arr, arr2<double> &sinweight, arr2<double> &cosweight);
    
  void conviqt_tod_loop_v4(arr<long> &lowerIndex, arr<long> &upperIndex, arr<double> &outpntarr, arr3<xcomplex<double> > &TODAsym, long thetaIndex, arr<long> &itheta0, long max_order, double inv_delta_theta, double theta0, double inv_delta_phi, double phioffset, long nphi, long ioffset, long npoints, long ntod, arr<double> &TODValue, long lat);
  
  void conviqt_tod_loop_pol_v5(arr<long> &lowerIndex, arr<long> &upperIndex, arr<double> &outpntarr, arr3<xcomplex<double> > &TODAsym, long thetaIndex, arr<long> &itheta0, long max_order, double inv_delta_theta, double theta0, double inv_delta_phi, double phioffset, long nphi, long ioffset, long npoints, long ntod, arr<double> &TODValue, long lat);

  void conviqt_hemiscm_pol_fast( arr2<xcomplex<double> > &tod1, arr2<xcomplex<double> > &tod2, double beta, long &countdlm, arr<long> &effM, arr<double> &dp1, int corenum );

  void conviqt_tod_loop_read(arr<long> &lowerIndex, arr<long> &upperIndex, arr2<xcomplex<double> > &TODAsym, long thetaIndex, long nphi, long npoints, arr<double> &TODValue, arr<int> &iphi0arr, arr2<double> &sinweight, arr2<double> &cosweight);

  void ithetacalc(arr<long> &itheta0, arr<double> &outpntarr, long ntod, double inv_delta_theta, double theta0, long ioffset, long ntheta, long npoints);

  void interpolTOD_arrTestcm_pol_fill( arr<double> &outpntarr1, arr<double> &outpntarr2, arr<double> &TODValue, arr<double> &TODValue2, long ntod, long ntod2, int corenum, int cores, arr<double> &betaIni_arr, long &NThetaIndex1, arr<long> &lowerIndex, arr<long> &upperIndex, arr<long> &lowerIndex2, arr<long> &upperIndex2, arr<int> &iphi0arr, arr2<double> &sinweight, arr2<double> &cosweight, arr<int> &iphi0arr2, arr2<double> &sinweight2, arr2<double> &cosweight2, arr<long> &effM, arr<double> &dp1 );

  void interpolTOD_arrTestcm_pol_fast( arr<double> &TODValue, arr<double> &TODValue2, int corenum, arr<double> &betaIni_arr, long &NThetaIndex1, arr<int> &iphi0arr, arr<int> &iphi0arr2, arr2<double> &sinweight, arr2<double> &cosweight, arr2<double> &sinweight2, arr2<double> &cosweight2, arr<long> &lowerIndex, arr<long> &upperIndex, arr<long> &lowerIndex2, arr<long> &upperIndex2, arr<long> &effM, arr<double> &dp1 );

  void interpolTOD_arrTestcm_v4( arr<double> &outpntarr1, arr<double> &outpntarr2, arr<double> &TODValue, arr<double> &TODValue2, long ntod, long ntod2, int corenum );

  void interpolTOD_arrTestcm_pol_v4( arr<double> &outpntarr1, arr<double> &outpntarr2, arr<double> &TODValue, arr<double> &TODValue2, long ntod, long ntod2, int corenum);

  void itheta0SetUp(int npoints, int ioffset, long ntheta, double theta0, double inv_delta_theta, long nphi, arr<double> outpntarr, long ntod, long &NThetaIndex, arr<long> &itheta0, arr<long> &lowerIndex, arr<long> &upperIndex, arr3<xcomplex<double> > &TODAsym);
  
  void fillingBetaSeg( arr < double > & pntarr, long & arrsize, double ratiodeltas, int cores, int corenum, arr<double> &corethetaarr, arr<int> &inBetaSeg);

  void todRedistribution5cm( arr<double> pntarr, arr<int> inBetaSeg, arr<int> outBetaSeg, arr<int> &inBetaSegAcc, arr<int> &outBetaSegAcc, long &outBetaSegSize, int cores, arr<double> &pntarr2, long totsize, arr<int> &inOffset, arr<int> &outOffset, double ratiodeltas, arr<double> &corethetaarr );

  void todgen_fill( long ntod, long ntod2, arr<double> &todTest_arr, arr<double> &timeTest_arr, arr<double> &todTest_arr2, arr<double> &timeTest_arr2, arr<double> &outpntarr, int corenum, int cores, arr<double> &betaIni_arr, long &NThetaIndex1, arr<long> &lowerIndex, arr<long> &upperIndex, arr<long> &lowerIndex2, arr<long> &upperIndex2, arr<int> &iphi0arr, arr2<double> &sinweight, arr2<double> &cosweight, arr<int> &iphi0arr2, arr2<double> &sinweight2, arr2<double> &cosweight2, arr<long> &effM, arr<double> &dp1 );

  void todgen_fast( long ntod, long ntod2, arr<double> &todTest_arr, arr<double> &todTest_arr2, int corenum, arr<double> &betaIni_arr, long &NThetaIndex1, arr<long> &lowerIndex, arr<long> &upperIndex, arr<long> &lowerIndex2, arr<long> &upperIndex2, arr<int> &iphi0arr, arr2<double> &sinweight, arr2<double> &cosweight, arr<int> &iphi0arr2, arr2<double> &sinweight2, arr2<double> &cosweight2, arr<long> &effM, arr<double> &dp1 );

  void todgen_v4( long ntod, long ntod2, arr<double> &todTest_arr, arr<double> &timeTest_arr, arr<double> &todTest_arr2, arr<double> &timeTest_arr2, arr<double> &outpntarr, int corenum );

  void preReorderingStep(long ntod, long ntod2, arr<double> &todAll, arr<double> &todTest_arr, arr<double> &todTest_arr2);

  void effMFiller( double beta, arr<long> &effM );

  void arrsizecounter( arr<long> &effM, long &numberofdlms );

  void arrFillingcm_v2( long ntod, arr<double> &timeTest_arr, arr<double> &outpntarrx, arr<double> &outpntarr, long offindex );
  
};

// wignercalc.cpp

double xpow( int expo, double val );
void wignerCalc( tsize n, tsize mmax, double theta, arr2<double> &d );
void wignerCalcGeneral( tsize n, tsize mmax, double theta, arr2<double> &d );
void wignerCalcHalfpi( tsize n, tsize mmax, arr2<double> &d );
double wignerCalc00_halfpi( tsize n );
double wignerCalc00( double theta, tsize n );

// conviqt_util.cpp

void sift_down_DL( arr<double> &ra, arr<long> &brr, const int l, const int r );
void hpsort_DL( arr<double> &ra, arr<long> &brr );
void sift_down_arrTheta( arr<double> &ra, const int l, const int r );
void hpsort_arrTheta( arr<double> &ra );
void sift_down_arrTOD( arr<double> &ra, const int l, const int r );
void hpsort_arrTOD( arr<double> &ra );
void sift_down_arrTime( arr<double> &ra, const int l, const int r );
void hpsort_arrTime( arr<double> &ra );
void sift_down_DDcm( arr<double> &ra, arr<double> &brr, const int l, const int r );
void hpsort_DDcm( arr<double> &ra, arr<double> &brr );

#endif
