#ifndef __CONVIQT_HPP__

// Include the autotools-provided configuration macros
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <mpi.h>

#include <new> // std::bad_alloc
#include <utility> // std::swap
#include <stdexcept> //std::runtime_error

// Healpix headers
#include "ls_alm.h"
#include "ls_alm_fitsio.h"
#include "ls_alm_powspec_tools.h"

// LevelS cxxsupport and c_utils headers
#include "ls_constants.h"
#include "ls_arr.h"
#include "ls_mpi_support.h"
#include "ls_wigner.h"
#include "ls_fftpack_support.h"

// convolver accuracy parameter (Wigner matrix elements smaller than
// conv_acc are ignored)
#define conv_acc 1.e-30

// convolver verbosity
#define CMULT_VERBOSITY 1

// The declarations go here

#define __CONVIQT_HPP__

using namespace levels;

namespace conviqt {

class beam {
public :
    int read( long beamlmax, long beammmax, bool pol, std::string infile_beam, MPI_Comm comm=MPI_COMM_WORLD );
    Alm< xcomplex<float> > & blmT( void );
    Alm< xcomplex<float> > & blmG( void );
    Alm< xcomplex<float> > & blmC( void );
    int get_lmax( void ) { return lmax; }
    int get_mmax( void ) { return mmax; }
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
    int get_lmax( void ) { return lmax; }
    void remove_monopole( void );
    void remove_dipole( void );
private :
    Alm< xcomplex<float> > slmT_, slmG_, slmC_;
    long lmax;
    bool pol;
    double fwhm_deconv;
    std::string fname;
};

class pointing : public levels::arr<double> {
    // The storage is 1D but there are 5 virtual columns:
    // row*5 + 0 = longitude (theta)
    // row*5 + 1 = latitude (phi)
    // row*5 + 2 = position angle (psi)
    // row*5 + 3 = signal
    // row*5 * 4 = time
public :

private :

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
    convolver( sky *s, beam *b, detector *d, bool pol=false, long lmax=5000, long beammax=14, long order=5, MPI_Comm comm=MPI_COMM_WORLD );
    int convolve( pointing & pnt, bool calibrate=true );
    int set_sky( sky *s );
    int set_beam( beam *b );
    int set_detector( detector *d );

    void thetaDeltaThetacm( double thetaini, int corenum, double &theta, double &deltatheta );
    void deltaTheta2( long iival, double thetaini, levels::arr<double> &dbeta );
    void ratiobetacalcgreatercm( long &ratiobeta, double theta, levels::arr<double> &corethetaarr );
    void ratiobetacalcsmallercm( long &ratiobeta, double theta, levels::arr<double> &corethetaarr );

private :
    MPI_Manager mpiMgr;
    int corenum, cores;
    sky *s;
    beam *b;
    detector *d;
    bool pol;
    long lmax, beammmax, order;
    std::vector<double> base_wgt;

    double t_convolve, t_todgen_v4, t_arrFillingcm_v2, t_interpolTOD_arrTestcm_v4,
        t_itheta0SetUp, t_ithetacalc, t_conviqt_tod_loop_v4, t_weight_ncm,
        t_conviqt_hemiscm_single, t_conviqt_hemiscm_v4,
        t_interpolTOD_arrTestcm_pol_v4, t_conviqt_hemiscm_pol_v4, t_conviqt_hemiscm_pol_single,
        t_conviqt_tod_loop_pol_v5, t_todAnnulus_v3,
        t_wigner_init, t_wigner_prepare, t_wigner_calc,
        t_lat_iter, t_sincos_iter;

    void weight_ncm( double x, levels::arr<double> &wgt );
    void weight_ncm( double x, std::vector<double> &wgt );

    void conviqt_hemiscm_v4( levels::arr3<xcomplex<double> > &tod1, levels::arr3<xcomplex<double> > &tod2, long NThetaIndex1, levels::arr<double> &rthetas );

    void conviqt_hemiscm_pol_v4( levels::arr3<xcomplex<double> > &tod1, levels::arr3<xcomplex<double> > &tod2, long NThetaIndex1, levels::arr<double> &rthetas );

    void conviqt_hemiscm_single( levels::arr3<xcomplex<double> > &tod1, long NThetaIndex1, levels::arr<double> &rthetas );

    void conviqt_hemiscm_pol_single( levels::arr3<xcomplex<double> > &tod1, long NThetaIndex1, levels::arr<double> &rthetas );

    void todAnnulus_v3(levels::arr3<xcomplex<double> > &tod1, levels::arr3<xcomplex<double> > &Cmm, levels::arr<double> &cs, levels::arr<double> &sn, levels::arr<double> &cs0, levels::arr<double> &sn0, long NThetaIndex1);

    void conviqt_tod_loop_v4(levels::arr<long> &lowerIndex, levels::arr<long> &upperIndex, levels::arr<double> &outpntarr, levels::arr3<xcomplex<double> > &TODAsym, long thetaIndex, levels::arr<long> &itheta0, long max_order, double inv_delta_theta, double theta0, double inv_delta_phi, double phioffset, long nphi, long ioffset, long npoints, long ntod, levels::arr<double> &TODValue, long lat);

    void conviqt_tod_loop_pol_v5(levels::arr<long> &lowerIndex, levels::arr<long> &upperIndex, levels::arr<double> &outpntarr, levels::arr3<xcomplex<double> > &TODAsym, long thetaIndex, levels::arr<long> &itheta0, long max_order, double inv_delta_theta, double theta0, double inv_delta_phi, double phioffset, long nphi, long ioffset, long npoints, long ntod, levels::arr<double> &TODValue, long lat);

    void ithetacalc(levels::arr<long> &itheta0, levels::arr<double> &outpntarr, long ntod, double inv_delta_theta, double theta0, long ioffset, long ntheta, long npoints);

    void interpolTOD_arrTestcm_v4( levels::arr<double> &outpntarr1, levels::arr<double> &outpntarr2, levels::arr<double> &TODValue, levels::arr<double> &TODValue2, long ntod, long ntod2 );

    void interpolTOD_arrTestcm_pol_v4( levels::arr<double> &outpntarr1, levels::arr<double> &outpntarr2, levels::arr<double> &TODValue, levels::arr<double> &TODValue2, long ntod, long ntod2 );

    void itheta0SetUp(int npoints, int ioffset, long ntheta, double theta0, double inv_delta_theta, long nphi, levels::arr<double> outpntarr, long ntod, long &NThetaIndex, levels::arr<long> &itheta0, levels::arr<long> &lowerIndex, levels::arr<long> &upperIndex, levels::arr3<xcomplex<double> > &TODAsym);

    void fillingBetaSeg( levels::arr<double> & pntarr, long & arrsize, double ratiodeltas, levels::arr<double> &corethetaarr, levels::arr<int> &inBetaSeg);

    void todRedistribution5cm( levels::arr<double> pntarr, levels::arr<int> inBetaSeg, levels::arr<int> outBetaSeg, levels::arr<int> &inBetaSegAcc, levels::arr<int> &outBetaSegAcc, long &outBetaSegSize, levels::arr<double> &pntarr2, long totsize, levels::arr<int> &inOffset, levels::arr<int> &outOffset, double ratiodeltas, levels::arr<double> &corethetaarr );

    void todgen_v4( long ntod, long ntod2, levels::arr<double> &todTest_arr, levels::arr<double> &timeTest_arr, levels::arr<double> &todTest_arr2, levels::arr<double> &timeTest_arr2, levels::arr<double> &outpntarr );

    void preReorderingStep(long ntod, long ntod2, levels::arr<double> &todAll, levels::arr<double> &todTest_arr, levels::arr<double> &todTest_arr2);

    void arrFillingcm_v2( long ntod, levels::arr<double> &timeTest_arr, levels::arr<double> &outpntarrx, levels::arr<double> &outpntarr, long offindex );

    void report_timing();
    void timing_line(std::string label, double timer);
};

// wignercalc.cpp

double xpow( int expo, double val );
void wignerCalc( tsize n, tsize mmax, double theta, levels::arr2<double> &d );
void wignerCalcGeneral( tsize n, tsize mmax, double theta, levels::arr2<double> &d );
void wignerCalcHalfpi( tsize n, tsize mmax, levels::arr2<double> &d );
double wignerCalc00_halfpi( tsize n );
double wignerCalc00( double theta, tsize n );

// conviqt_util.cpp

void sift_down_DL( levels::arr<double> &ra, levels::arr<long> &brr, const int l, const int r );
void hpsort_DL( levels::arr<double> &ra, levels::arr<long> &brr );
void sift_down_arrTheta( levels::arr<double> &ra, const int l, const int r );
void hpsort_arrTheta( levels::arr<double> &ra );
void sift_down_arrTOD( levels::arr<double> &ra, const int l, const int r );
void hpsort_arrTOD( levels::arr<double> &ra );
void sift_down_arrTime( levels::arr<double> &ra, const int l, const int r );
void hpsort_arrTime( levels::arr<double> &ra );
void sift_down_DDcm( levels::arr<double> &ra, levels::arr<double> &brr, const int l, const int r );
void hpsort_DDcm( levels::arr<double> &ra, levels::arr<double> &brr );

} // namespace conviqt

#endif
