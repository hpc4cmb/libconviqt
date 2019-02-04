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

// The declarations go here

#define __CONVIQT_HPP__

using namespace levels;

namespace conviqt {

class beam {
public :
    beam() {}
    beam(long beamlmax, long beammmax, bool pol,
         std::string infile_beam, MPI_Comm comm);
    int read(long beamlmax, long beammmax, bool pol,
             std::string infile_beam,
             MPI_Comm comm=MPI_COMM_WORLD);
    Alm< xcomplex<float> > & blmT(void);
    Alm< xcomplex<float> > & blmG(void);
    Alm< xcomplex<float> > & blmC(void);
    int get_lmax(void) { return lmax; }
    int get_mmax(void) { return mmax; }
    double normalize(void);
private :
    Alm< xcomplex<float> > blmT_, blmG_, blmC_;
    long lmax, mmax;
    bool pol;
    int verbosity = 0;
    std::string fname;
};

class sky {
public :
    sky() {}
    sky(long skylmax, bool pol,
        std::string infile_sky,
        double fwhm_deconv_sky,
        MPI_Comm comm);
    int read(long skylmax, bool pol,
             std::string infile_sky,
             double fwhm_deconv_sky=0,
             MPI_Comm comm=MPI_COMM_WORLD);
    Alm< xcomplex<float> > & slmT(void);
    Alm< xcomplex<float> > & slmG(void);
    Alm< xcomplex<float> > & slmC(void);
    int get_lmax(void) { return lmax; }
    void remove_monopole(void);
    void remove_dipole(void);
private :
    Alm< xcomplex<float> > slmT_, slmG_, slmC_;
    long lmax;
    bool pol;
    double fwhm_deconv;
    int verbosity = 0;
    std::string fname;
};

class pointing : public levels::arr<double> {
    // The storage is 1D but there are 5 virtual columns:
    // row * 5 + 0 = longitude (phi)
    // row * 5 + 1 = co-latitude (theta)
    // row * 5 + 2 = position angle (psi)
    // row * 5 + 3 = signal
    // row * 5 * 4 = time
public :

private :

};


class detector {

public :

    detector() {};
    detector(std::string det_id) : det_id_(det_id) {};
    void set_epsilon(double epsilon) { epsilon_ = epsilon; };
    double get_epsilon() {return epsilon_;};
    void set_id(std::string det_id) { det_id_ = det_id; };
    std::string get_id() { return det_id_; };

private :

    std::string det_id_;
    double epsilon_;
    int verbosity = 0;

};

class convolver {

public :

    convolver(sky *s, beam *b, detector *d, bool pol=false,
              long lmax=5000, long beammax=14,
              long order=5, int verbosity=1,
              MPI_Comm comm=MPI_COMM_WORLD);
    int convolve(pointing & pnt, bool calibrate=true);
    int set_sky(sky *s);
    int set_beam(beam *b);
    int set_detector(detector *d);

private :

    int verbosity = 0;

    // Variables related to the gridding of the data cube

    double phi0, dphi, inv_delta_phi, phioffset, dtheta, theta0, inv_delta_theta;
    long npsi, nphi, halfmargin, margin, ntheta;
    int npoints, ioffset;

    long beta_to_itheta(double beta);
    double itheta_to_beta1(const long itheta);
    double itheta_to_beta2(const long itheta);

    void distribute_colatitudes(levels::arr<double> &pntarr, const long totsize,
                                levels::arr<double> &corethetaarr);
    MPI_Manager mpiMgr;
    int corenum, cores;
    sky *s;
    beam *b;
    detector *d;
    bool pol;
    long lmax, beammmax, order;
    std::vector<double> base_wgt;

    double t_convolve, t_todgen, t_arrFillingcm, t_interpolTOD,
        t_itheta0SetUp, t_ithetacalc, t_conviqt_tod_loop,
        t_conviqt_hemiscm, t_todAnnulus,
        t_wigner_init, t_wigner_prepare,
        t_alltoall, t_todRedistribution5cm,
        t_distribute_colatitudes, t_conviqt_hemiscm_pol,
        t_alltoall_datacube, t_sort, t_fillingBetaSeg;

    long n_convolve, n_todgen, n_arrFillingcm, n_interpolTOD,
        n_itheta0SetUp, n_ithetacalc, n_conviqt_tod_loop,
        n_conviqt_hemiscm, n_todAnnulus,
        n_wigner_init, n_wigner_prepare,
        n_alltoall, n_todRedistribution5cm,
        n_distribute_colatitudes, n_conviqt_hemiscm_pol,
        n_alltoall_datacube, n_sort, n_fillingBetaSeg;

    void weight_ncm_initialize();
    void weight_ncm(double x, std::vector<double> &wgt);

    void get_latitude_tables(const long NThetaIndex1,
                             const long NThetaIndex2,
                             const int ithetaoffset1,
                             const int ithetaoffset2,
                             std::vector<char> &need_itheta1_core,
                             std::vector<char> &need_itheta2_core,
                             std::vector<long> &my_itheta,
                             std::vector<int> &itheta_core);

    void conviqt_hemiscm(levels::arr3< xcomplex<double> > &tod1,
                         levels::arr3< xcomplex<double> > &tod2,
                         const long NThetaIndex1,
                         const long NThetaIndex2,
                         const int ithetaoffset1,
                         const int ithetaoffset2);

    void conviqt_hemiscm_pol(levels::arr3< xcomplex<double> > &tod1,
                             levels::arr3< xcomplex<double> > &tod2,
                             const long NThetaIndex1,
                             const long NThetaIndex2,
                             const int ithetaoffset1,
                             const int ithetaoffset2);

    void alltoall_datacube(std::vector<char> &need_itheta_core,
                           std::vector<int> &itheta_core,
                           std::vector<long> &my_itheta,
                           levels::arr3< xcomplex<double> > &my_tod,
                           levels::arr3< xcomplex<double> > &tod,
                           const long ithetaoffset, const long NThetaIndex);

    void todAnnulus(levels::arr3< xcomplex<double> > &tod,
                    levels::arr3< xcomplex<double> > &Cmm,
                    levels::arr<double> &cs,
                    levels::arr<double> &sn,
                    levels::arr<double> &cs0,
                    levels::arr<double> &sn0,
                    const long NThetaIndex);

    void conviqt_tod_loop(levels::arr<long> &lowerIndex,
                          levels::arr<long> &upperIndex,
                          levels::arr<double> &outpntarr,
                          levels::arr3< xcomplex<double> > &TODAsym,
                          const long thetaIndex,
                          levels::arr<long> &itheta0,
                          const long ntod,
                          levels::arr<double> &TODValue,
                          const long lat);

    void ratiobetacalcgreatercm(int &corenum, double theta,
                                levels::arr<double> &corethetaarr);
    void ratiobetacalcsmallercm(int &corenum, double theta,
                                levels::arr<double> &corethetaarr);

    void interpolTOD(levels::arr<double> &outpntarr1,
                     levels::arr<double> &outpntarr2,
                     levels::arr<double> &TODValue1,
                     levels::arr<double> &TODValue2,
                     const long ntod1, const long ntod2);

    void itheta0SetUp(levels::arr<double> outpntarr,
                      const long ntod,
                      long &NThetaIndex,
                      levels::arr<long> &itheta0,
                      levels::arr<long> &lowerIndex,
                      levels::arr<long> &upperIndex,
                      levels::arr3<xcomplex<double> > &TODAsym);

    void ithetacalc(levels::arr<long> &itheta0,
                    levels::arr<double> &outpntarr,
                    const long ntod);

    void fillingBetaSeg(levels::arr<double> & pntarr, const long arrsize,
                        const double ratiodeltas,
                        levels::arr<double> &corethetaarr, levels::arr<int> &inBetaSeg);

    void todRedistribution5cm(levels::arr<double> pntarr, levels::arr<int> inBetaSeg,
                              levels::arr<int> outBetaSeg, levels::arr<int> &inBetaSegAcc,
                              levels::arr<int> &outBetaSegAcc, long &outBetaSegSize,
                              levels::arr<double> &pntarr2, const long totsize,
                              levels::arr<int> &inOffset, levels::arr<int> &outOffset,
                              const double ratiodeltas, levels::arr<double> &corethetaarr);

    void todgen(const long ntod1,
                const long ntod2,
                levels::arr<double> &outpntarr);

    void arrFillingcm(long ntod,
                      levels::arr<double> &outpntarrx,
                      levels::arr<double> &outpntarr,
                      long offindex);

    void report_timing();
    void timing_line(std::string label, double timer, long counter);
};

// wignercalc.cpp

double xpow(int expo, double val);
void wignerCalc(tsize n, tsize mmax, double theta, levels::arr2<double> &d);
void wignerCalcGeneral(tsize n, tsize mmax, double theta, levels::arr2<double> &d);
void wignerCalcHalfpi(tsize n, tsize mmax, levels::arr2<double> &d);
double wignerCalc00_halfpi(tsize n);
double wignerCalc00(double theta, tsize n);

// conviqt_util.cpp

void sift_down_DL(levels::arr<double> &ra, levels::arr<long> &brr, const int l, const int r);
void hpsort_DL(levels::arr<double> &ra, levels::arr<long> &brr);
void sift_down_arrTheta(levels::arr<double> &ra, const int l, const int r);
void hpsort_arrTheta(levels::arr<double> &ra);
void sift_down_arrTOD(levels::arr<double> &ra, const int l, const int r);
void hpsort_arrTOD(levels::arr<double> &ra );
void sift_down_arrTime(levels::arr<double> &ra, const int l, const int r);
void hpsort_arrTime(levels::arr<double> &ra);
void sift_down_DDcm(levels::arr<double> &ra, levels::arr<double> &brr, const int l, const int r);
void hpsort_DDcm(levels::arr<double> &ra, levels::arr<double> &brr);

} // namespace conviqt

#endif
