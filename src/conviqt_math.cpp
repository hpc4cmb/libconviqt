#include "conviqt.hpp"

#include <sstream>

#include <cstring>

#include "omp.h"

/*
  int convolve()
      void fillingBetaSeg()
          void ratiobetacalcsmallercm()
          void ratiobetacalcgreatercm()
      void todRedistribution5cm()
          void ratiobetacalcsmallercm()
          void ratiobetacalcgreatercm()
      mpiMgr.all2allv()
      hpsort_arrTheta()

      void todgen()
          void arrFillingcm()
          void interpolTOD()
              void itheta0SetUp()
                  void ithetacalc()
              EITHER
              void conviqt_hemiscm()
              OR
              void conviqt_hemiscm_pol()
                  get_latitude_tables()
                  wignergen()
                  wignergen.prepare()
                  wignergen.calc()
                  void todAnnulus()
                  alltoall_datacube()
              void conviqt_tod_loop()
                  void weight_ncm()

      hpsort_DL()
      hpsort_arrTOD()
      mpiMgr.all2allv()
      hpsort_arrTime
      hpsort_DDcm()
      mpiMgr.all2allv()
      hpsort_DDcm()

  int set_sky()
  int set_beam()
  int set_detector()

 */

// This file will contain the actual operations on the sky, no I/O

namespace conviqt {

beam::beam(long beamlmax, long beammmax, bool pol,
           std::string infile_beam,
           MPI_Comm comm=MPI_COMM_WORLD) {
    read(beamlmax, beammmax, pol, infile_beam, comm);
}

double beam::normalize() {
    /*
      Scale the beam expansion to integrate to 0.5 over the sphere
    */
    double scale = 0;
    if (blmT_.Lmax() >= 0) {
        double b00 = blmT_(0, 0).re;
        double current_norm = 2 * b00 / sqrt(1 / pi);
        scale = 0.5 / current_norm;
        if (CMULT_VERBOSITY > 1) {
            std::cerr << "Normalizing beam from " << current_norm << " to 0.5 with "
                      << scale << std::endl;
        }
        blmT_.Scale(scale);
        if (pol) {
            blmG_.Scale(scale);
            blmC_.Scale(scale);
        }
    }
    return scale;
}

sky::sky(long skylmax, bool pol, std::string infile_sky, double fwhm_deconv_sky=0,
         MPI_Comm comm=MPI_COMM_WORLD) {
    read(skylmax, pol, infile_sky, fwhm_deconv_sky, comm);
}

void sky::remove_monopole(void) {
  if (slmT_.Lmax() >= 0) {
    slmT_(0, 0) = 0;
  }
}

void sky::remove_dipole(void) {
  if (slmT_.Lmax() >= 1) {
    slmT_(1, 0) = 0;
    if (slmT_.Mmax() >= 1) {
      slmT_(1, 1) = 0;
    }
  }
}


convolver::convolver(sky *s, beam *b, detector *d, bool pol,
                     long lmax, long beammmax,
                     long order, MPI_Comm comm)
    : s(s), b(b), d(d), pol(pol), lmax(lmax), beammmax(beammmax), order(order) {
    mpiMgr = MPI_Manager(comm);
    cores = mpiMgr.num_ranks();
    corenum = mpiMgr.rank();

    const int lmax_sky = s->get_lmax();
    const int lmax_beam = b->get_lmax();
    if (lmax < 0) {
        this->lmax = std::max(lmax_beam, lmax_sky);
    } else if (lmax > lmax_sky || lmax > lmax_beam) {
        throw std::runtime_error("Convolver lmax exceeds input expansion order.");
    }

    const int mmax_beam = b->get_mmax();
    if (beammmax < 0) {
        this->beammmax = mmax_beam;
    } else if (beammmax > mmax_beam) {
        throw std::runtime_error("Convolver mmax exceeds input expansion order.");
    }

    // data cube gridding

    npoints = order + 1;
    ioffset = order / 2;
    phi0 = halfpi;
    npsi = this->beammmax + 1;
    nphi = 2 * this->lmax + 1;
    dphi = 2 * pi / nphi;
    inv_delta_phi = 1. / dphi;
    phioffset = phi0 / dphi;
    halfmargin = npoints / 2 + 1;
    margin = 2 * halfmargin;
    ntheta = this->lmax + 1 + margin;
    dtheta = -pi / (ntheta - margin);
    theta0 = pi - halfmargin * dtheta;
    inv_delta_theta = 1 / dtheta;

    // profiling

    t_convolve = 0;
    n_convolve = 0;
    t_alltoall = 0;
    n_alltoall = 0;
    t_sort = 0;
    n_sort = 0;
    t_fillingBetaSeg = 0;
    n_fillingBetaSeg = 0;
    t_todRedistribution5cm = 0;
    n_todRedistribution5cm = 0;
    t_distribute_colatitudes = 0;
    n_distribute_colatitudes = 0;

    t_todgen = 0;
    t_arrFillingcm = 0;
    t_interpolTOD = 0;
    t_itheta0SetUp = 0;
    t_ithetacalc = 0;
    t_conviqt_tod_loop = 0;
    t_conviqt_hemiscm = 0;
    t_todAnnulus = 0;
    t_conviqt_hemiscm_pol = 0;
    t_alltoall_datacube = 0;

    n_todgen = 0;
    n_arrFillingcm = 0;
    n_interpolTOD = 0;
    n_itheta0SetUp = 0;
    n_ithetacalc = 0;
    n_conviqt_tod_loop = 0;
    n_conviqt_hemiscm = 0;
    n_todAnnulus = 0;
    n_conviqt_hemiscm_pol = 0;
    n_alltoall_datacube = 0;

    t_wigner_init = 0;
    t_wigner_prepare = 0;

    n_wigner_init = 0;
    n_wigner_prepare = 0;
}


int convolver::set_sky(sky *s) {

  this->s = s;

  return 0;
}

int convolver::set_beam(beam *b) {

  this->b = b;

  return 0;
}

int convolver::set_detector(detector *d) {

  this->d = d;

  return 0;
}


void convolver::weight_ncm_initialize() {
    /*
      Allocate and initialize base_wgt to be used in weight_ncm
     */
    if (base_wgt.size() == 0) {
        // Initialize base_wgt only when needed
        base_wgt.resize(npoints, 1);
        for (int m = 0; m < npoints; ++m) {
            for (int n = 0; n < m; ++n) {
                base_wgt[m] *= m - n;
            }
            for (int n = m + 1; n < npoints; ++n) {
                base_wgt[m] *= m - n;
            }
        }
        for (int m = 0; m < npoints; ++m) {
            base_wgt[m] = 1. / base_wgt[m];
        }
    } else if (base_wgt.size() != npoints) {
        throw std::runtime_error("base_wgt changed dimensions");
    }
}


void convolver::weight_ncm(const double x, std::vector<double> &wgt) {
    /*
      Return the interpolation weights used in conviqt_tod_loop

      This method is called from multiple threads at once
    */
    memcpy(&(wgt[0]), &(base_wgt[0]), sizeof(double) * npoints);

    double mul1 = x;
    for (int m = 1; m < npoints; ++m) {
        wgt[m] *= mul1;
        mul1 *= x - m;
    }

    double mul2 = x - npoints + 1;
    const long loffset = npoints - 1;
    const double doffset = x - npoints + 1;
    for (int m = 1; m < npoints; ++m) {
        wgt[loffset - m] *= mul2;
        mul2 *= doffset + m;
    }
}


void convolver::distribute_colatitudes(levels::arr<double> &pntarr,
                                       const long totsize,
                                       levels::arr<double> &corethetaarr) {
    /*
      Distribute the latitudes between 0..90 degrees between
      processes. 90..180 deg follow from symmetry.  The lower limits are
      stored in corethetaarr.
     */

    double tstart = mpiMgr.Wtime();
    ++n_distribute_colatitudes;

    long totsize_tot = 0;
    mpiMgr.allreduceRaw<long>(&totsize, &totsize_tot, 1, mpiMgr.Sum);

    // Build a latitude hit map

    long nbin = 100000;
    double wbin = halfpi / nbin;
    arr<int> my_hits(nbin, 0), hits(nbin, 0);
    for (long ii = 0; ii < totsize; ++ii) {
        double beta = pntarr[5 * ii + 1];
        if (beta > halfpi) {
            beta = pi - beta;
        }
        long ibin = beta / wbin;
        ++my_hits[ibin];
    }

    mpiMgr.allreduce(my_hits, hits, mpiMgr.Sum);

    long nhitleft = totsize_tot;

    corethetaarr[0] = 0;
    long ibin = -1;
    for (int corenum = 1; corenum < cores; ++corenum) {
        // Width of the bin is based on remaining latitude
        // and samples, which ever we meet first
        double divisor = cores - corenum + 1;
        double nhit_target = nhitleft / divisor;
        long bin_hit = 0;
        while (ibin < nbin) {
            ++ibin;
            // Skip over isolatitudes that have no hits
            if (hits[ibin] != 0) {
                bin_hit += hits[ibin];
                nhitleft -= hits[ibin];
                if (bin_hit >= nhit_target) {
                    break;
                }
            }
        }
        double theta = ibin * wbin;
        corethetaarr[corenum] = theta;
    }

    t_distribute_colatitudes += mpiMgr.Wtime() - tstart;
}


void convolver::fillingBetaSeg(levels::arr<double> &pntarr,
                               const long arrsize,
                               const double ratiodeltas,
                               levels::arr<double> &corethetaarr,
                               levels::arr<int> &inBetaSeg) {
    /*
      This method determines the number of doubles to send from this
      process to all other processes.  The counts are recorded in inBetaSeg
    */

    double tstart = mpiMgr.Wtime();
    ++n_fillingBetaSeg;
    if (CMULT_VERBOSITY > 1) {
        std::cerr << corenum << " : Entering fillingBetaSeg" << std::endl;
    }

    inBetaSeg.alloc(cores);
    inBetaSeg.fill(0);

    int corenum = 0;
    for (long ii = 0; ii < arrsize; ++ii) {
        double theta = pntarr[5 * ii + 1];
        if (theta >= halfpi) {
            theta = pi - theta;
        }
        // Initialize corenum to a rough estimate of which task
        // should own this colatitude ..
        // corenum = theta / ratiodeltas;
        // .. and then perform a linear search to find the process
        // that owns it:
        //   corethetaarr[corenum] < theta < corethetaarr[corenum + 1]
        if (theta >= corethetaarr[corenum]) {
            // increment corenum
            ratiobetacalcgreatercm(corenum, theta, corethetaarr);
        } else {
            // decrement corenum
            ratiobetacalcsmallercm(corenum, theta, corethetaarr);
        }
        inBetaSeg[corenum] += 5; // Each entry comprises 5 double elements
    }

    if (CMULT_VERBOSITY > 1) {
        std::cerr << corenum << " : Leaving fillingBetaSeg" << std::endl;
    }

    t_fillingBetaSeg += mpiMgr.Wtime() - tstart;
}


void convolver::ratiobetacalcgreatercm(int &corenum, const double theta,
                                       levels::arr<double> &corethetaarr) {
    while (corenum < cores - 1) {
        if (theta < corethetaarr[corenum + 1]) {
            break;
        }
        ++corenum;
    }
}


void convolver::ratiobetacalcsmallercm(int &corenum, const double theta,
                                       levels::arr<double> &corethetaarr) {
    while (corenum > 0) {
        if (theta > corethetaarr[corenum]) {
            break;
        }
        --corenum;
    }
}

void convolver::todRedistribution5cm(levels::arr<double> pntarr,
                                     levels::arr<int> inBetaSeg,
                                     levels::arr<int> outBetaSeg,
                                     levels::arr<int> &inBetaSegAcc,
                                     levels::arr<int> &outBetaSegAcc,
                                     long &outBetaSegSize,
                                     levels::arr<double> &pntarr2,
                                     const long totsize,
                                     levels::arr<int> &inOffset,
                                     levels::arr<int> &outOffset,
                                     const double ratiodeltas,
                                     levels::arr<double> &corethetaarr) {
    /*
      Establish offsets for alltoallv and pack the send buffer.
     */
    double tstart = mpiMgr.Wtime();
    ++n_todRedistribution5cm;
    if (CMULT_VERBOSITY > 1) {
        std::cerr << "Entered todRedistribution5cm" << std::endl;
    }

    inBetaSegAcc.alloc(cores);
    inBetaSegAcc.fill(0);
    outBetaSegAcc.alloc(cores);
    outBetaSegAcc.fill(0);
    for (long ii = 0; ii < cores; ++ii) {
        for (long jj = 0; jj <= ii; ++jj) {
            inBetaSegAcc[ii] += inBetaSeg[jj];
            outBetaSegAcc[ii] += outBetaSeg[jj];
        }
    }
    outBetaSegSize = outBetaSegAcc[cores - 1];
    inOffset.alloc(cores);
    inOffset.fill(0);
    outOffset.alloc(cores);
    outOffset.fill(0);
    for (long ii = 1; ii < cores; ++ii) {
        outOffset[ii] = outBetaSegAcc[ii - 1];
        inOffset[ii] = inBetaSegAcc[ii - 1];
    }
    if (totsize > 0) {
        try {
            pntarr2.alloc(5 * totsize);
        } catch (std::bad_alloc &e) {
            std::cerr << " todRedistribution5cm : Out of memory allocating "
                      << (5 * totsize) * 8. / 1024 / 1024 << "MB for pntarr2"
                      << std::endl;
            throw;
        }
        // Pack the data in pntarr into pntarr2 for sending
        pntarr2.fill(0);
        long offset1 = 0;
        levels::arr<int> offset2(inOffset);
        int corenum = 0;
        for (long ii = 0; ii < totsize; ++ii) {
            double theta = pntarr[5 * ii + 1];
            if (theta < 0 || theta > pi) {
                throw std::runtime_error("ERROR: Illegal latitude in pntarr: theta = " +
                                         std::to_string(theta) + " not in 0..pi");
            }
            if (theta >= halfpi) {
                theta = pi - theta;
            }
            // corenum = theta / ratiodeltas;
            if (theta >= corethetaarr[corenum]) {
                ratiobetacalcgreatercm(corenum, theta, corethetaarr);
            } else {
                ratiobetacalcsmallercm(corenum, theta, corethetaarr);
            }
            memcpy(&(pntarr2[offset2[corenum]]), &(pntarr[offset1]), 5 * sizeof(double));
            offset1 += 5;
            offset2[corenum] += 5;
	}
    }

    if (CMULT_VERBOSITY > 1) {
        std::cerr << "Leaving todRedistribution5cm" << std::endl;
    }

    t_todRedistribution5cm += mpiMgr.Wtime() - tstart;
}


void convolver::todgen(const long ntod1, const long ntod2,
                       levels::arr<double> &outpntarr) {
    double tstart = mpiMgr.Wtime();
    ++n_todgen;

    levels::arr<double> todTest_arr, todTest_arr2;
    levels::arr<double> outpntarr1, outpntarr2;

    // Store the input ordering of outpntarr in the signal column
    for (long ii = 0; ii < ntod1 + ntod2; ++ii) {
        outpntarr[5 * ii + 3] = ii;
    }

    if (ntod1 != 0) {
        // Collect northern hemisphere data from outpntarr to outpntarr1
        arrFillingcm(ntod1, outpntarr1, outpntarr, 0);
        todTest_arr.alloc(ntod1);
        todTest_arr.fill(0);
    }
    if (ntod2 != 0) {
        // Collect southern hemisphere data from outpntarr to outpntarr2
        arrFillingcm(ntod2, outpntarr2, outpntarr, ntod1);
        todTest_arr2.alloc(ntod2);
        todTest_arr2.fill(0);
    }

    if (CMULT_VERBOSITY > 1) {
        std::cerr << "corenum = " << corenum << "  order = " << order
                  << "  ntod1 = " << ntod1 << "  ntod2 = " << ntod2 << std::endl;
    }

    mpiMgr.barrier();
    interpolTOD(outpntarr1, outpntarr2,
                todTest_arr, todTest_arr2,
                ntod1, ntod2);

    // Copy signal from todTest_arr into outpntarr

    for (long ii = 0; ii < ntod1; ++ii) {
        long old_row = outpntarr1[5 * ii + 3] + .5;
        outpntarr[5 * old_row + 3] = todTest_arr[ii];
    }
    todTest_arr.dealloc();
    outpntarr1.dealloc();

    for (long ii = 0; ii < ntod2; ++ii) {
        long old_row = outpntarr2[5 * ii + 3] + .5;
        outpntarr[5 * old_row + 3] = todTest_arr2[ii];
    }
    todTest_arr2.dealloc();
    outpntarr2.dealloc();

    t_todgen += mpiMgr.Wtime() - tstart;
}


void convolver::itheta0SetUp(levels::arr<double> outpntarr,
                             const long ntod,
                             long &NThetaIndex,
                             levels::arr<long> &itheta0,
                             levels::arr<long> &lowerIndex,
                             levels::arr<long> &upperIndex,
                             levels::arr3<xcomplex<double> > &TODAsym) {
    // Map each TOD sample onto isolatitude rings (itheta0)
    // and find consecutive ranges of samples that share a ring
    // (lowerIndex and upperIndex)
    double tstart = mpiMgr.Wtime();
    ++n_itheta0SetUp;
    itheta0.alloc(ntod);
    ithetacalc(itheta0, outpntarr, ntod);
    NThetaIndex = itheta0[0] - itheta0[ntod - 1] + npoints;
    lowerIndex.alloc(NThetaIndex);
    lowerIndex.fill(ntod);
    upperIndex.alloc(NThetaIndex);
    upperIndex.fill(0);
    for (long jj = 0; jj < NThetaIndex; ++jj) {
        for (long ii = ntod - 1; ii >= 0; --ii) {
            if (itheta0[ii] > itheta0[ntod - 1] + jj - npoints) {
                lowerIndex[jj] = ii;
                break;
            }
        }
    }

    for (long jj = 0; jj < NThetaIndex; ++jj) {
        for (long ii = 0; ii < ntod; ++ii) {
            if (itheta0[ii] <= itheta0[ntod - 1] + jj) {
                // this yields the smallest ii where
                // itheta0[ii] <= itheta0[ntod - 1] + jj
                upperIndex[jj] = ii;
                break;
            }
        }
    }

    try {
        TODAsym.alloc(nphi, npsi, NThetaIndex);
    } catch (std::bad_alloc &e) {
        std::cerr << "itheta0SetUp :  Out of memory allocating "
                  << nphi * npsi * NThetaIndex * 8. / 1024 / 1024
                  << "MB for TODAsym" << std::endl;
        throw;
    }

    TODAsym.fill(0.);
    t_itheta0SetUp += mpiMgr.Wtime() - tstart;
}


long convolver::beta_to_itheta(double beta) {
    /*
      Translate co-latitude into latitude index
    */
    if (beta > halfpi) {
        beta = pi - beta;
    }
    // Note that the larger the ii the smaller frac is and the smaller itheta0[ii] is.
    double frac = (beta - theta0) * inv_delta_theta;
    // Note also that itheta0[ii] is always positive.
    long itheta = int(frac) - ioffset;
    if (itheta > (ntheta - npoints)) {
        itheta = ntheta - npoints;
    }
    if (itheta < 0) {
        itheta = 0;
    }
    return itheta;
}


double convolver::itheta_to_beta1(const long itheta) {
    /*
      Translate latitude index into co-latitude (Northern hemisphere)
    */
    return theta0 + itheta * dtheta;
}


double convolver::itheta_to_beta2(const long itheta) {
    /*
      Translate latitude index into co-latitude (Southern hemisphere)
    */
    return pi - (theta0 + itheta * dtheta);
}


void convolver::ithetacalc(levels::arr<long> &itheta0,
                           levels::arr<double> &outpntarr,
                           const long ntod) {
    // Map each TOD sample to an isolatitude ring
    // The TOD is already ordered by latitude
    double tstart = mpiMgr.Wtime();
    ++n_ithetacalc;
    for (long ii = 0; ii < ntod; ++ii) {
        double beta = outpntarr[5 * ii + 1]; // co-latitude
        long itheta = beta_to_itheta(beta);
        itheta0[ii] = itheta;
    }
    t_ithetacalc += mpiMgr.Wtime() - tstart;
}


void convolver::conviqt_hemiscm(levels::arr3<xcomplex<double> > &tod1,
                                levels::arr3<xcomplex<double> > &tod2,
                                const long NThetaIndex1,
                                const long NThetaIndex2,
                                const int ithetaoffset1,
                                const int ithetaoffset2) {
    /*
      Build the beam-convolved theta/phi/psi data cube.

      This unpolarized version duplicates most of conviqt_hemiscm_pol
      but avoids conditionals in performance-critical code.
     */
    if (CMULT_VERBOSITY > 1) {
        mpiMgr.barrier();
        std::cerr << corenum << " : Entering conviqt_hemiscm" << std::endl;
    }
    double tstart = mpiMgr.Wtime();
    ++n_conviqt_hemiscm;

    std::vector<char> need_itheta1_core;
    std::vector<char> need_itheta2_core;
    std::vector<int> itheta_core;  // Cores assigned to each latitude
    std::vector<long> my_itheta;  // Locally assigned latitudes

    get_latitude_tables(NThetaIndex1, NThetaIndex2, ithetaoffset1, ithetaoffset2,
                        need_itheta1_core, need_itheta2_core, my_itheta, itheta_core);

    long my_ntheta = my_itheta.size();

    levels::arr<double> rthetas(my_ntheta);
    for (long itheta_local = 0; itheta_local < my_ntheta; ++itheta_local) {
        long itheta_global = my_itheta[itheta_local];
        rthetas[itheta_local] = itheta_to_beta1(itheta_global);
    }

    // Locally owned part of the data cube

    levels::arr3< xcomplex<double> > my_Cmm1, my_Cmm2, my_tod1, my_tod2;

    try {
        my_Cmm1.alloc(nphi, npsi, my_ntheta);
        my_Cmm2.alloc(nphi, npsi, my_ntheta);
        my_tod1.alloc(nphi, npsi, my_ntheta);
        my_tod2.alloc(nphi, npsi, my_ntheta);
    } catch (std::bad_alloc &e) {
        std::cerr << "conviqt_hemiscm_pol :  Out of memory allocating "
                  << nphi * npsi * my_ntheta * 4 * 2 * 8. / 1024 / 1024
                  << "MB for locally owned data cube" << std::endl;
        throw;
    }

    my_Cmm1.fill(0.);
    my_Cmm2.fill(0.);
    my_tod1.fill(0.);
    my_tod2.fill(0.);

    Alm< xcomplex<float> > &blmT = b->blmT();
    Alm< xcomplex<float> > &slmT = s->slmT();
    levels::arr<double> cs(nphi), sn(nphi);
    levels::arr<double> cs0(nphi), sn0(nphi);

#pragma omp parallel default(shared)
    {
        double t1 = mpiMgr.Wtime();
        wignergen wgen(lmax, rthetas, conv_acc), wgen_neg(lmax, rthetas, conv_acc);
        if (omp_get_thread_num() == 0) {
            t_wigner_init += mpiMgr.Wtime() - t1;
            ++n_wigner_init;
        }
        //wigner_estimator estimator(lmax,100);
        for(long mbeam = 0; mbeam < npsi; ++mbeam) {
            double dsignb = levels::xpow(mbeam, 1);
#pragma omp for schedule(static, 8)
            for (long msky = 0; msky <= lmax; ++msky) {
                const double dsign = levels::xpow(msky, 1);
                const double dsb = dsign * dsignb;
                // estimator.prepare_m(mbeam, msky);
                // if (estimator.canSkip(rthetas[my_ntheta - 1]))
                //     continue; // negligible dmm
                t1 = mpiMgr.Wtime();
                wgen.prepare(mbeam, msky);
                wgen_neg.prepare(mbeam, -msky);
                if (omp_get_thread_num() == 0) {
                    t_wigner_prepare += mpiMgr.Wtime() - t1;
                    ++n_wigner_prepare;
                }
                for (long lat = 0; lat < my_ntheta; ++lat) {
                    xcomplex<double> &Cmm1_pos = my_Cmm1(msky + lmax, mbeam, lat);
                    xcomplex<double> &Cmm1_neg = my_Cmm1(-msky + lmax, mbeam, lat);
                    xcomplex<double> &Cmm2_pos = my_Cmm2(msky + lmax, mbeam, lat);
                    xcomplex<double> &Cmm2_neg = my_Cmm2(-msky + lmax, mbeam, lat);
                    int firstl1, firstl2;
                    const levels::arr<double> &dmm = wgen.calc(lat, firstl1);
                    const levels::arr<double> &dmmneg = wgen_neg.calc(lat, firstl2);
                    const int firstl = std::max(firstl1, firstl2);
                    double dlb = -levels::xpow(firstl, dsignb);
                    for (long ii = firstl; ii <= lmax; ++ii) {
                        dlb = -dlb;
                        // Note that msky in dlm is located to the left of
                        // mb and they have to be interchanged in the convolution
                        const double dMatrixElementmskypos = dsb * dmm[ii];
                        const double dMatrixElementmskyneg = dsb * dmmneg[ii];
                        const double dMatrixElementmskypos2 = dlb * dMatrixElementmskyneg;
                        const double dMatrixElementmskyneg2 = dlb * dMatrixElementmskypos;
                        const xcomplex<float> sT = slmT(ii, msky);
                        const xcomplex<float> bT = blmT(ii, mbeam);
                        const double prod1 = sT.re * bT.re;
                        const double prod3 = sT.im * bT.re;
                        const double prod2 = sT.im * bT.im;
                        const double prod4 = sT.re * bT.im;
                        const double tmp_1 = prod1 + prod2;
                        const double tmp_2 = prod3 - prod4;
                        const double xtmp_1 = tmp_1 * dMatrixElementmskypos;
                        const double xtmp_2 = tmp_2 * dMatrixElementmskypos;
                        const double xtmp_5 = tmp_1 * dMatrixElementmskypos2;
                        const double xtmp_6 = tmp_2 * dMatrixElementmskypos2;
                        Cmm1_pos.re += xtmp_1;
                        Cmm1_pos.im += xtmp_2;
                        Cmm2_pos.re += xtmp_5;
                        Cmm2_pos.im += xtmp_6;
                        if (msky != 0) {
                            const double tmp_3 = dsign * (prod1 - prod2);
                            const double tmp_4 = -dsign * (prod3 + prod4);
                            const double xtmp_3 = tmp_3 * dMatrixElementmskyneg;
                            const double xtmp_4 = tmp_4 * dMatrixElementmskyneg;
                            const double xtmp_7 = tmp_3 * dMatrixElementmskyneg2;
                            const double xtmp_8 = tmp_4 * dMatrixElementmskyneg2;
                            Cmm1_neg.re += xtmp_3;
                            Cmm1_neg.im += xtmp_4;
                            Cmm2_neg.re += xtmp_7;
                            Cmm2_neg.im += xtmp_8;
                        }
                    }
                }
            }
        }

        const double lmaxinv = double(2 * pi * lmax) / double(2 * lmax + 1);
#pragma omp for schedule(static, 8)
        for (long msky = -lmax; msky <= lmax; ++msky) {
            double arg = -halfpi * msky;
            long ii = lmax + msky;
            cs[ii] = cos(arg);
            sn[ii] = sin(arg);
            arg = ii * lmaxinv;
            cs0[ii] = cos(arg);
            sn0[ii] = sin(arg);
        }
    } // End of parallel section

    todAnnulus(my_tod1, my_Cmm1, cs, sn, cs0, sn0, my_ntheta);
    todAnnulus(my_tod2, my_Cmm2, cs, sn, cs0, sn0, my_ntheta);

    alltoall_datacube(need_itheta1_core, itheta_core, my_itheta, my_tod1, tod1,
                      ithetaoffset1, NThetaIndex1);
    alltoall_datacube(need_itheta2_core, itheta_core, my_itheta, my_tod2, tod2,
                      ithetaoffset2, NThetaIndex2);

    t_conviqt_hemiscm += mpiMgr.Wtime() - tstart;

    if (CMULT_VERBOSITY > 1) {
        std::cerr << corenum << " : Exiting conviqt_hemiscm" << std::endl;
    }
}


void convolver::get_latitude_tables(const long NThetaIndex1,
                                    const long NThetaIndex2,
                                    const int ithetaoffset1,
                                    const int ithetaoffset2,
                                    std::vector<char> &need_itheta1_core,
                                    std::vector<char> &need_itheta2_core,
                                    std::vector<long> &my_itheta,
                                    std::vector<int> &itheta_core) {

    // Local tables of required latitudes

    std::vector<char> my_need_itheta1(ntheta, false);
    std::vector<char> my_need_itheta2(ntheta, false);

    for (int itheta = 0; itheta < NThetaIndex1; ++itheta) {
        my_need_itheta1[ithetaoffset1 + itheta] = true;
    }

    for (int itheta = 0; itheta < NThetaIndex2; ++itheta) {
        my_need_itheta2[ithetaoffset2 + itheta] = true;
    }

    // Global tables of required latitudes

    std::vector<char> need_itheta1(ntheta);
    std::vector<char> need_itheta2(ntheta);
    std::vector<char> need_itheta(ntheta);
    int err;

    err = MPI_Allreduce(my_need_itheta1.data(), need_itheta1.data(),
                        int(ntheta), MPI_CHAR, MPI_BOR, mpiMgr.comm());
    if (err) throw std::runtime_error("Failed to reduce need_itheta1");

    err = MPI_Allreduce(my_need_itheta2.data(), need_itheta2.data(),
                        int(ntheta), MPI_CHAR, MPI_BOR, mpiMgr.comm());
    if (err) throw std::runtime_error("Failed to reduce need_itheta2");

    need_itheta1_core.resize(ntheta * cores);
    err = MPI_Allgather(my_need_itheta1.data(), int(ntheta), MPI_CHAR,
                        need_itheta1_core.data(), int(ntheta), MPI_CHAR,
                        mpiMgr.comm());
    if (err) throw std::runtime_error("Failed to gather need_itheta1_core");

    need_itheta2_core.resize(ntheta * cores);
    err = MPI_Allgather(my_need_itheta2.data(), int(ntheta), MPI_CHAR,
                        need_itheta2_core.data(), int(ntheta), MPI_CHAR,
                        mpiMgr.comm());
    if (err) throw std::runtime_error("Failed to gather need_itheta2_core");

    // Collapse needed indices between North and South hemisphere

    long ntheta_need = 0;  // Total number of required latitudes
    for (long itheta = 0; itheta < ntheta; ++itheta) {
        if (need_itheta1[itheta] or need_itheta2[itheta]) {
            need_itheta[itheta] = true;
            ++ntheta_need;
        }
    }

    // Assign each latitude to exactly one core

    itheta_core.resize(ntheta, -1);
    // Number of latitudes per core
    long ntheta_core = ceil(ntheta_need / double(cores));
    long core = 0;
    long nassigned = 0;
    long nleft = ntheta_need;

    for (long itheta = 0; itheta < ntheta; ++itheta) {
        if (need_itheta[itheta]) {
            itheta_core[itheta] = core;
            if (core == corenum) {
                my_itheta.push_back(itheta);
            }
            ++nassigned;
            --nleft;
            if (nassigned == ntheta_core) {
                ++core;
                nassigned = 0;
                // Adjust latitudes per core to avoid cores not
                // having any latitudes
                ntheta_core = ceil(nleft / double(cores - core));
            }
        }
    }
}



void convolver::conviqt_hemiscm_pol(levels::arr3< xcomplex<double> > &tod1,
                                    levels::arr3< xcomplex<double> > &tod2,
                                    const long NThetaIndex1,
                                    const long NThetaIndex2,
                                    const int ithetaoffset1,
                                    const int ithetaoffset2) {
    /*
      Build the beam-convolved theta/phi/psi data cube
     */
    if (CMULT_VERBOSITY > 1) {
        mpiMgr.barrier();
        std::cerr << corenum << " : Entering conviqt_hemiscm_pol" << std::endl;
    }
    double tstart = mpiMgr.Wtime();
    ++n_conviqt_hemiscm_pol;

    std::vector<char> need_itheta1_core;
    std::vector<char> need_itheta2_core;
    std::vector<int> itheta_core;  // Cores assigned to each latitude
    std::vector<long> my_itheta;  // Locally assigned latitudes

    get_latitude_tables(NThetaIndex1, NThetaIndex2, ithetaoffset1, ithetaoffset2,
                        need_itheta1_core, need_itheta2_core, my_itheta, itheta_core);

    long my_ntheta = my_itheta.size();

    levels::arr<double> rthetas(my_ntheta);
    for (long itheta_local = 0; itheta_local < my_ntheta; ++itheta_local) {
        long itheta_global = my_itheta[itheta_local];
        rthetas[itheta_local] = itheta_to_beta1(itheta_global);
    }

    // Locally owned part of the data cube

    levels::arr3< xcomplex<double> > my_Cmm1, my_Cmm2, my_tod1, my_tod2;

    try {
        my_Cmm1.alloc(nphi, npsi, my_ntheta);
        my_Cmm2.alloc(nphi, npsi, my_ntheta);
        my_tod1.alloc(nphi, npsi, my_ntheta);
        my_tod2.alloc(nphi, npsi, my_ntheta);
    } catch (std::bad_alloc &e) {
        std::cerr << "conviqt_hemiscm_pol :  Out of memory allocating "
                  << nphi * npsi * my_ntheta * 4 * 2 * 8. / 1024 / 1024
                  << "MB for locally owned data cube" << std::endl;
        throw;
    }

    my_Cmm1.fill(0.);
    my_Cmm2.fill(0.);
    my_tod1.fill(0.);
    my_tod2.fill(0.);

    Alm< xcomplex<float> > &blmT = b->blmT();
    Alm< xcomplex<float> > &blmG = b->blmG();
    Alm< xcomplex<float> > &blmC = b->blmC();
    Alm< xcomplex<float> > &slmT = s->slmT();
    Alm< xcomplex<float> > &slmG = s->slmG();
    Alm< xcomplex<float> > &slmC = s->slmC();
    levels::arr<double> cs(nphi), sn(nphi);
    levels::arr<double> cs0(nphi), sn0(nphi);

#pragma omp parallel default(shared)
    {
        double t1 = mpiMgr.Wtime();
        wignergen wgen(lmax, rthetas, conv_acc), wgen_neg(lmax, rthetas, conv_acc);
        if (omp_get_thread_num() == 0) {
            t_wigner_init += mpiMgr.Wtime() - t1;
            ++n_wigner_init;
        }
        // wigner_estimator estimator(lmax, 100);
        for(long mbeam = 0; mbeam < npsi; ++mbeam) {
            const double dsignb = levels::xpow(mbeam, 1);
#pragma omp for schedule(static, 8)
            for (long msky = 0; msky <= lmax; ++msky) {
                const double dsign = levels::xpow(msky, 1);
                const double dsb = dsign * dsignb;
                //estimator.prepare_m(mbeam, msky);
                // if (estimator.canSkip(rthetas[my_ntheta - 1]))
                //     continue; // negligible dmm
                t1 = mpiMgr.Wtime();
                wgen.prepare(mbeam, msky);
                wgen_neg.prepare(mbeam, -msky);
                if (omp_get_thread_num() == 0) {
                    t_wigner_prepare += mpiMgr.Wtime() - t1;
                    ++n_wigner_prepare;
                }
                for (long itheta = 0; itheta < my_ntheta; ++itheta) {
                    int firstl1, firstl2;
                    const levels::arr<double> &dmm = wgen.calc(itheta, firstl1);
                    const levels::arr<double> &dmmneg = wgen_neg.calc(itheta, firstl2);
                    xcomplex<double> &Cmm1_pos = my_Cmm1(msky + lmax, mbeam, itheta);
                    xcomplex<double> &Cmm1_neg = my_Cmm1(-msky + lmax, mbeam, itheta);
                    xcomplex<double> &Cmm2_pos = my_Cmm2(msky + lmax, mbeam, itheta);
                    xcomplex<double> &Cmm2_neg = my_Cmm2(-msky + lmax, mbeam, itheta);
                    const int firstl = std::max(firstl1, firstl2);
                    double dlb = -levels::xpow(firstl, dsignb);
                    for (long ii = firstl; ii <= lmax; ++ii) {
                        dlb = -dlb;
                        // Note that msky in dlm is located to the
                        // left of mb and they have to be interchanged in the convolution
                        const double dMatrixElementmskypos = dsb * dmm[ii];
                        const double dMatrixElementmskyneg = dsb * dmmneg[ii];
                        const double dMatrixElementmskypos2 = dlb * dMatrixElementmskyneg;
                        const double dMatrixElementmskyneg2 = dlb * dMatrixElementmskypos;
                        const xcomplex<float> &sT = slmT(ii, msky);
                        const xcomplex<float> &sG = slmG(ii, msky);
                        const xcomplex<float> &sC = slmC(ii, msky);
                        const xcomplex<float> &bT = blmT(ii, mbeam);
                        const xcomplex<float> &bG = blmG(ii, mbeam);
                        const xcomplex<float> &bC = blmC(ii, mbeam);
                        const double prod1 = sT.re * bT.re + sG.re * bG.re + sC.re * bC.re;
                        const double prod2 = sT.im * bT.im + sG.im * bG.im + sC.im * bC.im;
                        const double prod3 = sT.im * bT.re + sG.im * bG.re + sC.im * bC.re;
                        const double prod4 = sT.re * bT.im + sG.re * bG.im + sC.re * bC.im;
                        const double tmp_1 = prod1 + prod2;
                        const double tmp_2 = prod3 - prod4;
                        const double xtmp_1 = tmp_1 * dMatrixElementmskypos;
                        const double xtmp_2 = tmp_2 * dMatrixElementmskypos;
                        const double xtmp_5 = tmp_1 * dMatrixElementmskypos2;
                        const double xtmp_6 = tmp_2 * dMatrixElementmskypos2;
                        Cmm1_pos.re += xtmp_1;
                        Cmm1_pos.im += xtmp_2;
                        Cmm2_pos.re += xtmp_5;
                        Cmm2_pos.im += xtmp_6;
                        if (msky != 0) {
                            double tmp_3 = dsign * (prod1 - prod2);
                            double tmp_4 = -dsign * (prod3 + prod4);
                            double xtmp_3 = tmp_3 * dMatrixElementmskyneg;
                            double xtmp_4 = tmp_4 * dMatrixElementmskyneg;
                            double xtmp_7 = tmp_3 * dMatrixElementmskyneg2;
                            double xtmp_8 = tmp_4 * dMatrixElementmskyneg2;
                            Cmm1_neg.re += xtmp_3;
                            Cmm1_neg.im += xtmp_4;
                            Cmm2_neg.re += xtmp_7;
                            Cmm2_neg.im += xtmp_8;
                        }
                    }
                }
            }
        }

        const double lmaxinv = double(2 * pi * lmax) / double(2 * lmax + 1);
#pragma omp for schedule(static, 8)
        for (long msky = -lmax; msky <= lmax; ++msky) {
            double arg = -halfpi * msky;
            long ii = lmax + msky;
            cs[ii] = cos(arg);
            sn[ii] = sin(arg);
            arg = ii * lmaxinv;
            cs0[ii] = cos(arg);
            sn0[ii] = sin(arg);
        }
    } // End of parallel section

    todAnnulus(my_tod1, my_Cmm1, cs, sn, cs0, sn0, my_ntheta);
    todAnnulus(my_tod2, my_Cmm2, cs, sn, cs0, sn0, my_ntheta);

    alltoall_datacube(need_itheta1_core, itheta_core, my_itheta, my_tod1, tod1,
                      ithetaoffset1, NThetaIndex1);
    alltoall_datacube(need_itheta2_core, itheta_core, my_itheta, my_tod2, tod2,
                      ithetaoffset2, NThetaIndex2);

    t_conviqt_hemiscm_pol += mpiMgr.Wtime() - tstart;

    if (CMULT_VERBOSITY > 1) {
        std::cerr << corenum << " : Exiting conviqt_hemiscm_pol" << std::endl;
    }
}


void convolver::alltoall_datacube(std::vector<char> &need_itheta_core,
                                  std::vector<int> &itheta_core,
                                  std::vector<long> &my_itheta,
                                  levels::arr3< xcomplex<double> > &my_tod,
                                  levels::arr3< xcomplex<double> > &tod,
                                  const long ithetaoffset, const long NThetaIndex) {
    if (CMULT_VERBOSITY > 1) {
        std::cerr << corenum << " : Entering alltoall_datacube" << std::endl;
    }
    double tstart = mpiMgr.Wtime();
    ++n_alltoall_datacube;

    // Calculate the latitudes we send and receive

    std::vector<int> sendcounts_lat(cores, 0);
    std::vector<int> recvcounts_lat(cores, 0);
    std::vector<int> itheta_send, itheta_recv;

    for (long core = 0; core < cores; ++core) {

        // How many to send to each core

        for (long pos = 0; pos < my_itheta.size(); ++pos) {
            long itheta = my_itheta[pos];
            if (need_itheta_core[core * ntheta + itheta]) {
                ++sendcounts_lat[core];
                // local rather than global latitude index
                itheta_send.push_back(pos);
            }
        }

        // How many to receive from each core

        for (long pos = 0; pos < NThetaIndex; ++pos) {
            long itheta = pos + ithetaoffset;
            if (itheta_core[itheta] == core) {
                ++recvcounts_lat[core];
                // local rather than global latitude index
                itheta_recv.push_back(pos);
            }
        }
    }

    long nsend = itheta_send.size();
    long nrecv = itheta_recv.size();

    // Number of doubles per latitude

    long bufsize = nphi * npsi * 2;
    std::vector<double> sendbuf, recvbuf;
    try {
        sendbuf.resize(nsend * bufsize);
        recvbuf.resize(nrecv * bufsize);
    } catch (std::bad_alloc &e) {
        std::cerr << "alltoall_datacube :  Out of memory allocating "
                  << (nsend + nrecv) * bufsize * 8. / 1024 / 1024
                  << "MB for receiving and sending" << std::endl;
        throw;
    }

    std::vector<int> sendcounts(cores), sdispls(cores);
    std::vector<int> recvcounts(cores), rdispls(cores);

    // Pack send buffer.  Unfortunately the TOD is stored row-major
    // so we are accessing my_tod at a stride of my_ntheta

    double *p = sendbuf.data();
    for (long isend = 0; isend < nsend; ++isend) {
        long pos = itheta_send[isend];
        for (long iphi = 0; iphi < nphi; ++iphi) {
            for (long ipsi = 0; ipsi < npsi; ++ipsi) {
                *(p++) = my_tod(iphi, ipsi, pos).re;
                *(p++) = my_tod(iphi, ipsi, pos).im;
            }
        }
    }

    for (int core = 0; core < cores; ++core) {
        sendcounts[core] = sendcounts_lat[core] * bufsize;
        recvcounts[core] = recvcounts_lat[core] * bufsize;
        if (core == 0) {
            sdispls[core] = 0;
            rdispls[core] = 0;
        } else {
            sdispls[core] = sdispls[core - 1] + sendcounts[core - 1];
            rdispls[core] = rdispls[core - 1] + recvcounts[core - 1];
        }
    }

    int err = MPI_Alltoallv(sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_DOUBLE,
                            recvbuf.data(), recvcounts.data(), rdispls.data(), MPI_DOUBLE,
                            mpiMgr.comm());
    if (err) throw std::runtime_error("Failed to alltoall datacube");

    // Unpack recv buffer

    p = recvbuf.data();
    for (long irecv = 0; irecv < nrecv; ++irecv) {
        long pos = itheta_recv[irecv];
        for (long iphi = 0; iphi < nphi; ++iphi) {
            for (long ipsi = 0; ipsi < npsi; ++ipsi) {
                tod(iphi, ipsi, pos).re = *(p++);
                tod(iphi, ipsi, pos).im = *(p++);
            }
        }
    }

    t_alltoall_datacube += mpiMgr.Wtime() - tstart;

    if (CMULT_VERBOSITY > 1) {
        std::cerr << corenum << " : Exiting alltoall_datacube" << std::endl;
    }
}


void convolver::todAnnulus(levels::arr3<xcomplex<double> > &tod,
                           levels::arr3<xcomplex<double> > &Cmm,
                           levels::arr<double> &cs,
                           levels::arr<double> &sn,
                           levels::arr<double> &cs0,
                           levels::arr<double> &sn0,
                           long NThetaIndex) {
    if (CMULT_VERBOSITY > 1) {
        std::cerr << corenum << " : Entering todAnnulus" << std::endl;
    }
    double tstart = mpiMgr.Wtime();
    ++n_todAnnulus;

#pragma omp parallel default(shared)
    {
        levels::arr< xcomplex<double> > Cmsky(nphi, 0);
        cfft p1(nphi);

#pragma omp for schedule(static, 8)
        for (long msky = -lmax; msky <= lmax; ++msky) {
            const long ii = msky + lmax;
            const double cos1 = cs[ii];
            const double sin1 = sn[ii];
            for (long mbeam = 0; mbeam < npsi; ++mbeam) {
                for (long lat = 0; lat < NThetaIndex; ++lat) {
                    // Get references to the complex parts to
                    // avoid overhead
                    xcomplex<double> &Cmmref = Cmm(ii, mbeam, lat);
                    double &CmmR = Cmmref.re;
                    double &CmmI = Cmmref.im;
                    const double tempR = CmmR;
                    const double tempI = CmmI;
                    CmmR = cos1 * tempR + sin1 * tempI;
                    CmmI = cos1 * tempI - sin1 * tempR;
                }
            }
        } // end of parallel for

#pragma omp for schedule(static, 1)
        for (long mbeam = 0; mbeam < npsi; ++mbeam) {
            for (long lat = 0; lat < NThetaIndex; ++lat) {
                for (long msky = -lmax; msky <= lmax; ++msky) {
                    const long ii = msky + lmax;
                    Cmsky[ii] = Cmm(ii, mbeam, lat);
                }
                p1.backward(Cmsky);
                for (long msky = 0; msky < nphi; msky++) {
                    tod(msky, mbeam, lat) = Cmsky[msky];
                }
            }
        } // end of parallel for

#pragma omp for schedule(static, 8)
        for (long msky = -lmax; msky <= lmax; ++msky) {
            const long ii = msky + lmax;
            const double cos0 = cs0[ii];
            const double sin0 = sn0[ii];
            for (long mbeam = 0; mbeam < npsi; ++mbeam) {
                const long jj = mbeam + lmax;
                const double cos1 = cs[jj];
                const double sin1 = sn[jj];
                double tempR, tempI;
                for (long lat = 0; lat < NThetaIndex; ++lat) {
                    // Get references to the complex parts to
                    // avoid overhead
                    xcomplex<double> &todref = tod(ii, mbeam, lat);
                    double &todR = todref.re;
                    double &todI = todref.im;
                    tempR = todR;  // Copy the value
                    tempI = todI;
                    todR = cos0 * tempR + sin0 * tempI;
                    todI = cos0 * tempI - sin0 * tempR;
                    // Rotate in psi space
                    tempR = todR;  // Copy updated value
                    tempI = todI;
                    todR = cos1 * tempR + sin1 * tempI;
                    todI = cos1 * tempI - sin1 * tempR;
                }
            }
        } // end of parallel for

    } // end of parallel region

    // Finished with convolution and FFT over msky.

    t_todAnnulus += mpiMgr.Wtime() - tstart;

    if (CMULT_VERBOSITY > 1) {
        std::cerr << corenum << " : Exiting todAnnulus" << std::endl;
    }
}


void convolver::interpolTOD(levels::arr<double> &outpntarr1,
                            levels::arr<double> &outpntarr2,
                            levels::arr<double> &TODValue1,
                            levels::arr<double> &TODValue2,
                            const long ntod1, const long ntod2) {
    double tstart = mpiMgr.Wtime();
    ++n_interpolTOD;

    if (CMULT_VERBOSITY > 1) {
        std::cerr << corenum << " : Entering interpolTOD" << std::endl;
    }

    levels::arr<long> itheta0_1, itheta0_2;
    long NThetaIndex1 = 0, NThetaIndex2 = 0;
    levels::arr<long> lowerIndex1, lowerIndex2;
    levels::arr<long> upperIndex1, upperIndex2;
    levels::arr3< xcomplex<double> > TODAsym1, TODAsym2;
    long ithetaoffset1 = 0, ithetaoffset2 = 0;

    if (ntod1 != 0) {
        itheta0SetUp(outpntarr1, ntod1, NThetaIndex1, itheta0_1,
                     lowerIndex1, upperIndex1, TODAsym1);
        ithetaoffset1 = itheta0_1[ntod1 - 1];
    }

    if (ntod2 != 0) {
        itheta0SetUp(outpntarr2, ntod2, NThetaIndex2, itheta0_2,
                     lowerIndex2, upperIndex2, TODAsym2);
        ithetaoffset2 = itheta0_2[ntod2 - 1];
    }

    if (CMULT_VERBOSITY > 1) {
        for (int core = 0; core < cores; ++core) {
            if (core == corenum) {
                std::cerr << corenum << " :"
                          << " NThetaIndex1 = " << NThetaIndex1
                          << " ntod1 = " << ntod1
                          << " NThetaIndex2 = " << NThetaIndex2
                          << " ntod2 = " << ntod2
                          << std::endl;
            }
            mpiMgr.barrier();
        }
    }

    if (pol) {
        conviqt_hemiscm_pol(TODAsym1, TODAsym2,
                            NThetaIndex1, NThetaIndex2,
                            ithetaoffset1, ithetaoffset2);
    } else {
        conviqt_hemiscm(TODAsym1, TODAsym2,
                        NThetaIndex1, NThetaIndex2,
                        ithetaoffset1, ithetaoffset2);
    }

    for (int thetaIndex = 0; thetaIndex < NThetaIndex1; ++thetaIndex) {
        conviqt_tod_loop(lowerIndex1, upperIndex1, outpntarr1, TODAsym1,
                         thetaIndex, itheta0_1, ntod1, TODValue1,
                         thetaIndex);
    }

    for (int thetaIndex = 0; thetaIndex < NThetaIndex2; ++thetaIndex) {
        conviqt_tod_loop(lowerIndex2, upperIndex2, outpntarr2, TODAsym2,
                         thetaIndex, itheta0_2, ntod2, TODValue2,
                         thetaIndex);
    }

    if (CMULT_VERBOSITY > 1) {
        std::cerr << corenum << " : Exiting interpolTOD" << std::endl;
    }

    t_interpolTOD += mpiMgr.Wtime() - tstart;
}


void convolver::conviqt_tod_loop(levels::arr<long> &lowerIndex,
                                 levels::arr<long> &upperIndex,
                                 levels::arr<double> &outpntarr,
                                 levels::arr3<xcomplex<double> > &TODAsym,
                                 const long thetaIndex,
                                 levels::arr<long> &itheta0,
                                 const long ntod,
                                 levels::arr<double> &TODValue,
                                 const long lat) {
    if (CMULT_VERBOSITY > 1) {
        std::cerr << corenum << " : Entering conviqt_tod_loop" << std::endl;
    }

    double tstart = mpiMgr.Wtime();
    ++n_conviqt_tod_loop;
    levels::arr2< xcomplex<double> > conviqtarr;
    try {
        conviqtarr.alloc(nphi, npsi);
    } catch (std::bad_alloc &e) {
        std::cerr << "conviqt_tod_loop : Out of memory allocating "
                  << nphi * npsi * 16. / 1024 / 1024
                  << "MB for conviqtarr" << std::endl;
        throw;
    }

    for (long ii = 0; ii < nphi; ++ii) {
        for (long jj = 0; jj < npsi; ++jj) {
            conviqtarr[ii][jj] = TODAsym(ii, jj, lat);
        }
    }

    weight_ncm_initialize();  // Ensure base_wgt is set-up

    const long thetaoffset = thetaIndex + itheta0[ntod - 1];

#pragma omp parallel default(shared)
    {
        std::vector<double> cosang(npsi), sinang(npsi);
        std::vector<double> wgt1(npoints, 0.);
        std::vector<double> wgt2(npoints, 0.);

#pragma omp for schedule(static, 1)
        for (int ii = lowerIndex[thetaIndex]; ii >= upperIndex[thetaIndex]; --ii) {
            // Co-latitude, theta
            double frac = (outpntarr[5 * ii + 1] - theta0) * inv_delta_theta;
            frac -= itheta0[ii];
            weight_ncm(frac, wgt1);
            // Longitude, phi
            frac = outpntarr[5 * ii] * inv_delta_phi - phioffset;
            frac = levels::fmodulo(frac, double(nphi));
            int iphi0 = int(frac) - ioffset;
            frac -= iphi0;
            if (iphi0 >= nphi) {
                iphi0 -= nphi;
            } else if (iphi0 < 0) {
                iphi0 += nphi;
            }
            weight_ncm(frac, wgt2);
            // Position angle (psi)
            const double omega = outpntarr[5 * ii + 2] + halfpi;
            const double sinomg = sin(omega);
            const double cosomg = cos(omega);
            cosang[0] = 1;
            sinang[0] = 0;
            for (long ipsi = 0; ipsi < npsi - 1; ++ipsi) {
                cosang[ipsi + 1] = cosang[ipsi] * cosomg - sinang[ipsi] * sinomg;
                sinang[ipsi + 1] = sinang[ipsi] * cosomg + cosang[ipsi] * sinomg;
            }
            cosang[0] = 0.5;
            const double weight1 = 2 * wgt1[thetaoffset - itheta0[ii]];
            for (int iphi = 0; iphi < npoints; ++iphi) {
                const double weight = wgt2[iphi] * weight1;
                long iphinew = iphi0 + iphi;
                if (iphinew >= nphi) {
                    iphinew -= nphi;
                }
                const xcomplex<double> *ca = &(conviqtarr[iphinew][0]);
                double x = 0;
                for (long ipsi = 0; ipsi < npsi; ++ipsi) {
                    x += cosang[ipsi] * ca[ipsi].re - sinang[ipsi] * ca[ipsi].im;
                }
                TODValue[ii] += weight * x;
            }
        }

    } // end parallel region

    t_conviqt_tod_loop += mpiMgr.Wtime() - tstart;
    if (CMULT_VERBOSITY > 1) {
        std::cerr << corenum << " : Exiting conviqt_tod_loop" << std::endl;
    }
}


void convolver::arrFillingcm(long ntod,
                             levels::arr<double> &outpntarrx,
                             levels::arr<double> &outpntarr,
                             long offindex) {
    /*
      Copy North or South hemisphere data into a separate pntarr
     */
    if (CMULT_VERBOSITY > 1) {
        std::cerr << corenum << " : Entering arrFillingcm" << std::endl;
    }
    double tstart = mpiMgr.Wtime();
    ++n_arrFillingcm;
    outpntarrx.alloc(5 * ntod);
    for (long ii = 0; ii < ntod; ++ii) {
        double theta = outpntarr[5 * (ii + offindex) + 1];
        if (theta > halfpi) {
            theta = pi - theta;
        }
        outpntarrx[5 * ii] = outpntarr[5 * (ii + offindex)];
        outpntarrx[5 * ii + 1] = theta;
        outpntarrx[5 * ii + 2] = outpntarr[5 * (ii + offindex) + 2];
        outpntarrx[5 * ii + 3] = outpntarr[5 * (ii + offindex) + 3];
        outpntarrx[5 * ii + 4] = outpntarr[5 * (ii + offindex) + 4];
    }
    hpsort_arrTheta(outpntarrx);

    t_arrFillingcm += mpiMgr.Wtime() - tstart;
    if (CMULT_VERBOSITY > 1) {
        std::cerr << corenum << " : Exiting arrFillingcm" << std::endl;
    }
}


int convolver::convolve(pointing &pntarr, bool calibrate) {

    mpiMgr.barrier();

    if (CMULT_VERBOSITY > 1) {
        std::cerr << corenum << " : Convolving : pol = " << pol
                  << " lmax = " << lmax
                  << " beammmax = " << beammmax
                  << " order = " << order
                  << std::endl;
    }

    double tstart = mpiMgr.Wtime(), t1;
    ++n_convolve;

    const long totsize = pntarr.size() / 5;

    // Assign a running index across the communicator to the last column of pntarr
    // It is needed to collect the convolved data.
    long my_offset = 0;
    int err = MPI_Scan(&totsize, &my_offset, 1, MPI_LONG, MPI_SUM, mpiMgr.comm());
    if (err) throw std::runtime_error("Error scanning total TOD size");
    my_offset -= totsize;
    for (long ii = 0; ii < totsize; ++ii) {
        pntarr[5 * ii + 4] = my_offset + ii;
    }

    // Check the angles for conformity

    for (long i = 0; i < totsize; ++i) {
        double theta = pntarr[5 * i + 1];
        if (theta < 0 || theta > pi) {
            throw std::runtime_error("convolver::convolve: Illegal theta: " +
                                     std::to_string(theta) + " not in 0..pi");
        }
    }

    // distribute the co-latitudes.

    levels::arr<double> corethetaarr(cores);
    distribute_colatitudes(pntarr, totsize, corethetaarr);

    long outBetaSegSize;
    levels::arr<int> inBetaSeg, outBetaSeg;
    levels::arr<int> inBetaSegAcc, outBetaSegAcc;
    levels::arr<double> pntarr2, outpntarr;
    levels::arr<int> inOffset, outOffset;
    // Make sure theta=pi doesn't break anything
    const double ratiodeltas = (halfpi + 1e-10) / cores;

    fillingBetaSeg(pntarr, totsize, ratiodeltas, corethetaarr, inBetaSeg);

    mpiMgr.barrier();
    t1 = mpiMgr.Wtime();
    mpiMgr.all2all(inBetaSeg, outBetaSeg);
    t_alltoall += mpiMgr.Wtime() - t1;
    ++n_alltoall;

    for (long ii = 0; ii < totsize; ++ii) {
        // The signal column in pntarr is replaced with a
        // local row number to allow reordering
        pntarr[5 * ii + 3] = ii;
    }

    todRedistribution5cm(pntarr, inBetaSeg, outBetaSeg, inBetaSegAcc, outBetaSegAcc,
                         outBetaSegSize, pntarr2, totsize, inOffset, outOffset,
                         ratiodeltas, corethetaarr);

    pntarr.dealloc();
    std::vector<long> orig_order(totsize);
    if (totsize != 0) {
        t1 = mpiMgr.Wtime();
        for(long ii = 0; ii < totsize; ++ii) {
            orig_order[ii] = long(pntarr2[5 * ii + 3] + .5);
        }
    }

    // Redistribute the data from pntarr2 to outpntarr

    if (totsize != 0 || outBetaSegSize != 0) {
        mpiMgr.barrier();
        t1 = mpiMgr.Wtime();
        mpiMgr.all2allv(pntarr2, inBetaSeg, inOffset,
                        outpntarr, outBetaSeg, outOffset, outBetaSegSize);
        t_alltoall += mpiMgr.Wtime() - t1;
        ++n_alltoall;
    }

    if (totsize > 0) {
        pntarr2.dealloc();
    }

    // Count elements on Northern (ntod1) and Southern (ntod2) hemisphere

    long ntod = outBetaSegSize / 5;
    long ntod1 = 0, ntod2 = 0;
    for (long ii = 0; ii < ntod; ++ii) {
        if (outpntarr[5 * ii + 1] < halfpi) {
            ++ntod1;
        } else {
            ++ntod2;
        }
        // The signal column in outpntarr is replaced with a
        // local row number to allow reordering
        outpntarr[5 * ii + 3] = ii;
    }

    std::vector<long> old_order(ntod);
    if (ntod != 0) {
        t1 = mpiMgr.Wtime();
        hpsort_arrTheta(outpntarr); // Sort according to latitude
        t_sort += mpiMgr.Wtime() - t1;
        ++n_sort;
        for(long ii = 0; ii < ntod; ++ii) {
            old_order[ii] = long(outpntarr[5 * ii + 3] + .5);
        }
    }

    // Run the convolution

    todgen(ntod1, ntod2, outpntarr);

    // Restore the time-based order in outpntarr so we can re-shuffle
    // the TOD back to the original owners.  We need to make a copy of
    // the array for the reorder.

    if (ntod != 0) {
        t1 = mpiMgr.Wtime();
        levels::arr<double> outpntarr_temp(5 * ntod);
        memcpy(&(outpntarr_temp[0]), &(outpntarr[0]), sizeof(double) * ntod * 5);
        for(long new_row = 0; new_row < ntod; ++new_row) {
            long old_row = old_order[new_row];
            memcpy(&(outpntarr[5 * old_row]),
                   &(outpntarr_temp[5 * new_row]),
                   sizeof(double) * 5);
        }
        outpntarr_temp.dealloc();
        old_order.clear();
        t_sort += mpiMgr.Wtime() - t1;
        ++n_sort;
    }

    // Communicate the pointing and TOD to the original owners

    if (totsize != 0 || ntod != 0) {
        mpiMgr.barrier();
        t1 = mpiMgr.Wtime();
        mpiMgr.all2allv(outpntarr, outBetaSeg, outOffset, pntarr,
                        inBetaSeg, inOffset, 5 * totsize);
        t_alltoall += mpiMgr.Wtime() - t1;
        ++n_alltoall;
    }

    if (ntod != 0) {
        outpntarr.dealloc();
    }

    // Now restore the original, time-based order in pntarr

    if (totsize != 0) {
        t1 = mpiMgr.Wtime();
        levels::arr<double> pntarr_temp(5 * totsize);
        memcpy(&(pntarr_temp[0]), &(pntarr[0]), sizeof(double) * totsize * 5);
        for(long new_row = 0; new_row < totsize; ++new_row) {
            long old_row = orig_order[new_row];
            memcpy(&(pntarr[5 * old_row]),
                   &(pntarr_temp[5 * new_row]),
                   sizeof(double) * 5);
        }
        pntarr_temp.dealloc();
        orig_order.clear();
        t_sort += mpiMgr.Wtime() - t1;
        ++n_sort;
    }

    if (totsize != 0 && calibrate) {
        double calibration = 2. / (1. + d->get_epsilon());
        for (long ii = 0; ii < totsize; ++ii) {
            // Insert convolved TOD into the output array
            pntarr[5 * ii + 3] *= calibration;
        }
    }

    t_convolve += mpiMgr.Wtime() - tstart;
    //if (CMULT_VERBOSITY > 1)
    report_timing();

    return 0;
}

void convolver::report_timing() {
    if (corenum == 0)
        std::cerr << "Conviqt::convolve timing:" << std::endl;
    timing_line(std::string("convolve"), t_convolve, n_convolve);
    timing_line(std::string("    distribute_colatitudes"),
                t_distribute_colatitudes, n_distribute_colatitudes);
    timing_line(std::string("    fillingBetaSeg"),
                t_fillingBetaSeg, n_fillingBetaSeg);
    timing_line(std::string("    todRedistribution5cm"),
                t_todRedistribution5cm, n_todRedistribution5cm);
    timing_line(std::string("    todgen"), t_todgen, n_todgen);
    timing_line(std::string("        arrFillingcm"), t_arrFillingcm, n_arrFillingcm);
    timing_line(std::string("        interpolTOD"), t_interpolTOD, n_interpolTOD);
    timing_line(std::string("            itheta0SetUp"), t_itheta0SetUp, n_itheta0SetUp);
    timing_line(std::string("                ithetacalc"), t_ithetacalc, n_ithetacalc);
    timing_line(std::string("            conviqt_hemiscm"),
                t_conviqt_hemiscm, n_conviqt_hemiscm);
    timing_line(std::string("            conviqt_hemiscm_pol"),
                t_conviqt_hemiscm_pol, n_conviqt_hemiscm_pol);
    timing_line(std::string("                wigner_init"), t_wigner_init, n_wigner_init);
    timing_line(std::string("                wigner_prepare"),
                t_wigner_prepare, n_wigner_prepare);
    timing_line(std::string("                todAnnulus"),  t_todAnnulus, n_todAnnulus);
    timing_line(std::string("                    alltoall_datacube"),
                t_alltoall_datacube, n_alltoall_datacube);
    timing_line(std::string("            conviqt_tod_loop"),
                t_conviqt_tod_loop, n_conviqt_tod_loop);
    timing_line(std::string("    alltoall"), t_alltoall, n_alltoall);
    timing_line(std::string("    sort TOD"), t_sort, n_sort);
}

void convolver::timing_line(std::string label, double timer, long counter) {
    double tmin, tmax;
    mpiMgr.reduceRaw<double>(&timer, &tmin, 1, mpiMgr.Min);
    mpiMgr.reduceRaw<double>(&timer, &tmax, 1, mpiMgr.Max);
    double tmean, tstd;
    mpiMgr.allreduceRaw<double>(&timer, &tmean, 1, mpiMgr.Sum);
    tmean /= cores;
    double tsquare = ((timer - tmean) * (timer - tmean));
    mpiMgr.reduceRaw<double>(&tsquare, &tstd, 1, mpiMgr.Sum);
    tstd = sqrt(tstd / cores);

    long cmin, cmax;
    mpiMgr.reduceRaw<long>(&counter, &cmin, 1, mpiMgr.Min);
    mpiMgr.reduceRaw<long>(&counter, &cmax, 1, mpiMgr.Max);
    double cmean, cstd;
    double dcounter = counter;
    mpiMgr.allreduceRaw<double>(&dcounter, &cmean, 1, mpiMgr.Sum);
    cmean /= cores;
    double csquare = ((counter - cmean) * (counter - cmean));
    mpiMgr.reduceRaw<double>(&csquare, &cstd, 1, mpiMgr.Sum);
    cstd = sqrt(cstd / cores);

    if (corenum == 0 && cmax > 0) {
        fprintf(stderr,
                "%-40s %8.2f < %8.2f +- %8.2f < %8.2f "
                "(%8ld < %10.1f +- %10.1f < %8ld)\n",
                label.c_str(),
                tmin, tmean, tstd, tmax,
                cmin, cmean, cstd, cmax);
        //std::cerr << label << " : " << tmin
        //          << " < " << tmean
        //          << " +- " << tstd
        //          << " < " << tmax << std::endl;
    }
}


} // namespace conviqt
