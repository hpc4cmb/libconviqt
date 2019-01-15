#include "conviqt.hpp"

#include <sstream>

#include <cstring>

/*
  int convolve()
      void thetaDeltaThetacm()
      void deltaTheta2()
      void fillingBetaSeg()
          void ratiobetacalcsmallercm()
          void ratiobetacalcgreatercm()
      void todRedistribution5cm()
          void ratiobetacalcsmallercm()
          void ratiobetacalcgreatercm()
      mpiMgr.all2allv()
      hpsort_arrTheta()

      void todgen_v4()
          void arrFillingcm_v2()
          void interpolTOD_arrTestcm_v4()
              void itheta0SetUp()
                  void ithetacalc()
              void conviqt_hemiscm_single()
                  void todAnnulus_v3()
              void conviqt_hemiscm_v4()
                  void todAnnulus_v3()
              void conviqt_tod_loop_v4()
                  void weight_ncm()
          void interpolTOD_arrTestcm_pol_v4()
              void itheta0SetUp()
                  void ithetacalc()
              EITHER
              void conviqt_hemiscm_pol_v4()
                  wignergen()
                  wignergen.prepare()
                  wignergen.calc()
                  void todAnnulus_v3()
              OR
              void conviqt_hemiscm_pol_single()
                  wignergen()
                  wignergen.prepare()
                  wignergen.calc()
                  void todAnnulus_v3()
              void conviqt_tod_loop_pol_v5()
                  void weight_ncm()

      hpsort_DL()
      hpsort_arrTOD()
      mpiMgr.all2allv()
      hpsort_arrTime
      void preReorderingStep()
      hpsort_DDcm()
      mpiMgr.all2allv()
      hpsort_DDcm()

  int set_sky()
  int set_beam()
  int set_detector()

 */

// This file will contain the actual operations on the sky, no I/O

namespace conviqt {

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

    // data cube gridding

    phi0 = halfpi;
    npsi = beammmax + 1;
    nphi = 2 * lmax + 1;
    dphi = 2 * pi / nphi;
    inv_delta_phi = 1. / dphi;
    phioffset = phi0 / dphi;
    halfmargin = 10;
    margin = halfmargin * 2 + 1;
    ntheta = lmax + 1 + margin;
    dtheta = -pi / (ntheta - margin);
    theta0 = pi - halfmargin * dtheta;
    inv_delta_theta = 1. / dtheta;
    max_order = 19;
    npoints = order + 1;
    ioffset = order / 2;

    // profiling

    t_convolve = 0;
    n_convolve = 0;
    t_alltoall = 0;
    n_alltoall = 0;
    t_todRedistribution5cm = 0;
    n_todRedistribution5cm = 0;

    t_todgen_v4 = 0;
    t_arrFillingcm_v2 = 0;
    t_interpolTOD_arrTestcm_v4 = 0;
    t_itheta0SetUp = 0;
    t_ithetacalc = 0;
    t_conviqt_tod_loop_v4 = 0;
    t_conviqt_hemiscm_single = 0;
    t_todAnnulus_v3 = 0;
    t_conviqt_hemiscm_v4 = 0;
    t_interpolTOD_arrTestcm_pol_v4 = 0;
    t_conviqt_hemiscm_pol_v4 = 0;
    t_conviqt_hemiscm_pol_single = 0;
    t_conviqt_tod_loop_pol_v5 = 0;

    n_todgen_v4 = 0;
    n_arrFillingcm_v2 = 0;
    n_interpolTOD_arrTestcm_v4 = 0;
    n_itheta0SetUp = 0;
    n_ithetacalc = 0;
    n_conviqt_tod_loop_v4 = 0;
    n_conviqt_hemiscm_single = 0;
    n_todAnnulus_v3 = 0;
    n_conviqt_hemiscm_v4 = 0;
    n_interpolTOD_arrTestcm_pol_v4 = 0;
    n_conviqt_hemiscm_pol_v4 = 0;
    n_conviqt_hemiscm_pol_single = 0;
    n_conviqt_tod_loop_pol_v5 = 0;

    t_wigner_init = 0;
    t_wigner_prepare = 0;
    t_wigner_calc = 0;

    n_wigner_init = 0;
    n_wigner_prepare = 0;
    n_wigner_calc = 0;

    t_lat_iter = 0;
    t_sincos_iter = 0;

    n_lat_iter = 0;
    n_sincos_iter = 0;

    t_weight_ncm = 0;
    n_weight_ncm = 0;
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


void convolver::weight_ncm(const double x, levels::arr<double> &wgt) {
    double tstart = mpiMgr.Wtime();
    ++n_weight_ncm;
    const int npoints = order + 1;

    if (base_wgt.size() == 0) {
        // Initialize base_wgt only when needed
        base_wgt.resize(npoints, 1);
        for (int m = 0; m < npoints; ++m) {
            for (int n = 0; n<m; ++n)
                base_wgt[m] *= m - n;
            for (int n = m + 1; n < npoints; ++n)
                base_wgt[m] *= m - n;
        }
        for (int m = 0; m < (int)npoints; ++m)
            base_wgt[m] = 1. / base_wgt[m];
    } else if (base_wgt.size() != npoints) {
        throw std::runtime_error("base_wgt changed dimensions");
    }

    std::vector<double> temp_wgt(base_wgt);
    double mul1 = x;
    for (int m = 1; m < npoints; ++m) {
        temp_wgt[m] *= mul1;
        mul1 *= x - m;
    }

    double mul2 = x - npoints + 1;
    for (int m = 1; m < npoints; ++m) {
        temp_wgt[npoints - m - 1] *= mul2;
        mul2 *= x - npoints + m + 1;
    }

    memcpy(&(wgt[0]), &(temp_wgt[0]), sizeof(double) * npoints);
    t_weight_ncm += mpiMgr.Wtime() - tstart;
}


void convolver::weight_ncm(const double x, std::vector<double> &wgt) {
    // This method is called from multiple threads at once
    double tstart = mpiMgr.Wtime();
    ++n_weight_ncm;
    const unsigned int npoints = order + 1;

    if (base_wgt.size() == 0) {
        // Initialize base_wgt only when needed
        base_wgt.resize(npoints, 1);
        for (int m = 0; m < (int)npoints; ++m) {
            for (int n = 0; n < m; ++n)
                base_wgt[m] *= m - n;
            for (int n = m + 1; n < (int)npoints; ++n)
                base_wgt[m] *= m - n;
        }
        for (int m = 0; m < (int)npoints; ++m)
            base_wgt[m] = 1. / base_wgt[m];
    } else if (base_wgt.size() != npoints) {
        throw std::runtime_error("base_wgt changed dimensions");
    }

    std::vector<double> temp_wgt(base_wgt);
    double mul1 = x;
    for (int m = 1; m < (int)npoints; ++m) {
        temp_wgt[m] *= mul1;
        mul1 *= x - m;
    }

    double mul2 = x - npoints + 1;
    for (int m = 1; m < (int)npoints; ++m) {
        temp_wgt[npoints - m - 1] *= mul2;
        mul2 *= x - npoints + m + 1;
    }

    memcpy(&(wgt[0]), &(temp_wgt[0]), sizeof(double) * npoints);
    t_weight_ncm += mpiMgr.Wtime() - tstart;
}


void convolver::thetaDeltaThetacm(const int corenum,
                                  const double thetaini,
                                  double &theta,
                                  double &deltatheta) {
    if (corenum == 0) {
        theta = 0;
        deltatheta = thetaini;
    } else if (corenum == cores - 1) {
        deltatheta = 1. / cores;
        theta = thetaini - deltatheta;
    } else {
        double a = -0.5 * cos(thetaini);
        double b = -sin(thetaini) + thetaini * cos(thetaini);
        double c = (thetaini * sin(thetaini) -
                    0.5 * cos(thetaini) * thetaini * thetaini -
                    1. / cores);
        theta = -(b + sqrt(b * b - 4 * a * c)) / (2 * a);
        deltatheta = thetaini - theta;
    }

    std::cerr << "corenum = " << corenum
              << " thetaini = " << thetaini * 180 / pi
              << " theta = " << theta * 180 / pi
              << " deltatheta = " << deltatheta * 180 / pi
              << std::endl;
}


void convolver::deltaTheta2(const int corenum,
                            const double thetaini,
                            levels::arr<double> &dbeta) {
    dbeta.alloc(corenum + 1);
    double dy = 0.5 * (1 - thetaini);
    double yshift = -0.7;
    double theta2 = pi - abs(asin(1 - dy)) + yshift;
    double theta1 = pi - abs(asin(sin(theta2) - thetaini));
    double dt = (theta2 - theta1) / (corenum + 1.);
    double dbeta0 = sin(theta2);

    for (long ii = corenum; ii >= 0; --ii) {
        dbeta[ii] = dbeta0 - sin(theta2 - (corenum + 1 - ii) * dt);
        dbeta0 = sin(theta2 - (corenum + 1 - ii) * dt);
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

    long totsize_tot = 0;
    mpiMgr.allreduceRaw<long>(&totsize, &totsize_tot, 1, mpiMgr.Sum);

    // Build a latitude hit map using the same gridding as ithetacalc
    // size_t nbin = 10000;
    // double wbin = halfpi / nbin;
    arr<int> my_hits(ntheta, 0.), hits(ntheta, 0.);
    for (long ii = 0; ii < totsize; ++ii) {
        double beta = pntarr[5 * ii + 1];
        long itheta = beta_to_itheta(beta);
        ++my_hits[itheta];
    }
    
    mpiMgr.allreduce(my_hits, hits, mpiMgr.Sum);

    long nhit_bin = 0;
    for (long itheta = 0; itheta < ntheta; ++itheta) {
        if (hits[itheta] > 0) {
            ++nhit_bin;
        }
    }

    // the bin width, dtheta, is negative for convenience elsewhere
    double thetaleft = -nhit_bin * dtheta;
    long nhitleft = totsize_tot, nbinleft = nhit_bin;

    corethetaarr[0] = 0;
    long itheta = 0;    
    for (int corenum = cores - 1; corenum > 0; --corenum) {
        // Width of the bin is based on remaining latitude
        // and samples, which ever we meet first
        double divisor = corenum + 1;
        double nbin_target = nbinleft / divisor;
        double dtheta_target = thetaleft / divisor;
        double nhit_target = nhitleft / divisor;
        long bin_hit = 0;
        double bin_width = 0;
        long bin_nbin = 0;
        while (itheta < ntheta) {
            // Skip over isolatitudes that have no hits
            if (hits[itheta] != 0) {
                bin_hit += hits[itheta];
                nhitleft -= hits[itheta];
                bin_width -= dtheta;  // dtheta is negative
                thetaleft += dtheta;                
                ++bin_nbin;
                --nbinleft;
                double score1 = bin_hit / nhit_target;
                double score2 = bin_width / dtheta_target;
                double score3 = bin_nbin / nbin_target;
                if (score1 + score2 + score3 > 4 ||
                    (score1 > 1 && score2 > 1 && score3 > 1)) {
                    break;
                }
            }
            ++itheta;
        }
        double theta = itheta_to_beta1(itheta);        
        corethetaarr[corenum] = theta;
        if (mpiMgr.master()) std::cerr
                                 << " corenum = " << corenum
                                 << " theta = " << theta * 180 / pi
                                 << " nbinleft = " << nbinleft << " / " << nhit_bin
                                 << " thetaleft = " << thetaleft * 180 / pi << " / " << -nhit_bin * dtheta * 180 / pi
                                 << " nhitleft = " << nhitleft << " / " << totsize_tot
                                 << std::endl;
    }
    corethetaarr[0] = 0;

    /*
    double dtheta, newtheta, thetaini = halfpi;
    for (int corenum = cores - 1; corenum >= 0; --corenum) {
        if (thetaini > 0.18) {
            thetaDeltaThetacm(corenum, thetaini, newtheta, dtheta);
        } else {
            levels::arr<double> dbeta;
            deltaTheta2(corenum, thetaini, dbeta);
            for (long jj = corenum; jj >= 0; --jj) {
                newtheta -= dbeta[jj];
                corethetaarr[jj] = newtheta;
                thetaini = newtheta;
            }
            corethetaarr[0] = 0;
            break;
        }
        corethetaarr[corenum] = newtheta;
        thetaini = newtheta;
    }
    */
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

    if (CMULT_VERBOSITY > 1)
        std::cout << "Entered fillingBetaSeg in core = " << corenum << std::endl;

    inBetaSeg.alloc(cores);
    inBetaSeg.fill(0);

    for (long ii = 0; ii < arrsize; ++ii) {
        double theta = pntarr[5 * ii + 1];
        if (theta >= halfpi) {
            theta = pi - theta;
        }
        // Initialize corenum to a rough estimate of which task
        // should own this colatitude ..
        int corenum = theta / ratiodeltas;
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

    if (CMULT_VERBOSITY > 1)
        std::cout << "Leaving fillingBetaSeg in core = " << corenum << std::endl;
}


void convolver::ratiobetacalcgreatercm(int &corenum, const double theta,
                                       levels::arr<double> &corethetaarr) {
    while (corenum < cores - 1) {
        if (theta < corethetaarr[corenum + 1])
            break;
        ++corenum;
    }
}


void convolver::ratiobetacalcsmallercm(int &corenum, const double theta,
                                       levels::arr<double> &corethetaarr) {
    while (corenum > 0) {
        if (theta > corethetaarr[corenum])
            break;
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
    if (CMULT_VERBOSITY > 1)
        std::cout << "Entered todRedistribution5cm" << std::endl;

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
        } catch (std::bad_alloc & e) {
            std::cerr << " todRedistribution5cm : Out of memory allocating "
                  << (5 * totsize) * 8. / 1024 / 1024 << "MB for pntarr2"
                      << std::endl;
            throw;
        }
        // Pack the data in pntarr into pntarr2 for sending
        pntarr2.fill(0);
        levels::arr<int> countbeta(cores, 0);
        for (long ii = 0; ii < totsize; ++ii) {
            double theta = pntarr[5 * ii + 1];
            if (theta < 0 || theta > pi) {
                throw std::runtime_error("ERROR: Illegal latitude in pntarr: theta = " +
                                         std::to_string(theta) + " not in 0..pi");
            }
            if (theta >= halfpi) {
                theta = pi - theta;
            }
            int corenum = theta / ratiodeltas;
            if (theta >= corethetaarr[corenum]){
                ratiobetacalcgreatercm(corenum, theta, corethetaarr);
            } else {
                ratiobetacalcsmallercm(corenum, theta, corethetaarr);
            }
            const long i1 = 5 * ii;
            const long i2 = inOffset[corenum] + countbeta[corenum];
            for (long k = 0; k < 5; ++k) {
                pntarr2[i2 + k] = pntarr[i1 + k];
            }
            countbeta[corenum] += 5;
	}
    }

    if (CMULT_VERBOSITY > 1) {
        std::cout << "Leaving todRedistribution5cm" << std::endl;
    }

    t_todRedistribution5cm = mpiMgr.Wtime() - tstart;
}


void convolver::preReorderingStep(const long ntod1, const long ntod2,
                                  levels::arr<double> &todAll,
                                  levels::arr<double> &todTest_arr,
                                  levels::arr<double> &todTest_arr2) {
    if (ntod1 + ntod2 != 0) {
        todAll.alloc(ntod1 + ntod2);
    }

    if (ntod1 != 0) {
        for (long ii = 0; ii < ntod1; ++ii) {
            todAll[ii] = todTest_arr[ii];
        }
        todTest_arr.dealloc();
    }

    if (ntod2 != 0) {
        for (long ii = 0; ii < ntod2; ++ii) {
            todAll[ntod1 + ii] = todTest_arr2[ii];
        }
        todTest_arr2.dealloc();
    }
}


void convolver::todgen_v4(const long ntod1, const long ntod2,
                          levels::arr<double> &todTest_arr,
                          levels::arr<double> &timeTest_arr,
                          levels::arr<double> &todTest_arr2,
                          levels::arr<double> &timeTest_arr2,
                          levels::arr<double> &outpntarr) {
    double tstart = mpiMgr.Wtime();
    ++n_todgen_v4;
    levels::arr<double> outpntarr1, outpntarr2;

    if (ntod1 != 0) {
        // Collect northern hemisphere data from outpntarr to outpntarr1
        arrFillingcm_v2(ntod1, timeTest_arr, outpntarr1, outpntarr, 0);
        todTest_arr.alloc(ntod1);
        todTest_arr.fill(0);
    }
    if (ntod2 != 0) {
        // Collect southern hemisphere data from outpntarr to outpntarr2
        arrFillingcm_v2(ntod2, timeTest_arr2, outpntarr2, outpntarr, ntod1);
        todTest_arr2.alloc(ntod2);
        todTest_arr2.fill(0);
    }

    if (CMULT_VERBOSITY > 1)
        std::cout << "corenum = " << corenum << "  order = " << order
                  << "  ntod1 = " << ntod1 << "  ntod2 = " << ntod2 << std::endl;

    mpiMgr.barrier();
    if (!pol) {
        interpolTOD_arrTestcm_v4(outpntarr1, outpntarr2,
                                 todTest_arr, todTest_arr2, ntod1, ntod2);
    } else {
        interpolTOD_arrTestcm_pol_v4(outpntarr1, outpntarr2,
                                     todTest_arr, todTest_arr2, ntod1, ntod2);
    }

    t_todgen_v4 += mpiMgr.Wtime() - tstart;
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
    } catch (std::bad_alloc & e) {
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


void convolver::conviqt_hemiscm_v4(levels::arr3<xcomplex<double> > &tod1,
                                   levels::arr3<xcomplex<double> > &tod2,
                                   long NThetaIndex1,
                                   levels::arr<double> &rthetas) {
    double tstart = mpiMgr.Wtime();
    ++n_conviqt_hemiscm_v4;
    levels::arr3<xcomplex<double> > Cmm(nphi, npsi, NThetaIndex1);
    Cmm.fill(0.);
    levels::arr3<xcomplex<double> > Cmm2(nphi, npsi, NThetaIndex1);
    Cmm2.fill(0.);

    Alm< xcomplex<float> > & blmT = b->blmT();
    Alm< xcomplex<float> > & slmT = s->slmT();
    levels::arr<double> cs(nphi), sn(nphi);
    levels::arr<double> cs0(nphi), sn0(nphi);

#pragma omp parallel default(shared)
    {
        wignergen wgen(lmax, rthetas, conv_acc), wgen_neg(lmax, rthetas, conv_acc);
        //wigner_estimator estimator(lmax,100);
        for(long beamIndex = 0; beamIndex <= beammmax; ++beamIndex) {
            double dsignb = levels::xpow(beamIndex, 1);
#pragma omp for schedule(static, 8)
            for (long msky = 0; msky <= lmax; ++msky) {
                const double dsign = levels::xpow(msky, 1);
                const double dsb = dsign * dsignb;
                // estimator.prepare_m(beamIndex, msky);
                // if (estimator.canSkip(rthetas[NThetaIndex1-1]))
                //     continue; // negligible dmm
                wgen.prepare(beamIndex, msky);
                wgen_neg.prepare(beamIndex, -msky);
                for (long lat = 0; lat < NThetaIndex1; ++lat) {
                    xcomplex<double> Cmm_pos = Cmm(msky + lmax, beamIndex, lat);
                    xcomplex<double> Cmm_neg = Cmm(-msky + lmax, beamIndex, lat);
                    xcomplex<double> Cmm2_pos = Cmm2(msky + lmax, beamIndex, lat);
                    xcomplex<double> Cmm2_neg = Cmm2(-msky + lmax, beamIndex, lat);
                    int firstl1, firstl2;
                    const levels::arr<double> &dmm = wgen.calc(lat, firstl1);
                    const levels::arr<double> &dmmneg = wgen_neg.calc(lat, firstl2);
                    const int firstl = (firstl1 > firstl2) ? firstl1 : firstl2;
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
                        const xcomplex<float> bT = blmT(ii, beamIndex);
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
                        Cmm_pos.re += xtmp_1;
                        Cmm_pos.im += xtmp_2;
                        Cmm2_pos.re += xtmp_5;
                        Cmm2_pos.im += xtmp_6;
                        if (msky != 0) {
                            const double tmp_3 = dsign * (prod1 - prod2);
                            const double tmp_4 = -dsign * (prod3 + prod4);
                            const double xtmp_3 = tmp_3 * dMatrixElementmskyneg;
                            const double xtmp_4 = tmp_4 * dMatrixElementmskyneg;
                            const double xtmp_7 = tmp_3 * dMatrixElementmskyneg2;
                            const double xtmp_8 = tmp_4 * dMatrixElementmskyneg2;
                            Cmm_neg.re += xtmp_3;
                            Cmm_neg.im += xtmp_4;
                            Cmm2_neg.re += xtmp_7;
                            Cmm2_neg.im += xtmp_8;
                        }
                    }
                }
            }
        }
        const double lmaxinv = 2 * pi * lmax / (2. * lmax + 1.);
#pragma omp for schedule(static, 8)
        for (long msky = -lmax; msky <= lmax; ++msky) {
            double arg = -halfpi * msky;
            long ii = lmax + msky;
            cs[ii] = cos(arg);
            sn[ii] = sin(arg);
            // arg = 2.*pi*(msky+lmax)*lmax/(2.*lmax+1.);
            arg = ii * lmaxinv;
            cs0[ii] = cos(arg);
            sn0[ii] = sin(arg);
        }
    }

    todAnnulus_v3(tod1, Cmm, cs, sn, cs0, sn0, NThetaIndex1);
    todAnnulus_v3(tod2, Cmm2, cs, sn, cs0, sn0, NThetaIndex1);

    t_conviqt_hemiscm_v4 += mpiMgr.Wtime() - tstart;
}


void convolver::conviqt_hemiscm_pol_v4(levels::arr3< xcomplex<double> > &tod1,
                                       levels::arr3< xcomplex<double> > &tod2,
                                       const long NThetaIndex,
                                       levels::arr<double> &rthetas) {
    double tstart = mpiMgr.Wtime();
    ++n_conviqt_hemiscm_pol_v4;
    levels::arr3< xcomplex<double> > Cmm(nphi, npsi, NThetaIndex);
    Cmm.fill(0.);
    levels::arr3< xcomplex<double> > Cmm2(nphi, npsi, NThetaIndex);
    Cmm2.fill(0.);

    Alm< xcomplex<float> > & blmT = b->blmT();
    Alm< xcomplex<float> > & blmG = b->blmG();
    Alm< xcomplex<float> > & blmC = b->blmC();
    Alm< xcomplex<float> > & slmT = s->slmT();
    Alm< xcomplex<float> > & slmG = s->slmG();
    Alm< xcomplex<float> > & slmC = s->slmC();
    levels::arr<double> cs(nphi), sn(nphi);
    levels::arr<double> cs0(nphi), sn0(nphi);

#pragma omp parallel default(shared)
    {
        double t1 = mpiMgr.Wtime();
        wignergen wgen(lmax, rthetas, conv_acc), wgen_neg(lmax, rthetas, conv_acc);
        t_wigner_init += mpiMgr.Wtime() - t1;
        ++n_wigner_init;
        // wigner_estimator estimator(lmax, 100);
        for(long beamIndex = 0; beamIndex <= beammmax; ++beamIndex) {
            const double dsignb = levels::xpow(beamIndex, 1);
#pragma omp for schedule(static, 8)
            for (long msky = 0; msky <= lmax; ++msky) {
                const double dsign = levels::xpow(msky, 1);
                const double dsb = dsign * dsignb;
                //estimator.prepare_m(beamIndex, msky);
                // if (estimator.canSkip(rthetas[NThetaIndex - 1]))
                //     continue; // negligible dmm
                t1 = mpiMgr.Wtime();
                wgen.prepare(beamIndex, msky);
                wgen_neg.prepare(beamIndex, -msky);
                t_wigner_prepare += mpiMgr.Wtime() - t1;
                ++n_wigner_prepare;
                for (long lat = 0; lat < NThetaIndex; ++lat) {
                    int firstl1, firstl2;
                    t1 = mpiMgr.Wtime();
                    const levels::arr<double> &dmm = wgen.calc(lat, firstl1);
                    const levels::arr<double> &dmmneg = wgen_neg.calc(lat, firstl2);
                    t_wigner_calc += mpiMgr.Wtime() - t1;
                    ++n_wigner_calc;
                    t1 = mpiMgr.Wtime();
                    xcomplex<double> Cmm_pos = Cmm(msky + lmax, beamIndex, lat);
                    xcomplex<double> Cmm_neg = Cmm(-msky + lmax, beamIndex, lat);
                    xcomplex<double> Cmm2_pos = Cmm2(msky + lmax, beamIndex, lat);
                    xcomplex<double> Cmm2_neg = Cmm2(-msky + lmax, beamIndex, lat);
                    const int firstl = (firstl1 > firstl2) ? firstl1 : firstl2;
                    double dlb = -levels::xpow(firstl, dsignb);
                    for (long ii = firstl; ii <= lmax; ++ii) {
                        dlb = -dlb;
                        // Note that msky in dlm is located to the
                        // left of mb and they have to be interchanged in the convolution
                        const double dMatrixElementmskypos = dsb * dmm[ii];
                        const double dMatrixElementmskyneg = dsb * dmmneg[ii];
                        const double dMatrixElementmskypos2 = dlb * dMatrixElementmskyneg;
                        const double dMatrixElementmskyneg2 = dlb * dMatrixElementmskypos;
                        const xcomplex<float> sT = slmT(ii, msky);
                        const xcomplex<float> sG = slmG(ii, msky);
                        const xcomplex<float> sC = slmC(ii, msky);
                        const xcomplex<float> bT = blmT(ii, beamIndex);
                        const xcomplex<float> bG = blmG(ii, beamIndex);
                        const xcomplex<float> bC = blmC(ii, beamIndex);
                        const double prod1 = sT.re * bT.re + sG.re * bG.re + sC.re * bC.re;
                        const double prod3 = sT.im * bT.re + sG.im * bG.re + sC.im * bC.re;
                        const double prod2 = sT.im * bT.im + sG.im * bG.im + sC.im * bC.im;
                        const double prod4 = sT.re * bT.im + sG.re * bG.im + sC.re * bC.im;
                        const double tmp_1 = prod1 + prod2;
                        const double tmp_2 = prod3 - prod4;
                        const double xtmp_1 = tmp_1 * dMatrixElementmskypos;
                        const double xtmp_2 = tmp_2 * dMatrixElementmskypos;
                        const double xtmp_5 = tmp_1 * dMatrixElementmskypos2;
                        const double xtmp_6 = tmp_2 * dMatrixElementmskypos2;
                        Cmm_pos.re += xtmp_1;
                        Cmm_pos.im += xtmp_2;
                        Cmm2_pos.re += xtmp_5;
                        Cmm2_pos.im += xtmp_6;
                        if (msky != 0) {
                            double tmp_3 = dsign * (prod1 - prod2);
                            double tmp_4 = -dsign * (prod3 + prod4);
                            double xtmp_3 = tmp_3 * dMatrixElementmskyneg;
                            double xtmp_4 = tmp_4 * dMatrixElementmskyneg;
                            double xtmp_7 = tmp_3 * dMatrixElementmskyneg2;
                            double xtmp_8 = tmp_4 * dMatrixElementmskyneg2;
                            Cmm_neg.re += xtmp_3;
                            Cmm_neg.im += xtmp_4;
                            Cmm2_neg.re += xtmp_7;
                            Cmm2_neg.im += xtmp_8;
                        }
                    }
                    t_lat_iter += mpiMgr.Wtime() - t1;
                    ++n_lat_iter;
                }
            }
        }
        t1 = mpiMgr.Wtime();
        const double lmaxinv = 2 * pi * lmax / (2. * lmax + 1.);
#pragma omp for schedule(static, 8)
        for (long msky = -lmax; msky <= lmax; ++msky) {
            double arg = -halfpi * msky;
            long ii = lmax + msky;
            cs[ii] = cos(arg);
            sn[ii] = sin(arg);
            // arg = 2.*pi*(msky+lmax)*lmax/(2.*lmax+1.);
            arg = ii * lmaxinv;
            cs0[ii] = cos(arg);
            sn0[ii] = sin(arg);
        }
        t_sincos_iter += mpiMgr.Wtime() - t1;
        ++n_sincos_iter;
    } // End of parallel section

  todAnnulus_v3(tod1, Cmm, cs, sn, cs0, sn0, NThetaIndex);
  todAnnulus_v3(tod2, Cmm2, cs, sn, cs0, sn0, NThetaIndex);

  t_conviqt_hemiscm_pol_v4 += mpiMgr.Wtime() - tstart;
}


void convolver::conviqt_hemiscm_single(levels::arr3<xcomplex<double> > &tod,
                                       long NThetaIndex,
                                       levels::arr<double> &rthetas) {
    double tstart = mpiMgr.Wtime();
    ++n_conviqt_hemiscm_single;
    levels::arr3<xcomplex<double> > Cmm(nphi, npsi, NThetaIndex);
    Cmm.fill(0.);

    Alm< xcomplex<float> > & blmT = b->blmT();
    Alm< xcomplex<float> > & slmT = s->slmT();
    levels::arr<double> cs(nphi), sn(nphi);
    levels::arr<double> cs0(nphi), sn0(nphi);

#pragma omp parallel default(shared)
    {
        wignergen wgen(lmax, rthetas, conv_acc), wgen_neg(lmax, rthetas, conv_acc);
        // wigner_estimator estimator(lmax, 100);

        for(long beamIndex = 0; beamIndex <= beammmax; ++beamIndex) {
            const double dsignb = levels::xpow(beamIndex, 1);
#pragma omp for schedule(static,8)
            for (long msky = 0; msky <= lmax; ++msky) {
                const double dsign = levels::xpow(msky, 1);
                const double dsb = dsign * dsignb;
                // estimator.prepare_m(beamIndex, msky);
                // if (estimator.canSkip(rthetas[NThetaIndex - 1]))
                //     continue; // negligible dmm
                wgen.prepare(beamIndex, msky);
                wgen_neg.prepare(beamIndex, -msky);
                for (long lat = 0; lat < NThetaIndex; ++lat) {
                    xcomplex<double> Cmm_pos = Cmm(msky + lmax, beamIndex, lat);
                    xcomplex<double> Cmm_neg = Cmm(-msky + lmax, beamIndex, lat);
                    int firstl1, firstl2;
                    const levels::arr<double> &dmm = wgen.calc(lat, firstl1);
                    const levels::arr<double> &dmmneg = wgen_neg.calc(lat, firstl2);
                    const int firstl = (firstl1 > firstl2) ? firstl1 : firstl2;
                    for (long ii = firstl; ii <= lmax; ++ii) {
                        // Note that msky in dlm is located to the left of
                        // mb and they have to be interchanged in the convolution
                        const double dMatrixElementmskypos = dsb * dmm[ii];
                        const double dMatrixElementmskyneg = dsb * dmmneg[ii];
                        const xcomplex<float> sT = slmT(ii, msky);
                        const xcomplex<float> bT = blmT(ii, beamIndex);
                        const double prod1 = sT.re * bT.re;
                        const double prod3 = sT.im * bT.re;
                        const double prod2 = sT.im * bT.im;
                        const double prod4 = sT.re * bT.im;
                        const double tmp_1 = prod1 + prod2;
                        const double tmp_2 = prod3 - prod4;
                        const double xtmp_1 = tmp_1 * dMatrixElementmskypos;
                        const double xtmp_2 = tmp_2 * dMatrixElementmskypos;
                        Cmm_pos.re += xtmp_1;
                        Cmm_pos.im += xtmp_2;
                        if (msky != 0) {
                            const double tmp_3 = dsign * (prod1 - prod2);
                            const double tmp_4 = -dsign * (prod3 + prod4);
                            const double xtmp_3 = tmp_3 * dMatrixElementmskyneg;
                            const double xtmp_4 = tmp_4 * dMatrixElementmskyneg;
                            Cmm_neg.re += xtmp_3;
                            Cmm_neg.im += xtmp_4;
                        }
                    }
                }
            }
        }
        const double lmaxinv = 2 * pi * lmax / (2. * lmax + 1.);
#pragma omp for schedule(static, 8)
        for (long msky = -lmax; msky <= lmax; ++msky) {
            double arg = -halfpi * msky;
            long ii = lmax + msky;
            cs[ii] = cos(arg);
            sn[ii] = sin(arg);
            // arg = 2 * pi * (msky + lmax) * lmax / (2. * lmax + 1.);
            arg = ii * lmaxinv;
            cs0[ii] = cos(arg);
            sn0[ii] = sin(arg);
        }
    } // end parallel region

    todAnnulus_v3(tod, Cmm, cs, sn, cs0, sn0, NThetaIndex);

    t_conviqt_hemiscm_single += mpiMgr.Wtime() - tstart;
}


void convolver::conviqt_hemiscm_pol_single(levels::arr3<xcomplex<double> > &tod,
                                           long NThetaIndex,
                                           levels::arr<double> &rthetas) {
    double tstart = mpiMgr.Wtime();
    ++n_conviqt_hemiscm_pol_single;
    levels::arr3<xcomplex<double> > Cmm(nphi, npsi, NThetaIndex);
    Cmm.fill(0.);

    Alm< xcomplex<float> > & blmT = b->blmT();
    Alm< xcomplex<float> > & blmG = b->blmG();
    Alm< xcomplex<float> > & blmC = b->blmC();
    Alm< xcomplex<float> > & slmT = s->slmT();
    Alm< xcomplex<float> > & slmG = s->slmG();
    Alm< xcomplex<float> > & slmC = s->slmC();
    levels::arr<double> cs(nphi), sn(nphi);
    levels::arr<double> cs0(nphi), sn0(nphi);

#pragma omp parallel default(shared)
    {
        double t1 = mpiMgr.Wtime();
        wignergen wgen(lmax, rthetas, conv_acc), wgen_neg(lmax, rthetas, conv_acc);
        t_wigner_init += mpiMgr.Wtime() - t1;
        ++n_wigner_init;
        // wigner_estimator estimator(lmax, 100);
        for (long beamIndex = 0; beamIndex <= beammmax; ++beamIndex) {
            const double dsignb = levels::xpow(beamIndex, 1);
#pragma omp for schedule(static, 8)
            for (long msky = 0; msky <= lmax; ++msky) {
                const double dsign = levels::xpow(msky, 1);
                const double dsb = dsign * dsignb;
                //estimator.prepare_m(beamIndex, msky);
                //if (estimator.canSkip(rthetas[NThetaIndex - 1]))
                //    continue; // negligible dmm
                t1 = mpiMgr.Wtime();
                wgen.prepare(beamIndex, msky);
                wgen_neg.prepare(beamIndex, -msky);
                t_wigner_prepare += mpiMgr.Wtime() - t1;
                ++n_wigner_prepare;
                for (long lat = 0; lat < NThetaIndex; ++lat) {
                    int firstl1, firstl2;
                    t1 = mpiMgr.Wtime();
                    const levels::arr<double> &dmm = wgen.calc(lat, firstl1);
                    const levels::arr<double> &dmmneg = wgen_neg.calc(lat, firstl2);
                    t_wigner_calc += mpiMgr.Wtime() - t1;
                    ++n_wigner_calc;
                    t1 = mpiMgr.Wtime();
                    xcomplex<double> Cmm_pos = Cmm(msky + lmax, beamIndex, lat);
                    xcomplex<double> Cmm_neg = Cmm(-msky + lmax, beamIndex, lat);
                    const int firstl = (firstl1 > firstl2) ? firstl1 : firstl2;
                    for (long ii = firstl; ii <= lmax; ++ii) {
                        // Note that msky in dlm is located to the left
                        // of mb and they have to be interchanged in the convolution
                        const double dMatrixElementmskypos = dsb * dmm[ii];
                        const xcomplex<float> sT = slmT(ii, msky);
                        const xcomplex<float> sG = slmG(ii, msky);
                        const xcomplex<float> sC = slmC(ii, msky);
                        const xcomplex<float> bT = blmT(ii, beamIndex);
                        const xcomplex<float> bG = blmG(ii, beamIndex);
                        const xcomplex<float> bC = blmC(ii, beamIndex);
                        const double prod1 = sT.re * bT.re + sG.re * bG.re + sC.re * bC.re;
                        const double prod3 = sT.im * bT.re + sG.im * bG.re + sC.im * bC.re;
                        const double prod2 = sT.im * bT.im + sG.im * bG.im + sC.im * bC.im;
                        const double prod4 = sT.re * bT.im + sG.re * bG.im + sC.re * bC.im;
                        const double tmp_1 = prod1 + prod2;
                        const double tmp_2 = prod3 - prod4;
                        const double xtmp_1 = tmp_1 * dMatrixElementmskypos;
                        const double xtmp_2 = tmp_2 * dMatrixElementmskypos;
                        Cmm_pos.re += xtmp_1;
                        Cmm_pos.im += xtmp_2;
                        if (msky != 0) {
                            const double tmp_3 = prod1 - prod2;
                            const double tmp_4 = -prod3 - prod4;
                            const double dMatrixElementmskyneg = dsb * dmmneg[ii];
                            const double xtmp_3 = dsign * tmp_3 * dMatrixElementmskyneg;
                            const double xtmp_4 = dsign * tmp_4 * dMatrixElementmskyneg;
                            Cmm_neg.re += xtmp_3;
                            Cmm_neg.im += xtmp_4;
                        }
                    }
                    t_lat_iter += mpiMgr.Wtime() - t1;
                    ++n_lat_iter;
                }
            }
        } // end of parallel for
        t1 = mpiMgr.Wtime();
        const double lmaxinv = 2 * pi * lmax / (2. * lmax + 1.);
#pragma omp for schedule(static, 8)
        for (long msky = -lmax; msky <= lmax; ++msky) {
            double arg = -halfpi * msky;
            long ii = lmax + msky;
            cs[ii] = cos(arg);
            sn[ii] = sin(arg);
            // arg = 2 * pi * (msky + lmax) * lmax / (2. * lmax + 1.);
            arg = ii * lmaxinv;
            cs0[ii] = cos(arg);
            sn0[ii] = sin(arg);
        }
        t_sincos_iter += mpiMgr.Wtime() - t1;
        ++n_sincos_iter;
    } // end of parallel region

    todAnnulus_v3(tod, Cmm, cs, sn, cs0, sn0, NThetaIndex);

    t_conviqt_hemiscm_pol_single += mpiMgr.Wtime() - tstart;
}


void convolver::todAnnulus_v3(levels::arr3<xcomplex<double> > &tod,
                              levels::arr3<xcomplex<double> > &Cmm,
                              levels::arr<double> &cs,
                              levels::arr<double> &sn,
                              levels::arr<double> &cs0,
                              levels::arr<double> &sn0,
                              long NThetaIndex) {
    double tstart = mpiMgr.Wtime();
    ++n_todAnnulus_v3;

#pragma omp parallel default(shared)
    {
        levels::arr< xcomplex<double> > Cmsky(nphi, 0);
        cfft p1(nphi);
#pragma omp for schedule(static, 8)
        for (long msky = -lmax; msky <= lmax; ++msky) {
            const long ii = msky + lmax;
            for (long beamIndex = 0; beamIndex <= beammmax; ++beamIndex) {
                for (long lat = 0; lat < NThetaIndex; ++lat) {
                    //double dPhi = -halfpi;
                    const double tmpR = Cmm(ii, beamIndex, lat).re;
                    const double tmpI = Cmm(ii, beamIndex, lat).im;
                    Cmm(ii, beamIndex, lat).re = cs[ii] * tmpR + sn[ii] * tmpI;
                    Cmm(ii, beamIndex, lat).im = cs[ii] * tmpI - sn[ii] * tmpR;
                }
            }
        } // end of parallel for
#pragma omp for schedule(static, 1)
        for (long beamIndex = 0; beamIndex <= beammmax; ++beamIndex) {
            for (long lat = 0; lat < NThetaIndex; ++lat) {
                for (long msky = -lmax; msky <= lmax; ++msky) {
                    const long ii = msky + lmax;
                    Cmsky[ii] = Cmm(ii, beamIndex, lat);
                }
                p1.backward(Cmsky);
                for (long msky = 0; msky < nphi; msky++)
                    tod(msky, beamIndex, lat) = Cmsky[msky];
            }
        } // end of parallel for
#pragma omp for schedule(static, 8)
        for (long msky = -lmax; msky <= lmax; ++msky) {
            const long ii = msky + lmax;
            for (long beamIndex = 0; beamIndex <= beammmax; ++beamIndex) {
                const long jj = beamIndex + lmax;
                for (long lat = 0; lat < NThetaIndex; ++lat) {
                    double tmpR = tod(ii, beamIndex, lat).re;
                    double tmpI = tod(ii, beamIndex, lat).im;
                    tod(ii, beamIndex, lat).re = cs0[ii] * tmpR + sn0[ii] * tmpI;
                    tod(ii, beamIndex, lat).im = cs0[ii] * tmpI - sn0[ii] * tmpR;
                    // Rotate in psi space
                    tmpR = tod(ii, beamIndex, lat).re;
                    tmpI = tod(ii, beamIndex, lat).im;
                    tod(ii, beamIndex, lat).re = cs[jj] * tmpR + sn[jj] * tmpI;
                    tod(ii, beamIndex, lat).im = cs[jj] * tmpI - sn[jj] * tmpR;
                }
            }
        } // end of parallel for
    } // end of parallel region
    // Note that now, -lmax <= msky <= lmax corresponds to the angle
    // 0 <= phi <= 2.*pi/2./(2.*lmax+1.) and -pi <= phi <= -2.pi/(2.*lmax+1.)
    // Finished with convolution and FFT over msky.

    t_todAnnulus_v3 += mpiMgr.Wtime() - tstart;
}


void convolver::interpolTOD_arrTestcm_v4(levels::arr<double> &outpntarr1,
                                         levels::arr<double> &outpntarr2,
                                         levels::arr<double> &TODValue1,
                                         levels::arr<double> &TODValue2,
                                         const long ntod1, const long ntod2) {
    double tstart = mpiMgr.Wtime();
    ++n_interpolTOD_arrTestcm_v4;
    if (CMULT_VERBOSITY > 1) {
        std::cout << "Entered interpolTOD_arrTestcm_v4" << std::endl;
    }
    levels::arr<long> itheta0_1, itheta0_2;
    long NThetaIndex1(0), NThetaIndex2(0);
    levels::arr<long> lowerIndex1, lowerIndex2;
    levels::arr<long> upperIndex1, upperIndex2;
    levels::arr3<xcomplex<double> > TODAsym1, TODAsym2;

    if (ntod1 != 0) {
        itheta0SetUp(outpntarr1, ntod1, NThetaIndex1, itheta0_1,
                     lowerIndex1, upperIndex1, TODAsym1);
    }

    if (ntod2 != 0) {
        itheta0SetUp(outpntarr2, ntod2, NThetaIndex2, itheta0_2,
                     lowerIndex2, upperIndex2, TODAsym2);
    }

    if (CMULT_VERBOSITY > 1) {
        std::cout << "DONE with NThetaIndex1 = " << NThetaIndex1
                  << "  and NThetaIndex2 = " << NThetaIndex2 << "  corenum = "
                  << corenum << "  ntod1 = " << ntod1 << "  ntod2 = " << ntod2 << std::endl;
    }

    bool fNThetaIndex = true;
    if (ntod1 != 0 && ntod2 != 0) {
        if (itheta0_1[ntod1 - 1] != itheta0_2[ntod2 - 1]) {
            fNThetaIndex = false;
        }
        if (NThetaIndex1 != NThetaIndex2) {
            fNThetaIndex = false;
        }
    }
    if (ntod1 == 0 || ntod2 == 0) {
        fNThetaIndex = false;
    }

    levels::arr<double> rthetas1, rthetas2;
    if (NThetaIndex1 != 0) {
        rthetas1.alloc(NThetaIndex1);
        for (int thetaIndex = 0; thetaIndex < NThetaIndex1; ++thetaIndex) {
            rthetas1[thetaIndex] = itheta_to_beta1(itheta0_1[ntod1 - 1] + thetaIndex);
        }
    }
    if (NThetaIndex2 != 0) {
        rthetas2.alloc(NThetaIndex2);
        for (int thetaIndex = 0; thetaIndex < NThetaIndex2; ++thetaIndex) {
            rthetas2[thetaIndex] = itheta_to_beta2(itheta0_2[ntod2 - 1] + thetaIndex);
        }
    }

    if (fNThetaIndex) {
        conviqt_hemiscm_v4(TODAsym1, TODAsym2, NThetaIndex1, rthetas1);
    } else {
        if (NThetaIndex1 != 0) {
            conviqt_hemiscm_single(TODAsym1, NThetaIndex1, rthetas1);
        }
        if (NThetaIndex2 != 0) {
            conviqt_hemiscm_single(TODAsym2, NThetaIndex2, rthetas2);
        }
    }

    if (CMULT_VERBOSITY > 1) {
        std::cout << "After conviqt if" << std::endl;
    }

    if (ntod1 > 0) {
        for (int thetaIndex = 0; thetaIndex < NThetaIndex1; ++thetaIndex) {
            conviqt_tod_loop_v4(lowerIndex1, upperIndex1, outpntarr1, TODAsym1,
                                thetaIndex, itheta0_1, ntod1, TODValue1, thetaIndex);
        }
    }

    if (ntod2 > 0) {
        for (int thetaIndex = 0; thetaIndex < NThetaIndex2; ++thetaIndex) {
            conviqt_tod_loop_v4(lowerIndex2, upperIndex2, outpntarr2, TODAsym2,
                                thetaIndex, itheta0_2, ntod2, TODValue2, thetaIndex);
        }
    }

    if (CMULT_VERBOSITY > 1) {
        std::cout << "Leaving interpolTOD_arrTestcm_v4" << std::endl;
    }

    t_interpolTOD_arrTestcm_v4 += mpiMgr.Wtime() - tstart;
}


void convolver::interpolTOD_arrTestcm_pol_v4(levels::arr<double> &outpntarr1,
                                             levels::arr<double> &outpntarr2,
                                             levels::arr<double> &TODValue1,
                                             levels::arr<double> &TODValue2,
                                             long ntod1, long ntod2) {
    double tstart = mpiMgr.Wtime();
    ++n_interpolTOD_arrTestcm_pol_v4;

    if (CMULT_VERBOSITY > 1) {
        std::cout << "Entered interpolTOD_arrTestcm_pol_v4" << std::endl;
    }

    levels::arr<long> itheta0_1, itheta0_2;
    long NThetaIndex1(0), NThetaIndex2(0);
    levels::arr<long> lowerIndex1, lowerIndex2;
    levels::arr<long> upperIndex1, upperIndex2;
    levels::arr3< xcomplex<double> > TODAsym1, TODAsym2;

    if (ntod1 != 0) {
        itheta0SetUp(outpntarr1, ntod1, NThetaIndex1, itheta0_1,
                     lowerIndex1, upperIndex1, TODAsym1);
    }
    if (ntod2 != 0) {
        itheta0SetUp(outpntarr2, ntod2, NThetaIndex2, itheta0_2,
                     lowerIndex2, upperIndex2, TODAsym2);
    }

    // DEBUG begin
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
    // DEBUG end

    if (CMULT_VERBOSITY > 1) {
        std::cout << "DONE with NThetaIndex1 = " << NThetaIndex1
                  << "  and NThetaIndex2 = " << NThetaIndex2
                  << "  corenum = " << corenum << "  ntod1 = " << ntod1
                  << "  ntod2 = " << ntod2 << std::endl;
    }

    bool fNThetaIndex = true;
    if (ntod1 != 0 && ntod2 != 0) {
        if (itheta0_1[ntod1 - 1] != itheta0_2[ntod2 - 1]) {
            fNThetaIndex = false;
        }
        if (NThetaIndex1 != NThetaIndex2) {
            fNThetaIndex = false;
        }
    }

    if (ntod1 == 0 || ntod2 == 0) {
        fNThetaIndex = false;
    }

    levels::arr<double> rthetas1, rthetas2;
    if (NThetaIndex1 != 0) {
        rthetas1.alloc(NThetaIndex1);
        for (int thetaIndex = 0; thetaIndex < NThetaIndex1; ++thetaIndex) {
            rthetas1[thetaIndex] = itheta_to_beta1(itheta0_1[ntod1 - 1] + thetaIndex);
        }
    }

    if (NThetaIndex2 != 0) {
        rthetas2.alloc(NThetaIndex2);
        for (int thetaIndex = 0; thetaIndex < NThetaIndex2; ++thetaIndex) {
            rthetas2[thetaIndex] = itheta_to_beta2(itheta0_2[ntod2 - 1] + thetaIndex);
        }
    }

    if (fNThetaIndex) {
        conviqt_hemiscm_pol_v4(TODAsym1, TODAsym2, NThetaIndex1, rthetas1);
    } else {
        if (NThetaIndex1 != 0) {
            conviqt_hemiscm_pol_single(TODAsym1, NThetaIndex1, rthetas1);
        }
        if (NThetaIndex2 != 0) {
            conviqt_hemiscm_pol_single(TODAsym2, NThetaIndex2, rthetas2);
        }
    }

    if (CMULT_VERBOSITY > 1) {
        std::cout << "After conviqt if" << std::endl;
    }

    if (ntod1 > 0) {
        for (int thetaIndex = 0; thetaIndex < NThetaIndex1; ++thetaIndex) {
            conviqt_tod_loop_pol_v5(lowerIndex1, upperIndex1, outpntarr1, TODAsym1,
                                    thetaIndex, itheta0_1, ntod1, TODValue1,
                                    thetaIndex);
        }
    }

    if (ntod2 > 0) {
        for (int thetaIndex = 0; thetaIndex < NThetaIndex2; ++thetaIndex) {
            conviqt_tod_loop_pol_v5(lowerIndex2, upperIndex2, outpntarr2, TODAsym2,
                                    thetaIndex, itheta0_2, ntod2, TODValue2,
                                    thetaIndex);
        }
    }

    if (CMULT_VERBOSITY > 1) {
        std::cout << "Leaving interpolTOD_arrTestcm_pol_v4" << std::endl;
    }

    t_interpolTOD_arrTestcm_pol_v4 += mpiMgr.Wtime() - tstart;
}


void convolver::conviqt_tod_loop_v4(levels::arr<long> &lowerIndex,
                                    levels::arr<long> &upperIndex,
                                    levels::arr<double> &outpntarr,
                                    levels::arr3<xcomplex<double> > &TODAsym,
                                    const long thetaIndex,
                                    levels::arr<long> &itheta0,
                                    const long ntod,
                                    levels::arr<double> &TODValue, const long lat) {
    double tstart = mpiMgr.Wtime();
    ++n_conviqt_tod_loop_v4;
    levels::arr2< xcomplex<double> > conviqtarr;
    try {
        conviqtarr.alloc(nphi, beammmax + 1);
    } catch (std::bad_alloc & e) {
        std::cerr << "conviqt_tod_loop_v4 : Out of memory allocating "
                  << nphi * (beammmax + 1) * 16. / 1024 / 1024
                  << "MB for conviqtarr" << std::endl;
        throw;
    }

    for (long ii = 0; ii < nphi; ++ii) {
        for (long jj = 0; jj < beammmax + 1; ++jj) {
            conviqtarr[ii][jj] = TODAsym(ii, jj, lat);
        }
    }

    // Call weight_ncm once unthreaded to make sure the auxiliary arrays are allocated

    std::vector<double> wgt_temp(max_order + 1, 0.);
    weight_ncm(0., wgt_temp);

#pragma omp parallel default(shared)
    {
        std::vector<double> cosang(beammmax + 1), sinang(beammmax + 1);
        std::vector<double> wgt1(max_order + 1, 0.);
        std::vector<double> wgt2(max_order + 1, 0.);
#pragma omp for schedule(static, 1)
        for (int ii = lowerIndex[thetaIndex]; ii >= upperIndex[thetaIndex]; --ii) {
            // Note that the larger the ii the smaller frac is
            // and the smaller itheta0[ii] is.            
            double frac = (outpntarr[5 * ii + 1] - theta0) * inv_delta_theta;
            frac -= itheta0[ii];
            weight_ncm(frac, wgt1);
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
            const double omega = outpntarr[5 * ii + 2] + halfpi;
            const double sinomg = sin(omega);
            const double cosomg = cos(omega);
            cosang[0] = 1;
            sinang[0] = 0;
            for (long ipsi = 0; ipsi < beammmax; ++ipsi) {
                cosang[ipsi + 1] = cosang[ipsi] * cosomg - sinang[ipsi] * sinomg;
                sinang[ipsi + 1] = sinang[ipsi] * cosomg + cosang[ipsi] * sinomg;
            }
            cosang[0] = 0.5;
            const double weight1 = 2 * wgt1[thetaIndex - (itheta0[ii] - itheta0[ntod - 1])];
            for (int phiIndex = 0; phiIndex < npoints; ++phiIndex) {
                const double weight = wgt2[phiIndex] * weight1;
                long newphiIndex = iphi0 + phiIndex;
                if (newphiIndex >= nphi)
                    newphiIndex -= nphi;
                xcomplex<double> *ca = &(conviqtarr[newphiIndex][0]);
                double *cang = &(cosang[0]);
                double *sang = &(sinang[0]);
                for (long ipsi = 0; ipsi <= beammmax; ++ipsi) {
                    TODValue[ii] += weight * ((*cang) * (*ca).re - (*sang) * (*ca).im);
                    ++ca;
                    ++cang;
                    ++sang;
                }
            }
        }
    } // end parallel region

    t_conviqt_tod_loop_v4 += mpiMgr.Wtime() - tstart;
}


void convolver::conviqt_tod_loop_pol_v5(levels::arr<long> &lowerIndex,
                                        levels::arr<long> &upperIndex,
                                        levels::arr<double> &outpntarr,
                                        levels::arr3<xcomplex<double> > &TODAsym,
                                        long thetaIndex,
                                        levels::arr<long> &itheta0,
                                        long ntod,
                                        levels::arr<double> &TODValue,
                                        long lat) {
    double tstart = mpiMgr.Wtime();
    ++n_conviqt_tod_loop_pol_v5;
    levels::arr2< xcomplex<double> > conviqtarr;
    try {
        conviqtarr.alloc(nphi, beammmax + 1);
    } catch (std::bad_alloc & e) {
        std::cerr << "conviqt_tod_loop_pol_v5 : Out of memory allocating "
                  << nphi * (beammmax + 1) * 16. / 1024 / 1024 << "MB for conviqtarr"
                  << std::endl;
        throw;
    }

    for (long ii = 0; ii < nphi; ++ii) {
        for (long jj = 0; jj < beammmax + 1; ++jj) {
            conviqtarr[ii][jj] = TODAsym(ii, jj, lat);
        }
    }

    // Call weight_ncm once unthreaded to make sure the auxiliary arrays are allocated

    std::vector<double> wgt_temp(max_order + 1, 0.);
    weight_ncm(0., wgt_temp);

#pragma omp parallel default(shared)
    {
        std::vector<double> cosang(beammmax + 1), sinang(beammmax + 1);
        std::vector<double> wgt1(max_order + 1, 0.);
        std::vector<double> wgt2(max_order + 1, 0.);
#pragma omp for schedule(static, 1)
        for (int ii = lowerIndex[thetaIndex]; ii >= upperIndex[thetaIndex]; --ii) {
            // Note that the larger the ii the smaller frac is
            // and the smaller itheta0[ii] is.
            double frac = (outpntarr[5 * ii + 1] - theta0) * inv_delta_theta;
            frac -= itheta0[ii];
            weight_ncm(frac, wgt1);
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
            const double omega = outpntarr[5 * ii + 2] + halfpi;
            const double sinomg = sin(omega);
            const double cosomg = cos(omega);
            cosang[0] = 1;
            sinang[0] = 0;
            for (long ipsi = 0; ipsi < beammmax; ++ipsi) {
                cosang[ipsi + 1] = cosang[ipsi] * cosomg - sinang[ipsi] * sinomg;
                sinang[ipsi + 1] = sinang[ipsi] * cosomg + cosang[ipsi] * sinomg;
            }
            cosang[0] = 0.5;
            double weight1 = 2 * wgt1[thetaIndex - (itheta0[ii] - itheta0[ntod - 1])];
            for (int iphi = 0; iphi < npoints; ++iphi) {
                double weight = wgt2[iphi] * weight1;
                long iphinew = iphi0 + iphi;
                if (iphinew >= nphi) {
                    iphinew -= nphi;
                }
                xcomplex<double> *ca = &(conviqtarr[iphinew][0]);
                double *cang = &(cosang[0]);
                double *sang = &(sinang[0]);
                for (long ipsi = 0; ipsi <= beammmax; ++ipsi) {
                    TODValue[ii] += weight * ((*cang) * (*ca).re - (*sang) * (*ca).im);
                    ++ca;
                    ++cang;
                    ++sang;
                }
            }
        }
    }  // end parallel region

    t_conviqt_tod_loop_pol_v5 += mpiMgr.Wtime() - tstart;
}


void convolver::arrFillingcm_v2(long ntod,
                                levels::arr<double> &timeTest_arr,
                                levels::arr<double> &outpntarrx,
                                levels::arr<double> &outpntarr,
                                long offindex) {
    double tstart = mpiMgr.Wtime();
    ++n_arrFillingcm_v2;
    timeTest_arr.alloc(ntod);
    outpntarrx.alloc(5 * ntod);
    for (long ii = 0; ii < ntod; ++ii) {
        double theta = (outpntarr[5 * (ii + offindex) + 1] > halfpi) ?
            pi - outpntarr[5 * (ii + offindex) + 1] : outpntarr[5 * (ii + offindex) + 1];
        outpntarrx[5 * ii] = outpntarr[5 * (ii + offindex)];
        outpntarrx[5 * ii + 1] = theta;
        outpntarrx[5 * ii + 2] = outpntarr[5 * (ii + offindex) + 2];
        outpntarrx[5 * ii + 3] = outpntarr[5 * (ii + offindex) + 3];
        outpntarrx[5 * ii + 4] = outpntarr[5 * (ii + offindex) + 4];
    }
    hpsort_arrTheta(outpntarrx);
    for (long ii = 0; ii < ntod; ++ii)
        timeTest_arr[ii] = outpntarrx[5 * ii + 4];

    t_arrFillingcm_v2 += mpiMgr.Wtime() - tstart;
}


int convolver::convolve(pointing & pntarr, bool calibrate) {

    double tstart = mpiMgr.Wtime();
    ++n_convolve;
    const int lmax_sky = s->get_lmax();
    const int lmax_beam = b->get_lmax();
    if (lmax < 0) {
        lmax = lmax_sky < lmax_beam ? lmax_sky : lmax_beam;
    } else if (lmax > lmax_sky || lmax > lmax_beam) {
        throw std::runtime_error("Convolver lmax exceeds input expansion order.");
    }

    const int mmax_beam = b->get_mmax();
    if (beammmax < 0) {
        beammmax = mmax_beam;
    } else if (beammmax > mmax_beam) {
        throw std::runtime_error("Convolver mmax exceeds input expansion order.");
    }

    const long totsize = pntarr.size() / 5;

    // Assign a running index across the communicator to the last column of pntarr
    // It is needed to collect the convolved data.
    long my_offset = 0;
    int err = MPI_Scan(&totsize, &my_offset, 1, MPI_LONG, MPI_SUM, mpiMgr.comm());
    if (err != 0) {
        throw std::runtime_error("Error scannign total TOD size");
    }
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

    // DEBUG begin
    if (mpiMgr.master()) {
        for (int corenum=0; corenum < cores; ++corenum) {
            double theta = corethetaarr[corenum];
            std::cerr << corenum << " : theta = " << theta * 180 / pi
                      << std::endl;
        }
    }
    // DEBUG end

    long outBetaSegSize;
    levels::arr<int> inBetaSeg, outBetaSeg;
    levels::arr<int> inBetaSegAcc, outBetaSegAcc;
    levels::arr<double> pntarr2, outpntarr;
    levels::arr<int> inOffset, outOffset;
    // Make sure theta=pi doesn't break anything
    const double ratiodeltas = (halfpi + 1e-10) / cores;

    fillingBetaSeg(pntarr, totsize, ratiodeltas, corethetaarr, inBetaSeg);
    mpiMgr.barrier();
    double t1 = mpiMgr.Wtime();
    mpiMgr.all2all(inBetaSeg, outBetaSeg);
    t_alltoall += mpiMgr.Wtime() - t1;
    ++n_alltoall;

    todRedistribution5cm(pntarr, inBetaSeg, outBetaSeg, inBetaSegAcc, outBetaSegAcc,
                         outBetaSegSize, pntarr2, totsize, inOffset, outOffset,
                         ratiodeltas, corethetaarr);

    // Store input TOD for comparison

    levels::arr<double> todtmp;
    if (totsize > 0) {
        todtmp.alloc(totsize);
        for (long ii = 0; ii < totsize; ++ii)
            todtmp[ii] = pntarr[5 * ii + 3];
        pntarr.dealloc();
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

    if (totsize > 0)
        pntarr2.dealloc();

    // Count elements on Northern (ntod1) and Southern (ntod2) hemisphere

    long ntod1 = 0, ntod2 = 0;
    for (long ii = 0; ii < outBetaSegSize / 5; ++ii) {
        if (outpntarr[5 * ii + 1] < halfpi) {
            ntod1++;
        } else {
            ntod2++;
        }
        // The signal column in outpntarr is replaced with a
        // local row number to allow reordering
        outpntarr[5 * ii + 3] = ii;
    }

    if (outBetaSegSize != 0)
        hpsort_arrTheta(outpntarr); // Sort according to latitude

    // Run the convolution

    levels::arr<double> todTest_arr, todTest_arr2, timeTest_arr, timeTest_arr2;

    todgen_v4(ntod1, ntod2, todTest_arr, timeTest_arr,
              todTest_arr2, timeTest_arr2, outpntarr);

    // Sort outpntarr according to the 4th column
    if (outBetaSegSize != 0) {
        hpsort_arrTOD(outpntarr);
    }

    if (totsize != 0 || outBetaSegSize != 0) {
        mpiMgr.barrier();
        t1 = mpiMgr.Wtime();
        mpiMgr.all2allv(outpntarr, outBetaSeg, outOffset, pntarr,
                        inBetaSeg, inOffset, 5 * totsize);
        t_alltoall += mpiMgr.Wtime() - t1;
        ++n_alltoall;
    }

    if (outBetaSegSize != 0) {
        outpntarr.dealloc();
    }

    // Sort according to the 5th column, time stamp
    if (totsize != 0) {
        hpsort_arrTime(pntarr);
    }

    levels::arr<double> todAll;
    preReorderingStep(ntod1, ntod2, todAll, todTest_arr, todTest_arr2);
    levels::arr<double> timeAll;
    preReorderingStep(ntod1, ntod2, timeAll, timeTest_arr, timeTest_arr2);

    // sort time and signal according to time
    if (outBetaSegSize != 0) {
        hpsort_DDcm(timeAll, todAll);
    }

    for (long ii = 0; ii < cores; ++ii) {
        inBetaSeg[ii] = inBetaSeg[ii] / 5;
        outBetaSeg[ii] = outBetaSeg[ii] / 5;
        inOffset[ii] = inOffset[ii] / 5;
        outOffset[ii] = outOffset[ii] / 5;
    }

    // collect the convolved TOD

    levels::arr<double> outtodarr, outnumarr;
    if (totsize != 0 || outBetaSegSize != 0) {
        mpiMgr.barrier();
        t1 = mpiMgr.Wtime();
        mpiMgr.all2allv(todAll, outBetaSeg, outOffset, outtodarr,
                        inBetaSeg, inOffset, totsize);
        mpiMgr.all2allv(timeAll, outBetaSeg, outOffset, outnumarr,
                        inBetaSeg, inOffset, totsize);
        t_alltoall += mpiMgr.Wtime() - t1;
        ++n_alltoall;
        ++n_alltoall;
    }

    if (outBetaSegSize != 0) {
        todAll.dealloc();
        timeAll.dealloc();
    }

    if (CMULT_VERBOSITY > 1) {
        double maxtodall(0), mintodall(1e10);
        for (long ii = 0; ii < totsize; ++ii) {
            maxtodall = (maxtodall < outtodarr[ii]) ? outtodarr[ii] : maxtodall;
            mintodall = (mintodall > outtodarr[ii]) ? outtodarr[ii] : mintodall;
        }

        if (totsize > 0)
            std::cout << "maxtodall = " << maxtodall  << "   mintodall = "
                      << mintodall << "   outtodarr.size()/totsize = "
                      << outtodarr.size() / totsize * 1. << std::endl;
    }

    double calibration = 1;
    if (calibrate) {
        calibration = 2. / (1. + d->get_epsilon());
    }

    if (totsize > 0) {
        hpsort_DDcm(outnumarr, outtodarr); // Sort time and signal according to time
        outnumarr.dealloc();
        double maxtoddiff = 0;
        for (long ii = 0; ii < totsize; ++ii) {
            // Insert convolved TOD into the output array
            pntarr[5 * ii + 3] = calibration * outtodarr[ii];
            if (CMULT_VERBOSITY > 1) {
                maxtoddiff = (abs(todtmp[ii]-pntarr[5 * ii + 3]) > maxtoddiff) ?
                    abs(todtmp[ii] - pntarr[5 * ii + 3]) : maxtoddiff;
                if (ii % 100000 == 0)
                    std::cout << "todtmp[ii] = " << todtmp[ii] << "  pntarr[5*ii+3] = "
                              << pntarr[5 * ii + 3] << "  difference = "
                              << abs(todtmp[ii] - pntarr[5 * ii + 3]) << std::endl;
            }
        }
        if (CMULT_VERBOSITY > 1)
            std::cout << "  corenum = " << corenum << "   maxtoddiff = "
                      << maxtoddiff << "   calibration = " << calibration << std::endl;
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
    timing_line(std::string("    todRedistribution5cm"),
                t_todRedistribution5cm, n_todRedistribution5cm);
    timing_line(std::string("    todgen_v4"), t_todgen_v4, n_todgen_v4);
    timing_line(std::string("        arrFillingcm_v2"), t_arrFillingcm_v2, n_arrFillingcm_v2);
    timing_line(std::string("        interpolTOD_arrTestcm_v4"),
                t_interpolTOD_arrTestcm_v4, n_interpolTOD_arrTestcm_v4);
    timing_line(std::string("        interpolTOD_arrTestcm_pol_v4"),
                t_interpolTOD_arrTestcm_pol_v4, n_interpolTOD_arrTestcm_pol_v4);
    timing_line(std::string("            itheta0SetUp"), t_itheta0SetUp, n_itheta0SetUp);
    timing_line(std::string("                ithetacalc"), t_ithetacalc, n_ithetacalc);
    timing_line(std::string("            conviqt_hemiscm_single"),
                t_conviqt_hemiscm_single, n_conviqt_hemiscm_single);
    timing_line(std::string("            conviqt_hemiscm_v4"),
                t_conviqt_hemiscm_v4, n_conviqt_hemiscm_v4);
    timing_line(std::string("            conviqt_hemiscm_pol_v4"),
                t_conviqt_hemiscm_pol_v4, n_conviqt_hemiscm_pol_v4);
    timing_line(std::string("            conviqt_hemiscm_pol_single"),
                t_conviqt_hemiscm_pol_single, n_conviqt_hemiscm_pol_single);
    timing_line(std::string("                wigner_init"), t_wigner_init, n_wigner_init);
    timing_line(std::string("                wigner_prepare"),
                t_wigner_prepare, n_wigner_prepare);
    timing_line(std::string("                wigner_calc"), t_wigner_calc, n_wigner_calc);
    timing_line(std::string("                lat_iter"), t_lat_iter, n_lat_iter);
    timing_line(std::string("                sincos_iter"), t_sincos_iter, n_sincos_iter);
    timing_line(std::string("                todAnnulus_v3"),
                t_todAnnulus_v3, n_todAnnulus_v3);
    timing_line(std::string("            conviqt_tod_loop_v4"),
                t_conviqt_tod_loop_v4, n_conviqt_tod_loop_v4);
    timing_line(std::string("            conviqt_tod_loop_pol_v5"),
                t_conviqt_tod_loop_pol_v5, n_conviqt_tod_loop_pol_v5);
    timing_line(std::string("                weight_ncm"), t_weight_ncm, n_weight_ncm);
    timing_line(std::string("    alltoall"), t_alltoall, n_alltoall);
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
