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

      void todgen()
          void arrFillingcm()
          void interpolTOD_arrTestcm()
              void itheta0SetUp()
                  void ithetacalc()
              void conviqt_hemiscm_alltoall()
                  void todAnnulus()
              void conviqt_tod_loop()
                  void weight_ncm()
          void interpolTOD_arrTestcm_pol()
              void itheta0SetUp()
                  void ithetacalc()
              void conviqt_hemiscm_pol_alltoall()
                  wignergen()
                  wignergen.prepare()
                  wignergen.calc()
                  void todAnnulus()
              void conviqt_tod_loop_pol()
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

    phi0 = halfpi;
    npsi = this->beammmax + 1;
    nphi = 2 * this->lmax + 1;
    dphi = 2 * pi / nphi;
    inv_delta_phi = 1. / dphi;
    phioffset = phi0 / dphi;
    halfmargin = 10;
    margin = halfmargin * 2 + 1;
    ntheta = this->lmax + 1 + margin;
    dtheta = -pi / (ntheta - margin);
    theta0 = pi - halfmargin * dtheta;
    inv_delta_theta = 1 / dtheta;
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
    t_distribute_colatitudes = 0;
    n_distribute_colatitudes = 0;

    t_todgen = 0;
    t_arrFillingcm = 0;
    t_interpolTOD_arrTestcm = 0;
    t_itheta0SetUp = 0;
    t_ithetacalc = 0;
    t_conviqt_tod_loop = 0;
    t_conviqt_hemiscm_alltoall = 0;
    t_todAnnulus = 0;
    t_interpolTOD_arrTestcm_pol = 0;
    t_conviqt_tod_loop_pol = 0;
    t_conviqt_hemiscm_pol_alltoall = 0;
    t_alltoall_datacube = 0;

    n_todgen = 0;
    n_arrFillingcm = 0;
    n_interpolTOD_arrTestcm = 0;
    n_itheta0SetUp = 0;
    n_ithetacalc = 0;
    n_conviqt_tod_loop = 0;
    n_conviqt_hemiscm_alltoall = 0;
    n_todAnnulus = 0;
    n_interpolTOD_arrTestcm_pol = 0;
    n_conviqt_tod_loop_pol = 0;
    n_conviqt_hemiscm_pol_alltoall = 0;
    n_alltoall_datacube = 0;

    t_wigner_init = 0;
    t_wigner_prepare = 0;
    t_wigner_calc = 0;

    n_wigner_init = 0;
    n_wigner_prepare = 0;
    n_wigner_calc = 0;

    t_lat_iter = 0;

    n_lat_iter = 0;

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
        if (mpiMgr.master()) std::cerr
                                 << " corenum = " << corenum
                                 << " ibin = " << ibin  << " / " << nbin
                                 << " theta = " << theta * 180 / pi
                                 << " bin_hit = " << bin_hit
                                 << " nhitleft = " << nhitleft << " / " << totsize_tot
                                 << std::endl; // DEBUG
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

    if (CMULT_VERBOSITY > 1) {
        std::cerr << "Entered fillingBetaSeg in core = " << corenum << std::endl;
    }

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

    if (CMULT_VERBOSITY > 1) {
        std::cerr << "Leaving fillingBetaSeg in core = " << corenum << std::endl;
    }
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
        std::cerr << "Leaving todRedistribution5cm" << std::endl;
    }

    t_todRedistribution5cm += mpiMgr.Wtime() - tstart;
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


void convolver::todgen(const long ntod1, const long ntod2,
                       levels::arr<double> &todTest_arr,
                       levels::arr<double> &timeTest_arr,
                       levels::arr<double> &todTest_arr2,
                       levels::arr<double> &timeTest_arr2,
                       levels::arr<double> &outpntarr) {
    double tstart = mpiMgr.Wtime();
    ++n_todgen;
    levels::arr<double> outpntarr1, outpntarr2;

    if (ntod1 != 0) {
        // Collect northern hemisphere data from outpntarr to outpntarr1
        arrFillingcm(ntod1, timeTest_arr, outpntarr1, outpntarr, 0);
        todTest_arr.alloc(ntod1);
        todTest_arr.fill(0);
    }
    if (ntod2 != 0) {
        // Collect southern hemisphere data from outpntarr to outpntarr2
        arrFillingcm(ntod2, timeTest_arr2, outpntarr2, outpntarr, ntod1);
        todTest_arr2.alloc(ntod2);
        todTest_arr2.fill(0);
    }

    if (CMULT_VERBOSITY > 1) {
        std::cerr << "corenum = " << corenum << "  order = " << order
                  << "  ntod1 = " << ntod1 << "  ntod2 = " << ntod2 << std::endl;
    }

    mpiMgr.barrier();
    if (!pol) {
        interpolTOD_arrTestcm(outpntarr1, outpntarr2,
                                 todTest_arr, todTest_arr2, ntod1, ntod2);
    } else {
        interpolTOD_arrTestcm_pol(outpntarr1, outpntarr2,
                                  todTest_arr, todTest_arr2, ntod1, ntod2);
    }

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


void convolver::conviqt_hemiscm_alltoall(levels::arr3<xcomplex<double> > &tod1,
                                         levels::arr3<xcomplex<double> > &tod2,
                                         long NThetaIndex1,
                                         long NThetaIndex2,
                                         int ithetaoffset1,
                                         int ithetaoffset2) {
    if (CMULT_VERBOSITY > 1) {
        std::cerr << corenum << " : Entering conviqt_hemiscm_alltoall" << std::endl;
    }
    double tstart = mpiMgr.Wtime();
    ++n_conviqt_hemiscm_alltoall;

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
        std::cerr << "conviqt_hemiscm_pol_alltoall :  Out of memory allocating "
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
        wignergen wgen(lmax, rthetas, conv_acc), wgen_neg(lmax, rthetas, conv_acc);
        //wigner_estimator estimator(lmax,100);
        for(long mbeam = 0; mbeam < npsi; ++mbeam) {
            double dsignb = levels::xpow(mbeam, 1);
#pragma omp for schedule(static, 8)
            for (long msky = 0; msky <= lmax; ++msky) {
                const double dsign = levels::xpow(msky, 1);
                const double dsb = dsign * dsignb;
                // estimator.prepare_m(mbeam, msky);
                // if (estimator.canSkip(rthetas[NThetaIndex1-1]))
                //     continue; // negligible dmm
                wgen.prepare(mbeam, msky);
                wgen_neg.prepare(mbeam, -msky);
                for (long lat = 0; lat < my_ntheta; ++lat) {
                    xcomplex<double> Cmm1_pos = my_Cmm1(msky + lmax, mbeam, lat);
                    xcomplex<double> Cmm1_neg = my_Cmm1(-msky + lmax, mbeam, lat);
                    xcomplex<double> Cmm2_pos = my_Cmm2(msky + lmax, mbeam, lat);
                    xcomplex<double> Cmm2_neg = my_Cmm2(-msky + lmax, mbeam, lat);
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

    t_conviqt_hemiscm_alltoall += mpiMgr.Wtime() - tstart;

    if (CMULT_VERBOSITY > 1) {
        std::cerr << corenum << " : Exiting conviqt_hemiscm_alltoall" << std::endl;
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



void convolver::conviqt_hemiscm_pol_alltoall(levels::arr3< xcomplex<double> > &tod1,
                                             levels::arr3< xcomplex<double> > &tod2,
                                             const long NThetaIndex1,
                                             const long NThetaIndex2,
                                             int ithetaoffset1,
                                             int ithetaoffset2) {
    if (CMULT_VERBOSITY > 1) {
        std::cerr << corenum << " : Entering conviqt_hemiscm_pol_alltoall" << std::endl;
    }
    double tstart = mpiMgr.Wtime();
    ++n_conviqt_hemiscm_pol_alltoall;

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
        std::cerr << "conviqt_hemiscm_pol_alltoall :  Out of memory allocating "
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
        t_wigner_init += mpiMgr.Wtime() - t1;
        ++n_wigner_init;
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
                t_wigner_prepare += mpiMgr.Wtime() - t1;
                ++n_wigner_prepare;
                for (long itheta = 0; itheta < my_ntheta; ++itheta) {
                    int firstl1, firstl2;
                    t1 = mpiMgr.Wtime();
                    const levels::arr<double> &dmm = wgen.calc(itheta, firstl1);
                    const levels::arr<double> &dmmneg = wgen_neg.calc(itheta, firstl2);
                    t_wigner_calc += mpiMgr.Wtime() - t1;
                    ++n_wigner_calc;
                    t1 = mpiMgr.Wtime();
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
                    t_lat_iter += mpiMgr.Wtime() - t1;
                    ++n_lat_iter;
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

    t_conviqt_hemiscm_pol_alltoall += mpiMgr.Wtime() - tstart;

    if (CMULT_VERBOSITY > 1) {
        std::cerr << corenum << " : Exiting conviqt_hemiscm_pol_alltoall" << std::endl;
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
            for (long mbeam = 0; mbeam < npsi; ++mbeam) {
                for (long lat = 0; lat < NThetaIndex; ++lat) {
                    const double tmpR = Cmm(ii, mbeam, lat).re;
                    const double tmpI = Cmm(ii, mbeam, lat).im;
                    Cmm(ii, mbeam, lat).re = cs[ii] * tmpR + sn[ii] * tmpI;
                    Cmm(ii, mbeam, lat).im = cs[ii] * tmpI - sn[ii] * tmpR;
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
            for (long mbeam = 0; mbeam < npsi; ++mbeam) {
                const long jj = mbeam + lmax;
                for (long lat = 0; lat < NThetaIndex; ++lat) {
                    double tmpR = tod(ii, mbeam, lat).re;
                    double tmpI = tod(ii, mbeam, lat).im;
                    tod(ii, mbeam, lat).re = cs0[ii] * tmpR + sn0[ii] * tmpI;
                    tod(ii, mbeam, lat).im = cs0[ii] * tmpI - sn0[ii] * tmpR;
                    // Rotate in psi space
                    tmpR = tod(ii, mbeam, lat).re;
                    tmpI = tod(ii, mbeam, lat).im;
                    tod(ii, mbeam, lat).re = cs[jj] * tmpR + sn[jj] * tmpI;
                    tod(ii, mbeam, lat).im = cs[jj] * tmpI - sn[jj] * tmpR;
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


void convolver::interpolTOD_arrTestcm(levels::arr<double> &outpntarr1,
                                         levels::arr<double> &outpntarr2,
                                         levels::arr<double> &TODValue1,
                                         levels::arr<double> &TODValue2,
                                         const long ntod1, const long ntod2) {
    if (CMULT_VERBOSITY > 1) {
        std::cerr << corenum << " : Entered interpolTOD_arrTestcm" << std::endl;
    }
    double tstart = mpiMgr.Wtime();
    ++n_interpolTOD_arrTestcm;
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

    conviqt_hemiscm_alltoall(TODAsym1, TODAsym2,
                             NThetaIndex1, NThetaIndex2,
                             itheta0_1[ntod1 - 1], itheta0_2[ntod2 - 1]);

    for (int thetaIndex = 0; thetaIndex < NThetaIndex1; ++thetaIndex) {
        conviqt_tod_loop(lowerIndex1, upperIndex1, outpntarr1, TODAsym1,
                         thetaIndex, itheta0_1, ntod1, TODValue1, thetaIndex);
    }

    for (int thetaIndex = 0; thetaIndex < NThetaIndex2; ++thetaIndex) {
        conviqt_tod_loop(lowerIndex2, upperIndex2, outpntarr2, TODAsym2,
                         thetaIndex, itheta0_2, ntod2, TODValue2, thetaIndex);
    }

    if (CMULT_VERBOSITY > 1) {
        std::cerr << corenum << " : Exiting interpolTOD_arrTestcm" << std::endl;
    }

    t_interpolTOD_arrTestcm += mpiMgr.Wtime() - tstart;
}


void convolver::interpolTOD_arrTestcm_pol(levels::arr<double> &outpntarr1,
                                          levels::arr<double> &outpntarr2,
                                          levels::arr<double> &TODValue1,
                                          levels::arr<double> &TODValue2,
                                          long ntod1, long ntod2) {
    double tstart = mpiMgr.Wtime();
    ++n_interpolTOD_arrTestcm_pol;

    if (CMULT_VERBOSITY > 1) {
        std::cerr << corenum << " : Entering interpolTOD_arrTestcm_pol" << std::endl;
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

    conviqt_hemiscm_pol_alltoall(TODAsym1, TODAsym2,
                                 NThetaIndex1, NThetaIndex2,
                                 itheta0_1[ntod1 - 1], itheta0_2[ntod2 - 1]);

    for (int thetaIndex = 0; thetaIndex < NThetaIndex1; ++thetaIndex) {
        conviqt_tod_loop_pol(lowerIndex1, upperIndex1, outpntarr1, TODAsym1,
                             thetaIndex, itheta0_1, ntod1, TODValue1,
                             thetaIndex);
    }

    for (int thetaIndex = 0; thetaIndex < NThetaIndex2; ++thetaIndex) {
        conviqt_tod_loop_pol(lowerIndex2, upperIndex2, outpntarr2, TODAsym2,
                             thetaIndex, itheta0_2, ntod2, TODValue2,
                             thetaIndex);
    }

    if (CMULT_VERBOSITY > 1) {
        std::cerr << corenum << " : Exiting interpolTOD_arrTestcm_pol" << std::endl;
    }

    t_interpolTOD_arrTestcm_pol += mpiMgr.Wtime() - tstart;
}


void convolver::conviqt_tod_loop(levels::arr<long> &lowerIndex,
                                 levels::arr<long> &upperIndex,
                                 levels::arr<double> &outpntarr,
                                 levels::arr3<xcomplex<double> > &TODAsym,
                                 const long thetaIndex,
                                 levels::arr<long> &itheta0,
                                 const long ntod,
                                 levels::arr<double> &TODValue, const long lat) {
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

    // Call weight_ncm once unthreaded to make sure the auxiliary arrays are allocated

    std::vector<double> wgt_temp(max_order + 1, 0.);
    weight_ncm(0., wgt_temp);

#pragma omp parallel default(shared)
    {
        std::vector<double> cosang(npsi), sinang(npsi);
        std::vector<double> wgt1(max_order + 1, 0.);
        std::vector<double> wgt2(max_order + 1, 0.);
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
            const double weight1 = 2 * wgt1[thetaIndex - (itheta0[ii] - itheta0[ntod - 1])];
            for (int phiIndex = 0; phiIndex < npoints; ++phiIndex) {
                const double weight = wgt2[phiIndex] * weight1;
                long newphiIndex = iphi0 + phiIndex;
                if (newphiIndex >= nphi) {
                    newphiIndex -= nphi;
                }
                xcomplex<double> *ca = &(conviqtarr[newphiIndex][0]);
                double x = 0;
                for (long ipsi = 0; ipsi < npsi; ++ipsi) {
                    x += cosang[ipsi] * (*ca).re - sinang[ipsi] * (*ca).im;
                    ++ca;
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


void convolver::conviqt_tod_loop_pol(levels::arr<long> &lowerIndex,
                                     levels::arr<long> &upperIndex,
                                     levels::arr<double> &outpntarr,
                                     levels::arr3<xcomplex<double> > &TODAsym,
                                     long thetaIndex,
                                     levels::arr<long> &itheta0,
                                     long ntod,
                                     levels::arr<double> &TODValue,
                                     long lat) {
    if (CMULT_VERBOSITY > 1) {
        std::cerr << corenum << " : Entering conviqt_tod_loop_pol" << std::endl;
    }
    double tstart = mpiMgr.Wtime();
    ++n_conviqt_tod_loop_pol;
    levels::arr2< xcomplex<double> > conviqtarr;
    try {
        conviqtarr.alloc(nphi, npsi);
    } catch (std::bad_alloc &e) {
        std::cerr << "conviqt_tod_loop_pol : Out of memory allocating "
                  << nphi * npsi * 16. / 1024 / 1024 << "MB for conviqtarr"
                  << std::endl;
        throw;
    }

    for (long ii = 0; ii < nphi; ++ii) {
        for (long jj = 0; jj < npsi; ++jj) {
            conviqtarr[ii][jj] = TODAsym(ii, jj, lat);
        }
    }

    // Call weight_ncm once unthreaded to make sure the auxiliary arrays are allocated

    std::vector<double> wgt_temp(max_order + 1, 0.);
    weight_ncm(0., wgt_temp);

#pragma omp parallel default(shared)
    {
        std::vector<double> cosang(npsi), sinang(npsi);
        std::vector<double> wgt1(max_order + 1, 0.);
        std::vector<double> wgt2(max_order + 1, 0.);
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
            double weight1 = 2 * wgt1[thetaIndex - (itheta0[ii] - itheta0[ntod - 1])];
            for (int iphi = 0; iphi < npoints; ++iphi) {
                double weight = wgt2[iphi] * weight1;
                long iphinew = iphi0 + iphi;
                if (iphinew >= nphi) {
                    iphinew -= nphi;
                }
                xcomplex<double> *ca = &(conviqtarr[iphinew][0]);
                double x = 0;
                for (long ipsi = 0; ipsi < npsi; ++ipsi) {
                    x += cosang[ipsi] * (*ca).re - sinang[ipsi] * (*ca).im;
                    ++ca;
                }
                TODValue[ii] += weight * x;
            }
        }
    }  // end parallel region

    t_conviqt_tod_loop_pol += mpiMgr.Wtime() - tstart;
    if (CMULT_VERBOSITY > 1) {
        std::cerr << corenum << " : Exiting conviqt_tod_loop_pol" << std::endl;
    }
}


void convolver::arrFillingcm(long ntod,
                             levels::arr<double> &timeTest_arr,
                             levels::arr<double> &outpntarrx,
                             levels::arr<double> &outpntarr,
                             long offindex) {
    double tstart = mpiMgr.Wtime();
    ++n_arrFillingcm;
    timeTest_arr.alloc(ntod);
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
    for (long ii = 0; ii < ntod; ++ii)
        timeTest_arr[ii] = outpntarrx[5 * ii + 4];

    t_arrFillingcm += mpiMgr.Wtime() - tstart;
}


int convolver::convolve(pointing &pntarr, bool calibrate) {

    double tstart = mpiMgr.Wtime();
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

    // DEBUG begin
    /*
    if (mpiMgr.master()) {
        for (int corenum=0; corenum < cores; ++corenum) {
            double theta = corethetaarr[corenum];
            std::cerr << corenum << " : theta = " << theta * 180 / pi
                      << std::endl;
        }
    }
    */
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

    todgen(ntod1, ntod2, todTest_arr, timeTest_arr,
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
            maxtodall = std::max(maxtodall, outtodarr[ii]);
            mintodall = std::min(mintodall, outtodarr[ii]);
        }

        if (totsize > 0)
            std::cerr << "maxtodall = " << maxtodall  << "   mintodall = "
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
                maxtoddiff = std::max(fabs(todtmp[ii]-pntarr[5 * ii + 3]), maxtoddiff);
                if (ii % 100000 == 0) {
                    std::cerr << "todtmp[ii] = " << todtmp[ii] << "  pntarr[5*ii+3] = "
                              << pntarr[5 * ii + 3] << "  difference = "
                              << abs(todtmp[ii] - pntarr[5 * ii + 3]) << std::endl;
                }
            }
        }
        if (CMULT_VERBOSITY > 1) {
            std::cerr << "  corenum = " << corenum << "   maxtoddiff = "
                      << maxtoddiff << "   calibration = " << calibration << std::endl;
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
    timing_line(std::string("    todRedistribution5cm"),
                t_todRedistribution5cm, n_todRedistribution5cm);
    timing_line(std::string("    todgen"), t_todgen, n_todgen);
    timing_line(std::string("        arrFillingcm"), t_arrFillingcm, n_arrFillingcm);
    timing_line(std::string("        interpolTOD_arrTestcm"),
                t_interpolTOD_arrTestcm, n_interpolTOD_arrTestcm);
    timing_line(std::string("        interpolTOD_arrTestcm_pol"),
                t_interpolTOD_arrTestcm_pol, n_interpolTOD_arrTestcm_pol);
    timing_line(std::string("            itheta0SetUp"), t_itheta0SetUp, n_itheta0SetUp);
    timing_line(std::string("                ithetacalc"), t_ithetacalc, n_ithetacalc);
    timing_line(std::string("            conviqt_hemiscm_alltoall"),
                t_conviqt_hemiscm_alltoall, n_conviqt_hemiscm_alltoall);
    timing_line(std::string("            conviqt_hemiscm_pol_alltoall"),
                t_conviqt_hemiscm_pol_alltoall, n_conviqt_hemiscm_pol_alltoall);
    timing_line(std::string("                wigner_init"), t_wigner_init, n_wigner_init);
    timing_line(std::string("                wigner_prepare"),
                t_wigner_prepare, n_wigner_prepare);
    timing_line(std::string("                wigner_calc"), t_wigner_calc, n_wigner_calc);
    timing_line(std::string("                lat_iter"), t_lat_iter, n_lat_iter);
    timing_line(std::string("                todAnnulus"),  t_todAnnulus, n_todAnnulus);
    timing_line(std::string("                    alltoall_datacube"),
                t_alltoall_datacube, n_alltoall_datacube);
    timing_line(std::string("            conviqt_tod_loop"),
                t_conviqt_tod_loop, n_conviqt_tod_loop);
    timing_line(std::string("            conviqt_tod_loop_pol"),
                t_conviqt_tod_loop_pol, n_conviqt_tod_loop_pol);
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
