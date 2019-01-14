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

    t_todgen_v4 = 0;
    t_arrFillingcm_v2 = 0;
    t_interpolTOD_arrTestcm_v4 = 0;
    t_itheta0SetUp = 0;
    t_ithetacalc = 0;
    t_conviqt_tod_loop_v4 = 0;
    t_weight_ncm = 0;
    t_conviqt_hemiscm_single = 0;
    t_todAnnulus_v3 = 0;
    t_conviqt_hemiscm_v4 = 0;
    t_interpolTOD_arrTestcm_pol_v4 = 0;
    t_conviqt_hemiscm_pol_v4 = 0;
    t_conviqt_hemiscm_pol_single = 0;
    t_conviqt_tod_loop_pol_v5 = 0;

    t_wigner_init = 0;
    t_wigner_prepare = 0;
    t_wigner_calc = 0;

    t_lat_iter = 0;
    t_sincos_iter = 0;    
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


void convolver::weight_ncm(double x, levels::arr<double> &wgt) {
    unsigned int npoints = order + 1;

    if (base_wgt.size() == 0) {
        // Initialize base_wgt only when needed
        base_wgt.resize(npoints, 1);
        for (int m = 0; m < (int)npoints; ++m) {
            for (int n = 0; n<m; ++n)
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
}


void convolver::weight_ncm(double x, std::vector<double> &wgt) {
    // This method is called from multiple threads at once
    unsigned int npoints = order + 1;

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
}


void convolver::ithetacalc(levels::arr<long> &itheta0,
                           levels::arr<double> &outpntarr,
                           long ntod, double inv_delta_theta,
                           double theta0, long ioffset,
                           long ntheta, long npoints) {
    double tstart = mpiMgr.Wtime();
    for (long ii = 0; ii < ntod; ++ii) {
        double beta = outpntarr[5 * ii + 1]; // latitude
        beta = (beta > halfpi) ? pi - beta : beta;
        // Note that the larger the ii the smaller frac is and the smaller itheta0[ii] is.
        double frac = (beta - theta0) * inv_delta_theta;
        //Note also that itheta0[ii] is always positive.
        itheta0[ii] = int(frac) - ioffset;
        if (itheta0[ii] > (ntheta - npoints))
            itheta0[ii] = ntheta - npoints;
        if (itheta0[ii] < 0)
            itheta0[ii] = 0;
    }
    t_ithetacalc = mpiMgr.Wtime() - tstart;
}


void convolver::ratiobetacalcgreatercm(long &ratiobeta, double theta,
                                       levels::arr<double> &corethetaarr) {
    for (long ii = 0; ii < cores; ++ii) {
        if (ratiobeta + 1 >= cores)
            break;
        if (theta < corethetaarr[ratiobeta + 1])
            break;
        ++ratiobeta;
    }
}


void convolver::ratiobetacalcsmallercm(long &ratiobeta, double theta,
                                       levels::arr<double> &corethetaarr) {
    for (long ii = 0; ii < 100000; ++ii) {
        if (ratiobeta < 0) break;
        if (theta > corethetaarr[ratiobeta])
            break;
        --ratiobeta;
    }
}

void convolver::thetaDeltaThetacm(double thetaini,
                                  int corenum,
                                  double &theta,
                                  double &deltatheta) {
    if (corenum == 0) {
        theta = 0;
        deltatheta = thetaini;
    } else if (corenum == cores-1) {
      deltatheta = 1. / cores;
      theta = thetaini - deltatheta;
    } else {
      double a = -0.5 * cos(thetaini);
      double b = -sin(thetaini) + thetaini * cos(thetaini);
      double c = (thetaini * sin(thetaini) -
                  cos(thetaini) * thetaini * thetaini / 2. -
                  1. / cores);
      theta = (-b - sqrt(b * b - 4. * a * c)) / (2. * a);
      deltatheta = thetaini - theta;
    }
}


void convolver::deltaTheta2(long iival,
                            double thetaini,
                            levels::arr<double> &dbeta) {
    dbeta.alloc(iival + 1);
    double dy = (1 - thetaini) / 2.;
    double yshift = -0.7;
    double theta2 = pi - abs(asin(1 - dy)) + yshift;
    double theta1 = pi - abs(asin(sin(theta2) - thetaini));
    double dt = (theta2 - theta1) / (iival + 1.);
    double dbeta0 = sin(theta2);

    if (CMULT_VERBOSITY > 1)
        std::cout << "dy = " << dy << "   theta1 = " << theta1
                  << "   theta2 = " << theta2
                  << "   dbeta0 = " << dbeta0
                  << "  dt = " << dt
                  << "   iival = " << iival
                  << "  sin(theta2-(iival+1-iival)*dt) = "
                  << sin(theta2 - (iival + 1 - iival) * dt)
                  << "  sin(theta2-(iival+1-(iival-1))*dt) = "
                  << sin(theta2 - (iival + 1 - (iival - 1)) * dt) << std::endl;

  for (long ii = iival; ii >= 0; --ii) {
      dbeta[ii] = dbeta0 - sin(theta2 - (iival + 1 - ii) * dt);
      dbeta0 = sin(theta2 - (iival + 1 - ii) * dt);
    }
}


void convolver::fillingBetaSeg (levels::arr<double> & pntarr,
                                long & arrsize, double ratiodeltas,
                                levels::arr<double> &corethetaarr,
                                levels::arr<int> &inBetaSeg) {
    if (CMULT_VERBOSITY > 1)
        std::cout << "Entered fillingBetaSeg in core = " << corenum << std::endl;

  for (long ii = 0; ii < arrsize; ++ii) {
      double theta = pntarr[5 * ii + 1];
      if (theta >= halfpi)
          theta = pi - theta;
      long ratiobeta;
      ratiobeta = theta / ratiodeltas;
      if (theta >= corethetaarr[ratiobeta])
          ratiobetacalcgreatercm(ratiobeta, theta, corethetaarr);
      else
          ratiobetacalcsmallercm(ratiobeta, theta, corethetaarr);
      inBetaSeg[ratiobeta] += 5; // Each entry comprises 5 elements
    }

  if (CMULT_VERBOSITY > 1)
      std::cout << "Leaving fillingBetaSeg in core = " << corenum << std::endl;
}


void convolver::todRedistribution5cm (levels::arr<double> pntarr,
                                      levels::arr<int> inBetaSeg,
                                      levels::arr<int> outBetaSeg,
                                      levels::arr<int> &inBetaSegAcc,
                                      levels::arr<int> &outBetaSegAcc,
                                      long &outBetaSegSize,
                                      levels::arr<double> &pntarr2,
                                      long totsize,
                                      levels::arr<int> &inOffset,
                                      levels::arr<int> &outOffset,
                                      double ratiodeltas,
                                      levels::arr<double> &corethetaarr) {
    if (CMULT_VERBOSITY > 1) std::cout << "Entered todRedistribution5cm" << std::endl;

    inBetaSegAcc.alloc(cores);
    inBetaSegAcc.fill(0);
    outBetaSegAcc.alloc(cores);
    outBetaSegAcc.fill(0);
    for (long ii = 0; ii < cores; ++ii)
        for (long jj = 0; jj <= ii; ++jj) {
            inBetaSegAcc[ii] += inBetaSeg[jj];
            outBetaSegAcc[ii] += outBetaSeg[jj];
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
      pntarr2.fill(0);
      levels::arr<int> countbeta(cores, 0);
      for (long ii = 0; ii < totsize; ++ii) {
	  long ratiobeta; //, nbeta;
	  double theta = pntarr[5 * ii + 1];
	  if (theta < 0 || theta > pi)
	      throw std::runtime_error("ERROR: Illegal latitude in pntarr: theta = " +
                                       std::to_string(theta) + " not in 0..pi");
	  if (theta >= halfpi)
              theta = pi - theta;
	  ratiobeta = theta / ratiodeltas;
	  if (theta >= corethetaarr[ratiobeta])
              ratiobetacalcgreatercm(ratiobeta, theta, corethetaarr);
	  else
              ratiobetacalcsmallercm(ratiobeta, theta, corethetaarr);

          long i1 = 5 * ii;
          long i2 = inOffset[ratiobeta] + countbeta[ratiobeta];
          for (long k = 0; k < 5; ++k)
              pntarr2[i2 + k] = pntarr[i1 + k];
	  countbeta[ratiobeta] += 5;
	}
    }

  if (CMULT_VERBOSITY > 1)
      std::cout << "Leaving todRedistribution5cm" << std::endl;
}


void convolver::preReorderingStep(long ntod, long ntod2,
                                  levels::arr<double> &todAll,
                                  levels::arr<double> &todTest_arr,
                                  levels::arr<double> &todTest_arr2) {
  if (ntod + ntod2 != 0) todAll.alloc(ntod + ntod2);

  if (ntod != 0) {
      for (long ii = 0; ii < ntod; ++ii)
          todAll[ii] = todTest_arr[ii];
      todTest_arr.dealloc();
  }

  if (ntod2 != 0) {
      for (long ii = 0; ii < ntod2; ++ii)
          todAll[ntod + ii] = todTest_arr2[ii];
      todTest_arr2.dealloc();
  }

}


void convolver::todgen_v4(long ntod, long ntod2,
                          levels::arr<double> &todTest_arr,
                          levels::arr<double> &timeTest_arr,
                          levels::arr<double> &todTest_arr2,
                          levels::arr<double> &timeTest_arr2,
                          levels::arr<double> &outpntarr) {
    double tstart = mpiMgr.Wtime();
    levels::arr<double> outpntarr1, outpntarr2;

    if (ntod != 0) {
        // Collect northern hemisphere data from outpntarr to outpntarr1
        arrFillingcm_v2(ntod, timeTest_arr, outpntarr1, outpntarr, 0);
        todTest_arr.alloc(ntod);
        todTest_arr.fill(0);
    }
    if (ntod2 != 0) {
        // Collect southern hemisphere data from outpntarr to outpntarr2
        arrFillingcm_v2(ntod2, timeTest_arr2, outpntarr2, outpntarr, ntod);
        todTest_arr2.alloc(ntod2);
        todTest_arr2.fill(0);
    }

    if (CMULT_VERBOSITY > 1)
        std::cout << "corenum = " << corenum << "  order = " << order
                  << "  ntod = " << ntod << "  ntod2 = " << ntod2 << std::endl;

    if (!pol)
        interpolTOD_arrTestcm_v4(outpntarr1, outpntarr2,
                                 todTest_arr, todTest_arr2, ntod, ntod2);
    else
        interpolTOD_arrTestcm_pol_v4(outpntarr1, outpntarr2,
                                     todTest_arr, todTest_arr2, ntod, ntod2);

    t_todgen_v4 = mpiMgr.Wtime() - tstart;
}


void convolver::itheta0SetUp(int npoints, int ioffset, long ntheta, double theta0,
                             double inv_delta_theta, long nphi,
                             levels::arr<double> outpntarr, long ntod,
                             long &NThetaIndex,
                             levels::arr<long> &itheta0,
                             levels::arr<long> &lowerIndex,
                             levels::arr<long> &upperIndex,
                             levels::arr3<xcomplex<double> > &TODAsym) {
    double tstart = mpiMgr.Wtime();
    itheta0.alloc(ntod);
    ithetacalc(itheta0, outpntarr, ntod, inv_delta_theta, theta0, ioffset,
               ntheta, npoints);
    NThetaIndex = itheta0[0] - itheta0[ntod - 1] + npoints;
    lowerIndex.alloc(NThetaIndex);
    lowerIndex.fill(ntod);
    upperIndex.alloc(NThetaIndex);
    upperIndex.fill(0);
    for (long jj = 0; jj < NThetaIndex; ++jj)
        for (long ii = ntod - 1; ii >= 0; --ii)
            if (itheta0[ii] > itheta0[ntod - 1] + jj - npoints) {
                lowerIndex[jj] = ii;
                break;
            }
    for (long jj = 0; jj < NThetaIndex; ++jj)
        for (long ii = 0; ii < ntod; ++ii)
            if (itheta0[ii] <= itheta0[ntod - 1] + jj) {
                // this yields the smallest ii where
                // itheta0[ii] <= itheta0[ntod - 1] + jj
                upperIndex[jj] = ii;
                break;
            }
    try {
        TODAsym.alloc(nphi, beammmax + 1, NThetaIndex);
    } catch (std::bad_alloc & e) {
        std::cerr << "itheta0SetUp :  Out of memory allocating "
                  << nphi * (beammmax + 1) * NThetaIndex * 8. / 1024 / 1024
                  << "MB for TODasym" << std::endl;
        throw;
    }

    TODAsym.fill(0.);
    t_itheta0SetUp = mpiMgr.Wtime() - tstart;
}


void convolver::conviqt_hemiscm_v4(levels::arr3<xcomplex<double> > &tod1,
                                   levels::arr3<xcomplex<double> > &tod2,
                                   long NThetaIndex1,
                                   levels::arr<double> &rthetas) {
    double tstart = mpiMgr.Wtime();
    long binElements = 2 * lmax + 1;
    long psiElements = beammmax + 1;
    levels::arr3<xcomplex<double> > Cmm(binElements, psiElements, NThetaIndex1);
    Cmm.fill(0.);
    levels::arr3<xcomplex<double> > Cmm2(binElements, psiElements, NThetaIndex1);
    Cmm2.fill(0.);

    Alm< xcomplex<float> > & blmT = b->blmT();
    Alm< xcomplex<float> > & slmT = s->slmT();
    levels::arr<double> cs(2 * lmax + 1), sn(2 * lmax + 1);
    levels::arr<double> cs0(binElements), sn0(binElements);

#pragma omp parallel default(shared)
    {
        wignergen wgen(lmax, rthetas, conv_acc), wgen_neg(lmax, rthetas, conv_acc);
        //wigner_estimator estimator(lmax,100);
        for(long beamIndex = 0; beamIndex <= beammmax; ++beamIndex) {
            double dsignb = levels::xpow(beamIndex, 1);
#pragma omp for schedule(static, 8)
            for (long msky = 0; msky <= lmax; ++msky) {
                double dsign = levels::xpow(msky, 1);
                double dsb = dsign * dsignb;
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
                    int firstl = (firstl1 > firstl2) ? firstl1 : firstl2;
                    double dlb = -levels::xpow(firstl, dsignb);
                    for (long ii = firstl; ii <= lmax; ++ii) {
                        dlb = -dlb;
                        // Note that msky in dlm is located to the left of
                        // mb and they have to be interchanged in the convolution
                        double dMatrixElementmskypos = dsb * dmm[ii];
                        double dMatrixElementmskyneg = dsb * dmmneg[ii];
                        double dMatrixElementmskypos2 = dlb * dMatrixElementmskyneg;
                        double dMatrixElementmskyneg2 = dlb * dMatrixElementmskypos;
                        xcomplex<float> sT = slmT(ii, msky);
                        xcomplex<float> bT = blmT(ii, beamIndex);
                        double prod1 = sT.re * bT.re;
                        double prod3 = sT.im * bT.re;
                        double prod2 = sT.im * bT.im;
                        double prod4 = sT.re * bT.im;
                        double tmp_1 = prod1 + prod2;
                        double tmp_2 = prod3 - prod4;
                        double xtmp_1 = tmp_1 * dMatrixElementmskypos;
                        double xtmp_2 = tmp_2 * dMatrixElementmskypos;
                        double xtmp_5 = tmp_1 * dMatrixElementmskypos2;
                        double xtmp_6 = tmp_2 * dMatrixElementmskypos2;
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

    t_conviqt_hemiscm_v4 = mpiMgr.Wtime() - tstart;
}


void convolver::conviqt_hemiscm_pol_v4(levels::arr3< xcomplex<double> > &tod1,
                                       levels::arr3< xcomplex<double> > &tod2,
                                       long NThetaIndex1,
                                       levels::arr<double> &rthetas) {
    double tstart = mpiMgr.Wtime();
    long binElements = 2 * lmax + 1;
    long psiElements = beammmax + 1;
    levels::arr3< xcomplex<double> > Cmm(binElements, psiElements, NThetaIndex1);
    Cmm.fill(0.);
    levels::arr3< xcomplex<double> > Cmm2(binElements, psiElements, NThetaIndex1);
    Cmm2.fill(0.);

    Alm< xcomplex<float> > & blmT = b->blmT();
    Alm< xcomplex<float> > & blmG = b->blmG();
    Alm< xcomplex<float> > & blmC = b->blmC();
    Alm< xcomplex<float> > & slmT = s->slmT();
    Alm< xcomplex<float> > & slmG = s->slmG();
    Alm< xcomplex<float> > & slmC = s->slmC();
    levels::arr<double> cs(2 * lmax + 1), sn(2 * lmax + 1);
    levels::arr<double> cs0(binElements), sn0(binElements);

#pragma omp parallel default(shared)
    {
        double t1 = mpiMgr.Wtime();
        wignergen wgen(lmax, rthetas, conv_acc), wgen_neg(lmax, rthetas, conv_acc);
        t_wigner_init += mpiMgr.Wtime() - t1;
        // wigner_estimator estimator(lmax, 100);
        for(long beamIndex = 0; beamIndex <= beammmax; ++beamIndex) {
            double dsignb = levels::xpow(beamIndex, 1);
#pragma omp for schedule(static, 8)
            for (long msky = 0; msky <= lmax; ++msky) {
                double dsign = levels::xpow(msky, 1);
                double dsb = dsign * dsignb;
                //estimator.prepare_m(beamIndex, msky);
                // if (estimator.canSkip(rthetas[NThetaIndex1-1]))
                //     continue; // negligible dmm
                t1 = mpiMgr.Wtime();        
                wgen.prepare(beamIndex, msky);
                wgen_neg.prepare(beamIndex, -msky);
                t_wigner_prepare += mpiMgr.Wtime() - t1;
                for (long lat = 0; lat < NThetaIndex1; ++lat) {
                    int firstl1, firstl2;
                    t1 = mpiMgr.Wtime();        
                    const levels::arr<double> &dmm = wgen.calc(lat, firstl1);
                    const levels::arr<double> &dmmneg = wgen_neg.calc(lat, firstl2);
                    t_wigner_calc += mpiMgr.Wtime() - t1;
                    t1 = mpiMgr.Wtime();
                    xcomplex<double> Cmm_pos = Cmm(msky + lmax, beamIndex, lat);
                    xcomplex<double> Cmm_neg = Cmm(-msky + lmax, beamIndex, lat);
                    xcomplex<double> Cmm2_pos = Cmm2(msky + lmax, beamIndex, lat);
                    xcomplex<double> Cmm2_neg = Cmm2(-msky + lmax, beamIndex, lat);
                    int firstl = (firstl1 > firstl2) ? firstl1 : firstl2;
                    double dlb = -levels::xpow(firstl, dsignb);
                    for (long ii = firstl; ii <= lmax; ++ii) {
                        dlb = -dlb;
                        // Note that msky in dlm is located to the
                        // left of mb and they have to be interchanged in the convolution
                        double dMatrixElementmskypos = dsb * dmm[ii];
                        double dMatrixElementmskyneg = dsb * dmmneg[ii];
                        double dMatrixElementmskypos2 = dlb * dMatrixElementmskyneg;
                        double dMatrixElementmskyneg2 = dlb * dMatrixElementmskypos;
                        xcomplex<float> sT = slmT(ii, msky);
                        xcomplex<float> sG = slmG(ii, msky);
                        xcomplex<float> sC = slmC(ii, msky);
                        xcomplex<float> bT = blmT(ii, beamIndex);
                        xcomplex<float> bG = blmG(ii, beamIndex);
                        xcomplex<float> bC = blmC(ii, beamIndex);
                        double prod1 = sT.re * bT.re + sG.re * bG.re + sC.re * bC.re;
                        double prod3 = sT.im * bT.re + sG.im * bG.re + sC.im * bC.re;
                        double prod2 = sT.im * bT.im + sG.im * bG.im + sC.im * bC.im;
                        double prod4 = sT.re * bT.im + sG.re * bG.im + sC.re * bC.im;
                        double tmp_1 = prod1 + prod2;
                        double tmp_2 = prod3 - prod4;
                        double xtmp_1 = tmp_1 * dMatrixElementmskypos;
                        double xtmp_2 = tmp_2 * dMatrixElementmskypos;
                        double xtmp_5 = tmp_1 * dMatrixElementmskypos2;
                        double xtmp_6 = tmp_2 * dMatrixElementmskypos2;
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
    } // End of parallel section

  todAnnulus_v3(tod1, Cmm, cs, sn, cs0, sn0, NThetaIndex1);
  todAnnulus_v3(tod2, Cmm2, cs, sn, cs0, sn0, NThetaIndex1);

  t_conviqt_hemiscm_pol_v4 = mpiMgr.Wtime() - tstart;
}


void convolver::conviqt_hemiscm_single(levels::arr3<xcomplex<double> > &tod1,
                                       long NThetaIndex1,
                                       levels::arr<double> &rthetas) {
    double tstart = mpiMgr.Wtime();
    long binElements = 2 * lmax + 1;
    long psiElements = beammmax + 1;
    levels::arr3<xcomplex<double> > Cmm(binElements, psiElements, NThetaIndex1);
    Cmm.fill(0.);

    Alm< xcomplex<float> > & blmT = b->blmT();
    Alm< xcomplex<float> > & slmT = s->slmT();
    levels::arr<double> cs(2 * lmax + 1), sn(2 *lmax + 1);
    levels::arr<double> cs0(binElements), sn0(binElements);

#pragma omp parallel default(shared)
    {
        wignergen wgen(lmax, rthetas, conv_acc), wgen_neg(lmax, rthetas, conv_acc);
        // wigner_estimator estimator(lmax, 100);

        for(long beamIndex = 0; beamIndex <= beammmax; ++beamIndex) {
            double dsignb = levels::xpow(beamIndex, 1);
#pragma omp for schedule(static,8)
            for (long msky = 0; msky <= lmax; ++msky) {
                double dsign = levels::xpow(msky, 1);
                double dsb = dsign * dsignb;
                // estimator.prepare_m(beamIndex, msky);
                // if (estimator.canSkip(rthetas[NThetaIndex1 - 1]))
                //     continue; // negligible dmm
                wgen.prepare(beamIndex, msky);
                wgen_neg.prepare(beamIndex, -msky);
                for (long lat = 0; lat < NThetaIndex1; ++lat) {
                    xcomplex<double> Cmm_pos = Cmm(msky + lmax, beamIndex, lat);
                    xcomplex<double> Cmm_neg = Cmm(-msky + lmax, beamIndex, lat);
                    int firstl1, firstl2;
                    const levels::arr<double> &dmm = wgen.calc(lat, firstl1);
                    const levels::arr<double> &dmmneg = wgen_neg.calc(lat, firstl2);
                    int firstl = (firstl1 > firstl2) ? firstl1 : firstl2;
                    for (long ii = firstl; ii <= lmax; ++ii) {
                        // Note that msky in dlm is located to the left of
                        // mb and they have to be interchanged in the convolution
                        double dMatrixElementmskypos = dsb * dmm[ii];
                        double dMatrixElementmskyneg = dsb * dmmneg[ii];
                        xcomplex<float> sT = slmT(ii, msky);
                        xcomplex<float> bT = blmT(ii, beamIndex);
                        double prod1 = sT.re * bT.re;
                        double prod3 = sT.im * bT.re;
                        double prod2 = sT.im * bT.im;
                        double prod4 = sT.re * bT.im;
                        double tmp_1 = prod1 + prod2;
                        double tmp_2 = prod3 - prod4;
                        double xtmp_1 = tmp_1 * dMatrixElementmskypos;
                        double xtmp_2 = tmp_2 * dMatrixElementmskypos;
                        Cmm_pos.re += xtmp_1;
                        Cmm_pos.im += xtmp_2;
                        if (msky != 0) {
                            double tmp_3 = dsign * (prod1 - prod2);
                            double tmp_4 = -dsign * (prod3 + prod4);
                            double xtmp_3 = tmp_3 * dMatrixElementmskyneg;
                            double xtmp_4 = tmp_4 * dMatrixElementmskyneg;
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

    todAnnulus_v3(tod1, Cmm, cs, sn, cs0, sn0, NThetaIndex1);

    t_conviqt_hemiscm_single = mpiMgr.Wtime() - tstart;
}


void convolver::conviqt_hemiscm_pol_single(levels::arr3<xcomplex<double> > &tod1,
                                           long NThetaIndex1,
                                           levels::arr<double> &rthetas) {
    double tstart = mpiMgr.Wtime();
    long binElements = 2 * lmax + 1;
    long psiElements = beammmax + 1;
    levels::arr3<xcomplex<double> > Cmm(binElements, psiElements, NThetaIndex1);
    Cmm.fill(0.);

    Alm< xcomplex<float> > & blmT = b->blmT();
    Alm< xcomplex<float> > & blmG = b->blmG();
    Alm< xcomplex<float> > & blmC = b->blmC();
    Alm< xcomplex<float> > & slmT = s->slmT();
    Alm< xcomplex<float> > & slmG = s->slmG();
    Alm< xcomplex<float> > & slmC = s->slmC();
    levels::arr<double> cs(2 * lmax + 1), sn(2 * lmax + 1);
    levels::arr<double> cs0(binElements), sn0(binElements);

#pragma omp parallel default(shared)
    {
        double t1 = mpiMgr.Wtime();
        wignergen wgen(lmax, rthetas, conv_acc), wgen_neg(lmax, rthetas, conv_acc);
        t_wigner_init += mpiMgr.Wtime() - t1;
        // wigner_estimator estimator(lmax, 100);
        for (long beamIndex = 0; beamIndex <= beammmax; ++beamIndex) {
            double dsignb = levels::xpow(beamIndex, 1);
#pragma omp for schedule(static, 8)
            for (long msky = 0; msky <= lmax; ++msky) {
                double dsign = levels::xpow(msky, 1);
                double dsb = dsign * dsignb;
                //estimator.prepare_m(beamIndex, msky);
                //if (estimator.canSkip(rthetas[NThetaIndex1 - 1]))
                //    continue; // negligible dmm
                t1 = mpiMgr.Wtime();
                wgen.prepare(beamIndex, msky);
                wgen_neg.prepare(beamIndex, -msky);
                t_wigner_prepare += mpiMgr.Wtime() - t1;
                for (long lat = 0; lat < NThetaIndex1; ++lat) {
                    int firstl1, firstl2;
                    t1 = mpiMgr.Wtime();
                    const levels::arr<double> &dmm = wgen.calc(lat, firstl1);
                    const levels::arr<double> &dmmneg = wgen_neg.calc(lat, firstl2);
                    t_wigner_calc += mpiMgr.Wtime() - t1;
                    t1 = mpiMgr.Wtime();
                    xcomplex<double> Cmm_pos = Cmm(msky + lmax, beamIndex, lat);
                    xcomplex<double> Cmm_neg = Cmm(-msky + lmax, beamIndex, lat);
                    int firstl = (firstl1 > firstl2) ? firstl1 : firstl2;
                    for (long ii = firstl; ii <= lmax; ++ii) {
                        // Note that msky in dlm is located to the left
                        // of mb and they have to be interchanged in the convolution
                        double dMatrixElementmskypos = dsb * dmm[ii];
                        xcomplex<float> sT = slmT(ii, msky);
                        xcomplex<float> sG = slmG(ii, msky);
                        xcomplex<float> sC = slmC(ii, msky);
                        xcomplex<float> bT = blmT(ii, beamIndex);
                        xcomplex<float> bG = blmG(ii, beamIndex);
                        xcomplex<float> bC = blmC(ii, beamIndex);
                        double prod1 = sT.re * bT.re + sG.re * bG.re + sC.re * bC.re;
                        double prod3 = sT.im * bT.re + sG.im * bG.re + sC.im * bC.re;
                        double prod2 = sT.im * bT.im + sG.im * bG.im + sC.im * bC.im;
                        double prod4 = sT.re * bT.im + sG.re * bG.im + sC.re * bC.im;
                        double tmp_1 = prod1 + prod2;
                        double tmp_2 = prod3 - prod4;
                        double xtmp_1 = tmp_1 * dMatrixElementmskypos;
                        double xtmp_2 = tmp_2 * dMatrixElementmskypos;
                        Cmm_pos.re += xtmp_1;
                        Cmm_pos.im += xtmp_2;
                        if (msky != 0) {
                            double tmp_3 = prod1 - prod2;
                            double tmp_4 = -prod3 - prod4;
                            double dMatrixElementmskyneg = dsb * dmmneg[ii];
                            double xtmp_3 = dsign * tmp_3 * dMatrixElementmskyneg;
                            double xtmp_4 = dsign * tmp_4 * dMatrixElementmskyneg;
                            Cmm_neg.re += xtmp_3;
                            Cmm_neg.im += xtmp_4;
                        }
                    }
                    t_lat_iter += mpiMgr.Wtime() - t1;
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
    } // end of parallel region

    todAnnulus_v3(tod1, Cmm, cs, sn, cs0, sn0, NThetaIndex1);

    t_conviqt_hemiscm_pol_single = mpiMgr.Wtime() - tstart;
}


void convolver::todAnnulus_v3(levels::arr3<xcomplex<double> > &tod1,
                              levels::arr3<xcomplex<double> > &Cmm,
                              levels::arr<double> &cs,
                              levels::arr<double> &sn,
                              levels::arr<double> &cs0,
                              levels::arr<double> &sn0,
                              long NThetaIndex1) {
    double tstart = mpiMgr.Wtime();
    long binElements = 2 * lmax + 1;

#pragma omp parallel default(shared)
    {
#pragma omp for schedule(static,8)
        for (long msky = -lmax; msky <= lmax; ++msky) {
            for (long beamIndex = 0; beamIndex <= beammmax; ++beamIndex) {
                for (long lat = 0; lat < NThetaIndex1; ++lat) {
                    //double dPhi = -halfpi;
                    double tmpR = Cmm(msky + lmax, beamIndex, lat).real();
                    double tmpI = Cmm(msky + lmax, beamIndex, lat).imag();
                    Cmm(msky + lmax, beamIndex, lat).real()
                        = (cs[msky + lmax] * tmpR + sn[msky + lmax] * tmpI);
                    Cmm(msky + lmax, beamIndex, lat).imag()
                        = (cs[msky + lmax] * tmpI - sn[msky + lmax] * tmpR);
                }
            }
        } // end of parallel for
#pragma omp for schedule(static, 1)
        for (long beamIndex = 0; beamIndex <= beammmax; ++beamIndex) {
            for (long lat = 0; lat < NThetaIndex1; ++lat) {
                levels::arr<xcomplex<double> > Cmsky(binElements, 0);
                for (long msky = -lmax; msky <= lmax; ++msky)
                    Cmsky[msky + lmax] = Cmm(msky + lmax, beamIndex, lat);
                cfft p1(binElements);
                p1.backward(Cmsky);
                for (long msky = 0; msky < binElements; msky++)
                    tod1(msky, beamIndex, lat) = Cmsky[msky];
            }
        } // end of parallel for
#pragma omp for schedule(static, 8)
        for (long msky = -lmax; msky <= lmax; ++msky) {
            for (long beamIndex = 0; beamIndex <= beammmax; ++beamIndex) {
                for (long lat = 0; lat < NThetaIndex1; ++lat) {
                    double tmpR = tod1(msky + lmax, beamIndex, lat).real();
                    double tmpI = tod1(msky+lmax, beamIndex, lat).imag();
                    tod1(msky + lmax, beamIndex, lat).real()
                        = (cs0[lmax + msky] * tmpR + sn0[lmax + msky] * tmpI);
                    tod1(msky + lmax, beamIndex, lat).imag()
                        = (cs0[lmax + msky] * tmpI - sn0[lmax + msky] * tmpR);
                    //Rotate in psi space
                    tmpR = tod1(msky + lmax, beamIndex, lat).real();
                    tmpI = tod1(msky + lmax, beamIndex, lat).imag();
                    tod1(msky + lmax, beamIndex, lat).real()
                        = (cs[beamIndex + lmax] * tmpR + sn[beamIndex + lmax] * tmpI);
                    tod1(msky + lmax, beamIndex, lat).imag()
                        = (cs[beamIndex + lmax] * tmpI - sn[beamIndex + lmax] * tmpR);
                }
            }
        } // end of parallel for
    } // end of parallel region
    // Note that now, -lmax <= msky <= lmax corresponds to the angle
    // 0 <= phi <= 2.*pi/2./(2.*lmax+1.) and -pi <= phi <= -2.pi/(2.*lmax+1.)
    // Finished with convolution and FFT over msky.
    t_todAnnulus_v3 = mpiMgr.Wtime() - tstart;
}


void convolver::interpolTOD_arrTestcm_v4(levels::arr<double> &outpntarr1,
                                         levels::arr<double> &outpntarr2,
                                         levels::arr<double> &TODValue,
                                         levels::arr<double> &TODValue2,
                                         long ntod, long ntod2) {
    double tstart = mpiMgr.Wtime();
    if (CMULT_VERBOSITY > 1)
        std::cout << "Entered interpolTOD_arrTestcm_v4" << std::endl;
    double phi0 = halfpi;
    //long npsi = 2*beammmax+1;
    long nphi = 2 * lmax + 1;
    double dphi = 2 * pi / nphi;
    double inv_delta_phi = 1. / dphi;
    double phioffset = phi0 / dphi;
    long ntheta = lmax + 1 + 21;
    double dtheta = -pi / (ntheta - 21);
    double theta0 = pi - 10 * dtheta;
    double inv_delta_theta = 1. / dtheta;
    long max_order = 19;
    int npoints = order + 1;
    int ioffset = order / 2;
    levels::arr<long> itheta0, itheta0_2;
    long NThetaIndex1(0), NThetaIndex2(0);
    levels::arr<long> lowerIndex, lowerIndex2;
    levels::arr<long> upperIndex, upperIndex2;
    levels::arr3<xcomplex<double> > TODAsym, TODAsym2;

    if (ntod != 0)
        itheta0SetUp(npoints, ioffset, ntheta, theta0, inv_delta_theta, nphi,
                     outpntarr1, ntod, NThetaIndex1, itheta0, lowerIndex,
                     upperIndex, TODAsym);

    if (ntod2 != 0)
        itheta0SetUp(npoints, ioffset, ntheta, theta0, inv_delta_theta, nphi,
                     outpntarr2, ntod2, NThetaIndex2, itheta0_2, lowerIndex2,
                     upperIndex2, TODAsym2);

    if (CMULT_VERBOSITY > 1)
        std::cout << "DONE with NThetaIndex1 = " << NThetaIndex1
                  << "  and NThetaIndex2 = " << NThetaIndex2 << "  corenum = "
                  << corenum << "  ntod = " << ntod << "  ntod2 = " << ntod2 << std::endl;

    bool fNThetaIndex=true;
    if (ntod != 0 && ntod2 != 0) {
        if (itheta0[ntod-1] != itheta0_2[ntod2-1]) fNThetaIndex = false;
        if (NThetaIndex1 != NThetaIndex2) fNThetaIndex = false;
    }
    if (ntod == 0 || ntod2 == 0) fNThetaIndex = false;

    levels::arr<double> rthetas1, rthetas2;
    if (NThetaIndex1 != 0) {
        rthetas1.alloc(NThetaIndex1);
        for (int thetaIndex = 0; thetaIndex < NThetaIndex1; ++thetaIndex)
            rthetas1[thetaIndex] = theta0 + (itheta0[ntod-1] + thetaIndex) * dtheta;
    }
    if (NThetaIndex2 != 0) {
        rthetas2.alloc(NThetaIndex2);
        for (int thetaIndex = 0; thetaIndex < NThetaIndex2; ++thetaIndex)
            rthetas2[thetaIndex] = pi - (theta0 + (itheta0_2[ntod2-1] +
                                                   thetaIndex) * dtheta);
    }

    double elapsed_secs = 0;
    clock_t begin = clock();

    if (fNThetaIndex)
        conviqt_hemiscm_v4(TODAsym, TODAsym2, NThetaIndex1, rthetas1);
    else {
        if (NThetaIndex1 != 0)
            conviqt_hemiscm_single(TODAsym, NThetaIndex1, rthetas1);
        if (NThetaIndex2 != 0)
            conviqt_hemiscm_single(TODAsym2, NThetaIndex2, rthetas2);
    }

    if (CMULT_VERBOSITY > 1) std::cout << "After conviqt if" << std::endl;

    clock_t end = clock();
    elapsed_secs += double(end - begin) / CLOCKS_PER_SEC;

    if (ntod > 0)
        for (int thetaIndex = 0; thetaIndex < NThetaIndex1; ++thetaIndex)
            conviqt_tod_loop_v4(lowerIndex, upperIndex, outpntarr1, TODAsym,
                                thetaIndex, itheta0, max_order, inv_delta_theta, theta0,
                                inv_delta_phi, phioffset, nphi, ioffset, npoints, ntod,
                                TODValue, thetaIndex);
    if (ntod2 > 0)
        for (int thetaIndex = 0; thetaIndex < NThetaIndex2; ++thetaIndex)
            conviqt_tod_loop_v4(lowerIndex2, upperIndex2, outpntarr2, TODAsym2,
                                thetaIndex, itheta0_2, max_order, inv_delta_theta, theta0,
                                inv_delta_phi, phioffset, nphi, ioffset, npoints, ntod2,
                                TODValue2, thetaIndex);

  if (CMULT_VERBOSITY > 2) {
      double maxtod = 0, maxtod2 = 0;
      if (ntod > 0)
          for (long ii = 0; ii < ntod; ++ii)
              if (maxtod < abs(TODValue[ii])) maxtod = abs(TODValue[ii]);
      if (ntod2 > 0)
          for (long ii = 0; ii < ntod2; ++ii)
              if (maxtod2 < abs(TODValue2[ii])) maxtod2 = abs(TODValue2[ii]);

      std::cout << "corenum = " << corenum << "  elapsed_secs for conviqt = "
                << elapsed_secs << "   average = " << elapsed_secs
                << "  NThetaIndex1 = " << NThetaIndex1 << "  NThetaIndex2 = "
                << NThetaIndex2 << "  maxtod = " << maxtod << "  maxtod2 = "
                << maxtod2 << std::endl;
  }

  if (CMULT_VERBOSITY > 1) std::cout << "Leaving interpolTOD_arrTestcm_v4" << std::endl;

  t_interpolTOD_arrTestcm_v4 =  mpiMgr.Wtime() - tstart;
}


void convolver::interpolTOD_arrTestcm_pol_v4(levels::arr<double> &outpntarr1,
                                             levels::arr<double> &outpntarr2,
                                             levels::arr<double> &TODValue,
                                             levels::arr<double> &TODValue2,
                                             long ntod, long ntod2) {
    double tstart = mpiMgr.Wtime();
    if (CMULT_VERBOSITY > 1)
        std::cout << "Entered interpolTOD_arrTestcm_pol_v4" << std::endl;
    double phi0 = halfpi;
    //long npsi = 2*beammmax + 1;
    long nphi = 2 * lmax + 1;
    double dphi = 2 * pi / nphi;
    double inv_delta_phi = 1. / dphi;
    double phioffset = phi0 / dphi;
    long ntheta = lmax + 1 + 21;
    double dtheta = -pi / (ntheta - 21);
    double theta0 = pi - 10 * dtheta;
    double inv_delta_theta = 1. / dtheta;
    long max_order = 19;
    int npoints = order + 1;
    int ioffset = order / 2;
    levels::arr<long> itheta0, itheta0_2;
    long NThetaIndex1(0), NThetaIndex2(0);
    levels::arr<long> lowerIndex, lowerIndex2;
    levels::arr<long> upperIndex, upperIndex2;
    levels::arr3<xcomplex<double> > TODAsym, TODAsym2;

    if (ntod != 0)
        itheta0SetUp(npoints, ioffset, ntheta, theta0, inv_delta_theta, nphi,
                     outpntarr1, ntod, NThetaIndex1, itheta0, lowerIndex,
                     upperIndex, TODAsym);

    if (ntod2 != 0)
        itheta0SetUp(npoints, ioffset, ntheta, theta0, inv_delta_theta, nphi,
                     outpntarr2, ntod2, NThetaIndex2, itheta0_2, lowerIndex2,
                     upperIndex2, TODAsym2);

    if (CMULT_VERBOSITY > 1) {
        std::cout << "DONE with NThetaIndex1 = " << NThetaIndex1
                  << "  and NThetaIndex2 = " << NThetaIndex2
                  << "  corenum = " << corenum << "  ntod = " << ntod
                  << "  ntod2 = " << ntod2 << std::endl;
    }

    bool fNThetaIndex = true;
    if (ntod != 0 && ntod2 != 0) {
        if (itheta0[ntod - 1] != itheta0_2[ntod2 - 1]) fNThetaIndex = false;
        if (NThetaIndex1 != NThetaIndex2) fNThetaIndex = false;
    }

    if (ntod == 0 || ntod2 == 0) fNThetaIndex = false;

    levels::arr<double> rthetas1, rthetas2;
    if (NThetaIndex1 != 0) {
        rthetas1.alloc(NThetaIndex1);
        for (int thetaIndex = 0; thetaIndex < NThetaIndex1; ++thetaIndex)
            rthetas1[thetaIndex]
                = theta0 + (itheta0[ntod - 1] + thetaIndex) * dtheta;
    }
    if (NThetaIndex2 != 0) {
        rthetas2.alloc(NThetaIndex2);
        for (int thetaIndex = 0; thetaIndex < NThetaIndex2; ++thetaIndex)
            rthetas2[thetaIndex]
                = pi - (theta0 + (itheta0_2[ntod2 - 1] + thetaIndex) * dtheta);
    }

    double elapsed_secs = 0;
    clock_t begin = clock();

    if (fNThetaIndex)
        conviqt_hemiscm_pol_v4(TODAsym, TODAsym2, NThetaIndex1, rthetas1);
    else {
        if (NThetaIndex1 != 0)
            conviqt_hemiscm_pol_single(TODAsym, NThetaIndex1, rthetas1);
        if (NThetaIndex2 != 0)
            conviqt_hemiscm_pol_single(TODAsym2, NThetaIndex2, rthetas2);
    }

    if (CMULT_VERBOSITY > 1) std::cout << "After conviqt if" << std::endl;

    clock_t end = clock();
    elapsed_secs += double(end - begin) / CLOCKS_PER_SEC;

    if (ntod > 0)
        for (int thetaIndex = 0; thetaIndex < NThetaIndex1; ++thetaIndex)
            conviqt_tod_loop_pol_v5(lowerIndex, upperIndex, outpntarr1, TODAsym,
                                    thetaIndex, itheta0, max_order, inv_delta_theta,
                                    theta0, inv_delta_phi, phioffset, nphi, ioffset,
                                    npoints, ntod, TODValue, thetaIndex);

    if (ntod2 > 0)
        for (int thetaIndex = 0; thetaIndex < NThetaIndex2; ++thetaIndex)
            conviqt_tod_loop_pol_v5(lowerIndex2, upperIndex2, outpntarr2, TODAsym2,
                                    thetaIndex, itheta0_2, max_order, inv_delta_theta,
                                    theta0, inv_delta_phi, phioffset, nphi, ioffset,
                                    npoints, ntod2, TODValue2, thetaIndex);

    if (CMULT_VERBOSITY > 1) {
        double maxtod = 0, maxtod2 = 0;
        if (ntod > 0)
          for (long ii = 0; ii < ntod; ++ii)
              if (maxtod < abs(TODValue[ii]))
                  maxtod = abs(TODValue[ii]);

        if (ntod2 > 0)
            for (long ii = 0; ii < ntod2; ii++)
              if (maxtod2 < abs(TODValue2[ii]))
                  maxtod2 = abs(TODValue2[ii]);

        std::cout << "corenum = " << corenum << "  elapsed_secs for conviqt = "
                  << elapsed_secs << "   average = " << elapsed_secs
                  << "  NThetaIndex1 = " << NThetaIndex1 << "  NThetaIndex2 = "
                  << NThetaIndex2 << "  maxtod = " << maxtod << "  maxtod2 = "
                  << maxtod2 << std::endl;

        std::cout << "Leaving interpolTOD_arrTestcm_pol_v4" << std::endl;
    }

    t_interpolTOD_arrTestcm_pol_v4 = mpiMgr.Wtime() - tstart;
}


void convolver::conviqt_tod_loop_v4(levels::arr<long> &lowerIndex,
                                    levels::arr<long> &upperIndex,
                                    levels::arr<double> &outpntarr,
                                    levels::arr3<xcomplex<double> > &TODAsym,
                                    long thetaIndex,
                                    levels::arr<long> &itheta0,
                                    long max_order, double inv_delta_theta,
                                    double theta0, double inv_delta_phi,
                                    double phioffset, long nphi, long ioffset,
                                    long npoints, long ntod,
                                    levels::arr<double> &TODValue, long lat) {
    double tstart = mpiMgr.Wtime();
    levels::arr2<std::complex<double> > conviqtarr;
    try {
        conviqtarr.alloc(nphi, beammmax + 1);
    } catch (std::bad_alloc & e) {
        std::cerr << "conviqt_tod_loop_v4 : Out of memory allocating "
                  << nphi * (beammmax + 1) * 16. / 1024 / 1024
                  << "MB for conviqtarr" << std::endl;
        throw;
    }

    for (long ii = 0; ii < nphi; ++ii)
        for (long jj = 0; jj < beammmax + 1; ++jj)
            conviqtarr[ii][jj] = TODAsym(ii, jj, lat);

    // Call weight_ncm once unthreaded to make sure the auxiliary arrays are allocated

    std::vector<double> wgt_temp(max_order + 1, 0.);
    weight_ncm(0., wgt_temp);

#pragma omp parallel for default(shared) schedule(static, 1)
    for (int ii = lowerIndex[thetaIndex]; ii >= upperIndex[thetaIndex]; --ii) {
        std::vector<double> wgt1(max_order + 1, 0.);
        // Note that the larger the ii the smaller frac is and the smaller itheta0[ii] is.
        double frac = (outpntarr[5 * ii + 1] - theta0) * inv_delta_theta;
        frac -= itheta0[ii];
        weight_ncm(frac, wgt1);

        std::vector<double> wgt2(max_order + 1, 0.);
        frac = outpntarr[5 * ii] * inv_delta_phi - phioffset;
        frac = levels::fmodulo(frac, double(nphi));
        int iphi0 = int (frac) - ioffset;
        frac -= iphi0;
        if (iphi0 >= nphi)
            iphi0 -= nphi;
        if (iphi0 < 0)
            iphi0 += nphi;
        weight_ncm(frac, wgt2);

        double omega = outpntarr[5 * ii + 2] + halfpi;
        double sinomg = sin(omega), cosomg = cos(omega);
        std::vector<double> cosang(beammmax + 1), sinang(beammmax + 1);
        cosang[0] = 1;
        sinang[0] = 0;
        for (long psiIndex = 1; psiIndex <= beammmax; ++psiIndex) {
            cosang[psiIndex] =
                cosang[psiIndex - 1] * cosomg -
                sinang[psiIndex - 1] * sinomg;
            sinang[psiIndex] =
                sinang[psiIndex - 1] * cosomg +
                cosang[psiIndex - 1] * sinomg;
        }
        cosang[0] = 0.5;
        double weight1 = wgt1[thetaIndex - (itheta0[ii] - itheta0[ntod - 1])];
        for (int phiIndex = 0; phiIndex < npoints; ++phiIndex) {
            double weight = wgt2[phiIndex] * weight1;
            long newphiIndex = iphi0 + phiIndex;
            if (newphiIndex >= nphi) newphiIndex -= nphi;
            std::complex<double> *ca = &(conviqtarr[newphiIndex][0]);
            double *cang = &(cosang[0]);
            double *sang = &(sinang[0]);
            for (long psiIndex = 0; psiIndex <= beammmax; ++psiIndex) {
                TODValue[ii] += 2 * weight * ((*cang) * (*ca).real() -
                                              (*sang) * (*ca).imag());
                ++ca;
                ++cang;
                ++sang;
            }
        }
    }

    t_conviqt_tod_loop_v4 = mpiMgr.Wtime() - tstart;
}


void convolver::conviqt_tod_loop_pol_v5(levels::arr<long> &lowerIndex,
                                        levels::arr<long> &upperIndex,
                                        levels::arr<double> &outpntarr,
                                        levels::arr3<xcomplex<double> > &TODAsym,
                                        long thetaIndex, levels::arr<long> &itheta0,
                                        long max_order, double inv_delta_theta,
                                        double theta0, double inv_delta_phi,
                                        double phioffset, long nphi, long ioffset,
                                        long npoints, long ntod,
                                        levels::arr<double> &TODValue, long lat) {
    double tstart = mpiMgr.Wtime();
    levels::arr2<std::complex<double> > conviqtarr;
    try {
        conviqtarr.alloc(nphi,beammmax+1);
    } catch (std::bad_alloc & e) {
        std::cerr << "conviqt_tod_loop_pol_v5 : Out of memory allocating "
                  << nphi*(beammmax+1)*16./1024/1024 << "MB for conviqtarr"
                  << std::endl;
        throw;
    }

    for (long ii = 0; ii < nphi; ++ii)
        for (long jj = 0; jj < beammmax + 1; ++jj)
            conviqtarr[ii][jj] = TODAsym(ii, jj, lat);

    // Call weight_ncm once unthreaded to make sure the auxiliary arrays are allocated

    std::vector<double> wgt_temp(max_order + 1, 0.);
    weight_ncm(0., wgt_temp);

#pragma omp parallel for default(shared) schedule(static, 1)
    for (int ii = lowerIndex[thetaIndex]; ii >= upperIndex[thetaIndex]; --ii) {
        std::vector<double> wgt1(max_order + 1, 0.);
        // Note that the larger the ii the smaller frac is and the smaller itheta0[ii] is.
        double frac = (outpntarr[5 * ii + 1] - theta0) * inv_delta_theta;
        frac -= itheta0[ii];
        weight_ncm(frac, wgt1);

        std::vector<double> wgt2(max_order + 1, 0.);
        frac = outpntarr[5 * ii] * inv_delta_phi - phioffset;
        frac = levels::fmodulo(frac, double(nphi));
        int iphi0 = int(frac) - ioffset;
        frac -= iphi0;
        if (iphi0 >= nphi)
            iphi0 -= nphi;
        if (iphi0 < 0)
            iphi0 += nphi;
        weight_ncm(frac, wgt2);

        double omega = outpntarr[5 * ii + 2] + halfpi;
        double sinomg = sin(omega), cosomg = cos(omega);
        std::vector<double> cosang(beammmax + 1), sinang(beammmax + 1);
        cosang[0] = 1;
        sinang[0] = 0;
        for (long psiIndex = 1; psiIndex <= beammmax; ++psiIndex) {
            cosang[psiIndex] = cosang[psiIndex - 1] * cosomg - sinang[psiIndex - 1] * sinomg;
            sinang[psiIndex] = sinang[psiIndex - 1] * cosomg + cosang[psiIndex - 1] * sinomg;
        }
        cosang[0] = 0.5;
        double weight1 = wgt1[thetaIndex - (itheta0[ii] - itheta0[ntod - 1])];
        for (int phiIndex = 0; phiIndex < npoints; ++phiIndex) {
            double weight = wgt2[phiIndex] * weight1;
            long newphiIndex = iphi0 + phiIndex;
            if (newphiIndex >= nphi) newphiIndex -= nphi;
            std::complex<double> *ca = &(conviqtarr[newphiIndex][0]);
            double *cang = &(cosang[0]);
            double *sang = &(sinang[0]);
            for (long psiIndex = 0; psiIndex <= beammmax; ++psiIndex) {
                TODValue[ii] += 2 * weight * ((*cang) * (*ca).real() -
                                              (*sang) * (*ca).imag());
                ++ca;
                ++cang;
                ++sang;
            }
        }
    }

    t_conviqt_tod_loop_pol_v5 = mpiMgr.Wtime() - tstart;
}


void convolver::arrFillingcm_v2(long ntod, levels::arr<double> &timeTest_arr,
                                levels::arr<double> &outpntarrx,
                                levels::arr<double> &outpntarr, long offindex) {
    double tstart = mpiMgr.Wtime();
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

    t_arrFillingcm_v2 = mpiMgr.Wtime() - tstart;
}


int convolver::convolve(pointing & pntarr, bool calibrate) {
    double tstart = mpiMgr.Wtime();
    int lmax_sky = s->get_lmax();
    int lmax_beam = b->get_lmax();
    if (lmax < 0) {
        lmax = lmax_sky < lmax_beam ? lmax_sky : lmax_beam;
    } else if (lmax > lmax_sky || lmax > lmax_beam) {
        throw std::runtime_error("Convolver lmax exceeds input expansion order.");
    }

    int mmax_beam = b->get_mmax();
    if (beammmax < 0) {
        beammmax = mmax_beam;
    } else if (beammmax > mmax_beam) {
        throw std::runtime_error("Convolver mmax exceeds input expansion order.");
    }

    long totsize = pntarr.size() / 5;

    // Assign a running index across the communicator to the last column of pntarr
    // It is needed to collect the convolved data.
    long my_offset = 0;
    int err = MPI_Scan(&totsize, &my_offset, 1, MPI_LONG, MPI_SUM, mpiMgr.comm());
    if (err != 0) throw std::runtime_error("Error scannign total TOD size");
    my_offset -= totsize;
    for (long ii = 0; ii < totsize; ++ii) pntarr[5 * ii + 4] = my_offset + ii;

    long ntod = 0, ntod2 = 0, outBetaSegSize;
    levels::arr<int> inBetaSeg, inBetaSegFir, inBetaSegSec, outBetaSeg;
    inBetaSeg.alloc(cores);
    inBetaSeg.fill(0);
    inBetaSegFir.alloc(cores);
    inBetaSegFir.fill(0);
    inBetaSegSec.alloc(cores);
    inBetaSegSec.fill(0);
    levels::arr<int> inBetaSegAcc, outBetaSegAcc;
    levels::arr<double> pntarr2, outpntarr;
    levels::arr<int> inOffset, outOffset;
    // Make sure theta=pi doesn't break anything
    double ratiodeltas = (halfpi + 1e-10) / cores;
    double dtheta, newtheta, thetaini = halfpi;
    double countdeltatheta = 0;
    levels::arr<double> corethetaarr(cores), coredthetaarr(cores);
    levels::arr<double> dbeta;

    // Check the angles for conformity

    for (long i = 0; i < totsize; ++i) {
        double theta = pntarr[5 * i + 1];
        if (theta < 0 || theta > pi)
            throw std::runtime_error("convolver::convolve: Illegal theta: "
                                     + std::to_string(theta) + " not in 0..pi");
    }

    // distribute the latitudes

    for (long ii = cores - 1; ii >= 0; --ii) {
        if (thetaini > 0.18) {
            thetaDeltaThetacm(thetaini, ii, newtheta, dtheta);
        } else {
            long iival = ii;
            deltaTheta2(iival, thetaini, dbeta);
            for (long jj = ii; jj >= 0; --jj) {
                newtheta -= dbeta[jj];
                corethetaarr[jj] = newtheta;
                coredthetaarr[jj] = dbeta[jj];
                thetaini = newtheta;
            }
            newtheta = 0;
            corethetaarr[0] = 0;
            break;
        }
        corethetaarr[ii] = newtheta;
        coredthetaarr[ii] = dtheta;
        thetaini = newtheta;
        countdeltatheta += sin(newtheta + dtheta / 2.) * dtheta;
    }
    fillingBetaSeg(pntarr, totsize, ratiodeltas, corethetaarr, inBetaSeg);
    mpiMgr.all2all(inBetaSeg, outBetaSeg);

    todRedistribution5cm(pntarr, inBetaSeg, outBetaSeg, inBetaSegAcc, outBetaSegAcc,
                         outBetaSegSize, pntarr2, totsize, inOffset, outOffset,
                         ratiodeltas, corethetaarr);

    // Store input TOD for comparison

    levels::arr<double> todtmp;
    if (totsize > 0) {
        todtmp.alloc(totsize);
        for (long ii = 0; ii < totsize; ii++) todtmp[ii] = pntarr[5 * ii + 3];
        pntarr.dealloc();
    }

    // Redistribute the data from pntarr2 to outpntarr

    if (totsize != 0 || outBetaSegSize != 0)
        mpiMgr.all2allv(pntarr2, inBetaSeg, inOffset, outpntarr, outBetaSeg, outOffset,
                        outBetaSegSize);

    if (totsize > 0) pntarr2.dealloc();

    // Count elements on Northern (ntod) and Southern (ntod2) hemisphere

    long nthetaCount(0);
    long nthetaCount2(0);
    long csmall(0), cbig(0);
    for (long ii = 0; ii < outBetaSegSize / 5; ++ii) {
        if (outpntarr[5 * ii + 1] < halfpi) {
            nthetaCount++;
            if (outpntarr[5 * ii + 1] <= 0)
                csmall++;
        } else {
            nthetaCount2++;
            if (outpntarr[5 * ii + 1] >= pi)
                cbig++;
        }
        // The signal column in outpntarr is replaced with a local row number to allow reordering
        outpntarr[5 * ii + 3] = ii;
    }
    ntod = nthetaCount;
    ntod2 = nthetaCount2;

    if (outBetaSegSize != 0)
        hpsort_arrTheta(outpntarr); // Sort according to latitude

    // Run the convolution

    clock_t begin = clock();

    levels::arr<double> todTest_arr, todTest_arr2, timeTest_arr, timeTest_arr2;

    todgen_v4(ntod, ntod2, todTest_arr,
              timeTest_arr, todTest_arr2, timeTest_arr2, outpntarr);

    if (CMULT_VERBOSITY > 1) {
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        levels::arr<double> tarr(cores, elapsed_secs), toutarr;
        levels::arr<long> arrcores(cores);
        for (long ii = 0; ii < cores; ++ii)
            arrcores[ii] = ii;
        levels::arr<long> outbetasegsizearr(cores, outBetaSegSize), obt;
        mpiMgr.all2all(tarr, toutarr);
        mpiMgr.all2all(outbetasegsizearr, obt);
        double maxtint = 0;
        double mintint = 1.e30;
        long corevaluemax, corevaluemin;
        for (long ii = 0; ii < cores; ++ii) {
            if (toutarr[ii] > maxtint) {
                maxtint = toutarr[ii];
                corevaluemax = ii;
            }
            if (toutarr[ii] < mintint) {
                mintint = toutarr[ii];
                corevaluemin = ii;
            }
        }
        std::cout << "corenum = " << corenum
                  << "   elapsed_secs = " << elapsed_secs << std::endl;
        if (corenum == 0)
            std::cout << "maxtint = " << maxtint << "   mintint = " << mintint
                      << "   corevaluemax = " << corevaluemax
                      << "   corevaluemin = " << corevaluemin << std::endl;

        // Sort task IDs in arrcores according to execution time in toutarr

        hpsort_DL(toutarr, arrcores);

        if (corenum == 0)
            for (long ii = 0; ii < cores; ii++)
                std::cout << "toutarr = " << toutarr[ii] << "   corenum = "
                          << arrcores[ii] << "   outBetaSegSize = "
                          << obt[arrcores[ii]] << std::endl;
    }

    // Sort outpntarr according to the 4th column
    if (outBetaSegSize != 0)
        hpsort_arrTOD(outpntarr);
    if (totsize != 0 || outBetaSegSize != 0)
        mpiMgr.all2allv(outpntarr, outBetaSeg, outOffset, pntarr,
                        inBetaSeg, inOffset, 5 * totsize);

    if (outBetaSegSize != 0) outpntarr.dealloc();
    // Sort according to the 5th column, time stamp
    if (totsize != 0) hpsort_arrTime(pntarr);

    levels::arr<double> todAll;
    preReorderingStep(ntod, ntod2, todAll, todTest_arr, todTest_arr2);
    levels::arr<double> timeAll;
    preReorderingStep(ntod, ntod2, timeAll, timeTest_arr, timeTest_arr2);

    // sort time and signal according to time
    if (outBetaSegSize != 0) hpsort_DDcm(timeAll, todAll);

    for (long ii = 0; ii < cores; ii++) {
        inBetaSeg[ii] = inBetaSeg[ii] / 5;
        outBetaSeg[ii] = outBetaSeg[ii] / 5;
        inOffset[ii] = inOffset[ii] / 5;
        outOffset[ii] = outOffset[ii] / 5;
    }

    // collect the convolved TOD

    levels::arr<double> outtodarr, outnumarr;
    if (totsize != 0 || outBetaSegSize != 0)
        mpiMgr.all2allv(todAll, outBetaSeg, outOffset, outtodarr, inBetaSeg, inOffset, totsize);
    if (totsize != 0 || outBetaSegSize != 0)
        mpiMgr.all2allv(timeAll, outBetaSeg, outOffset, outnumarr, inBetaSeg, inOffset, totsize);

    if (outBetaSegSize != 0) todAll.dealloc();
    if (outBetaSegSize != 0) timeAll.dealloc();

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
    if (calibrate)
        calibration = 2. / (1. + d->get_epsilon());

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

    t_convolve = mpiMgr.Wtime() - tstart;
    //if (CMULT_VERBOSITY > 1)
    report_timing();

    return 0;
}

void convolver::report_timing() {
    if (corenum == 0)
        std::cerr << "Conviqt::convolve timing:" << std::endl;
    timing_line(std::string("convolve"), t_convolve);
    timing_line(std::string("  todgen_v4"), t_todgen_v4);
    timing_line(std::string("    arrFillingcm_v2"), t_arrFillingcm_v2);
    timing_line(std::string("    interpolTOD_arrTestcm_v4"),
                t_interpolTOD_arrTestcm_v4);
    timing_line(std::string("    interpolTOD_arrTestcm_pol_v4"),
                t_interpolTOD_arrTestcm_pol_v4);
    timing_line(std::string("        itheta0SetUp"), t_itheta0SetUp);
    timing_line(std::string("            ithetacalc"), t_ithetacalc);
    timing_line(std::string("        conviqt_hemiscm_single"),
                t_conviqt_hemiscm_single);
    timing_line(std::string("        conviqt_hemiscm_v4"),
                t_conviqt_hemiscm_v4);
    timing_line(std::string("        conviqt_hemiscm_pol_v4"),
                t_conviqt_hemiscm_pol_v4);
    timing_line(std::string("        conviqt_hemiscm_pol_single"),
                t_conviqt_hemiscm_pol_single);
    timing_line(std::string("            wigner_init"), t_wigner_init);
    timing_line(std::string("            wigner_prepare"), t_wigner_prepare);
    timing_line(std::string("            wigner_calc"), t_wigner_calc);
    timing_line(std::string("            lat_iter"), t_lat_iter);
    timing_line(std::string("            sincos_iter"), t_sincos_iter);
    timing_line(std::string("            todAnnulus_v3"), t_todAnnulus_v3);
    timing_line(std::string("        conviqt_tod_loop_v4"), t_conviqt_tod_loop_v4);
    timing_line(std::string("        conviqt_tod_loop_pol_v5"),
                t_conviqt_tod_loop_pol_v5);
    timing_line(std::string("            weight_ncm"), t_weight_ncm);
}

void convolver::timing_line(std::string label, double timer) {
    double tmin, tmax, tmean, tstd;
    mpiMgr.reduceRaw<double>(&timer, &tmin, 1, mpiMgr.Min);
    mpiMgr.reduceRaw<double>(&timer, &tmax, 1, mpiMgr.Max);
    mpiMgr.reduceRaw<double>(&timer, &tmean, 1, mpiMgr.Sum);
    tmean /= cores;
    double tsquare = ((timer - tmean) * (timer - tmean));
    mpiMgr.reduceRaw<double>(&tsquare, &tstd, 1, mpiMgr.Sum);
    tstd = sqrt(tstd / cores);
    if (corenum == 0 && tmax > 1e-1) {
        fprintf(stderr, "%-35s %8.2f < %8.2f +- %8.2f < %8.2f\n",
                label.c_str(), tmin, tmean, tstd, tmax);
        //std::cerr << label << " : " << tmin
        //          << " < " << tmean
        //          << " +- " << tstd
        //          << " < " << tmax << std::endl;
    }
}


} // namespace conviqt
