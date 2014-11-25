/*
 *  This file is part of libc_utils.
 *
 *  libc_utils is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libc_utils is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libc_utils; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libc_utils is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*! \file wrec3a.c
 *  Computation of Wigner-3j symbols
 *  Algorithm implemented according to Schulten & Gordon:
 *  J. Math. Phys. 16, p. 10 (1975)
 *
 *  Copyright (C) 2009 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "c_utils.h"
#include "wrec3a.h"

static inline double max (double a, double b)
  { return (a>=b) ? a : b; }

static inline int nearest_int (double arg)
  { return (int)(floor(arg+0.5)); }

static inline int abs_approx (double a, double b, double epsilon)
  { return fabs(a-b) < epsilon; }

static inline int intcheck (double val)
  { return abs_approx(val,(double)nearest_int(val),1e-13); }

static inline double xpow (int m, double val)
  { return (m&1) ? -val : val; }

void wrec3jj (double l2, double l3, double m2, double m3, double *res, int sz)
  {
  const int expo=250;
  const double srhuge=ldexp(1.,expo),
               tiny=ldexp(1.,-2*expo), srtiny=ldexp(1.,-expo);

  const double m1 = -m2 -m3;
  const double l1min = max(fabs(l2-l3),fabs(m1)),
               l1max = l2 + l3;
  const int ncoef = nearest_int(l1max-l1min)+1;

  const double l2ml3sq = (l2-l3)*(l2-l3),
               pre1 = (l2+l3+1.)*(l2+l3+1.),
               m1sq = m1*m1,
               pre2 = m1*(l2*(l2+1.)-l3*(l3+1.)),
               m3mm2 = m3-m2;

  int i=0;
  double sumfor,sumbac;
  double c1=1e300;
  double oldfac=0.;

  UTIL_ASSERT (l2>=fabs(m2),"l2<fabs(m2)");
  UTIL_ASSERT (l3>=fabs(m3),"l3<fabs(m3)");
  UTIL_ASSERT (intcheck(l2+fabs(m2)),"l2+fabs(m2) is not integer");
  UTIL_ASSERT (intcheck(l3+fabs(m3)),"l3+fabs(m3) is not integer");

  UTIL_ASSERT (intcheck(l1max-l1min), "l1max-l1min is not integer");
  UTIL_ASSERT (l1max>=l1min, "l1max is smaller than l1min");

  UTIL_ASSERT(ncoef<=sz,"array too small");

  res[i] = srtiny;
  sumfor = (2.*l1min+1.) * res[i]*res[i];

  while(1)
    {
    if (i==ncoef-1) break; /* all done */
    ++i;
    {
    const double l1 = l1min+i,
                 l1sq = l1*l1,
                 c1old=fabs(c1),
                 newfac = sqrt((l1sq-l2ml3sq)*(pre1-l1sq)*(l1sq-m1sq));

    if (i>1)
      {
      const double tmp1 = 1./((l1-1.)*newfac);
      c1 = (2.*l1-1.)*(pre2-(l1sq-l1)*m3mm2) * tmp1;
      res[i] = res[i-1]*c1 - res[i-2]*l1*oldfac*tmp1;
      }
    else
      {
      c1 = (l1>1.000001) ? (2.*l1-1.)*(pre2-(l1sq-l1)*m3mm2)/((l1-1.)*newfac)
                         : -(2.*l1-1.)*l1*(m3mm2)/newfac;
      res[i] = res[i-1]*c1;
      }

    oldfac=newfac;

    sumfor += (2.*l1+1.)*res[i]*res[i];
    if (fabs(res[i])>=srhuge)
      {
      int k;
      for (k=0; k<=i; ++k)
        res[k]*=srtiny;
      sumfor*=tiny;
      }
    if (c1old<=fabs(c1)) break;
    }
    }

  sumbac=0.;
  if (i<ncoef-1) /* we have to iterate from the other side */
    {
    const double x1=res[i-2], x2=res[i-1], x3=res[i];
    const int nstep2 = i-2;

    i=ncoef-1;
    res[i] = srtiny;
    sumbac = (2.*l1max+1.) * res[i]*res[i];

    do
      {
      --i;
      {
      const double l1 = l1min+i,
                   l1p1sq = (l1+1.)*(l1+1.),
                   newfac = sqrt((l1p1sq-l2ml3sq)*(pre1-l1p1sq)*(l1p1sq-m1sq));

      if (i<ncoef-2)
        res[i] = (res[i+1] * (2.*l1+3.)*(pre2-(l1p1sq+l1+1.)*m3mm2)
                 -res[i+2] * (l1+1.)*oldfac)
                 / ((l1+2.)*newfac);
      else
        res[i] = res[i+1]*(2.*l1+3.)*(pre2-(l1p1sq+l1+1.)*m3mm2)
                 /((l1+2.)*newfac);

      oldfac=newfac;

      if (i>nstep2+2)
        sumbac += (2.*l1+1.)*res[i]*res[i];
      if (fabs(res[i])>=srhuge)
        {
        int k;
        for (k=i; k<ncoef; ++k)
          res[k]*=srtiny;
        sumbac*=tiny;
        }
      }
      }
    while (i>nstep2);

    {
    const double ratio = (x1*res[i]+x2*res[i+1]+x3*res[i+2])
                         /(x1*x1+x2*x2+x3*x3);
    int k;
    for (k=0; k<nstep2; ++k)
      res[k]*=ratio;
    sumfor*=ratio*ratio;
    }
    }

  {
  double cnorm=1./sqrt(sumfor+sumbac);
  int k;
  if (xpow(nearest_int(l2-l3-m1),res[ncoef-1])<0.)
    cnorm = -cnorm;

  for (k=0; k<ncoef; ++k)
    res[k]*=cnorm;
  }
  }
