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

/*
 *  Functions for singular value decomposition
 *
 *  Copyright (C) 2008, 2009 Max-Planck-Society
 *  Author: Martin Reinecke
 */

/* Port of the Golub & Reinsch algorithm for singular value decomposition.
   This code has not been extensively tested yet, so use at your own risk! */

#include <stdio.h>
#include <math.h>
#include "c_utils.h"
#include "svd_c.h"

static double pythag (double a, double b)
  {
  double at=fabs(a), bt=fabs(b), ct;
  if (at>bt)
    {
    ct = bt/at;
    return at*sqrt(1.+ct*ct);
    }
  else
    {
    if (bt==0) return 0.;
    ct = at/bt;
    return bt*sqrt(1.+ct*ct);
    }
  }

void svd_destroy (svd_obj *obj)
  {
  DEALLOC2D(obj->u);
  DEALLOC2D(obj->v);
  DEALLOC(obj->w);
  }

void svd_init (double **matrix, double epsilon, int m, int n,
  svd_obj *obj)
  {
  int flag, l = 0, nm = 0, mn;
  double c, f, h, s, x, y, z;
  double *stmp, *rv1=RALLOC(double,n);
  double **u, **v, *w;
  int i, j, k, its, jj;
  double anorm, g, scale, wmin, wmax;

  ALLOC2D(u,double,n,m);
  ALLOC2D(v,double,n,n);
  ALLOC(w,double,n);

  obj->m=m;
  obj->n=n;
  obj->u=u;
  obj->v=v;
  obj->w=w;

  for (i=0; i<n; ++i)
    for (j=0; j<m; ++j)
      u[i][j]=matrix[j][i];

  /* Householder reduction to bidiagonal form. */
  anorm = g = scale = 0.0;
  stmp=RALLOC(double,m);
  for (i=0; i<n; ++i)
    {
    l = i + 1;
    rv1[i] = scale*g;
    g = scale = 0.0;
    if (i<m)
      {
      for (k=i; k<m; ++k) scale += fabs(u[i][k]);
      if (scale)
        {
        s = 0;
        for (k=i; k<m; ++k)
          {
          u[i][k] /= scale;
          s += u[i][k]*u[i][k];
          }
        f = u[i][i];
        g = (f>=0.) ? -sqrt(s) : sqrt(s);
        h = f*g - s;
        u[i][i] = f - g;
        if (i!=(n-1))
          {
          for (j=l; j<n; ++j)
            {
            s = 0.;
            for (k=i; k<m; ++k) s += u[i][k]*u[j][k];
            f = s/h;
            for (k=i; k<m; ++k) u[j][k] += f*u[i][k];
            }
          }
        for (k=i; k<m; ++k) u[i][k] *= scale;
        }
      }
    w[i] = scale*g;
    g = scale = 0.0;
    if ((i<m) && (i!=(n-1)))
      {
      for (k=l; k<n; ++k)
        scale += fabs(u[k][i]);
      if (scale)
        {
        s=0;
        for (k=l; k<n; ++k)
          {
          u[k][i] /= scale;
          s += u[k][i]*u[k][i];
          }
        f = u[l][i];
        g = (f>=0.) ? -sqrt(s) : sqrt(s);
        h = f*g - s;
        u[l][i] = f - g;
        for (k=l; k<n; ++k) rv1[k] = u[k][i]/h;
        if (i!=(m-1))
          {
          for (j=l; j<m; ++j)
            stmp[j]=0;
          for (k=l; k<n; ++k)
            for (j=l; j<m; ++j)
              stmp[j]+=u[k][j]*u[k][i];
          for (k=l; k<n; ++k)
            for (j=l; j<m; ++j)
              u[k][j]+=stmp[j]*rv1[k];
          }
        for (k=l; k<n; ++k) u[k][i] *= scale;
        }
      }
    anorm = IMAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
    }
  DEALLOC(stmp);

  /* Accumulation of right-hand transformations. */
  stmp=RALLOC(double,n);
  for (i=n-1; i>=0; --i)
    {
    if (i<(n-1))
      {
      if (g)
        {
        for (j=l; j<n; ++j)
          {
          /* Double division to avoid possible underflow: */
          v[i][j] = (u[j][i]/u[l][i])/g;
          stmp[j] = u[j][i];
          }
        for (j=l; j<n; ++j)
          {
          s=0;
          for (k=l; k<n; ++k) s += stmp[k]*v[j][k];
          for (k=l; k<n; ++k) v[j][k] += s*v[i][k];
          }
        }
      for (j=l; j<n; ++j) v[i][j] = v[j][i] = 0.0;
      }
    v[i][i] = 1.0;
    g = rv1[i];
    l = i;
    }
  DEALLOC(stmp);

  /* Accumulation of left-hand transformations. */
  mn=n;
  if (m<n) mn=m;
  for (i=mn-1; 0<=i; --i)
    {
    l = i + 1;
    g = w[i];
    if (i<(n-1))
      for (j=l; j<n; ++j) u[j][i] = 0.0;
    if (g)
      {
      g = 1.0/g;
      if (i!=(n-1))
        {
        for (j=l; j<n; ++j)
          {
          s=0;
          for (k=l; k<m; ++k) s += u[i][k]*u[j][k];
          f = (s/u[i][i])*g;
          for (k=i; k<m; ++k) u[j][k] += f*u[i][k];
          }
        }
      for (j=i; j<m; ++j) u[i][j] *= g;
      }
    else
      for (j=i; j<m; ++j) u[i][j] = 0.0;
    ++u[i][i];
    }

  /* Diagonalization of the bidiagonal form. */
  for (k=n-1; 0<=k; --k) /* Loop over singular values. */
    {
    for (its=0; its<30; ++its) /* Loop over allowed iterations. */
      {
      flag = 1;
      for (l=k; 0<=l; --l) /* Test for splitting: */
        {
        nm = l - 1; /* Note that rv1[0] is always zero. */
        if (fabs(rv1[l]) + anorm == anorm)
          {
          flag = 0;
          break;
          }
        if (fabs(w[nm]) + anorm == anorm)
          break;
        }
      if (flag)
        {
        c = 0.0; /* Cancellation of rv1[l], if l>0: */
        s = 1.0;
        for (i=l; i<=k; ++i)
          {
          f = s*rv1[i];
          if (fabs(f) + anorm != anorm)
            {
            g = w[i];
            h = pythag(f, g);
            w[i] = h;
            h = 1.0/h;
            c = g*h;
            s = (-f*h);
            for (j=0; j<m; ++j)
              {
              y = u[nm][j];
              z = u[i][j];
              u[nm][j] = y*c + z*s;
              u[i][j]  = z*c - y*s;
              }
            }
          }
        }
      z = w[k];
      if (l==k) /* Convergence. */
        {
        if (z<0.0) /* Singular value is made non-negative. */
          {
          w[k] = -z;
          for (j=0; j<n; ++j)
            v[k][j] = -v[k][j];
          }
        break;
        };
      if (its>=29)
        {
        fprintf(stderr,"No convergence in 30 svdcmp iterations.");
        exit(1);
        }
      x = w[l]; /* Shift from bottom 2-by-2 minor. */
      nm = k - 1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y - z)*(y + z) + (g - h)*(g + h))/(2.0*h*y);
      g = (f>=0) ? pythag(f, 1.0) : -pythag(f,1.0);
      f = ((x - z)*(x + z) + h*((y/(f + g)) - h))/x;
      /* Next QR transformation: */
      c = s = 1.0;
      for (j=l; j<=nm; ++j)
        {
        i = j + 1;
        g = rv1[i];
        y = w[i];
        h = s*g;
        g = c*g;
        z = pythag(f, h);
        rv1[j] = z;
        c = f/z;
        s = h/z;
        f = x*c + g*s;
        g = g*c - x*s;
        h = y*s;
        y = y*c;
        for (jj=0; jj<n; ++jj)
          {
          x = v[j][jj];
          z = v[i][jj];
          v[j][jj] = x*c + z*s;
          v[i][jj] = z*c - x*s;
          }
        z = pythag(f, h);
        w[j] = z; /* Rotation can be arbitrary if z = 0. */
        if (z)
          {
          z = 1.0/z;
          c = f*z;
          s = h*z;
          }
        f = (c*g) + (s*y);
        x = (c*y) - (s*g);
        for (jj=0; jj<m; ++jj)
          {
          y = u[j][jj];
          z = u[i][jj];
          u[j][jj] = y*c + z*s;
          u[i][jj] = z*c - y*s;
          }
        }
      rv1[l] = 0.0;
      rv1[k] = f;
      w[k] = x;
      }
    }

  DEALLOC(rv1);

  wmax = 0.0;
  for (j=0; j<n; ++j)
    if (w[j]>wmax) wmax=w[j];
  wmin = wmax*epsilon;
  for (j=0; j<n; ++j)
    if (w[j]<wmin) w[j]=0.0;
  }

void svd_solve (const svd_obj *obj, double *b)
  {
  int i,j,jj,n=obj->n,m=obj->m;
  double **u=obj->u, **v=obj->v, *w=obj->w;
  double *tmp = RALLOC(double,n);

  for (j=0; j<n; ++j)
    {
    double s=0.0;
    if (w[j])
      {
      for (i=0; i<m; ++i) s += u[j][i]*b[i];
      s /= w[j];
      }
    tmp[j]=s;
    }
  for (j=0; j<n; ++j)
    {
    double s=0.0;
    for (jj=0; jj<n; ++jj) s += v[jj][j]*tmp[jj];
    b[j]=s;
    }
  DEALLOC(tmp);
  }
