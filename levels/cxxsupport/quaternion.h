/*
 *  This file is part of libcxxsupport.
 *
 *  libcxxsupport is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libcxxsupport is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libcxxsupport; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libcxxsupport is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*! \file quaternion.h
 *  Quaternion class for rotation transforms in 3D space
 *
 *  Copyright (C) 2005-2011 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_QUATERNION_H
#define PLANCK_QUATERNION_H

#include "vec3.h"
#include "rotmatrix.h"
#include <cmath>

/*! \defgroup quaterniongroup Quaternions */
/*! \{ */

/*! Quaternion class for rotation transforms in 3D space */
class quaternion
  {
  public:
    double w,x,y,z;

    quaternion() {}
    quaternion(double W, double X, double Y, double Z)
      : w(W), x(X), y(Y), z(Z) {}
    quaternion(const vec3 &v)
      : w(0), x(v.x), y(v.y), z(v.z) {}
    quaternion(const vec3 &axis, double angle);
    quaternion(const rotmatrix &mat);

    void Set (double W, double X, double Y, double Z)
      { w=W; x=X; y=Y; z=Z; }

    double Norm() const
      { return w*w+x*x+y*y+z*z; }

    void Flip()
      { w=-w; x=-x; y=-y; z=-z; }
    quaternion operator- () const
      { return quaternion(-w,-x,-y,-z); }

    void Conjugate()
      { x=-x; y=-y; z=-z; }
    quaternion conj() const
      { return quaternion(w,-x,-y,-z); }

    quaternion &operator*= (const quaternion &b)
      {
      quaternion a=*this;
      w = a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z;
      x = a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y;
      y = a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x;
      z = a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w;
      return *this;
      }
    quaternion operator* (const quaternion &b) const
      {
      return quaternion (w*b.w - x*b.x - y*b.y - z*b.z,
                         w*b.x + x*b.w + y*b.z - z*b.y,
                         w*b.y - x*b.z + y*b.w + z*b.x,
                         w*b.z + x*b.y - y*b.x + z*b.w);
      }

    void toAxisAngle (vec3 &axis, double &angle) const
      {
      double norm = x*x + y*y + z*z;
      if (norm==0.)
        { axis.Set(0,0,1); angle=0; return; }
      norm = sqrt(norm);
      angle = 2*atan2(norm,w);
      double inorm = 1./norm;
      axis.Set(x*inorm, y*inorm, z*inorm);
      }

    void make_rotmatrix (rotmatrix &mat) const
      {
      double s=2./(x*x+y*y+z*z+w*w);
      double xs=s*x, ys=s*y, zs=s*z;
      {
      double xx=x*xs, yy=y*ys, zz=z*zs;
      mat.entry[0][0] = 1.-yy-zz;
      mat.entry[1][1] = 1.-xx-zz;
      mat.entry[2][2] = 1.-xx-yy;
      }
      {
      double xy=x*ys, zw=zs*w;
      mat.entry[0][1] = xy-zw; mat.entry[1][0] = xy+zw;
      }
      {
      double xz=x*zs, yw=ys*w;
      mat.entry[0][2] = xz+yw; mat.entry[2][0] = xz-yw;
      }
      {
      double yz=y*zs, xw=xs*w;
      mat.entry[1][2] = yz-xw; mat.entry[2][1] = yz+xw;
      }
      }

    operator rotmatrix() const
      {
      rotmatrix result;
      make_rotmatrix(result);
      return result;
      }
  };

inline double inner_product (const quaternion &a, const quaternion &b)
  { return a.w*b.w + a.x*b.x + a.y*b.y + a.z*b.z; }

std::ostream &operator<< (std::ostream &os, const quaternion &q);

/*! \} */

#endif
