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

/*! \file quaternion.cc
 *  Quaternion class for rotation transforms in 3D space
 *
 *  Copyright (C) 2005-2011 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include "quaternion.h"

using namespace std;

quaternion::quaternion(const vec3 &axis, double angle)
  {
  angle*=0.5;
  double sa=sin(angle);
  w = cos(angle);
  x = sa*axis.x;
  y = sa*axis.y;
  z = sa*axis.z;
  }

quaternion::quaternion(const rotmatrix &mat)
  {
  double trace = mat.entry[0][0]+mat.entry[1][1]+mat.entry[2][2];
  if (trace>0.)
    {
    double s = 2.*sqrt(1.+trace);
    w = 0.25*s;
    s = 1./s;
    x = s * (mat.entry[2][1]-mat.entry[1][2]);
    y = s * (mat.entry[0][2]-mat.entry[2][0]);
    z = s * (mat.entry[1][0]-mat.entry[0][1]);
    }
  else if (mat.entry[0][0]>max(mat.entry[1][1],mat.entry[2][2]))
    {
    double s = 2.*sqrt(1.+mat.entry[0][0]-mat.entry[1][1]-mat.entry[2][2]);
    x = 0.25*s;
    s = 1./s;
    y = s * (mat.entry[1][0]+mat.entry[0][1]);
    z = s * (mat.entry[0][2]+mat.entry[2][0]);
    w = s * (mat.entry[2][1]-mat.entry[1][2]);
    }
  else if (mat.entry[1][1]>mat.entry[2][2])
    {
    double s = 2.*sqrt(1.+mat.entry[1][1]-mat.entry[0][0]-mat.entry[2][2]);
    y = 0.25*s;
    s = 1./s;
    x = s * (mat.entry[1][0]+mat.entry[0][1]);
    z = s * (mat.entry[2][1]+mat.entry[1][2]);
    w = s * (mat.entry[0][2]-mat.entry[2][0]);
    }
  else
    {
    double s = 2.*sqrt(1.+mat.entry[2][2]-mat.entry[0][0]-mat.entry[1][1]);
    z = 0.25*s;
    s = 1./s;
    x = s * (mat.entry[0][2]+mat.entry[2][0]);
    y = s * (mat.entry[2][1]+mat.entry[1][2]);
    w = s * (mat.entry[1][0]-mat.entry[0][1]);
    }
  }

ostream &operator<< (ostream &os, const quaternion &q)
  {
  os << "(" << q.w << "," << q.x << "," << q.y << "," << q.z << ")" << endl;
  return os;
  }
