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

/*! \file io_utils.h
 *  Convenience functions for I/O with iohandles.
 *
 *  Copyright (C) 2009, 2010 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_IO_UTILS_H
#define PLANCK_IO_UTILS_H

#include "iohandle.h"
#include "pointing.h"
#include "quaternion.h"
#include "arr.h"

/*! Reads colatitudes from the column with the number \a ctheta and longitudes
    from the column with the number \a cphi from \a inp, starting at the column
    entry with the number \a offset, and returns them in \a pt. The size of
    \a pt determines the number of values read. */
void readPointing (iohandle &inp, int ctheta, int cphi,
  arr<pointing> &pt, uint64 offset);

/*! Reads colatitudes from the column with the name \a ntheta and longitudes
    from the column with the name \a nphi from \a inp, starting at the column
    entry with the number \a offset, and returns them in \a pt. The size of
    \a pt determines the number of values read. */
inline void readPointing (iohandle &inp, const std::string &ntheta,
  const std::string &nphi, arr<pointing> &pt, uint64 offset)
  {
  readPointing (inp, inp.columnNumber(ntheta), inp.columnNumber(nphi),
                pt, offset);
  }

/*! Reads three sets of pointings from 6 columns of \a inp, starting at the
    column entry with the number \a offset. */
void readPointing (iohandle &inp, int ctx, int cpx, int cty, int cpy,
  int ctz, int cpz, arr<pointing> &px, arr<pointing> &py, arr<pointing> &pz,
  uint64 offset);


/*! Reads all colatitudes from the column with the number \a ctheta and
    all longitudes from the column with the number \a cphi from \a inp, and
    returns them in \a pt (which is resized accordingly). */
void readEntirePointing (iohandle &inp, int ctheta, int cphi,
  arr<pointing> &pt);

/*! Reads all colatitudes from the column with the name \a ntheta and
    all longitudes from the column with the name \a nphi from \a inp, and
    returns them in \a pt (which is resized accordingly). */
inline void readEntirePointing (iohandle &inp, const std::string &ntheta,
  const std::string &nphi, arr<pointing> &pt)
  {
  readEntirePointing (inp,inp.columnNumber(ntheta),inp.columnNumber(nphi),pt);
  }

/*! Appends the colatitudes and longitudes from \a pt to the columns with the
    numbers \a ctheta and \a cphi, respectively, in \a out. */
void appendPointing (iohandle &out, int ctheta, int cphi,
  const arr<pointing> &pt);

/*! Appends the colatitudes and longitudes from \a pt to the columns with the
    numbers \a ctheta and \a cphi, respectively, in \a out. Also appends the
    orientations in \a psi to the column with the number \a cpsi. */
void appendPointing (iohandle &out, int ctheta, int cphi, int cpsi,
  const arr<pointing> &pt, const arr<double> &psi);

/*! Reads the x,y,z, and w quaternion components from the columns with the
    numbers \a cx, \a cy, \a cz and \a cw in \a inp, starting at the column
    entry with the number \a offset, and returns them in \a quat. The size of
    \a quat determines the number of values read.   */
void readQuaternions (iohandle &inp, int cx, int cy, int cz, int cw,
  arr<quaternion> &quat, uint64 offset);

/*! Appends the x,y,z, and w components of \a quat to the columns with the
    numbers \a cx, \a cy, \a cz and \a cw in \a out. */
void appendQuaternions (iohandle &out, int cx, int cy, int cz, int cw,
  const arr<quaternion> &quat);

void readPointSources (const std::string &infile, const std::string &suffix,
  bool polarisation, arr<pointing> &pos, arr<double> &flux,
  arr<double> &polang, arr<double> &polfrac, arr<std::string> &name,
  int nshares, int myshare);

void readPointSourcesNew (const std::string &infile,
  bool polarisation, arr<pointing> &pos, arr<double> &flux,
  arr<double> &polang, arr<double> &polfrac, arr<std::string> &name,
  int nshares, int myshare);

#endif
