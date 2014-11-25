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

/*
 *  Utilities for accessing ephemeris data
 *
 *  Copyright (C) 2009 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#ifndef PLANCK_EPHEMERIDES_H
#define PLANCK_EPHEMERIDES_H

#include "iohandle.h"
#include "arr.h"
#include "vec3.h"

/*! Class providing ephemeris information. All times are in seconds since
    Jan 1, 1958, 00:00:00. All quantities are relative to the Planck satellite
    at the requested time. */
class ephemerides
  {
  public:
//    ephemerides (const std::string &obj);
    virtual ~ephemerides () {}

    /*! Returns a list of all bodies stored in the ephemeris object. */
    virtual const arr<std::string> &getBodies () const = 0;
    /*! All data requests following this call will refer to \a name. */
    virtual void loadBody (const std::string &name) = 0;
    /*! Returns a vector from the satellite to the current object
        (ecliptic coordinates, measured in metres). */
    virtual vec3 posRelSat_m (double time) const = 0;
    /*! Returns the distance of the current object to the Sun
        (measured in metres). */
    virtual double dSun_m (double time) const = 0;
    /*! Returns the angular diameter of the current object as seen from
        the satellite (measured in arcmin). */
    virtual double angdiam_arcmin (double time) const = 0;
    /*! Returns the angle between the satellite and the Sun, as seen from the
        current object (measured in radians). */
    virtual double angle_S_T_O_rad (double time) const = 0;
  };

/*! Factory function generating an \a ephemerides object associated with
    \a obj.*/
ephemerides *getEphemerides (const std::string &obj);

#endif
