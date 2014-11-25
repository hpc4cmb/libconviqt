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

/*! \file focalplane_db.h
 *  Class for holding information on the detectors of the Planck satellite
 *
 *  Copyright (C) 2003-2011 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_FOCALPLANE_DB_H
#define PLANCK_FOCALPLANE_DB_H

#include <string>
#include <map>
#include "iohandle.h"

class paramfile;

/*! Class for reading focal plane database info from a FITS file and making
    it accessible for modules. */
class focalplane_db
  {
  private:
    std::map<std::string,int> det2line;
    safe_ptr<iohandle> handle;
    std::string version_;
    double theta_b_;

    void init(const std::string &fname);

  public:
    /*! Reads the focal plane database from the file specified under
        \c focalplane_db in \a parfile. */
    focalplane_db (const paramfile &parfile);
    focalplane_db (const std::string &fname);
    ~focalplane_db();

    /*! Returns the boresight angle (in degrees). */
    double theta_b() const { return theta_b_; }
    /*! Returns the version string of the focal plane database. */
    const std::string &version() const { return version_; }

    /*! Returns \a true if the requested quantity is present, else \a false. */
    bool quantityPresent (const std::string &quantity) const
      { return handle->columnPresent(quantity); }

    /*! Returns the quantity \a quantity for the detector \a detname in
        \a result. */
    template<typename T> void getValue (const std::string &detname,
      const std::string &quantity, T &result) const;
    /*! Returns the quantity \a quantity for the detector \a detname. */
    template<typename T> T getValue (const std::string &detname,
      const std::string &quantity) const
      {
      T result;
      getValue(detname, quantity, result);
      return result;
      }
  };

#endif
