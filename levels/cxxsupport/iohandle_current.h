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
 *  Copyright (C) 2004, 2005, 2006, 2009, 2010 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#ifndef PLANCK_IOHANDLE_CURRENT_H
#define PLANCK_IOHANDLE_CURRENT_H

#include "iohandle_fits.h"
#include "iohandle_dmc.h"

class iohandle_current
  {
  public:
    class Manager
      {
      public:
        Manager (const std::string &properties);
        Manager (int argc, const char ** argv);
        ~Manager();
        paramfile getParams (bool verbose=true) const;
      };
    static std::string Type();
  };

#endif
