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

#include "iohandle_current.h"

using namespace std;

iohandleManager HandleManager;

namespace {

class handleRegistrator
  {
  public:
    handleRegistrator()
      {
      HandleManager.registerType("FITS",iohandle_fits::makeHandle);
#if (defined(USE_TOODI) || defined(USE_HFIDMC))
      HandleManager.registerDefaultType
        (iohandle_dmc::Type(),iohandle_dmc::makeHandle);
#else
      HandleManager.registerDefaultType("FITS",iohandle_fits::makeHandle);
#endif
      }
  };

handleRegistrator localReg;

#if (defined(USE_TOODI) || defined(USE_HFIDMC))
#include "iohandle_dmc.h"
typedef iohandle_dmc iohc;
#else
typedef iohandle_fits iohc;
#endif

iohc::Manager *iohcmng=0;

} // unnamed namespace

iohandle_current::Manager::Manager(const string &properties)
  {
  planck_assert(iohcmng==0,"iohandle manager already initialised");
  iohcmng = new iohc::Manager(properties);
  }

iohandle_current::Manager::Manager(int argc, const char **argv)
  {
  planck_assert(iohcmng==0,"iohandle manager already initialised");
  iohcmng = new iohc::Manager(argc,argv);
  }

iohandle_current::Manager::~Manager()
  { delete iohcmng; }

paramfile iohandle_current::Manager::getParams (bool verbose) const
  { return iohcmng->getParams(verbose); }

//static
string iohandle_current::Type()
  { return iohc::Type(); }
