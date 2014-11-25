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
 *  Class for holding information on the detectors of the Planck satellite
 *
 *  Copyright (C) 2003-2011 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "focalplane_db.h"
#include "iohandle.h"
#include "paramfile.h"

using namespace std;

void focalplane_db::init(const string &fname)
  {
// FIXME: in principle this should better be opened without explicitly
//        specifying a DDL type. HFIDMC currently does not allow this,
//        however ...
  handle = HandleManager.openObject(fname,"focalplane.LS_focalplanedb");

  handle->getKey("fpdb_version", version_);
//  cout << "Focalplane DB: found version " << version_ << endl;
  theta_b_ = handle->getKey<double>("theta_b");

  int c1  = handle->columnNumber("detector");
  arr<string> name;
  handle->readEntireColumn(c1,name);
  for (tsize i=0; i<name.size(); ++i)
    det2line[name[i]] = i;
  }

//FIXME: temporary
focalplane_db::focalplane_db (const paramfile &parfile)
  {
  string fname = parfile.find<string>("focalplane_db");
  init(fname);
  }

focalplane_db::focalplane_db (const string &fname)
  {
  init(fname);
  }

focalplane_db::~focalplane_db () {}

template<typename T> void focalplane_db::getValue (const string &detname,
  const string &quantity, T &result) const
  {
  map<string,int>::const_iterator line = det2line.find(detname);
  planck_assert (line!=det2line.end(), "detector " + detname + " not found.");
  handle->readColumn(quantity,result,line->second);
  }

template void focalplane_db::getValue (const string &detname,
  const string &quantity, int &result) const;
template void focalplane_db::getValue (const string &detname,
  const string &quantity, double &result) const;
