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
 *  Copyright (C) 2004, 2005, 2006, 2009 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "iohandle.h"

using namespace std;

void iohandleManager::getTypeAndName (const string &fullname, string &type,
  string &name) const
  {
  for (maptype::const_iterator iter=dict.begin(); iter!=dict.end(); ++iter)
    {
    string tofind = iter->first + ":";
    if (fullname.find(tofind)==0) //found at beginning
      {
      type = iter->first;
      name = fullname;
      name.erase(0,tofind.size());
      return;
      }
    }
  name = fullname;
  type = deflt;
  }

iohandle *iohandleManager::makeHandle(const string &type) const
  {
  planck_assert (dict.size()>0, "makeHandle: no types registered");
  maptype::const_iterator iter = dict.find(type);
  planck_assert (iter!=dict.end(),
    "makeHandle(): could not find type '"+type+"'");
  return (*(iter->second))();
  }
iohandle *iohandleManager::makeHandle() const
  {
  planck_assert (deflt!="", "makeHandle: no default type provided");
  return makeHandle(deflt);
  }
iohandle *iohandleManager::createObject(const string &fullname,
  const string &ddltype) const
  {
  string type, name;
  getTypeAndName (fullname, type, name);
  iohandle *hnd = makeHandle(type);
  hnd->create (name, ddltype);
  return hnd;
  }
iohandle *iohandleManager::openObject(const string &fullname,
  const string &ddltype) const
  {
  string type, name;
  getTypeAndName (fullname, type, name);
  iohandle *hnd = makeHandle(type);
  hnd->open (name, ddltype);
  return hnd;
  }
iohandle *iohandleManager::openObject(const string &fullname) const
  {
  string type, name;
  getTypeAndName (fullname, type, name);
  iohandle *hnd = makeHandle(type);
  hnd->open (name);
  return hnd;
  }
void iohandleManager::deleteObject(const string &fullname) const
  {
  string type, name;
  getTypeAndName (fullname, type, name);
  safe_ptr<iohandle> hnd(makeHandle(type));
  hnd->deleteObject (name);
  }
string iohandleManager::getID(string &fullname, const string &ddltype)
  {
  string type, name;
  getTypeAndName (fullname, type, name);
  safe_ptr<iohandle> hnd(makeHandle(type));
  return hnd->getID(name, ddltype);
  }
