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
 *  This file contains the implementation of the FITS I/O interface
 *  used by the Planck LevelS package.
 *
 *  Copyright (C) 2005-2011 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include <sstream>
#include "iohandle_fits.h"
#include "safe_cast.h"
#include "string_utils.h"

using namespace std;

namespace {

static map<string,string> paramHolder;

bool keyReserved (const string &key)
  {
  const char *exclist[] = { "objType" };
  const tsize n_exc = sizeof(exclist)/sizeof(char *);

  for (tsize i=0; i<n_exc; ++i)
    if (equal_nocase(key,exclist[i])) return true;
  return false;
  }

bool conversionAcceptable (PDT from, PDT to)
  {
  if (from==to) return true;
  switch (to)
    {
    case PLANCK_FLOAT32:
    case PLANCK_FLOAT64:
      switch (from)
        {
        case PLANCK_FLOAT32:
        case PLANCK_FLOAT64:
        case PLANCK_INT32:
        case PLANCK_INT64:
          return true;
        default:
          return false;
        }
    case PLANCK_UINT8:
    case PLANCK_INT8:
      switch (from)
        {
        case PLANCK_INT8:
        case PLANCK_UINT8:
          return true;
        default:
          return false;
        }
    default:
      return false;
    }
  return false;
  }

} // unnamed namespace

iohandle_fits::coltype::coltype() {}
iohandle_fits::coltype::coltype (const string &nm, const string &un,
  PDT tp, int hdu_, int idx_, int repc, uint64 sz)
  : name(nm), unit(un), type(tp), hdu(hdu_), idx(idx_), width(repc), size(sz) {}

iohandle_fits::coltype::~coltype() {}

iohandle_fits::Manager::Manager (const string &properties)
  { if (properties!="") parse_file(properties,paramMap); }
iohandle_fits::Manager::Manager (int argc, const char **argv)
  {
  planck_assert(argc>=2,"incorrect command line format");
  if ((argc==2) && (string(argv[1]).find("=")==string::npos))
    {
    if (string(argv[1])!="")
      parse_file(argv[1],paramMap);
    }
  else
    parse_cmdline_equalsign(argc,argv,paramMap);
  }
iohandle_fits::Manager::~Manager ()
  {}

paramfile iohandle_fits::Manager::getParams (bool verbose) const
  {
  paramfile pfile(paramMap,verbose);
  paramHolder = pfile.getParams();
  return pfile;
  }

//static
iohandle_fits::ddlmap iohandle_fits::ddltable;
//static
map<string,string> iohandle_fits::ddl_hierarchy;

//static
void iohandle_fits::init_ddltables()
  {
#include "ddl.cc"
  istringstream ddl(ddlstring);
  while(!ddl.eof())
    {
    string tn, cn, un, tp;
    int ncol, hdu, coln, repc;
    ddl >> tn >> ncol;
    vector<coltype> vcol;
    vcol.reserve(ncol);
    if (!ddl) continue;
    for (int i=0; i<ncol; ++i)
      {
      ddl >> cn >> un >> tp >> hdu >> coln >> repc;
      planck_assert(ddl,"read error during DDL parsing");
      PDT itype = (tp=="BYTE") ? PLANCK_INT8 : string2type(tp);
      vcol.push_back(coltype(cn,un,itype,hdu,coln,repc));
      }
    ddltable[tn].swap(vcol);
    }
  istringstream hier(ddl_inheritance);
  while(!hier.eof())
    {
    string child, parent;
    hier >> child >> parent;
    if (!hier) continue;
    ddl_hierarchy[child]=parent;
    }
  }

void iohandle_fits::assert_connected() const
  { planck_assert (connected, "not connected to an object"); }
void iohandle_fits::assert_write() const
  {
  assert_connected();
  planck_assert (!readonly, "trying to write to a read-only object");
  }
void iohandle_fits::assert_read() const
  {
  assert_connected();
  planck_assert (readonly, "reading not allowed from this object");
  }
void iohandle_fits::clean_all()
  {
  if (!connected) return;
  handle.close();
  col.clear();
  type_ = "INVALID";
  connected=false;
  readonly=true;
  }

void iohandle_fits::fill_cache()
  {
  assert_read();

  handle.goto_hdu(2);
  if (handle.key_present("objType"))
    type_ = handle.get_key<string> ("objType");
  else
    type_ = "UNKNOWN";

  for (int m=2; m<=handle.num_hdus(); ++m)
    {
    handle.goto_hdu(m);
    for (int n=1; n<=handle.ncols(); ++n)
      {
      coltype tcol;
      tcol.name = handle.colname(n);
      tcol.unit = handle.colunit(n);
      tcol.width = safe_cast<int>(handle.repcount(n));
      tcol.type = handle.coltype(n);
      tcol.hdu = m;
      tcol.idx = n;
      if (handle.coltype(n)==PLANCK_STRING)
        tcol.size = uint64(handle.nrows());
      else
        tcol.size = uint64(handle.nrows())*handle.repcount(n);
      col.push_back(tcol);
      }
    }
  }

namespace {

void consistency_warning (int hdu, int idx, const char *thing,
  const string &expect, const string &found)
  {
  cerr << "Warning: HDU " << hdu << ", column " << idx << ": expected "
       << thing << " '" << expect << "', but found '" << found << "'." << endl;

  }

} // unnamed namespace

void iohandle_fits::check_consistency (const string &fname) const
  {
  assert_read();

  handle.goto_hdu(2);
  if (handle.key_present("objType"))
    {
    string tp = handle.get_key<string> ("objType");
    string tp2=tp;
    while ((tp2!=type_) && (tp2!=""))
      {
      map<string,string>::const_iterator it=ddl_hierarchy.find(tp2);
      tp2 = (it==ddl_hierarchy.end()) ? "" : it->second;
      }
    if (type_!=tp2)
      cerr << "\n*\n"
           << "* WARNING: File '" << fname <<"':\n"
           << "* Expected object type '" << type_ << "', but found\n"
           << "* '" << tp << "', which doesn't seem to be a parent type.\n"
           << "* The module would crash here when running in DMC mode.\n*\n\n";
    }
  else
    cerr << "WARNING: FITS file does not contain type information" << endl;

  for (tsize m=0; m<col.size(); ++m)
    {
    int hdu=col[m].hdu, idx=col[m].idx;
    planck_assert(handle.num_hdus()>=hdu,"FITS file has too few HDUs");
    handle.goto_hdu(hdu);
    planck_assert(handle.ncols()>=idx,"FITS HDU has too few columns");
    if (col[m].name!=handle.colname(idx))
      consistency_warning (hdu,idx,"name",col[m].name,handle.colname(idx));
    if (col[m].type!=handle.coltype(idx))
      {
      consistency_warning (hdu,idx,"type",type2string(col[m].type),
                           type2string(handle.coltype(idx)));
      planck_assert(conversionAcceptable(col[m].type,handle.coltype(idx))
                  &&conversionAcceptable(handle.coltype(idx),col[m].type),
        "incompatible column types");
      }
    if (col[m].unit!="[arbitrary]" && col[m].unit!=handle.colunit(idx))
      consistency_warning (hdu,idx,"unit",col[m].unit, handle.colunit(idx));
    if (col[m].width!=handle.repcount(idx))
      consistency_warning (hdu,idx,"repcount",dataToString(col[m].width),
        dataToString(handle.repcount(idx)));
    }
  }

void iohandle_fits::create_hdus_up_to (int hdu) const
  {
  assert_write();

  for (int m=maxhdu_written+1; m<=hdu; ++m)
    {
    handle.goto_hdu(m-1);
    int maxidx = 0;
    for (tsize n=0; n<col.size(); ++n)
      if ((col[n].hdu==m)&&(col[n].idx>maxidx)) maxidx=col[n].idx;

    vector<fitscolumn> cols(maxidx);
    for (tsize n=0; n<col.size(); ++n)
      {
      if (col[n].hdu==m)
        cols[col[n].idx-1]=fitscolumn (col[n].name,col[n].unit,col[n].width,
                                       col[n].type);
      }
    handle.insert_bintab(cols);
    }
  if (hdu>maxhdu_written) maxhdu_written = hdu;
  }

void iohandle_fits::populate_fits_file() const
  {
  assert_write();

  //FIXME Doesn't work in mixed scenario with the DMC providing the parameters!
  handle.goto_hdu(1);
  for (map<string,string>::const_iterator iter=paramHolder.begin();
       iter!=paramHolder.end(); ++iter)
    handle.set_key(iter->first, iter->second);

  maxhdu_written = 1;
  create_hdus_up_to(2);
  handle.goto_hdu(2);
  handle.set_key("objType",type_);
  }

void iohandle_fits::fill_cache_from_ddl (const string &type)
  {
#pragma omp critical (fits_ddl_init)
{
  if (ddltable.size()==0)
    {
    init_ddltables();
    planck_assert (ddltable.size()>0, "no types in DDL");
    }
}
  ddlmap::const_iterator iter = ddltable.find(type);
  planck_assert (iter!=ddltable.end(), "type '"+type+"' not found in DDL");
  type_ = type;
  col = iter->second;
  }

void iohandle_fits::update_colsizes()
  {
  assert_read();

  for (tsize m=0; m<col.size(); ++m)
    {
    handle.goto_hdu(col[m].hdu);
    if (handle.coltype(col[m].idx)==PLANCK_STRING)
      col[m].size = uint64(handle.nrows());
    else
      col[m].size = uint64(handle.nrows())*handle.repcount(col[m].idx);
    }
  }

//static
iohandle *iohandle_fits::makeHandle()
  { return new iohandle_fits; }

iohandle_fits::iohandle_fits ()
  : connected(false) {}
// virtual
iohandle_fits::~iohandle_fits()
  { clean_all(); }

// virtual
void iohandle_fits::create (const string &name, const string &type)
  {
  clean_all();
  fill_cache_from_ddl(type);
  handle.create(name);
  connected = true;
  readonly = false;
  populate_fits_file();
  }
// virtual
void iohandle_fits::open (const string &name)
  {
  clean_all();
  handle.open(name);
  connected = true;
  readonly = true;

// Do not perform a check for objType in the file.
// Simply take the column information from the FITS headers.
  fill_cache();
  update_colsizes();
  }
// virtual
void iohandle_fits::open (const string &name, const string &type)
  {
  clean_all();
  handle.open(name);
  connected = true;
  readonly = true;

  fill_cache_from_ddl(type);
  check_consistency(name);
  update_colsizes();
  }
// virtual
void iohandle_fits::close()
  {
  if (!readonly) // if writing, create potentially empty trailing HDUs
    {
    int maxhdu=2;
    for (tsize n=0; n<col.size(); ++n)
      maxhdu=max(maxhdu,col[n].hdu);
    create_hdus_up_to (maxhdu);
    }
  clean_all();
  }
// virtual
void iohandle_fits::flush()
  { /* do nothing (at least for now) */ }

//virtual
string iohandle_fits::getID (const string &, const string &) const
  {
  planck_fail ("iohandle_fits::getID() called");
  }

//virtual
void iohandle_fits::deleteObject(const string &name) const
  {
  fitshandle::delete_file (name);
  }

// virtual
string iohandle_fits::objectName() const
  {
  assert_connected();
  return handle.fileName();
  }
// virtual
string iohandle_fits::objectType() const
  {
  assert_connected();
  return type_;
  }
// virtual
int iohandle_fits::numColumns() const
  {
  assert_connected();
  return col.size();
  }

// virtual
bool iohandle_fits::keyPresent (const string &key) const
  {
  assert_read();
  handle.goto_hdu(2);
  return handle.key_present(key);
  }
//virtual
void iohandle_fits::setKeyVoid (const string &key, const void *value, PDT type)
  {
  assert_write();

  if (keyReserved(key))
    {
    cerr << "Warning: key '"+key+"' is reserved and will not be written!"
         << endl;
    return;
    }

  handle.goto_hdu(2);
  handle.set_key_void(key, value, type);
  }
//virtual
void iohandle_fits::setKeyTypeless (const string &key, const string &value)
  {
  if (keyReserved(key))
    {
    cerr << "Warning: key '"+key+"' is reserved and will not be written!"
         << endl;
    return;
    }
  setKey(key,value);
  }
//virtual
void iohandle_fits::getKeyVoid (const string &key, void *value, PDT type) const
  {
  assert_read();
  handle.goto_hdu(2);
  handle.get_key_void(key, value, type);
  }
//virtual
string iohandle_fits::getKeyTypeless (const string &key) const
  {
  return getKey<string>(key);
  }

// virtual
void iohandle_fits::deleteKey (const string &key)
  {
  assert_write();
  handle.goto_hdu(2);
  handle.delete_key(key);
  }

//virtual
void iohandle_fits::getAllKeys (vector<string> &keys) const
  {
  assert_read();
  keys.clear();
  vector<string> keys2;
  handle.goto_hdu(2);
  handle.get_all_keys(keys2);

  for (tsize i=0; i<keys2.size(); ++i)
    if (!keyReserved(keys2[i]))
      keys.push_back(keys2[i]);
  }

//virtual
void iohandle_fits::query (const string &, vector<string> &) const
  {
  planck_fail ("no queries supported for FITS handles");
  }

//virtual
bool iohandle_fits::columnPresent (const string &colname) const
  {
  assert_connected();
  for (int m=0; m<numColumns(); ++m)
    if (col[m].name==colname) return true;
  return false;
  }

//virtual
int iohandle_fits::columnNumber (const string &colname) const
  {
  assert_connected();
  for (int m=0; m<numColumns(); ++m)
    if (col[m].name==colname) return m;
  planck_fail ("Column name '" + colname + "' not found.");
  }
// virtual
string iohandle_fits::columnName (int colnum) const
  {
  assert_connected();
  return col[colnum].name;
  }
// virtual
PDT iohandle_fits::columnType (int colnum) const
  {
  assert_connected();
  return col[colnum].type;
  }
// virtual
string iohandle_fits::columnUnit (int colnum) const
  {
  assert_connected();
  return col[colnum].unit;
  }
//  virtual
uint64 iohandle_fits::columnLength (int colnum) const
  {
  assert_read();
  return col[colnum].size;
  }
//  virtual
uint64 iohandle_fits::efficientChunkSize (int colnum) const
  {
  assert_connected();
  if (!readonly) create_hdus_up_to(col[colnum].hdu);
  handle.goto_hdu(col[colnum].hdu);
  return handle.efficientChunkSize(col[colnum].idx);
  }

// virtual
void iohandle_fits::readColumnRawVoid
  (int colnum, void *data, PDT type, long num, uint64 offset) const
  {
  assert_read();
  planck_assert (conversionAcceptable(col[colnum].type,type),
    "invalid data type conversion requested");
  handle.goto_hdu(col[colnum].hdu);
  handle.read_column_raw_void(col[colnum].idx, data, type, num, offset);
  }

// virtual
void iohandle_fits::writeColumnRawVoid
  (int colnum, const void *data, PDT type, long num, uint64 offset)
  {
  assert_write();
  planck_assert (conversionAcceptable(type,col[colnum].type),
    "invalid data type conversion requested");
  create_hdus_up_to(col[colnum].hdu);
  handle.goto_hdu(col[colnum].hdu);
  handle.write_column_raw_void(col[colnum].idx, data, type, num, offset);
  if ((offset+num)>col[colnum].size) col[colnum].size = offset+num;
  }

// virtual
void iohandle_fits::appendColumnRawVoid
  (int colnum, const void *data, PDT type, long num)
  {
  assert_write();
  planck_assert (conversionAcceptable(type,col[colnum].type),
    "invalid data type conversion requested");
  create_hdus_up_to(col[colnum].hdu);
  handle.goto_hdu(col[colnum].hdu);
  handle.write_column_raw_void(col[colnum].idx, data, type, num,
    col[colnum].size);
  col[colnum].size += num;
  }
