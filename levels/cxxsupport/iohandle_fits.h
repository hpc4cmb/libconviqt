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
 *  This file contains the declaration of the FITS I/O interface
 *  used by the Planck LevelS package.
 *
 *  Copyright (C) 2005-2011 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#ifndef PLANCK_IOHANDLE_FITS_H
#define PLANCK_IOHANDLE_FITS_H

#include <string>
#include "iohandle.h"
#include "paramfile.h"
#include "fitshandle.h"

class iohandle_fits: public iohandle
  {
  private:
    struct coltype
      {
      std::string name, unit;
      PDT type;
      int hdu, idx, width;
      uint64 size;

      coltype();
      coltype (const std::string &nm, const std::string &un,
               PDT tp, int hdu_, int idx_, int repc, uint64 sz=0);
      ~coltype();
      };

    typedef std::map<std::string, std::vector<coltype> > ddlmap;
    static ddlmap ddltable;
    static std::map<std::string,std::string> ddl_hierarchy;
    static void init_ddltables();

    bool connected;
    bool readonly;
    mutable fitshandle handle;
    mutable int maxhdu_written;

    std::vector<coltype> col;
    std::string type_;

    void assert_connected() const;
    void assert_write() const;
    void assert_read() const;
    void clean_all();
    void check_consistency(const std::string &fname) const;
    void fill_cache();
    void create_hdus_up_to(int hdu) const;
    void populate_fits_file() const;
    void fill_cache_from_ddl (const std::string &type);
    void update_colsizes();

    virtual void setKeyVoid
      (const std::string &key, const void *value, PDT type);
    virtual void getKeyVoid
      (const std::string &key, void *value, PDT type) const;

    virtual void readColumnRawVoid
      (int colnum, void *data, PDT type, long num, uint64 offset=0) const;
    virtual void writeColumnRawVoid
      (int colnum, const void *data, PDT type, long num, uint64 offset=0);
    virtual void appendColumnRawVoid
      (int colnum, const void *data, PDT type, long num);

  public:
    class Manager
      {
      private:
        std::map<std::string,std::string> paramMap;

      public:
        Manager (const std::string &properties);
        Manager (int argc, const char **argv);
        ~Manager();
        paramfile getParams (bool verbose=true) const;
      };

    static iohandle *makeHandle();
    static std::string Type() { return "FITS"; }

    iohandle_fits ();
    virtual ~iohandle_fits();

    virtual void create (const std::string &name, const std::string &type);
    virtual void open (const std::string &name);
    virtual void open (const std::string &name, const std::string &type);
    virtual void close();
    virtual void flush();
    virtual std::string getID (const std::string &name, const std::string &type)
      const;
    virtual void deleteObject (const std::string &name) const;

    virtual std::string objectName() const;
    virtual std::string objectType() const;
    virtual int numColumns() const;

    virtual bool keyPresent (const std::string &key) const;
    virtual void setKeyTypeless (const std::string &key,
      const std::string &value);
    virtual std::string getKeyTypeless (const std::string &key) const;
    virtual void deleteKey (const std::string &key);
    virtual void getAllKeys (std::vector<std::string> &keys) const;
    virtual void query (const std::string &, std::vector<std::string> &) const;

    virtual bool columnPresent (const std::string &colname) const;
    virtual int columnNumber (const std::string &colname) const;
    virtual std::string columnName (int colnum) const;
    virtual PDT columnType (int colnum) const;
    virtual std::string columnUnit (int colnum) const;
    virtual uint64 columnLength (int colnum) const;
    virtual uint64 efficientChunkSize (int colnum) const;
  };

#endif
