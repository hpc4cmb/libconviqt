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

/*! \file mmap_handle.h
 *
 *  Copyright (C) 2012 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_MMAP_HANDLE_H
#define PLANCK_MMAP_HANDLE_H

#include "datatypes.h"
#include "error_handling.h"

class mmap_handle
  {
  private:
    tsize mapsize;
    int fd;
    char *dataptr;

    void clean_all();

    void assert_connected() const
      { planck_assert(fd!=-1,"trying to access an unconnected mmap_handle"); }

  public:
    mmap_handle()
      : mapsize(0), fd(-1), dataptr(0) {}
    ~mmap_handle()
      { clean_all(); }

    void create (const std::string &fname, tsize nbytes);
    void openRead (const std::string &fname);
    void openReadWrite (const std::string &fname);
    void close() { clean_all(); }

    char *data() const { assert_connected(); return dataptr; }
    tsize size() const { assert_connected(); return mapsize; }
  };

#endif
