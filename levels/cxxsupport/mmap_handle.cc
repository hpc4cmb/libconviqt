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

/*! \file mmap_handle.cc
 *
 *  Copyright (C) 2012 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>

#include "mmap_handle.h"

using namespace std;

const char *MapFailed=reinterpret_cast<char *>(-1);

void mmap_handle::clean_all()
  {
  if (fd==-1) return;
  if (dataptr!=MapFailed)
    planck_assert(munmap(dataptr,mapsize)!=-1,"Error during munmap()");
  planck_assert(::close(fd)!=-1,"Error during close()");
  fd=-1;
  mapsize=0;
  dataptr=0;
  }

void mmap_handle::create (const string &fname, tsize nbytes)
  {
  clean_all();
  mapsize=nbytes;
  fd = open(fname.c_str(), O_RDWR | O_CREAT | O_TRUNC, mode_t(0600));
  planck_assert(fd!=-1,"Error during open()");
  planck_assert(lseek(fd, mapsize-1, SEEK_SET)!=-1,"Error during lseek()");
  planck_assert(write(fd, "", 1)!=-1,"Error during write()");
  dataptr = static_cast<char *>(mmap(0, mapsize, PROT_READ | PROT_WRITE,
    MAP_SHARED, fd, 0));
  planck_assert(dataptr!=MapFailed,"Error during mmap()");
  }

void mmap_handle::openRead (const string &fname)
  {
  clean_all();
  fd = open(fname.c_str(), O_RDONLY);
  planck_assert(fd!=-1,"Error during open()");
  struct stat sttmp;
  fstat(fd,&sttmp);
  mapsize=sttmp.st_size;
  dataptr = static_cast<char *>(mmap(0, mapsize, PROT_READ, MAP_SHARED, fd, 0));
  planck_assert(dataptr!=MapFailed,"Error during mmap()");
  }

void mmap_handle::openReadWrite (const string &fname)
  {
  clean_all();
  fd = open(fname.c_str(), O_RDWR);
  planck_assert(fd!=-1,"Error during open()");
  struct stat sttmp;
  fstat(fd,&sttmp);
  mapsize=sttmp.st_size;
  dataptr = static_cast<char *>(mmap(0, mapsize, PROT_READ | PROT_WRITE,
    MAP_SHARED, fd, 0));
  planck_assert(dataptr!=MapFailed,"Error during mmap()");
  }
