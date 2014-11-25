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

/*! \file iohandle_helpers.h
 *
 *  Convenience functions for copying information between iohandles.
 *
 *  Copyright (C) 2009-2011 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#ifndef PLANCK_IOHANDLE_HELPERS
#define PLANCK_IOHANDLE_HELPERS

#include <algorithm>
#include "iohandle.h"
#include "datatypes.h"
#include "share_utils.h"

namespace {

template<typename T> void copyColumn (const iohandle &inp, iohandle &out,
  const std::string &name_in, const std::string &name_out, uint64 startpos,
  uint64 len)
  {
  arr<T> data;
  tsize icol = inp.columnNumber(name_in),
        ocol = out.columnNumber(name_out);
  if (len==0) len=inp.columnLength(icol);

  chunkMaker cm (len, std::min(inp.efficientChunkSize(icol),
                               out.efficientChunkSize(ocol)));
  uint64 offset, nsamples;
  while (cm.getNext(offset,nsamples))
    {
    data.alloc(nsamples);
    inp.readColumn(icol,data,startpos+offset);
    out.appendColumn(ocol,data);
    }
  }

template<typename T> void copyColumn (const iohandle &inp, iohandle &out,
  const std::string &name, uint64 startpos, uint64 len)
  { copyColumn<T> (inp,out,name,name,startpos,len); }

template<typename T> void copyEntireColumn (const iohandle &inp, iohandle &out,
  const std::string &name_in, const std::string &name_out)
  { copyColumn<T> (inp,out,name_in,name_out,0,0); }

template<typename T> void copyEntireColumn (const iohandle &inp, iohandle &out,
  const std::string &name)
  { copyColumn<T> (inp,out,name,name,0,0); }

template<typename T> void copyKey(iohandle &inp, iohandle &out,
  const std::string &key_in, const std::string &key_out)
  { out.setKey(key_out,inp.getKey<T>(key_in)); }

template<typename T> void copyKey(iohandle &inp, iohandle &out,
  const std::string &key)
  { copyKey<T>(inp,out,key,key); }

template<typename T> void copyKeyIfPresent(iohandle &inp, iohandle &out,
  const std::string &key_in, const std::string &key_out)
  { if (inp.keyPresent(key_in)) copyKey<T>(inp,out,key_in,key_out); }

template<typename T> void copyKeyIfPresent(iohandle &inp, iohandle &out,
  const std::string &key)
  { copyKeyIfPresent<T>(inp,out,key,key); }

} // unnamed namespace

#endif
