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

/*! \file io_utils.cc
 *  Convenience functions for I/O with iohandles.
 *
 *  Copyright (C) 2009-2012 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include "io_utils.h"
#include "lsconstants.h"
#include "share_utils.h"

using namespace std;

void readPointing (iohandle &inp, int ctheta, int cphi,
  arr<pointing> &pt, uint64 offset)
  {
  chunkMaker cm(pt.size(),inp.efficientChunkSize(ctheta));
  uint64 start, size;
  arr<double> tmp;
  while (cm.getNext(start,size))
    {
    tmp.alloc(size);
    inp.readColumn(ctheta,tmp,start+offset);
    for (tsize i=0; i<size; ++i) pt[start+i].theta=tmp[i];
    inp.readColumn(cphi,tmp,start+offset);
    for (tsize i=0; i<size; ++i) pt[start+i].phi=tmp[i];
    }
  }

void readPointing (iohandle &inp, int ctx, int cpx, int cty, int cpy,
  int ctz, int cpz, arr<pointing> &px, arr<pointing> &py, arr<pointing> &pz,
  uint64 offset)
  {
  planck_assert (multiequal(px.size(),py.size(),pz.size()),
    "pointing arrays have different length");
  chunkMaker cm(px.size(),inp.efficientChunkSize(ctx));
  uint64 start, size;
  arr<double> tmp;
  while (cm.getNext(start,size))
    {
    tmp.alloc(size);
    inp.readColumn(ctx,tmp,start+offset);
    for (tsize i=0; i<size; ++i) px[start+i].theta=tmp[i];
    inp.readColumn(cpx,tmp,start+offset);
    for (tsize i=0; i<size; ++i) px[start+i].phi=tmp[i];
    inp.readColumn(cty,tmp,start+offset);
    for (tsize i=0; i<size; ++i) py[start+i].theta=tmp[i];
    inp.readColumn(cpy,tmp,start+offset);
    for (tsize i=0; i<size; ++i) py[start+i].phi=tmp[i];
    inp.readColumn(ctz,tmp,start+offset);
    for (tsize i=0; i<size; ++i) pz[start+i].theta=tmp[i];
    inp.readColumn(cpz,tmp,start+offset);
    for (tsize i=0; i<size; ++i) pz[start+i].phi=tmp[i];
    }
  }


void readEntirePointing (iohandle &inp, int ctheta, int cphi,
  arr<pointing> &pt)
  {
  planck_assert(inp.columnLength(ctheta)==inp.columnLength(cphi),
    "readEntirePointing(): columns have different length");
  pt.alloc(inp.columnLength(ctheta));
  readPointing (inp,ctheta,cphi,pt,0);
  }

void appendPointing (iohandle &out, int ctheta, int cphi,
  const arr<pointing> &pt)
  {
  chunkMaker cm(pt.size(),out.efficientChunkSize(ctheta));
  uint64 start, size;
  arr<double> tmp;
  while (cm.getNext(start,size))
    {
    tmp.alloc(size);
    for (tsize i=0; i<size; ++i) tmp[i]=pt[start+i].theta;
    out.appendColumn(ctheta,tmp);
    for (tsize i=0; i<size; ++i) tmp[i]=pt[start+i].phi;
    out.appendColumn(cphi,tmp);
    }
  }

void appendPointing (iohandle &out, int ctheta, int cphi, int cpsi,
  const arr<pointing> &pt, const arr<double> &psi)
  {
  planck_assert (pt.size()==psi.size(), "array size mismatch");
  chunkMaker cm(pt.size(),out.efficientChunkSize(ctheta));
  uint64 start, size;
  arr<double> tmp;
  while (cm.getNext(start,size))
    {
    tmp.alloc(size);
    for (tsize i=0; i<size; ++i) tmp[i]=pt[start+i].theta;
    out.appendColumn(ctheta,tmp);
    for (tsize i=0; i<size; ++i) tmp[i]=pt[start+i].phi;
    out.appendColumn(cphi,tmp);
    out.appendColumnRaw(cpsi,&psi[start],size);
    }
  }

void readQuaternions (iohandle &inp, int cx, int cy, int cz, int cw,
  arr<quaternion> &quat, uint64 offset)
  {
  chunkMaker cm(quat.size(),inp.efficientChunkSize(cx));
  uint64 start, size;
  arr<double> tmp;
  while (cm.getNext(start,size))
    {
    tmp.alloc(size);
    inp.readColumn(cx,tmp,start+offset);
    for (tsize i=0; i<size; ++i) quat[start+i].x=tmp[i];
    inp.readColumn(cy,tmp,start+offset);
    for (tsize i=0; i<size; ++i) quat[start+i].y=tmp[i];
    inp.readColumn(cz,tmp,start+offset);
    for (tsize i=0; i<size; ++i) quat[start+i].z=tmp[i];
    inp.readColumn(cw,tmp,start+offset);
    for (tsize i=0; i<size; ++i) quat[start+i].w=tmp[i];
    }
  }

void appendQuaternions (iohandle &out, int cx, int cy, int cz, int cw,
  const arr<quaternion> &quat)
  {
  chunkMaker cm(quat.size(),out.efficientChunkSize(cx));
  uint64 start, size;
  arr<double> tmp;
  while (cm.getNext(start,size))
    {
    tmp.alloc(size);
    for (tsize i=0; i<size; ++i) tmp[i]=quat[start+i].x;
    out.appendColumn(cx,tmp);
    for (tsize i=0; i<size; ++i) tmp[i]=quat[start+i].y;
    out.appendColumn(cy,tmp);
    for (tsize i=0; i<size; ++i) tmp[i]=quat[start+i].z;
    out.appendColumn(cz,tmp);
    for (tsize i=0; i<size; ++i) tmp[i]=quat[start+i].w;
    out.appendColumn(cw,tmp);
    }
  }

void readPointSources (const string &infile, const string &suffix,
  bool polarisation, arr<pointing> &pos, arr<double> &flux,
  arr<double> &polang, arr<double> &polfrac, arr<string> &name,
  int nshares, int myshare)
  {
  safe_ptr<iohandle> inp;
  if (HandleManager.defaultType()=="HFIDMC")
    {
    string datatype = polarisation ?
      "catalog.LS_pointsource_catalog_pol" :
      "catalog.LS_pointsource_catalog";
    inp = HandleManager.openObject(infile,datatype);
    }
  else
    inp = HandleManager.openObject(infile);

  int thetacol = inp->columnNumber("theta_ecl");
  int phicol = inp->columnNumber("phi_ecl");
  int namecol = inp->columnNumber("source_name");
  int fluxcol = inp->columnNumber(string("flux_") + suffix);
  int anglecol = polarisation ?
    inp->columnNumber(string("polangle_") + suffix) : 0;
  int percentcol = polarisation ?
    inp->columnNumber(string("polpercent_") + suffix) : 0;

  int64 lo,hi;
  calcShareGeneral(0,inp->columnLength(thetacol),nshares,myshare,lo,hi);

  tsize nsrc=hi-lo;
  flux.alloc(nsrc);
  polang.alloc(nsrc);
  polfrac.alloc(nsrc);
  name.alloc(nsrc);
  pos.alloc(nsrc);

  chunkMaker cm(nsrc,inp->efficientChunkSize(thetacol));
  uint64 start, size;
  arr<double> tmp;
  while (cm.getNext(start,size))
    {
    tmp.alloc(size);

    inp->readColumn(thetacol,tmp,lo+start);
    for (tsize m=0; m<size; ++m) pos[m+start].theta=tmp[m];
    inp->readColumn(phicol,tmp,lo+start);
    for (tsize m=0; m<size; ++m) pos[m+start].phi=tmp[m];
    inp->readColumnRaw(fluxcol,&flux[start],size,lo+start);
    inp->readColumnRaw(namecol,&name[start],size,lo+start);

    if (polarisation)
      {
      inp->readColumnRaw(anglecol,&polang[start],size,lo+start);
      inp->readColumnRaw(percentcol,&polfrac[start],size,lo+start);
      }
    }

  if (polarisation)
    for (tsize m=0; m<nsrc; ++m)
      { polang[m]*=degr2rad; polfrac[m]*=0.01; }
  else
    { polang.fill(0); polfrac.fill(0); }
  }

void readPointSourcesNew (const string &infile,
  bool polarisation, arr<pointing> &pos, arr<double> &flux,
  arr<double> &polang, arr<double> &polfrac, arr<string> &name,
  int nshares, int myshare)
  {
  safe_ptr<iohandle> inp (HandleManager.openObject
    (infile,polarisation ? "catalog.LS_pointsource_catalog_pol_new"
                         : "catalog.LS_pointsource_catalog_new"));

  int thetacol = inp->columnNumber("theta_ecl");
  int phicol = inp->columnNumber("phi_ecl");
  int namecol = inp->columnNumber("source_name");
  int fluxcol = inp->columnNumber("flux");
  int anglecol = polarisation ? inp->columnNumber("polangle") : 0;
  int percentcol = polarisation ? inp->columnNumber("polpercent") : 0;

  int64 lo,hi;
  calcShareGeneral(0,inp->columnLength(thetacol),nshares,myshare,lo,hi);

  tsize nsrc=hi-lo;
  flux.alloc(nsrc);
  polang.alloc(nsrc);
  polfrac.alloc(nsrc);
  name.alloc(nsrc);
  pos.alloc(nsrc);

  chunkMaker cm(nsrc,inp->efficientChunkSize(thetacol));
  uint64 start, size;
  arr<double> tmp;
  while (cm.getNext(start,size))
    {
    tmp.alloc(size);

    inp->readColumn(thetacol,tmp,lo+start);
    for (tsize m=0; m<size; ++m) pos[m+start].theta=tmp[m];
    inp->readColumn(phicol,tmp,lo+start);
    for (tsize m=0; m<size; ++m) pos[m+start].phi=tmp[m];
    inp->readColumnRaw(fluxcol,&flux[start],size,lo+start);
    inp->readColumnRaw(namecol,&name[start],size,lo+start);

    if (polarisation)
      {
      inp->readColumnRaw(anglecol,&polang[start],size,lo+start);
      inp->readColumnRaw(percentcol,&polfrac[start],size,lo+start);
      }
    }

  if (polarisation)
    for (tsize m=0; m<nsrc; ++m)
      { polang[m]*=degr2rad; polfrac[m]*=0.01; }
  else
    { polang.fill(0); polfrac.fill(0); }
  }
