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
 * ls_mpi_support.h, mpi_support.cc and error_handling.cc have been modified to support libconviqt:
 *   - there is no longer a static instance, extern MPI_Manager mpiMgr, instead, 
 *     calling codes must instantiate their own managers and optionally supply the
 *     communicator
 *   - destroying an MPI_Manager instance no longer calls MPI_Finalize
 * 2014-12-01 - Reijo Keskitalo 
 */

/*
 *  libcxxsupport is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Copyright (C) 2009-2011 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef LEVELS_MPI_SUPPORT_H
#define LEVELS_MPI_SUPPORT_H

#ifdef USE_MPI
#include "mpi.h"
#else
typedef int MPI_Comm;
#define MPI_COMM_WORLD 0
#endif

#include "ls_datatypes.h"
#include "ls_arr.h"
#include "ls_share_utils.h"

namespace levels {

class MPI_Manager
  {
  public:
    enum redOp { Sum, Min, Max, Prod };

  private:
    void gatherv_helper1_m (int nval_loc, arr<int> &nval, arr<int> &offset,
      int &nval_tot) const;
    void gatherRawVoid (const void *in, tsize num, void *out, NDT type,
      int root=0) const;
    void gathervRawVoid (const void *in, tsize num, void *out,
      const int *nval, const int *offset, NDT type) const;
    void all2allv_easy_prep (tsize insz, const arr<int> &numin, arr<int> &disin,
      arr<int> &numout, arr<int> &disout) const;

    int num_ranks_, rank_;
    MPI_Comm comm_;

  public:
    MPI_Manager( MPI_Comm comm=MPI_COMM_WORLD );
    ~MPI_Manager();

    void abort() const;

    int num_ranks() const { return num_ranks_; }
    int rank() const { return rank_; }
    bool master() const { return rank_==0; }
    MPI_Comm comm() const { return comm_; }

    void barrier() const;
    double Wtime();

    void calcShare (int64 glo, int64 ghi, int64 &lo, int64 &hi) const
      { calcShareGeneral(glo,ghi,num_ranks_,rank_,lo,hi); }

    void sendRawVoid (const void *data, NDT type, tsize num, tsize dest) const;
    template<typename T> void sendRaw (const T *data, tsize num, tsize dest)
      const
      { sendRawVoid(data, nativeType<T>(), num, dest); }
    template<typename T> void send (const arr<T> &data, tsize dest) const
      { sendRaw(&data[0], data.size(), dest); }
    template<typename T> void send (const T &data, tsize dest) const
      { sendRaw(&data, 1, dest); }

    void recvRawVoid (void *data, NDT type, tsize num, tsize src) const;
    template<typename T> void recvRaw (T *data, tsize num, tsize src) const
      { recvRawVoid(data, nativeType<T>(), num, src); }
    template<typename T> void recv (arr<T> &data, tsize src) const
      { recvRaw(&data[0], data.size(), src); }
    template<typename T> void recv (T &data, tsize src) const
      { recvRaw(&data, 1, src); }

    void sendrecvRawVoid (const void *sendbuf, tsize sendcnt,
      tsize dest, void *recvbuf, tsize recvcnt, tsize src, NDT type) const;
    template<typename T> void sendrecvRaw (const T *sendbuf, tsize sendcnt,
      tsize dest, T *recvbuf, tsize recvcnt, tsize src) const
      {
      sendrecvRawVoid(sendbuf, sendcnt, dest, recvbuf, recvcnt, src,
                      nativeType<T>());
      }
    template<typename T> void sendrecv (const arr<T> &sendbuf, tsize dest,
      arr<T> &recvbuf, tsize src) const
      {
      sendrecvRaw(&sendbuf[0], sendbuf.size(), dest,
                  &recvbuf[0], recvbuf.size(), src);
      }
    template<typename T> void sendrecv (const T &sendval, tsize dest,
      T &recvval, tsize src) const
      { sendrecvRaw(&sendval, 1, dest, &recvval, 1, src); }

    void sendrecv_replaceRawVoid (void *data, NDT type, tsize num,
      tsize dest, tsize src) const;
    template<typename T> void sendrecv_replaceRaw (T *data, tsize num,
      tsize dest, tsize src) const
      { sendrecv_replaceRawVoid(data, nativeType<T>(), num, dest, src); }
    template<typename T> void sendrecv_replace (arr<T> &data, tsize dest,
      tsize src) const
      { sendrecv_replaceRaw(&data[0], data.size(), dest, src); }
    template<typename T> void sendrecv_replace (T &data, tsize dest,
      tsize src) const
      { sendrecv_replaceRaw(&data, 1, dest, src); }

    template<typename T> void gather_m (const T &in, arr<T> &out, int root=0)
      const
      {
      out.alloc(num_ranks_);
      gatherRawVoid (&in,1,&out[0],nativeType<T>(),root);
      }
    template<typename T> void gather_s (const T &in, int root=0) const
      { gatherRawVoid (&in,1,0,nativeType<T>(),root); }

    template<typename T> void gatherv_m (const arr<T> &in, arr<T> &out) const
      {
      int nval_loc = in.size(), nval_tot;
      arr<int> nval, offset;
      gatherv_helper1_m (nval_loc,nval,offset,nval_tot);
      out.alloc(nval_tot);
      gathervRawVoid (&in[0],nval_loc,&out[0],&nval[0],&offset[0],
        nativeType<T>());
      }
    template<typename T> void gatherv_s (const arr<T> &in) const
      {
      int nval_loc = in.size();
      gather_s (nval_loc);
      gathervRawVoid (&in[0],nval_loc,0,0,0,nativeType<T>());
      }

    template<typename T> void gatherv (const arr<T> &in, arr<T> &out) const
      { master() ? gatherv_m(in,out) : gatherv_s(in); }

    template<typename T> void gatherv_m (const arr2<T> &in, arr2<T> &out) const
      {
      int nval_loc = in.size(), nval_tot;
      arr<int> nval, offset;
      gatherv_helper1_m (nval_loc, nval, offset, nval_tot);
      out.alloc(nval_tot/in.size2(),in.size2());
      gathervRawVoid (&in[0][0],nval_loc,&out[0][0],&nval[0],&offset[0],
        nativeType<T>());
      }

    template<typename T> void gatherv_s (const arr2<T> &in) const
      {
      int nval_loc = in.size();
      gather_s (nval_loc);
      gathervRawVoid (&in[0][0],nval_loc,0,0,0,nativeType<T>());
      }

    void reduceRawVoid (const void *in, void *out, NDT type, tsize num,
      redOp op, int root=0) const;
    template<typename T> void reduceRaw (const T *in, T *out, tsize num,
      redOp op, int root=0) const
      { reduceRawVoid (in, out, nativeType<T>(), num, op, root); }
    template<typename T> void reduce (const arr<T> &in, arr<T> &out, redOp op,
      int root=0) const
      { (rank_==root) ? reduce_m (in, out, op) : reduce_s (in, op, root); }
    template<typename T> void reduce (const T &in, T &out, redOp op,
      int root=0) const
      { reduceRaw (&in,&out,1,op,root); }
    template<typename T> void reduce_m (const arr<T> &in, arr<T> &out,
      redOp op) const
      {
      out.alloc(in.size());
      reduceRaw (&in[0], &out[0], in.size(), op, rank_);
      }
    template<typename T> void reduce_s (const arr<T> &in, redOp op, int root=0)
      const
      { reduceRaw (&in[0], 0, in.size(), op, root); }

    void allgatherRawVoid (const void *in, void *out, NDT type, tsize num)
      const;
    template<typename T> void allgatherRaw (const T *in, T *out, tsize num)
      const
      { allgatherRawVoid (in, out, nativeType<T>(), num); }
    template<typename T> void allgather (const arr<T> &in, arr<T> &out) const
      {
      out.alloc(num_ranks_*in.size());
      allgatherRaw (&in[0], &out[0], in.size());
      }
    template<typename T> void allgather (const T &in, arr<T> &out) const
      {
      out.alloc(num_ranks_);
      allgatherRaw (&in, &out[0], 1);
      }

    void allreduceRawVoid (const void *in, void *out, NDT type, tsize num,
      redOp op) const;
    template<typename T> void allreduceRaw (const T *in, T *out, tsize num,
      redOp op) const
      { allreduceRawVoid (in, out, nativeType<T>(), num, op); }
    template<typename T> void allreduce (const arr<T> &in, arr<T> &out,
      redOp op) const
      {
      out.alloc(in.size());
      allreduceRaw (&in[0], &out[0], in.size(), op);
      }
    template<typename T> void allreduce (const T &in, T &out, redOp op) const
      { allreduceRaw (&in, &out, 1, op); }

    void allreduceRawVoid (void *data, NDT type, tsize num, redOp op) const;
    template<typename T> void allreduceRaw (T *data, tsize num,
      redOp op) const
      { allreduceRawVoid (data, nativeType<T>(), num, op); }
    template<typename T> void allreduce (arr<T> &data, redOp op) const
      { allreduceRaw (&data[0], data.size(), op); }
    template<typename T> void allreduce (T &data, redOp op) const
      { allreduceRaw (&data, 1, op); }

    void bcastRawVoid (void *data, NDT type, tsize num, int root=0) const;
    template<typename T> void bcastRaw (T *data, tsize num, int root=0) const
      { bcastRawVoid (data, nativeType<T>(), num, root); }
    template<typename T> void bcast (arr<T> &data, int root=0) const
      { bcastRaw (&data[0], data.size(), root); }
    template<typename T> void bcast (T &data, int root=0) const
      { bcastRaw (&data, 1, root); }

    /*! NB: \a num refers to the <i>total</i> number of items in the arrays;
        the individual message size is \a num/num_ranks(). */
    void all2allRawVoid (const void *in, void *out, NDT type, tsize num) const;
    /*! NB: \a num refers to the <i>total</i> number of items in the arrays;
        the individual message size is \a num/num_ranks(). */
    template<typename T> void all2allRaw (const T *in, T *out, tsize num) const
      { all2allRawVoid (in, out, nativeType<T>(), num); }
    template<typename T> void all2all (const arr<T> &in, arr<T> &out) const
      {
      out.alloc(in.size());
      all2allRaw (&in[0], &out[0], in.size());
      }

    void all2allvRawVoid (const void *in, const int *numin, const int *disin,
      void *out, const int *numout, const int *disout, NDT type) const;
    template<typename T> void all2allvRaw (const T *in, const int *numin,
      const int *disin, T *out, const int *numout, const int *disout) const
      { all2allvRawVoid (in,numin,disin,out,numout,disout,nativeType<T>()); }
    /*!\deprecated */
    template<typename T> void all2allv (const arr<T> &in, const arr<int> &numin,
      const arr<int> &disin, arr<T> &out, const arr<int> &numout,
      const arr<int> &disout, tsize outsize) const
      {
      out.alloc(outsize);
      all2allvRaw (&in[0],&numin[0],&disin[0],&out[0],&numout[0],&disout[0]);
      }
    template<typename T> void all2allv_easy (const arr<T> &in,
      const arr<int> &numin, arr<T> &out, arr<int> &numout) const
      {
      arr<int> disin,disout;
      all2allv_easy_prep (in.size(),numin,disin,numout,disout);
      out.alloc(disout[num_ranks_-1]+numout[num_ranks_-1]);
      all2allvRawVoid (&in[0], &numin[0], &disin[0], &out[0], &numout[0],
        &disout[0], nativeType<T>());
      }
  };

//extern MPI_Manager mpiMgr;

} // namespace levels

#endif
