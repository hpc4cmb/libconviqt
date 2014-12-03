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
 * 2014-12-01 - Reijo Keskitalo 
 */

/*
 *  libcxxsupport is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Utilities for error reporting
 *
 *  Copyright (C) 2003-2011 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "ls_error_handling.h"
#include "ls_mpi_support.h"

using namespace std;

namespace levels {

LEVELSError::LEVELSError(const string &message) : msg (message) {}
LEVELSError::LEVELSError(const char *message) : msg (message) {}

//virtual
LEVELSError::~LEVELSError() {}

void levels_failure__(const char *file, int line, const char *func,
  const string &msg)
  {
  cerr << "Error encountered at " << file << ", line " << line << endl;
  if (func) cerr << "(function " << func << ")" << endl;
  if (msg!="") cerr << endl << msg << endl;
  cerr << endl;
  }

void levels_failure__(const char *file, int line, const char *func,
  const char *msg)
  { levels_failure__ (file,line,func,string(msg)); }

void killjob__()
  {
    MPI_Manager mpiMgr;
    mpiMgr.abort();
  }

} // namespace levels
