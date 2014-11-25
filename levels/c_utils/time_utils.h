/*
 *  This file is part of libc_utils.
 *
 *  libc_utils is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libc_utils is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libc_utils; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libc_utils is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*! \file time_utils.h
 *  Functions for conversion of time formats
 *
 *  Copyright (C) 2009 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_TIME_UTILS_H
#define PLANCK_TIME_UTILS_H

#ifdef __cplusplus
extern "C" {
#endif

/*! Converts a given date and time into seconds since Jan 1, 1958 00:00 UTC. */
long long dateandtime2planck (int yr, int mth, int day, int hr, int min,
  int sec);
/*! Converts seconds since Jan 1, 1958 00:00 UTC into date and time. */
void planck2dateandtime (long long tplanck, int *yr, int *mth, int *day,
  int *hr, int *min, int *sec);

#ifdef __cplusplus
}
#endif

#endif
