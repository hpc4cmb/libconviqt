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

/*
 *  Functions for conversion of time formats
 *
 *  Copyright (C) 2009 Max-Planck-Society
 *  Author: Martin Reinecke
 */

/* get some POSIX functions from system headers */
#define _POSIX_C_SOURCE 200112L

#include <time.h>
#include <stdlib.h>
#include "c_utils.h"
#include "time_utils.h"

static time_t my_timegm(struct tm *tm)
  {
  time_t ret;
  char *tz;

  tz = getenv("TZ");
  setenv("TZ", "", 1);
  tzset();
  ret = mktime(tm);
  if (tz)
    setenv("TZ", tz, 1);
  else
    unsetenv("TZ");
  tzset();
  return ret;
  }

static time_t getunixzero(void)
  {
  struct tm unix_zero;
  time_t zero;
  unix_zero.tm_sec=0;
  unix_zero.tm_min=0;
  unix_zero.tm_hour=0;
  unix_zero.tm_mday=1;
  unix_zero.tm_mon=0;
  unix_zero.tm_year=70;
  unix_zero.tm_isdst=0;
  zero=my_timegm(&unix_zero);
  UTIL_ASSERT(zero!=-1,"error in mktime()");
  return zero;
  }

static time_t unixzero;
static int init=0;

static void initcheck()
  {
  if (init!=0) return;
  unixzero = getunixzero();
  init=1;
  }

#define SECS_1958_TO_1970 378691200

long long dateandtime2planck (int yr, int mth, int day, int hr, int min,
  int sec)
  {
  long long tplanck;
  struct tm mytime;
  initcheck();
  mytime.tm_sec = sec;
  mytime.tm_min = min;
  mytime.tm_hour = hr;
  mytime.tm_mday = day;
  mytime.tm_mon = mth-1;
  mytime.tm_year = yr-1900;
  mytime.tm_isdst=0;
  tplanck = my_timegm(&mytime);
  UTIL_ASSERT(tplanck!=-1,"error in mktime()");
  return tplanck-unixzero+SECS_1958_TO_1970;
  }

void planck2dateandtime (long long tplanck, int *yr, int *mth, int *day,
  int *hr, int *min, int *sec)
  {
  time_t mytimet;
  struct tm *mytime;
  initcheck();
  mytimet = tplanck+unixzero-SECS_1958_TO_1970;
  mytime = gmtime(&mytimet);
  UTIL_ASSERT(mytime!=NULL,"error in gmtime()");
  *sec = mytime->tm_sec;
  *min = mytime->tm_min;
  *hr = mytime->tm_hour;
  *day = mytime->tm_mday;
  *mth = mytime->tm_mon+1;
  *yr = mytime->tm_year+1900;
  }
