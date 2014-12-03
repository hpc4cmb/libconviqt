#
# SYNOPSIS
#
#   ACX_HEALPIX_CXX([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro looks for a version of the Healpix C++ library.  The HEALPIX_CXX_CPPFLAGS
#   and HEALPIX_CXX output variables hold the compile and link flags.
#
#   To link an application with Healpix, you should link with:
#
#   	$HEALPIX_CXX
#
#   The user may use:
# 
#       --with-healpix_cxx=<PATH>
#
#   to manually specify the Healpix C++ installation prefix.
#
#   ACTION-IF-FOUND is a list of shell commands to run if a Healpix C++ library is
#   found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it is
#   not found. If ACTION-IF-FOUND is not specified, the default action will
#   define HAVE_HEALPIX_CXX and set output variables above.
#
#   This macro requires autoconf 2.50 or later.
#
# LAST MODIFICATION
#
#   2013-09-30
#
# COPYING
#
#   Copyright (c) 2013 Reijo Keskitalo <rtkeskitalo@lbl.gov>
#
#   All rights reserved.
#
#   Redistribution and use in source and binary forms, with or without modification,
#   are permitted provided that the following conditions are met:
#
#   o  Redistributions of source code must retain the above copyright notice, 
#      this list of conditions and the following disclaimer.
#
#   o  Redistributions in binary form must reproduce the above copyright notice, 
#      this list of conditions and the following disclaimer in the documentation
#      and/or other materials provided with the distribution.
#
#   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
#   ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
#   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
#   IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
#   INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
#   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
#   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
#   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
#   OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
#   OF THE POSSIBILITY OF SUCH DAMAGE.
#

AC_DEFUN([ACX_HEALPIX_CXX], [
AC_PREREQ(2.50)
#AC_REQUIRE([ACX_MPI])
AC_REQUIRE([ACX_CFITSIO])
#AC_REQUIRE([AX_OPENMP])

acx_healpix_cxx_ok=no
acx_healpix_cxx_default="-lhealpix_cxx -lcxxsupport -lfftpack"
#acx_healpix_cxx_default="-lhealpix_cxx"

HEALPIX_CXX_CPPFLAGS=""
HEALPIX_CXX=""

AC_ARG_WITH(healpix_cxx, [AC_HELP_STRING([--with-healpix_cxx=<PATH>], [use the Healpix C++ installed in <PATH>.])])

if test x"$with_healpix_cxx" != x; then
   if test x"$with_healpix_cxx" != xno; then
      if test -d "$with_healpix_cxx/include/healpix_cxx"; then
        HEALPIX_CXX_CPPFLAGS="-I$with_healpix_cxx/include/healpix_cxx"
      else
        HEALPIX_CXX_CPPFLAGS="-I$with_healpix_cxx/include"
      fi
      HEALPIX_CXX="-L$with_healpix_cxx/lib $acx_healpix_cxx_default"
   else
      acx_healpix_cxx_ok=disable
   fi
fi

if test $acx_healpix_cxx_ok = disable; then
   echo "**** Healpix C++ explicitly disabled by configure."
else

   # Save environment

   acx_healpix_cxx_save_CXX="$CXX"
   acx_healpix_cxx_save_CXXCPP="$CXXCPP"
   acx_healpix_cxx_save_CPPFLAGS="$CPPFLAGS"
   acx_healpix_cxx_save_LIBS="$LIBS"

   # Test compile and linking.

   #CXX="$MPICXX"
   #CXXCPP="$MPICXX -E"
   CPPFLAGS="$CPPFLAGS $HEALPIX_CXX_CPPFLAGS"
   LIBS="$HEALPIX_CXX $acx_healpix_cxx_save_LIBS $LIBS -lm $OPENMP_CXXFLAGS $CFITSIO"

   AC_MSG_CHECKING([for Healpix_Map functionality])
   AC_LINK_IFELSE([AC_LANG_PROGRAM([
     [#include <healpix_map.h>]
   ],[
      int nside=8;
      Healpix_Ordering_Scheme scheme=NEST;
      Healpix_Map<float> map;
      map.SetNside( nside, scheme );
      double d = map.max_pixrad();
   ])],[acx_healpix_cxx_ok=yes;AC_DEFINE(HAVE_HEALPIX_CXX,1,[Define if you have the Healpix C++ library.])])

   AC_MSG_RESULT($acx_healpix_cxx_ok)

   if test $acx_healpix_cxx_ok = no; then
      HEALPIX_CXX="$acx_healpix_cxx_default"
      LIBS="$HEALPIX_CXX $acx_healpix_cxx_save_LIBS -lm $OPENMP_CXXFLAGS $CFITSIO"

      AC_MSG_CHECKING([for Healpix C++ in default location])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([
        [#include <healpix_map.h>]
      ],[
         int nside=8;
         Healpix_Ordering_Scheme scheme=NEST;
         Healpix_Map<float> map;
         map.SetNside( nside, scheme );
         double d = map.max_pixrad();
      ])],[acx_healpix_cxx_ok=yes;AC_DEFINE(HAVE_HEALPIX_CXX,1,[Define if you have the Healpix C++ library.])])

      AC_MSG_RESULT($acx_healpix_cxx_ok)
   fi

   if test $acx_healpix_cxx_ok = no; then
      HEALPIX_CXX=""
   fi

   # Restore environment

   CXX="$acx_healpix_cxx_save_CXX"
   CXXCPP="$acx_healpix_cxx_save_CXXCPP"
   CPPFLAGS="$acx_healpix_Cxx_save_CPPFLAGS"
   LIBS="$acx_healpix_cxx_save_LIBS"

fi

# Define exported variables

AC_SUBST(HEALPIX_CXX_CPPFLAGS)
AC_SUBST(HEALPIX_CXX)

# Execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
   
if test x"$acx_healpix_cxx_ok" = xyes; then
   ifelse([$1],,[echo "**** Enabling support for Healpix C++."],[$1])
else
   ifelse([$2],,[echo "**** Healpix C++ not found - disabling support."],[$2])
fi

])
