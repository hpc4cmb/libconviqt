#
# SYNOPSIS
#
#   ACX_CFITSIO([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro looks for a version of the CFITSIO library.  The CFITSIO_CPPFLAGS
#   and CFITSIO output variables hold the compile and link flags.
#
#   To link an application with CFITSIO, you should link with:
#
#   	$CFITSIO
#
#   The user may use:
# 
#       --with-cfitsio=<path> 
#
#   to manually specify the path to the cfitsio installation.
#
#   ACTION-IF-FOUND is a list of shell commands to run if an CFITSIO library is
#   found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it is
#   not found. If ACTION-IF-FOUND is not specified, the default action will
#   define HAVE_CFITSIO and set output variables above.
#
#   This macro requires autoconf 2.50 or later.
#
# LAST MODIFICATION
#
#   2013-03-09
#
# COPYING
#
#   Copyright (c) 2013 Theodore Kisner <tskisner@lbl.gov>
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


AC_DEFUN([ACX_CFITSIO], [
AC_PREREQ(2.50)

acx_cfitsio_ok=no
acx_cfitsio_default="-lcfitsio"

CFITSIO_CPPFLAGS=""
CFITSIO=""

AC_ARG_WITH(cfitsio, [AC_HELP_STRING([--with-cfitsio=<PATH>], [use CFITSIO installed in <path>.])])

if test x"$with_cfitsio" != x; then
   if test x"$with_cfitsio" != xno; then
      CFITSIO_CPPFLAGS="-I$with_cfitsio/include"
      CFITSIO="-L$with_cfitsio/lib -lcfitsio"
   else
      acx_cfitsio_ok=disable
   fi
fi

if test $acx_cfitsio_ok = disable; then
   echo "**** CFITSIO explicitly disabled by configure."
else

   # Save environment

   acx_cfitsio_save_CC="$CC"
   acx_cfitsio_save_CPP="$CPP"
   acx_cfitsio_save_CPPFLAGS="$CPPFLAGS"
   acx_cfitsio_save_LIBS="$LIBS"

   # Test serial compile and linking

   CPPFLAGS="$CPPFLAGS $CFITSIO_CPPFLAGS"
   LIBS="$CFITSIO $acx_cfitsio_save_LIBS -lm"

   AC_CHECK_HEADERS([fitsio.h])

   AC_MSG_CHECKING([for ffopen in user specified location])
   AC_TRY_LINK_FUNC(ffopen, [acx_cfitsio_ok=yes;AC_DEFINE(HAVE_CFITSIO,1,[Define if you have the CFITSIO library.])], [])
   AC_MSG_RESULT($acx_cfitsio_ok)

   if test $acx_cfitsio_ok = no; then
      CFITSIO="$acx_cfitsio_default"
      LIBS="$acx_cfitsio_default $acx_cfitsio_save_LIBS -lm"
      AC_MSG_CHECKING([for ffopen in default location])
      AC_TRY_LINK_FUNC(ffopen, [acx_cfitsio_ok=yes;AC_DEFINE(HAVE_CFITSIO,1,[Define if you have the CFITSIO library.])], [])
      AC_MSG_RESULT($acx_cfitsio_ok)
   fi

   if test $acx_cfitsio_ok = no; then
      CFITSIO=""
   fi

   # Restore environment

   CC="$acx_cfitsio_save_CC"
   CPP="$acx_cfitsio_save_CPP"
   LIBS="$acx_cfitsio_save_LIBS"
   CPPFLAGS="$acx_cfitsio_save_CPPFLAGS"

fi

# Define exported variables

AC_SUBST(CFITSIO_CPPFLAGS)
AC_SUBST(CFITSIO)

# Execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
   
if test x"$acx_cfitsio_ok" = xyes; then
   ifelse([$1],,[echo "**** Enabling support for CFITSIO."],[$1])
else
   ifelse([$2],,[echo "**** CFITSIO not found - disabling support."],[$2])
fi

])
