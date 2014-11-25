AC_INIT([libconviqt], [0.1], [bug-report@address])

AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([foreign -Wall -Werror])

# Ensure we are in the right directory
AC_CONFIG_SRCDIR([src/conviqt.hpp])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile levels/Makefile levels/cxxsupport/Makefile levels/c_utils/Makefile])

# Ensure the necessary compilers work. Libconviqt is an MPI C++ code
AC_PROG_CC_C99
AC_LANG_PUSH([C])
CPPFLAGS="$CPPFLAGS $OPENMP_CFLAGS"
CFLAGS="$CFLAGS $OPENMP_CFLAGS"

AC_PROG_CXX
AC_LANG_PUSH([C++])
AX_CXX_COMPILE_STDCXX_11
AC_OPENMP
CXXFLAGS="$CXXFLAGS $OPENMP_CFLAGS"
#AX_OPENMP
#AC_SUBST(OPENMP_CXXFLAGS)

ACX_MPI
AM_CONDITIONAL([HAVE_AM_MPI], [test x"$MPICXX" != x])
if test x"$MPICXX" != x; then
  CXX=$MPICXX
fi

AC_REQUIRE_CPP

# We need cfitsio and Healpix
ACX_CFITSIO([], [AC_MSG_ERROR([Could not find the CFITSIO library!])])
ACX_HEALPIX_CXX([], [AC_MSG_ERROR([could not find the Healpix C++ library])])

# We use libtool to build the libraries
AM_PROG_AR
LT_INIT

# Let the user know the results of ./configure. These lines are copied from TOAST/configure.ac
AC_MSG_NOTICE([                                           ])
AC_MSG_NOTICE([=========== Build Configuration ===========])
AC_MSG_NOTICE([                                           ])
AC_MSG_NOTICE([  C Compiler                : $CC])
AC_MSG_NOTICE([  C Compile Flags           : $CFLAGS])
AC_MSG_NOTICE([  C++ Compiler              : $CXX])
AC_MSG_NOTICE([  C++ Compile Flags         : $CXXFLAGS])
#AC_MSG_NOTICE([  C++ OpenMP Flags          : $OPENMP_CXXFLAGS])
#if test x"${MPICXX}" = x; then
#AC_MSG_NOTICE([  MPI C++ Compiler          : Not found])
#else
#AC_MSG_NOTICE([  MPI C++ Compiler          : $MPICXX])
#fi
AC_MSG_NOTICE([  CFITSIO Compile           : $CFITSIO_CPPFLAGS])
AC_MSG_NOTICE([  CFITSIO Link              : $CFITSIO])
AC_MSG_NOTICE([  Healpix C++ Compile       : $HEALPIX_CXX_CPPFLAGS])
AC_MSG_NOTICE([  Healpix C++ Link          : $HEALPIX_CXX])

AC_OUTPUT

