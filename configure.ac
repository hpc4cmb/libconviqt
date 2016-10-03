dnl
dnl This is the Autoconf file for libconviqt
dnl
dnl +------------------------
dnl | Initialize package info
dnl +------------------------
AC_INIT([libconviqt], [2.0], [https://github.com/hpc4cmb/libconviqt/issues])
AC_CONFIG_SRCDIR([Makefile.am])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_HEADERS(config.h)
AC_CONFIG_MACRO_DIR([m4])

dnl +-------------------------
dnl | If the library source code has changed at all since the last update, 
dnl | then increment revision (c:r:a becomes c:r+1:a).  If any interfaces 
dnl | have been added, removed, or changed since the last update, increment
dnl | current, and set revision to 0.  If any interfaces have been added 
dnl | since the last public release, then increment age.  If any interfaces 
dnl | have been removed since the last public release, then set age to 0.
dnl +-------------------------
LIBCONVIQT_LT_VERSION="2:0:0"
AC_SUBST(LIBCONVIQT_LT_VERSION)

dnl +---------------------
dnl | Compiler features
dnl +---------------------
dnl Set default optimization to '-O3' instead of the default '-g -O2'
if test x"${CFLAGS}" = x; then
   CFLAGS="-O3"
fi
if test x"${CXXFLAGS}" = x; then
   CXXFLAGS="-O3"
fi
AC_CANONICAL_HOST
AC_PROG_INSTALL

dnl +-------------------------
dnl | C++ support -- use MPI compiler throughout
dnl +-------------------------
AC_PROG_CXX
AX_PROG_CXX_MPI([],[],[
  AC_MSG_FAILURE([MPI C++ compiler is required but not found.])
])
AX_CXX_COMPILE_STDCXX([11])

dnl +-------------------------
dnl | C support -- use MPI compiler throughout
dnl +-------------------------
dnl AC_PROG_CC
AX_PROG_CC_MPI([],[],[
  AC_MSG_FAILURE([MPI CC compiler is required but not found.])
])

dnl libtool must only be called after the c/c++ compilers are overridden
dnl by MPI versions using AX_PROG_CC_MPI and AX_PROG_CXX_MPI

LT_INIT

dnl +------------------------------------------------
dnl | This is a C++ project...
dnl +------------------------------------------------
AC_LANG([C++])
AC_LANG_PUSH([C])
AX_OPENMP
AC_LANG_POP([C])

dnl +------------------------------------------------
dnl | ... that will one day be OpenMP parallel
dnl +------------------------------------------------
AX_OPENMP
AC_SUBST([CXXFLAGS],["$CXXFLAGS $OPENMP_CXXFLAGS"])
AC_SUBST([LDFLAGS],["$LDFLAGS $OPENMP_CXXFLAGS"])

AC_LANG_PUSH([C])
AX_OPENMP
AC_LANG_POP([C])
AC_SUBST([CFLAGS],["$CFLAGS $OPENMP_CFLAGS"])

dnl +------------------------------------------------
dnl | Check for CFITSIO
dnl +------------------------------------------------
AC_LANG_PUSH([C])
AX_OPENMP
ACX_CFITSIO([], [AC_MSG_ERROR([Could not find the CFITSIO library!])])
AC_LANG_POP([C])

dnl +------------------------------------------------
dnl | Set outputs
dnl +------------------------------------------------
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([levels/Makefile])
AC_CONFIG_FILES([levels/cxxsupport/Makefile])
AC_CONFIG_FILES([levels/c_utils/Makefile])
AC_CONFIG_FILES([levels/libfftpack/Makefile])
AC_CONFIG_FILES([levels/Healpix_cxx/Makefile])

dnl +------------------------------------------------
dnl | Generate Makefiles
dnl +------------------------------------------------
AC_OUTPUT

dnl +------------------------------------------------
dnl | Print out detected build options
dnl +------------------------------------------------
AC_MSG_NOTICE([                                           ])
AC_MSG_NOTICE([=========== Build Configuration ===========])
AC_MSG_NOTICE([                                           ])
AC_MSG_NOTICE([  C++ Compiler              : $CXX])
AC_MSG_NOTICE([  C++ Compile Flags         : $CXXFLAGS])
AC_MSG_NOTICE([  C++ Linking Flags         : $LDFLAGS])
AC_MSG_NOTICE([  C Compiler                : $CC])
AC_MSG_NOTICE([  C Compile Flags           : $CFLAGS])
AC_MSG_NOTICE([  CFITSIO Compile           : $CFITSIO_CPPFLAGS])
AC_MSG_NOTICE([  CFITSIO Link              : $CFITSIO])
AC_MSG_NOTICE([                                           ])
