
CC = @MPICC@
CXX = @MPICXX@

# Initialize compile and linking flags

AM_CPPFLAGS = \
-I$(top_srcdir)/levels/c_utils \
-I$(top_srcdir)/levels/libfftpack \
-I$(top_srcdir)/levels/cxxsupport \
-I$(top_srcdir)/levels/Healpix_cxx \
-I$(top_srcdir)/src

AM_CFLAGS =
AM_CXXFLAGS =
AM_LDFLAGS =
AM_LIBS = 

# Append CFITSIO to compile and link flags
AM_CPPFLAGS += $(CFITSIO_CPPFLAGS)
AM_LIBS += $(CFITSIO)

# Always need libm
AM_LIBS += -lm


include_HEADERS = \
    conviqt.hpp \
    cconviqt.h

lib_LTLIBRARIES = libconviqt.la

libconviqt_la_SOURCES = \
    conviqt_math.cpp \
    conviqt_io.cpp \
    conviqt_util.cpp \
    wignercalc.cpp \
    cconviqt.cpp

libconviqt_la_LIBADD = \
    $(top_builddir)/levels/Healpix_cxx/libhealpix_cxx.la \
    $(top_builddir)/levels/cxxsupport/libls_cxxsupport.la \
    $(top_builddir)/levels/libfftpack/libfftpack.la \
    $(top_builddir)/levels/c_utils/libc_utils.la


noinst_PROGRAMS = test_libconviqt test_libcconviqt

test_libconviqt_SOURCES = test_libconviqt.cpp
test_libconviqt_LDADD = libconviqt.la $(AM_LIBS)

test_libcconviqt_SOURCES = test_libcconviqt.c
test_libcconviqt_LDADD = libconviqt.la $(AM_LIBS)

TESTS = test_libconviqt test_libcconviqt test_libconviqt.py
