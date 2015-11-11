# Initialize compile and linking flags

AM_CPPFLAGS = \
-I$(top_srcdir)/levels/c_utils \
-I$(top_srcdir)/levels/libfftpack \
-I$(top_srcdir)/levels/cxxsupport \
-I$(top_srcdir)/levels/Healpix_cxx

AM_CFLAGS =
AM_CXXFLAGS =
AM_LDFLAGS =
AM_LIBS = 

# Append CFITSIO to compile and link flags
AM_CPPFLAGS += $(CFITSIO_CPPFLAGS)
AM_LIBS += $(CFITSIO)

# Always need libm
AM_LIBS += -lm


noinst_HEADERS = \
    ls_alm.h \
    ls_alm_fitsio.h \
    ls_alm_powspec_tools.h \
    ls_powspec.h

noinst_LTLIBRARIES = libhealpix_cxx.la

libhealpix_cxx_la_SOURCES = \
    alm.cc \
    alm_fitsio.cc \
    alm_powspec_tools.cc \
    powspec.cc
