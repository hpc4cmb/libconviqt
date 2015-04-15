AM_CPPFLAGS = $(CFITSIO_CPPFLAGS) -I$(top_srcdir)/levels/c_utils -I$(top_srcdir)/levels/cxxsupport
LDADD = $(CFITSIO)

noinst_HEADERS = ls_alm.h ls_alm_fitsio.h ls_alm_powspec_tools.h ls_powspec.h

noinst_LTLIBRARIES = libhealpix_cxx.la

libhealpix_cxx_la_SOURCES = ls_alm.h alm.cc ls_alm_fitsio.h alm_fitsio.cc ls_alm_powspec_tools.h alm_powspec_tools.cc ls_powspec.h powspec.cc
