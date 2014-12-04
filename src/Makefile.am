AM_CPPFLAGS = -I$(top_srcdir)/levels/c_utils -I$(top_srcdir)/levels/cxxsupport -I$(top_srcdir)/levels/libfftpack $(CFITSIO_CPPFLAGS) $(HEALPIX_CXX_CPPFLAGS)
LDADD = $(HEALPIX_CXX) $(CFITSIO)

include_HEADERS = conviqt.hpp

lib_LTLIBRARIES = libconviqt.la
libconviqt_la_LIBADD = $(top_srcdir)/levels/c_utils/libc_utils.la $(top_srcdir)/levels/cxxsupport/libls_cxxsupport.la $(top_srcdir)/levels/libfftpack/libfftpack.la
libconviqt_la_SOURCES = conviqt_math.cpp conviqt_io.cpp conviqt_util.cpp wignercalc.cpp conviqt.hpp

noinst_PROGRAMS = test_libconviqt
test_libconviqt_SOURCES = test_libconviqt.cpp
test_libconviqt_LDADD = libconviqt.la $(LDADD)

TESTS = test_libconviqt
