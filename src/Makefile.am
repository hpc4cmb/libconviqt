AM_CPPFLAGS = $(HEALPIX_CXX_CPPFLAGS)
#AM_CXXFLAGS = $(OPENMP_CXXFLAGS)
#LDADD = $(HEALPIX_CXX) $(OPENMP_CXXFLAGS)
LDADD = $(HEALPIX_CXX) $(OPENMP_CXXFLAGS)

include_HEADERS = conviqt.hpp

lib_LTLIBRARIES = libconviqt.la
libconviqt_la_SOURCES = conviqt_math.cpp conviqt_io.cpp conviqt.hpp

bin_PROGRAMS = test_libconviqt
test_libconviqt_SOURCES = test_libconviqt.cpp
test_libconviqt_LDADD = libconviqt.la $(LDADD)
