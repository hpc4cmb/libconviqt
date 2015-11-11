CC = @MPICC@
CXX = @MPICXX@

# Initialize compile and linking flags

AM_CPPFLAGS = \
-I$(top_srcdir)/levels/c_utils \
-I$(top_srcdir)/levels/cxxsupport

AM_CFLAGS =
AM_CXXFLAGS =
AM_LDFLAGS =
AM_LIBS = 

# Append CFITSIO to compile and link flags
AM_CPPFLAGS += $(CFITSIO_CPPFLAGS)
AM_LIBS += $(CFITSIO)

# Always need libm
AM_LIBS += -lm

# Level-S source expects this macro for use with MPI:
AM_CPPFLAGS += -DUSE_MPI


noinst_HEADERS = \
    ls_alloc_utils.h \
    ls_arr.h \
    ls_datatypes.h \
    ls_error_handling.h \
    ls_math_utils.h \
    ls_mpi_support.h \
    ls_share_utils.h \
    ls_sse_utils_cxx.h \
    ls_wigner.h \
    ls_fftpack_support.h \
    ls_constants.h \
    ls_xcomplex.h \
    ls_fitshandle.h \
    ls_rng.h \
    ls_string_utils.h \
    ls_rotmatrix.h \
    ls_openmp_support.h \
    ls_safe_cast.h \
    ls_vec3.h

noinst_LTLIBRARIES = libls_cxxsupport.la

libls_cxxsupport_la_SOURCES = \
    error_handling.cc \
    font_data.inc \
    mpi_support.cc \
    wigner.cc \
    fitshandle.cc \
    string_utils.cc \
    rotmatrix.cc
