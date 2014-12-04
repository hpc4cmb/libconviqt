AM_CPPFLAGS = $(CFITSIO_CPPFLAGS) -I$(top_srcdir)/levels/c_utils -DUSE_MPI
LDADD = $(CFITSIO)

noinst_HEADERS = ls_alloc_utils.h ls_arr.h ls_datatypes.h ls_error_handling.h ls_math_utils.h ls_mpi_support.h ls_share_utils.h ls_sse_utils_cxx.h ls_wigner.h ls_fftpack_support.h ls_constants.h

noinst_LTLIBRARIES = libls_cxxsupport.la

libls_cxxsupport_la_SOURCES = ls_error_handling.h error_handling.cc font_data.inc ls_alloc_utils.h ls_arr.h ls_datatypes.h ls_math_utils.h ls_mpi_support.h mpi_support.cc ls_share_utils.h ls_sse_utils_cxx.h ls_wigner.h wigner.cc ls_fftpack_support.h ls_constants.h
