AM_CPPFLAGS = $(CFITSIO_CPPFLAGS) -I$(top_srcdir)/levels/c_utils -DUSE_MPI
LDADD = $(CFITSIO)

#include_HEADERS = alloc_utils.h bstream.h ephemerides.h fitshandle.h io_utils.h iohandle_dmc.h levels_facilities.h ls_math_utils.h openmp_support.h pointing.h rangeset.h safe_ptr.h se_utils_cxx.h ls_vec3.h xcomplex.h announce.h ls_cxxutils.h error_handling.h focalplane_db.h iohandle.h iohandle_fits.h ls_image.h mmap_handle.h paramfile.h psort.h rotmatrix.h share_utils.h string_utils.h walltimer.h ls_arr.h datatypes.h ls_fftpack_support.h geom_utils.h iohandle_current.h iohandle_helpers.h lsconstants.h ls_mpi_support.h planck_rng.h quaternion.h safe_cast.h sort_utils.h trafos.h ls_wigner.h

noinst_HEADERS = ls_alloc_utils.h ls_arr.h ls_datatypes.h ls_error_handling.h ls_math_utils.h ls_mpi_support.h ls_share_utils.h ls_sse_utils_cxx.h ls_wigner.h ls_fftpack_support.h

noinst_LTLIBRARIES = libls_cxxsupport.la

#libcxxsupport_la_SOURCES = alloc_utils.h ls_cxxutils.h error_handling.cc focalplane_db.cc io_utils.cc iohandle_current.h iohandle_helpers.h ls_math_utils.h openmp_support.h pointing.h rotmatrix.cc sort_utils.h trafos.h wigner.h datatypes.h error_handling.h focalplane_db.h io_utils.h levels_facilities.h mmap_handle.cc paramfile.cc psort.h rotmatrix.h sse_utils_cxx.h ls_vec3.h xcomplex.h ls_fftpack_support.h font_data.inc iohandle.cc ls_image.cc mmap_handle.h paramfile.h quaternion.cc safe_cast.h string_utils.cc walltimer.cc ls_arr.h ephemerides.cc fitshandle.cc geom_utils.cc iohandle.h iohandle_fits.cc ls_image.h mpi_support.cc planck_rng.h quaternion.h safe_ptr.h string_utils.h walltimer.h bstream.h ephemerides.h fitshandle.h geom_utils.h iohandle_current.cc iohandle_fits.h lsconstants.h ls_mpi_support.h pointing.cc rangeset.h share_utils.h trafos.cc wigner.cc
libls_cxxsupport_la_SOURCES = ls_error_handling.h error_handling.cc font_data.inc ls_alloc_utils.h ls_arr.h ls_datatypes.h ls_math_utils.h ls_mpi_support.h mpi_support.cc ls_share_utils.h ls_sse_utils_cxx.h ls_wigner.h wigner.cc ls_fftpack_support.h

#EXTRA_DIST = font_data.inc

