AM_CPPFLAGS = -I$(top_srcdir)/levels/c_utils

#include_HEADERS = bluestein.h fftpack.h ls_fft.h

noinst_LTLIBRARIES = libfftpack.la
libfftpack_la_SOURCES = bluestein.c bluestein.h fftpack.c fftpack.h ls_fft.c ls_fft.h

EXTRA_DIST = fftpack_inc.c

bin_PROGRAMS = ffttest
ffttest_SOURCES = ffttest.c
ffttest_LDADD = libfftpack.la $(top_builddir)/levels/c_utils/libc_utils.la $(LDADD)
