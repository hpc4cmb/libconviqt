# Initialize compile and linking flags

AM_CPPFLAGS = \
-I$(top_srcdir)/levels/c_utils \
-I$(top_srcdir)/levels/libfftpack

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
    bluestein.h \
    fftpack.h \
    ls_fft.h

EXTRA_DIST = fftpack_inc.c

noinst_LTLIBRARIES = libfftpack.la

libfftpack_la_SOURCES = \
    bluestein.c \
    fftpack.c \
    ls_fft.c

noinst_PROGRAMS = ffttest

ffttest_SOURCES = ffttest.c
ffttest_LDADD = \
    libfftpack.la \
    $(top_builddir)/levels/c_utils/libc_utils.la \
    $(AM_LIBS)
