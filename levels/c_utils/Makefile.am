# Initialize compile and linking flags

AM_CPPFLAGS = \
-I$(top_srcdir)/levels/c_utils

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
    c_utils.h \
    memusage.h \
    sse_utils.h \
    svd_c.h \
    time_utils.h \
    walltime_c.h \
    wrec3a.h

noinst_LTLIBRARIES = libc_utils.la

libc_utils_la_SOURCES = \
    c_utils.c \
    memusage.c \
    svd_c.c \
    time_utils.c \
    walltime_c.c \
    wrec3a.c
