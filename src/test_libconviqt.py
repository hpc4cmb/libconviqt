import ctypes

# Make sure to check your LD_LIBRARY_PATH

libmpi = ctypes.CDLL( 'libmpi.so' )
libcfitsio = ctypes.CDLL( 'libcfitsio.so' )
#libconviqt = ctypes.CDLL( '../lib/libconviqt.so' )
libconviqt = ctypes.CDLL( 'libconviqt.so' )

