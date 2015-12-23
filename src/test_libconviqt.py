#!/usr/bin/env python

from mpi4py import MPI
import ctypes
import numpy as np
import sys

# Make sure to check your LD_LIBRARY_PATH

comm = MPI.COMM_WORLD
itask = comm.Get_rank()
ntask = comm.Get_size()

#comm_ptr = MPI._addressof( comm )
#fcomm = comm.py2f()

if itask == 0: print 'Running with ', ntask, ' MPI tasks'

#libmpi = ctypes.CDLL( 'libmpi.so' )
#libmpi_cxx = ctypes.CDLL( 'libmpi_cxx.so' )
#libcfitsio = ctypes.CDLL( 'libcfitsio.so' )
try:
    libconviqt = ctypes.CDLL( '.libs/libconviqt.so' )
except:
    try:
        libconviqt = ctypes.CDLL( 'libconviqt.so' )
    except Exception, e:
        raise Exception('Unable to load libconviqt.so. You may need to install it first with "make install": {}'.format(e))

lmax = 32
beamlmax = lmax
beammmax = 4
pol = 0
fwhm = 4.
beamfile = '../data/mb_lfi_30_27_x_rescaled.alm'
skyfile = '../data/slm.fits'
det_id = 'LFITEST'

ntheta = 3
nphi = 3
npsi = 3
nsamp = ntheta*nphi*npsi

nbetafac = 10
mcsamples = 0
lmaxout = 32
order = 3


beam = libconviqt.conviqt_beam_new()
err = libconviqt.conviqt_beam_read( beam, ctypes.c_long(beamlmax), ctypes.c_long(beammmax), ctypes.c_byte(pol), beamfile, comm.py2f()  )
if err != 0: raise Exception( 'Failed to load ' + beamfile )

sky = libconviqt.conviqt_sky_new()
err = libconviqt.conviqt_sky_read( sky, ctypes.c_long(lmax), ctypes.c_byte(pol), skyfile, ctypes.c_double(fwhm), comm.py2f()  )
if err != 0: raise Exception( 'Failed to load ' + skyfile )

detector = libconviqt.conviqt_detector_new_with_id( det_id )

s = ctypes.create_string_buffer( 255 )
libconviqt.conviqt_detector_get_id( detector, s )
print 'detector name = ' + s.value
epsilon = 1.32495160e-04
libconviqt.conviqt_detector_set_epsilon( detector, ctypes.c_double(epsilon) )
eps = ctypes.c_double( 0. )
libconviqt.conviqt_detector_get_epsilon( detector, ctypes.byref(eps) )
print 'epsilon = ',eps.value

print 'Allocating pointing array'

pnt = libconviqt.conviqt_pointing_new()
err = libconviqt.conviqt_pointing_alloc( pnt, ctypes.c_long(nsamp*5) )
if err != 0: raise Exception( 'Failed to allocate pointing array' )
libconviqt.conviqt_pointing_data.restype = ctypes.POINTER( ctypes.c_double )
ppnt = libconviqt.conviqt_pointing_data( pnt )

print 'Populating pointing array'

row = 0
for theta in np.arange(ntheta) * np.pi / ntheta:
    for phi in np.arange(nphi) * 2 * np.pi / nphi:
        for psi in np.arange(npsi) * np.pi / npsi:
            ppnt[ row*5 + 0 ] = phi
            ppnt[ row*5 + 1 ] = theta
            ppnt[ row*5 + 2 ] = psi
            ppnt[ row*5 + 3 ] = 0
            ppnt[ row*5 + 4 ] = row
            row += 1

for i in range(10):
    print 'ppnt[',i,'] = ',ppnt[i]

print 'Creating convolver'

convolver = libconviqt.conviqt_convolver_new( sky, beam, detector, ctypes.c_byte(pol), ctypes.c_long(lmax), ctypes.c_long(beammmax), ctypes.c_long(nbetafac), ctypes.c_long(mcsamples), ctypes.c_long(lmaxout), ctypes.c_long(order), comm.py2f() )
if convolver == 0: raise Exception( "Failed to instantiate convolver" );

print 'Convolving data'

err = libconviqt.conviqt_convolver_convolve( convolver, pnt )
if err != 0: raise Exception( 'Convolution FAILED!' )

# The pointer to the data will have changed during the convolution call ...

ppnt = libconviqt.conviqt_pointing_data( pnt )

if itask == 0:

    print 'Convolved TOD: '
    for row in np.arange( nsamp ):
        print ppnt[row*5+4], ' ', ppnt[row*5+3]

    if pol:
        if np.abs( ppnt[ 0*5+3] -  0.8546349819096275 ) > 1e-6: raise Exception( "Row 0 should be 0.8546349819096275, not " + str(ppnt[ 0*5+3]) )
        if np.abs( ppnt[10*5+3] +  25.53734467183137 )  > 1e-6: raise Exception( "Row 10 should be -25.53734467183137, not " + str(ppnt[10*5+3]) )
        if np.abs( ppnt[15*5+3] +  76.04945574990082 )  > 1e-6: raise Exception( "Row 15 should be -76.04945574990082, not " + str(ppnt[15*5+3]) )
    else:
        if np.abs( ppnt[ 0*5+3] -  0.8545846415739397 ) > 1e-6: raise Exception( "Row 0 should be 0.8545846415739397, not " + str(ppnt[ 0*5+3]) )
        if np.abs( ppnt[10*5+3] +  25.20150061107036 )  > 1e-6: raise Exception( "Row 10 should be -25.20150061107036, not " + str(ppnt[ 10*5+3]) )
        if np.abs( ppnt[15*5+3] +  76.14723911261254 )  > 1e-6: raise Exception( "Row 15 should be -76.14723911261254, not " + str(ppnt[ 15*5+3]) )

libconviqt.conviqt_convolver_del( convolver )
libconviqt.conviqt_pointing_del( pnt )
libconviqt.conviqt_detector_del( detector )
libconviqt.conviqt_sky_del( sky )
libconviqt.conviqt_beam_del( beam )
