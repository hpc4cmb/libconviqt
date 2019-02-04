#!/usr/bin/env python

from __future__ import print_function
from mpi4py import MPI
import ctypes as ct
import numpy as np
import os
import sys

_libdir = os.path.dirname(__file__)
if not os.path.isdir(os.path.join(_libdir, ".libs")):
    _libdir = "."
if not os.path.isdir(os.path.join(_libdir, ".libs")):
    raise Exception("Cannot find the libtool .libs directory")

path = os.path.join(_libdir, ".libs", "libconviqt.so")
try:
    print("Trying to import libconviqt from", path, end="")
    sys.stdout.flush()
    libconviqt = ct.CDLL(path)
except:
    try:
        path2 = path.replace(".so", ".dylib")
        print(" FAILED.\nTrying to import libconviqt from", path2, end="")
        sys.stdout.flush()
        libconviqt = ct.CDLL(path2)
    except:
        print(" FAILED")
        sys.stdout.flush()
        raise
print(" SUCCESS")
sys.stdout.flush()

try:
    if MPI._sizeof(MPI.Comm) == ct.sizeof(ct.c_int):
        MPI_Comm = ct.c_int
    else:
        MPI_Comm = ct.c_void_p
except Exception as e:
    raise Exception(
        "Failed to set the portable MPI communicator datatype. MPI4py is probably too old. You may need to install from a git checkout. ({})".format(
            e
        )
    )

comm = MPI.COMM_WORLD
itask = comm.Get_rank()
ntask = comm.Get_size()

comm_ptr = MPI._addressof(comm)
comm = MPI_Comm.from_address(comm_ptr)

if itask == 0:
    print("Running with ", ntask, " MPI tasks")

# Beam functions

libconviqt.conviqt_beam_new.restype = ct.c_void_p
libconviqt.conviqt_beam_new.argtypes = []

libconviqt.conviqt_beam_del.restype = ct.c_int
libconviqt.conviqt_beam_del.argtypes = [ct.c_void_p]

libconviqt.conviqt_beam_read.restype = ct.c_int
libconviqt.conviqt_beam_read.argtypes = [
    ct.c_void_p,
    ct.c_long,
    ct.c_long,
    ct.c_byte,
    ct.c_char_p,
    MPI_Comm,
]

libconviqt.conviqt_beam_lmax.restype = ct.c_int
libconviqt.conviqt_beam_lmax.argtypes = [ct.c_void_p]

libconviqt.conviqt_beam_mmax.restype = ct.c_int
libconviqt.conviqt_beam_mmax.argtypes = [ct.c_void_p]

libconviqt.conviqt_beam_normalize.restype = ct.c_double
libconviqt.conviqt_beam_normalize.argtypes = [ct.c_void_p]

# Sky functions

libconviqt.conviqt_sky_new.restype = ct.c_void_p
libconviqt.conviqt_sky_new.argtypes = []

libconviqt.conviqt_sky_del.restype = ct.c_int
libconviqt.conviqt_sky_del.argtypes = [ct.c_void_p]

libconviqt.conviqt_sky_read.restype = ct.c_int
libconviqt.conviqt_sky_read.argtypes = [
    ct.c_void_p,
    ct.c_long,
    ct.c_byte,
    ct.c_char_p,
    ct.c_double,
    MPI_Comm,
]

libconviqt.conviqt_sky_lmax.restype = ct.c_int
libconviqt.conviqt_sky_lmax.argtypes = [ct.c_void_p]

libconviqt.conviqt_sky_remove_monopole.restype = ct.c_int
libconviqt.conviqt_sky_remove_monopole.argtypes = [ct.c_void_p]

libconviqt.conviqt_sky_remove_dipole.restype = ct.c_int
libconviqt.conviqt_sky_remove_dipole.argtypes = [ct.c_void_p]

# Detector functions

libconviqt.conviqt_detector_new.restype = ct.c_void_p
libconviqt.conviqt_detector_new.argtypes = []

libconviqt.conviqt_detector_new_with_id.restype = ct.c_void_p
libconviqt.conviqt_detector_new_with_id.argtypes = [ct.c_char_p]

libconviqt.conviqt_detector_del.restype = ct.c_int
libconviqt.conviqt_detector_del.argtypes = [ct.c_void_p]

libconviqt.conviqt_detector_set_epsilon.restype = ct.c_int
libconviqt.conviqt_detector_set_epsilon.argtypes = [ct.c_void_p, ct.c_double]

libconviqt.conviqt_detector_get_epsilon.restype = ct.c_int
libconviqt.conviqt_detector_get_epsilon.argtypes = [
    ct.c_void_p,
    ct.POINTER(ct.c_double),
]

libconviqt.conviqt_detector_get_id.restype = ct.c_int
libconviqt.conviqt_detector_get_id.argtypes = [ct.c_void_p, ct.c_char_p]

# Pointing functions

libconviqt.conviqt_pointing_new.restype = ct.c_void_p
libconviqt.conviqt_pointing_new.argtypes = []

libconviqt.conviqt_pointing_del.restype = ct.c_int
libconviqt.conviqt_pointing_del.argtypes = [ct.c_void_p]

libconviqt.conviqt_pointing_alloc.restype = ct.c_int
libconviqt.conviqt_pointing_alloc.argtypes = [ct.c_void_p, ct.c_long]

libconviqt.conviqt_pointing_data.restype = ct.POINTER(ct.c_double)
libconviqt.conviqt_pointing_data.argtypes = [ct.c_void_p]

# Convolver functions

libconviqt.conviqt_convolver_new.restype = ct.c_void_p
libconviqt.conviqt_convolver_new.argtypes = [
    ct.c_void_p,
    ct.c_void_p,
    ct.c_void_p,
    ct.c_byte,
    ct.c_long,
    ct.c_long,
    ct.c_long,
    ct.c_int,
    MPI_Comm,
]

libconviqt.conviqt_convolver_convolve.restype = ct.c_int
libconviqt.conviqt_convolver_convolve.argtypes = [ct.c_void_p, ct.c_void_p, ct.c_byte]

libconviqt.conviqt_convolver_del.restype = ct.c_int
libconviqt.conviqt_convolver_del.argtypes = [ct.c_void_p]

verbosity = 3

lmax = 32
beamlmax = lmax
beammmax = 4
pol = 0
fwhm = 4.0
beamfile = "../data/mb_lfi_30_27_x_rescaled.alm"
skyfile = "../data/slm.fits"
det_id = "LFITEST".encode()

ntheta = 3
nphi = 3
npsi = 3
nsamp = ntheta * nphi * npsi

order = 3

beam = libconviqt.conviqt_beam_new()
beam_auto = libconviqt.conviqt_beam_new()
err = libconviqt.conviqt_beam_read(
    beam, beamlmax, beammmax, pol, beamfile.encode(), comm
)
err = libconviqt.conviqt_beam_read(beam_auto, -1, -1, pol, beamfile.encode(), comm)
if err != 0:
    raise Exception("Failed to load " + beamfile)

beam_lmax = libconviqt.conviqt_beam_lmax(beam)
beam_mmax = libconviqt.conviqt_beam_mmax(beam)
beam_auto_lmax = libconviqt.conviqt_beam_lmax(beam_auto)
beam_auto_mmax = libconviqt.conviqt_beam_mmax(beam_auto)
print(
    "Read beam. lmax = {}, lmax_auto = {}, mmax = {}, mmax_auto = {}".format(
        beam_lmax, beam_auto_lmax, beam_mmax, beam_auto_mmax
    )
)

sky = libconviqt.conviqt_sky_new()
sky_auto = libconviqt.conviqt_sky_new()
err = libconviqt.conviqt_sky_read(sky, lmax, pol, skyfile.encode(), fwhm, comm)
err = libconviqt.conviqt_sky_read(sky_auto, -1, pol, skyfile.encode(), fwhm, comm)
err = libconviqt.conviqt_sky_remove_monopole(sky_auto)
err = libconviqt.conviqt_sky_remove_dipole(sky_auto)
if err != 0:
    raise Exception("Failed to load " + skyfile)

sky_lmax = libconviqt.conviqt_sky_lmax(sky)
sky_auto_lmax = libconviqt.conviqt_sky_lmax(sky_auto)
print("Read sky. lmax = {}, lmax_auto = {}".format(sky_lmax, sky_auto_lmax))

detector = libconviqt.conviqt_detector_new_with_id(det_id)

s = ct.create_string_buffer(255)
libconviqt.conviqt_detector_get_id(detector, s)
print("detector name = ", s.value)
epsilon = 1.32495160e-04
libconviqt.conviqt_detector_set_epsilon(detector, epsilon)
eps = ct.c_double(0.0)
libconviqt.conviqt_detector_get_epsilon(detector, ct.byref(eps))
print("epsilon = ", eps.value)

print("Allocating pointing array")

pnt = libconviqt.conviqt_pointing_new()
pnt_auto = libconviqt.conviqt_pointing_new()
err = libconviqt.conviqt_pointing_alloc(pnt, nsamp * 5)
err = libconviqt.conviqt_pointing_alloc(pnt_auto, nsamp * 5)
if err != 0:
    raise Exception("Failed to allocate pointing array")
ppnt = libconviqt.conviqt_pointing_data(pnt)
ppnt_auto = libconviqt.conviqt_pointing_data(pnt_auto)

print("Populating pointing array")

for p in ppnt, ppnt_auto:
    row = 0
    for theta in np.arange(ntheta) * np.pi / ntheta:
        for phi in np.arange(nphi) * 2 * np.pi / nphi:
            for psi in np.arange(npsi) * np.pi / npsi:
                p[row * 5 + 0] = phi
                p[row * 5 + 1] = theta
                p[row * 5 + 2] = psi
                p[row * 5 + 3] = 0
                p[row * 5 + 4] = row
                row += 1

for i in range(10):
    print("ppnt[", i, "] = ", ppnt[i])

print("Creating convolver")

convolver = libconviqt.conviqt_convolver_new(
    sky_auto, beam_auto, detector, pol, lmax, beammmax, order, verbosity, comm
)
convolver_auto = libconviqt.conviqt_convolver_new(
    sky, beam, detector, pol, -1, -1, order, verbosity, comm
)
if convolver == 0:
    raise Exception("Failed to instantiate convolver")

print("Convolving data")

for (cnv, p) in zip([convolver, convolver_auto], [pnt, pnt_auto]):

    calibrate = True
    err = libconviqt.conviqt_convolver_convolve(cnv, p, calibrate)
    if err != 0:
        raise Exception("Convolution FAILED!")

    # The pointer to the data will have changed during the convolution call ...

    ppnt = libconviqt.conviqt_pointing_data(p)

    if itask == 0:

        print("Convolved TOD: ")
        for row in np.arange(nsamp):
            print(ppnt[row * 5 + 4], " ", ppnt[row * 5 + 3])

        if pol:
            if np.abs(ppnt[0 * 5 + 3] - 0.8546349819096275) > 1e-6:
                raise Exception(
                    "Row 0 should be 0.8546349819096275, not " + str(ppnt[0 * 5 + 3])
                )
            if np.abs(ppnt[10 * 5 + 3] + 25.53734467183137) > 1e-6:
                raise Exception(
                    "Row 10 should be -25.53734467183137, not " + str(ppnt[10 * 5 + 3])
                )
            if np.abs(ppnt[15 * 5 + 3] + 76.04945574990082) > 1e-6:
                raise Exception(
                    "Row 15 should be -76.04945574990082, not " + str(ppnt[15 * 5 + 3])
                )
        else:
            if np.abs(ppnt[0 * 5 + 3] - 0.8545846415739397) > 1e-6:
                raise Exception(
                    "Row 0 should be 0.8545846415739397, not " + str(ppnt[0 * 5 + 3])
                )
            if np.abs(ppnt[10 * 5 + 3] + 25.20150061107036) > 1e-6:
                raise Exception(
                    "Row 10 should be -25.20150061107036, not " + str(ppnt[10 * 5 + 3])
                )
            if np.abs(ppnt[15 * 5 + 3] + 76.14723911261254) > 1e-6:
                raise Exception(
                    "Row 15 should be -76.14723911261254, not " + str(ppnt[15 * 5 + 3])
                )

libconviqt.conviqt_convolver_del(convolver)
libconviqt.conviqt_convolver_del(convolver_auto)
libconviqt.conviqt_pointing_del(pnt)
libconviqt.conviqt_detector_del(detector)
libconviqt.conviqt_sky_del(sky)
libconviqt.conviqt_sky_del(sky_auto)
libconviqt.conviqt_beam_del(beam)
libconviqt.conviqt_beam_del(beam_auto)
