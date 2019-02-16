# Copyright (c) 2019 by the parties listed in the AUTHORS
# file.  All rights reserved.  Use of this source code is governed
# by a BSD-style license that can be found in the LICENSE file.

from __future__ import division
from __future__ import print_function

import ctypes as ct
import ctypes.util as ctu
import os
import sys
import glob
import shutil

import numpy as np
import numpy.ctypeslib as npc

from mpi4py import MPI


try:
    _conviqt = ct.CDLL("_conviqt.so")
except OSError:
    path = ctu.find_library("conviqt")
    if path is not None:
        _conviqt = ct.CDLL(path)

available = _conviqt is not None

try:
    if MPI._sizeof(MPI.Comm) == ct.sizeof(ct.c_int):
        MPI_Comm = ct.c_int
    else:
        MPI_Comm = ct.c_void_p
except Exception as e:
    raise Exception(
        'Failed to set the portable MPI communicator datatype: "{}". '
        "MPI4py is probably too old. ".format(e)
    )


def encode_comm(comm):
    comm_ptr = MPI._addressof(comm)
    return MPI_Comm.from_address(comm_ptr)


# Beam functions

_conviqt.conviqt_beam_new.restype = ct.c_void_p
_conviqt.conviqt_beam_new.argtypes = []

_conviqt.conviqt_beam_del.restype = ct.c_int
_conviqt.conviqt_beam_del.argtypes = [ct.c_void_p]

_conviqt.conviqt_beam_read.restype = ct.c_int
_conviqt.conviqt_beam_read.argtypes = [
    ct.c_void_p,
    ct.c_long,
    ct.c_long,
    ct.c_byte,
    ct.c_char_p,
    MPI_Comm,
]

_conviqt.conviqt_beam_lmax.restype = ct.c_int
_conviqt.conviqt_beam_lmax.argtypes = [ct.c_void_p]

_conviqt.conviqt_beam_mmax.restype = ct.c_int
_conviqt.conviqt_beam_mmax.argtypes = [ct.c_void_p]

_conviqt.conviqt_beam_normalize.restype = ct.c_double
_conviqt.conviqt_beam_normalize.argtypes = [ct.c_void_p]


class Beam(object):
    """
    Conviqt beam expansion object
    """

    def __init__(self, lmax, mmax, pol, beamfile, comm):
        self._beam = _conviqt.conviqt_beam_new()
        err = _conviqt.conviqt_beam_read(
            self._beam, lmax, mmax, pol, beamfile.encode(), encode_comm(comm)
        )
        if err != 0:
            raise RuntimeError("Failed to load {}".format(beamfile))

    def __del__(self):
        try:
            err = _conviqt.conviqt_beam_del(self._beam)
            if err != 0:
                raise RuntimeError("Unknown error")
        except Exception as e:
            print('WARNING: failed to free conviqt beam: "{}"'.format(e), flush=True)

    def lmax(self):
        return _conviqt.conviqt_beam_lmax(self._beam)

    def mmax(self):
        return _conviqt.conviqt_beam_mmax(self._beam)

    def normalize(self):
        scale = _conviqt.conviqt_beam_normalize(self._beam)
        if scale <= 0:
            raise RuntimeError("Failed to normalize beam. scale = {}".format(scale))
        return scale


# Sky functions

_conviqt.conviqt_sky_new.restype = ct.c_void_p
_conviqt.conviqt_sky_new.argtypes = []

_conviqt.conviqt_sky_del.restype = ct.c_int
_conviqt.conviqt_sky_del.argtypes = [ct.c_void_p]

_conviqt.conviqt_sky_read.restype = ct.c_int
_conviqt.conviqt_sky_read.argtypes = [
    ct.c_void_p,
    ct.c_long,
    ct.c_byte,
    ct.c_char_p,
    ct.c_double,
    MPI_Comm,
]

_conviqt.conviqt_sky_lmax.restype = ct.c_int
_conviqt.conviqt_sky_lmax.argtypes = [ct.c_void_p]

_conviqt.conviqt_sky_remove_monopole.restype = ct.c_int
_conviqt.conviqt_sky_remove_monopole.argtypes = [ct.c_void_p]

_conviqt.conviqt_sky_remove_dipole.restype = ct.c_int
_conviqt.conviqt_sky_remove_dipole.argtypes = [ct.c_void_p]


class Sky(object):
    """
    Conviqt sky expansion object
    """

    def __init__(self, lmax, pol, skyfile, fwhm, comm):
        self._sky = _conviqt.conviqt_sky_new()
        err = _conviqt.conviqt_sky_read(
            self._sky, lmax, pol, skyfile.encode(), fwhm, encode_comm(comm)
        )
        if err != 0:
            raise RuntimeError("Failed to load {}".format(skyfile))

    def __del__(self):
        try:
            err = _conviqt.conviqt_sky_del(self._sky)
            if err != 0:
                raise RuntimeError("Unknown error")
        except Exception as e:
            print('WARNING: failed to free conviqt sky: "{}"'.format(e), flush=True)

    def lmax(self):
        return _conviqt.conviqt_sky_lmax(self._sky)

    def remove_monopole(self):
        err = _conviqt.conviqt_sky_remove_monopole(self._sky)

    def remove_dipole(self):
        err = _conviqt.conviqt_sky_remove_dipole(self._sky)


# Detector functions

_conviqt.conviqt_detector_new.restype = ct.c_void_p
_conviqt.conviqt_detector_new.argtypes = []

_conviqt.conviqt_detector_new_with_id.restype = ct.c_void_p
_conviqt.conviqt_detector_new_with_id.argtypes = [ct.c_char_p]

_conviqt.conviqt_detector_del.restype = ct.c_int
_conviqt.conviqt_detector_del.argtypes = [ct.c_void_p]

_conviqt.conviqt_detector_set_epsilon.restype = ct.c_int
_conviqt.conviqt_detector_set_epsilon.argtypes = [ct.c_void_p, ct.c_double]

_conviqt.conviqt_detector_get_epsilon.restype = ct.c_int
_conviqt.conviqt_detector_get_epsilon.argtypes = [ct.c_void_p, ct.POINTER(ct.c_double)]

_conviqt.conviqt_detector_get_id.restype = ct.c_int
_conviqt.conviqt_detector_get_id.argtypes = [ct.c_void_p, ct.c_char_p]


class Detector(object):
    """
    Conviqt detector object
    """

    def __init__(self, name=None, epsilon=None):
        if name is None:
            self._det = _conviqt.conviqt_detector_new()
        else:
            self._det = _conviqt.conviqt_detector_new_with_id(name.encode())
        if epsilon is not None:
            self.set_epsilon(epsilon)

    def __del__(self):
        try:
            err = _conviqt.conviqt_detector_del.restype(self._det)
            if err != 0:
                raise RuntimeError("Unknown error")
        except Exception as e:
            print(
                'WARNING: failed to free conviqt detector: "{}"'.format(e), flush=True
            )

    def set_epsilon(self, epsilon):
        err = _conviqt.conviqt_detector_set_epsilon(self._det, epsilon)
        if err != 0:
            raise RuntimeError("Failed to set detector epsilon")

    def get_epsilon(self):
        eps = ct.c_double(0.0)
        err = _conviqt.conviqt_detector_get_epsilon(self._det, ct.byref(eps))
        if err != 0:
            raise RuntimeError("Failed to get detector epsilon")
        return eps.value

    def get_id(self):
        s = ct.create_string_buffer(1024)
        _conviqt.conviqt_detector_get_id(self._det, s)
        return s.value


# Pointing functions

_conviqt.conviqt_pointing_new.restype = ct.c_void_p
_conviqt.conviqt_pointing_new.argtypes = []

_conviqt.conviqt_pointing_del.restype = ct.c_int
_conviqt.conviqt_pointing_del.argtypes = [ct.c_void_p]

_conviqt.conviqt_pointing_alloc.restype = ct.c_int
_conviqt.conviqt_pointing_alloc.argtypes = [ct.c_void_p, ct.c_long]

_conviqt.conviqt_pointing_data.restype = ct.POINTER(ct.c_double)
_conviqt.conviqt_pointing_data.argtypes = [ct.c_void_p]


class Pointing(object):
    """
    Conviqt pointing object
    """

    _nrow = 0
    _ncol = 5

    def __init__(self, nrow=None):
        self._pnt = _conviqt.conviqt_pointing_new()
        if nrow is not None:
            self.alloc(nrow)

    def __del__(self):
        try:
            err = _conviqt.conviqt_pointing_del(self._pnt)
            if err != 0:
                raise RuntimeError("Unknown error")
        except Exception as e:
            print(
                'WARNING: failed to free conviqt pointing: "{}"'.format(e), flush=True
            )

    def alloc(self, nrow):
        self._nrow = nrow
        err = _conviqt.conviqt_pointing_alloc(self._pnt, self._nrow * self._ncol)
        if err != 0:
            raise RuntimeError("Failed to allocate {} rows to pointing".format(nrow))

    def data(self):
        buf = _conviqt.conviqt_pointing_data(self._pnt)
        arr = np.ctypeslib.as_array(buf, shape=(self._nrow, self._ncol))
        return arr


# Convolver functions

_conviqt.conviqt_convolver_new.restype = ct.c_void_p
_conviqt.conviqt_convolver_new.argtypes = [
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

_conviqt.conviqt_convolver_convolve.restype = ct.c_int
_conviqt.conviqt_convolver_convolve.argtypes = [ct.c_void_p, ct.c_void_p, ct.c_byte]

_conviqt.conviqt_convolver_del.restype = ct.c_int
_conviqt.conviqt_convolver_del.argtypes = [ct.c_void_p]


class Convolver(object):
    """
    Conviqt convolver object
    """

    def __init__(self, sky, beam, detector, pol, lmax, mmax, order, verbosity, comm):
        self._conv = _conviqt.conviqt_convolver_new(
            sky._sky,
            beam._beam,
            detector._det,
            pol,
            lmax,
            mmax,
            order,
            verbosity,
            encode_comm(comm),
        )
        if self._conv == 0:
            raise RuntimeError("Failed to instantiate convolver")

    def __del__(self):
        try:
            err = _conviqt.conviqt_convolver_del(self._conv)
            if err != 0:
                raise RuntimeError("Unknown error")
        except Exception as e:
            print(
                'WARNING: failed to free conviqt convolver: "{}"'.format(e), flush=True
            )

    def convolve(self, pointing, calibrate=True):
        err = _conviqt.conviqt_convolver_convolve(self._conv, pointing._pnt, calibrate)
        if err != 0:
            raise RuntimeError("Failed to convolve")
