# Copyright (c) 2019 by the parties listed in the AUTHORS
# file.  All rights reserved.  Use of this source code is governed
# by a BSD-style license that can be found in the LICENSE file.

from __future__ import division
from __future__ import print_function

from mpi4py import MPI
import ctypes as ct
import numpy as np
import os
import sys

from unittest import TestCase

import libconviqt_wrapper as conviqt


class ConviqtTest(TestCase):

    comm = MPI.COMM_WORLD
    itask = comm.Get_rank()
    ntask = comm.Get_size()

    if itask == 0:
        print("Running with ", ntask, " MPI tasks")

    verbosity = 3

    lmax = 32
    beamlmax = lmax
    beammmax = 4
    pol = 0
    fwhm = 4.0
    beamfile = "../data/mb_lfi_30_27_x_rescaled.alm"
    skyfile = "../data/slm.fits"
    det_id = "LFITEST"

    ntheta = 3
    nphi = 3
    npsi = 3
    nsamp = ntheta * nphi * npsi

    order = 3

    beam = conviqt.Beam(beamlmax, beammmax, pol, beamfile, comm)
    beam_auto = conviqt.Beam(-1, -1, pol, beamfile, comm)

    beam_lmax = beam.lmax()
    beam_mmax = beam.mmax()
    beam_auto_lmax = beam_auto.lmax()
    beam_auto_mmax = beam_auto.mmax()
    print(
        "Read beam. lmax = {}, lmax_auto = {}, mmax = {}, mmax_auto = {}".format(
            beam_lmax, beam_auto_lmax, beam_mmax, beam_auto_mmax
        )
    )

    sky = conviqt.Sky(lmax, pol, skyfile, fwhm, comm)
    sky_auto = conviqt.Sky(-1, pol, skyfile, fwhm, comm)
    sky_auto.remove_monopole()
    sky_auto.remove_dipole()

    sky_lmax = sky.lmax()
    sky_auto_lmax = sky_auto.lmax()
    print("Read sky. lmax = {}, lmax_auto = {}".format(sky_lmax, sky_auto_lmax))

    detector = conviqt.Detector(det_id)

    print("detector name = ", detector.get_id())
    epsilon = 1.32495160e-04
    detector.set_epsilon(epsilon)
    eps = detector.get_epsilon()
    print("epsilon = ", eps)

    print("Allocating pointing array")

    pnt = conviqt.Pointing()
    pnt.alloc(nsamp)
    pnt_auto = conviqt.Pointing(nsamp)
    ppnt = pnt.data()
    ppnt_auto = pnt_auto.data()

    print("Populating pointing array")

    for p in ppnt, ppnt_auto:
        row = 0
        for theta in np.arange(ntheta) * np.pi / ntheta:
            for phi in np.arange(nphi) * 2 * np.pi / nphi:
                for psi in np.arange(npsi) * np.pi / npsi:
                    p[row, 0] = phi
                    p[row, 1] = theta
                    p[row, 2] = psi
                    p[row, 3] = 0
                    p[row, 4] = row
                    row += 1

    for i in range(10):
        print("ppnt[", i, "] = ", ppnt[i])

    print("Creating convolver")

    convolver = conviqt.Convolver(
        sky_auto, beam_auto, detector, pol, lmax, beammmax, order, verbosity, comm
    )
    convolver_auto = conviqt.Convolver(
        sky, beam, detector, pol, -1, -1, order, verbosity, comm
    )

    print("Convolving data")

    for (cnv, p) in zip([convolver, convolver_auto], [pnt, pnt_auto]):

        calibrate = True
        cnv.convolve(p, calibrate)

        # The pointer to the data will have changed during the convolution call ...

        ppnt = p.data()

        if itask == 0:

            print("Convolved TOD: ")
            for row in np.arange(nsamp):
                print(ppnt[row, 4], " ", ppnt[row, 3])

            if pol:
                if np.abs(ppnt[0, 3] - 0.8546349819096275) > 1e-6:
                    raise Exception(
                        "Row 0 should be 0.8546349819096275, not " + str(ppnt[0, 3])
                    )
                if np.abs(ppnt[10, 3] + 25.53734467183137) > 1e-6:
                    raise Exception(
                        "Row 10 should be -25.53734467183137, not " + str(ppnt[10, 3])
                    )
                if np.abs(ppnt[15, 3] + 76.04945574990082) > 1e-6:
                    raise Exception(
                        "Row 15 should be -76.04945574990082, not " + str(ppnt[15, 3])
                    )
            else:
                if np.abs(ppnt[0, 3] - 0.8545846415739397) > 1e-6:
                    raise Exception(
                        "Row 0 should be 0.8545846415739397, not " + str(ppnt[0, 3])
                    )
                if np.abs(ppnt[10, 3] + 25.20150061107036) > 1e-6:
                    raise Exception(
                        "Row 10 should be -25.20150061107036, not " + str(ppnt[10, 3])
                    )
                if np.abs(ppnt[15, 3] + 76.14723911261254) > 1e-6:
                    raise Exception(
                        "Row 15 should be -76.14723911261254, not " + str(ppnt[15, 3])
                    )

    del convolver
    del convolver_auto
    del pnt
    del detector
    del sky
    del sky_auto
    del beam
    del beam_auto
