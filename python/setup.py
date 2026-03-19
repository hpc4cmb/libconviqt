import glob
import os
import re
import unittest

from setuptools import setup, Extension
from setuptools.command.test import test as TestCommand

import numpy as np


def get_version():
    with open("../configure.ac", "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            if "AC_INIT" in line:
                ver = line.split()[1][1:-2]
                break
    return ver


current_version = get_version()

setup(
    name="libconviqt_wrapper",
    provides="libconviqt_wrapper",
    version=current_version,
    description="Python wrapper to libConviqt",
    author="Reijo Keskitalo",
    author_email="rtkeskitalo@lbl.gov",
    url="https://github.com/hpc4cmb/libconviqt",
    packages=["libconviqt_wrapper"],
    license="BSD",
    requires=["Python (>3.9.0)"],
)
