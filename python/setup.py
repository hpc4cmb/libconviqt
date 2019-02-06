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


# run unit tests


class PTestCommand(TestCommand):
    def __init__(self, *args, **kwargs):
        super(PTestCommand, self).__init__(*args, **kwargs)

    def initialize_options(self):
        TestCommand.initialize_options(self)

    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_suite = True

    def run(self):
        loader = unittest.TestLoader()
        runner = unittest.TextTestRunner(verbosity=2)
        suite = loader.discover("tests", pattern="test_*.py", top_level_dir=".")
        runner.run(suite)


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
    requires=["Python (>3.4.0)"],
    cmdclass={"test": PTestCommand},
)
