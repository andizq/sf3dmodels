# Enforce Python version check during package import.
# This is the same check as the one at the top of setup.py
import sys

class UnsupportedPythonError(Exception):
    pass

if sys.version_info < tuple((int(val) for val in "2.7".split('.'))):
    raise UnsupportedPythonError("sf3dmodels does not support Python < {}".format(2.7))
