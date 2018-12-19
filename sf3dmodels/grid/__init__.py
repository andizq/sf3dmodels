# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
SF3dmodels : Star-Forming regions 3d modelling package
"""
# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from .._astropy_init import *
# ----------------------------------------------------------------------------

class GridInit(object):
    def __init__(self):
        pass

if not _ASTROPY_SETUP_:
    # For egg_info test builds to pass, put package imports here.
    #from .example_mod import *

    #from .core import GridSet
    from .spherical import Build_r, Build_theta, Build_phi

__all__ = ['Build_r', 'Build_theta', 'Build_phi']

