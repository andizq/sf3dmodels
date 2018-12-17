# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
SF3dmodels : Star-Forming regions 3d modelling package
"""
# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from .._astropy_init import *
# ----------------------------------------------------------------------------

if not _ASTROPY_SETUP_:
    # For egg_info test builds to pass, put package imports here.
    #from .example_mod import *

    from .core import Build_r, Build_theta, Build_phi

