# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package contains useful tools to operate on data 
from the hybrid (SPH+AMR) code `AREPO`_ (Under development).
"""
# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from .._astropy_init import *
# ----------------------------------------------------------------------------

if not _ASTROPY_SETUP_:
    # For egg_info test builds to pass, put package imports here.
    from .core import UniqueCells, ReadSnapshot, ArepoTags

__all__ = ['UniqueCells', 'ReadSnapshot', 'ArepoTags']

