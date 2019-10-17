# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package initializes the model hosting grid (Under development).
"""
# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from .._astropy_init import *
# ----------------------------------------------------------------------------

class GridInit(object):
    def __init__(self):
        pass

class GridSet(object):
    """
    Base class for all classes invoking grids to operate with.
    """
    _pars = """
    Parameters
    ----------
    GRID : `~sf3dmodels.Model.Struct`
       Grid structure to compute.
    """
    __doc__ += _pars
    def __init__(self, GRID): self.GRID = GRID
    def _set_flag(self, a): 
        try: 
            self.GRID.flag[a] = True
        except AttributeError: #To do: Call flag setter class
            self.GRID.flag = {}
            self.GRID.flag[a] = True

if not _ASTROPY_SETUP_:
    # For egg_info test builds to pass, put package imports here.
    #from .example_mod import *

    from . import fillgrid
    from .core import Grid, Overlap, RandomGridAroundAxis, Build_r, Build_theta, Build_phi, NeighbourRegularGrid
    
__all__ = ['Grid',
           'NeighbourRegularGrid',
           'Overlap',
           'RandomGridAroundAxis',
           'Build_r', 'Build_theta', 'Build_phi',
           'GridSet'] 

