import numpy as np
from . import GridInit

__all__ = ['GridSet']
class GridSet(object):
    """
    Base class for all classes invoking non-new grids.
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
        except AttributeError: #Call flag setter class
            self.GRID.flag = {}
            self.GRID.flag[a] = True
