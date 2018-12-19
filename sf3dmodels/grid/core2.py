import numpy as np
from . import GridInit

class GridSet(object):
    """
    Base class for calling the grid and setting new flags.
    """
    _pars = """
    Parameters
    ----------
    GRID : `~sf3dmodels.Model.Struct`
       Grid structure to compute.
    """
    __doc__ += _pars
    def __init__(self, GRID): self.GRID = GRID
    def _set_flag(self, a): self.GRID.flag[a] = True
