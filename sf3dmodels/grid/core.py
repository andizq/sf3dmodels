import numpy as np

#*************
#GRID BUILDER
#*************

class Build_r(object):
    _base = """
    Computes the r - spherical coordinate on each grid point.

    .. math:: r = \\sqrt{x^2 + y^2 + z^2}
    """
    _pars = """
    Parameters
    ----------
    GRID : `~sf3dmodels.Model.Struct`
       Grid structure to compute.
   
    get_r : bool
       If True, stores the computed r coordinate in a new attribute ``GRID.r``

       Defaults to %r
    """       
    _def0 = True
    _returns = """
    Returns
    -------
    New attributes:
    
    GRID.r : `numpy.ndarray`, shape (GRID.NPoints,)
       Computed r - spherical coordinate.
    """
    __doc__ = _base + _pars%_def0 + _returns

    def __init__(self, GRID, get_r=_def0):
        self.GRID = GRID
        self._r_grid = np.linalg.norm(self.GRID.XYZ, axis = 0)
        self._r_max = np.max(self._r_grid)
        if get_r: self.GRID.r = self._r_grid
    
class Build_theta(object):
    def __init__(self, GRID):
        self.GRID = GRID
        self.r_grid = np.linalg.norm(self.GRID.XYZ, axis = 0)
        self.r_max = np.max(self.r_grid)
        
class Build_phi(object):
    def __init__(self, GRID):
        self.GRID = GRID
        self.r_grid = np.linalg.norm(self.GRID.XYZ, axis = 0)
        self.r_max = np.max(self.r_grid)



