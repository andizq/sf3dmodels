import numpy as np
from .core import GridSet
#*************
#GRID BUILDER
#*************

class Build_r(GridSet):
    _base = """
    Computes the spherical coordinate :math:`r` on each grid point.

       .. math:: r = \\sqrt{x^2 + y^2 + z^2}
    """
    _pars = GridSet._pars
    _returns = """
    Returns
    -------
    New attributes:
    
    GRID.r : `numpy.ndarray`, shape ``(GRID.NPoints,)``
       Array of :math:`r` in `input units`.
    """
    __doc__ = _base + _pars + _returns

    def __init__(self, GRID):
        super(Build_r,self).__init__(GRID)
        self._r_grid = np.linalg.norm(self.GRID.XYZ, axis = 0)
        self._r_max = np.max(self._r_grid)
        self.GRID.r = self._r_grid
        self._set_flag('r')

class Build_theta(GridSet):
    _base = """
    Computes the polar angle :math:`\\theta` on each grid point.
       .. math:: 
          \\theta = \\cos^{-1}\\left(\\frac{z}{r}\\right) 
    """
    _pars = GridSet._pars
    _returns = """
    Returns
    -------
    New attributes:
    
    GRID.theta : `numpy.ndarray`, shape ``(GRID.NPoints,)``
       Array of angles in radians, in the range `[0, pi/2]` if :math:`z>=0` and 
       `[pi/2, 0]` if :math:`z<0`; where :math:`\\theta=\\pi/2` for the plane :math:`z=0`.
    """
    __doc__ = _base + _pars + _returns

    def __init__(self, GRID):
        self.GRID = GRID
        self.r_grid = np.linalg.norm(self.GRID.XYZ, axis = 0)
        self.r_max = np.max(self.r_grid)
        
class Build_phi(GridSet):
    _base = """
    Computes the azimuthal angle :math:`\\phi` on each grid point. Uses `numpy.arctan2`.

       .. math:: \\phi = \\tan^{-1}\\left(\\frac{y}{x}\\right) + 2\\pi
    """
    _pars = GridSet._pars
    _returns = """
    Returns
    -------
    New attributes:
    
    GRID.phi : `numpy.ndarray`, shape ``(GRID.NPoints,)``
       Array of angles in `radians`, in the range `[0, pi/2]`.
    """
    __doc__ = _base + _pars + _returns

    def __init__(self, GRID):
        self.GRID = GRID
        self.r_grid = np.linalg.norm(self.GRID.XYZ, axis = 0)
        self.r_max = np.max(self.r_grid)



