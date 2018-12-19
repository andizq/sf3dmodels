import numpy as np

#*************
#GRID BUILDER
#*************

class Spherical(object):
    def __init__(self, GRID):
        self.GRID = GRID
    def set_flag(self, a):
        self.GRID.flag[a] = True

class Build_r(Spherical):
    _base = """
    Computes the :math:`r` - spherical coordinate on each grid point.

       .. math:: r = \\sqrt{x^2 + y^2 + z^2}
    """
    _pars = """
    Parameters
    ----------
    GRID : `~sf3dmodels.Model.Struct`
       Grid structure to compute.
   
    get_r : bool
       If True, stores the computed :math:`r` - coordinate as a new attribute ``GRID.r``.

       Defaults to %r.
    """       
    _def0 = True
    _returns = """
    Returns
    -------
    New attributes:
    
    GRID.r : `numpy.ndarray`, shape (GRID.NPoints,)
       Array of :math:`r` in the input units.
    """
    __doc__ = _base + _pars%_def0 + _returns

    def __init__(self, GRID, get_r=_def0):
        super(Build_r,self).__init__(GRID)
        self._r_grid = np.linalg.norm(self.GRID.XYZ, axis = 0)
        self._r_max = np.max(self._r_grid)
        if get_r: 
            self.GRID.r = self._r_grid
            self.set_flag('r')

class Build_theta(Spherical):
    _base = """
    Computes :math:`\\theta` - polar angle - on each grid point.
       .. math:: 
          \\theta = \\cos^{-1}\\left(\\frac{z}{r}\\right) 
    """
    _pars = """
    Parameters
    ----------
    GRID : `~sf3dmodels.Model.Struct`
       Grid structure to compute.
   
    get_theta : bool
       If True, stores the computed :math:`\\theta` - coordinate as a new attribute ``GRID.theta``.

       Defaults to %r.
    """       
    _def0 = True
    _returns = """
    Returns
    -------
    New attributes:
    
    GRID.theta : `numpy.ndarray`, shape (GRID.NPoints,)
       Array of angles in radians, in the range `[0, pi/2]` if :math:`z>=0` and 
       `[pi/2, 0]` if :math:`z<0`; where :math:`\\theta=\\pi/2` for the plane :math:`z=0`.
    """
    __doc__ = _base + _pars%_def0 + _returns

    def __init__(self, GRID):
        self.GRID = GRID
        self.r_grid = np.linalg.norm(self.GRID.XYZ, axis = 0)
        self.r_max = np.max(self.r_grid)
        
class Build_phi(Spherical):
    _base = """
    Computes :math:`\\phi` - azimuthal angle - on each grid point. Uses `numpy.arctan2`.

       .. math:: \\phi = \\tan^{-1}\\left(\\frac{y}{x}\\right) + 2\\pi
    """
    _pars = """
    Parameters
    ----------
    GRID : `~sf3dmodels.Model.Struct`
       Grid structure to compute.
   
    get_theta : bool
       If True, stores the computed :math:`\\phi` - coordinate as a new attribute ``GRID.phi``.

       Defaults to %r.
    """       
    _def0 = True
    _returns = """
    Returns
    -------
    New attributes:
    
    GRID.phi : `numpy.ndarray`, shape (GRID.NPoints,)
       Array of angles in radians, in the range `[0, pi/2]`.
    """
    __doc__ = _base + _pars%_def0 + _returns

    def __init__(self, GRID):
        self.GRID = GRID
        self.r_grid = np.linalg.norm(self.GRID.XYZ, axis = 0)
        self.r_max = np.max(self.r_grid)



