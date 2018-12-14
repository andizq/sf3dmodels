import numpy as np

#*******************
#Useful TOOLS
#*******************

def _list2array(a):
    """
    Converts the input list (if it's the case) to an np.ndarray 

    Parameters
    ----------
    a : list
       list to be converted

    Returns
    -------
    out: np.array(a)
    """

    if isinstance(a,list): return np.array(a)
    else: return a 

def spherical2cartesian(r = None, theta = None, phi = None):
    """
    Converts spherical coordinates to rectangular cartesian coordinates.
    
    Parameters
    ----------
    r : scalar or array_like, shape (n,) 
       Radial distance from the origin of coordinates to the point in question.
    
    theta : scalar or array_like, shape (n,) 
       Polar angle referred to the positive `z` - axis.
    
    phi : scalar or array_like, shape (n,) 
       Azimuthal angle referred to the positive `x` - axis.
       
    Returns
    -------
    out: numpy.ndarray, shape (3,n) 
       Array-like object with the resulting x,y,z transformation.
    
    Notes
    -----
    Transformation equations:
    
    .. math:: x &= r\\sin(\\theta)\\cos(\\phi) \n 
              y &= r\\sin(\\theta)\\sin(\\phi) \n
              z &= r\\cos(\\theta)
    """

    r = _list2array(r)
    theta = _list2array(theta)
    phi = _list2array(phi)
    
    return np.array([
            r * np.sin(theta) * np.cos(phi),
            r * np.sin(theta) * np.sin(phi),
            r * np.cos(theta)
            ])

def cartesian2spherical(x=None, y=None, z=None):
    return 0

class Fill_grid(object):
    r"""
    Fills the grid in with dummy (empty) points.

    This is particularly useful for irregular grids containing large voids between its
    real cells and the border of the eventual radiative transfer domain.   
    """

    def __init__(self, GRID):
        self.GRID = GRID
        self.r_grid = np.linalg.norm(self.GRID.XYZ, axis = 0)
        self.r_max = np.max(self.r_grid)
        
    def random(self, r_dummy = False, n_dummy = False):
        r"""
        Fills the grid in with uniformly distributed random points.
        
        Parameters
        ----------
        r_dummy : float
           Radius of the sphere enclosing the random dummy points.

           Defaults to False. In that case the radius takes the distance value of the farthest point in the grid.
        n_dummy : float
           Number of dummy points to be generated.

           Defaults to False. In that case n_dummy = GRID.NPoints / 100.
        
        Returns
        -------
        n_dummy : float
        
        See Also
        --------
        random_by_density, random_shell
        """

        if not r_dummy: r_dummy = self.r_max
        if not n_dummy: n_dummy = int(round(self.GRID.NPoints / 100.))
        r_rand = np.random.uniform(0, r_dummy, size = n_dummy)
        th_rand = np.random.uniform(0, np.pi, size = n_dummy)
        phi_rand = np.random.uniform(0, 2*np.pi, size = n_dummy)
        x_rand, y_rand, z_rand = spherical2cartesian(r = r_rand,
                                                     theta = th_rand,
                                                     phi = phi_rand)  
        self.GRID.XYZ = np.hstack((self.GRID.XYZ,[x_rand,y_rand,z_rand]))
        self.GRID.NPoints = self.GRID.NPoints + n_dummy
        print("New number of grid points:", self.GRID.NPoints)
        return n_dummy 
    
    def random_by_density(self):
        r"""
        Under development.

        Fills the grid in with random points weighted by the inverse of the density field.
        """
        t = 4+4
        return t

    def random_shell(self):
        r"""
        Under development.
        
        Fills the grid with random points in a spherical shell. 
        """
        pass
