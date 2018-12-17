import numpy as np

#******************
#Useful TOOLS
#******************

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
    out: `numpy.ndarray`, shape (3,n) 
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
    real cells and the border of the eventual radiative transfer domain. Have a look at 
    the methods available below.
    
    Initialize this class providing the grid to fill in.

    Parameters
    ----------
    GRID : `~sf3dmodels.Model.Struct`
    """

    def __init__(self, GRID):
        self.GRID = GRID
        self.r_grid = np.linalg.norm(self.GRID.XYZ, axis = 0)
        self.r_max = np.max(self.r_grid)
        
    def random_spherical(self, r_min = 0., r_max = None, n_dummy = None):
        r"""
        Fills the grid with uniformly distributed random dummy points in a spherical section.
        
        Parameters
        ----------
        r_min : float
           Inner radius of the spherical section enclosing the random dummy points.

           Defaults to zero.

        r_max : float
           Outer radius of the spherical section enclosing the random dummy points.

           Defaults to `None`. In that case the maximum radius takes the distance of the farthest point in the grid.

        n_dummy : int
           Number of dummy points to be generated.

           Defaults to `None`. In that case n_dummy = GRID.NPoints / 100
        
        Returns
        -------
        out : dict

        Returns a dictionary with the following keys:

        r_rand : `numpy.ndarray`
           Radial coordinate of the computed random points.

        n_dummy : int
           Number of dummy points generated.
        
        See Also
        --------
        random_by_density, random_shell
        """

        if r_max == None: r_max = self.r_max
        if n_dummy == None: n_dummy = int(round(self.GRID.NPoints / 100.))
        r_rand = np.random.uniform(r_min, r_max, size = n_dummy)
        th_rand = np.random.uniform(0, np.pi, size = n_dummy)
        phi_rand = np.random.uniform(0, 2*np.pi, size = n_dummy)
        x_rand, y_rand, z_rand = spherical2cartesian(r = r_rand,
                                                     theta = th_rand,
                                                     phi = phi_rand)  
        self.GRID.XYZ = np.hstack((self.GRID.XYZ,[x_rand,y_rand,z_rand]))
        self.GRID.NPoints = self.GRID.NPoints + n_dummy
        print("New number of grid points:", self.GRID.NPoints)
        return {"r_rand": r_rand, "n_dummy": n_dummy}
    
    def random_by_mass(self, mass_field, mass_fraction = 0.5, r_max = None, n_dummy = None, r_steps = 100):
        r"""
        Fills the grid with random dummy points in a spherical section according to a mass criterion.
        
        The inner radius of the spherical section will be equal to the radius of a sphere which encloses the input ``mass_fraction``.
        The function iterates over the radial coordinate until the enclosed mass fraction surpass the input ``mass_fraction``.
        
        Parameters
        ----------
        mass_field: array_like
           `list` or `numpy.ndarray` containing the mass of each cell, where the i-th mass corresponds to the i-th cell of the GRID.
        
        fraction : float
           The mass fraction to be enclosed by the inner radius of the spherical section.
           Sets the inner radius r_min at `random_spherical`.

           Defaults to 0.5, i.e, by default the computed inner radius will enclose (approx.) the 50% of the total mass.
        
        r_max : float
           Outer radius of the spherical section enclosing the random dummy points.

           Defaults to `None`. In that case it takes the distance of the farthest point in the grid.

        n_dummy : int
           Number of dummy points to be generated.

           Defaults to `None`. In that case n_dummy = GRID.NPoints / 100

        r_steps : int
           Number of divisions along the radial coordinate to compute the enclosed mass.
           
           Defaults to 100.

        Returns
        -------
        out : dict
        
        Returns a dictionary with the following keys:
        
        r_rand : `numpy.ndarray`
           Radial coordinate of the computed random points.
        
        r_min : float
           Computed inner radius of the spherical section which encloses (approximately) the input mass fraction and above which the random points are generated.

        r_max : float
           Outer radius of the spherical section enclosing the random dummy points.

        comp_fraction : float
           Computed mass fraction enclosed by r_min. The higher r_steps the closer this value to the wished mass fraction. 
           
        n_dummy : int
           Number of dummy points generated.
        
        See Also
        --------
        random_spherical
        """
        if r_max == None: r_max = self.r_max

        mass_field = _list2array(mass_field)
        total_mass = np.sum(mass_field)
        min_mass = mass_fraction * total_mass
        enc_mass = 0
        r_min = None
        for r in np.linspace(0,r_max,r_steps):
            inds = np.where(self.r_grid < r)
            enc_mass = np.sum(mass_field[inds])
            if enc_mass > min_mass: 
                r_min = r
                break
        print r_min    
        comp_fraction = enc_mass / total_mass
        rand_spher = self.random_spherical(r_min=r_min, r_max=r_max, n_dummy=n_dummy)
        return {"r_rand": rand_spher["r_rand"], "r_min": r_min, "r_max": r_max, 
                "comp_fraction": comp_fraction, "n_dummy": rand_spher["n_dummy"]}


    def random_by_density(self):
        r"""
        Under development.

        Fills the grid with random points weighted by the inverse of the density field.
        """
        t = 4+4
        return t

    def random_shell(self):
        r"""
        Under development.
        
        Fills the grid with random points in a spherical shell. 
        """
        pass
