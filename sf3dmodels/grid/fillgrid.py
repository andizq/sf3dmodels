"""
Fills up the input grid with empty points (dummy points).

This is particularly useful for irregular grids that have large voids between the
real cells and the border of an eventual radiative transfer domain.
"""
import numpy as np
from ..tools.transform import spherical2cartesian
from .core import Build_r
from copy import copy, deepcopy

__all__ = ['Random']
#******************
#Useful TOOLS
#******************
class SmartRejectionDummies(object):
    def __init__(self, pregrid, n_real, n_dummy, delaunay=False):
        tri = delaunay
        if not delaunay:
            from scipy.spatial import Delaunay
            tri = Delaunay(pregrid.T)

        self.n_real = n_real
        self.n_dummy = n_dummy
        self.indices, self.indptr = tri.vertex_neighbor_vertices
        self.accepted_dummies = [] #np.zeros(n_dummy, dtype=bool)
        self.inspect_dummies()
        self.accepted_dummies = np.asarray(self.accepted_dummies)
        self.n_accepted_dummies = len(self.accepted_dummies)

    def find_dummy_neighs(self, k):
        ind_neighs = self.indptr[self.indices[k]:self.indices[k+1]]
        if np.sum(ind_neighs >= self.n_real) < 2: pass #If there is just one neighbour dummy, then reject the k-th dummy.
        else: self.accepted_dummies.append(k)
        
    def inspect_dummies(self):
        for k in range(self.n_real, self.n_real+self.n_dummy):
            self.find_dummy_neighs(k)
        
class Random(Build_r, SmartRejectionDummies): 
    """    
    Methods for filling in the grid randomly with (non-physical) dummy points.
    """
    __doc__ += Build_r._pars + Build_r._returns
    
    def __init__(self, grid, smart=True):
        self._smart = smart
        self.grid_orig = copy(grid)
        super(Random, self).__init__(grid)

    def by_density(self, density, mass_fraction = 0.5, r_max = None, n_dummy = None, r_steps = 100):
        r"""
        Under development.

        Fills up the grid with uniformly-distributed random points according to a mass criterion, given a density field.

        Useful when the model does provide the density but not the mass. However the aim is identical to the method `by_mass`.

        The mass is calculated via the mean density enclosed in a sphere growing on each iteration until the ``mass_fraction`` 
        is reached, which sets the inner radius of the spherical section where the dummy points will be generated.  
        
        Parameters
        ----------
        density : array_like
           `list` or `numpy.ndarray` with the density of each cell, where the i-th density corresponds to the i-th cell of the GRID.
        
        mass_fraction : float
           The mass fraction to be enclosed by the inner radius of the spherical section.
           Sets the inner radius ``r_min`` at `spherical`.

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
        In prep.
        
        See Also
        --------
        by_mass, spherical
        """
        t = 4+4
        return t

    def spherical(self, prop, prop_fill = {}, r_min = 0., r_max = None, n_dummy = None):
        r"""
        Fills up the grid with uniformly-distributed random dummy points in a spherical section.
        
        Parameters
        ----------
        prop : dict
           Dictionary of arrays of physical properties to be filled up by dummy values. 

        prop_fill : dict, optional
           Dictionary specifying the dummy value (scalar) with which each physical property will be filled up.
           
           Defaults to {}. In that case, the filling dummy value is zero for all the properties. 

        r_min : float, optional
           Inner radius of the spherical section enclosing the random dummy points.

           Defaults to zero.

        r_max : float, optional
           Outer radius of the spherical section enclosing the random dummy points.

           Defaults to `None`. In that case, the maximum radius takes the distance of the farthest point in the grid from the coordinates centre.

        n_dummy : int, optional
           Number of dummy points to be generated.

           Defaults to `None`. In that case, n_dummy = GRID.NPoints / 100

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
        box
        """

        if r_max == None: r_max = self._r_max
        if n_dummy == None: n_dummy = int(round(self.GRID.NPoints / 100.))
        else: n_dummy = int(n_dummy)
        r_rand = np.random.uniform(r_min, r_max, size = n_dummy)
        th_rand = np.random.uniform(0, np.pi, size = n_dummy)
        phi_rand = np.random.uniform(0, 2*np.pi, size = n_dummy)
        x_rand, y_rand, z_rand = spherical2cartesian(r = r_rand,
                                                     theta = th_rand,
                                                     phi = phi_rand)  
        
        self.pregrid = np.hstack((self.GRID.XYZ,[x_rand,y_rand,z_rand]))
        if self._smart: 
            SmartRejectionDummies.__init__(self, self.pregrid, self.GRID.NPoints, n_dummy) 
            n_dummy = self.n_accepted_dummies
            accepted_dummies_ids = self.accepted_dummies - self.GRID.NPoints
            self.GRID.XYZ = np.hstack((self.GRID.XYZ,[x_rand[accepted_dummies_ids],
                                                      y_rand[accepted_dummies_ids],
                                                      z_rand[accepted_dummies_ids]]))
        else: self.GRID.XYZ = self.pregrid

        self.GRID.NPoints = self.GRID.NPoints + n_dummy

        print("New number of grid points:", self.GRID.NPoints)

        dummies = np.zeros(n_dummy)
        
        fill = {}
        for p in prop: fill[p] = 0.0
        fill.update(prop_fill)
        for p in prop: prop[p] = np.append(prop[p], dummies+fill[p]) 
                
        return {"r_rand": r_rand, "n_dummy": n_dummy}


    def by_mass(self, mass, prop, prop_fill = {}, mass_fraction = 0.5, r_max = None, n_dummy = None, r_steps = 100):
        r"""
        Fills up the grid with uniformly-distributed random dummy points in a spherical section based on a mass threshold.
        
        The inner radius of the spherical section will be equal to the radius of a sphere enclosing the input ``mass_fraction``.

        The function iterates over the radial coordinate until it finds an r=r0 where the enclosed mass fraction exceeds the input ``mass_fraction``. 
        Afterwards, it invokes spherical(r_min = r0, r_max = r_max, n_dummy = n_dummy), see `spherical`.
        
        Parameters
        ----------
        mass : array_like
           `list` or `numpy.ndarray` with the mass of each cell, where the i-th mass corresponds to the i-th cell of the GRID.

        prop : dict
           Dictionary of arrays of physical properties to be filled up by dummy values. 

        prop_fill : dict, optional
           Dictionary specifying the dummy value (scalar) with which each physical property will be filled up.
           
           Defaults to {}. In that case, the filling dummy value is zero for all the properties. 
        
        mass_fraction : float, optional
           The mass fraction enclosed by the inner radius of the spherical section. \n
           - Note that this sets the inner radius ``r_min`` at `spherical`.

           Defaults to 0.5, i.e, by default the computed inner radius will enclose (approx.) the 50% of the total mass.
        
        r_max : float, optional
           Outer radius of the spherical section enclosing the random dummy points.

           Defaults to `None`. In that case it takes the distance of the farthest point in the grid.

        n_dummy : int, optional
           Number of dummy points to be generated.

           Defaults to `None`. In that case n_dummy = GRID.NPoints / 100

        r_steps : int, optional
           Number of divisions along the radial coordinate to compute the enclosed mass.
           
           Defaults to 100.

        Returns
        -------
        out : dict
        
        Returns a dictionary with the following keys:
        
        r_rand : `numpy.ndarray`
           Radial coordinate of the computed random points.
        
        r_min : float
           Computed inner radius of the spherical section which encloses (approximately) the input ``mass fraction`` and above which the random`dummy points are generated.

        r_max : float
           Outer radius of the spherical section enclosing the random dummy points.

        comp_fraction : float
           Computed mass fraction enclosed by r_min. The higher r_steps the closer this value will be to the input ``mass fraction``. 
           
        n_dummy : int
           Number of dummy points generated.
        
        See Also
        --------
        by_density, spherical
        """
        if r_max == None: r_max = self._r_max

        mass = np.asarray(mass)
        total_mass = np.sum(mass)
        min_mass = mass_fraction * total_mass
        enc_mass = 0
        r_min = None
        for r in np.linspace(0,r_max,r_steps):
            inds = np.where(self._r_grid < r)
            enc_mass = np.sum(mass[inds])
            if enc_mass > min_mass: 
                r_min = r
                break
        print("Inner radius:", r_min)
        print("Outer radius:", r_max)
        comp_fraction = enc_mass / total_mass
        rand_spher = self.spherical(prop, prop_fill=prop_fill, r_min=r_min, r_max=r_max, n_dummy=n_dummy)
        return {"r_rand": rand_spher["r_rand"], "r_min": r_min, "r_max": r_max, 
                "comp_fraction": comp_fraction, "n_dummy": rand_spher["n_dummy"]}



