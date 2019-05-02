from __future__ import print_function
from ..Model import Struct
import numpy as np
import inspect

class WriteProp(object):
    def __init__(self, prop):
        self._write(prop)
    def _write(self, prop):
        pass
        #print (prop)

class RandomGrid(object):
    """
    Base class for Random grids
    """
    def __init__(self, pos_c, axis, z_min, z_max, dx):
        self.pos_c = np.asarray(pos_c)
        self.axis = np.asarray(axis)
        self.pos_f = self.pos_c + self.axis
        self.z_min = z_min
        self.z_max = z_max
        self.dx = dx        
        
    def _axis(self):
        r_seg = self.axis
        r_seg_mag = np.linalg.norm(r_seg) #Outflow axis magnitude (preliminary)

        self.r_seg_dir = r_seg/r_seg_mag #Axis direction
        if self.z_max is not None: self.r_seg = (self.z_max - self.z_min) * self.r_seg_dir #If set z_max, then the segment length is zmax-zmin  
        else: self.r_seg = r_seg - self.z_min * self.r_seg_dir #New axis, only taking into account z_min
        self.r_seg_mag = np.linalg.norm(self.r_seg) #Final axis magnitude
        self.dr = self.dx * 2.**-1 #3.**-1#2.5**-1
        return self.r_seg_dir

    def _grid(self, func_width):
        
        #Guess the number of random grid points to generate: (approximation to a rectangular region)
        # Number of divisions along the main segment, 
        #  times the number of divisions along an axis perpendicular to the segment,
        #   times the number of divisions along an axis perpendicular to both of the segments above.
        #    times an exageration factor (in dr) to increase the number of grid points.
        #     Twice to accoutn for the other branch of the jet.

        mean_w = func_width(0.5 * self.r_seg_mag)
        self.NPoints = 2 * int(self.r_seg_mag/self.dr * (mean_w/self.dr)**2 ) 
        print ('Number of grid points:', self.NPoints)
        
        self.grid = np.zeros((self.NPoints, 4)) #x,y,z,r
        half_points = int(self.NPoints/2)
        r_vec = np.zeros((half_points, 3))
        rand_vec_plane = np.zeros((half_points, 3))
                
        r = np.random.uniform(self.z_min, self.z_min+self.r_seg_mag, size = half_points) 
        width = func_width(r)
        for i in range(half_points):
            r_vec[i] = r[i] * self.r_seg_dir #Random vector along the long axis
            R = np.random.uniform(0, width[i]) #Random R from the long axis

            rand_vec = np.random.uniform(-1,1,3) #Random vector to generate the perpendicular vector to the outflow axis
            rand_vec = rand_vec / np.linalg.norm(rand_vec)
            cross_unit = np.cross(self.r_seg_dir, rand_vec)
            cross_unit = cross_unit / np.linalg.norm(cross_unit) #Perpendicular (random) unitary vector to the outflow axis
            rand_vec_plane[i] = R * cross_unit #Perpendicular (random) vector to the outflow axis

        r_c = r_vec + rand_vec_plane #Vector from the outflow center to the generated point
        r_real = self.pos_c + r_c #Real position from the origin of coordinates
        r_real_n = self.pos_c - r_c #Mirror point to real position from the origin of coordinates

        r = np.expand_dims(r, axis=0) # or r[np.newaxis]
        self.grid[0:half_points] = np.append(r_real, r.T, axis=1)
        self.grid[half_points:] = np.append(r_real_n, r.T, axis=1)
        r_c_dir = r_c / np.linalg.norm(r_c, axis = 1)[np.newaxis].T
        self.r_c_dir = np.append(r_c_dir, -1*r_c_dir, axis = 0) 
        """
        self.r_c_dir = np.append(np.repeat([self.r_seg_dir], half_points, axis=0),
                                 np.repeat([-1*self.r_seg_dir], half_points, axis=0),
                                 axis = 0)
        """
        self.width = np.append(width, width)
                
class OutflowModel(RandomGrid):
   
    def __init__(self, pos_c, axis, z_min, z_max, dx):
        """
        Host class for outflow models. The grid points are generated randomly taking into account the Jet model geometry.

        Available models:

           - 'reynolds86': Jet Model (`Reynolds+1986`_)

        Input units: SI (metres, kilograms, seconds).

        Parameters
        ----------
        pos_c : array_like, shape (3,) 
           Origin of coordinates on the axis.

        axis : array_like, shape (3,) 
           Long axis direction, of arbitrary length. For example, an axis pointing to the z-direction: axis = [0,0,1]  

        z_min : scalar
           Axis lower limit to compute the physical properties.

        z_max : scalar
           Axis upper limit to compute the physical properties. The axis extent will be [z_min, z_max].

        dx : scalar
           Maximum separation between two adjacent grid points. This will prevent void holes when merging this random grid into a regular grid of nodes separation = dx.       

        mirror : bool, optional
           If True, it is assumed that the model is symmetric to the reference position ``pos_c``. Defaults to True.
        """
        RandomGrid.__init__(self, pos_c, axis, z_min, z_max, dx)
        print ('Invoked %s'%self._get_classname())

    @classmethod
    def _get_classname(cls):
        return cls.__name__
    
    def reynolds86(self, width_pars, dens_pars=None, ionfrac_pars=None, temp_pars=None, z0=None, v0=None, abund_pars=None, gtdratio=None, vsys=[0,0,0]):
        """
        Outflow Model from `Reynolds+1986`_.
        
           - See the **Table 1** on `Reynolds+1986`_ for examples on the combination of parameters and their physical interpretation.
           - See the model equations and a sketch of the jet geometry in the **Notes** section below.
           - The resulting **Attributes** depend on which parameters were set different to None.
        
        Input and output units: SI (metres, kilograms, seconds).

        Parameters
        ----------
        width_pars : array_like, shape (2,) 
           [:math:`w_0`, :math:`\\epsilon`]: parameters to compute the Jet half-width at each :math:`z`.

        dens_pars : array_like, shape (2,) 
           [:math:`n_0`, :math:`q_n`]: parameters to compute the Jet number density.

        ionfrac_pars : array_like, shape (2,) 
           [:math:`x_0`, :math:`q_x`]: parameters to compute the ionized fraction of gas.

        temp_pars : array_like, shape (2,) 
           [:math:`T_0`, :math:`q_T`]: parameters to compute the Jet temperature.

        abund_pars : array_like, shape (2,) 
           [:math:`a_0`, :math:`q_a`]: parameters to compute the molecular abundance.

        z0 : scalar
           Normalization distance. Must be different to 0 (zero) to avoid divergences. \n
           If None, ``z_max/100`` is used. \n
           Defaults to None.

        v0 : scalar
           Gas speed at `z0`. Speed normalization factor. For mass conservation: :math:`q_v=-q_n-2\\epsilon`. 

        gtdratio : scalar
           Gas-to-dust ratio value (mass ratio). Uniform at every point of the grid.
        
        vsys : array_like, shape (3,)
           [:math:`v_{0x}, v_{0y}, v_{0z}`]: Systemic velocity to be added to the final Jet velocity field. \n
           Defaults to [0,0,0].

        Attributes
        ----------
        density : `numpy.ndarray`, shape (n,) 
           Atomic hydrogen number density.

        temperature : `numpy.ndarray`, shape (n,) 
           Gas temperature.
        
        ionfraction : `numpy.ndarray`, shape (n,) 
           Ionized fraction of hydrogen.

        density_ion : `numpy.ndarray`, shape (n,) 
           Ionized hydrogen number density: density * ionfraction.

        speed : `numpy.ndarray`, shape (n,) 
           Speed.

        vel : `~sf3dmodels.Model.Struct`, 
           The velocity field. Attributes: vel.x, vel.y, vel.z, shape (n, ) each.

        abundance : `numpy.ndarray`, shape (n,) 
           Molecular abundance.

        gtdratio : `numpy.ndarray`, shape (n,)
           Gas-to-dust ratio.

        GRID : `~sf3dmodels.Model.Struct`, 
           Grid point coordinates and number of grid points. \n
           Attributes: GRID.XYZ, shape(3,n): [x,y,z]; GRID.NPoints, scalar.

        r : `numpy.ndarray`, shape (n,) 
           Grid points distance to the outflow centre.

        Examples
        --------
        Have a look at the `examples/outflows/ <https://github.com/andizq/star-forming-regions/tree/master/examples/outflows>`_ folder on the GitHub repository. 
        
        The example shows how to compute the free-free emission from a couple of outflows using `RADMC-3D`_. 
        Figures below, from left to right: Outflows distribution of grid points; Spectral Energy Distribution; Free-free continuum image at :math:`\\lambda=` 1000 microns. 

        .. image:: ../../examples/outflows/global_grid_dens.png
           :width: 220px
           :height: 170px

        .. image:: ../../examples/outflows/sed_outflows.png
           :width: 200px
           :height: 170px

        .. image:: ../../examples/outflows/img_outflows.png
           :width: 220px
           :height: 170px
        
        Notes
        -----
        Model equations:
        
        \t With :math:`r = \sqrt{x^2 + y^2+ z^2}`,

        .. math:: w(z) &= w_0(z/z_0)^\\epsilon \n 
                  n_{H}(z) &= n_0(z/z_0)^{q_n} \n
                  T(z) &= T_0(z/z_0)^{q_T} \n
                  x(z) &= x_0(z/z_0)^{q_x} \n
                  n_{ion}(z) &= n_{H}(z) x(z)\n
                  a(z) &= a_0(z/z_0)^{q_a} \n
                  v(z) &= v_0(z/z_0)^{q_v} \n
                  \\vec{v}(r) &= v(z) \\hat{r}, \n
        

        where :math:`q_v = -q_n-2\epsilon` for mass conservation.

        .. figure:: ../../images/reynolds86_outflow.png
           :width: 350px
           :align: center
           :height: 400px
           :alt: reynolds86 model
           :figclass: align-center

           Figure 1 from `Reynolds+1986`_. Jet geometry. Note that the :math:`r` from the sketch is :math:`z` in our equations.                           
        """
             
        if z0 is None: r0 = 1e-2*self.z_max
        else: r0 = float(z0)
        
        cw = width_pars[0] * r0**-width_pars[1]
        def func_width(r): #Jet half-width 
            return cw * r**width_pars[1]
        
        if dens_pars is not None: cd = dens_pars[0] * r0**-dens_pars[1]
        def func_density(r): #Jet density
            return cd * r**dens_pars[1]
        
        if temp_pars is not None: ct = temp_pars[0] * r0**-temp_pars[1]
        def func_temperature(r): #Jet temperature
            return ct * r**temp_pars[1]
                
        if ionfrac_pars is not None: ci = ionfrac_pars[0] * r0**-ionfrac_pars[1]
        def func_ionfraction(r): #Jet ionization fraction
            return ci * r**ionfrac_pars[1]
        
        if abund_pars is not None: ca = abund_pars[0] * r0**-abund_pars[1]
        def func_abundance(r):
            return ca * r**abund_pars[1]

        if v0 is not None and dens_pars is not None:
            qv = -2 * width_pars[1] - dens_pars[1]
            cv = v0 * r0**-qv
            print ('Velocity powerlaw, qv = -2eps - qn :', qv)
        def func_speed(r): #Jet velocity
            return cv * r**qv
        
        self._axis()
        self._grid(func_width)
        self.x, self.y, self.z, r = self.grid.T
        self.r = r
        self.GRID = Struct( **{'XYZ': [self.x,self.y,self.z], 'NPoints': self.NPoints})
        
        if v0 is not None and dens_pars is not None: 
            self.speed = func_speed(r)
            self.vx, self.vy, self.vz = (self.speed[np.newaxis].T * self.r_c_dir + np.asarray(vsys)).T
            self.vel = Struct( **{'x': self.vx, 'y': self.vy, 'z': self.vz})
        if dens_pars is not None: self.density = func_density(r)
        if temp_pars is not None: self.temperature = func_temperature(r)
        if ionfrac_pars is not None: self.ionfraction = func_ionfraction(r)
        if dens_pars is not None and ionfrac_pars is not None: self.density_ion = self.density * self.ionfraction
        if abund_pars is not None: self.abundance = func_abundance(r)
        if gtdratio is not None: self.gtdratio = gtdratio + np.zeros(self.NPoints)
        
        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------\n-------------------------------------------------')
