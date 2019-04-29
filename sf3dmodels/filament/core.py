from __future__ import print_function
from ..Model import Struct
from ..grid import RandomGridAroundAxis
from ..utils.units import au, pc
import numpy as np
import inspect

class DefaultFunctions(object):
    """
    Default functions for filament models.
    
    :math:`(R,\\theta,z)` are cylindrical coordinates referred to the model long axis.
    """
    """
    def _func_width(self,z,*pars):
        half_width = pars
        return half_width + np.zeros(shape=np.shape(z))
    @property
    def func_width(self): return self._func_width
    @func_width.setter
    def func_width(self, func): self._func_width = func
    """
    def func_width(self,z,*pars):
        """
        Default half-width function.
        
        Parameters
        ----------
        z : array_like, shape(n,)
           Array of z's on the model long axis.

        *pars : scalar(s)
           Function parameters. Default expected parameters: ``w0`` 

        Returns
        -------
        half_width : array_like, shape(n,)
           Half-width at each z. 
        
        Notes
        -----
        Default model: constant width
        
        .. math:: w(z) = w0
        """
        w0 = pars
        return w0 + np.zeros(shape=np.shape(z))

    def func_dens(self,R,theta,z,*pars):
        """
        Default temperature function.
        
        Parameters
        ----------
        R : array_like, shape(n,)
           Array of R's with respect to the model long axis.

        theta : array_like, shape(n,)
           Array of theta's with respect to the model long axis.

        z : array_like, shape(n,)
           Array of z's with respect to the model long axis.

        *pars : scalar(s)
           Function parameters. Default expected parameters: ``nc``, ``Rflat``, ``p`` 

        Returns
        -------
        density : array_like, shape(n,)
           Number density at each (R,theta,z) point.
        
        Notes
        -----
        Default model: Plummer-like function. See the section 2.4 of `Smith+2014b`_
        
        .. math:: n(R,\\theta,z) = \\frac{n_c}{\\big[1+(R/R_{\\rm flat})^2]^{\\rm p/2}}
        """
        nc, Rflat, p = pars
        return nc/(1.+(R/Rflat)**2)**(0.5*p)     

    def func_temp(self,R,theta,z,*pars):
        """
        Default number density function.
        
        Parameters
        ----------
        R : array_like, shape(n,)
           Array of R's with respect to the model long axis.

        theta : array_like, shape(n,)
           Array of theta's with respect to the model long axis.

        z : array_like, shape(n,)
           Array of z's with respect to the model long axis.

        *pars : scalar(s)
           Function parameters. Default expected parameters: ``T0``, ``R0``, ``p`` 

        Returns
        -------
        density : array_like, shape(n,)
            Temperature at each (R,theta,z) point.
        
        Notes
        -----
        Default model: Power-law on R.
        
        .. math:: T(R,\\theta,z) = T_0(R/R_0)^{\\rm p}
        """
        T0, R0, p = pars 
        ct = T0*R0**-p
        return ct*R**p

    def func_abund(self,R,theta,z,*pars):
        """
        Default abundance function.
        
        Parameters
        ----------
        R : array_like, shape(n,)
           Array of R's with respect to the model long axis.

        theta : array_like, shape(n,)
           Array of theta's with respect to the model long axis.

        z : array_like, shape(n,)
           Array of z's with respect to the model long axis.

        *pars : scalar(s)
           Function parameters. Default expected parameters: ``a0`` 

        Returns
        -------
        abundance : array_like, shape(n,)
           Abundance at each (R,theta,z) point.
        
        Notes
        -----
        Default model: constant abundance
        
        .. math:: a(R,\\theta,z) = a_0
        """    
        a0 = pars
        return a0 + np.zeros(shape=np.shape(z))

    def func_gtdratio(self,R,theta,z,*pars):
        """
        Default gas-to-dust ratio function.
        
        Parameters
        ----------
        R : array_like, shape(n,)
           Array of R's with respect to the model long axis.

        theta : array_like, shape(n,)
           Array of theta's with respect to the model long axis.

        z : array_like, shape(n,)
           Array of z's with respect to the model long axis.

        *pars : scalar(s)
           Function parameters. Default expected parameters: ``gtd0`` 

        Returns
        -------
        gtdratio : array_like, shape(n,)
            gas-to-dust ratio at each (R,theta,z) point.
        
        Notes
        -----
        Default model: constant gas-to-dust ratio
        
        .. math:: {\\rm gtd}(R,\\theta,z) = {\\rm gtd}0
        """       
        gtd0 = pars
        return gtd0 + np.zeros(shape=np.shape(z))

    def func_vel(self,
                 R, R_dir,
                 theta, theta_dir,
                 z, z_dir, 
                 *pars):
        """
        Default velocity function.
        
        Parameters
        ----------
        R : array_like, shape(n,1)
           Array of R's with respect to the model long axis.

        R_dir : array_like, shape(n,3)
           Unit radial vectors with respect to the model long axis.

        theta : array_like, shape(n,1)
           Array of theta's with respect to the model long axis.

        theta_dir : array_like, shape(n,3)
           Unit angular vectors with respect to the model long axis.

        z : array_like, shape(n,1)
           Array of z's with respect to the model long axis.

        z_dir : array_like, shape(n,3)
           Unit vectors parallel to the model long axis.

        *pars : scalar(s)
           Function parameters. Default expected parameters: ``(vR0,R0,p), vth0, vz0`` 

        Returns
        -------
        gtdratio : array_like, shape(n,)
            gas-to-dust ratio at each (R,theta,z) point.
        
        Notes
        -----
        Default model: (in :math:`\\hat{R}`) infalling material as a power-law; (in :math:`\\hat{\\theta}`) rotating at a constant rate; (in :math:`\\hat{z}`) flowing uniformly towards the centre along the long axis.
        
        .. math:: {\\vec{v}}(R,\\theta,z) = -v_{\\rm \\small{R_0}}(R/R_0)^{\\rm p}\\hat{R} + v_{\\theta_0}\\hat{\\theta} - v_{\\rm z_0}\\hat{z}
        """
        print(R.shape, z_dir.shape, R_dir.shape)
        (vR0,R0,p), vth0, vz0 = pars
        cvR = vR0*R0**-p
        vR = cvR*R**p 
        vth = vth0 + np.zeros(shape=np.shape(theta))
        vz = vz0 + np.zeros(shape=np.shape(z))
        return vR * (-R_dir) + vth * theta_dir + vz * (-np.sign(z)*z_dir)


class FilamentModel(RandomGridAroundAxis, DefaultFunctions):

    def __init__(self, pos_c, axis, z_min, z_max, dx, mirror=False):
        """
        Host class for filament models. The grid points are generated randomly taking into account the filament geometry.

        Available models:

           - `Smith+2014b`_: See the model functions on `DefaultFunctions`

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
           If True, it is assumed that the model is symmetric to the reference position ``pos_c``. Defaults to False.
        """
        RandomGridAroundAxis.__init__(self, pos_c, axis, z_min, z_max, dx, mirror=mirror)
        print ('Invoked %s'%self._get_classname())

    @classmethod
    def _get_classname(cls):
        return cls.__name__

    def cylinder(self, width_pars, R_min,
                 dens_pars=[1e10,0.1*pc,1.0], 
                 temp_pars=[100.,0.01*pc,-0.5], 
                 vel_pars=[(1e3,0.01*pc,-0.5), 0.5e3, 0.8e3], 
                 abund_pars=1e-8, gtdratio_pars=100., 
                 vsys=[0,0,0]):

        """
        Filament Model from `Smith+2014b`_.
        
           - See the **Section 2.4** and the **Tables 2,5 and 6** on `Smith+2014b`_ for examples on the combination of parameters.
           - See the model equations and a sketch of the cylinder geometry in the **Notes** section below.        
           - See the model functions on `DefaultFunctions`.
           - It is possible to re-define the model functions in order to modify the default filament geometry, density, etc. See the **example** section below.
           - The resulting **Attributes** depend on which parameters were set on.

        Input and output units: SI (metres, kilograms, seconds).

        Parameters
        ----------
        width_pars : scalar or array_like, shape (npars,) 
           Parameters to compute the filament half-width at each :math:`z`. \n
           Expecting: :math:`w_0`

        R_min : scalar,
           Minimum radius from the long axis to compute the physical properties.\n 
           Must be different to zero to avoid divergences in the model functions.
        
        dens_pars : scalar or array_like, shape (npars,) 
           Parameters to compute the filament number density. \n
           Expecting: [:math:`n_c`, :math:`R_{\\rm flat}`, :math:`p`]\n
           Default values: [1e10, 0.1pc, 1.0]

        temp_pars : scalar or array_like, shape (npars,) 
           Parameters to compute the filament temperature. \n
           Expecting: [:math:`T_0`, :math:`R_0`, :math:`p`]\n
           Default values: [100., 0.01pc, -0.5]

        vel_pars : scalar or array_like, shape (npars,) 
           Parameters to compute the filament temperature. \n
           Expecting: [(:math:`v_{\\small{\\rm R_0}}`, :math:`R_0`, :math:`p`), :math:`v_{\\theta_0}`, :math:`v_{z_0}`]\n
           Default values: [(1000., 0.01pc, -0.5), 500., 800.]

        abund_pars : scalar or array_like, shape (npars,) 
           Parameters to compute the filament molecular abundance. \n
           Expecting: :math:`a_0`\n
           Default values: 1e-8

        gtdratio_pars : scalar or array_like, shape (npars,) 
           Parameters to compute the filament gas-to-dust ratio. \n
           Expecting: :math:`\\rm{gtd}0`\n
           Default values: 100.

        vsys : array_like, shape (3,)
           [:math:`v_{0x}, v_{0y}, v_{0z}`]: Systemic velocity to be added to the final cylinder velocity field. \n
           Defaults to [0,0,0].

        Attributes
        ----------
        density : `numpy.ndarray`, shape (n,) 
           Hydrogen number density.

        temperature : `numpy.ndarray`, shape (n,) 
           Jet temperature.

        vel : `~sf3dmodels.Model.Struct`, 
           Jet velocity. Attributes: vel.x, vel.y, vel.z, shape (n, ) each.

        abundance : `numpy.ndarray`, shape (n,) 
           Molecular abundance.

        gtdratio : `numpy.ndarray`, shape (n,)
           Uniform gas-to-dust ratio.

        GRID : `~sf3dmodels.Model.Struct`, 
           Grid point coordinates and number of grid points. \n
           Attributes: GRID.XYZ, shape(3,n): [x,y,z]; GRID.NPoints, scalar.

        r : `numpy.ndarray`, shape (n,) 
           Grid points distance to the outflow centre.

        Raises
        ------
        ValueError : If R_min == 0.0

        Examples
        --------
        Have a look at the `examples/outflows/ <https://github.com/andizq/star-forming-regions/tree/master/examples/outflows>`_ folder on the GitHub repository. 
        
        The example shows how to compute the free-free emission from a couple of outflows using `RADMC-3D`_. 
        Figures below, from left to right: Outflows grid points distribution; Spectral Energy Distribution; Free-free continuum image at :math:`\\lambda=` 1000 microns. 

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

        if R_min == 0.0: raise ValueError('R_min must be different to zero to avoid divergences in the model functions') 

        pars = {'width': width_pars, 
                'dens': dens_pars, 
                'temp': temp_pars, 
                'vel': vel_pars, 
                'abund': abund_pars, 
                'gtdratio': gtdratio_pars}

        for key in pars:
            if isinstance(pars[key], (float,int,long,complex)):
                pars[key] = [pars[key]]

        
        self._grid(self.func_width, pars['width'], R_min = R_min)
        self.GRID = Struct( **{'XYZ': self.grid.T, 'NPoints': self.NPoints})
    
        if dens_pars is not None: self.dens_gas = self.func_dens(self.R, self.theta, self.z, *pars['dens'])
        if temp_pars is not None: self.temp_gas = self.func_temp(self.R, self.theta, self.z, *pars['temp'])
        if abund_pars is not None: self.abund = self.func_abund(self.R, self.theta, self.z, *pars['abund'])
        if gtdratio_pars is not None: self.gtdratio = self.func_gtdratio(self.R, self.theta, self.z, *pars['gtdratio'])
        if vel_pars is not None: self.vel = self.func_vel(self.R[:,None], self.R_dir,
                                                          self.theta[:,None], self.theta_dir,
                                                          self.z[:,None], self.z_dir, *pars['vel'])
    
        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------\n-------------------------------------------------')
