from __future__ import print_function
from ..Model import Struct
from ..grid import RandomGridAroundAxis
from ..utils.units import au, pc
from ..utils.constants import temp_cmb
import numpy as np
import inspect

class DefaultFilamentFunctions(object):
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
        (vR0,R0,p), vth0, vz0 = pars
        cvR = vR0*R0**-p
        vR = cvR*R**p 
        vth = vth0 + np.zeros(shape=np.shape(theta))
        vz = vz0 + np.zeros(shape=np.shape(z))
        return vR * (-1*R_dir) + vth * theta_dir + vz * (-np.sign(z)*z_dir)


class FilamentModel(RandomGridAroundAxis, DefaultFilamentFunctions):

    def __init__(self, pos_c, axis, z_min, z_max, dx, mirror=False):
        """
        Host class for filament models. 

        The grid points are generated randomly taking into account the filament width as a function of the Z-coordinate.   

        Available models:

           - `Smith+2014b`_: See the model functions on `DefaultFilamentFunctions`

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
           Axis upper limit to compute the physical properties. The axis extent will be [``z_min``, ``z_max``].

        dx : scalar
           Maximum separation between two adjacent grid points. This will prevent void holes when merging this random grid into a regular grid of nodes separation = ``dx``.       

        mirror : bool, optional
           If True, it is assumed that the model is symmetric to the reference position ``pos_c``. Defaults to False.
        
        Notes
        -----
        The reference axis to compute :math:`\\theta` is gotten as follows::
        
           >>> cart = np.zeros(3) #Initializing an empty cartesian vector.
           >>> cart[np.argmin(axis)] = 1. #Make 1 the component where the axis vector is shorter. 
           >>> axis_theta = np.cross(axis, cart) #The reference axis for theta is the cross product of axis and cart.
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
                 vsys=[0,0,0], dummy_frac=0.0, 
                 dummy_values = {'density': 1e3,
                                 'temperature': temp_cmb,
                                 'abundance': 1e-12,
                                 'gtdratio': 100.,
                                 'vx': 0.0,
                                 'vy': 0.0,
                                 'vz': 0.0} ):

        """
        Filament Model from `Smith+2014b`_.
        
           - See the **Section 2.4** and the **Tables 2,5 and 6** on `Smith+2014b`_ for examples on the combination of parameters.
           - See the model equations and a sketch of the cylinder geometry in the **Notes** section below.        
           - See the model functions on `DefaultFilamentFunctions`.
           - It is possible to customise the model functions in order to control the filament geometry, density, etc. See the **example** section below.
           - The resulting **Attributes** depend on which parameters were set different to None.

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

        dummy_frac : scalar
           Fraction of additional points (dummy points) to fill the grid borders.
           For radiative transfer purposes it is strongly recommended to set on this parameter (to at least 0.3)
           in case the width function of the filament is not constant, 
           so that the RT code(s) can know explicitly where the empty spaces are and does not perform 
           any kind of emission interpolation.           
           Defaults to 0.0\n
           ndummies = dummy_frac*npoints
        
        dummy_values : dict
           Dictionary containing the values of the properties of the dummy points. It will only be used if dummy_frac > 0.0

        Warnings
        --------
        For radiative transfer purposes, if you modify the width function `~DefaultFilamentFunctions.func_width` 
        it is strongly recommended to set on the ``dummy_frac`` parameter (to at least 0.3), 
        so that the RT code(s) can know explicitily where the empty spaces are and does 
        not perform any kind of emission interpolation.\n 

        Attributes
        ----------
        density : `numpy.ndarray`, shape (n,) 
           Gas number density.

        temperature : `numpy.ndarray`, shape (n,) 
           Gas temperature.

        vel : `~sf3dmodels.Model.Struct`, 
           The velocity field. Attributes: vel.x, vel.y, vel.z, shape (n, ) each.

        abundance : `numpy.ndarray`, shape (n,) 
           Molecular abundance.

        gtdratio : `numpy.ndarray`, shape (n,)
           Gas-to-dust ratio.

        GRID : `~sf3dmodels.Model.Struct`, 
           Coordinates of the grid points and number of grid points. \n
           Attributes: GRID.XYZ, shape(3,n): [x,y,z]; GRID.NPoints, scalar.

        R : `numpy.ndarray`, shape (n,) 
           Radial distance of the grid points, referred to the filament frame of reference.

        theta : `numpy.ndarray`, shape (n,) 
           Azimuthal coordinate of the grid points, referred to the filament frame of reference.

        z : `numpy.ndarray`, shape (n,) 
           Height of the grid points, referred to the filament frame of reference.
        
        R_dir : `numpy.ndarray`, shape (n,3) 
           Radial unit vector of the grid points, referred to the filament frame of reference.

        theta_dir : `numpy.ndarray`, shape (n,3) 
           Azimuthal unit vector of the grid points, referred to the filament frame of reference.

        z_dir : `numpy.ndarray`, shape (n,3) 
           Z unit vector of the grid points, referred to the filament frame of reference.

        Raises
        ------
        ValueError : If R_min == 0.0
        
        Notes
        -----
        Default model equations: `DefaultFilamentFunctions`
                
        .. image:: https://github.com/andizq/andizq.github.io/blob/master/star-forming-regions/examples/filament/filament_edgeon.png?raw=true
           :width: 320px
           :height: 320px
           :alt: Filament sketch - edgeon

        .. image:: https://github.com/andizq/andizq.github.io/blob/master/star-forming-regions/examples/filament/filament_faceon.png?raw=true
           :width: 320px
           :height: 270px
           :alt: Filament sketch - faceon

        **Figure 1**. Sketches showing the basic input geometrical components that shape the filament. Left: edge-on view, Right: face-on view. The solid vectors and lengths are referred to the **filament** frame of reference and the dashed vectors to the **GRID** frame of reference.

        Examples
        --------

        Let's get started with the simplest model.
        The default filament is used except for the temperature and abundance parameters. 
        It's centred at (0,0,0); pointing to :math:`z`, i.e (0,0,1); with 0.4 pc of length and 0.1 pc of width.  

        .. plot:: 
           :include-source: True

           import numpy as np
           import sf3dmodels.filament as sf
           import sf3dmodels.Plot_model as pm
           from sf3dmodels.utils.units import pc

           f1 = sf.FilamentModel([0,0,0], [0,0,1], -0.2*pc, 0.2*pc, 0.01*pc) # Invoke the class with the filament axis geometry
           f1.cylinder(0.1*pc, 1e-3*pc, temp_pars = [500, 0.02*pc, -0.3], abund_pars = 1e-4) # Specify the method and the physical parameters
           
           pm.scatter3D(f1.GRID, f1.density, np.mean(f1.density), axisunit = pc,
                        colordim = f1.temperature, 
                        colorlabel = 'T [K]',
                        NRand = 10000, show=True)
                        
        **Figure 2**. The plot shows 10000 grid points, randomly picked based on their density. The colours represent the temperature of the grid points. 

        For the radiative transfer we need to construct the ``prop`` dictionary and write the necessary files depending on the code that will be used (see `~sf3dmodels.rt.Lime` or `~sf3dmodels.rt.Radmc3d`):

        .. code-block:: python
           
           import sf3dmodels.rt as rt
           prop = {'dens_H': f1.density,
                   'temp_gas': f1.temperature,
                   'abundance': f1.abundance,
                   'gtdratio': f1.gtdratio,
                   'vel_x': f1.vel.x,
                   'vel_y': f1.vel.y,
                   'vel_z': f1.vel.z,
                   }

           lime = rt.Lime(f1.GRID)
           lime.submodel(prop, output='datatab.dat', folder='./', 
                         lime_header=True, lime_npoints=True)

        Since the grid of our filament is non-regular the command to run Lime in sf3dmodels mode will have an additional -G flag: 
        
        .. code-block:: bash
        
           $ lime -nSG -p 4 model.c
        
        It will generate a set of cubes containing the :math:`^{12}\\rm CO` J=1-0 line emission and the dust emission for this model. 

        .. image:: https://github.com/andizq/andizq.github.io/blob/master/star-forming-regions/examples/filament/filament_limedust_edgeon.png?raw=true
           :width: 335px
           :height: 275px
           :alt: Filament dust emission - edgeon

        .. image:: https://github.com/andizq/andizq.github.io/blob/master/star-forming-regions/examples/filament/filament_limedust_faceon.png?raw=true
           :width: 335px
           :height: 275px
           :alt: Filament dust emission - faceon
        
        **Figure 3**. 115 GHz dust emission from the default cylindrical filament for two different orientations.
        
        Let's now customise some model functions and compute again the radiative transfer. 
        This time the **temperature** will also be a function of :math:`z`, the **abundance** is no longer constant, 
        and the **width** is a sine function of :math:`z`.
        
        .. plot::
           :include-source: True
           
           import numpy as np
           import sf3dmodels.filament as sf
           import sf3dmodels.Plot_model as pm
           from sf3dmodels.utils.units import pc


           f1 = sf.FilamentModel([0,0,0], [0,0,1], -0.2*pc, 0.2*pc, 8e-3*pc)
           
           def new_width(z, *width_pars):
               w0, period = width_pars
               return w0*((0.5*np.sin(z*2*np.pi/period)**2) + 0.5)
    
           def new_abund(R,theta,z, *abund_pars):
               a0, R0, p = abund_pars
               ca = a0*R0**-p
               return ca*R**p

           def new_temp(R,theta,z, *temp_pars):
               TR, R0, pR, zh = temp_pars
               cR = TR*R0**-pR
               return cR*R**pR * np.exp(np.abs(z)/zh)
    
           f1.func_width = new_width
           f1.func_abund = new_abund
           f1.func_temp = new_temp

           f1.cylinder([0.1*pc, 0.3*pc], 1e-4*pc, 
                       abund_pars = [1e-5, 0.05*pc, -0.15],
                       temp_pars = [200, 0.02*pc, -0.15, -0.17*pc],
                       dummy_frac = 0.5)

           pm.scatter3D(f1.GRID, f1.density, np.mean(f1.density), axisunit = pc,
                        colordim = f1.temperature,
                        colorlabel = 'T [K]',
                        NRand = 10000, 
                        cmap = 'nipy_spectral_r',
                        azim=45, elev=15, show=True)
        
           pm.scatter3D(f1.GRID, f1.density, np.min(f1.density), axisunit = pc,
                        colordim = f1.abundance,
                        colorlabel = 'Molec. abund.',
                        NRand = 10000, 
                        cmap = 'nipy_spectral_r',
                        azim=45, elev=15, show=True)

        **Figure 4**. 10000 grid points, randomly chosen according to their density. 
        In the top plot the colours represent the temperature of the grid points and in the bottom one the abundance.
        For the bottom plot the weighting value was reduced in order to make the dummy points visible (in gray colour), 
        which were activated here for radiative transfer purposes.

        And the block concerning to the radiative transfer:

        .. code-block:: python
           
           import sf3dmodels.rt as rt
           prop = {'dens_H': f1.density,
                   'temp_gas': f1.temperature,
                   'abundance': f1.abundance,
                   'gtdratio': f1.gtdratio,
                   'vel_x': f1.vel.x,
                   'vel_y': f1.vel.y,
                   'vel_z': f1.vel.z,
                   }

           lime = rt.Lime(f1.GRID)
           lime.submodel(prop, output='datatab.dat', folder='./', 
                         lime_header=True, lime_npoints=True)

        .. code-block:: bash
        
           $ lime -nSG -p 4 model.c

        The dust emission: 

        .. image:: https://github.com/andizq/andizq.github.io/blob/master/star-forming-regions/examples/filament/filamentcustom_limedust_edgeon.png?raw=true
           :width: 335px
           :height: 275px
           :alt: Filament dust emission - edgeon

        .. image:: https://github.com/andizq/andizq.github.io/blob/master/star-forming-regions/examples/filament/filamentcustom_limedust_faceon.png?raw=true
           :width: 335px
           :height: 275px
           :alt: Filament dust emission - faceon
        
        **Figure 5**. 115 GHz dust emission from the customised filament for two different orientations.

        Finally let's compute the first two moment maps of :math:`^{12}\\rm CO J=1-0` from the Lime's output::
        
           >>> from astropy.io import fits
           >>> from astropy.wcs import WCS
           >>> from spectral_cube import SpectralCube

           >>> fn = 'img_custom_CO_J1-0_LTE_jypxl_edgeon.fits'
           >>> hdu = fits.open(fn)[0]
           >>> hdu.header['CUNIT3'] = 'm/s' 

           >>> w = WCS(hdu.header)
           >>> cube = SpectralCube(data=hdu.data.squeeze(), wcs=w.dropaxis(3))
           
           >>> m0 = cube.moment(order=0) # Moment 0 (Jy/pxl m/s)
           >>> m1 = cube.moment(order=1) # Moment 1 (m/s)
           >>> m0.write('moment0_'+fn, overwrite=True)
           >>> m1.write('moment1_'+fn, overwrite=True)

           >>> fig, ax = plt.subplots(ncols=2, subplot_kw = {'projection': w.celestial}, figsize = (12,6))
           
           >>> im0 = ax[0].imshow(m0.array, cmap='hot')
           >>> im1 = ax[1].imshow(m1.array, cmap='nipy_spectral')           
           >>> cb0 = fig.colorbar(im0, ax=ax[0], orientation='horizontal', format='%.2f', pad=0.1)
           >>> cb1 = fig.colorbar(im1, ax=ax[1], orientation='horizontal', pad=0.1)
           >>> ax[0].set_title('Moment 0 - edgeon', fontsize = 14)
           >>> ax[1].set_title('Moment 1 - edgeon', fontsize = 14)
           >>> cb0.set_label(r'Jy pxl$^{-1}$ m s$^{-1}$', fontsize = 14)
           >>> cb1.set_label(r'm s$^{-1}$', fontsize = 14)
           
           >>> plt.tight_layout(pad=2.0)
           >>> fig.savefig('customfilament_moments_edgeon.png')
           
        .. image:: https://github.com/andizq/andizq.github.io/blob/master/star-forming-regions/examples/filament/customfilament_moments_edgeon.png?raw=true
           :width: 770px
           :height: 350px
           :alt: Filament moment maps - edgeon

        Find the source codes on LINKTOFOLDER.

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

        self._grid(self.func_width, pars['width'], R_min=R_min, dummy_frac=dummy_frac)
        if self.ndummies > 0: self.GRID = Struct( **{'XYZ': np.append(self.grid, self.grid_dummy, axis=0).T, 'NPoints': self.NPoints})
        else: self.GRID = Struct( **{'XYZ': self.grid.T, 'NPoints': self.NPoints})

        append_dummies = lambda prop, tag, n: np.append(prop, np.zeros(n)+dummy_values[tag]) 

        if dens_pars is not None: 
            self.density = self.func_dens(self.R, self.theta, self.z, *pars['dens'])
            if self.ndummies > 0: self.density = append_dummies(self.density, 'density', self.ndummies)
        if temp_pars is not None: 
            self.temperature = self.func_temp(self.R, self.theta, self.z, *pars['temp'])
            if self.ndummies > 0: self.temperature = append_dummies(self.temperature, 'temperature', self.ndummies)
        if abund_pars is not None: 
            self.abundance = self.func_abund(self.R, self.theta, self.z, *pars['abund'])
            if self.ndummies > 0: self.abundance = append_dummies(self.abundance, 'abundance', self.ndummies)
        if gtdratio_pars is not None: 
            self.gtdratio = self.func_gtdratio(self.R, self.theta, self.z, *pars['gtdratio'])
            if self.ndummies > 0: self.gtdratio = append_dummies(self.gtdratio, 'gtdratio', self.ndummies)
        if vel_pars is not None: 
            vx, vy, vz = self.func_vel(self.R[:,None], self.R_dir,
                                       self.theta[:,None], self.theta_dir,
                                       self.z[:,None], self.z_dir, *pars['vel']).T
            if self.ndummies > 0: 
                vx = append_dummies(vx, 'vx', self.ndummies)
                vy = append_dummies(vy, 'vy', self.ndummies)
                vz = append_dummies(vz, 'vz', self.ndummies)
            self.vel = Struct( **{'x': vx, 'y': vy, 'z': vz})
    
        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------\n-------------------------------------------------')
