"""
2D Disc models
==============
Classes: Rosenfeld2d, General2d, Velocity, Intensity
"""

from ..utils.constants import G, kb
from ..utils import units as u
import numpy as np
import numbers
import warnings
from scipy.interpolate import griddata
import os

#warnings.filterwarnings("error")

class Tools:
    @staticmethod
    def _rotate_sky_plane(x, y, ang):
        xy = np.array([x,y])
        cos_ang = np.cos(ang)
        sin_ang = np.sin(ang)
        rot = np.array([[cos_ang, -sin_ang],
                        [sin_ang, cos_ang]])
        return np.dot(rot, xy)

    @staticmethod
    def _rotate_sky_plane3d(x, y, z, ang, axis='z'):
        xyz = np.array([x,y,z])
        cos_ang = np.cos(ang)
        sin_ang = np.sin(ang)
        if axis == 'x':
            rot = np.array([[1, 0, 0],
                            [0, cos_ang, -sin_ang],
                            [0, sin_ang, cos_ang]])
        if axis == 'y':
            rot = np.array([[cos_ang, 0, -sin_ang],
                            [0, 1, 0],
                            [sin_ang, 0, cos_ang]])
            
        if axis == 'z':
            rot = np.array([[cos_ang, -sin_ang , 0],
                            [sin_ang, cos_ang, 0], 
                            [0, 0, 1]])
        return np.dot(rot, xyz)

    @staticmethod
    def _project_on_skyplane(x, y, z, cos_incl, sin_incl):
        x_pro = x
        y_pro = y * cos_incl - z * sin_incl
        z_pro = y * sin_incl + z * cos_incl
        return x_pro, y_pro, z_pro

    @staticmethod
    def _compute_prop(grid, prop_funcs, prop_kwargs):
        n_funcs = len(prop_funcs)
        props = [{} for i in range(n_funcs)]
        for side in ['near', 'far']:
            x, y, z, R, phi = grid[side]
            coord = {'x': x, 'y': y, 'z': z, 'phi': phi, 'R': R}
            for i in range(n_funcs): props[i][side] = prop_funcs[i](coord, **prop_kwargs[i])
        return props

    @staticmethod
    def _line_profile(v_chan, v, T, v_turb=0, mmol=2*u.amu):
        v_sigma = np.sqrt(2*kb*T/mmol + v_turb**2) * 1e-3 #in km/s
        return 1/(np.sqrt(np.pi)*v_sigma) * np.exp(-((v-v_chan)/v_sigma)**2)

    @staticmethod
    def get_channel(velocity2d, intensity2d, temperature2d, v_chan, 
                    v_turb=0.0, mmol=2*u.amu):

        vel2d, temp2d, int2d = velocity2d, {}, {}
        if isinstance(temperature2d, numbers.Number): temp2d['near'] = temp2d['far'] = temperature2d
        else: temp2d = temperature2d

        if isinstance(intensity2d, numbers.Number): int2d['near'] = int2d['far'] = intensity2d
        else: int2d = intensity2d
    
        v_near = Tools._line_profile(v_chan, vel2d['near'], temp2d['near'], mmol=28.0101*u.amu) #12CO
        v_far = Tools._line_profile(v_chan, vel2d['far'], temp2d['far'], mmol=28.0101*u.amu) #12CO

        v_near_clean = np.where(np.isnan(v_near), -np.inf, v_near)
        v_far_clean = np.where(np.isnan(v_far), -np.inf, v_far)
        
        #int2d_near = int2d['near'] * v_near_clean / v_near_clean.max()
        int2d_near = np.where(np.isnan(int2d['near']), -np.inf, int2d['near'] * v_near_clean / v_near_clean.max())
        int2d_far = np.where(np.isnan(int2d['far']), -np.inf, int2d['far'] * v_far_clean / v_far_clean.max())
        
        vmap_full = np.array([v_near_clean, v_far_clean]).max(axis=0)
        int2d_full = np.array([int2d_near, int2d_far]).max(axis=0)

        return int2d_full

    @staticmethod
    def make_channels_movie(vchan0, vchan1, velocity2d, intensity2d, temperature2d, nchans=30, folder='./movie_channels/', **kwargs):
        import matplotlib.pyplot as plt
        vchannels = np.linspace(vchan0, vchan1, num=nchans)
        int2d_cube = []
        for i, vchan in enumerate(vchannels):
            int2d = Tools.get_channel(velocity2d, intensity2d, temperature2d, vchan, **kwargs)
            int2d_cube.append(int2d)
            extent = [-600, 600, -600, 600]
            plt.imshow(int2d, cmap='binary', extent=extent, origin='lower left', vmax=np.max(temperature2d['near']))
            plt.xlabel('au')
            plt.ylabel('au')
            plt.text(200, 500, '%.1f km/s'%vchan)
            cbar = plt.colorbar()
            cbar.set_label(r'Brightness Temperature [K]')
            plt.contour(velocity2d['near'], levels=[vchan], colors='red', linestyles='--', linewidths=1.3, extent = extent)
            plt.contour(velocity2d['far'], levels=[vchan], colors='red', linestyles=':', linewidths=1.3, extent = extent)
            plt.plot([None],[None], color='red', linestyle='--', linewidth=2, label='Near side') 
            plt.plot([None],[None], color='red', linestyle=':', linewidth=2, label='Far side') 
            plt.legend(loc='upper left')
            plt.savefig(folder+'int2d_chan%04d'%i)
            print ('Saved channel %d'%i)
            plt.close()

        os.chdir(folder)
        print ('Making channels movie...')
        os.system('convert -delay 10 *int2d* cube_channels.gif')
        return np.array(int2d_cube)

class Velocity:
    @property
    def velocity_func(self): 
        return self._velocity_func
          
    @velocity_func.setter 
    def velocity_func(self, vel): 
        print('Setting velocity function to', vel) 
        self._velocity_func = vel

    @velocity_func.deleter 
    def velocity_func(self): 
        print('Deleting velocity function') 
        del self._velocity_func

    @staticmethod
    def keplerian(coord, Mstar=u.Msun):
        if 'R' not in coord.keys(): R = np.hypot(coord['x'], coord['y'])
        else: R = coord['R'] 
        return np.sqrt(G*Mstar/R) 
    
    @staticmethod
    def keplerian_vertical(coord, Mstar=u.Msun):
        if 'R' not in coord.keys(): R = np.hypot(coord['x'], coord['y'])
        else: R = coord['R'] 
        if 'r' not in coord.keys(): r = np.hypot(R, coord['z'])
        else: r = coord['r']
        return np.sqrt(G*Mstar/r**3)*R


class Intensity:
    @property
    def intensity_func(self): 
        return self._intensity_func
          
    @intensity_func.setter 
    def intensity_func(self, intensity): 
        print('Setting intensity function to', intensity) 
        self._intensity_func = intensity

    @intensity_func.deleter 
    def intensity_func(self): 
        print('Deleting intensity function') 
        del self._intensity_func

    @staticmethod
    def powerlaw(coord, I0=30.0, R0=100*u.au, p=-0.4, z0=100*u.au, q=0.3):
        if 'R' not in coord.keys(): R = np.hypot(coord['x'], coord['y'])
        else: R = coord['R'] 
        z = coord['z']        
        A = I0*R0**-p*z0**-q
        return A*R**p*abs(z)**q
        
    @staticmethod
    def nuker(coord, I0=30.0, Rt=100*u.au, alpha=-0.5, gamma=0.1, beta=0.2):
        if 'R' not in coord.keys(): R = np.hypot(coord['x'], coord['y'])
        else: R = coord['R'] 
        A = I0*Rt**gamma
        return A*(R**-gamma) * (1+(R/Rt)**alpha)**((gamma-beta)/alpha)

   
class General2d(Velocity, Intensity, Tools):
    def __init__(self, grid):
        self.flags = {'disc': True, 'env': False}
        self.grid = grid
        self._velocity_func = General2d.keplerian
        self._intensity_func = General2d.powerlaw
        
    def make_model(self, incl, z_func, PA=0.0, get_2d=True, z_far=None, int_kwargs={}, vel_kwargs={}):

        #*************************************
        #MAKE TRUE GRID FOR NEAR AND FAR SIDES
        cos_incl, sin_incl = np.cos(incl), np.sin(incl)

        x_init, y_init = self.grid.XYZ[:2]
        x_true, y_true = x_init, y_init
        phi_true = np.arctan2(y_true, x_true)        
        R_true = np.hypot(x_init, y_init)
        z_true = z_func(R_true)

        if z_far is not None: z_true_far = z_far(R_true) 
        else: z_true_far = -z_true
            
        grid_true = {'near': [x_true, y_true, z_true, R_true, phi_true], 
                     'far': [x_true, y_true, z_true_far, R_true, phi_true]}

        #*******************************
        #COMPUTE PROPERTIES ON TRUE GRID
        avai_kwargs = [vel_kwargs, int_kwargs]
        avai_funcs = [self.velocity_func, self.intensity_func]
        true_kwargs = [isinstance(kwarg, dict) for kwarg in avai_kwargs]
        prop_kwargs = [kwarg for i, kwarg in enumerate(avai_kwargs) if true_kwargs[i]]
        prop_funcs = [func for i, func in enumerate(avai_funcs) if true_kwargs[i]]
        props = self._compute_prop(grid_true, prop_funcs, prop_kwargs)
        #Positive vel is positive along z, i.e. pointing to the observer, for that reason imposed a (-) factor to convert to the standard convention: (+) receding  
        if true_kwargs[0]:
            ang_fac = -sin_incl * np.cos(phi_true)
            props[0]['near'] *= ang_fac 
            props[0]['far'] *= ang_fac

        #***********************************
        #PROJECT PROPERTIES ON THE SKY PLANE
        #grid_axes = np.meshgrid(*self.grid.XYZgrid[:2]) #Avoiding this to keep backward compatibility with python2
        grid_axes = np.meshgrid(self.grid.XYZgrid[0], self.grid.XYZgrid[1]) 
        
        for side in ['near', 'far']:
            xt, yt, zt = grid_true[side][:3]
            x_pro, y_pro, z_pro = self._project_on_skyplane(xt, yt, zt, cos_incl, sin_incl)
            if PA: x_pro, y_pro = self._rotate_sky_plane(x_pro, y_pro, PA)             
            for prop in props: 
                prop[side] = griddata((x_pro, y_pro), prop[side], (grid_axes[0], grid_axes[1]))
        #*************************************

        return props
    
class Rosenfeld2d(Velocity, Intensity, Tools):
    """
    Host class for Rosenfeld+2013 model to describe the velocity field of a flared disc in 2D. 
    This model assumes a (Keplerian) double cone to account for the near and far sides of the disc 
    and solves analytical equations to find the line-of-sight velocity v_obs projected on the sky-plane. 
    
    Parameters
    ----------
    grid : array_like, shape (nrows, ncols)
       (x', y') map of the sky-plane onto which the disc velocity field will be projected.

    Attributes
    ----------
    velocity_func : function(coord, **kwargs) 
       Velocity function describing the kinematics of the disc. The argument coord is a dictionary
       of coordinates (e.g. 'x', 'y', 'z', 'r', 'R', 'theta', 'phi') at which the function will be evaluated. 
       Additional arguments are optional and depend upon the function definition, e.g. Mstar=1.0*u.Msun
    """

    def __init__(self, grid):
        self.flags = {'disc': True, 'env': False}
        self.grid = grid
        self._velocity_func = Rosenfeld2d.keplerian
        self._intensity_func = Rosenfeld2d.powerlaw

    def _get_t(self, A, B, C):
        t = []
        for i in range(self.grid.NPoints):
            p = [A, B[i], C[i]]
            t.append(np.sort(np.roots(p)))
        return np.array(t)

    def make_model(self, incl, psi, PA=0.0, get_2d=True, int_kwargs={}, vel_kwargs={}):
        """
        Executes the Rosenfeld+2013 model.

        Parameters
        ----------
        incl : scalar
           Inclination of the disc midplane with respect to the x'y' plane; pi/2 radians is edge-on.
    
        psi : scalar
           Opening angle of the cone describing the velocity field of the gas emitting layer in the disc; 
           0 radians returns the projected velocity field of the disc midplane (i.e no conic emission). 

        PA : scalar, optional
           Position angle in radians. Measured from North (+y) to East (-x).

        get_2d : bool, optional
           If True regrids the resulting velocity field into a 2D map and stores it in the attribute 'velocity2d'. 

        Attributes
        ----------
        velocity : array_like, size (n,)
           Velocity field computed using the Rosenfeld+2013 model.

        velocity2d : array_like, size (nx, ny)
           If set get_2d=True: Velocity field computed using the Rosenfeld+2013 model, reshaped to 2D to facilitate plotting.
        """
        if PA: x_plane, y_plane = Rosenfeld2d._rotate_sky_plane(self.grid.XYZ[0], self.grid.XYZ[1], -PA)
        else: x_plane, y_plane = self.grid.XYZ[:2]

        cos_incl = np.cos(incl)
        sin_incl = np.sin(incl)
        y_plane_cos_incl = y_plane/cos_incl

        #**********************
        #ROSENFELD COEFFICIENTS
        fac = -2*np.sin(psi)**2
        A = np.cos(2*incl) + np.cos(2*psi)
        B = fac * 2*(sin_incl/cos_incl) * y_plane
        C = fac * (x_plane**2 + (y_plane_cos_incl)**2)
        t = self._get_t(A,B,C).T

        #****************************
        #ROSENFELD CONVERSION X<-->X'
        x_true_near = x_plane
        y_true_near = y_plane_cos_incl + t[1]*sin_incl
            
        x_true_far = x_plane
        y_true_far = y_plane_cos_incl + t[0]*sin_incl
        
        #np.hypot 2x faster than np.linalg.norm([x,y], axis=0)
        R_true_near = np.hypot(x_true_near, y_true_near) 
        R_true_far = np.hypot(x_true_far, y_true_far)

        z_true_near = t[1] * cos_incl
        z_true_far = t[0] * cos_incl 

        phi_true_near = np.arctan2(y_true_near, x_true_near)        
        phi_true_far = np.arctan2(y_true_far, x_true_far)        

        #****************************
            
        grid_true =  {'near': [x_true_near, y_true_near, z_true_near, R_true_near, phi_true_near], 
                      'far': [x_true_far, y_true_far, z_true_far, R_true_far, phi_true_far]}

        #*******************************
        #COMPUTE PROPERTIES ON TRUE GRID
        avai_kwargs = [vel_kwargs, int_kwargs]
        avai_funcs = [self.velocity_func, self.intensity_func]
        true_kwargs = [isinstance(kwarg, dict) for kwarg in avai_kwargs]
        prop_kwargs = [kwarg for i, kwarg in enumerate(avai_kwargs) if true_kwargs[i]]
        prop_funcs = [func for i, func in enumerate(avai_funcs) if true_kwargs[i]]
        props = self._compute_prop(grid_true, prop_funcs, prop_kwargs)
        #Positive vel is positive along z, i.e. pointing to the observer, for that reason imposed a (-) factor to convert to the standard convention: (+) receding  
        if true_kwargs[0]:
            ang_fac_near = -sin_incl * np.cos(phi_true_near)
            ang_fac_far = -sin_incl * np.cos(phi_true_far)
            props[0]['near'] *= ang_fac_near 
            props[0]['far'] *= ang_fac_far
                
        #*************************************

        return [{side: prop[side].reshape(self.grid.Nodes[:2]) for side in ['near', 'far']} for prop in props]

