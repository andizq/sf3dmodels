"""
2D Disc models
==============
Classes: Rosenfeld2d, General2d, Velocity, Intensity
"""

from ..utils.constants import G
from ..utils import units as u
import numpy as np
import numbers
from scipy.interpolate import griddata

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

    def keplerian(coord, Mstar=u.Msun):
        if 'R' not in coord.keys(): R = np.hypot(coord['x'], coord['y'])
        else: R = coord['R'] 
        return np.sqrt(G*Mstar/R) 
    
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
    def intensity_func(self, vel): 
        print('Setting intensity function to', vel) 
        self._intensity_func = vel

    @intensity_func.deleter 
    def intensity_func(self): 
        print('Deleting intensity function') 
        del self._intensity_func

    def powerlaw(coord, I0=20.0, R0=10*u.au, p=-1.0):
        if 'R' not in coord.keys(): R = np.hypot(coord['x'], coord['y'])
        else: R = coord['R'] 
        A = I0*R0**-p
        return A*R**p
    
    def nuker(coord, I0=20.0, Rt=10*u.au, alpha=1.0, gamma=1.0, beta=2.0):
        if 'R' not in coord.keys(): R = np.hypot(coord['x'], coord['y'])
        else: R = coord['R'] 
        A = I0*Rt**gamma
        return A*(R**-gamma) * (1+(R/Rt)**alpha)**((gamma-beta)/alpha)
   
class General2d(Velocity):
    def __init__(self, grid):
        self.flags = {'disc': True, 'env': False}
        self.grid = grid
        self._velocity_func = General2d.keplerian

    @staticmethod
    def _rotate_sky_plane(x, y, ang):
        xy = np.array([x,y])
        cos_ang = np.cos(ang)
        sin_ang = np.sin(ang)
        rot = np.array([[cos_ang, -sin_ang],
                        [sin_ang, cos_ang]])
        return np.dot(rot, xy)
    
    @staticmethod
    def _project_on_skyplane(x, y, z, cos_incl, sin_incl):
        x_pro = x
        y_pro = y * cos_incl - z * sin_incl
        z_pro = y * sin_incl + z * cos_incl
        return x_pro, y_pro, z_pro

    def _compute_velocity(self, incl, grid, get_vel2d, PA, **vel_kwargs):
        vel = {}
        cos_incl, sin_incl = np.cos(incl), np.sin(incl)
        grid_axes = np.meshgrid(*np.array(self.grid.XYZgrid[:2])) 
        for side in ['near', 'far']:
            x, y, z, R, phi = grid[side]
            ang_fac = sin_incl * np.cos(phi)
            coord = {'x': x, 'y': y, 'z': z, 'phi': phi, 'R': R}
            #Positive vel means positive along z, which means approaching to the observer, for that reason imposed a (-) factor
            vel[side] = -self.velocity_func(coord, **vel_kwargs)*ang_fac 

            #The following lines should be outside to allow using the same for the intensity. 
            x_pro, y_pro, z_pro = self._project_on_skyplane(x, y, z, cos_incl, sin_incl)
            if PA: x_pro, y_pro = self._rotate_sky_plane(x_pro, y_pro, PA) 

            vel[side] = griddata((x_pro, y_pro), vel[side], (grid_axes[0], grid_axes[1]))

        if get_vel2d: self.velocity2d = vel
        self.velocity = {side: vel[side].flatten() for side in ['near', 'far']}


    def make_model(self, incl, z_func, PA=0.0, get_vel2d=True, z_far=None, int_kwargs={}, **vel_kwargs):
        x_init, y_init = self.grid.XYZ[:2]
        x_true, y_true = x_init, y_init
        phi_true = np.arctan2(y_true, x_true)        

        R_true = np.hypot(x_init, y_init)
        z_true = z_func(R_true)
        if z_far is not None: z_true_far = z_far(R_true) 
        else: z_true_far = -z_true
            
        grid_true = {'near': [x_true, y_true, z_true, R_true, phi_true], 
                     'far': [x_true, y_true, z_true_far, R_true, phi_true]}

        self._compute_velocity(incl, grid_true, get_vel2d, PA, **vel_kwargs)
        #self.get_channel = self._get_channel  
        #if int_kwargs: self._compute_intensity(incl, grid_true, PA, **int_kwargs)

    def get_channel(self, v_chan, T, v_turb=0.0, mmol=2*u.amu):
        #Use here the velocity and intensity fields. Double check first whether they exist and if not raise an error message
        pass
    
class Rosenfeld2d(Velocity):
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

    @staticmethod
    def _rotate_sky_plane(x, y, ang):
        xy = np.array([x,y])
        cos_ang = np.cos(ang)
        sin_ang = np.sin(ang)
        rot = np.array([[cos_ang, -sin_ang],
                        [sin_ang, cos_ang]])
        return np.dot(rot, xy)

    def _get_t(self, A, B, C):
        t = []
        for i in range(self.grid.NPoints):
            p = [A, B[i], C[i]]
            t.append(np.sort(np.roots(p)))
        return np.array(t)

    def _compute_velocity(self, sin_incl, grid, get_vel2d, **vel_kwargs):
        vel = {}
        for side in ['near', 'far']:
            x, y, z, R = grid[side]
            phi = np.arctan2(y, x)
            ang_fac = sin_incl * np.cos(phi)
            coord = {'x': x, 'y': y, 'z': z, 'phi': phi, 'R': R}
            vel[side] = -self.velocity_func(coord, **vel_kwargs)*ang_fac #Positive vel means positive along z, which means approaching to the observer, for that reason imposed a (-) factor.

        if get_vel2d: self.velocity2d = {side: vel[side].reshape((*self.grid.Nodes[:2])) for side in ['near', 'far']}
        self.velocity = vel

    def make_model(self, incl, psi, PA=0.0, get_vel2d=True, **vel_kwargs):
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

        get_vel2d : bool, optional
           If True regrids the resulting velocity field into a 2D map and stores it in the attribute 'velocity2d'. 

        Attributes
        ----------
        velocity : array_like, size (n,)
           Velocity field computed using the Rosenfeld+2013 model.

        velocity2d : array_like, size (nx, ny)
           If set get_vel2d=True: Velocity field computed using the Rosenfeld+2013 model, reshaped to 2D to facilitate plotting.
        """
        if PA: 
            x_plane, y_plane = Rosenfeld2d._rotate_sky_plane(self.grid.XYZ[0], self.grid.XYZ[1], -PA)
        else:
            x_plane, y_plane = self.grid.XYZ[:2]

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
        #**********************
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
        #****************************
            
        grid_true =  {'near': [x_true_near, y_true_near, z_true_near, R_true_near], 
                      'far': [x_true_far, y_true_far, z_true_far, R_true_far]}
        
        self._compute_velocity(sin_incl, grid_true, get_vel2d, **vel_kwargs)

    #*********************************

    def get_channel(self, v_chan, v_turb=0): pass

            
    def rotate_sky_plane3d(self, x, y, z, ang, axis='z'):
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
