"""
Disc models collection
======================
Classes: Rosenfeld2d
"""

from ..utils.constants import G
from ..utils import units as u
import numpy as np
import numbers
from scipy.interpolate import griddata


class Rosenfeld2d(object):
    """
    Host class for Rosenfeld+2013 toy model to describe the velocity field of a flared disc in 2D. 
    This model assumes a (Keplerian) double cone to account for the near and far sides of a flared disc 
    and projects their line-of-sight velocity v_obs on the sky-plane. 
    
    Parameters
    ----------
    grid : array_like, shape (nrows, ncols)
       (x', y') map of the sky-plane onto which the disc velocity field will be projected.
       
    Mstar : scalar
       Mass of the star to compute keplerian rotation.
    
    incl : scalar
       Inclination of the disc midplane with respect to the x'y' plane; pi/2 radians is edge-on.
    
    psi : scalar
       Opening angle of the cone describing the velocity field of the gas emitting layer in the disc; 
       0 radians returns the projected velocity field of the disc midplane (i.e no conic emission). 

    get2d : bool
       If True returns the velocity field regrided in a 2D map.
       
    velocity : str from ['keplerian', 'keplerian_vertical']
       Orbital velocity function to compute the velocity field. 

    """

    def __init__(self, grid, Mstar, incl, psi, get2d = True, velocity='keplerian', z_func=False, PA=0.0):
        self.flags = {'disc': True, 'env': False}
        self.grid = grid
        self.n = self.grid.NPoints
        self.Mstar = Mstar
        self.incl = incl
        if not z_func: self.psi = psi
        else:            
            x_plane, y_plane = self.grid.XYZ[:2]
            x_mid = x_plane
            self.y_mid = y_plane * np.cos(incl)**-1 # doing y_mid*=np.cos(ang)**-1 changes also the grid.XYZ object
            R_mid = np.linalg.norm([x_mid, self.y_mid], axis=0)
            self.R_mid = np.linalg.norm([x_mid, self.y_mid], axis=0)
            self.z_val = z_func(R_mid)
            #self.psi = np.arctan(z_func(self.grid.rRTP[1]) / self.grid.rRTP[1])
            self.psi = np.arctan(z_func(R_mid) / R_mid)
            #self.psi_far = np.arctan(-z_func(R_mid) / R_mid)

        self.z_func = z_func
        self.PA = PA
        self.grid_true = self.cone_to_2d() #(x,y,z) grid as a function of the x'y' plane coordinates.
        velocity_func = {'keplerian': self.velocity_keplerian, 'keplerian_vertical': self.velocity_keplerian_vertical}
        self.velocity = velocity_func[velocity]()
        if get2d: 
            self.velocity2d = {}
            for side in ['near', 'far']:
                self.velocity2d[side] = self.convert_array_to_matrix(self.velocity[side])

    def make_cone(self, Mstar, psi):
        pass
    def get_channel(self, v_chan, v_turb=0): pass
    def convert_array_to_matrix(self, vec):
        matrix = np.zeros(self.grid.Nodes[self.grid.Nodes>1])
        k = 0
        for j in range(self.grid.Nodes[1]):
            for i in range(self.grid.Nodes[0]):
                matrix[j,i] = vec[k]
                k+=1
        return matrix
        
    def _get_t(self, A, B, C):
        t = []
        if isinstance(self.psi, numbers.Number): 
            for i in range(self.n):
                p = [A, B[i], C[i]]
                t.append(np.sort(np.roots(p)))
        else: 
            for i in range(self.n):
                p = [A[i], B[i], C[i]]
                t.append(np.sort(np.roots(p)))
        return np.array(t)
            
    def solve_quadratic(self, x, y):
        fac = -2*np.sin(self.psi)**2
        A = np.cos(2*self.incl) + np.cos(2*self.psi)
        B = fac * 2*np.tan(self.incl) * y
        C = fac * (x**2 + (y / np.cos(self.incl))**2)
        return self._get_t(A,B,C)

    def velocity_keplerian(self):
        vel = {}
        for side in ['near', 'far']:
            x, y, z, R = self.grid_true[side]
            phi = np.arctan2(y, x)
            ang_fac = np.sin(self.incl) * np.cos(phi)
            vel[side] = -np.sqrt(G * self.Mstar/R) * ang_fac #Positive vel means positive along z, which means approaching to the observer, for that reason imposed a (-) factor.
        return vel

    def velocity_keplerian_vertical(self):
        vel = {}
        for side in ['far', 'near']:
            x, y, z, R = self.grid_true[side]
            r = np.linalg.norm([R,z], axis=0)
            phi = np.arctan2(y, x) 
            ang_fac = np.sin(self.incl) * np.cos(phi)
            vel[side] = -np.sqrt(G * self.Mstar/r**3) * R * ang_fac #Positive vel means positive along z, which means approaching to the observer, for that reason imposed a (-) factor.

        x_near, y_near, z_near = self.grid_true['near'][:3]
        x_far, y_far, z_far = self.grid_true['far'][:3]
        
        """
        vel['near'] = np.ma.masked_where(z_near>80*u.au, vel['near'])
        vel['far'] = np.ma.masked_where(np.logical_or(abs(z_far)>80*u.au, ~vel['near'].mask), vel['far'])
        """
        x_dep_near = x_near
        y_dep_near = y_near * np.cos(self.incl) - z_near * np.sin(self.incl)
        z_dep_near = y_near * np.sin(self.incl) + z_near * np.cos(self.incl)
        x_dep_near, y_dep_near = self.rotate_sky_plane(x_dep_near, y_dep_near, self.PA)

        x_dep_far = x_far
        y_dep_far = y_far * np.cos(self.incl) - z_far * np.sin(self.incl)
        z_dep_far = y_far * np.sin(self.incl) + z_far * np.cos(self.incl)
        x_dep_far, y_dep_far = self.rotate_sky_plane(x_dep_far, y_dep_far, self.PA)

        grid_axes = np.meshgrid(*np.array(self.grid.XYZgrid[:2])) 
        
        z_near = griddata((x_dep_near, y_dep_near), z_near, (grid_axes[0], grid_axes[1])).flatten()
        z_far = griddata((x_dep_far, y_dep_far), z_far, (grid_axes[0], grid_axes[1])).flatten()

        vel['near'] = np.ma.masked_where(np.logical_or(z_near>150*u.au, np.isnan(z_near)), vel['near'])
        vel['far'] = np.ma.masked_where(np.logical_or(-z_far>150*u.au, ~vel['near'].mask), vel['far'])
        
        vel['near'] = np.where(~vel['near'].mask, griddata((x_dep_near, y_dep_near), vel['near'], (grid_axes[0], grid_axes[1])).flatten(), np.nan)
        vel['far'] = np.where(~vel['far'].mask, griddata((x_dep_far, y_dep_far), vel['far'], (grid_axes[0], grid_axes[1])).flatten(), np.nan)
        #vel['near'] = griddata((x_dep, y_dep), vel['near'], (grid_axes[0], grid_axes[1])).flatten()
        
        return vel
            
    def rotate_sky_plane(self, x, y, ang):
        xy = np.array([x,y])
        cos_ang = np.cos(ang)
        sin_ang = np.sin(ang)
        rot = np.array([[cos_ang, -sin_ang],
                        [sin_ang, cos_ang]])
        return np.dot(rot, xy)

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

    def cone_to_2d(self):
        if self.PA != 0.0: 
            x_plane, y_plane = self.rotate_sky_plane(self.grid.XYZ[0], self.grid.XYZ[1], -self.PA)
            if self.z_func:
                x_init, y_init = self.grid.XYZ[:2]
                R_init = np.hypot(x_init, y_init)
                z_init = self.z_func(R_init)
                x_incl, y_incl, z_incl = self.rotate_sky_plane3d(x_init, y_init, z_init, self.incl, axis='x')
                x_inclf, y_inclf, z_inclf = self.rotate_sky_plane3d(x_init, y_init, -z_init, self.incl, axis='x')
                
                """
                x_pa, y_pa, z_pa = self.rotate_sky_plane3d(x_incl, y_incl, z_incl, -self.PA)
                x_paf, y_paf, z_paf = self.rotate_sky_plane3d(x_inclf, y_inclf, z_inclf, -self.PA)
                """
                x_pa, y_pa, z_pa = self.rotate_sky_plane3d(x_init, y_init, z_init, -self.PA)
                x_paf, y_paf, z_paf = self.rotate_sky_plane3d(x_init, y_init, -z_init, -self.PA)

                #xx, self.y_mid, self.z_val = self.rotate_sky_plane3d(self.grid.XYZ[0], self.y_mid, self.z_val, -self.PA)
                #self.R_mid = np.hypot(xx, self.y_mid)
                
        else:
            x_plane, y_plane = self.grid.XYZ[:2]

        t = self.solve_quadratic(x_plane,y_plane).T        

        x_true_near = x_plane
        y_true_near = y_plane / np.cos(self.incl) + t[1]*np.sin(self.incl)
            
        x_true_far = x_plane
        y_true_far = y_plane / np.cos(self.incl) + t[0]*np.sin(self.incl)
        
        R_true_near = np.linalg.norm([x_true_near, y_true_near], axis=0)
        R_true_far = np.linalg.norm([x_true_far, y_true_far], axis=0)

        if not self.z_func: 
            z_true_near = t[1] * np.cos(self.incl) 
            z_true_far = t[0] * np.cos(self.incl) 
        else: 
            x_pa, y_pa = x_init, y_init
            z_true_far = self.z_func(R_true_far) 
            self.R_mid = np.hypot(x_pa, y_pa)
            R_true_near = self.R_mid
            z_true_near = self.z_func(self.R_mid)
            x_true_near =  x_pa
            y_true_near =  y_pa

            self.R_midf = np.hypot(x_paf, y_paf)
            R_true_far = self.R_mid
            z_true_far = -z_true_near #-self.z_func(self.R_midf)
            x_true_far =  x_pa #x_paf
            y_true_far =  y_pa #y_paf

            #z_true_far = -z_true_near
            #y_true_near = self.y_mid 

            #z_true_near = z_pa
            
            #z_true_far = -self.z_func(self.R_mid)
            #z_true_near = self.z_func(self.grid.rRTP[1])
            #z_true_far = -z_true_near
            
        return {'near': [x_true_near, y_true_near, z_true_near, R_true_near], 
                'far': [x_true_far, y_true_far, z_true_far, R_true_far]}
        
        
