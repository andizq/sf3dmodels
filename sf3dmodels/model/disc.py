"""
Disc models collection
======================
Classes: Transition, (missing: Keplerian --> Burgers, Pringle+1981, Keto+2010)
"""

from ..utils.units import au
import numpy as np

#****************
#TRANSITION DISCS
#****************
class Transition(object):
    """
    Host class for transition disc models.
    Included: van Dishoeck+2015
    """
    #Missing: include default models, e.g vandishoeck2015 = {'dens': Transition.func1, 'temp': Transition.func2} 
    
    def __init__(self):
        self.flags = {'disc': True, 'env': False}
        func_1d = {'powerlaw_cavity': self._powerlaw_cavity1d}

    def constant_profile(self, x): return x
    def gaussian_profile(self, x, x_mean, stddev):
        return np.exp(-0.5*((x-x_mean)/stddev)**2)

    def _powerlaw_cavity1d(self, n_cav=1e16, power=-1.0, dn_cav=1e-4,
                           R_cav=30*au, #Radial cavity, must be > 0 
                           z_mean=0, z_stddev=5*au, #Gaussian mean and standard deviation for the disc scale-height
                           phi_mean=0, phi_stddev=None,
                           grid=None, coord=None, func_1d=False):
        R = coord['R']
        z = coord['z']
        a_cav = n_cav
        if R < R_cav: a_cav *= dn_cav 
        if phi_stddev is not None: 
            phi = coord['phi'] - phi_mean
            phi = (phi>np.pi)*(phi-2*np.pi) + (phi<=np.pi)*phi #Making the grid symmetric with respect to the gaussian val phi_mean
            phi_val = self.gaussian_profile(phi, 0, phi_stddev)
        else: phi_val = 1.0
        val = a_cav*(R/R_cav)**power * self.gaussian_profile(z, z_mean, z_stddev) * phi_val
        return val

    def powerlaw_cavity(self, n_cav=1e16, power=-1.0, dn_cav=1e-4,
                        R_cav=30*au, #Radial cavity, must be > 0 
                        z_mean=0, z_stddev=5*au, #Gaussian mean and standard deviation for the disc scale-height
                        phi_mean=0, phi_stddev=None,
                        grid=None, coord=None, func_1d=False):
        if func_1d: return self._powerlaw_cavity1d #If the coord input is scalar
        if coord is not None:
            R = np.asarray(coord['R'])
            z = np.asarray(coord['z'])
            if phi_stddev is not None: 
                phi = np.asarray(coord['phi']) - phi_mean
                phi = np.where(phi > np.pi, phi-2*np.pi, phi)
                phi_val = self.gaussian_profile(phi, 0, phi_stddev)
            else: phi_val = 1.0
            val = n_cav*(R/R_cav)**power * self.gaussian_profile(z, z_mean, z_stddev) * phi_val
            #val = np.zeros(R.shape)
            #cav = R < R_cav
            #val[cav] *= dn_cav
            val = np.where(R < R_cav, dn_cav*val, val)
        if grid is not None:
            profile = (grid.rRTP[1]/R_cav)**power
            if phi_stddev is not None: 
                phi = grid.rRTP[3] - phi_mean
                phi = np.where(phi > np.pi, phi-2*np.pi, phi)
                print (phi.max()*180/np.pi, phi.min()*180/np.pi, grid.rRTP[3].max(), grid.rRTP[3].min())
                phi_val = self.gaussian_profile(phi, 0, phi_stddev)
            else: phi_val = 1.0
            val = np.where(grid.rRTP[1] > R_cav, n_cav*profile, dn_cav*n_cav*profile) * self.gaussian_profile(grid.XYZ[2], z_mean, z_stddev) * phi_val
        return val

    
