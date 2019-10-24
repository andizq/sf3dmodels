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

    def gaussian_profile(self, z, stddev):
        return np.exp(-0.5*(z/stddev)**2)

    def _powerlaw_cavity1d(self, R_cav=30*au, h=5*au, n_cav=1e16, power=-1.0, dn_cav=1e-4, coord={'R': None, 'z': None}):
        R = coord['R']
        z = coord['z']
        a_cav = n_cav
        if R < R_cav: a_cav *= dn_cav 
        val = a_cav*(R/R_cav)**power * self.gaussian_profile(z, h)
        return val

    def powerlaw_cavity(self, grid=None, R_cav=30*au, h=5*au, n_cav=1e16, power=-1.0, dn_cav=1e-4, coord={'R': None, 'z': None}, func_1d=False):
        if func_1d: return self._powerlaw_cavity1d
        if coord['R'] is not None:
            R = np.asarray(coord['R'])
            z = np.asarray(coord['z'])
            val = np.zeros(R.shape)
            cav = R < R_cav
            val = n_cav*(R/R_cav)**power * self.gaussian_profile(z, h)
            #val[cav] *= dn_cav
            val = np.where(cav,dn_cav*val,val)
        if grid is not None:
            #print (coord)
            profile = (grid.rRTP[1]/R_cav)**power
            val = np.where(grid.rRTP[1] > R_cav, n_cav*profile, dn_cav*n_cav*profile) * self.gaussian_profile(grid.XYZ[2], h)
        return val

    
