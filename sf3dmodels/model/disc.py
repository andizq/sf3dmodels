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
    
    def gaussian_profile(self, z, stddev):
        return np.exp(-0.5*(z/stddev)**2)

    def powerlaw_cavity(self, grid=None, R_cav=30*au, h=5*au, n_cav=1e16, power=-1.0, dn_cav=1e-4, loc={'R': None, 'z': None}):
        if loc['R'] is not None:
            R = np.asarray(loc['R'])
            z = np.asarray(loc['z'])
            val = np.zeros(R.shape)
            cav = R < R_cav
            val = n_cav*(R/R_cav)**power * self.gaussian_profile(z, h)
            val[cav] *= dn_cav
        if grid is not None:
            print (loc)
            profile = (grid.rRTP[1]/R_cav)**power
            val = np.where(grid.rRTP[1] > R_cav, n_cav*profile, dn_cav*n_cav*profile) * self.gaussian_profile(grid.XYZ[2], h)
        return val

    
