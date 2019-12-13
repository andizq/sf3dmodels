"""
Envelope models collection
==========================
Classes: Powerlaw, (missing: Ulrich, Uniform, Shells)
"""

from ..utils.units import au, pc
import numpy as np

#****************
#TRANSITION DISCS
#****************
class Powerlaw(object):
    """
    Host class for powerlaw models in envelopes.
    """
    def __init__(self):
        self.flags = {'disc': False, 'env': True}
        self.func_scalar = {'monotonic': self._monotonic_scalar}
        self.background = 1.0
        
    def _shells_scalar(self): pass
    def shells(self): pass

    def _monotonic_scalar(self, 
                          val_min = None, r_min = None, 
                          r_max = None, q = None, 
                          grid=None, coord=None):
        r = coord['r']
        if r >= r_min and r <= r_max: 
            rq = r**q
            val = val_min/r_min**q * rq
        else: val = self.background
        return val

    def monotonic(self,
                  val_min = 1e15, r_min = au, r_max = 100*au, q = -1.0, 
                  grid=None, coord=None, func_scalar=False):

        if func_scalar: return self.func_scalar['monotonic'] #If the coord input is scalar
        if coord is not None:
            print ('Computing Envelope property using 1D power-law...')
            r = np.asarray(coord['r'])
        if grid is not None:
            r = grid.rRTP[0]
            print ('Computing Envelope property using 3D power-law...')
        rq = np.where((r >= r_min) & (r <= r_max), r**q, 0.)
        val = val_min/r_min**q * rq 
        val[val==0.] = self.background
        return val

    
