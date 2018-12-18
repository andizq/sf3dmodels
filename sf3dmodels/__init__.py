# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
SF3dmodels : Star-Forming regions 3d modelling package
"""
# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
import version_error
# ----------------------------------------------------------------------------

if not _ASTROPY_SETUP_:
    # For egg_info test builds to pass, put package imports here.
    #from .example_mod import *
    
    from . import Model
    from . import Resolution
    from . import Plot_model
#    from . import Tools_v0
#    from . import grid
    
    from . import Utils
    from . import create_cylinder
    from . import create_parabola
    Res = Resolution
    U = Utils
    pass

#from .Model import grid, sphe_cart, streamline, density_Env_Disc

