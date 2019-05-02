import sf3dmodels.Plot_model as pm
import sf3dmodels.filament as sf
import numpy as np
from sf3dmodels.utils.units import pc

f1 = sf.FilamentModel([0,0,0],
                      [0,0,1],
                      -5e15,
                      5e15,
                      2e14, mirror=True*0)

f1.func_width = lambda z, *width_pars: width_pars[0]*((0.5*np.sin(z*2*np.pi/1e16)**2) + 0.5)
f1.cylinder(3e15, 1e-3*pc)

pm.scatter3D(f1.GRID, f1.density, f1.density.min(),
             axisunit = pc,
             colordim = f1.theta, 
             NRand = 10000, show=True*1)

