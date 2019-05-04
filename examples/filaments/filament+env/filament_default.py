import numpy as np
import sf3dmodels.filament as sf
import sf3dmodels.Plot_model as pm
from sf3dmodels.utils.units import pc
import sf3dmodels.rt as rt

f1 = sf.FilamentModel([0,0,0], [1,0,0], -0.2*pc, 0.2*pc, 4e-3*pc)

def new_temp(R,theta,z, *temp_pars):
    TR, R0, pR, zh = temp_pars
    cR = TR*R0**-pR
    return cR*R**pR * np.exp(np.abs(z)/zh)

f1.func_temp = new_temp

f1.cylinder(0.1*pc, 1e-4*pc,
            dens_pars = [7e10, 0.03*pc, 2.0],
            temp_pars = [200, 0.02*pc, -0.15, -0.17*pc],
            abund_pars = 1e-4)

pm.scatter3D(f1.GRID, f1.density, np.min(f1.density), axisunit = pc,
             colordim = f1.density,
             colorlabel = 'T [K]',
             NRand = 10000, show=True)

prop = {'dens_H': f1.density,
        'temp_gas': f1.temperature,
        'abundance': f1.abundance,
        'gtdratio': f1.gtdratio,
        'vel_x': f1.vel.x,
        'vel_y': f1.vel.y,
        'vel_z': f1.vel.z,
        }

lime = rt.Lime(f1.GRID)
lime.submodel(prop, output='filament.dat')
