import numpy as np
import sf3dmodels.filament as sf
import sf3dmodels.Plot_model as pm
from sf3dmodels.utils.units import pc
import sf3dmodels.rt as rt

f1 = sf.FilamentModel([0,0,0], [0,0,1], -0.2*pc, 0.2*pc, 8e-3*pc)

def new_width(z, *width_pars):
    w0, period = width_pars
    return w0*((0.5*np.sin(z*2*np.pi/period)**2) + 0.5)
    
def new_abund(R,theta,z, *abund_pars):
    a0, R0, p = abund_pars
    ca = a0*R0**-p
    return ca*R**p

def new_temp(R,theta,z, *temp_pars):
    TR, R0, pR, zh = temp_pars
    cR = TR*R0**-pR
    return cR*R**pR * np.exp(np.abs(z)/zh)
    
f1.func_width = new_width
f1.func_abund = new_abund
f1.func_temp = new_temp

f1.cylinder([0.1*pc, 0.3*pc], 1e-4*pc, 
            abund_pars = [1e-6, 0.05*pc, -0.15],
            temp_pars = [200, 0.02*pc, -0.15, -0.17*pc],
            dummy_frac = 0.5)

lims=np.array([-0.3,0.3])*pc
pm.scatter3D(f1.GRID, f1.density, np.mean(f1.density), axisunit = pc,
             colordim = f1.temperature,
             colorlabel = 'T [K]',
             NRand = 10000, 
             cmap = 'nipy_spectral_r',
             xlim=lims, ylim=lims, zlim=lims,
             azim=45, elev=15, output='fig_filament_temp.png', show=True)

pm.scatter3D(f1.GRID, f1.density, np.min(f1.density), axisunit = pc,
             colordim = f1.abundance,
             colorlabel = 'Molec. abund.',
             NRand = 10000, 
             cmap = 'nipy_spectral_r',
             xlim=lims, ylim=lims, zlim=lims,
             azim=45, elev=15, output='fig_filament_abund.png', show=True)

prop = {'dens_H': f1.density,
        'temp_gas': f1.temperature,
        'abundance': f1.abundance,
        'gtdratio': f1.gtdratio,
        'vel_x': f1.vel.x,
        'vel_y': f1.vel.y,
        'vel_z': f1.vel.z,
        }

lime = rt.Lime(f1.GRID)
lime.submodel(prop, output='datatab.dat', folder='./', lime_header=True, lime_npoints=True)
