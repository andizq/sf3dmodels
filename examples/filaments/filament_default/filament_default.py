import sf3dmodels.filament as sf
import sf3dmodels.Plot_model as pm
from sf3dmodels.utils.units import pc
import sf3dmodels.rt as rt

f1 = sf.FilamentModel([0,0,0], [0,0,1], -0.2*pc, 0.2*pc, 0.01*pc)
f1.cylinder(0.1*pc, 1e-3*pc, temp_pars = [500, 0.02*pc, -0.3], abund_pars = 1e-4)

lims=[-0.3*pc,0.3*pc]
pm.scatter3D(f1.GRID, f1.density, f1.density.min(), axisunit = pc,
             colordim = f1.temperature, 
             colorlabel = 'T [K]',
             xlim=lims, ylim=lims, zlim=lims,
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
lime.submodel(prop, output='datatab.dat', folder='./', lime_header=True, lime_npoints=True)
