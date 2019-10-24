"""
Transition Disc Example
=======================

This example uses the gas and dust density profiles inferred for the transitional disc HD135344B.
"""
#******************
#sf3dmodels modules
#******************
from sf3dmodels import Model, Plot_model
import sf3dmodels.utils.units as u
import sf3dmodels.rt as rt
from sf3dmodels.model import disc
from sf3dmodels.grid import Grid
#******************
#External libraries
#******************
import numpy as np
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import time

t0 = time.time()
#******************
#Model parameters
#******************
kwargs_dens={'dn_cav': 1e-5, 'n_cav': 1.5e17, 'R_cav': 25*u.au, 'power':-0.9}
kwargs_dtg={'dn_cav': 1e-2, 'n_cav': 0.01, 'R_cav': 35*u.au, 'power': 0}

#******************
#Grid creation
#******************
transition = disc.Transition()
init_grid = Grid()
grid = init_grid.random(function=transition.powerlaw_cavity, power=0.4, r_size=70*u.au, kwargs_func=kwargs_dens, normalization=1e16, npoints=100000)

#*****************************************
#Computing relevant physical distributions
#*****************************************
density = transition.powerlaw_cavity(grid=grid, **kwargs_dens)
gtdratio = 1./transition.powerlaw_cavity(grid=grid, **kwargs_dtg) #Converting to gas-to-dust ratio
temperature = np.zeros(grid.NPoints) + 100. #Constant temperature

prop = {'dens_H2': density,
        'gtdratio': gtdratio,
        'temp_gas': temperature}

#**********************
#Writing files for LIME 
#**********************
lime = rt.Lime(grid)
lime.submodel(prop, output='datatab.dat', lime_npoints = True, lime_header = True)
print('Output columns', lime.columns)

print ('Ellapsed time:', time.time()-t0)
#******************
#1D plotting
#****************
#fig, ax = plt.subplots(ncols=2, figsize=(6.4,4.8))
#r_plot = np.arange(1,70)*u.au
#transition.powerlaw_cavity(**kwargs_dens)


#******************
#3D plotting
#****************
lims = np.array([-70,70])
canvas3d = Plot_model.Canvas3d(ax_kw={'xlim': lims, 'ylim': lims, 'zlim': lims, 'azim': -50, 'elev': 30})
ax = canvas3d.ax 
sp = canvas3d.scatter_random(grid, density/1e6, 1e10, GRID_unit=u.au, power=0, NRand=10000, prop_min=1.0, #function arguments
                             marker = 'o', cmap = 'jet', s = 3, edgecolors = 'none', vmin=1e6, norm = colors.LogNorm()) #Scatter kwargs
cbar = plt.colorbar(sp)
cbar.ax.set_ylabel(r'$n_{gas}$ cm$^{-3}$')
ax.set_xlabel('au')

plt.savefig('transition_disc_3d.png')
plt.show()
