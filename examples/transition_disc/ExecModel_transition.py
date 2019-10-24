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
import matplotlib
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
grid = init_grid.random(func=transition.powerlaw_cavity, power=0.4, r_size=70*u.au, kwargs_func=kwargs_dens, normalization=1e16, npoints=100000)

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
#******************
matplotlib.rc('font', size=15)
fig = plt.figure(figsize=(12,8.8))
ax0 = fig.add_subplot(222)
r_plot = np.arange(1,60)*u.au
dens1d = transition.powerlaw_cavity(**kwargs_dens, coord={'R': r_plot, 'z': 0})
ax0.plot(r_plot/u.au, dens1d/1e6, lw=3, color='green')
#ax0.set_xscale('log')
ax0.set_yscale('log')
ax0.set_xlabel('au', labelpad=-40)
ax0.set_ylabel('H$_2$ density [cm$^{-3}$]')
ax0.tick_params(axis='y', labelcolor='green')

ax0b = ax0.twinx()
gtd1d = 1/transition.powerlaw_cavity(**kwargs_dtg, coord={'R': r_plot, 'z': 0})
ax0b.plot(r_plot/u.au, gtd1d, lw=3, ls='--', color='red')
ax0b.tick_params(axis='y', labelcolor='red')
ax0b.set_yscale('log')
ax0b.set_ylabel('gas-to-dust ratio')

#******************
#2D plotting
#******************
ax1 = fig.add_subplot(224, projection='polar')
theta = np.linspace(0,2*np.pi,100)
r_gas = kwargs_dens['R_cav']/u.au
r_dust = kwargs_dtg['R_cav']/u.au
ax1.plot(theta, np.zeros(len(theta))+r_gas, lw=3, color='green')
ax1.plot(theta, np.zeros(len(theta))+r_dust, lw=3, ls='--', color='red')
ax1.set_xlabel('Radial cavities')
ax1.set_rmax(r_plot[-1]/u.au)

#******************
#2D plotting
#******************
from astropy.io import fits
ax2 = fig.add_subplot(223)
data = fits.getdata('./Subgrids/img_cont.fits')
Tb = ax2.imshow(data.squeeze(), cmap='nipy_spectral', extent = [-70,70,-70,70], norm=colors.Normalize(vmax=47))
ax2.set_xlabel('au')
cbar_dat = plt.colorbar(Tb)
cbar_dat.ax.set_ylabel(r'T$_b$ [K]')

#******************
#3D plotting
#******************
lims = np.array([-70,70])
ax_kw = {'projection': '3d', 'xlim': lims, 'ylim': lims, 'zlim': lims, 'azim': -50, 'elev': 30}
ax3 = fig.add_axes([0.03,0.53,0.43,0.43], **ax_kw)
canvas3d = Plot_model.Canvas3d(fig=fig, ax=ax3)
sp = canvas3d.scatter_random(grid, density/1e6, 1e10, GRID_unit=u.au, power=0, NRand=10000, prop_min=1.0, #function arguments
                             marker = 'o', cmap = 'jet', s = 3, edgecolors = 'none', vmin=1e6, norm = colors.LogNorm()) #Scatter kwargs
cbar = plt.colorbar(sp)
cbar.ax.set_ylabel(r'$n_{gas}$ [cm$^{-3}$]')
ax3.set_xlabel('au')

plt.savefig('transition_disc_3d.png', bbox_inches='tight')
plt.show()
