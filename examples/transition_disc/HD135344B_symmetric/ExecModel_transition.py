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
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import AutoMinorLocator
import time

t0 = time.time()
#******************
#Model parameters
#******************
kwargs_dens={'dn_cav': 1e-5, 'n_cav': 2e17, 'R_cav': 30*u.au, 'power':-0.9}
kwargs_dtg={'dn_cav': 1e-2, 'n_cav': 0.01, 'R_cav': 40*u.au, 'power': 0}
R_disc = 70*u.au

#******************
#Grid creation
#******************
transition = disc.Transition()
init_grid = Grid()
grid = init_grid.random(func=transition.powerlaw_cavity, power=0.4, r_size=R_disc, kwargs_func=kwargs_dens, normalization=1e16, npoints=100000)

#*****************************************
#Computing relevant physical distributions
#*****************************************
density = transition.powerlaw_cavity(grid=grid, **kwargs_dens)
gtdratio = 1./transition.powerlaw_cavity(grid=grid, **kwargs_dtg) #Converting to gas-to-dust ratio for LIME
temperature = np.zeros(grid.NPoints) + 100. #Constant temperature 100 K

prop = {'dens_H2': density,
        'gtdratio': gtdratio,
        'temp_gas': temperature}

#**********************
#Writing files for LIME 
#**********************
lime = rt.Lime(grid)
lime.submodel(prop, output='datatab.dat', folder='./', lime_npoints = True, lime_header = True)
print('Output columns', lime.columns)

print ('Ellapsed time:', time.time()-t0)

#******************
#PLOTTING BLOCK
#******************
SMALL_SIZE = 11
MEDIUM_SIZE = 15
BIGGER_SIZE = 18
matplotlib.rc('axes', labelsize=MEDIUM_SIZE)
matplotlib.rc('xtick', labelsize=SMALL_SIZE) 
matplotlib.rc('ytick', labelsize=SMALL_SIZE)

fig = plt.figure(figsize=(11,8))
R_plot = np.arange(u.au, R_disc, u.au)
z_plot = 0
phi_plot = np.pi/4

#******************
#1D plotting
#******************
ax0 = fig.add_subplot(222)
dens1d = transition.powerlaw_cavity(**kwargs_dens, coord={'R': R_plot, 'z': z_plot, 'phi': phi_plot})
ax0.plot(R_plot/u.au, dens1d/1e6, lw=3, color='green')
#ax0.set_xscale('log')
ax0.set_xlabel('au')
ax0.xaxis.set_label_position('top') 
ax0.xaxis.set_minor_locator(AutoMinorLocator())
ax0.set_yscale('log')
ax0.set_ylabel(r'H$_2$ density [cm$^{-3}$]')
ax0.tick_params(axis='y', labelcolor='green')
ax0.tick_params(axis='x', labelbottom=False, labeltop=True, top=True) #By default bottom=True
ax0.tick_params(axis='x', which='minor', top=True)

ax0b = ax0.twinx()
gtd1d = 1/transition.powerlaw_cavity(**kwargs_dtg, coord={'R': R_plot, 'z': z_plot, 'phi': phi_plot})
ax0b.plot(R_plot/u.au, gtd1d, lw=3, ls='--', color='red')
ax0b.tick_params(axis='y', labelcolor='red')
ax0b.set_yscale('log')
ax0b.set_ylabel('gas-to-dust ratio')

#******************
#2D plotting
#******************
ax1 = fig.add_subplot(224, projection='polar')
theta = np.linspace(0,2*np.pi,100)
R_gas = kwargs_dens['R_cav']/u.au
R_dust = kwargs_dtg['R_cav']/u.au
ax1.plot(theta, np.zeros(len(theta))+R_gas, lw=3, color='green')
ax1.plot(theta, np.zeros(len(theta))+R_dust, lw=3, ls='--', color='red')
ax1.plot(np.zeros(R_plot.shape)+phi_plot, R_plot/u.au, lw=3, ls='-', color='blue')
ax1.set_xlabel('Radial cavities')
ax1.set_rmax(R_plot[-1]/u.au)

#******************
#2D plotting
#******************
from astropy.io import fits
ax2 = fig.add_subplot(223)
data = fits.getdata('img_cont_faceon.fits')
Tb = ax2.imshow(data.squeeze(), cmap='nipy_spectral', extent = [-70,70,-70,70], origin='lower left', norm=colors.Normalize(vmax=47))
ax2.set_xlabel('au')
cbar_dat = plt.colorbar(Tb)
cbar_dat.ax.set_ylabel(r'T$_b$ [K]')

#******************
#3D plotting
#******************
lims = np.array([-70,70])
weight = 1e10 
ax_kw = {'projection': '3d', 'xlim': lims, 'ylim': lims, 'zlim': lims, 'azim': -50, 'elev': 30}
ax3 = fig.add_axes([0.03,0.53,0.43,0.43], **ax_kw)
canvas3d = Plot_model.Canvas3d(fig=fig, ax=ax3)
sp = canvas3d.scatter_random(grid, density/1e6, weight, GRID_unit=u.au, power=0, NRand=10000, prop_min=1.0, #function arguments
                             marker = '+', cmap = 'jet', s = 3, edgecolors = 'none', vmin=1e6, norm = colors.LogNorm()) #Scatter kwargs
cbar = plt.colorbar(sp)
#cbar.ax.set_ylabel(r'H$_2$ density [cm$^{-3}$]')
ax3.set_xlabel('au')

plt.savefig('transition_disc_3d.png', bbox_inches='tight')
plt.show()
