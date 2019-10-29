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
from astropy.io import fits
import time

t0 = time.time()
#******************
#Model parameters
#******************
kwargs_dens={'dn_cav': 1e-5, 'n_cav': 2e17, 'R_cav': 30*u.au, 'power':-0.9}#, 'phi_stddev': 2*np.pi/3, 'phi_mean': -np.pi/2}
kwargs_dtg={'dn_cav': 1e-2, 'n_cav': 0.01, 'R_cav': 40*u.au, 'power': 0, 'phi_stddev': 2*np.pi/3, 'phi_mean': -np.pi/2}#, 'phi_stddev': np.pi/2}
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

fig = plt.figure(figsize=(11,12))
R_plot = np.linspace(u.au, R_disc, 100)
R_const = np.zeros(R_plot.shape) 
R_plot_0 = 45*u.au

phi = np.linspace(0,2*np.pi,100)
phi_const = np.zeros(R_plot.shape) 
phi_sine = 0.1*np.pi * np.sin(2*np.pi*R_plot/(0.5*R_disc)) + np.pi
phi_plot_0 = np.pi/3
phi_plot_1 = -np.pi/2

z_plot = 0
#******************
#2D plotting
#******************
ax1 = fig.add_subplot(324, projection='polar')

ax1.plot(phi_const+phi_plot_0, R_plot/u.au, lw=3, ls='--', color='blue',  alpha=1.0)
ax1.plot(phi_const+phi_plot_1, R_plot/u.au, lw=3, ls='-', color='blue', alpha=0.7)
ax1.plot(phi_sine, R_plot/u.au, lw=3, ls=':', color='blue', alpha=0.7)
ax1.plot(phi_sine, R_plot/u.au, lw=3, ls=':', color='blue', alpha=0.7)
ax1.plot(phi, R_const+R_plot_0/u.au, lw=3, ls='-', color='darkorange', alpha=1.0)

ax1.set_ylabel('Disc polar map', labelpad=30)
ax1.set_rmax(R_plot[-1]/u.au)

ax1.text(1.2, 0.0, 'Paths for Radial 1d profiles', transform=ax1.transAxes, rotation=90,
         size=SMALL_SIZE, ha='left', va='bottom', color='blue', fontweight='bold')
ax1.text(1.3, 0.0, 'Paths for Azimuthal 1d profiles', transform=ax1.transAxes, rotation=90,
         size=SMALL_SIZE, ha='left', va='bottom', color='darkorange', fontweight='bold')

#******************
#1D plotting
#******************
ax0 = fig.add_subplot(322)
dens_0 = transition.powerlaw_cavity(**kwargs_dens, coord={'R': R_plot, 'z': z_plot, 'phi': phi_plot_0})
dens_1 = transition.powerlaw_cavity(**kwargs_dens, coord={'R': R_plot, 'z': z_plot, 'phi': phi_plot_1})
dens_2 = transition.powerlaw_cavity(**kwargs_dens, coord={'R': R_plot, 'z': z_plot, 'phi': phi_sine})
ax0.plot(R_plot/u.au, dens_0/1e6, lw=3, ls='--', color='green')
ax0.plot(R_plot/u.au, dens_1/1e6, lw=3, ls='-', color='green', alpha=0.7)
ax0.plot(R_plot/u.au, dens_2/1e6, lw=3, ls=':', color='green', alpha=0.7)
ax0.set_xlabel('au')
ax0.xaxis.set_label_position('top') 
ax0.xaxis.set_minor_locator(AutoMinorLocator())
ax0.set_yscale('log')
ax0.set_ylabel(r'H$_2$ density [cm$^{-3}$]')
ax0.tick_params(axis='y', labelcolor='green')
ax0.tick_params(axis='x', labelbottom=False, labeltop=True, top=True) #By default bottom=True
ax0.tick_params(axis='x', which='minor', top=True)

ax0b = ax0.twinx()
gtd_0 = 1/transition.powerlaw_cavity(**kwargs_dtg, coord={'R': R_plot, 'z': z_plot, 'phi': phi_plot_0})
gtd_1 = 1/transition.powerlaw_cavity(**kwargs_dtg, coord={'R': R_plot, 'z': z_plot, 'phi': phi_plot_1})
gtd_2 = 1/transition.powerlaw_cavity(**kwargs_dtg, coord={'R': R_plot, 'z': z_plot, 'phi': phi_sine})
ax0b.plot(R_plot/u.au, gtd_0, lw=3, ls='--', color='red')
ax0b.plot(R_plot/u.au, gtd_1, lw=3, ls='-', color='red', alpha=0.7)
ax0b.plot(R_plot/u.au, gtd_2, lw=3, ls=':', color='red')
ax0b.tick_params(axis='y', labelcolor='red')
ax0b.set_yscale('log')
ax0b.set_ylabel('gas-to-dust ratio')

#******************
#1D plotting
#******************
ax2 = fig.add_subplot(326)
gtd_0 = 1/transition.powerlaw_cavity(**kwargs_dtg, coord={'R': R_plot_0, 'z': z_plot, 'phi': phi})
ax2.plot(phi*180/np.pi, gtd_0, lw=3, ls='-', color='red')
ax2.tick_params(axis='y', which='both', labelcolor='red', left=True, labelleft=False, right=True, labelright=True)
ax2.set_ylabel('gas-to-dust ratio')
ax2.set_xlabel('deg.')
#ax2.set_yscale('log')
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_label_position('right') 

#******************
#2D plotting
#******************
ax3 = fig.add_subplot(323)
data = fits.getdata('img_cont_faceon_band7.fits')
Tb = ax3.imshow(data.squeeze(), cmap='nipy_spectral', extent = [-70,70,-70,70], origin='lower left', norm=colors.Normalize(vmax=47))
ax3.set_xlabel('au')
ax3.set_ylabel('face-on')
cbar_dat = plt.colorbar(Tb)
cbar_dat.ax.set_ylabel(r'T$_b$ [K] - Band 7')

#******************
#2D plotting
#******************
ax4 = fig.add_subplot(325)
data = fits.getdata('img_cont_edgeon_phi90_band7.fits')
Tb = ax4.imshow(data.squeeze(), cmap='nipy_spectral', extent = [-70,70,-70,70], origin='lower left')#, norm=colors.Normalize(vmax=47))
ax4.set_xlabel('au')
ax4.set_ylabel(r'edge-on $\phi$=90$^o$')
cbar_dat = plt.colorbar(Tb)
cbar_dat.ax.set_ylabel(r'T$_b$ [K] - Band 7')

#******************
#3D plotting
#******************
lims = np.array([-70,70])
weight = 1e10 
ax_kw = {'projection': '3d', 'xlim': lims, 'ylim': lims, 'zlim': lims, 'azim': -50, 'elev': 30}
ax5 = fig.add_axes([0.03,0.65,0.43,0.3], **ax_kw)
canvas3d = Plot_model.Canvas3d(fig=fig, ax=ax5)
sp = canvas3d.scatter_random(grid, density/1e6, weight, GRID_unit=u.au, power=0, NRand=10000, prop_min=1.0, #function arguments
                             marker = '+', cmap = 'jet', s = 3, edgecolors = 'none', vmin=1e6, norm = colors.LogNorm()) #Scatter kwargs
cbar = plt.colorbar(sp)
#cbar.ax.set_ylabel(r'H$_2$ density [cm$^{-3}$]')
ax5.set_xlabel('au')

plt.savefig('transition_disc_3d.png', bbox_inches='tight')
plt.show()
