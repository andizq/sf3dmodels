"""
Powerlaw Sphere Example
=======================

This is an example of a powerlaw-density  spherical HII region
"""
#------------------
#Import the package
#------------------
from sf3dmodels import Model, Plot_model
import sf3dmodels.utils.units as u
import sf3dmodels.rt as rt
#-----------------
#Extra libraries
#-----------------
import numpy as np
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import time

t0 = time.time()

#------------------
#General Parameters
#------------------
#from Galvan-Madrid et al. 2009, Table 3:

MStar = 34 * u.MSun
r_max = 2530 * u.au #H II sphere size
r_min = r_max / 200 #Minimum distance (!= 0 to avoid indeterminations)
r_s = r_max #Normalization distance
rho_s = 1.4e5 * 1e6 #from cgs to SI. Density at r_s
q = 1.3 #Density powerlaw
t_e = 1.e4 #K

#---------------
#GRID Definition
#---------------

sizex = sizey = sizez = 2600 * u.au 
Nx = Ny = Nz = 63 #Number of divisions for each axis
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz], rt_code = 'radmc3d')
NPoints = GRID.NPoints #Final number of nodes in the grid

#-------------------
#PHYSICAL PROPERTIES
#-------------------
density = Model.density_Powerlaw_HII(r_min, r_max, r_s, rho_s, q, GRID)
temperature = Model.temperature_Constant(density, GRID, envTemp = t_e, backTemp=2.725)

Model.PrintProperties(density, temperature, GRID)

#-------
#TIMING
#-------
print ('Ellapsed time: %.3fs' % (time.time() - t0))
print ('-------------------------------------------------\n-------------------------------------------------\n')

#----------------------
#WRITING RADMC-3D FILES
#----------------------
prop = {'dens_e': density.total,
        'dens_ion': density.total,
        'temp_gas': temperature.total}

Rad = rt.Radmc3dDefaults(GRID)
Rad.freefree(prop)

#------------------------------------
#3D PLOTTING (weighting with density)
#------------------------------------
tag = 'plsphere_HII'
weight = 10*rho_s

Plot_model.scatter3D(GRID, density.total, weight, NRand = 4000, colordim = density.total / 1e6, axisunit = u.au, cmap = 'jet', 
                     colorscale = 'log', colorlabel = r'${\rm log}_{10}$($n_{\rm e}$ [cm$^{-3}$])', output = '3Ddens_%s.png'%tag, show = True)

Plot_model.scatter3D(GRID, density.total, weight, NRand = 4000, colordim = temperature.total, axisunit = u.au, cmap = 'binary', marker = 'o', s = 4,
                     colorscale = 'uniform', colorlabel = r'$T_{\rm e}$ [Kelvin]', output = '3Dtemp_%s.png'%tag, show = True)



fig = plt.figure(figsize=(8,6))
#ax = plt.axes(projection='3d')
lims = (-3000,3000)
canvas3d = Plot_model.Canvas3d(fig=fig, axes_kw={'xlim': lims, 'ylim': lims, 'zlim': lims, 'azim': 50, 'elev': 10})
ax = canvas3d.ax #generated with fig.add_axes from matplotlib, hence all the matplotlib functions are available
sp = canvas3d.scatter_random(GRID, density.total/1e6, weight/1e6, GRID_norm=u.au, power=0., NRand=4000, prop_color=temperature.total, prop_min=3.0, #function arguments
                             marker = 'o', cmap = 'jet', s = 4, edgecolors = 'none', vmin = None, vmax = None, norm = colors.LogNorm())#marker='o', s=10, cmap='hot') #Scatter kwargs
#ax.scatter(3000,0,0, c=[dens_e*1.2/1e6], norm=colors.Normalize(density.total.min()/1e6, vmax=density.total.max()/1e6), cmap=sp.cmap)

cbar = plt.colorbar(sp)
cbar.ax.set_ylabel(r'$n_{e^-}$ [cm$^{-3}$]')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

#plt.tight_layout()
plt.savefig('3D_%s.png'%tag, dpi = 500, bbox_inches='tight')
plt.show()



