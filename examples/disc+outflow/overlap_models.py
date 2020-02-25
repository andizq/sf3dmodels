"""
Basic docstring explaining example
"""
from __future__ import print_function
#********************
#sf3dmodels libraries
#********************
from sf3dmodels.outflow import OutflowModel   #Model functions
import sf3dmodels.utils.units as u            #Units
import sf3dmodels.rt as rt                    #Writing functions for radiative transfer
import sf3dmodels.Plot_model as Pm            #Plotting model
import sf3dmodels.Model as Model              #Grid
from sf3dmodels.grid import Overlap           #Overlap submodels
#********************
#Extra libraries
#********************
import numpy as np
import time
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors


t0 = time.time()

#********
#GRIDDING
#********
sizex = 100 * u.au
sizey = sizez = 100 * u.au 
Nx = Ny = Nz = 100
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz], rt_code='lime')

#********
#MERGING
#********
files = ['disc.dat', 'outflow.dat']
#Old style
#outflows = BGG.overlap(GRID, submodels = data2merge, rho_min = 1e6)

#New style
columns = ['id', 'x', 'y', 'z', 'dens_H2', 'dens_Hplus', 'temp_gas', 'vel_x', 'vel_y', 'vel_z', 'abundance', 'gtdratio']
overlap = Overlap(GRID)
finalprop = overlap.fromfiles(columns, submodels = files, rt_code = 'lime')

#**********
#WRITING
#**********
lime = rt.Lime(GRID)
lime.finalmodel(finalprop) 

#********
#TIMING
#********
print ('Ellapsed time: %.3fs' % (time.time() - t0))
print ('-------------------------------------------------\n-------------------------------------------------\n')

#********
#PLOTTING
#********
density = finalprop['dens_H2'] / 1e6 #dens. in cm^-3
temperature = finalprop['temp_gas']

weight = 1.0#100 * np.mean(density)

"""
#-----------------
#Plot for DENSITY
#-----------------
Pm.scatter3D(GRID, density, weight, NRand = 4000, axisunit = u.au, colorscale = 'log', cmap = 'cool',
             colorlabel = r'${\rm log}_{10}(n [cm^{-3}])$', output = 'global_grid_dens.png', vmin = 5)

#--------------------
#Plot for TEMPERATURE
#--------------------
Pm.scatter3D(GRID, density, weight, colordim = temperature, NRand = 4000, axisunit = u.au, colorscale = 'log',
             cmap = 'brg', colorlabel = r'${\rm log}_{10}(T$ $[K])$', output = 'global_grid_temp.png', vmin = 2)
"""

#******************
#3D plotting
#******************
lims = np.array([-100,100])
weight = 1.0 
ax_kw = {'projection': '3d'}#, 'xlim': lims, 'ylim': lims, 'zlim': lims, 'azim': -50, 'elev': 30}
canvas3d = Pm.Canvas3d(ax_kw=ax_kw)
sp = canvas3d.scatter_random(GRID, density, weight, GRID_unit=u.au, power=0, NRand=10000, prop_min=1.0, #function arguments
                             marker = '+', cmap = 'jet', s = 3, edgecolors = 'none', vmin=1, norm = colors.LogNorm()) #Scatter kwargs
cbar = plt.colorbar(sp)
cbar.ax.set_ylabel(r'H$_2$ density [cm$^{-3}$]')
canvas3d.ax.set_xlabel('au')
plt.savefig('grid_dens3d.png', bbox_inches='tight')


canvas3d = Pm.Canvas3d(ax_kw=ax_kw)
sp = canvas3d.scatter_random(GRID, density, weight, prop_color = temperature, GRID_unit=u.au, power=0, NRand=10000, prop_min=1.0, #function arguments
                             marker = '+', cmap = 'jet', s = 3, edgecolors = 'none', vmin=1, norm = colors.LogNorm()) #Scatter kwargs
cbar = plt.colorbar(sp)
cbar.ax.set_ylabel(r'T [K]')
canvas3d.ax.set_xlabel('au')
plt.savefig('grid_temp3d.png', bbox_inches='tight')

plt.show()
