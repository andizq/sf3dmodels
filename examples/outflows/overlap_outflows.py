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

t0 = time.time()

#********
#GRIDDING
#********
sizex = 500 * u.au
sizey = sizez = 500 * u.au 
Nx = Ny = Nz = 150
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz], rt_code='radmc3d')

#********
#MERGING
#********
files = ['datatab_outflow1.dat', 'datatab_outflow2.dat']
#Old style
#outflows = BGG.overlap(GRID, submodels = data2merge, rho_min = 1e6)

#New style
columns = ['id', 'x', 'y', 'z', 'dens_e', 'dens_ion', 'temp_gas']
outflows = Overlap(GRID)
finalprop = outflows.fromfiles(columns, submodels = files, rt_code = 'radmc3d')

radmc = rt.Radmc3dDefaults(GRID)
radmc.freefree(finalprop)
#********
#TIMING
#********
print ('Ellapsed time: %.3fs' % (time.time() - t0))
print ('-------------------------------------------------\n-------------------------------------------------\n')

#********
#PLOTTING
#********
density = finalprop['dens_e'] / 1e6 #dens. in cm^-3
temperature = finalprop['temp_gas']

weight = 100 * np.mean(density)

#-----------------
#Plot for DENSITY
#-----------------
Pm.scatter3D(GRID, density, weight, NRand = 4000, axisunit = u.au, colorscale = 'log', cmap = 'cool',
             colorlabel = r'${\rm log}_{10}(n [cm^{-3}])$', output = 'global_grid_dens.png', vmin = 5, show=False)

#--------------------
#Plot for TEMPERATURE
#--------------------
Pm.scatter3D(GRID, density, weight, colordim = temperature, NRand = 4000, axisunit = u.au, colorscale = 'log',
             cmap = 'brg', colorlabel = r'${\rm log}_{10}(T$ $[K])$', output = 'global_grid_temp.png', vmin = 2, show=False)
