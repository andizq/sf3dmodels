"""
Building global model from 2 submodels
======================================

Basic docstring explaining example
"""
#------------------
#Import the package
#------------------
from sf3dmodels import Model, Plot_model as Pm
import sf3dmodels.utils.units as u
import sf3dmodels.rt as rt
from sf3dmodels.grid import Overlap
#-----------------
#Extra libraries
#-----------------
import numpy as np

#********
#GRIDDING
#********
sizex = sizey = sizez = 1000 * u.au
Nx = Ny = Nz = 200
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz])

#**********
#OVERLAPING
#**********
columns = ['id', 'x', 'y', 'z', 'dens_H2', 'temp_gas', 'vel_x', 'vel_y', 'vel_z', 'abundance_0', 'gtdratio']
data2merge = ['Main.dat', 'Burger.dat']
overlap = Overlap(GRID)
finalprop = overlap.fromfiles(columns, submodels = data2merge)

#**********
#WRITING
#**********
lime = rt.Lime(GRID)
lime.finalmodel(finalprop) 

#--------
#PLOTTING
#--------

density = finalprop['dens_H2'] / 1e6 #1e6 to convert from m^-3 to cm^-3
temperature = finalprop['temp_gas']

weight = 400 * np.mean(density)

#-----------------
#Plot for DENSITY
#-----------------
Pm.scatter3D(GRID, density, weight, NRand = 7000, axisunit = u.au, colorscale = 'log', cmap = 'hot',
  	     colorlabel = r'${\rm log}_{10}(\rho [cm^{-3}])$', output = 'global_grid_dens.png')

#--------------------
#Plot for TEMPERATURE
#--------------------
Pm.scatter3D(GRID, density, weight, colordim = temperature, NRand = 7000, axisunit = u.au, colorscale = 'log',
             cmap = 'brg', colorlabel = r'${\rm log}_{10}(T$ $[K])$', output = 'global_grid_temp.png')
