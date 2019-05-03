"""
Global model from 2 submodels
=============================

Example
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
sizex = sizey = sizez = 0.3 * u.pc
Nx = Ny = Nz = 100
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz])

#**********
#OVERLAPING
#**********
columns = ['id', 'x', 'y', 'z', 'dens_H', 'temp_gas', 'vel_x', 'vel_y', 'vel_z', 'abundance', 'gtdratio']
data2merge = ['filament.dat', 'envelope.dat']
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
density = finalprop['dens_H'] / 1e6 #1e6 to convert from m^-3 to cm^-3
temperature = finalprop['temp_gas']

weight = 50*np.mean(density)

#-----------------
#Plot for DENSITY
#-----------------
from matplotlib.colors import LogNorm
lims=np.array([-0.3,0.3])*u.pc

Pm.scatter3D(GRID, density, weight, NRand = 7000, axisunit = u.pc, 
             norm=LogNorm(vmin=1e4, vmax=1e7), cmap = 'nipy_spectral',
             colorlabel = r'${\rm log}_{10}(\rho [cm^{-3}])$', output = 'global_grid_dens.png',
             xlim=lims, ylim=lims, zlim=lims)

#--------------------
#Plot for TEMPERATURE
#--------------------
Pm.scatter3D(GRID, density, weight, colordim = temperature, NRand = 7000, axisunit = u.au, colorscale = 'log',
             cmap = 'brg', colorlabel = r'${\rm log}_{10}(T$ $[K])$', output = 'global_grid_temp.png',
             xlim=lims, ylim=lims, zlim=lims)
