import BuildGlobalGrid as BGG
import Model
import Plot_model as Pm
import Utils as U
import numpy as np

sizex = sizey = sizez = 1000 * U.AU
Nx = Ny = Nz = 120
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz])
list_sub = ['datatab_Main.dat', 'datatab_Burger.dat']
global_prop = BGG.overlap(GRID, submodels = list_sub)

#--------
#PLOTTING
#--------

GRID = global_prop.GRID 
density = global_prop.density / 1e6 #1e6 to convert from m^-3 to cm^-3
temperature = global_prop.temperature

weight = 400 * np.mean(density)

#-----------------
#Plot for DENSITY
#-----------------
Pm.scatter3D(GRID, density, weight, NRand = 7000, axisunit = U.AU, colorscale = 'log', palette = 'hot',
  	     colorlabel = r'${\rm log}_{10}(\rho [cm^{-3}])$', output = 'global_grid_dens.png')

#--------------------
#Plot for TEMPERATURE
#--------------------
Pm.scatter3D(GRID, density, weight, colordim = temperature, NRand = 7000, axisunit = U.AU, colorscale = 'log',
             palette = 'brg', colorlabel = r'${\rm log}_{10}(T$ $[K])$', output = 'global_grid_temp.png')
