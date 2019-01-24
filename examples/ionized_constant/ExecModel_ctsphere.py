"""
Constant Sphere Example
=======================

This is an exmaple of a uniform-density (constant-density) spherical HII region
"""
#------------------
#Import the package
#------------------
from sf3dmodels import *
#-----------------
#Extra libraries
#-----------------
import numpy as np
import time

#------------------
#General Parameters
#------------------
r_max = 2530 * U.AU #H II sphere size
dens_e = 1.4e5 * 1e6 #Electron number density, from cgs to SI
t_e = 1.e4 #K

#---------------
#GRID Definition
#---------------
sizex = sizey = sizez = 2600 * U.AU 
Nx = Ny = Nz = 63 #Number of divisions for each axis
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz], radmc3d = True)
NPoints = GRID.NPoints #Final number of nodes in the grid

#-------------------
#PHYSICAL PROPERTIES
#-------------------
density = Model.density_Constant(r_max, GRID, envDens = dens_e)
temperature = Model.temperature_Constant(density, GRID, envTemp = t_e, backTemp = 2.725)

Model.PrintProperties(density, temperature, GRID) #Printing resultant properties (mass, mean temperature, etc)

#----------------------
#WRITING RADMC-3D FILES
#----------------------
Rad = Model.Radmc3dRT(GRID)
prop = {'dens_elect': density.total,
        'dens_ion': density.total,
        'tgas': temperature.total}
Rad.freefree(prop)

#------------------------------------
#3D PLOTTING (weighting with density)
#------------------------------------
tag = 'ctsphere_HII'
weight = dens_e

Plot_model.scatter3D(GRID, density.total, weight, NRand = 4000, colordim = density.total / 1e6 / 1e5, axisunit = U.AU, cmap = 'winter', 
                     marker = 'o', colorlabel = r'$n_{\rm e}$ [cm$^{-3}$] x $10^5$', output = '3Ddens_%s.png'%tag, show = True)

Plot_model.scatter3D(GRID, density.total, weight, NRand = 4000, colordim = temperature.total, axisunit = U.AU, cmap = 'winter', 
                     marker = 'o', colorlabel = r'$T_{\rm e}$ [Kelvin]', output = '3Dtemp_%s.png'%tag, show = True)
