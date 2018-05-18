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
dens_e = 1.4e5 * 1e6 #Electronic numerical density, from cgs to SI
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
temperature = Model.temperature_Constant(density, GRID, envTemp = t_e, backTemp=2.725)

Model.PrintProperties(density, temperature, GRID) #Printing resultant properties (mass, mean temperature, etc)

#---------------------------------
#WRITING DATA with RADMC-3D FORMAT
#---------------------------------
Model.Datatab_RADMC3D_FreeFree(density.total, temperature.total, GRID)

#------------------------------------
#3D PLOTTING (weighting with density)
#------------------------------------
tag = 'ctsphere_HII'
weight = dens_e
#Plot_model.scatter3D(GRID, density.total, weight, NRand = 4000, colordim = density.total / 1e6, axisunit = U.AU, palette = 'jet', 
#                     colorscale = 'log', colorlabel = r'$n_{\rm e}$ [cm$^{-3}$]', output = '%s.png'%tag, show = True)
