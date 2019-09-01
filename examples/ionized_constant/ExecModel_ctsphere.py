"""
Constant Sphere Example
=======================

This is an exmaple of a uniform-density (constant-density) spherical HII region
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
import time

#------------------
#General Parameters
#------------------
r_max = 2530 * u.au #H II sphere size
dens_e = 1.4e5 * 1e6 #Electron number density, from cgs to SI
t_e = 1.e4 #K

#---------------
#GRID Definition
#---------------
sizex = sizey = sizez = 2600 * u.au 
Nx = Ny = Nz = 63 #Number of divisions for each axis
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz], rt_code='radmc3d')
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
prop = {'dens_e': density.total,
        'dens_ion': density.total,
        'temp_gas': temperature.total}

Rad = rt.Radmc3dDefaults(GRID)
Rad.freefree(prop)

#------------------------------------
#3D PLOTTING (weighting with density)
#------------------------------------
tag = 'ctsphere_HII'
weight = dens_e

Plot_model.scatter3D(GRID, density.total, weight, power=0, NRand = 4000, colordim = density.total / 1e6 / 1e5, axisunit = u.au, cmap = 'winter', 
                     marker = 'o', colorlabel = r'$n_{\rm e}$ [cm$^{-3}$] x $10^5$', output = '3Ddens_%s.png'%tag, show = True)

Plot_model.scatter3D(GRID, density.total, weight, power=0, NRand = 4000, colordim = temperature.total, axisunit = u.au, cmap = 'winter', 
                     marker = 'o', colorlabel = r'$T_{\rm e}$ [Kelvin]', output = '3Dtemp_%s.png'%tag, show = True)
