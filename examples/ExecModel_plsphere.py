#------------------
#Import the package
#------------------
from sf3dmodels import *
#-----------------
#Extra libraries
#-----------------
import numpy as np
import time

t0 = time.time()

#------------------
#General Parameters
#------------------
#from Galvan-Madrid et al. 2009, Table 3:

MStar = 34 * U.MSun
r_max = 2530 * U.AU #1000 * U.AU #H II sphere size
r_min = r_max / 200
r_s = r_max
rho_s = 1.4e5 * 1e6 #from cgs to SI
q = 1.3 #Powerlaw
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
density = Model.density_Powerlaw_HII(r_min, r_max, r_s, rho_s, q, GRID)
temperature = Model.temperature_Constant(density, GRID, envTemp = t_e, backTemp=2.725)

Model.PrintProperties(density, temperature, GRID)

#-------
#TIMING
#-------
print ('Ellapsed time: %.3fs' % (time.time() - t0))
print ('-------------------------------------------------\n-------------------------------------------------\n')

#---------------------------------
#WRITING DATA with RADMC-3D FORMAT
#---------------------------------

Model.Datatab_RADMC3D_FreeFree(density.total, temperature.total, GRID)

#------------------------------------
#3D PLOTTING (weighting with density)
#------------------------------------
tag = 'plsphere_HII'
weight = 10*rho_s
Plot_model.scatter3D(GRID, density.total, weight, NRand = 4000, colordim = density.total / 1e6, axisunit = U.AU, palette = 'jet', 
                     colorscale = 'log', colorlabel = r'$n_{\rm e}$ [cm$^{-3}$]', output = '%s.png'%tag, show = True)
