"""
Basic docstring explaining example
"""
#------------------
#Import the package
#------------------
from sf3dmodels import *
#-----------------
#Extra libraries
#-----------------
from matplotlib import colors
import numpy as np
import time

t0 = time.time()

#------------------
#General Parameters
#------------------
#from Galvan-Madrid et al. 2009, Table 3:

MStar = 34 * U.MSun
r_max = 2530 * U.AU #1000 * U.AU #H II sphere size
r_min = r_max / 200 #Minimum distance (!= 0 to avoid indeterminations)
rho_s = 2 * 1.5e6 * 1e6 / 5 #from cgs to SI. Density at sonic radius
q = 1.3 #Density powerlaw
t_e = 1.e4 #K

#-------------------------------
#Parameters for the Pringle disc
#-------------------------------
MRate = 3e-4 * U.MSun_yr
RStar = U.RSun * ( MStar/U.MSun )**0.8
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

#--------
#ENVELOPE
#--------
densEnv = Model.density_Keto_HII(MStar, r_min, r_max, rho_s, t_e, GRID, q = 1.5)

#-------
#DISC
#-------
Rd = 10*densEnv.rs #10 times the sonic radius, just to make it visible
Rho0 = Res.Rho0(MRate, Rd, MStar)
Arho = 2 * 2.0 / 80 # 60.0
densDisc = Model.density_Env_Disc(RStar, Rd, Rho0, Arho, GRID, discFlag = True, envFlag = False, 
                                  rdisc_max = Rd)

density = Model.Struct( **{ 'total': densEnv.total + densDisc.total,
                            'disc': densDisc.total, 
                            'env': densEnv.total,
                            'discFlag': True,
                            'envFlag': True,
                            'r_disc': densDisc.r_disc, 
                            'r_env': densEnv.r_env} )

temperature = Model.temperature_Constant(density, GRID, discTemp=t_e, envTemp=t_e, backTemp=2.725)

Model.PrintProperties(density, temperature, GRID)

#-------
#TIMING
#-------
print ('Ellapsed time: %.3fs' % (time.time() - t0))
print ('-------------------------------------------------\n-------------------------------------------------\n')

#------------------------------
#WRITING DATA in RADMC3D FORMAT
#------------------------------

Model.Datatab_RADMC3D_FreeFree(density.total, temperature.total, GRID)

#------------------------------------
#3D PLOTTING (weighting with density)
#------------------------------------
tag = 'keto+disc_HII'
weight = rho_s

norm = colors.LogNorm()
Plot_model.scatter3D(GRID, density.total, weight, NRand = 4000, colordim = density.total / 1e6, axisunit = U.AU, cmap = 'jet', 
                     colorscale = 'log', colorlabel = r'${\rm log}_{10}$($n_{\rm e}$ [cm$^{-3}$])', output = '3Ddens_%s.png'%tag, show = True)

Plot_model.scatter3D(GRID, density.total, weight, NRand = 4000, colordim = temperature.total, axisunit = U.AU, cmap = 'jet', marker = 'o',
                     colorscale = 'uniform', colorlabel = r'$T_{\rm e}$ [Kelvin]', output = '3Dtemp_%s.png'%tag, show = True)


