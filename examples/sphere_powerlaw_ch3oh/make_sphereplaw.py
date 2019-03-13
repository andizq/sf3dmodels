"""
Powerlaw Sphere Example
=======================

This is an example of a spherical powerlaw density and temperature profiles.
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
import os 

t0 = time.time()

#------------------
#General Parameters
#------------------
r_max = 2000*u.au
r_min = r_max / 256 #Minimum distance (!= 0 to avoid indeterminations)

rho_mean = 1.e5*1e6 #Mean number density in part/m3
q = -2.0 #Density powerlaw

T10Env = 1000
Tmin_env = 25. #lower threshold for envelope temperature
p = 0.5 #temperature powerlaw

#The following parameters are only relevant when disc component is invoked.

MStar = None #30.0 * U.MSun 
MRate = None #4e-4 * U.MSun_yr 

RStar = None #26 * U.RSun * ( MStar/U.MSun )**0.27 * ( MRate / (1e-3*U.MSun_yr) )**0.41 
LStar = None #8.6e4 * U.LSun # 4
TStar = None #U.TSun * ( (LStar/U.LSun) / (RStar/U.RSun)**2 )**0.25

BT = None
Rd = None #centrifugal radius
Tmin_disc = None

#---------------
#GRID Definition
#---------------

sizex = sizey = sizez = r_max 
Nx = Ny = Nz = 128 #Number of divisions for each axis 
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz], rt_code = 'lime')
NPoints = GRID.NPoints #Final number of nodes in the grid

#-------------------
#PHYSICAL PROPERTIES
#-------------------
density = Model.density_Powerlaw(r_max, rho_mean, q, GRID, rho_min = 1.0e4)

temperature = Model.temperature(TStar, Rd, T10Env, RStar, MStar, MRate, BT, density, GRID, 
                                Tmin_disc = Tmin_disc, Tmin_env = Tmin_env , p = p, ang_cavity = False)

#--------
#VELOCITY
#--------
vel = Model.velocity_random(1800,GRID.NPoints)

#-------------------------------
#ABUNDANCE and GAS-to-DUST RATIO
#-------------------------------
ab0 = 1.e-8 #CH3OH abundance
abundance = Model.abundance(ab0, NPoints)

gtd0 = 100. #Gas to dust ratio
gtdratio = Model.gastodust(gtd0, NPoints)

#--------------------
#PRINTING and WRITING
#--------------------
prop = {'dens_p_H2': density.total,
        'temp_gas': temperature.total,
        'vel_x': vel.x,
        'vel_y': vel.y,
        'vel_z': vel.z,
        'abundance_0': abundance,
        'gtdratio': gtdratio}
lime = rt.Lime(GRID)
lime.finalmodel(prop)

Model.PrintProperties(density, temperature, GRID)                
#-------
#TIMING
#-------
print ('Ellapsed time: %.3fs' % (time.time() - t0))
print ('-------------------------------------------------\n-------------------------------------------------\n')
