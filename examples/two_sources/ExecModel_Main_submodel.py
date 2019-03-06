"""
Basic docstring explaining example
"""
from __future__ import print_function
#------------------
#Import the package
#------------------
from sf3dmodels import Model, Plot_model
from sf3dmodels import Resolution as Res
import sf3dmodels.utils.units as u            
import sf3dmodels.rt as rt                   
#-----------------
#Extra libraries
#-----------------
from matplotlib import colors
import numpy as np
import os
import time

t0 = time.time()

#------------------
#General Parameters
#------------------
MStar = 7.0 * u.MSun 
MRate = 4e-4 * u.MSun_yr 
RStar = 26 * u.RSun * ( MStar/u.MSun )**0.27 * ( MRate / (1e-3*u.MSun_yr) )**0.41 
LStar = 3.2e4 * u.LSun  
TStar = u.TSun * ( (LStar/u.LSun) / (RStar/u.RSun)**2 )**0.25
Rd = 152. * u.au

print ('RStar:', RStar/u.RSun,', LStar:', LStar/u.LSun, ', TStar:', TStar)

#---------------
#GRID Definition
#---------------
#Cubic grid, each edge ranges [-500, 500] au.

sizex = sizey = sizez = 500 * u.au
Nx = Ny = Nz = 100 #Number of divisions for each axis
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz])
NPoints = GRID.NPoints #Number of nodes in the grid

#--------
#DENSITY
#--------
Rho0 = Res.Rho0(MRate, Rd, MStar)
Arho = 24.1 #Disc-envelope density factor
Renv = 500 * u.au #Envelope radius
Cavity = 40 * np.pi/180 #Cavity opening angle 
density = Model.density_Env_Disc(RStar, Rd, Rho0, Arho, GRID, discFlag = True, envFlag = True, 
                                 renv_max = Renv, ang_cavity = Cavity)

#-----------
#TEMPERATURE
#-----------
T10Env = 375. #Envelope temperature at 10 au
BT = 5. #Adjustable factor for disc temperature. Extra, or less, disc heating.
temperature = Model.temperature(TStar, Rd, T10Env, RStar, MStar, MRate, BT, density, GRID,
                                ang_cavity = Cavity)

#--------
#VELOCITY
#--------
vel = Model.velocity(RStar, MStar, Rd, density, GRID)

#-------------------------------
#ABUNDANCE and GAS-to-DUST RATIO
#-------------------------------
ab0 = 1.8e-7 #CH3CN abundance
abundance = Model.abundance(ab0, NPoints)

gtd0 = 100. #Gas to dust ratio
gtdratio = Model.gastodust(gtd0, NPoints)

#-------------------------
#ROTATION, VSYS, CENTERING
#-------------------------
xc, yc, zc = [-250*u.au, 350*u.au, 300*u.au]
CENTER = [xc, yc, zc] #New center of the region in the global grid
v_sys = 3320. #Systemic velocity (vz) of the region (in m/s)
newProperties = Model.ChangeGeometry(GRID, center = CENTER, vsys = v_sys,  vel = vel,
	      	 	             rot_dict = { 'angles': [np.pi/4, 1.87*np.pi], 'axis': ['x','z'] })

#At the minute, the Model library only modifies the XYZ lists. 
 #This is enough information for LIME
GRID.XYZ = newProperties.newXYZ #Redefinition of the XYZ grid

#The velocity should inherit the new velocity distribution 
 #(because we rotated the system and added a systemic velocity)
vel.x, vel.y, vel.z = newProperties.newVEL 

#-----------------------------
#WRITING DATA for LIME
#-----------------------------
tag = 'Main.dat' 
prop = {'dens_H2': density.total,
        'temp_gas': temperature.total,
        'vel_x': vel.x,
        'vel_y': vel.y,
        'vel_z': vel.z,
        'abundance': abundance,
        'gtdratio': gtdratio}
lime = rt.Lime(GRID)
lime.submodel(prop, output=tag)
print('Output columns', lime.columns)

#-----------------------------
#PRINTING resultant PROPERTIES
#-----------------------------
Model.PrintProperties(density, temperature, GRID)

#-------
#TIMING
#-------
print ('Ellapsed time: %.3fs' % (time.time() - t0))
print ('-------------------------------------------------\n-------------------------------------------------\n')

#------------------------------------
#3D PLOTTING (weighting with density)
#------------------------------------
tag = 'Main'
weight = 10*Rho0
r = GRID.rRTP[0] / u.au #GRID.rRTP hosts [r, R, Theta, Phi] --> Polar GRID
Plot_model.scatter3D(GRID, density.total, weight, NRand = 4000, colordim = r, axisunit = u.au, cmap = 'jet', 
                     colorscale = 'log', colorlabel = r'${\rm log}_{10}(r$ $[au])$', output = '3Dpoints%s.png'%tag, show = True)
