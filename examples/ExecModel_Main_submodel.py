from __future__ import print_function
#-----------------
#Package libraries
#-----------------
import Model
import Resolution as Res
import Plot_model
import Utils as U
#-----------------
#Extra libraries
#-----------------
import numpy as np
import os
import time

t0 = time.time()

#------------------
#General Parameters
#------------------
MStar = 7.0 * U.MSun 
MRate = 4e-4 * U.MSun_yr 
RStar = 26 * U.RSun * ( MStar/U.MSun )**0.27 * ( MRate / (1e-3*U.MSun_yr) )**0.41 
LStar = 3.2e4 * U.LSun  
TStar = U.TSun * ( (LStar/U.LSun) / (RStar/U.RSun)**2 )**0.25
Rd = 152. * U.AU

print ('RStar:', RStar/U.RSun,', LStar:', LStar/U.LSun, ', TStar:', TStar)

#---------------
#GRID Definition
#---------------
#Cubic grid, each edge ranges [-500, 500] AU.

sizex = sizey = sizez = 500 * U.AU
Nx = Ny = Nz = 150 #Number of divisions for each axis
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz])
NPoints = GRID.NPoints #Number of nodes in the grid

#--------
#DENSITY
#--------
Rho0 = Res.Rho0(MRate, Rd, MStar)
Arho = 24.1 #Disc-envelope density factor
Renv = 500 * U.AU #Envelope radius
Cavity = 40 * np.pi/180 #Cavity opening angle 
density = Model.density_Ulrich(RStar, Rd, Rho0, Arho, GRID, discFlag = True, envFlag = True, 
                               renv_max = Renv, ang_cavity = Cavity)

#-----------
#TEMPERATURE
#-----------
T10Env = 375. #Envelope temperature at 10 AU
BT = 5. #Adjustable factor for disc temperature. Extra, or less, disc heating.
temperature = Model.temperature(TStar, Rd, T10Env, RStar, MStar, MRate, BT, None, density, GRID,
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
xc, yc, zc = [-250*U.AU, 350*U.AU, 300*U.AU]
CENTER = [xc, yc, zc] #Center of the region in the global grid
v_sys = 3320. #m/s
newProperties = Model.ChangeGeometry(GRID, center = CENTER, vsys = v_sys,  vel = vel,
	      	 	             rot_dict = { 'angles': [np.pi/4, 1.87*np.pi], 'axis': ['x','z'] })

#At the minute, the Model library only modifies the XYZ lists. 
 #It is enough information for LIME
GRID.XYZ = newProperties.newXYZ

#The velocity should inherit the new velocity distribution 
 #(as we rotated the system and added a systemic velocity)
vel.x, vel.y, vel.z = newProperties.newVEL 

#-----------------------------
#WRITING DATA with LIME format
#-----------------------------
tag = '_Main' #A tag to identify the final files from others
Model.DataTab_LIME(density.total, temperature.total, vel, abundance, gtdratio, GRID, is_submodel = True, tag = tag)

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
r = GRID.rRTP[0] / U.AU #GRID.rRTP hosts [r, R, Theta, Phi] --> Polar GRID
Plot_model.scatter3D(GRID, density.total, weight, NRand = 4000, colordim = r, axisunit = U.AU, palette = 'jet', 
                     colorscale = 'log', colorlabel = r'${\rm log}_{10}(r [au])$', output = 'totalPoints%s.png'%tag, show = True)
