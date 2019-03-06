"""
Basic docstring explaining example
"""
from __future__ import print_function
#********************
#sf3dmodels libraries
#********************
from sf3dmodels.outflow import OutflowModel   #Model functions
import sf3dmodels.utils.units as u            #Units
import sf3dmodels.rt as rt                    #Writing functions for radiative transfer
#********************
#Extra libraries
#********************
import numpy as np
import time

t0 = time.time()
#--------------------

dx_grid = 7*u.au
#*********************
#OUTFLOW 1
#*********************
tag = '_outflow1'

#---------------------
#GEOMETRIC PARAMETERS
#---------------------
pos_c = np.array([-200*u.au, 50*u.au, 0])
axis = np.array([0,1,0]) 
z_min = 5*u.au
z_max = 350*u.au

#---------------------
#PHYSICAL PROPERTIES
#---------------------
w0 = 12*u.au
eps = 0.45
w = [w0, eps]

T0 = 10000.
qT = -0.6
temp = [T0, qT]

dens0 = 5.e13
qn = -0.9
dens = [dens0, qn]

ionfrac = [1,-0.5]

abund = [1e-4, 0]
gtd = 100

v0 = 200 * 1e3 #km/s

#---------------------
#COMPUTING MODEL
#---------------------
Outf1 = OutflowModel(pos_c, axis, z_min, z_max, dx_grid) #Initializing Class with grid parameters
Outf1.reynolds86(w, dens, ionfrac, temp, v0, abund, gtd) #Invoking the outflow model from Reynolds et al. 1986

#---------------------
#WRITING FOR RADMC3D
#---------------------
#Using the Radmc3d class
prop1 = {'dens_e': Outf1.density,
         'dens_ion': Outf1.density,
         'temp_gas': Outf1.temperature}
radmc1 = rt.Radmc3d(Outf1.GRID)
radmc1.submodel(prop1, output='datatab'+tag+'.dat')
print('Output columns', radmc1.columns)

#*********************
#OUTFLOW 2
#*********************
tag = '_outflow2'

#---------------------
#GEOMETRY
#---------------------
pos_c = np.array([200*u.au, -50*u.au, 0])
axis = np.array([1,1,-1]) 
z_min = 5*u.au
z_max = 300*u.au

#---------------------
#PHYSICAL PROPERTIES
#---------------------
w0 = 10*u.au
eps = 3/4.
w = [w0, eps]

T0 = 10000.
qT = -1.
temp = [T0, qT]

dens0 = 7.e13
qn = -3/2.
dens = [dens0, qn]

ionfrac = [1,0]

abund = [1e-4, 0]
gtd = 100

v0 = 150 * 1e3 #km/s

#---------------------
#COMPUTING MODEL
#---------------------
Outf2 = OutflowModel(pos_c, axis, z_min, z_max, dx_grid) #Initializing Class with grid parameters
Outf2.reynolds86(w, dens, ionfrac, temp, v0, abund, gtd) #Invoking the outflow model from Reynolds et al. 1986

#---------------------
#WRITING FOR RADMC3D
#---------------------
#Using the Radmc3d class
prop2 = {'dens_e': Outf2.density,
         'dens_ion': Outf2.density,
         'temp_gas': Outf2.temperature}
radmc2 = rt.Radmc3d(Outf2.GRID)
radmc2.submodel(prop2, output='datatab'+tag+'.dat')
print('Output columns', radmc2.columns)

#-------
#TIMING
#-------
print ('Ellapsed time: %.3fs' % (time.time() - t0))
print ('-------------------------------------------------\n-------------------------------------------------\n')
