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
import sf3dmodels.Plot_model as Pm
#********************
#Extra libraries
#********************
import numpy as np
import time

t0 = time.time()
dx_grid = 1*u.au

#*********************
#OUTFLOW
#*********************

#---------------------
#GEOMETRIC PARAMETERS
#---------------------
pos_c = np.array([0*u.au, 0*u.au, 0])
axis = np.array([0,0,1]) 
z_min = 1*u.au
z_max = 100*u.au

#---------------------
#PHYSICAL PROPERTIES
#---------------------
w0 = 2*u.au
eps = 0.45
w = [w0, eps]

T0 = 10000.
qT = -0.5
temp = [T0, qT]

dens0 = 1.e14
qn = -0.5
dens = [dens0, qn]

ionfrac = [1,-0.5]

abund = [1e-4, 0]
gtd = 100

v0 = 20 * 1e3 #km/s

#---------------------
#COMPUTING MODEL
#---------------------
Outf = OutflowModel(pos_c, axis, z_min, z_max, dx_grid) #Initializing Class with grid parameters
Outf.reynolds86(w, dens, ionfrac, temp, v0, abund, gtd) #Invoking the outflow model from Reynolds et al. 1986

#---------------------
#WRITING FOR LIME
#---------------------
prop = {'dens_Hplus': Outf.density_ion,
        'dens_H2': np.ones(Outf.GRID.NPoints), #THERE IS NO H2 BUT THIS IS NEEDED FOR OVERLAPING PURPOSES
        'temp_gas': Outf.temperature,
        'vel_x': Outf.vel.x,
        'vel_y': Outf.vel.y,
        'vel_z': Outf.vel.z,
        'abundance': Outf.abundance,
        'gtdratio': Outf.gtdratio
        }

lime = rt.Lime(Outf.GRID)
lime.submodel(prop, output='outflow.dat')
print('Output columns', lime.columns)

#-----------------
#Plot for DENSITY
#-----------------
Pm.scatter3D(Outf.GRID, prop['dens_Hplus']/1e6, weight=1.0, NRand = 4000, axisunit = u.au, colorscale = 'log', cmap = 'cool',
             colorlabel = r'${\rm log}_{10}(n [cm^{-3}])$', output = 'global_grid_dens.png', vmin = 5)

#-------
#TIMING
#-------
print ('Ellapsed time: %.3fs' % (time.time() - t0))
print ('-------------------------------------------------\n-------------------------------------------------\n')
