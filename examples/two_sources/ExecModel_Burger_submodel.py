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
MStar = 0.86 * u.MSun 
MRate = 5.e-6 * u.MSun_yr 
RStar = u.RSun * ( MStar/u.MSun )**0.8 
LStar = u.LSun * ( MStar/u.MSun )**4 
TStar = u.TSun * ( (LStar/u.LSun) / (RStar/u.RSun)**2 )**0.25 
Rd = 264. * u.au

print ('RStar:', RStar/u.RSun, ', LStar:', LStar/u.LSun, ', TStar:', TStar)

#---------------
#GRID Definition
#---------------
#Cubic grid, each edge ranges [-500, 500] au.

sizex = sizey = sizez = 500 * u.au
Nx = Ny = Nz = 100 #Number of divisions for each axis
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz])
NPoints = GRID.NPoints #Number of nodes in the grid

#-------------
#DENSITY
#-------------

#--------
#ENVELOPE
#--------
Rho0 = Res.Rho0(MRate, Rd, MStar)
Arho = None
Renv = 2.5 * Rd
densEnv = Model.density_Env_Disc(RStar, Rd, Rho0, Arho, GRID, discFlag = False, envFlag = True, 
                                 renv_max = Renv)
#-------
#DISC
#-------
H0sf = 0.03 #Disc scale height factor (H0 = 0.03 * RStar)
Arho = 5.25
Rdisc = 1.5 * Rd
densDisc = Model.density_Hamburgers(RStar, H0sf, Rd, Rho0, Arho, GRID, discFlag = True, 
                                    rdisc_max = Rdisc)
#---------------------
#The COMPOSITE DENSITY
#---------------------
density = Model.Struct( **{ 'total': densEnv.total + densDisc.total,
                            'disc': densDisc.total, 
                            'env': densEnv.total,
                            'H': densDisc.H,
                            'discFlag': True,
                            'envFlag': True,
                            'Rt': densDisc.Rt,
                            'r_disc': densDisc.r_disc, 
                            'r_env': densEnv.r_env,
                            'streamline': densEnv.streamline} )

#-----------
#TEMPERATURE
#-----------
T10Env = 250. #Envelope temperature at 10 au
Tmin = 10. #Minimum possible temperature. Every node with T<Tmin will inherit Tmin. 
BT = 60. #Adjustable factor for disc temperature. Extra, or less, disc heating.
temperature = Model.temperature_Hamburgers(TStar, RStar, MStar, MRate, Rd, T10Env, BT, 
                                           density, GRID, Tmin_disc = Tmin, inverted = False)

#--------
#VELOCITY
#--------
vel = Model.velocity(RStar, MStar, Rd, density, GRID)

#-------------------------------
#ABUNDANCE and GAS-to-DUST RATIO
#-------------------------------
ab0 = 5e-8 #CH3CN abundance vs H2
abundance = Model.abundance(ab0, NPoints)

gtd0 = 100. #Gas to dust ratio (H2 vs Dust)
gtdratio = Model.gastodust(gtd0, NPoints)

#-------------------------
#ROTATION, VSYS, CENTERING
#-------------------------
xc, yc, zc = [350*u.au, -150*u.au, -200*u.au]
CENTER = [xc, yc, zc] #Center of the region in the global grid
v_sys = -2000. #Systemic velocity (vz) of the region (in m/s)
newProperties = Model.ChangeGeometry(GRID, center = CENTER, vsys = v_sys,  vel = vel,
	      	 	             rot_dict = { 'angles': [np.pi/2, np.pi/3], 'axis': ['x','z'] })

#At the minute, the Model library only modifies the XYZ lists. 
 #This is enough information for LIME
GRID.XYZ = newProperties.newXYZ

#The velocity should inherit the new velocity distribution 
 #(as we rotated the system and added a systemic velocity)
vel.x, vel.y, vel.z = newProperties.newVEL 

#-----------------------------
#WRITING DATA for LIME
#-----------------------------
tag = 'Burger.dat'
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

#----------------------------------------
#3D PLOTTING (weighting with temperature)
#----------------------------------------
tag = 'Burger'
weight = 10*T10Env

vmin, vmax = np.array([5e11, 5e15]) / 1e6
norm = colors.LogNorm(vmin=vmin, vmax=vmax)

Plot_model.scatter3D(GRID, temperature.total, weight, NRand = 4000, colordim = density.total / 1e6 , axisunit = u.au, cmap = 'hot', 
                     norm = norm, colorlabel = r'$\rho$ $[cm^{-3}]$', output = '3Dpoints%s.png'%tag, show = True)


