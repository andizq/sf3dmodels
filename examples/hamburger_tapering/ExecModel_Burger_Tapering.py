"""
Basic docstring explaining example
"""
from __future__ import print_function
#-----------------
#Package libraries
#-----------------
from sf3dmodels import *
#-----------------
#Extra libraries
#-----------------
from matplotlib import colors
import numpy as np
import os
import time
from argparse import ArgumentParser

t0 = time.time()

parser = ArgumentParser(prog='Burgers', description='Self-obscuration')
parser.add_argument('-A', '--Arho', help='Disc density factor')
args = parser.parse_args()

#------------------
#General Parameters
#------------------
MStar = 0.3 * U.MSun #0.3 Msun from paper by Li+2017
MRate = 1.e-5 * U.MSun_yr 
RStar = 10 * U.RSun * ( MStar/U.MSun )**0.8 #Multiply by 10 Main Sequence relationship to account for bloated protostar 
LStar = 9 * U.LSun #U.LSun * ( MStar/U.MSun )**4 # 9Lsun from Li+2017
TStar = U.TSun * ( (LStar/U.LSun) / (RStar/U.RSun)**2 )**0.25 
Rd = 68. * U.AU

print ('RStar:', RStar/U.RSun, ', LStar:', LStar/U.LSun, ', TStar:', TStar)

#---------------
#GRID Definition
#---------------
#Cubic grid, each edge ranges [-100, 100] AU.

sizex = sizey = sizez = 100 * U.AU
Nx = Ny = Nz = 200 #Number of divisions for each axis
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz])
NPoints = GRID.NPoints #Number of nodes in the grid

#-------------
#DENSITY
#-------------

#-------
#DISC
#-------
if args.Arho: Arho = float(args.Arho)
else: Arho = 150.0  # (Mdisk=0.106*1.4=Msun=0.149Msun)

H0sf = 0.04 #Disc scale-height factor (H0 = H0sf * RStar)
Rdisc = 1.0 * Rd
Rho0 = Res.Rho0(MRate, Rd, MStar) #Normalization density
Rtap = 36 * U.AU #Radius where the tapering begins
density = Model.density_Hamburgers(RStar, H0sf, Rd, Rho0, Arho, GRID, 
                                   discFlag = True, Rt = Rtap, rdisc_max = Rdisc)

#-----------
#TEMPERATURE
#-----------
T10Env = None #There is no envelope
Tmin = 10. #Minimum possible temperature. Every node with T<Tmin will inherit Tmin. 
BT = 2.7 #Adjustable factor for disc temperature. Extra, or less, disc heating.
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

#-----------------------------
#WRITING DATA with LIME format
#-----------------------------
Model.DataTab_LIME(density.total, temperature.total, vel, abundance, gtdratio, GRID) 

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
tag = 'Burger_Tapering'
dens_plot = density.total / 1e6

weight = 5000. #Kelvin
norm = colors.LogNorm()
Plot_model.scatter3D(GRID, temperature.total, weight, NRand = 4000, colordim = density.total / 1e6 , axisunit = U.AU, cmap = 'ocean_r', 
                     norm = norm, power = 0.8, colorlabel = r'$\rho [cm^{-3}]$', output = '3Dpoints%s.png'%tag, show = False)

#----------------------------------------
#2D PLOTTING (Density and Temperature)
#----------------------------------------

vmin, vmax = np.array([9e15, 5e19]) / 1e6
norm = colors.LogNorm(vmin=vmin, vmax=vmax)

Plot_model.plane2D(GRID, dens_plot, axisunit = U.AU, cmap = 'ocean_r', plane = {'z': 0*U.AU},
                   norm = norm, colorlabel = r'$[\rm cm^{-3}]$', output = 'DensMidplane_%s.png'%tag, show = False)

vmin, vmax = np.array([5e14, 1e18]) / 1e6
norm = colors.LogNorm(vmin=vmin, vmax=vmax)

Plot_model.plane2D(GRID, dens_plot, axisunit = U.AU, cmap = 'ocean_r', plane = {'y': 0*U.AU},
                   norm = norm, colorlabel = r'$[\rm cm^{-3}]$', output = 'DensVertical_%s.png'%tag, show = False)

vmin, vmax = np.array([27, 1e3])
norm = colors.LogNorm(vmin=vmin, vmax=vmax)

Plot_model.plane2D(GRID, temperature.total, axisunit = U.AU, cmap = 'ocean_r', plane = {'z': 0*U.AU},
                   norm = norm, colorlabel = r'[Kelvin]', output = 'TempMidplane_%s.png'%tag, show = False)

vmin, vmax = np.array([1e1, 1e3])
norm = colors.LogNorm(vmin=vmin, vmax=vmax)

Plot_model.plane2D(GRID, temperature.total, axisunit = U.AU, cmap = 'ocean_r', plane = {'y': 0*U.AU},
                   norm = norm, colorlabel = r'[Kelvin]', output = 'TempVertical_%s.png'%tag, show = False)


vmin, vmax = np.array([1e10, 1e14])
norm = colors.LogNorm(vmin=vmin, vmax=vmax)

Plot_model.plane2D(GRID, temperature.total * dens_plot, axisunit = U.AU, cmap = 'ocean_r', plane = {'y': 0*U.AU},
                   norm = norm, colorlabel = r'[$\rho$ T]', output = 'Emissivity_%s.png'%tag, show = False)
