"""
Powerlaw Shells Example
=======================

This is an example of spherical shells with non-constant power-laws for density and temperature.
"""
#------------------
#Import the package
#------------------
from sf3dmodels import Model, Plot_model
import sf3dmodels.utils.units as u
import sf3dmodels.utils.constants as ct
import sf3dmodels.rt as rt                   
#-----------------
#Extra libraries
#-----------------
from matplotlib.colors import LogNorm
import numpy as np
import time
import os 

t0 = time.time()

#------------------
#General Parameters
#------------------
r_max = 500*u.au
r_min = r_max / 256 #Minimum distance to the centre (!= 0 to avoid indeterminations)

rho0 = 1e8*1e6 #[part/m3] Number density at r_min
r_rho = [r_min, 200*u.au, 300*u.au, r_max] #r frontiers for density
q_rho = [-2.0, -1.5, -2.5] #Powerlaws for density

T0 = 1000. #[K] Temperature at r_min
r_T = [r_min, 250*u.au, r_max] #r frontiers for temperature
q_T = [-0.5, -1.0] #Powerlaws for temperature

#---------------
#GRID Definition
#---------------
sizex = sizey = sizez = r_max 
Nx = Ny = Nz = 95 #Number of divisions for each axis 
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz], rt_code = 'radmc3d')
NPoints = GRID.NPoints #Final number of nodes in the grid

#-------------------
#PHYSICAL PROPERTIES
#-------------------
density = Model.density_PowerlawShells(r_rho, q_rho, rho0, GRID, rho_min = 1.0e4)
temperature = Model.temperature_PowerlawShells(r_T, q_T, T0, GRID, T_min = 25.)

gtd0 = 100. #Gas to dust ratio
gtdratio = Model.gastodust(gtd0, NPoints)

#--------------------
#PRINTING and WRITING
#--------------------
prop = {'dens_H2': density.total,
        'temp_gas': temperature.total,
        'dens_dust': 2*ct.mH * density.total * 1/gtdratio, #dust mass density
        'gtdratio': gtdratio}

radmc = rt.Radmc3d(GRID, nphot=1000000)
wavelength_intervals = [1e-1,5e2,1e4]
wavelength_divisions = [20,20]
radmc.write_radmc3d_control(incl_dust=1, setthreads=8, incl_freefree=0)
radmc.write_amr_grid()
radmc.write_dust_density(prop['dens_dust']) #Needs mass density (not number density)
radmc.write_stars(nstars=1, pos=[[0,0,0]], rstars = [u.RSun], mstars = [u.MSun], flux = [[-u.TSun]], #flux --> if negative, radmc assumes the input number as the blackbody temperature of the star 
                  lam = wavelength_intervals, nxx = wavelength_divisions) 
radmc.write_wavelength_micron(lam = wavelength_intervals, nxx = wavelength_divisions) #lam --> wavelengths in microns, nxx --> number of divisions in between wavelengths

Model.PrintProperties(density, temperature, GRID)                

#---------------
#PLOTTING 2D z=0
#---------------
Plot_model.plane2D(GRID, prop['dens_H2'], axisunit=u.au, plane = {'z': 0}, norm=LogNorm(vmin=1e9, vmax=1e13), output='densH2_2D.png')
Plot_model.plane2D(GRID, prop['temp_gas'], axisunit=u.au, plane = {'z': 0}, norm=LogNorm(vmin=1e1, vmax=1e3), output='tempgas_2D.png')

#-------
#TIMING
#-------
print ('Ellapsed time: %.3fs' % (time.time() - t0))
print ('-------------------------------------------------\n-------------------------------------------------\n')
