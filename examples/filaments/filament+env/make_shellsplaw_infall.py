"""
Powerlaw Shells Example
=======================

This is an example of spherical shells with non-constant power-laws for density and temperature.
The velocity model is from `Murray+2017`_
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
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
import time
import os 

t0 = time.time()

#------------------
#General Parameters
#------------------
MStar = 1*u.MSun 
r_max = 10000*u.au
r_min = r_max / 256 #Minimum distance to the centre (!= 0 to avoid indeterminations)
r_stellar = 2000*u.au #Star influence radius

rho0 = 5e8*1e6 #[part/m3] Number density at r_min
r_rho = [r_min, r_stellar, r_max] #r frontiers for density
q_rho = [-1.5, -1.7] #Powerlaws for density

T0 = 1000. #[K] Temperature at r_min
r_T = [r_min, r_stellar, r_max] #r frontiers for temperature
q_T = [-0.5, -1.0] #Powerlaws for temperature

#---------------
#GRID Definition
#---------------

sizex = sizey = sizez = r_max 
Nx = Ny = Nz = 95 #Number of divisions for each axis 
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz], rt_code = 'lime')
NPoints = GRID.NPoints #Final number of nodes in the grid

#-------------------
#PHYSICAL PROPERTIES
#-------------------
density = Model.density_PowerlawShells(r_rho, q_rho, rho0, GRID, rho_min = 1.0e4)
temperature = Model.temperature_PowerlawShells(r_T, q_T, T0, GRID, T_min = 25.)

#--------
#VELOCITY
#--------
dens_dict = {'dens_H': density.total}
ff_factor = 1.
vel = Model.velocity_infall(dens_dict,ff_factor, MStar, r_stellar, GRID)

#-------------------------------
#ABUNDANCE and GAS-to-DUST RATIO
#-------------------------------
ab0 = 1e-8 #Molecule abundance
abundance = Model.abundance(ab0, NPoints)

gtd0 = 100. #Gas to dust ratio
gtdratio = Model.gastodust(gtd0, NPoints)

#---------
#PROP dict
#---------
prop = {'dens_H': density.total,
        'temp_gas': temperature.total,
        'vel_x': vel.x,
        'vel_y': vel.y,
        'vel_z': vel.z,
        'abundance_0': abundance,
        'gtdratio': gtdratio}

Model.PrintProperties(density, temperature, GRID)
print ('Ellapsed time: %.3fs' % (time.time() - t0))
print ('-------------------------------------------------\n-------------------------------------------------\n')

#---------------
#PLOTTING 2D z=0
#---------------
Plot_model.plane2D(GRID, prop['dens_H'], axisunit=u.au, plane = {'z': 0}, norm=LogNorm(vmin=1e9, vmax=1e14), 
                   colorlabel=r'$n_{H}$ $\rm{[m^{-3}]}$', cmap='cubehelix', output='densH_2D.png', show=True)

Plot_model.plane2D(GRID, prop['temp_gas'], axisunit=u.au, plane = {'z': 0}, norm=LogNorm(vmin=1e1, vmax=1e3), 
                   colorlabel=r'T [K]', cmap='hot', output='tempgas_2D.png', show=True)

Plot_model.plane2D(GRID, prop['vel_x']/1e3, axisunit=u.au, plane = {'z': 0}, 
                   colorlabel=r'$v_x$ [km/s]', cmap='gnuplot', output='velx_2D.png', show=True)

#---------------
#1D PLOT along r
#---------------
indz = GRID.XYZ[2] == 0.0 #indices where z==0
indphi = GRID.rRTP[3][indz] == np.pi/4 #indices where phi, on z=0, equals a certain value from the list

fig, ax = plt.subplots(ncols=2, figsize=(10,5))
ax[0].plot(GRID.rRTP[0][indz][indphi]/u.au, abs(prop['vel_x'][indz][indphi]), '+') #Plotting the prop along the 1D profile
ax[0].set_xlabel('au'); #ax[0].set_xscale('log')
ax[0].set_ylabel(r'$|v_x|$ [km/s]'); #ax[0].set_yscale('log')
prop2D = prop['temp_gas'][indz].reshape(GRID.Nodes[:2]) #Reshaping to Nx,Ny matrix
prop1D = np.ma.array(np.ones(len(indphi)), mask=~indphi).reshape(GRID.Nodes[:2]) #Array masked where ind != indphi
ax[1].imshow(prop2D, norm=LogNorm(vmin=1e1, vmax=1e3), cmap='gnuplot') #Illustrating where the 1D profile is coming from
ax[1].imshow(prop1D, cmap='binary_r') #Illustrating where the 1D profile is coming from
ax[1].set_title('Temperature on z=0')
fig.savefig('1Dprofile.png')
plt.show()

#------------------------
#RE-CENTERING and WRITING
#------------------------
shift = Model.ChangeGeometry(GRID, center=[15000*u.au,0,0])
GRID.XYZ = shift.newXYZ

lime = rt.Lime(GRID)
lime.submodel(prop, output='envelope.dat')
