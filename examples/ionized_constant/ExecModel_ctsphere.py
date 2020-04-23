"""
Constant Sphere Example
=======================

This is an exmaple of a uniform-density (constant-density) spherical HII region
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
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import time

#------------------
#General Parameters
#------------------
r_max = 2530 * u.au #H II sphere size
dens_e = 1.4e5 * 1e6 #Electron number density, from cgs to SI
t_e = 1.e4 #K

#---------------
#GRID Definition
#---------------
sizex = sizey = sizez = 2600 * u.au 
Nx = Ny = Nz = 63 #Number of divisions for each axis
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz], rt_code='radmc3d')
NPoints = GRID.NPoints #Final number of nodes in the grid

#-------------------
#PHYSICAL PROPERTIES
#-------------------
density = Model.density_Constant(r_max, GRID, envDens = dens_e)
temperature = Model.temperature_Constant(density, GRID, envTemp = t_e, backTemp = 2.725)

Model.PrintProperties(density, temperature, GRID) #Printing resultant properties (mass, mean temperature, etc)

#----------------------
#WRITING RADMC-3D FILES
#----------------------
prop = {'dens_e': density.total,
        'dens_ion': density.total,
        'temp_gas': temperature.total}

Rad = rt.Radmc3dDefaults(GRID)
Rad.freefree(prop)

#------------------------------------
#3D PLOTTING (weighting with density)
#------------------------------------
tag = 'ctsphere_HII'
weight = dens_e

"""
Plot_model.scatter3D(GRID, density.total, weight, power=0, NRand = 4000, colordim = density.total / 1e6 / 1e5, axisunit = u.au, cmap = 'winter', 
                     marker = 'o', colorlabel = r'$n_{\rm e}$ [cm$^{-3}$] x $10^5$', output = '3Ddens_%s.png'%tag, show = True)

Plot_model.scatter3D(GRID, density.total, weight, power=0, NRand = 4000, colordim = temperature.total, axisunit = u.au, cmap = 'winter', 
                     marker = 'o', colorlabel = r'$T_{\rm e}$ [Kelvin]', output = '3Dtemp_%s.png'%tag, show = True)
"""

fig = plt.figure(figsize=(6.4,4.8))
#ax = plt.axes(projection='3d')
lims = np.array([-sizex,sizex]) / u.au 
canvas3d = Plot_model.Canvas3d(fig=fig, ax_kw={'xlim': lims, 'ylim': lims, 'zlim': lims, 'azim': -50, 'elev': 30})
ax = canvas3d.ax #generated with fig.add_axes from matplotlib. All the matplotlib functions are therefore available on ax.
sp = canvas3d.scatter_random(GRID, density.total/1e6/1e5, weight/1e6/1e5, GRID_unit=u.au, power=0, NRand=4000, prop_min=1.0, #function arguments
                             marker = 'o', cmap = 'winter', s = 3, edgecolors = 'none', vmin = None, vmax = None, norm = None) #Scatter kwargs
#ax.scatter(3000,0,0, c=[dens_e*1.2/1e6], norm=colors.Normalize(density.total.min()/1e6, vmax=density.total.max()/1e6), cmap=sp.cmap)

x,y,z = GRID.XYZ
ind_rand = np.random.choice(np.arange(len(x)), size=10)
vel = np.zeros((3,10)) + 1000*np.sin([np.linspace(0,2*np.pi,10) for i in range(3)])#np.random.random(size=(3,10))
vp = canvas3d.vectors(x[ind_rand], y[ind_rand], z[ind_rand], vel[0], vel[1], vel[2], GRID_unit=u.au, length=1.0, linewidths=2, color='k')

cbar = plt.colorbar(sp)
cbar.ax.set_ylabel(r'$n_{e^-}$ ($\times$ 10$^5$) cm$^{-3}$')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.savefig('3Dtemp_%s.png'%tag, dpi = 500)#, bbox_inches='tight') 
# when a single density is used, bbox_inches='tight' associates a wrong colour to the scatter points (uses the min val of the colorbar instead of the mid value).
plt.show()



