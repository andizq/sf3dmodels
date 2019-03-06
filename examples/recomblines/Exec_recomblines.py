"""
Basic docstring explaining example
"""
from __future__ import print_function
#********************
#sf3dmodels libraries
#********************
from sf3dmodels import Model
import sf3dmodels.utils.units as u
from sf3dmodels.grid import Overlap
import sf3dmodels.rt as rt
#********************
#Extra libraries
#********************
import numpy as np
import time

t0 = time.time()
#--------------------

#**************
#GRID AND MODEL
#**************
r_max = 0.1 * u.pc  #H II sphere size
r_min = 0.02 * u.pc
rho0 = 28 * 1e6 * u.pc**2
q = -2

sizex = sizey = sizez = 1.2*r_max
Nx = Ny = Nz = 100 #Number of divisions along each axis
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz], rt_code = 'radmc3d')
NPoints = GRID.NPoints #Number of nodes in the grid

density = Model.density_Powerlaw2(r_max,r_min,rho0,q,GRID,rho_min=1.)
temperature = Model.temperature_Constant(density, GRID, envTemp=10000., backTemp=2.725)
vel = Model.velocity_random(15000,GRID.NPoints)

#********************
#WRITING FOR RADMC-3D
#********************
prop = {'vel_x' : vel.x, 'vel_y' : vel.y, 'vel_z' : vel.z, 
        'dens_e' : density.total, 'dens_ion' : density.total, 
        'temp_gas' : temperature.total}
A = rt.Radmc3dDefaults(GRID)
A.recomblines(prop, [110,109])
