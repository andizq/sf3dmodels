#******************
#sf3dmodels modules
#******************
from sf3dmodels import Model, Plot_model
import sf3dmodels.utils.units as u
import sf3dmodels.utils.constants as ct
import sf3dmodels.rt as rt
from sf3dmodels.model import disc
from sf3dmodels.grid import Grid, fillgrid
from sf3dmodels.arepo import UniqueCells
#******************
#External libraries
#******************
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import AutoMinorLocator
from astropy.io import fits
import time
import sys
#********************************
#Gridding and writing for Polaris
#********************************
from scipy.spatial import Delaunay
from collections import defaultdict
import struct
import bisect

t0 = time.time()
nfrac = 300000*1

file = 'snap_00108.ascii'
data = np.loadtxt(file) 
#Quantities in code units: i.e length referred to semi-major axis between stars a0, mass referred to mass of stars, time referred to orbital period between binaries

a0 = 5*u.au
Ms = 2*u.MSun
Rs = 2.5*u.RSun
Teff = 9500.

mu = ct.G * 2*Ms #Reduced mass
P = 2*np.pi*a0**1.5/mu**0.5 #Orbital period

star1_pos = np.array([3.1349969E-01,   3.9121187E-01,  -8.4040770E-03]) * a0 
star2_pos = np.array([-3.1147322E-01,  -3.8484186E-01,   1.1820099E-02]) * a0

data_gas = data[data[:,12]==1]
ngas = len(data_gas)

if nfrac: 
    ids_rand = np.random.choice(np.arange(0,ngas), size=nfrac)
    data_gas = data_gas[ids_rand]

data_gas[:,0:3] *= a0
data_r = np.linalg.norm(data_gas[:,0:3], axis=1)
data_gas = data_gas[data_r <= 100*u.au]
data_gas[:,5] *= 0.01*Ms/a0**3/ct.mH
data_gas[:,5] = np.where(data_gas[:,5]<1.0, 1.0, data_gas[:,5])
data_gas[:,6:9] *= a0/P


grid = Model.Struct(XYZ = data_gas[:,0:3].T, NPoints=len(data_gas))
prop = {'dens_H2': data_gas[:,5],
        'vel_x': data_gas[:,6],
        'vel_y': data_gas[:,7],
        'vel_z': data_gas[:,8]}

unique_cells = UniqueCells(prop, None)
id_pos = unique_cells.getuniques(grid).astype(np.int)

grid.XYZ = grid.XYZ.T[id_pos].T
grid.NPoints = len(id_pos)
for key in prop: prop[key] = prop[key][id_pos]

fill = fillgrid.Random(grid, smart=False)
fill_rand = fill.spherical(prop,
                           prop_fill = {'dens_H2': 1.0, 'temp_dust': 3.0, 'temp_gas': 3.0},
                           r_min = 1*u.au, r_max = 100*u.au,
                           n_dummy = grid.NPoints/10.)

lime = rt.Lime(grid)
lime.submodel(prop, output='datatab.dat', folder='./', lime_npoints = True, lime_header = True)
print('Output columns', lime.columns)

print ('Ellapsed time:', time.time()-t0)

#******************
#3D plotting
#******************
lims = np.array([-100,100])
weight = 1. 
ax_kw = {'projection': '3d', 'xlim': lims, 'ylim': lims, 'zlim': lims, 'azim': 30, 'elev': 13}
canvas3d = Plot_model.Canvas3d(ax_kw=ax_kw)
ax = canvas3d.ax
sp = canvas3d.scatter_random(grid, prop['dens_H2']/1e6, weight, GRID_unit=u.au, power=0, NRand=10000, prop_min=1.0, #function arguments
                             marker = '+', cmap = 'jet', s = 3, edgecolors = 'none', vmin=1e7, norm = colors.LogNorm()) #Scatter kwargs

cbar = plt.colorbar(sp)
cbar.ax.set_ylabel(r'H$_2$ density [cm$^{-3}$]')
ax.set_xlabel('au')
plt.savefig('img_H2_density_sph.png')
plt.show()


