import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors

from sf3dmodels.Model import Struct
from sf3dmodels import Plot_model
from sf3dmodels.grid import fillgrid
from sf3dmodels.arepo import UniqueCells
from sf3dmodels.utils import (constants as sfc, units as sfu)
from sf3dmodels.arepo.units import *
import sf3dmodels.rt as rt

import binary_read as rsnap
import cgs_constants as cgs

from copy import copy, deepcopy
import os
import sys
import time

t0 = time.time()

plot_hist = True #plot diagnostic 2D histograms?
random_choice = None #If integer, considers only a (random) smaller portion from the total grid points for the radiative transfer.

#******************************
#PATH TO SNAPSHOT
#******************************
snap = 'Restest_extracted_001_240'

#******************************
#AREPO 3D-grids (not projected)
#******************************
rsnap.io_flags['mc_tracer']=True
rsnap.io_flags['time_steps']=True
rsnap.io_flags['sgchem']=True
rsnap.io_flags['variable_metallicity']=False
#rsnap.io_flags['MHD']=True #This simulation is not magnetised
xHe = 0.1 #Helium abundance

#******************************
#READING SNAPSHOT
#******************************
data, header = rsnap.read_snapshot(snap)

#yn_real = data['rho']*udensity / ((1. + 4.0 * xHe + 2 * NH2_ + 1 * NHp_) * cgs.mp)
yn = data['rho'] * udensity / ((1. + 4.0 * xHe) * cgs.mp) #Number density of H-species 
energy = data['u_therm'] * data['rho'] * uenergy/ulength/ulength/ulength
yntot = (1. + xHe - data['chem'][:,0] + data['chem'][:,1]) * yn
data['tgas'] = 2.*energy / (3.*yntot*cgs.kb)

#*************************************
#DICTS TO OBJECTS AND SOME DEFINITIONS
#*************************************
nparts = header['num_particles']
ngas = nparts[0]
ntracer = nparts[2]
nsinks = nparts[5]
total_part = ngas + nsinks

pos = data['pos'][:ngas]
gasmass = data['mass'][:ngas]
sinkmass = gasmass[ngas:ngas+nsinks]

print('Mass in gas', np.sum(gasmass))
print('Mass in sinks', np.sum(sinkmass))

#******************************
#2D HISTROGRAMS
#******************************
if plot_hist:
    from plot_histograms import plot_all
    rr = np.linalg.norm(pos, axis=1)
    iref = np.where(rr!=None)
    vol = gasmass/data['rho']
    dcell = vol**0.333
    plot_all(data, yn, gasmass, dcell, rr, iref)

#***********************************************************
#MEAN DENSITY - weighted by mass
#***********************************************************
mean_by_mass = lambda prop, mass: np.sum(prop * mass) / mass.sum()  
mean_dens = mean_by_mass(data['rho'], data['mass'][0:ngas])
mean_dens *= udensity/((1. + 4.0 * xHe) * cgs.mp)

print ('mean number density (snapshot) in cm^-3:', mean_dens)

#*******************************************
#PRE-PROCESSING FOR LIME
#Note: LIME's input properties must be in SI
#*******************************************

#Microturbulence
def get_sound_speed(T):
    return np.sqrt(5/3. * sfc.kb * T / sfc.mH)

#Function to prepare the physical data:
def get_prop_dicts(id_pos):

    dens_3d = data['rho'][id_pos] * udensity / ((1. + 4.0 * xHe) * cgs.mp) * 1e6 #kg/m3 --> Number density of H nucleons 
    dens_H2 = dens_3d*data['chem'][:,0][id_pos]
    dens_Hplus = dens_3d*data['chem'][:,1][id_pos]
    coords = ['x','y','z']
    mean_vel = [np.mean(data['vel'][:,i][id_pos]) for i in range(3)]
    vel_3d = Struct( **{coords[i]: (data['vel'][:,i][id_pos] - mean_vel[i]) * uvel / 100. for i in range(3)} )
    
    prop = {'dens_H2': dens_H2,
            'dens_Hplus': dens_Hplus,
            'dens_H': dens_3d - 2*dens_H2 - dens_Hplus,
            'temp_gas': data['tgas'][id_pos],
            'temp_dust': data['tdust'][id_pos],
            'doppler': get_sound_speed(data['tgas'][id_pos]),
            'abundance': data['chem'][:,2][id_pos],
            'vel_x': vel_3d.x,
            'vel_y': vel_3d.y,
            'vel_z': vel_3d.z,
            'gtdratio': 100 * np.ones(GRID.NPoints)
            }
    
    prop_fill = {'dens_H2': 0.1,
                 'dens_Hplus': 0.1,
                 'temp_gas': 2.725,
                 'temp_dust': 2.725,
                 'abundance': 1e-8,
                 'vel_x': 1e6,
                 'vel_y': 1e6,
                 'vel_z': 1e6,
                 'gtdratio': 100,
                 }
    
    return prop, prop_fill

#******************
#MERGING TWIN CELLS
#******************
#Merging twin cell masses into a single cell 
data_orig = deepcopy(data) 
unique_cells = UniqueCells(data, header)
id_pos = unique_cells.mergemass() #This function returns the indices from non-twin cells 

if random_choice is not None:
    id_pos_0 = id_pos
    id_pos = np.random.choice(id_pos, size=random_choice, replace=False)

#************
#GRID OBJECT
#************
pos_max = np.max(data['pos'], axis=0)
pos_min = np.min(data['pos'], axis=0)
pos_mean = 0.5 * (pos_max + pos_min) * ulength / 100. #Shifting mean position to (0,0,0)

GRID = Struct( XYZ = np.array([ data['pos'][:,i][id_pos] * ulength / 100 - pos_mean[i] for i in range(3)]) ) #Making GRID.XYZ object
GRID.rRTP = [np.linalg.norm(GRID.XYZ, axis=0), None, None, None] #Only interested in the coordinate r (for fillgrid() later).
GRID.NPoints = len(id_pos)
    
#********************************************
#FILLING OUTER EMPTY SPACES WITH DUMMY CELLS
#********************************************
prop, prop_fill = get_prop_dicts(id_pos)
fill = fillgrid.Random(GRID)
fill_rand = fill.by_mass(data['mass'][id_pos], prop,
                         prop_fill = prop_fill,
                         mass_fraction = 0.9, 
                         n_dummy = GRID.NPoints/10, 
                         r_steps = 1000)

#*******************************
#WRITING FORMATTED DATA FOR LIME
#*******************************
lime = rt.Lime(GRID)
lime.submodel(prop, output='datatab.dat', lime_npoints = True, lime_header = True)
print('Output columns', lime.columns)

#******-
#TIMING
#******-
print ('Ellapsed time: %.3fs' % (time.time() - t0))
print ('-------------------------------------------------\n-------------------------------------------------\n')

#*******************************************
#3D PLOTTING
# The scatter points are randomly
# tracing the hydrogen (H) density.
#   The colormap is associated to the
#   molecular hydrogen (H2) density in cgs
#*******************************************
weight = 10*mean_dens

canvas3d = Plot_model.Canvas3d(ax_kw={'azim': -50, 'elev': 50})
ax = canvas3d.ax
sp = canvas3d.scatter_random(GRID, prop['dens_H'], weight, prop_color = prop['dens_H2']/1e6, GRID_unit=sfu.pc, power=0.5, NRand=20000, prop_min=None, #function arguments
                             marker = 'o', s = 3, cmap = 'nipy_spectral_r', edgecolors = 'none', vmin = 1.0, vmax = None, norm = colors.LogNorm()) #Scatter kwargs

cbar = plt.colorbar(sp)
cbar.ax.set_ylabel(r'$\rm n_{H_2} [\rm cm^{-3}]$')

ax.set_xlabel('X (pc)')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.savefig('3Dpoints_snap.png')#, dpi=200)
plt.show()
