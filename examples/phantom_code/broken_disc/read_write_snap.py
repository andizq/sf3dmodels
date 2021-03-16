"""
Author: Andres F. Izquierdo
This is a customised version of the converter script AREPOtoPOLARIS
provided by Reissl,S and Tress,R on the Polaris website
"""
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
from scipy.spatial import Delaunay, Voronoi
from collections import defaultdict
import struct
import bisect

t0 = time.time()
nfrac = 300000

file = './snap_00108.ascii'
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

data_gas = data[data[:,12]==1] #Gas particles id=1, others are sinks or swallowed particles
ngas = len(data_gas)

if nfrac: 
    ids_rand = np.random.choice(np.arange(0,ngas), size=nfrac, replace=False)
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

unique_cells = UniqueCells(prop,None)
id_pos = unique_cells.getuniques(grid).astype(np.int) #Removing twin particles

grid.XYZ = grid.XYZ.T[id_pos].T
grid.NPoints = len(id_pos)
for key in prop: prop[key] = prop[key][id_pos]

fill = fillgrid.Random(grid)
fill_rand = fill.spherical(prop,
                           prop_fill = {'dens_H2': 1.0},
                           r_min = 1*u.au, r_max = 100*u.au,
                           n_dummy = grid.NPoints/10.) #Smart distribution of dummies

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

#*****************************************************
#Writing file to compute dust temperature with POLARIS
#*****************************************************
output_file = 'voronoi_grid.dat'
grid_id = 50 #polaris grid ID (50 = voronoi)
data_ids = [0]

def is_in_hull(a, x):
    i = bisect.bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return -1

    return 1

def find_neighbours(tri):
    neighbours = defaultdict(set)

    tri_len=len(tri.simplices)
    counter = 0

    for simplex in tri.simplices:
        if counter % 100000 == 0:
            sys.stdout.write('Determining neighbours: ' + str(100.0 * counter / tri_len) + '%    \r')
            sys.stdout.flush()

        counter+=1

        for idx in simplex:
            other = set(simplex)
            other.remove(idx)
            neighbours[idx] = neighbours[idx].union(other)

    return neighbours

def create_convex_hull_list(hull):
    neighbours = defaultdict(set)

    for simplex in hull:
        for idx in simplex:
            other = set(simplex)
            other.remove(idx)
            neighbours[0] = neighbours[0].union(other)

    return sorted(list(neighbours[0]))

def tetrahedron_volume(a, b, c, d):
    return np.abs(np.einsum('ij,ij->i', a-d, np.cross(b-d, c-d))) / 6


CLR_LINE =   "                                                               \r"
print ("creating delaunay")

points = np.array(grid.XYZ).T
tri = Delaunay(points)
tets = tri.points[tri.simplices]
#vol = tetrahedron_volume(tets[:, 0], tets[:, 1], 
#                         tets[:, 2], tets[:, 3])

vor = Voronoi(points)

print ("finding all neighbours")

neighbours = find_neighbours(tri)
sys.stdout.write(CLR_LINE)

print ("creating convex hull")
convex_hull = create_convex_hull_list(tri.convex_hull)
sys.stdout.write(CLR_LINE)
tri.close()

print ("length of convex hull: ", len(convex_hull))
print ("Delaunay: done")

num_cols = len(data_ids)
l_max = 2*100*u.au
file = open(output_file, "wb")

file.write(struct.pack("H", grid_id))
file.write(struct.pack("H", num_cols))

for d_ids in data_ids:
    file.write(struct.pack("H", d_ids))

n_neighbours = []
has_neigh = []
for i in range(grid.NPoints):
    p_list = list(neighbours[i])
    sign = int(is_in_hull(convex_hull, i))
    n_neigh = int(len(p_list))
    n_neigh *= sign
    if n_neigh==0: 
        has_neigh.append(False) #Should be False but rejecting points here causes segm.fault in Polaris... I think that after removing cells with no neighbours, the Delaunay grid must be recalculated.
        print ('cell %d has no neighbours, omiting it'%i)
    else: has_neigh.append(True)
    n_neighbours.append(n_neigh)
 
has_neigh = np.array(has_neigh)
id_list = np.arange(grid.NPoints)[has_neigh]
grid.NPoints -= np.sum(~has_neigh)
grid.XYZ = grid.XYZ.T[has_neigh].T
for key in prop: prop[key] = prop[key][id_list]

print('Omitted cells with no neighbours:', np.sum(~has_neigh))


#**************************************************
#REPEATING VORONOI WITHOUT CELLS WITH NO NEIGHBOURS
#**************************************************
CLR_LINE =   "                                                               \r"
print ("creating delaunay")

points = np.array(grid.XYZ).T
tri = Delaunay(points)
tets = tri.points[tri.simplices]
#vol = tetrahedron_volume(tets[:, 0], tets[:, 1], 
#                         tets[:, 2], tets[:, 3])

vor = Voronoi(points)

print ("finding all neighbours")

neighbours = find_neighbours(tri)
sys.stdout.write(CLR_LINE)

print ("creating convex hull")
convex_hull = create_convex_hull_list(tri.convex_hull)
sys.stdout.write(CLR_LINE)
tri.close()

print ("length of convex hull: ", len(convex_hull))
print ("Delaunay: done")

n_neighbours = []
has_neigh = []
for i in range(grid.NPoints):
    p_list = list(neighbours[i])
    sign = int(is_in_hull(convex_hull, i))
    n_neigh = int(len(p_list))
    n_neigh *= sign
    if n_neigh==0: 
        has_neigh.append(False) #Should be False but rejecting points here causes segm.fault in Polaris... I think that after removing cells with no neighbours, the Delaunay grid must be recalculated.
        print ('cell %d has no neighbours, omiting it'%i)
    else: has_neigh.append(True)
    n_neighbours.append(n_neigh)
 
has_neigh = np.array(has_neigh)
id_list = np.arange(grid.NPoints)[has_neigh]
grid.NPoints -= np.sum(~has_neigh)

print('Omitted cells with no neighbours:', np.sum(~has_neigh))


#****************************
np.savetxt('pre_datatab.dat', np.array([prop['dens_H2'], prop['vel_x'], prop['vel_y'], prop['vel_z']]).T, fmt='%.6e')

file.write(struct.pack("d", grid.NPoints))
file.write(struct.pack("d", l_max))

pos_counter = 0
volume = np.zeros(grid.NPoints)+1.0

for i in id_list:
    x = grid.XYZ[0][i]
    y = grid.XYZ[1][i]
    z = grid.XYZ[2][i]
    n_gas = prop['dens_H2'][i]
    
    vor_region_i = np.array(vor.regions[vor.point_region[i]]) #Indices of vertices forming region i.
    vor_region_valid = vor_region_i[vor_region_i!=-1] #-1 indicates vertex outside the Voronoi diagram.
    verts_region_valid = vor.vertices[vor_region_valid]
    delau_i = Delaunay(verts_region_valid) #Delaunay triangulation within voronoi cell i
    tets_i = delau_i.points[delau_i.simplices] #Vertices of Delaunay tetrahedra within vor cell i
    vol = tetrahedron_volume(tets_i[:, 0], tets_i[:, 1], 
                             tets_i[:, 2], tets_i[:, 3]).sum()
    volume[pos_counter] = vol

    file.write(struct.pack("f", x))
    file.write(struct.pack("f", y))
    file.write(struct.pack("f", z))

    file.write(struct.pack("d", vol))
    file.write(struct.pack("f", n_gas))
    file.write(struct.pack("i", int(n_neighbours[i])))

    p_list = list(neighbours[i])
    for j in range(0, int(abs(n_neighbours[i]))):
        tmp_n=int(p_list[j])
        file.write(struct.pack("i", tmp_n))

    pos_counter+=1
    if pos_counter % 10000 == 0:
        sys.stdout.write('Writing voronoi grid: ' + str(100.0 * pos_counter / grid.NPoints) + '%    \r')
        sys.stdout.flush()

sys.stdout.write(CLR_LINE)
print ("DONE")

print ('Ellapsed time:', time.time()-t0)
#******************

