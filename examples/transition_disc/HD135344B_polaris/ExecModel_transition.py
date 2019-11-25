"""
Transition Disc Example
=======================

This example uses the gas and dust density profiles inferred for the transitional disc HD135344B.
"""
#******************
#sf3dmodels modules
#******************
from sf3dmodels import Model, Plot_model
import sf3dmodels.utils.units as u
import sf3dmodels.rt as rt
from sf3dmodels.model import disc
from sf3dmodels.grid import Grid
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
from scipy.spatial import Delaunay, ConvexHull, Voronoi
from collections import defaultdict
import struct
import bisect

t0 = time.time()
#******************
#Model parameters
#******************
kwargs_dens={'dn_cav': 1e-5, 'n_cav': 2e17, 'R_cav': 30*u.au, 'power':-0.9, 'phi_stddev': 2*np.pi/3, 'phi_mean': -np.pi/2}
kwargs_dtg={'dn_cav': 1e-2, 'n_cav': 0.01, 'R_cav': 40*u.au, 'power': 0} #, 'phi_stddev': 2*np.pi/3, 'phi_mean': -np.pi/2}
R_disc = 70*u.au

#******************
#Grid creation
#******************
transition = disc.Transition()
init_grid = Grid()
grid = init_grid.random(func=transition.powerlaw_cavity, power=0.4, r_size=R_disc, kwargs_func=kwargs_dens, normalization=1e16, npoints=100000)

#*****************************************
#Computing relevant physical distributions
#*****************************************
density = transition.powerlaw_cavity(grid=grid, **kwargs_dens)
gtdratio = 1./transition.powerlaw_cavity(grid=grid, **kwargs_dtg) #Converting to gas-to-dust ratio for LIME
temperature = np.zeros(grid.NPoints) + 100. #Constant temperature 100 K

prop = {'dens_H2': density,
        'gtdratio': gtdratio,
        'temp_gas': temperature}

#*******************************************************************************
#Writing file for self-consistently calculation of dust temperature with POLARIS
#*******************************************************************************
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

def convex_hull_volume(pts):
    ch = ConvexHull(pts)
    dt = Delaunay(pts[ch.vertices])
    tets = dt.points[dt.simplices]
    return np.sum(tetrahedron_volume(tets[:, 0], tets[:, 1],
                                     tets[:, 2], tets[:, 3]))
def tetrahedron_volume(a, b, c, d):
    return np.abs(np.einsum('ij,ij->i', a-d, np.cross(b-d, c-d))) / 6

CLR_LINE =   "                                                               \r"
print ("creating delaunay")

points = np.array(grid.XYZ).T
tri = Delaunay(points)
tets = tri.points[tri.simplices]
#vol = tetrahedron_volume(tets[:, 0], tets[:, 1], 
#                         tets[:, 2], tets[:, 3])
#convHull = ConvexHull(points) #This computes the smallest region enclosing the whole set of points
#volume = convHull.volume/grid.NPoints #Averaged volume

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
l_max = 2*R_disc
file = open(output_file, "wb")

file.write(struct.pack("H", grid_id))
file.write(struct.pack("H", num_cols))

for d_ids in data_ids:
    file.write(struct.pack("H", d_ids))

file.write(struct.pack("d", grid.NPoints))
file.write(struct.pack("d", l_max))

pos_counter = 0
volume = np.zeros(grid.NPoints)+1.0
for i in range(grid.NPoints):
    x = grid.XYZ[0][i]
    y = grid.XYZ[1][i]
    z = grid.XYZ[2][i]
    r = grid.rRTP[0][i]
    #dens = data['rho'][i]*arepoDensity #in g/ccm
    n_gas = prop['dens_H2'][i]

    vor_region_i = np.array(vor.regions[vor.point_region[i]]) #Indices of vertices forming region i.
    vor_region_valid = vor_region_i[vor_region_i!=-1] #-1 indicates vertex outside the Voronoi diagram.
    verts_region_valid = vor.vertices[vor_region_valid]
    #convHull_i = ConvexHull(verts_region_valid) #Failing to get convexhull when vertices are close to a common plane. Computing the delaunay and summing the volume of tetrahedrons leads to same results apparently.
    #volume_2[i] = convHull_i.volume
    delau_i = Delaunay(verts_region_valid)
    tets_i = delau_i.points[delau_i.simplices]
    vol = tetrahedron_volume(tets_i[:, 0], tets_i[:, 1], 
                             tets_i[:, 2], tets_i[:, 3]).sum()
      
    if r > R_disc-u.au or abs(z)/u.au > 20: vol = 1.0
    else:
        delau_i = Delaunay(verts_region_valid)
        tets_i = delau_i.points[delau_i.simplices]
        vol = tetrahedron_volume(tets_i[:, 0], tets_i[:, 1], 
                             tets_i[:, 2], tets_i[:, 3]).sum()
    volume[i] = vol
        
    file.write(struct.pack("f", x))
    file.write(struct.pack("f", y))
    file.write(struct.pack("f", z))

    file.write(struct.pack("d", vol))

    file.write(struct.pack("f", n_gas))
    
    p_list = list(neighbours[pos_counter])
    sign = int(is_in_hull(convex_hull, pos_counter))
    n_neighbours = int(len(p_list))
    n_neighbours *= sign

    file.write(struct.pack("i", int(n_neighbours)))

    for j in range(0, int(abs(n_neighbours))):
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
