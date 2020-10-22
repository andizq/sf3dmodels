#******************
#sf3dmodels modules
#******************
from sf3dmodels import Model, Plot_model
import sf3dmodels.utils.units as u
import sf3dmodels.rt as rt
from sf3dmodels.model import disc
from sf3dmodels.grid import Grid, fillgrid
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

data_ids_dict = {0: 'dens_gas', 2: 'temp_dust', 3: 'temp_gas'}
struct_dtypes = {'H': np.ushort, 'i': np.int32, 'f': np.float32, 'd': np.float64}

def progress(percent=0, width=30):
    left = width * percent // 100 
    right = width - left
    print('\r[', '#' * left, ' ' * right, ']',
          f' {percent:.0f}%',
          sep='', end='', flush=True)

def read_block(f, dtype=np.float64, count=None):
    return np.fromfile(f, dtype, count)
    
def read_temp(filename):

    f = open(filename, mode='rb')
    print ('Loading file %s' %filename)
    
    grid_id, = np.fromfile(f, struct_dtypes['H'], 1)
    num_cols, = np.fromfile(f, struct_dtypes['H'], 1)
    data_ids = []
    for i in range(num_cols):
        data_ids.append(np.fromfile(f, struct_dtypes['H'], 1)[0])
    npoints, = np.fromfile(f, struct_dtypes['d'], 1).astype('int')
    l_max, = np.fromfile(f, struct_dtypes['d'], 1) #Max extent along longest axis
    
    data = {'xyz': [], 'props': [], 'volume': [], 'n_neighbours': [], 'neighbours': []}
    prop_tags = [data_ids_dict[id] for id in data_ids]
    data.update({data_ids_dict[id]: [] for id in data_ids})
    
    offset=0
    for n in range(npoints):
        data['xyz'].append(np.fromfile(f, struct_dtypes['f'], 3, offset=offset))
        data['volume'].append(np.fromfile(f, struct_dtypes['d'], 1)[0])
        data['props'].append(np.fromfile(f, struct_dtypes['f'], num_cols))
        n_neighbours = abs(np.fromfile(f, struct_dtypes['i'], 1)[0])
        data['n_neighbours'].append(n_neighbours)
        #data['neighbours'].append(np.fromfile(f, struct_dtypes['i'], n_neighbours))
        offset = n_neighbours*np.dtype(struct_dtypes['i']).itemsize
        if n%1000==0: progress(int(100*n/npoints))
    progress(100)
    sys.stdout.write('\n')
    data['xyz'] = np.array(data['xyz'])
    data['props'] = np.array(data['props'])

    return prop_tags, npoints, data
    

t0 = time.time()

filename = "grid_temp.dat"
prop_tags, npoints, data = read_temp(filename)

npoints = npoints
data['xyz'] = data['xyz'][0:npoints]
data['props'] = data['props'][0:npoints]
abund_co = np.where(data['props'][:,2]>19.0, 1e-5, 1e-8) #CO Freeze-out

pre_data = np.loadtxt('pre_datatab.dat')

grid = Model.Struct(XYZ = data['xyz'].T, NPoints= npoints)
prop = {'dens_H2': data['props'][:,0],
        'temp_dust': data['props'][:,1],
        'temp_gas': data['props'][:,2],
        'vel_x': pre_data[:,1],
        'vel_y': pre_data[:,2],
        'vel_z': pre_data[:,3],
        'abundance_0': abund_co,
        }

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
sp = canvas3d.scatter_random(grid, prop['dens_H2']/1e6,  weight, prop_color=prop['temp_dust'], GRID_unit=u.au, power=0, NRand=10000, prop_min=1.0, #function arguments
                             marker = '+', cmap = 'jet', s = 3, edgecolors = 'none', vmax=300, norm = colors.Normalize()) #Scatter kwargs
cbar = plt.colorbar(sp)
cbar.ax.set_ylabel(r'Dust temperature [K]')
ax.set_xlabel('au')
plt.savefig('img_H2_temperature_sph.png')
plt.show()
