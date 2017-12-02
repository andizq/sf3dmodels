"""Plotting library, by AFIC"""

from __future__ import print_function
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import random
import inspect

AU = 1.495978707e11 #1AU to meters  


def scatter3D_rand(XYZ, prop, NRand = 10000, unit = 1.0, palette = 'hot', scale = 'uniform', output = 'figscatter.png', show = True):

    print ('Plotting 3D model with %d random-weighted points...'%NRand)
    x,y,z = XYZ
    NTotal = len(x)
    
    palette_c = getattr(cm , palette)

#    population = range(NTotal) #All the population!
#    population = list( np.where(abs(prop) > 1.e9)[0] ) #Where the property is not in the "vacuum"
    population = list( np.where(abs(prop) > 2)[0] ) #Where the property is not in the "vacuum"
#    population = list( np.where(prop != 0.)[0] ) #Where the property is not in the "vacuum"

    indices = []
    power = 0.
    weight = np.max(prop)
    print(weight)
    count = 0
    k = 0
    for i in range(NRand):
        rand = random.random()
        while 1:
            index = random.sample(population, 1) #Selects 1 point from the given list
        
            if (prop[index]/weight)**power > rand:
                indices.append(index)
                break
            count += 1

            if count == 50: 
                k += 1
                print ("warning", k)
                count = 0

                rand = random.random()

    colors = palette_c(np.linspace(0, 1, NRand))

    x,y,z = x[indices], y[indices], z[indices]
    
    if scale == 'uniform': prop2plot = np.sort(prop[indices])
    elif scale == 'log': prop2plot = np.sort(np.log10(prop[indices]))

    ind2plot = np.argsort(prop[indices])
    x,y,z = x[ind2plot]/unit, y[ind2plot]/unit, z[ind2plot]/unit


    fig = plt.figure()
    ax = fig.gca(projection='3d')#fig.add_subplot(111, projection = '3d')
    sp = ax.scatter(x,y,z, s = 20, c = prop2plot, cmap = palette)
    
    ax.set_xlabel('X (AU)')
    ax.set_ylabel('Y (AU)')
    ax.set_zlabel('Z (AU)')
        
    cbar = plt.colorbar(sp)
    cbar.ax.set_ylabel('%s scale'%scale)
#    cbar.set_clim(1e16,1e17)

    if output == '':
        print ('No output to save image...')
    else:
        print ('Saving image in %s'%output)
        plt.savefig(output)
    
    if show: plt.show()
    else:
        print ('The show-image mode is off!')
        plt.close()

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')
    

def scatter3D(XYZ, prop, NRand = 10000, unit = 1.0, palette = 'hot', scale = 'uniform', output = 'figscatter.png', show = True):

    print ('Plotting 3D model with %d random points...'%NRand)
    
    x,y,z = XYZ
    NTotal = len(x)
    
    palette_c = getattr(cm , palette)

#    population = range(NTotal) #All the population!
#    population = list( np.where(abs(prop) > 1.e9)[0] ) #Where the property is not in the "vacuum"
    population = list( np.where(abs(prop) > 2)[0] ) #Where the property is not in the "vacuum"
#    population = list( np.where(prop != 0.)[0] ) #Where the property is not in the "vacuum"
    indices = random.sample(population, NRand) #Selects NRand points from the given list
    indices = np.array(indices)

    colors = palette_c(np.linspace(0, 1, NRand))
    x,y,z = x[indices], y[indices], z[indices]
    
    if scale == 'uniform': prop2plot = np.sort(prop[indices])
    elif scale == 'log': prop2plot = np.sort(np.log10(prop[indices]))

    ind2plot = np.argsort(prop[indices])
    x,y,z = x[ind2plot]/unit, y[ind2plot]/unit, z[ind2plot]/unit


    fig = plt.figure()
    ax = fig.gca(projection='3d')#fig.add_subplot(111, projection = '3d')
    sp = ax.scatter(x,y,z, s = 20, c = prop2plot, cmap = palette)
    
    ax.set_xlabel('X (AU)')
    ax.set_ylabel('Y (AU)')
    ax.set_zlabel('Z (AU)')
        
    cbar = plt.colorbar(sp)
    cbar.ax.set_ylabel('%s scale'%scale)
#    cbar.set_clim(1e16,1e17)

    if output == '':
        print ('No output to save image...')
    else:
        print ('Saving image in %s'%output)
        plt.savefig(output)
    
    if show: plt.show()
    else:
        print ('The show-image mode is off!')
        plt.close()

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')
    
    
