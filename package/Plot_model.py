from __future__ import print_function
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import random
import inspect


def scatter3D(GRID, prop, weight, colordim = [False], NRand = 1000, axisunit = 1.0, palette = 'hot',
              colorscale = 'uniform', colorlabel = '', output = 'figscatter.png', show = True):

    print ('Plotting 3D model with %d random-weighted points...'%NRand)
    x,y,z = GRID.XYZ
    r = GRID.rRTP[0]
    NTotal = GRID.NPoints
    unit = axisunit
    scale = colorscale

    palette_c = getattr(cm , palette)

    #population = range(NTotal) #All the population
    population = list( np.where(abs(prop) > 2)[0] ) #Rejecting zero cells 

    indices = []
    power = 0.6
    count = 0

    for i in range(NRand):
        rand = random.random()
        while 1:
            index = random.sample(population, 1) #Selects 1 point from the given list
            val = (prop[index]/weight)**power

            if val > rand:
                indices.append(index)
                count = 0
                k = 0
                break
            count += 1

            if count == 50:
                count = 0
                rand = random.random()

    colors = palette_c(np.linspace(0, 1, NRand))
    indices = np.array(indices).T[0]
    x,y,z = x[indices], y[indices], z[indices]
    r = r[indices]
    
    if not np.array(colordim).any(): colordim = prop
    
    if scale == 'uniform': prop2plot = np.sort(colordim[indices])
    elif scale == 'log': prop2plot = np.sort(np.log10(colordim[indices]))

    ind2plot = np.argsort(colordim[indices])


    x,y,z = x[ind2plot]/unit, y[ind2plot]/unit, z[ind2plot]/unit

    fig = plt.figure()
    ax = fig.gca(projection='3d') 
    sp = ax.scatter(x,y,z, s = 5, c = prop2plot, cmap = palette, marker = '+')
    
    ax.set_xlabel('X (AU)')
    ax.set_ylabel('Y (AU)')
    ax.set_zlabel('Z (AU)')
        
    cbar = plt.colorbar(sp)
    cbar.ax.set_ylabel('%s'%colorlabel)

    if output == '': output = 'figure.png'

    print ('Saving image in %s'%output)
    plt.savefig(output, dpi = 1000)
    
    if show: plt.show()
    else:
        print ('The show-image mode is off!')
        plt.close()

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')
