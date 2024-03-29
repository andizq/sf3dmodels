from __future__ import print_function
import sys
import time
import random
import inspect
import warnings
import itertools
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#Needed to handle 3d projections, even if not called:
from mpl_toolkits.mplot3d import Axes3D 
from . import BuildGlobalGrid as BGG
from . import Model

"""
TINY_SIZE = 5
SMALL_SIZE = 10
MEDIUM_SIZE = 15
BIGGER_SIZE = 20

matplotlib.rcParams['font.family'] = 'monospace'
matplotlib.rcParams['font.weight'] = 'normal'
matplotlib.rcParams['lines.linewidth'] = 1.5
matplotlib.rcParams['axes.linewidth'] = 1.8

matplotlib.rcParams['xtick.minor.width']=1.
matplotlib.rcParams['ytick.minor.width']=1.
matplotlib.rcParams['xtick.major.width']=1.3
matplotlib.rcParams['ytick.major.width']=1.3

matplotlib.rcParams['xtick.minor.size']=2.8
matplotlib.rcParams['ytick.minor.size']=2.8
matplotlib.rcParams['xtick.major.size']=5.5
matplotlib.rcParams['ytick.major.size']=5.5

matplotlib.rc('font', size=MEDIUM_SIZE) 
matplotlib.rc('axes', titlesize=MEDIUM_SIZE, labelsize=MEDIUM_SIZE)    
matplotlib.rc('xtick', labelsize=SMALL_SIZE)   
matplotlib.rc('ytick', labelsize=SMALL_SIZE)   
matplotlib.rc('legend', fontsize=SMALL_SIZE)   
matplotlib.rc('figure', titlesize=BIGGER_SIZE)  
"""

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class InputError(Error):
    """Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message        
    def __str__(self):
        return self.expression+self.message

class MakeCanvas(object):
    """
    Base class for 1d,2d,3d canvases
    """
    def __init__(self, fig=None, ax=None, ax_rect=[0.05,0.05,0.9,0.9], fig_kw={}, ax_kw={}):
        if fig is None: fig = plt.figure(**fig_kw)
        self.fig= fig
        if ax is None: ax = fig.add_axes(ax_rect, **ax_kw)
        self.ax = ax

class Canvas3d(MakeCanvas):
    def __init__(self, fig=None, ax=None, ax_rect=[0.05,0.05,0.9,0.9], fig_kw={}, ax_kw={}):
        ax_kw.update({'projection': '3d'})
        if ax is not None and ax.name != '3d':
            raise InputError("__init__():\t", "input Axes instance is not a 3d projection, you can do for example ax=plt.axes(projection='3d')")                   
        super(Canvas3d, self).__init__(fig, ax, ax_rect, fig_kw, ax_kw)
        self.fig_kw = fig_kw
        self.ax_kw = ax_kw
        self.ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        self.ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        self.ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))    

    def scatter_random(self, GRID, prop, weight, NRand = 1000, prop_color = None, prop_min = None, GRID_unit=1.0, power = 0.5, count_max=50, **scatter_kw):
        t0 = time.time()
        print ('Plotting 3D model with %d random grid points...'%NRand)
        print ("Using '(prop_i/weight)**power > random[0,1] ?' to reject/accept the i-th point of \n  the randomly-selected sample. You set power=%.1f and weight=%.1e"%(power, weight))
        
        """
        defaults = dict(marker = '+', cmap = 'hot', s = 3, edgecolors = 'none', vmin = None, vmax = None, norm = None) #scatter3d kwargs
        for key_def in defaults.keys():
            if key_def in kwargs.keys(): continue
            else: kwargs[key_def] = defaults[key_def]
        """
        N = GRID.NPoints
        prop = np.asarray(prop)
        
        index_pop = np.arange(N) #All the population
        if prop_min is not None: index_pop = index_pop[prop > prop_min] #Can reject e.g. zero cells 

        indices = []
        count = 0

        for _ in itertools.repeat(None, NRand):
            rand = random.random()
            while 1:
                index = np.random.choice(index_pop) #Selects 1 point from the given list
                val = (prop[index]/weight)**power
                if val > rand:
                    indices.append(index)
                    count = 0
                    break
                count += 1
                if count == count_max: #If after 50 points the algorithm has not selected any, pick another point randomly.
                    warnings.warn("The algorithm is struggling to accept grid points. If it's taking too long try increasing the weighting 'weight' or reducing the exponent 'power'")
                    count = 0
                    rand = random.random()

        indices = np.array(indices)
        x,y,z = [xi[indices]/GRID_unit for xi in GRID.XYZ]
    
        if prop_color is None: prop_color = prop

        sp = self.ax.scatter(x,y,z, 
                             c = prop_color[indices], 
                             **scatter_kw)
    
        print ("Ellapsed time from %s: %.2f s"%(inspect.stack()[0][3], time.time()-t0))
        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------\n-------------------------------------------------')
        return sp 

    def vectors(self, x,y,z, vx,vy,vz, GRID_unit=1, length=1, arrow_length_ratio=0.1, **quiver_kw):
        """
        quiver_kw: Axes3D.quiver kwargs: https://matplotlib.org/mpl_toolkits/mplot3d/tutorial.html#quiver
        """
        t0 = time.time()
        print ('Plotting 3D vectorial field...')        
        quiver_kw.update({'length': length, 'arrow_length_ratio': arrow_length_ratio})
        x,y,z = [np.asarray(xi)/GRID_unit for xi in [x,y,z]]
        vp = self.ax.quiver(x,y,z, vx,vy,vz, **quiver_kw)
        print ("Ellapsed time from %s: %.2f s"%(inspect.stack()[0][3], time.time()-t0))
        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------\n-------------------------------------------------')
        return vp 
    

"""
class Canvas2d(MakeCanvas):
    def __init__():
        self._vectors_arrow = mpatches.ArrowStyle.CurveFilledB(head_length=0.25, head_width=0.1)

    @property
    def vectors_arrow(self):
        return self._vectors_arrow
    @vectors_arrow.setter
    def vectors_arrow(self, val):
        if not isinstance(val, mpatches.ArrowStyle):
            raise ValueError("Please use a matplotlib mpatches.ArrowStyle instance")
        self._vectors_arrow = val
    
    def vectors(x,y,z, vx,vy,vz, max_length=10., GRID_unit=1, **annotate_kw):
        t0 = time.time()
        print ('Plotting 2D vectorial field...')
        arrow = self._vectors_arrow
        ax[i].annotate("", xy=(xf[n], yf[n]), xytext=(x[n], y[n]), arrowprops=dict(arrowstyle=arrow, color='k',linewidth=1.0))
        print ("Ellapsed time from %s: %.2f s"%(inspect.stack()[0][3], time.time()-t0))
        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------\n-------------------------------------------------')
        return vp 
"""        

def scatter3D(GRID, prop, weight, 
              colordim = [False], NRand = 1000, axisunit = 1.0, power = 0.6,
              colorlabel = '', output = 'figscatter.png', show = True, 
              colorscale = 'uniform', 
              azim = None, elev = None, xlim=None, ylim=None, zlim=None,
              **kwargs):
    
    t0 = time.time()
    print ('Plotting 3D model with %d random-weighted points...'%NRand)

    defaults = dict(marker = '+', cmap = 'hot', s = 3, edgecolors = 'none', vmin = None, vmax = None, norm = None) #scatter3d kwargs
    for key_def in defaults.keys():
        if key_def in kwargs.keys(): continue
        else: kwargs[key_def] = defaults[key_def]

    x,y,z = GRID.XYZ
    #r = GRID.rRTP[0]
    NTotal = GRID.NPoints
    unit = axisunit
    scale = colorscale

    palette_c = cm.get_cmap(kwargs['cmap']) #getattr(cm , kwargs['cmap'])

    #population = range(NTotal) #All the population
    population = list( np.where(abs(prop) > 2.725)[0] ) # > 1e3. #Rejecting zero cells 

    indices = []
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

            if count == 50: #If after 50 points the algorithm has not selected any, pick another point randomly.
                count = 0
                rand = random.random()

    #colors = palette_c(np.linspace(0, 1, NRand))
    indices = np.array(indices).T[0]
    x,y,z = x[indices], y[indices], z[indices]
    #r = r[indices]
    
    if not np.array(colordim).any(): colordim = prop

    if kwargs['norm'] is not None or scale == 'uniform': prop2plot = np.sort(colordim[indices]) 
    elif scale == 'log': prop2plot = np.sort(np.log10(colordim[indices]))

    ind2plot = np.argsort(colordim[indices])

    x,y,z = x[ind2plot]/unit, y[ind2plot]/unit, z[ind2plot]/unit

    fig = plt.figure()
    #fig.clf()
    ax = plt.gca(projection='3d') 
    sp = ax.scatter(x,y,z, 
                    c = prop2plot, 
                    **kwargs)
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
    if xlim is not None: 
        xlim=np.asarray(xlim)
        xlim=np.where(xlim != None, xlim/unit, None)
    if ylim is not None: 
        ylim=np.asarray(ylim)
        ylim=np.where(ylim != None, ylim/unit, None)
    if zlim is not None: 
        zlim=np.asarray(zlim)
        zlim=np.where(zlim != None, zlim/unit, None)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_zlim(zlim)

    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

    ax.view_init(azim = azim, elev = elev)
    print ('3D camera azimuth: %.1f deg'%ax.azim)
    print ('3D camera elevation: %.1f deg'%ax.elev)

    cbar = plt.colorbar(sp)
    cbar.ax.set_ylabel('%s'%colorlabel)

    if output == '': output = 'figure.png'
    
    plt.tight_layout()
    plt.savefig(output, dpi = 1000)    
    print ('Image saved as %s'%output)
    print ("Ellapsed time from %s: %.2f s"%(inspect.stack()[0][3], time.time()-t0))

    if show:
        print ('Showing computed image...')
        plt.show()
    else:
        print ('The show-image mode is off!')
        plt.close()

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')
    return prop2plot

def plane2D(GRID, prop, axisunit = 1.0, plane = {'z': 0}, rot_dict = False, 
            colorlabel = '', output = 'figplane.png', show = True, auto = True, **kwargs):
    
    unit = axisunit
    i,j,k = 0,0,0

    key = list(plane.keys())
    if len(key) > 1: sys.exit('ERROR: Please define a single plane to plot')
    key = key[0]
    if key == 'x': ind = 0
    elif key == 'y': ind = 1
    else: ind = 2
    
    coord = GRID.XYZ[ind]
    dim2plot, = np.where(np.array([0,1,2]) != ind) #perpendicular dimensions to plane ind

    value = plane[key]
    indices, = np.where(coord == value)
    title = '%s = %d'%(key,value)

    lims = np.array([np.array( [min(GRID.XYZcentres[dd]), max(GRID.XYZcentres[dd])] ) / unit for dd in dim2plot]).flatten()
    defaults = dict(cmap = 'hot', extent = lims)
    for key_def in defaults.keys():
        if key_def in kwargs.keys(): continue
        else: kwargs[key_def] = defaults[key_def]

    if len(indices) == 0:
        domain_sorted = np.sort(GRID.XYZcentres[ind])        
        tmp, = np.where(domain_sorted < value)
        ind2 = tmp[-1] #closest below
        ind1 = ind2 + 1 #closest above
        value2 = min([domain_sorted[i] for i in [ind1,ind2]], 
                     key = lambda x: abs(x-value))
        indices, = np.where(coord == value2)
        title = 'SHIFTED to %s = %d'%(key,round(value2/unit))
    
    if rot_dict:
        GRID_rot = Model.Struct(**{'XYZ': None})
        GRID_rot.XYZ = [GRID.XYZ[i][indices] for i in range(3)]
        GRID_rot.NPoints = len(GRID_rot.XYZ[0])
        newProperties = Model.ChangeGeometry(GRID_rot, rot_dict = rot_dict)
        GRID_rot.XYZ = newProperties.newXYZ
        X, Y, Z = GRID_rot.XYZ #Rotated plane
        Xgrid, Ygrid, Zgrid = GRID.XYZcentres #Original grid
        Nx, Ny, Nz = GRID.Nodes
        Num = []
        for x,y,z in zip(X,Y,Z):
            i = BGG.mindistance(x,Xgrid,Nx)
            j = BGG.mindistance(y,Ygrid,Ny)
            k = BGG.mindistance(z,Zgrid,Nz)
            Num.append(i*(Ny)*(Nz)+j*(Nz)+k) #ID for the Global Grid
            
        Num = sorted(set(Num), key=lambda x: Num.index(x)) #Keep the order and delete duplicates
        Num = np.array(Num)
        prop_rot = prop[Num]

        dim0grid, dim1grid = [GRID.XYZ[i][Num] for i in dim2plot]#Final grid lists in each dimension
        perpdim0, perpdim1 = [sorted(list(set(dim))) for dim in [dim0grid, dim1grid]]
        prop2plot = prop_rot
        cells = len(prop2plot)
        prop_2d = np.zeros([len(perpdim1), len(perpdim0)])
        grid_2d = [[0 for i in range(len(perpdim0))] for i in range(len(perpdim1))]
        i,j,k,i0 = 0,0,0,0
        tmp_u, tmp_v, tmp_prop, list_prop, current_prop = None, None, prop2plot[k], [], None

        for u in perpdim0:
            for v in perpdim1:
                while (dim0grid[k] == tmp_u and dim1grid[k] == tmp_v): #Repeated X,Y
                    if tmp_v == perpdim1[-1]: i0 = 1 #If the repeated cell is in the last Y, turn on i0
                    else: i0 = 0
                    list_prop.append(prop2plot[k])
                    prop_2d[j-1,i-i0] = np.mean(np.append(current_prop, list_prop)) #Fixing the j-1 cell
                    k+=1
                    #tmp_prop = prop2plot[k]
                    #tmp_u, tmp_v = u, v

                prop_2d[j,i] = prop2plot[k] #tmp_prop
                grid_2d[j][i] = (u,v) 
            
                j+=1
                k+=1
                #if k == len(prop2plot): continue
                #tmp_prop = prop2plot[k]
                tmp_u, tmp_v = u, v
                list_prop, current_prop = [], prop_2d[j-1,i]
            
            j=0
            i+=1

        kwargs['extent'] = np.array([perpdim0[0], perpdim0[-1], perpdim1[0], perpdim1[-1]]) / unit
    
    else:
        
        perpdim0, perpdim1 = np.array(GRID.XYZcentres)[dim2plot]
        prop_2d = np.zeros([len(perpdim1), len(perpdim0)])
        grid_2d = [[0 for i in range(len(perpdim0))] for i in range(len(perpdim1))]
        i,j,k = 0,0,0
        prop2plot = prop[indices]
        print (np.shape(prop2plot), np.shape(prop_2d), dim2plot)
        for u in perpdim0:
            for v in perpdim1:
                prop_2d[j,i] = prop2plot[k] 
                grid_2d[j][i] = (u,v) 
            
                j+=1
                k+=1
            j=0
            i+=1

    if auto:    
        fig = plt.figure()
        ax = plt.axes()
        im = ax.imshow(prop_2d, 
                       #norm=colors.LogNorm(vmin = 5e8, vmax = 5e12), 
                       #cmap = palette,
                       **kwargs)
                            
        ax.set_title(title)
        xyz = ['X', 'Y', 'Z']
        #ax.set_xlabel(xyz[dim2plot[0]])
        #ax.set_ylabel(xyz[dim2plot[1]])
        #ax.set_ylabel('au')
        #ax.set_xlabel('au')
        
        min_cbar, max_cbar = im.get_clim()
        min_prop, max_prop = np.min(prop_2d), np.max(prop_2d) 

        if min_prop < min_cbar and max_prop > max_cbar: cbar_extend = 'both'
        elif min_prop < min_cbar: cbar_extend = 'min'
        elif max_prop > max_cbar: cbar_extend = 'max'
        else: cbar_extend = 'neither'
        
        cbar = plt.colorbar(im, pad = 0.01, extend = cbar_extend)
        cbar.ax.set_ylabel('%s'%colorlabel, rotation = 270, labelpad = 20)

        #plt.autoscale(tight=True)
        #plt.tight_layout()
        
        if output == '': output = 'figure.png'

        print ('Saving image in %s'%output)
        plt.savefig(output, dpi = 1000, bbox_inches='tight')
    
        if show: plt.show()
        else:
            print ('The show-image mode is off!')
            plt.close()
    
    else: return prop_2d
            
    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')
