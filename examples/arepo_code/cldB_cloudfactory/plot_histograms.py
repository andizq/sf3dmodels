import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arrow
import matplotlib.colors as mc
from copy import copy

#**********************
#PLOTTING 2D HISTOGRAMS
#**********************

def plot_all(data, yn, gasmass, dcell, rr, iref):
    #dcell vs yn
    plt.hist2d(np.log10(yn[iref]),np.log10(dcell[iref]*100.)
               ,bins=400
               ,norm = mc.PowerNorm(0.3)
               #,range=[[-2, 4], [-0.5, 2]]
               )
    plt.title('Cell size vs Number density')
    cbar = plt.colorbar()
    cbar.set_label('Counts', labelpad = 15)
    plt.grid(b=True, which='major', linestyle='--')
    #plt.grid(b=True, which='minor', linestyle=':')
    plt.minorticks_on()
    plt.xlabel('log yn [cm$^{-3}$]',fontsize=12)
    plt.ylabel('log dcell [pc]',fontsize=12)
    plt.savefig('cellsize_numdens-AREPOgrid.png')
    plt.close()

    #dcell vs xCO
    iref2 = np.where((rr!=None) & (data['chem'][:,2] > 0.0))
    plt.hist2d(np.log10(data['chem'][:,2][iref2]),np.log10(dcell[iref2]*100.)
               ,bins=400
               ,norm = mc.PowerNorm(0.3)
               )
    plt.xlim(-15,None)
    plt.title('Cell size vs Number density')
    cbar = plt.colorbar()
    cbar.set_label('Counts', labelpad = 15)
    plt.grid(b=True, which='major', linestyle='--')
    plt.minorticks_on()
    plt.xlabel('log xCO [nCO/nH]',fontsize=12)
    plt.ylabel('log dcell [pc]',fontsize=12)
    plt.savefig('cellsize_xco-AREPOgrid.png')
    plt.close()

    #u_therm vs yn
    colormm = plt.cm.cool
    palette = copy(colormm)
    palette.set_under('darkgray',0.8)

    x_, y_ = np.log10(yn[iref]),np.log10(data['u_therm'][iref])
    x_0,x_1, y_0,y_1 = x_.min(),x_.max(), y_.min()-2,y_.max()

    fig, ax = plt.subplots()
    ax.set_facecolor((0.15,0.15,0.15))

    im = ax.hist2d(x_,y_
                   ,bins=400
                   ,range=[[x_0,x_1], [y_0,y_1]]
                   ,norm = mc.LogNorm(vmin=14.0)
                   ,cmap = palette
                   )

    ax.add_patch(Arrow(4.5,-2.7,0,1, width=0.5, ec='red', fc='yellow'))

    plt.title('Thermal energy per unit mass vs Number density', fontsize=10)

    cbar = fig.colorbar(im[3], extend='min') 
    #hist2d returns (counts, xedges, yedges, Image). Then the (mappable) image is on the 4th entry (index 3)
    min_tick = min(im[3].get_clim()) #cmin
    cbar.set_ticks([min_tick] + list(cbar.get_ticks()))
    cbar.set_label('Counts', labelpad = 15)

    plt.grid(b=True, which='major', linestyle='--')
    plt.minorticks_on()
    plt.xlabel('log yn [cm$^{-3}$]',fontsize=12)
    plt.ylabel('log $\mu_{th}$ [?]',fontsize=12)
    plt.savefig('utherm_numdens-AREPOgrid.png')
    plt.close()

    #gasmass vs yn
    plt.hist2d(np.log10(yn[iref]),np.log10(gasmass[iref])
               ,bins = 400
               ,norm = mc.PowerNorm(0.4)           
               )
    plt.title('Gas mass vs Number density')
    cbar = plt.colorbar()
    cbar.set_label('Counts', labelpad = 15)
    plt.grid(b=True, which='major', linestyle='--')
    plt.minorticks_on()
    plt.xlabel('log yn [cm$^{-3}$]',fontsize=12)
    plt.ylabel('log gasmass [M$\odot$]',fontsize=12)
    plt.savefig('gasmass_numdens-AREPOgrid.png')
    plt.close()

