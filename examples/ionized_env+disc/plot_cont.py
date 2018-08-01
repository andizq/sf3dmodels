import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as patch
from matplotlib import ticker
import numpy as np
from astropy.io import fits
from copy import copy
import sys

id = "CONV_noise"
tag = "img_HII_env+disc_"+id

#------------------
#LOADING DATA
#------------------
if id: image_cont_349 = "img_HII_env+disc-"+id+".fits"
else: image_cont_349 = "img_HII_env+disc.fits"
data_349, header_349 = fits.getdata(image_cont_349, header=True)
if id == "CONV_noise": pass
else: data_349 = data_349[0]

#-------------------------
#RESOLUTION AND IMAGE SIZE
#-------------------------
distance = 4000 # parsecs
x_lim = 4000 #AUs
extent = [-x_lim, x_lim, -x_lim, x_lim]

#-----------
#CONT. DATA
#-----------
Npix_349 = header_349['NAXIS1']
resolution_349 = abs(header_349['CDELT1']) * 3600 # in arcsecs
#x_lim = (round (Npix / 2.) - x_init) * resolution * distance 
x_inits = int(round (round (Npix_349 / 2.) - (x_lim / (resolution_349 * distance)) ))
x_ends = Npix_349 - x_inits

#-----
#BEAM
#-----
header = header_349
i = 0
ellipse = 0
flag_beam = False
try:
    if header['BMAJ']: pass
    a, b = np.array([header['BMAJ'] , header['BMIN']]) * 3600 * distance
    if id == "CONV": f_c, e_c = 'gray', 'white'
    elif id == "CONV_noise": f_c, e_c = 'black', 'white'
    ellipse = patch.Ellipse(xy = (-3500,-3500), angle = 90 + header['BPA'], 
                            width = a, height = b, linewidth = 1, 
                            fill = True, facecolor = f_c,  edgecolor = e_c)
    flag_beam = True

except KeyError: pass
#-----
#-----

#-----------------------------
#COLORS IN LIMITS OF COLORBARS
#-----------------------------

colormm = plt.cm.hot
palette_hot = copy(colormm)
palette_hot.set_over('red', 0.8)#colormm(255), 0.8)
if id == "CONV_noise": palette_hot.set_under('gray', 0.9)#colormm(0), 0.8)
else: palette_hot.set_under('black', 1.0)#colormm(0), 0.8)

colormm = plt.cm.cubehelix_r
palette_helix = copy(colormm)
palette_helix.set_over('black', 1)#colormm(255), 0.8)
palette_helix.set_under('white', 1)#colormm(0), 0.8)
#-----------------------------
#-----------------------------

#---------
#PLOTTING
#---------
images = data_349

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 4.5))
if id == "CONV_noise": ax.set_facecolor('darkgray')
else: ax.set_facecolor('black')

im = np.ones((3,3))

#plt.subplots_adjust(bottom = 0.01, top = 0.99, left = 0.1, right = 0.9, hspace = 0.1, wspace = 0.3)


vmax = images.max()
vmin = vmax * 1e-3 

im = ax.imshow(images[::-1][x_inits : x_ends, x_inits : x_ends], 
               cmap = palette_hot,
               norm=colors.LogNorm(vmin = vmin, vmax = vmax), 
               extent = extent)
if flag_beam: ax.add_artist(ellipse)

#-----------------
#LABELS AND TITLES
#-----------------

if flag_beam: cbar_label = r'log$_{10}$($I_\nu$ [Jy beam$^{-1}$])'
else: cbar_label = r'log$_{10}$($I_\nu$ [Jy pixel$^{-1}$])'
cbar_extend = 'min'

ax_title = r'$\nu$ = 33 GHz, $S_\nu$ = 103 $\rm{mJy}$'            
if id == "CONV_noise": ax_title += r', $\sigma$ = 6 $\rm{\mu Jy}/bm$'

        
cbar = fig.colorbar(im, extend = cbar_extend, ax = ax, aspect = 15,
                    orientation = "vertical", pad = 0.02, shrink = 1.0)
    
cbar.set_label(cbar_label, labelpad = 2, fontsize = 8)
cbar.ax.tick_params(labelsize = 8)

ax.set_title(ax_title, fontsize = 10)
ax.set_ylabel('AU', fontsize = 10)
ax.set_xlabel('AU', fontsize = 10)

ax.tick_params(axis='both', labelsize = 7)

#-----------------
#-----------------
plt.savefig("%s.png"%tag, facecolor = "white", dpi = 500)
plt.show()
#---------
#---------
