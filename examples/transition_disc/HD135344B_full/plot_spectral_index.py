import numpy as np
from astropy.io import fits
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib
import matplotlib.pyplot as plt
import sys

func_line = lambda x,b,m: b+m*x

bands = ['band6', 'band7', 'band8']
nbands = len(bands)
bands_dict = {'band6': 230.0, 'band7': 340.0, 'band8': 450.0}
img_base = 'img_cont_faceon_%s.fits'
data = np.array([fits.getdata(img_base%band).squeeze() for band in bands])

spec_index = np.loadtxt('spectral_index.txt')
opacity = np.loadtxt('../../opacities_k05_230GHz_B_1_7.tab')
freq = 3e8/(opacity[:,0] * 1e-6) / 1e9 #in GHz

xdata = np.log10([bands_dict[band] for band in bands])

#*********************
#PLOTTING
#*********************
TINY_SIZE = 8
SMALL_SIZE = 10
MEDIUM_SIZE = 15
BIGGER_SIZE = 20

matplotlib.rcParams['font.family'] = 'monospace'
matplotlib.rcParams['axes.linewidth'] = 2.5  
matplotlib.rcParams['xtick.major.width']=2.0
matplotlib.rcParams['ytick.major.width']=2.0
matplotlib.rcParams['xtick.minor.width']=1.5
matplotlib.rcParams['ytick.minor.width']=1.5
matplotlib.rcParams['xtick.major.size']=5.0
matplotlib.rcParams['ytick.major.size']=5.0
matplotlib.rcParams['xtick.minor.size']=2.8
matplotlib.rcParams['ytick.minor.size']=2.8
matplotlib.rc('axes', labelsize=MEDIUM_SIZE+4)
matplotlib.rc('xtick', labelsize=MEDIUM_SIZE)
matplotlib.rc('ytick', labelsize=MEDIUM_SIZE)

fig = plt.figure(figsize=(16,7))
ax0 = fig.add_subplot(121)
im = ax0.imshow(spec_index, origin='lower left', extent=[-70,70,-70,70], vmin=1.0, cmap='nipy_spectral_r')
ax0.set_xlabel('au')
#CS = ax0.contour(spec_index, levels=[1.3,1.4,1.5,1.6], origin='lower', extent=[-70,70,-70,70], lw=3)
#ax0.clabel(CS, inline=1, fontsize=10)
cbar = fig.colorbar(im, ax=ax0, extend='min')
cbar.ax.set_ylabel(r'Spectral index')

ax1 = fig.add_subplot(122)
ax1.plot(freq, opacity[:,1], color = 'k', lw=5)
ax1.text(0.1,0.7,r'$\kappa_{\nu} = \kappa_0 \left(\frac{\nu}{\nu_0}\right)^{\beta}$'+'\n'+r'$\rightarrow\kappa_0=0.5$ cm$^2$ g$^{-1}$'+'\n'+r'$\rightarrow\nu_0=230$ GHz' +'\n'+r'$\rightarrow\beta=1.7$', ha='left', transform=ax1.transAxes, fontsize=MEDIUM_SIZE+2)
ax1.tick_params(axis='y', which='both', labelcolor='k', left=True, labelleft=False, right=True, labelright=True)
ax1.set_xlabel(r'$\nu$ [GHz]')
ax1.yaxis.set_label_position('right')
ax1.set_ylabel(r'$\kappa_{\nu}$')
ax1.grid()
 
plt.savefig('img_spectral_index.png')
plt.show()
