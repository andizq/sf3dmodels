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

xdata = np.log10([bands_dict[band] for band in bands])


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
matplotlib.rc('axes', labelsize=MEDIUM_SIZE)
matplotlib.rc('xtick', labelsize=MEDIUM_SIZE-2)
matplotlib.rc('ytick', labelsize=MEDIUM_SIZE-2)



fig = plt.figure(figsize=(16,7))
ax0 = fig.add_subplot(121)
im = ax0.imshow(spec_index, origin='lower left', extent=[-70,70,-70,70], cmap='nipy_spectral')
ax0.set_xlabel('au')
cbar = fig.colorbar(im, ax=ax0)
cbar.ax.set_ylabel(r'Spectral index')
plt.savefig('img_spectral_index.png')
plt.show()
