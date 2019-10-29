import numpy as np
from astropy.io import fits
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import sys

bands = ['band6', 'band7', 'band8']
nbands = len(bands)
bands_dict = {'band6': 230.0, 'band7': 340.0, 'band8': 450.0}
img_base = 'img_cont_faceon_%s.fits'
data = np.array([np.log10(fits.getdata(img_base%band).squeeze()) for band in bands])

func_line = lambda x,b,m: b+m*x

nrows, ncols = data[0].shape    
xdata = np.log10([bands_dict[band] for band in bands])
spec_index = np.zeros((nrows,ncols))

for row in range(nrows):
    for col in range(ncols):
        ydata = data[:,row,col]#np.array([data[i][row][col] for i in range(nbands)])
        if np.all(ydata == ydata[0]): continue
        slope, intercept, r_value, p_value, std_err = stats.linregress(xdata,ydata) #pt, pv = curve_fit(func_line, xdata, ydata)
        spec_index[row][col] = slope #pt[1]

np.savetxt('spectral_index.txt', spec_index, fmt='%.6f')
