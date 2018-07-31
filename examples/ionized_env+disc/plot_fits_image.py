import radmc3dPy.image as image
import matplotlib.pyplot as plt
import numpy as np

i = 7
dist = 5410.
#lambdas = np.array([ 37500., 20000., 13043.47826087, 9375., 7142.85714286, 
#                     5454.54545455, 4285.71428571, 3529.41176471])
#en GHz: np.array([8,15,23,32,42,55,70,85])*1. GHz.


im = image.readImage()

#image.plotImage(im, au=True, log=True, maxlog=10, saturate=1e-5, cmap=plt.cm.gist_heat)
"""
image.plotImage(im, dpc=140, arcsec=True, log=True, maxlog=10, saturate=1e-5, cmap=plt.cm.gist_heat)
image.plotImage(im, dpc=dist, arcsec=True, log=True, maxlog=10, saturate=1e-5, cmap=plt.cm.gist_heat)
cim = im.imConv(fwhm=[1.0, 0.6], pa=120., dpc=dist)
image.plotImage(cim, arcsec=True, dpc=dist, log=True, maxlog=10, bunit='snu', cmap=plt.cm.gist_heat)
"""
im.writeFits('img_HII_env+disc.fits', dpc=5410.)#, coord='03h10m05s -10d05m30s')


#Mirar lo de inputs desde terminal
