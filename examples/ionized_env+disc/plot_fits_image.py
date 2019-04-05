import radmc3dPy.image as image
import matplotlib.pyplot as plt
import numpy as np

dist = 4000. #pc
im = image.readImage()

#image.plotImage(im, au=True, log=True, maxlog=10, saturate=1e-5, cmap=plt.cm.gist_heat)
"""
image.plotImage(im, dpc=140, arcsec=True, log=True, maxlog=10, saturate=1e-5, cmap=plt.cm.gist_heat)
image.plotImage(im, dpc=dist, arcsec=True, log=True, maxlog=10, saturate=1e-5, cmap=plt.cm.gist_heat)
cim = im.imConv(fwhm=[1.0, 0.6], pa=120., dpc=dist)
image.plotImage(cim, arcsec=True, dpc=dist, log=True, maxlog=10, bunit='snu', cmap=plt.cm.gist_heat)
"""
im.writeFits('img_HII_env+disc.fits', dpc=dist)#, coord='03h10m05s -10d05m30s')

