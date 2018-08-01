from radmc3dPy.image import *
from matplotlib import cm
a=readImage()
plotImage(a,log=True,maxlog=2,cmap=cm.hot,bunit='snu',dpc=4000,au=True) #or arcsec=True
