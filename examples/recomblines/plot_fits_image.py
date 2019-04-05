import radmc3dPy.image as image
import matplotlib.pyplot as plt
import numpy as np
import subprocess

output = 'img_rl_powerlaw.fits'
subprocess.call(['rm',output])

dist = 2000.
im = image.readImage()

data = im.image #Accessing image data
plt.imshow(data[:,:,29], origin='lower', cmap='cool') #plotting a single channel
plt.show()

im.writeFits(fname=output, dpc=dist, coord='03h10m05s -10d05m30s') #writting 3d fits cube
