import radmc3dPy.image as image
import subprocess

output = 'img_rl_powerlaw.fits'
subprocess.call(['rm',output])

dist = 2000.
im = image.readImage()

data = im.image #Image data, channels lie on axis 2. 

image.plotImage(im, dpc=dist, au=True, bunit='jy/pixel', ifreq=29, cmap='gnuplot') #plotting channel 29
im.writeFits(fname=output, dpc=dist, coord='03h10m05s -10d05m30s') #writting 3d fits cube
