import numpy as np
import matplotlib.pyplot as plt

s = np.loadtxt('spectrum.out', skiprows = 3) #column 0: wavelength in microns; column 1: Flux density in cgs.
distance = 2000. #in pc. The spectrum.out file is normalized to a distance of 1 pc (see radmc3d docs)
F_nu = s[:,1] * distance**-2 * 1e23 #to Jy at the new distance
nu = 3e8 * s[:,0]**-1 * 1e6 * 1e-9 #microns to GHz

plt.plot(nu, F_nu)
plt.title('recomb.line - distance: %d pc'%distance)
plt.xlabel('Frequency [GHz]'); plt.ylabel('Flux density [Jy]')
plt.savefig('img_spectrum.png')
plt.show()
