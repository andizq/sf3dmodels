from radmc3dPy.analyze import *
import matplotlib.pyplot as plt

tag = 'ctsphere'

s = readSpectrum(fname = 'spectrum.out') #column 0: wavelength in microns; column 1: Flux in cgs.
distance = 4000. #in pc. The spectrum.out file is still normalized to a distance of 1 pc (see radmc3d docs)
F_nu = s[:,1] * distance**-2 * 1e23 #to Jy at the set distance
nu = 3e8 * s[:,0]**-1 * 1e6 * 1e-9 #microns to GHz
plt.plot(nu, F_nu)
plt.title('%s - distance: %d pc'%(tag,distance))
plt.xlabel('Frequency [GHz]'); plt.ylabel('Flux [Jy]')
plt.xscale('log'); plt.yscale('log')
plt.savefig('sed_'+tag+'.png')
plt.show()
