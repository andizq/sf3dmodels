import numpy as np
import matplotlib.pyplot as plt

spec = np.loadtxt('spectrum.out', skiprows = 3)
plt.plot(spec[:,0], spec[:,1])
plt.savefig('img_spectrum.png')
#plt.show()
