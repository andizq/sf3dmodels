from astropy.io import fits
import numpy as np

data, header = fits.getdata('img_HII_env+disc-CONV.fits', header = True)
result = data[0] + np.random.normal(scale = 6.0e-6, size = np.shape(data[0])) #sigma = 6 microJy/beam
fits.writeto('img_HII_env+disc-CONV_noise.fits', result, header, overwrite=True)


