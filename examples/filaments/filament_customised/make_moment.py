from astropy.io import fits
from astropy.wcs import WCS
from spectral_cube import SpectralCube
import matplotlib.pyplot as plt

fn = 'img_custom_CO_J1-0_LTE_jypxl_edgeon.fits'
hdu = fits.open(fn)[0]
hdu.header['CUNIT3'] = 'm/s' 
w = WCS(hdu.header)
cube = SpectralCube(data=hdu.data.squeeze(), wcs=w.dropaxis(3))

m0 = cube.moment(order=0)
m1 = cube.moment(order=1)
m0.write('moment0_'+fn, overwrite=True)
m1.write('moment1_'+fn, overwrite=True)

fig, ax = plt.subplots(ncols=2, subplot_kw = {'projection': w.celestial}, figsize = (12,6))

im0 = ax[0].imshow(m0.array, cmap='hot', origin='lower')
im1 = ax[1].imshow(m1.array, cmap='nipy_spectral', origin='lower')

cb0 = fig.colorbar(im0, ax=ax[0], orientation='horizontal', format='%.2f', pad=0.1)
cb1 = fig.colorbar(im1, ax=ax[1], orientation='horizontal', pad=0.1)

ax[0].set_title('Moment 0 - edgeon', fontsize = 14)
ax[1].set_title('Moment 1 - edgeon', fontsize = 14)
cb0.set_label(r'Jy pxl$^{-1}$ m s$^{-1}$', fontsize = 14)
cb1.set_label(r'm s$^{-1}$', fontsize = 14)

plt.tight_layout(pad=2.0)
fig.savefig('customfilament_moments_edgeon.png')
