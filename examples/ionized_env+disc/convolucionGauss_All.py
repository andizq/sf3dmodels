from __future__ import print_function
from astropy.convolution import Gaussian2DKernel, convolve, convolve_fft
from astropy.io import fits
import numpy as np

#FWHM/2 = Radio equivalente = (a*b)**0.5  -> radio que corresponde a la misma area
# de una elipse con semieje mayor a y semieje menor b.
#FWHM->Desviacion estandar: FWHM=2*raiz(2*log(2))*stddev = 2.35482*stddev

a_beam = 0.19 / 2 #in arcsecs
b_beam = 0.19 / 2


sizeau = 8000. #Size along each (x,y) direction
npix = 100. #Number of pixels along each (x,y) direction
dpc = 4000. #Distance in parsecs

Resolucion = (sizeau / npix) / dpc #in arcsecs 

a_pix = a_beam / Resolucion #in pxls
b_pix = b_beam / Resolucion

R_gauss = np.sqrt(a_pix * b_pix)
print ('Rgauss:',R_gauss)

HWHM=R_gauss
FWHM=2*HWHM
stddev=FWHM/(2.35482)
gaussian_2D_kernel = Gaussian2DKernel(stddev,stddev) #astropy requiere stddev

areaBeamPix=1.442*np.pi*HWHM**2
#2*np.log(2)=1.386 para Func.Gaussiana
#1.442 para Func.Bessel


#----------------
#READING PROCESS
#----------------

data_cont, header_cont = fits.getdata('img_HII_env+disc.fits', header = True)
result_cont = []
ch = 0

print ("Convolving continuum image...")
for image_channel in data_cont:
    result_cont.append(areaBeamPix*convolve(image_channel,gaussian_2D_kernel))
    print (ch)
    ch+=1

semiaxis_deg = 2*HWHM*Resolucion/3600    

header_cont['BUNIT'] = 'Jy/beam'
header_cont['BPA'] = 0
header_cont['BMIN'] = semiaxis_deg #Eje menor en degrees
header_cont['BMAJ'] = semiaxis_deg #Eje mayor en degrees

fits.writeto('img_HII_env+disc-CONV.fits',result_cont,header_cont,overwrite=True) 
