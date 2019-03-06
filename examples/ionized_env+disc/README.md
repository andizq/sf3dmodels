# Free-free emission and SED from HII envelope + disc

The radmc3d output image will be convolved by a given beam. Also, a gaussian noise will be added. 

Type the following commands:

```bash
$ python ExecModel_keto+disc.py
$ radmc3d sed
$ python plot_sed.py
$ radmc3d image lambda 9090 phi 90 incl 45 sizeau 8000 dpc 4000 npix 100
$ python plot_fits_image.py
$ python plot_cont.py
$ python convolucionGauss_All.py
$ python add_noise_fits.py
$ python plot_cont.py -i CONV_noise
```