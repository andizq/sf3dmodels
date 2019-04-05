# Recombination lines emission from a powerlaw-density HII sphere

Type the following commands:

```bash
$ python Exec_recomblines.py
```

## Spectrum of the whole region 
```bash
$ radmc3d spectrum lambdarange 59805.09 59873.11 nlam 150
$ python plot_spectrum.py
```

## Image cube with 50 velocity channels
```bash
$ radmc3d image lambdarange 59805.09 59873.11 nlam 50
$ python plot_fits_image.py
```

## External requirements (both need to be requested to the author): 

* The modified version of Radmc3d by [Peters+2012](http://adsabs.harvard.edu/abs/2012MNRAS.425.2352P) with the recombination lines extension.
* The computed table of fractional departure coefficients in your working directory (bn_cube_gs.dat).

## Model reference:

* Appendix B, model 1 from [Peters+2012](http://adsabs.harvard.edu/abs/2012MNRAS.425.2352P)