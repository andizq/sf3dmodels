# Recombination lines emission from a powerlaw-density HII sphere

Type the following commands:

```bash
$ python Exec_recomblines.py
$ radmc3d spectrum lambdarange 59805.09 59873.11 nlam 150
$ python plot_spectrum.py
```

## External requirements (both need to be requested to the author): 

* The modified version of Radmc3d by [Peters+2012](http://adsabs.harvard.edu/abs/2012MNRAS.425.2352P) with the recombination lines extension.
* The computed table of fractional departure coefficients.

## Model reference:

* Appendix of [Peters+2012](http://adsabs.harvard.edu/abs/2012MNRAS.425.2352P)