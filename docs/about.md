About SF3dmodels
------------

SF3dmodels is a star-forming-region(s) modelling package that brings together analytical physical models from different authors (Ulrich, Keto, Pringle, Whitney) to compute density, velocity, temperature, abundance and gas-to-dust ratio 3D distributions.

SF3dmodels makes it possible to construct diverse distributions by mixing the available standard models within it. However, new models can be added to the package if needed. Users willing to contribute and nourish the package with new models are very welcome!

In addition, SF3dmodels can couple different star forming regions together to recreate complex star forming systems as those being revealed by recent telescopes and interferometers. This feature was quite used by Izquierdo et al, 2018 to model the complex star forming system W33A MM1.

We made the SF3dmodels output data to be compatible with [LIME](https://github.com/lime-rt/lime): To simulate real observations you need first to perform radiative transfer calculations and LIME does this for you in far-infrared and (sub-)millimeter wavelengths.


Requirements
--------------------

* [Matplotlib](https://matplotlib.org/users/installing.html)
* [Numpy](https://www.scipy.org/install.html)
* [Pandas](http://pandas.pydata.org/pandas-docs/stable/install.html)
* [Astropy](http://docs.astropy.org/en/stable/install.html)
* [IPython](https://ipython.org/install.html) (optional, but recommended)

Installation process
--------------------

Once you have installed the required packages please clone or download the star-forming-regions repository from [https://github.com/andizq/star-forming-regions](https://github.com/andizq/star-forming-regions) (unzip it if necessary) and go into the folder using the Terminal:

```bash
cd /path/to/local/repository/
```

Finally, run the following command:

```bash
python setup.py install
```
