About SF3dmodels
------------

SF3dmodels is a star-forming-region(s) modelling package that brings together analytical physical models from different authors (Ulrich, Keto, Pringle, Whitney) to compute density, velocity, temperature, abundance and gas-to-dust ratio 3D distributions.

SF3dmodels makes it possible to construct diverse distributions by mixing the available standard models within it. However, new models can be added to the package if needed. Users willing to contribute and nourish the package with new models are very welcome!

In addition, SF3dmodels can couple different star forming regions together to recreate complex star forming systems as those being revealed by recent telescopes and interferometers. This feature was quite used by Izquierdo et al, 2018 to model the complex star forming system W33A MM1.

We made the SF3dmodels output data to be compatible with [LIME](https://github.com/lime-rt/lime): To simulate real observations you need first to perform radiative transfer calculations and LIME does this for you in far-infrared and (sub-)millimeter wavelengths.


Installation process
--------------------

Download the package folder and place it in your prefered location. Then, modify the Python Module Search Path in your shell startup file as follows:

#### for Bash shell (sh):

In your startup file (**`~/.bash_profile`** , **`~/.bash_login`** or **`~/.bashrc`** [Linux users], and **`~/.profile`** [Mac OS users]) include:

```bash
export PYTHONPATH=${PYTHONPATH}:/path/to/package
```

#### for C Shell (csh; tcsh):

In your startup file (**`~/.cshrc`** , **`~/.tcshrc`** or **`~/.login`**) include:

```csh
setenv PYTHONPATH $PYTHONPATH\:/path/to/package
```

### External required packages

* [Matplotlib](https://matplotlib.org/users/installing.html)
* [Numpy](https://www.scipy.org/install.html)
* [Pandas](http://pandas.pydata.org/pandas-docs/stable/install.html)
* [Astropy](http://docs.astropy.org/en/stable/install.html) (optional, but recommended)
* [IPython](https://ipython.org/install.html) (optional, but recommended)
