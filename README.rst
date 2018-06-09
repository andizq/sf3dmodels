`Documentation`_ | `View on Github`_ | `Download Stable ZIP`_

Star Forming regions 3D Modelling package
-----------------------------------------

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
    :alt: Powered by Astropy Badge

SF3dmodels is a star forming region(s) modelling package that brings together
analytical models in order to compute their physical properties in a 3D grid. The
package can couple different models in a single grid to recreate complex star
forming systems as those being revealed by current instruments. 
The output data can be read with `LIME <https://lime.readthedocs.io/en/latest/>`_ 
or `RADMC-3D <http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/>`_ 
to carry out radiative transfer calculations of the model.


Requirements
------------

* [Matplotlib](https://matplotlib.org/users/installing.html)
* [Numpy](https://www.scipy.org/install.html)
* [Pandas](http://pandas.pydata.org/pandas-docs/stable/install.html)
* [Astropy](http://docs.astropy.org/en/stable/install.html)
* [IPython](https://ipython.org/install.html) (optional, but recommended)

Installation
------------

Once you have installed the required packages please clone or download the star-forming-regions repository from [`GitHub https://github.com/andizq/star-forming-regions`](https://github.com/andizq/star-forming-regions) and go into the folder using the command line:

```bash
cd /path/to/local/repository/
```

Finally, run the following command:

```bash
python setup.py develop
```



License
-------

This project is Copyright (c) Andres Izquierdo and licensed under
the terms of the BSD 3-Clause license. This package is based upon
the `Astropy package template <https://github.com/astropy/package-template>`_
which is licensed under the BSD 3-clause licence. See the licenses folder for
more information.


.. _Download Stable ZIP: https://github.com/andizq/star-forming-regions/archive/master.zip
.. _View on Github: https://github.com/andizq/star-forming-regions/
.. _docs: http://star-forming-regions.readthedocs.io
.. _Documentation: http://star-forming-regions.readthedocs.io