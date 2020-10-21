.. _Download Stable ZIP: https://github.com/andizq/star-forming-regions/archive/master.zip
.. _Download: https://github.com/andizq/star-forming-regions/archive/master.zip
.. _View on Github: https://github.com/andizq/star-forming-regions/
.. _docs: http://star-forming-regions.readthedocs.io
.. _Full Documentation: http://star-forming-regions.readthedocs.io

`Full Documentation`_ | `View on Github`_ | `Download Stable ZIP`_

sf3dmodels
----------

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
    :alt: Powered by Astropy Badge

.. image:: https://img.shields.io/badge/ascl-2001.003-blue.svg?colorB=262255
    :target: http://ascl.net/2001.003
    :alt: ASCL Code Record

sf3dmodels is a 3D modelling package that collects analytical prescriptions of 
gas/dust **envelopes**, **discs**, **outflows** and **filaments** in order to reproduce complex star-forming 
systems such as those being revealed by state-of-the-art telescopes. The package can (i) model 
individual/isolated objects or, alternatively, (ii) couple multiple models in a single grid 
to recreate composite regions.
 
The output model can then be read with `LIME <https://lime.readthedocs.io/en/latest/>`_,
`RADMC-3D <http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/>`_ or `Polaris <http://www1.astrophysik.uni-kiel.de/~polaris/downloads.html>`_
to compute the radiative transfer of the modelled region, which is a key step prior comparison to real observations. 

The sf3models *grid* and *rt* (radiative transfer) modules can also be used as wrappers between hydrodynamical simulations and 
radiative transfer codes. These modules are especially dedicated to treat irregular meshes (e.g. Voronoi meshes or SPH grid particles).


Take a look at the following examples linking `AREPO <https://arepo-code.org>`_ and `Phantom <https://phantomsph.readthedocs.io>`_ hydrodynamical simulations with Polaris and LIME.

Documentation
-------------

Find the full documentation and tutorials on the sf3dmodels `website <http://star-forming-regions.readthedocs.io>`_.

Requirements
------------

* Python 2.7.x 3.5.x (or later)
* `Astropy <http://docs.astropy.org/en/stable/install.html>`__
* `Numpy <https://www.scipy.org/install.html>`_
* `Matplotlib <https://matplotlib.org/users/installing.html>`_
* `IPython <https://ipython.org/install.html>`_ (recommended)

Installation
------------

From source
***********

Clone the star-forming-regions repository from `GitHub <https://github.com/andizq/star-forming-regions>`_:

If you have a github account, type in a terminal:

.. code-block:: bash

   git clone git@github.com:andizq/star-forming-regions.git

if you don't have one:

.. code-block:: bash

   git clone https://github.com/andizq/star-forming-regions.git

Get into the star-forming-regions folder and run the ``setup.py`` script in installation mode:

.. code-block:: bash

   cd star-forming-regions
   python setup.py install

You can run any example from **star-forming-regions/examples** to check if the installation was succesful.

Upgrade the package
*******************
   
.. code-block:: bash
   
   cd star-forming-regions
   git fetch --all
   git reset --hard origin/master
   git submodule update --force --remote -- lime
   python setup.py install

Uninstall
*********

.. code-block:: bash
   
   pip uninstall sf3dmodels


Developers
----------

* `Andres Izquierdo <https://github.com/andizq>`_
* `Roberto Galvan-Madrid <https://github.com/rgalvanmadrid>`_
* `Adam Ginsburg <https://github.com/keflavich>`_
* `Luke Maud <https://local.strw.leidenuniv.nl/people/touchscreen2/persinline.php?id=1716>`_   

Special thanks to those who have reported bugs or whose ideas and discussions helped improve sf3dmodels, 

- Rowan Smith
- Yuxin Lin
- Antonio Hernandez
- Jonathan Henshaw
- Qizhou Zhang
- Leonardo Testi
- Stefano Facchini
- Ewine van Dishoeck
- Pietro Curone

Papers using sf3dmodels
-----------------------

- `Izquierdo et al. (2018) <https://ui.adsabs.harvard.edu/abs/2018MNRAS.478.2505I/abstract>`_
- `Galvan-Madrid et al. (2018) <https://ui.adsabs.harvard.edu/abs/2018ApJ...868...39G/abstract>`_
- `Soler et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020arXiv200707285S/abstract>`_
- Izquierdo et al. (2020, subm.) The Cloud Factory II

License
-------

This project is Copyright (c) Andres Izquierdo and licensed under
the terms of the BSD 3-Clause license. This package is based upon
the `Astropy package template <https://github.com/astropy/package-template>`_
which is licensed under the BSD 3-clause licence. See the licenses folder for
more information.

Citing sf3dmodels
-----------------

If you find **sf3dmodels** useful for your research please cite the work of `Izquierdo et al. (2018) <https://ui.adsabs.harvard.edu/abs/2018MNRAS.478.2505I/abstract>`_::

   @ARTICLE{2018MNRAS.478.2505I,
      author = {{Izquierdo}, Andr{\'e}s F. and {Galv{\'a}n-Madrid}, Roberto and
                {Maud}, Luke T. and {Hoare}, Melvin G. and {Johnston}, Katharine G. and
         	{Keto}, Eric R. and {Zhang}, Qizhou and {de Wit}, Willem-Jan},
      title = "{Radiative transfer modelling of W33A MM1: 3D structure and dynamics of a complex massive star-forming region}",
      journal = {\mnras},
      keywords = {radiative transfer, stars: formation, stars: massive, stars: protostars, Astrophysics - Astrophysics of Galaxies, Astrophysics - Solar and Stellar Astrophysics},
      year = "2018",
      month = "Aug",
      volume = {478},
      number = {2},
      pages = {2505-2525},
      doi = {10.1093/mnras/sty1096},
      archivePrefix = {arXiv},
      eprint = {1804.09204},
      primaryClass = {astro-ph.GA},
      adsurl = {https://ui.adsabs.harvard.edu/abs/2018MNRAS.478.2505I},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
   }