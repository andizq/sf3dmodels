**********
sf3dmodels
**********
  
This is the documentation for the Star-Forming regions 3D Modelling package.

Find the source code on `GitHub`_.

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

..
 .. image:: ../images/burger_face.png
    :width: 112.5

 .. image:: ../images/burger_edge.png
    :width: 112.5

 .. image:: ../images/burger_tapering.png
    :width: 112.5

 .. image:: ../images/global_objects.png
    :width: 112.5

 .. image:: ../images/compact_sources.png
    :width: 112.5

 .. image:: ../images/pIntensitypoints.png
    :width: 112.5

 
The output model can then be read with `LIME <https://lime.readthedocs.io/en/latest/>`_,
`RADMC-3D <http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/>`_ or `Polaris <http://www1.astrophysik.uni-kiel.de/~polaris/downloads.html>`_
to compute the radiative transfer of the modelled region, which is a key step prior comparison to real observations. 

.. image:: https://raw.githubusercontent.com/andizq/star-forming-regions/master/images/intro_sf3dmodels0_png.png 
   :width: 680

The sf3models *grid* and *rt* (radiative transfer) modules can also be used as wrappers between hydrodynamical simulations and 
radiative transfer codes. These modules are especially dedicated to treat irregular meshes (e.g. Voronoi meshes or SPH grid particles).

.. image:: https://raw.githubusercontent.com/andizq/star-forming-regions/master/images/intro_sf3dmodels1_png.png 
   :width: 680

Take a look at the following examples linking `AREPO <https://arepo-code.org>`_ and `Phantom <https://phantomsph.readthedocs.io>`_ hydrodynamical simulations with Polaris and LIME,

| `AREPO examples <https://github.com/andizq/star-forming-regions/tree/master/examples/arepo_code>`_
| `Phantom examples <https://github.com/andizq/star-forming-regions/tree/master/examples/phantom_code>`_

Requirements
============

* Python 2.7.x or 3.5.x (or later)
* `Astropy <http://docs.astropy.org/en/stable/install.html>`__
* `Numpy <https://www.scipy.org/install.html>`_
* `Matplotlib <https://matplotlib.org/users/installing.html>`_
* `IPython <https://ipython.org/install.html>`_ (recommended)

Installation
============

From source
-----------

Clone the star-forming-regions repository from `GitHub`_:

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

Update the package
------------------

.. code-block:: bash
   
   cd star-forming-regions
   git pull
   python setup.py install

Uninstall
---------

.. code-block:: bash
   
   pip uninstall sf3dmodels

Using sf3dmodels
================

.. toctree::
  :maxdepth: 1
  
  examples.rst
  single_source/single_source.rst	
  two_sources/two_sources.rst
  ionized_sources/ionized_sources.rst


LIME in *sf3dmodels* mode
=========================

The star-forming-regions repository includes a customised version of LIME v1.9.5, with additional routines for **sf3dmodels** users. 
This customisation allows using the **sf3dmodels** output as input for LIME in a user-friendly way. 

Installation
------------

LIME must be configured separately after installing the **sf3dmodels** package with success.

Have a look at the `LIME installation notes <https://github.com/andizq/lime/tree/sf3dmodels>`_ section. 
If you need to install third-party libraries (e.g. qhull, cfitsio, gsl), make sure to configure their paths to point to the LIME directory included in this repository. 

**We strongly recommend the user to install the 2010 version of Qhull, which can be done as follows:**

.. code-block:: bash

   wget https://github.com/qhull/qhull/archive/2010.1.tar.gz #download qhull v2010.1
   #curl -O https://github.com/qhull/qhull/archive/2010.1.tar.gz #you can also use curl
   gunzip 2010.1.tar.gz
   tar -xvf 2010.1.tar
   cd qhull-2010.1
   sh config/bootstrap.sh
   ./configure --prefix=/Users/andizq/star-forming-regions/lime
   emacs Makefile #if defined, delete flag -Wno-sign-conversion
   make
   make install

Running LIME
------------

LIME will look for the ``.dat`` files (with your model data) generated by **sf3dmodels** and load them into dedicated **sf3d** structures.
To invoke this option, a ``-S`` flag (capital S) must be added to the usual execution command as follows:

.. code-block:: bash

   lime -S model.c

The standard LIME command line options may also be invoked. For instance, to set LIME to 
(1.) produce normal output rather than the default ncurses output style, (2.)  
read the **sf3dmodels** output and (3.) run in parallel mode with 4 threads 
you should execute:

.. code-block:: bash

   lime -nS -p 4 model.c


*Note* that if the ``-S`` option *is not set* you will get back the 'default' operation of LIME.

LIME Examples
^^^^^^^^^^^^^

Take a look at the folders **lime/example** and **lime/example_sf3dmodels/** included in this repository.


Reference/API
=============

.. automodapi:: sf3dmodels
   :no-inheritance-diagram:

.. toctree::
   :maxdepth: 1
   
   Model.rst
   Plot_model.rst
   utils/utils.rst

New on sf3dmodels
=================

.. toctree::
   :maxdepth: 1   
   
   model/model.rst
   filament/filament.rst
   outflow/outflow.rst
   radiativetransfer/rt.rst
   grid/overlap.rst
   
Working on...
=============

.. toctree::
   :maxdepth: 1   
   
   tools/tools.rst
   arepo/arepo.rst
   grid/grid.rst

Developers
==========

* `Andres Izquierdo <https://github.com/andizq>`_
* `Roberto Galvan-Madrid <https://github.com/rgalvanmadrid>`_
* `Adam Ginsburg <https://github.com/keflavich>`_
* `Luke Maud <https://local.strw.leidenuniv.nl/people/touchscreen2/persinline.php?id=1716>`_   

Special thanks to those who have reported bugs or whose ideas and discussions have helped improve sf3dmodels, 

- Rowan Smith
- Yuxin Lin
- Antonio Hernandez
- Jonathan Henshaw
- Qizhou Zhang
- Leonardo Testi
- Stefano Facchini
- Ewine van Dishoeck
- Pietro Curone
- Carlos Carrasco-González
- Adriana Rodríguez-Kamenetzky

Papers using sf3dmodels
=======================

- `Izquierdo et al. (2018) <https://ui.adsabs.harvard.edu/abs/2018MNRAS.478.2505I/abstract>`_
- `Galvan-Madrid et al. (2018) <https://ui.adsabs.harvard.edu/abs/2018ApJ...868...39G/abstract>`_
- `Soler et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020arXiv200707285S/abstract>`_
- `Izquierdo et al. (2021a) <https://ui.adsabs.harvard.edu/abs/2021MNRAS.500.5268I/abstract>`_
- `Izquierdo et al. (2021b), <https://ui.adsabs.harvard.edu/abs/2021arXiv210409596I/abstract>`_
- `Carrasco-González et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021arXiv210601235C/abstract>`_
- Lin et al. (subm.)
- Curone et al. (subm.)

License
=======

This project is Copyright (c) Andres Izquierdo and licensed under
the terms of the BSD 3-Clause license. This package is based upon
the `Astropy package template <https://github.com/astropy/package-template>`_
which is licensed under the BSD 3-clause licence. See the licenses folder for
more information.

Citing sf3dmodels
=================

If you find **sf3dmodels** useful for your research please cite the work of `Izquierdo+2018`_::

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

   
.. warning:: This package, as well as its documentation, are currently under development and may often undergo modifications.
