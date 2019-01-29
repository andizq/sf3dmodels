Modelling HII regions. *RT with RADMC-3D*
=========================================

With the `~sf3dmodels.Model` module you can create electronic density distributions and use the `~sf3dmodels.Model.Radmc3d` and/or 
`~sf3dmodels.Model.Radmc3dRT` classes to obtain formatted data files that may be used later to predict the emission of your model with `RADMC-3D`_.

Example 1
---------

Source codes and figures on GitHub: `ionized_constant <https://github.com/andizq/star-forming-regions/tree/master/examples/ionized_constant>`_

.. note::
   `Uniform spherical HII region`
   
   `Model`: *Uniform density and constant temperature*
   
   `Radiative transfer`: free-free emission
   
   `Useful references`: `Pratap+1992`_, `Keto+2008`_


**The preamble**:  

.. code-block:: python

   #------------------
   #Import the package
   #------------------
   from sf3dmodels import *
   #-----------------
   #Extra libraries
   #-----------------
   import numpy as np
   import time


**a.** Define the general parameters and the GRID:

.. note:: The ``radmc3d`` flag must be turned on at the GRID definition (see below)

.. code-block:: python

   #------------------
   #General Parameters
   #------------------
   r_max = 2530 * U.AU #HII sphere size
   dens_e = 1.4e5 * 1e6 #Electronic numerical density, from cgs to SI
   t_e = 1.e4 #K

   #---------------
   #GRID Definition
   #---------------
   sizex = sizey = sizez = 2600 * U.AU 
   Nx = Ny = Nz = 63 #Number of divisions for each axis
   GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz], radmc3d = True)
   NPoints = GRID.NPoints #Final number of nodes in the grid


**b.** Invoke the `~sf3dmodels.Model` module to compute the physical properties at each ``GRID`` node:

.. code-block:: python

   #-------------------
   #PHYSICAL PROPERTIES
   #-------------------
   density = Model.density_Constant(r_max, GRID, envDens = dens_e)
   temperature = Model.temperature_Constant(density, GRID, envTemp = t_e, backTemp=2.725)


**c.** Write the data into a file with the RADMC-3D format:

.. code-block:: python

   #----------------------
   #WRITING RADMC-3D FILES
   #----------------------
   Rad = Model.Radmc3dRT(GRID)
   prop = {'dens_elect': density.total,
           'dens_ion': density.total,
	   'tgas': temperature.total}
   Rad.freefree(prop)


**d.** Plot the 3D spatial points distribution:

.. code-block:: python

   #------------------------------------
   #3D PLOTTING (weighting with density)
   #------------------------------------
   tag = 'HII'
   weight = dens_e
   Plot_model.scatter3D(GRID, density.total, weight, NRand = 4000, 
   			colordim = density.total / 1e6 / 1e5, axisunit = U.AU, 
			cmap = 'winter', marker = 'o', 
			colorlabel = r'$n_{\rm e}$ [cm$^{-3}$]', 
			output = '3Ddens_%s.png'%tag, show = True)

   Plot_model.scatter3D(GRID, density.total, weight, NRand = 4000, 
   			colordim = temperature.total, axisunit = U.AU, 
			cmap = 'winter', marker = 'o', 
			colorlabel = r'$T_{\rm e}$ [Kelvin]', 
			output = '3Dtemp_%s.png'%tag, show = True)


.. image:: https://github.com/andizq/andizq.github.io/blob/master/star-forming-regions/examples/ionized_constant/3Ddens_ctsphere_HII.png?raw=true
   :width: 49.5%

.. image:: https://github.com/andizq/andizq.github.io/blob/master/star-forming-regions/examples/ionized_constant/3Dtemp_ctsphere_HII.png?raw=true
   :width: 49.5%


Running RADMC-3D
^^^^^^^^^^^^^^^^

Making SEDs: In the folder where you stored the ``sf3dmodels`` output data files (**.inp**'s)
you should run the following command:

.. code-block:: bash

   $ radmc3d sed

This command writes the file ``spectrum.out`` in your working directory; 
it has two columns: 1. flux in cgs units and 2. wavelength in microns. 
Let's use that information to construct the Spectral Energy Distribution (SED) 
of the region, at a distance of 4 kpc:

.. code-block:: python

   from radmc3dPy.analyze import *
   import matplotlib.pyplot as plt

   tag = 'ctsphere'

   s = readSpectrum(fname = 'spectrum.out') #column 0: wavelength in microns; column 1: Flux in cgs. 
   distance = 4000. #in pc. The spectrum.out file is still normalized to a distance of 1 pc (see radmc3d docs)
   F_nu = s[:,1] * distance**-2 * 1e23 #to Jy at the set distance
   nu = 3e8 * s[:,0]**-1 * 1e6 * 1e-9 #microns to GHz
   plt.plot(nu, F_nu)
   plt.title('%s - distance: %d pc'%(tag,distance))
   plt.xlabel('Frequency [GHz]'); plt.ylabel('Flux Density [Jy]')
   plt.xscale('log'); plt.yscale('log')
   plt.savefig('sed_'+tag+'.png')
   plt.show()

.. image:: https://github.com/andizq/andizq.github.io/blob/master/star-forming-regions/examples/ionized_constant/sed_ctsphere.png?raw=true
   :width: 59.5%
   :align: center
   :alt: sed for constant-density sphere

Now let's have a look at the emission of this region at 300 GHz (or 1000 microns). The -simple- command for radmc3d would be:

.. code-block:: bash

   $ radmc3d image lambda 1000

which writes a file named ``image.out``. The following commands make a simple 2D plot from it:

.. code-block:: python
   
   from radmc3dPy.image import readImage, plotImage
   from matplotlib import cm
   a=readImage()
   plotImage(a,log=True,maxlog=4,cmap=cm.hot,bunit='snu',dpc=4000,arcsec=True) #or au=True

.. image:: https://github.com/andizq/andizq.github.io/blob/master/star-forming-regions/examples/ionized_constant/image_ctsphere.png?raw=true
   :width: 69.5%
   :align: center
   :alt: 2D emission for constant-density sphere


Example 2
---------

Source codes and figures on GitHub: `ionized_powerlaw <https://github.com/andizq/star-forming-regions/tree/master/examples/ionized_constant>`_

.. note::
   `Power-law spherical HII region`
   
   `Model`: *Power-law density distribution and constant temperature*

   `Radiative transfer`: free-free emission
   
   `Useful references`: `Keto+2008`_, `Galvan-Madrid+2009`_


The only difference with respect to the Example 1 is the density model (`~sf3dmodels.Model.density_Powerlaw_HII()`):

.. code-block:: python

   #------------------
   #General Parameters
   #------------------
   #from Galvan-Madrid et al. 2009, Table 3:

   MStar = 34 * U.MSun
   r_max = 2530 * U.AU #H II sphere size
   r_min = r_max / 200 #Minimum distance (!= 0 to avoid indeterminations).
   r_s = r_max #Normalization distance
   rho_s = 1.4e5 * 1e6 #from cgs to SI. Density at r_s
   q = 1.3 #Density powerlaw  
   t_e = 1.e4 #K

   #-------------------
   #PHYSICAL PROPERTIES
   #-------------------
   density = Model.density_Powerlaw_HII(r_min, r_max, r_s, rho_s, q, GRID)
   temperature = Model.temperature_Constant(density, GRID, envTemp = t_e, backTemp=2.725)

And the resulting plots:

.. image:: https://github.com/andizq/andizq.github.io/blob/master/star-forming-regions/examples/ionized_powerlaw/3Ddens_plsphere_HII.png?raw=true
   :width: 49.5%

.. image:: https://github.com/andizq/andizq.github.io/blob/master/star-forming-regions/examples/ionized_powerlaw/sed_plawsphere.png?raw=true
   :width: 45.5%
   :alt: sed for powerlaw-density HII sphere

.. image:: https://github.com/andizq/andizq.github.io/blob/master/star-forming-regions/examples/ionized_powerlaw/image_plawsphere.png?raw=true
   :width: 59.5%
   :align: center
   :alt: 2D emission for powerlaw-density HII sphere
