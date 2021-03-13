
Cloud Complex description
-------------------------

* In this example we use a self-gravitating cloud complex influenced by a large-scale Milky Way-like galactic potential.
* Local gravitational forces are turned on for 2 Myr in a high resolution portion of the simulation, from which the cloud complex is extracted. 
* There are also supernova explosions randomly distributed all over the Galaxy.
* There is also chemical evolution of CO and Hydrogen species; sink particles representing star-forming regions; radiative heating and cooling and galactic differential rotation.

Find full details of the simulation setup of this and other types of cloud complexes in the Cloud Factory series of papers (`I <https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.1594S/abstract>`_, `II <https://ui.adsabs.harvard.edu/abs/2021MNRAS.500.5268I/abstract>`_)

These simulations are powered by a customised version of the AREPO code (`Springel, 2010 <https://ui.adsabs.harvard.edu/abs/2010MNRAS.401..791S/abstract>`_)

pcafactory-data repository
--------------------------

The data required for this example can be downloaded `here <https://girder.hub.yt/#user/5da06b5868085e00016c2dee/folder/5da06ef668085e00016c2df3>`_,
where you will find the following files:
 
* Simulation snapshot of the cloud complex.
* 12CO J=1-0 intensity cubes: 3 line intensity cubes for different cloud orientations (face-on, edge-on phi=0, edge-on phi=90) generated with sf3dmodels and LIME.
* Optical depth cube for an edge-on view of the complex.

Quick Tutorial
--------------

1. Download the simulation snapshot 
   
.. code-block:: bash

   curl https://girder.hub.yt/api/v1/item/5da0777768085e00016c2e01/download -o Restest_extracted_001_240

2. Read and clean the snapshot, and save the formatted physical information of the cloud for radiative transfer with LIME.

.. code-block:: bash
      
   python make_arepo_lime.py

.. image:: https://github.com/andizq/andizq.github.io/blob/master/pcafactory-data/examples-data/cldB_cloudfactory/cellsize_numdens-AREPOgrid.png?raw=true
   :width: 30%

.. image:: https://github.com/andizq/andizq.github.io/blob/master/pcafactory-data/examples-data/cldB_cloudfactory/3Dpoints_snap.png?raw=true
   :width: 30%

3. The output files are stored in the folder ./Subgrids by default

.. code-block:: bash
   
   cd Subgrids

4. Download the CO excitation information from the LAMDA database. 

.. code-block:: bash
   
   curl https://home.strw.leidenuniv.nl/~moldata/datafiles/co.dat -o co.dat 

5. Run the sf3dmodels version of LIME. The flag -S means that the grid was generated with `sf3dmodels <https://github.com/andizq/star-forming-regions>`_, and the flag -G indicates that the input grid is not uniform. The flag -n is to show log messages in terminal. We use 8 cores by setting -p 8 (LIME uses openmp for parallel processing). 

.. code-block:: bash

   lime -nSG -p 8 rt-lime.c 

The resulting synthetic cubes (.fits) can be found in the repository prepared for this example.
