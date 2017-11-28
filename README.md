# star-forming-regions

About "Name" 
------------

"Name" is a star-forming-region(s) modelling package that brings together analytical physical models from different authors (Ulrich, Keto, Pringle, Whitney) to compute density, velocity, temperature, abundance and gas-to-dust ratio 3D distributions. 

"Name" makes it possible to construct diverse distributions by mixing the available standard models within it. However, new models can be added to the package if needed. Users willing to contribute and nourish the package with new models are very welcome! 

In addition, "Name" can couple different star-forming-regions together to recreate complex star forming systems as those being revealed by recent telescopes and interferometers. This feature was quite used by Izquierdo et al, 2018 to model the complex star forming system W33A MM1.

We made the "Name" output data to be compatible with [LIME](https://github.com/lime-rt/lime): To simulate real observations you need first to perform radiative transfer calculations and LIME does this for you in far-infrared and (sub-)millimeter wavelength.

Installation process
--------------------

Download the package folder and place it in your prefered location :open_file_folder:. Then, modify the Python Module Search Path in your shell startup file as follows:

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
* [Pandas](http://pandas.pydata.org/pandas-docs/stable/install.html) :panda_face:
* [Astropy](http://docs.astropy.org/en/stable/install.html) (optional, but recommended)
* [IPython](https://ipython.org/install.html) (optional, but recommended)


Examples
--------

In this section I will show some examples to illustrate the main features of the package. Detailed information about modules, functions and parameters of specific models can be found in the help page of the package. For example, to see the help of the module `Model` and its function `density_Ulrich`, type in your Python or IPython Command-Line the following commands: 

```python
>>> import Model
>>> help(Model)
>>> help(Model.density_Ulrich)
```

## Modelling a single star forming region :high_brightness:

**Example 1.** Creating a massive star forming region (with *Ulrich envelope + Pringle disc*):
```python
#-----------------
#Package libraries
#-----------------
import Model
import Resolution as Res
import Plot_model
import Utils as U #Module with useful units
#-----------------
#Extra libraries
#-----------------
import numpy as np
import os
import time
```
**a.** Define general parameters :page_with_curl::
```python
MStar = 7.0 * U.MSun 
MRate = 4e-4 * U.MSun_yr #Mass accretion rate                                                                                                         
RStar = 26 * U.RSun * ( MStar/U.MSun )**0.27 * ( MRate / (1e-3*U.MSun_yr) )**0.41                                                                                                               
LStar = 3.2e4 * U.LSun
TStar = U.TSun * ( (LStar/U.LSun) / (RStar/U.RSun)**2 )**0.25                                                                                       
Rd = 152. * U.AU #Centrifugal radius  
```

**b.** Create the grid that will host the region :house_with_garden::
```python
# Cubic grid, each edge ranges [-500, 500] AU. 

sizex = sizey = sizez = 500 * U.AU
Nx = Ny = Nz = 150 #Number of divisions for each axis
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz])
NPoints = GRID.NPoints #Number of nodes in the grid
```

**c.** Invoke the physical properties from a desired model(s) :pager::
```python
#--------
#DENSITY
#--------
Rho0 = Res.Rho0(MRate, Rd, MStar) #Base density for Ulrich model
Arho = 24.1 #Disc-envelope density factor 
Renv = 500 * U.AU #Envelope radius
Cavity = 40 * np.pi/180 #Cavity opening angle
density = Model.density_Ulrich(RStar, Rd, Rho0, Arho, GRID, discFlag = True, envFlag = True, 
                               renv_max = Renv, ang_cavity = Cavity)
                             
#-----------
#TEMPERATURE
#-----------
T10Env = 375. #Envelope temperature at 10 AU                                                                                                              
BT = 5. #Adjustable factor for disc temperature. Extra, or less, disc heating.
temperature = Model.temperature_Keto(TStar, Rd, T10Env, RStar, MStar, MRate, BT, density, GRID, 
                                     ang_cavity = Cavity)

#--------
#VELOCITY
#--------
vel = Model.velocity_Ulrich(RStar, MStar, Rd, density, GRID)

#-------------------------------
#ABUNDANCE and GAS-to-DUST RATIO
#-------------------------------
ab0 = 1.8e-7 #CH3CN abundance                                                                                                           
abundance = Model.abundance(ab0, NPoints) #Constant abundance

gtd0 = 100. #Gas to dust ratio
gtdratio = Model.gastodust(gtd0, NPoints) #Constant gtd ratio
```

**d.** Write the data into a file :memo::
```python
#-----------------------------
#WRITING DATA with LIME format
#-----------------------------
Model.DataTab_LIME(density.total, temperature.total, vel, abundance, gtdratio, GRID)
```

**e.** Plot the results :fireworks::
```python
#---------------------
#3D Plotting (density)
#---------------------
Plot_model.scatter3D(GRID.XYZ, density.total,  NRand = 2000,  unit=AU, palette='Blues', 
                     scale='log', output = 'density.png', show = True)

#-------------------------------------
#2D Plotting (density and temperature)
#-------------------------------------

#FACE-ON and EDGE-ON profiles:

#Density: colormap
#Temperature: contours

Plot_model.profile2D(GRID.XYZ, density.total, contours = temperature.total, unit=AU, 
                     palette='jet', output = 'density_faceon.png', tag = 'Main', show = True)
```

The resulting 3D distribution and 2D profiles: 


<p align="left">
  <img src="/images/totalPoints_Main.png" width="500"/>
  <img src="/images/Density_Temp_Main.png" width="250"/>
</p>

