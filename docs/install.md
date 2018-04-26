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


Examples
--------

In this section I will show some examples to illustrate the main features of the package. Detailed information about modules, functions and parameters of specific models can be found in the help page of the package. For example, to see the help of the module `Model` and its function `density_Env_Disc`, type in your Python or IPython Command-Line the following commands:

```python
>>> import Model
>>> help(Model)
>>> help(Model.density_Env_Disc)
```

## Modelling a single star forming region

**Example 1.** Creating a massive star forming region (with *Ulrich envelope + Pringle disc*):
```python
#------------------
#Import the package
#------------------
from sf3dmodels import *
#-----------------
#Extra libraries
#-----------------
import numpy as np
import os
import time
```
**a.** Define general parameters:
```python
MStar = 7.0 * U.MSun
MRate = 4e-4 * U.MSun_yr #Mass accretion rate                                                                                                         
RStar = 26 * U.RSun * ( MStar/U.MSun )**0.27 * ( MRate / (1e-3*U.MSun_yr) )**0.41                                                                                                               
LStar = 3.2e4 * U.LSun
TStar = U.TSun * ( (LStar/U.LSun) / (RStar/U.RSun)**2 )**0.25                                                                                       
Rd = 152. * U.AU #Centrifugal radius  
```

**b.** Create the grid that will host the region:
```python
# Cubic grid, each edge ranges [-500, 500] AU.

sizex = sizey = sizez = 500 * U.AU
Nx = Ny = Nz = 150 #Number of divisions for each axis
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz])
NPoints = GRID.NPoints #Number of nodes in the grid
```

**c.** Invoke the physical properties from a desired model(s):
```python
#--------
#DENSITY
#--------
Rho0 = Res.Rho0(MRate, Rd, MStar) #Base density for Ulrich model
Arho = 24.1 #Disc-envelope density factor
Renv = 500 * U.AU #Envelope radius
Cavity = 40 * np.pi/180 #Cavity opening angle
density = Model.density_Env_Disc(RStar, Rd, Rho0, Arho, GRID, discFlag = True, envFlag = True,
                                 renv_max = Renv, ang_cavity = Cavity)

#-----------
#TEMPERATURE
#-----------
p = 0
T10Env = 375. #Envelope temperature at 10 AU                                                                                                              
BT = 5. #Adjustable factor for disc temperature. Extra, or less, disc heating.
temperature = Model.temperature(TStar, Rd, T10Env, RStar, MStar, MRate, BT, p, density, GRID,
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

**d.** Write the data into a file with the LIME format:
```python
#-----------------------------
#WRITING DATA with LIME format
#-----------------------------
Model.DataTab_LIME(density.total, temperature.total, vel, abundance, gtdratio, GRID)
```

**e.** Plot the results:
```python
#------------------------------------
#3D PLOTTING (weighting with density)
#------------------------------------
tag = 'Main'
weight = 10*Rho0
r = GRID.rRTP[0] / U.AU #GRID.rRTP hosts [r, R, Theta, Phi] --> Polar GRID
Plot_model.scatter3D(GRID, density.total, weight, NRand = 4000, colordim = r, axisunit = U.AU,
	             palette = 'jet', colorscale = 'log', colorlabel = r'${\rm log}_{10}(r [au])$',
		     output = 'Points%s.png'%tag, show = True)

#-------------------------------------
#2D Plotting (density and temperature) - UNDER DEVELOPMENT
#-------------------------------------

#FACE-ON and EDGE-ON profiles:

#Density: colormap
#Temperature: contours

#Plot_model.profile2D(GRID.XYZ, density.total, contours = temperature.total, unit=U.AU,
#                     palette='jet', output = 'density_profiles.png', tag = 'Main', show = True)
```

The resulting 3D distribution and 2D profiles:

<p align="left">
  <img src="/Users/andrespipecar42/SFRegions/SF3dmodels/docs/_build/html/_images/totalPointsMain.png" width="450"/>
  <img src="/Users/andrespipecar42/SFRegions/SF3dmodels/docs/_build/html/_images/Density_Temp_Main.png" width="200"/>
</p>

Edge-on and Face-on 3D distribution:

<p align="left">
  <img src="/Users/andrespipecar42/SFRegions/SF3dmodels/docs/_build/html/_images/totalPointsMain_a.png" width="325"/>
  <img src="/Users/andrespipecar42/SFRegions/SF3dmodels/docs/_build/html/_images/totalPointsMain_b.png" width="325"/>
</p>

<br>

**Example 2.** Creating a low-mass star forming region with a composite model for the density (*Envelope: Ulrich density + Disc: Burger density*) and a different model for the temperature.

**a.** The general parameters:
```python
MStar = 0.86 * U.MSun
MRate = 5.e-6 * U.MSun_yr
RStar = U.RSun * ( MStar/U.MSun )**0.8
LStar = U.LSun * ( MStar/U.MSun )**4
TStar = U.TSun * ( (LStar/U.LSun) / (RStar/U.RSun)**2 )**0.25
Rd = 264. * U.AU
```

**b.** The grid:
```python
# Cubic grid, each edge ranges [-500, 500] AU.

sizex = sizey = sizez = 500 * U.AU
Nx = Ny = Nz = 200 #Number of divisions for each axis
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz])
NPoints = GRID.NPoints #Number of nodes in the grid
```

**c.** The physical properties. Note how the final density Structure should be defined joining both, the envelope density and the disc density that were calculated separately from 2 different models:
```python
#-------------
#DENSITY
#-------------

#--------
#ENVELOPE
#--------
Rho0 = Res.Rho0(MRate, Rd, MStar)
Arho = None
Renv = 2.5 * Rd
densEnv = Model.density_Env_Disc(RStar, Rd, Rho0, Arho, GRID, discFlag = False, envFlag = True,
                                 renv_max = Renv)
#-------
#DISC
#-------
H0sf = 0.03 #Disc scale height factor (H0 = 0.03 * RStar)
Arho = 5.25
Rdisc = 1.5 * Rd
densDisc = Model.density_Hamburgers(RStar, H0sf, Rd, Rho0, Arho, GRID, discFlag = True, envFlag = False,
                                    rdisc_max = Rdisc)
#---------------------
#The COMPOSITE DENSITY
#---------------------
density = Model.Struct( **{ 'total': densEnv.total + densDisc.total,
                            'disc': densDisc.total,
                            'env': densEnv.total,
                            'discFlag': True,
                            'envFlag': True,
                            'r_disc': densDisc.r_disc,
                            'r_env': densEnv.r_env,
                            'streamline': densEnv.streamline} )

#-----------
#TEMPERATURE
#-----------
T10Env = 250. #Envelope temperature at 10 AU
Tmin = 10. #Minimum possible temperature. Every node with T<Tmin will inherit Tmin.
BT = 60. #Adjustable factor for disc temperature. Extra, or less, disc heating.
temperature = Model.temperature_Hamburgers(TStar, RStar, MStar, MRate, Rd, T10Env, H0sf, Tmin,
                                           BT, None, density, GRID, inverted = False)

#--------
#VELOCITY
#--------
vel = Model.velocity(RStar, MStar, Rd, density, GRID)

#-------------------------------
#ABUNDANCE and GAS-to-DUST RATIO
#-------------------------------
ab0 = 5e-8 #CH3CN abundance vs H2
abundance = Model.abundance(ab0, NPoints)

gtd0 = 100. #Gas to dust ratio (H2 vs Dust)
gtdratio = Model.gastodust(gtd0, NPoints)
```

**d.** Write the data into a file:
```python
#-----------------------------
#WRITING DATA with LIME format
#-----------------------------
Model.DataTab_LIME(density.total, temperature.total, vel, abundance, gtdratio, GRID)
```

**e.** Plot the results:
```python
#----------------------------------------
#3D PLOTTING (weighting with temperature)
#----------------------------------------
tag = 'Burger'
weight = 10*T10Env
Plot_model.scatter3D(GRID, temperature.total, weight, NRand = 4000, colordim = density.total/1e6, axisunit = U.AU,
	             palette = 'hot', colorscale = 'log', colorlabel = r'${\rm log}_{10}(\rho [cm^{-3}])$',
	             output = 'Points%s.png'%tag, show = True)

#-------------------------------------
#2D Plotting (density and temperature) - UNDER DEVELOPMENT
#-------------------------------------

#FACE-ON and EDGE-ON profiles:

#Density: colormap
#Temperature: contours

#Plot_model.profile2D(GRID.XYZ, density.total, contours = temperature.total, unit=U.AU,
#                     palette='jet', output = 'density_profiles.png', tag = 'Burger', show = True)
```

The resulting 3D distribution and 2D profiles:

<p align="left">
  <img src="/Users/andrespipecar42/SFRegions/SF3dmodels/docs/_build/html/_images/totalPointsBurger.png" width="450"/>
  <img src="/Users/andrespipecar42/SFRegions/SF3dmodels/docs/_build/html/_images/Density_Temp_Burger.png" width="200"/>
</p>

Edge-on and Face-on 3D distribution:

<p align="center">
  <img src="/Users/andrespipecar42/SFRegions/SF3dmodels/docs/_build/html/_images/totalPointsBurger_a.png" width="325"/>
  <img src="/Users/andrespipecar42/SFRegions/SF3dmodels/docs/_build/html/_images/totalPointsBurger_b.png" width="325"/>
</p>

<br>

## Modelling multiple star forming regions

**Example 1.** I will use the last two examples to illustrate how to join them in a *global grid*. The spatial region that is shared by two or more *sub-models* will inherit physical properties by weighting them with the local density, as explained in section 3.2 of Izquierdo et al (2018).

**The execution codes for both star forming regions are identical until just before the "writing" section.**

As there is no longer a single model, geometric changes will probably be required in each sub-model to better reproduce real scenarios. Let's add a couple of lines in the latter codes to account for the centering, inclination and systemic velocity of each region.

In the first example:

```python
#-------------------------
#ROTATION, VSYS, CENTERING
#-------------------------
xc, yc, zc = [-250*U.AU, 350*U.AU, 300*U.AU]
CENTER = [xc, yc, zc] #Center of the region in the global grid
v_sys = 3320. #m/s
newProperties = Model.ChangeGeometry(GRID, center = CENTER, vsys = v_sys,  vel = vel,
	      	 	             rot_dict = { 'angles': [np.pi/4, 1.87*np.pi], 'axis': ['x','z'] })
```

The **`GRID`** and **`vel`** objects should inherit the modified properties that **`newProperties`** currently hosts:

```python
#At the minute, the Model library only modifies the XYZ lists.
 #It is enough information for LIME
GRID.XYZ = newProperties.newXYZ

#The velocity should inherit the new velocity distribution
 #(as we rotated the system and added a systemic velocity)
vel.x, vel.y, vel.z = newProperties.newVEL
```

Finally, the writing process. We have to specify that the model is actually a **sub-model**:

```python
#-----------------------------
#WRITING DATA with LIME format
#-----------------------------
tag = '_Main' #A tag to identify the final files from others
Model.DataTab_LIME(density.total, temperature.total, vel, abundance, gtdratio, GRID,
		   is_submodel = True, tag = tag)
```

**The second example region include the same additions, only differ in the first specific definitions**:

```python
#-------------------------
#ROTATION, VSYS, CENTERING
#-------------------------
xc, yc, zc = [350*U.AU, -150*U.AU, -200*U.AU]
CENTER = [xc, yc, zc] #Center of the region in the global grid
v_sys = -2000. #m/s
newProperties = Model.ChangeGeometry(GRID, center = CENTER, vsys = v_sys,  vel = vel,
	      	 	             rot_dict = { 'angles': [-np.pi/2, np.pi/4], 'axis': ['y','z'] })
```

PS: once a sub-model is defined, a new folder named "**Subgrids**" will be created in the current working directory. All the sub-model data files will be saved there automatically. This order must be respected so that other features of the package work well

<br>

### Overlapping sub-models
Now that we have the data of each sub-model separately, we should invoke a new library in order to properly overlap their physical properties in a single grid we call **global grid**.

You can overlap all the sub-models available in the "Subgrids" folder, or tell to the module explicitly the list of sub-models to overlap:

```python
#------------------
#Import the package
#------------------
from sf3dmodels import BuildGlobalGrid as BGG, Model, Plot_model as Pm, Utils as U 

#---------------
#DEFINE THE GRID
#---------------
sizex = sizey = sizez = 1000 * U.AU
Nx = Ny = Nz = 120
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz])

#---------------
#INVOKE BGG LIB
#---------------
global_prop = BGG.overlap(GRID, all = True)
```

The next block is equivalent to the latter:
```python
#------------------
#Import the package
#------------------
from sf3dmodels import BuildGlobalGrid as BGG, Model, Plot_model as Pm, Utils as U 

#---------------
#DEFINE THE GRID
#---------------
sizex = sizey = sizez = 1000 * U.AU
Nx = Ny = Nz = 120
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz])

#---------------
#INVOKE BGG LIB
#---------------
list_sub = ['datatab_Main.dat', 'datatab_Burger.dat']
global_prop = BGG.overlap(GRID, submodels = list_sub)
```

Plotting the result: the 3D points distribution follow the density field (see the `weight` parameter) in both plots. The colormaps represent the density in one plot and the temperature in the other.
```python
GRID = global_prop.GRID 
density = global_prop.density / 1e6 #1e6 to convert from m^-3 to cm^-3
temperature = global_prop.temperature

weight = 400 * np.mean(density)

#-----------------
#Plot for DENSITY
#-----------------
Pm.scatter3D(GRID, density, weight, NRand = 7000, axisunit = U.AU, colorscale = 'log', palette = 'hot',
  	     colorlabel = r'${\rm log}_{10}(\rho [cm^{-3}])$', output = 'global_grid_dens.png')

#--------------------
#Plot for TEMPERATURE
#--------------------
Pm.scatter3D(GRID, density, weight, colordim = temperature, NRand = 7000, axisunit = U.AU, colorscale = 'log',
             palette = 'brg', colorlabel = r'${\rm log}_{10}(T$ $[K])$', output = 'global_grid_temp.png')
```

***Left***: density colormap. ***Right***: temperature colormap.

<p align="center">
  <img src="/Users/andrespipecar42/SFRegions/SF3dmodels/docs/_build/html/_images/global_grid_dens.png" width="325"/>
  <img src="/Users/andrespipecar42/SFRegions/SF3dmodels/docs/_build/html/_images/global_grid_temp.png" width="325"/>
</p>

<br>

### Modelling HII regions + Radiative Transfer with RADMC-3D

You can create electronic density distributions with the library `Model` and use the module `Datatab_RADMC3D_FreeFree` to obtain the necessary formatted files to predict the Free-Free emission of the region with **RADMC-3D**.  

**Example 1.** Ionized spherical region with constant density and temperature.
```python
#------------------
#Import the package
#------------------
from sf3dmodels import *
#-----------------
#Extra libraries
#-----------------
import numpy as np
import time
```  
**a.** The general parameters and the GRID definition (the `radmc3d` flag should be turned on):

```python
#------------------
#General Parameters
#------------------
r_max = 2530 * U.AU #H II sphere size
dens_e = 1.4e5 * 1e6 #Electronic numerical density, from cgs to SI
t_e = 1.e4 #K

#---------------
#GRID Definition
#---------------
sizex = sizey = sizez = 2600 * U.AU 
Nx = Ny = Nz = 63 #Number of divisions for each axis
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz], radmc3d = True)
NPoints = GRID.NPoints #Final number of nodes in the grid
```

**b.** Invoke the library `Model` to assign the physical properties to each node in the `GRID`:

```python
#-------------------
#PHYSICAL PROPERTIES
#-------------------
density = Model.density_Constant(r_max, GRID, envDens = dens_e)
temperature = Model.temperature_Constant(density, GRID, envTemp = t_e, backTemp=2.725)

Model.PrintProperties(density, temperature, GRID) #Printing resultant properties (mass, mean temperature, etc)
```

**c.** Write the data into a file with the RADMC-3D format:
```python
#---------------------------------
#WRITING DATA with RADMC-3D FORMAT
#---------------------------------
Model.Datatab_RADMC3D_FreeFree(density.total, temperature.total, GRID)
```

**d.** Plot a random 3D distribution of points based on the physical properties of the model
```python
#------------------------------------
#3D PLOTTING (weighting with density)
#------------------------------------
tag = 'HII'
weight = dens_e
Plot_model.scatter3D(GRID, density.total, weight, NRand = 4000, colordim = density.total / 1e6, axisunit = U.AU, palette = 'jet', 
                     colorscale = 'log', colorlabel = r'$n_{\rm e}$ [cm$^{-3}$]', output = 'totalPoints%s.png'%tag, show = True)
```

<p align="center">
  <img src="/Users/andrespipecar42/SFRegions/SF3dmodels/docs/_build/html/_images/ctsphere_HII.png" width="325"/>
</p>

<br>

**e.** Now let's execute RADMC-3D. For a SED:

```console
radmc3d sed dpc 4000
```

And its plot:

```python
from radmc3dPy.analyze import *
import matplotlib.pyplot as plt

tag = 'ctsphere'

s = readSpectrum(fname = 'spectrum.out') #column 0: wavelength in microns; column 1: Flux in cgs. 
distance = 4000. #in pc. The spectrum.out file is still normalized to a distance of 1 pc (see radmc3d docs)
F_nu = s[:,1] * distance**-2 * 1e23 #to Jy at the set distance
nu = 3e8 * s[:,0]**-1 * 1e6 * 1e-9 #microns to GHz
plt.plot(nu, F_nu)
plt.title('%s - distance: %d pc'%(tag,distance))
plt.xlabel('Frequency [GHz]'); plt.ylabel('Flux [Jy]')
plt.xscale('log'); plt.yscale('log')
plt.savefig('sed_'+tag+'.png')
plt.show()
```

<p align="center">
  <img src="/Users/andrespipecar42/SFRegions/SF3dmodels/docs/_build/html/_images/sed_ctsphere.png" width="325"/>
</p>

<br>

**f.** Now let's calculate a 2D-Image at 300 GHz (1000 microns):

```console
radmc3d image lambda 1000
```

And its plot:

```python
from radmc3dPy.image import *
from matplotlib import cm
a=readImage()
plotImage(a,log=True,maxlog=4,cmap=cm.hot,bunit='snu',dpc=140,arcsec=True) #or au=True
```

<p align="center">
  <img src="/Users/andrespipecar42/SFRegions/SF3dmodels/docs/_build/html/_images/image_ctsphere.png" width="325"/>
</p>

<br>
  
**Example 2.** Ionized spherical region with a power-law density and constant temperature.

Here the only difference with the example 1 are the general parameters and the invocation of a different model for the density distribution.

```python
#------------------
#General Parameters
#------------------
#from Galvan-Madrid et al. 2009, Table 3:

MStar = 34 * U.MSun
r_max = 2530 * U.AU #H II sphere size
r_min = r_max / 200 #Minimum distance (!= 0 to avoid indeterminations)
r_s = r_max #Normalization distance
rho_s = 1.4e5 * 1e6 #from cgs to SI. Density at r_s
q = 1.3 #Density powerlaw  
t_e = 1.e4 #K

#---------------
#GRID Definition
#---------------

sizex = sizey = sizez = 2600 * U.AU 
Nx = Ny = Nz = 63 #Number of divisions for each axis
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz], radmc3d = True)
NPoints = GRID.NPoints #Final number of nodes in the grid

#-------------------
#PHYSICAL PROPERTIES
#-------------------
density = Model.density_Powerlaw_HII(r_min, r_max, r_s, rho_s, q, GRID)
temperature = Model.temperature_Constant(density, GRID, envTemp = t_e, backTemp=2.725)

Model.PrintProperties(density, temperature, GRID) #Printing resultant properties (mass, mean temperature, etc)
```

The resultant plots:

<p align="center">
  <img src="/Users/andrespipecar42/SFRegions/SF3dmodels/docs/_build/html/_images/plsphere_HII.png" width="325"/>
  <img src="/Users/andrespipecar42/SFRegions/SF3dmodels/docs/_build/html/_images/sed_plsphere.png" width="325"/>
</p>

<br>

<p align="center">
  <img src="/Users/andrespipecar42/SFRegions/SF3dmodels/docs/_build/html/_images/image_plsphere.png" width="325"/>
</p>

<br>


**Example 3.** Keto+03 spherical region envolving a Pringle disc; constant temperature.

```python
#------------------
#General Parameters
#------------------
#from Galvan-Madrid et al. 2009, Table 3:

MStar = 34 * U.MSun
r_max = 2530 * U.AU #1000 * U.AU #H II sphere size
r_min = r_max / 200 #Minimum distance (!= 0 to avoid indeterminations)
rho_s = 1.5e6 * 1e6 #from cgs to SI. Density at sonic radius
q = 1.3 #Density powerlaw
t_e = 1.e4 #K

#-------------------------------
#Parameters for the Pringle disc
#-------------------------------
MRate = 3e-4 * U.MSun_yr
RStar = U.RSun * ( MStar/U.MSun )**0.8

#---------------
#GRID Definition
#---------------

sizex = sizey = sizez = 2600 * U.AU 
Nx = Ny = Nz = 63 #Number of divisions for each axis
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz], radmc3d = True)
NPoints = GRID.NPoints #Final number of nodes in the grid

#-------------------
#PHYSICAL PROPERTIES
#-------------------

#--------
#ENVELOPE
#--------
densEnv = Model.density_Keto_HII(MStar, r_min, r_max, rho_s, t_e, GRID, q = 1.5)

#-------
#DISC
#-------
Rd = 10*densEnv.rs #10 times the sonic radius, just to make it visible
Rho0 = Res.Rho0(MRate, Rd, MStar)
Arho = 60.0 #/ 500 
densDisc = Model.density_Env_Disc(RStar, Rd, Rho0, Arho, GRID, discFlag = True, envFlag = False, 
                                rdisc_max = Rd)

density = Model.Struct( **{ 'total': densEnv.total + densDisc.total,
                            'disc': densDisc.total, 
                            'env': densEnv.total,
                            'discFlag': True,
                            'envFlag': True,
                            'r_disc': densDisc.r_disc, 
                            'r_env': densEnv.r_env} )

temperature = Model.temperature_Constant(density, GRID, discTemp=t_e, envTemp=t_e, backTemp=2.725)

Model.PrintProperties(density, temperature, GRID)
```

The resultant plots:

<p align="center">
  <img src="/Users/andrespipecar42/SFRegions/SF3dmodels/docs/_build/html/_images/keto+disc_HII.png" width="325"/>
  <img src="/Users/andrespipecar42/SFRegions/SF3dmodels/docs/_build/html/_images/sed_keto+disc.png" width="325"/>
</p>

<br>

<p align="center">
  <img src="/Users/andrespipecar42/SFRegions/SF3dmodels/docs/_build/html/_images/image_keto+disc.png" width="325"/>
</p>

<br>



