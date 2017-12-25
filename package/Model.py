from __future__ import print_function
from Utils import *

import numpy as np
import random
import pandas as pd
import inspect
import time
import sys
import os

class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)


#------------------------
#SPATIAL (Spherical-)GRID
#------------------------
 
def grid(XYZmax, NP, artist = False):

    """
    XYZmax: Maximum physical domain of the grid
    Please provide XYZmax as a 3D-list. XYZmax=[Xmax,Ymax,Zmax]
    NP: Number of grid points for each coordinate
    Please provide NP as a 3D-list. NP=[Nx,Ny,Nz]
    """

    print ('-------------------------------------------------\n-------------------------------------------------')
    
    if artist: artist = 1
    XYZmax = np.array(XYZmax)

    #----------------
    #NUMBER OF POINTS
    #----------------
    #Each NP becomes odd if it's not. This is done to force the grid to contain the midplane too.
    NP_dum = [ NP[i] + 1 if NP[i]%2 == 0 else NP[i] for i in range(3) ] 
    print('Changing the number of grid points to force it to contain the midplane...' 
          if NP != NP_dum
          else '... ... ...')
    NP = NP_dum
    
    #--------
    #XYZ GRID
    #--------
    print ('Calculating Grid...')
    
    Xmax,Ymax,Zmax = XYZmax    
    step = 2. * XYZmax / NP
    #epsilon = [RSun / 1.] * 3
    
    #In case of artist: a dummy point is created at the end of each coordinate. Artist won't read them but they are necessary for it to work fine! 
    XYZgrid = [np.linspace(-XYZmax[i], XYZmax[i] + artist * step[i], NP[i] + artist) for i in range(3)]

    X, Y, Z = XYZgrid
    
    #--------------------------------------
    #Extended Lists of distance coordinates
    #--------------------------------------
    rRxyzList = np.array([ ((x**2 + y**2 + z**2)**0.5, (x**2 + y**2)**0.5, x,y,z) for x in X for y in Y for z in Z])
    
    rList = rRxyzList[:,0] ; RList = rRxyzList[:,1]; xList = rRxyzList[:,2]; yList = rRxyzList[:,3]; zList = rRxyzList[:,4]
    rList = np.where(rList < 1., sorted(set(rList))[1] / 2 , rList ) # If r == 0: use the second minimum value of r divided by 2
    
    #"set" eliminates duplicates and "sorted" sorts values upward 
    RList = np.where(RList < 1., sorted(set(RList))[1] / 2 , RList )
    
    #-----
    #THETA
    #-----
    thetaList = np.arccos(zList / rList)
    halfPi = 0.5 * np.pi
    theta4vel = thetaList

    #If theta is under the xy plane, take its image above the plane.
    # This is done to avoid negative cosines in func streamline(), 
    #  taking advantage on the symmetry of discs and envelopes along xy.         
    deltaTh = np.where(thetaList > np.pi/2, thetaList - halfPi, 0.)
    thetaList = thetaList - 2 * deltaTh

    #-----
    #PHI
    #-----
    phiList = np.arctan2(yList,xList)
    twoPi =  2 * np.pi
    phiList = np.where(phiList < 0, phiList + twoPi, phiList )
    #-----
    #-----

    XYZ = [xList,yList,zList]
    rRTP = [rList,RList,thetaList,phiList]


    print ('Number of grid nodes for x,y,z:', NP)
    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

        
    return Struct( **{'XYZgrid': XYZgrid, 'XYZ': XYZ, 'rRTP': rRTP, 'theta4vel': theta4vel, 'NPoints': len(rList), 'Nodes': NP})


"""
#Not tested

#How the grid_FIXED function should work:                               
# 1: Provide the mapping function f(x)                            
# 2: Provide the cleared function x(f)                            
# 3: Provide the Maximum range of the coordinate                  
# 4: Provide the Number of points to make the map                 
# 5: (optional) Provide the coordinate datafile:                  
#    In this case must be done MU.grid_FIXED(0,0,0,0,"x.dat")     
#    Remember that the code is built to doesn't take the last     
#    value of the grid                                            
                                                                  
                                                                  
x=MU.grid_FIXED(lambda x: x*x, lambda x: x**0.5,8.*1900,15)       
y=MU.grid_FIXED(lambda x: x, lambda x: x,8.*1900,15)              
z=MU.grid_FIXED(lambda x: x, lambda x: x,8.*1900,15)              
                                                                  
                                                                  
#x=MU.grid_FIXED(0,0,0,0,"x.dat")                                 
#y=MU.grid_FIXED(0,0,0,0,"y.dat")                                 
#z=MU.grid_FIXED(0,0,0,0,"z.dat")                                 
                                                                  
                                                                  
GRID=[x,y,z]                                                      

def grid_FIXED(Fx=False,xF=False,Fmax=False,NP=False,File=False):

#Fx->F(x) -> Mapping function
#xF->x(F) -> Inverse function
#Fmax -> Real maximum in the domain

    
    xList=[]
    fList=[]
    xmax=0

    if Fx:
        
        xmax = xF(Fmax) #Finding the virtual maximum of x, 
                        # such that the domain remains between
                        #  [-Fmax,Fmax] after calculating F(x)
        
        
        xList=np.arange(-xmax,xmax+2*xmax/NP+2./NP,2*xmax/NP)
        fList=Fx(xList)
    
    else: fList=np.loadtxt(File)/AU
    '''
        if (fList[-1]-Fmax)<1e-6*Fmax:
            fList[0]=Fmax
            fList[-1]=Fmax
    '''
        
        
    return fList
"""

#---------------------------------
#Spherical base --> Cartesian base
#---------------------------------

def sphe_cart(V_sphe, theta, phi):
    
    #--------------------------------------------------------------
    #Wikipedia: Coordenadas esfericas, subseccion 'Base coordenada'
    #--------------------------------------------------------------
   
    if hasattr(theta, '__iter__') and hasattr(phi, '__iter__'): #A list of vectors
    #np.multiply and operator * are esentially the same even in terms of time
        vx = np.multiply( V_sphe, np.array([np.sin(theta) * np.cos(phi), np.cos(theta) * np.cos(phi), -np.sin(phi)]).T ).sum(1) 
        vy = np.multiply( V_sphe, np.array([np.sin(theta) * np.sin(phi), np.cos(theta) * np.sin(phi), np.cos(phi)]).T ).sum(1) 
        vz = np.multiply( V_sphe, np.array([np.cos(theta) , -np.sin(theta), np.zeros(len(theta))]).T ).sum(1)
     
    else: #Just one vector
     
        vx = np.dot( V_sphe, np.array([np.sin(theta) * np.cos(phi), np.cos(theta) * np.cos(phi), -np.sin(phi)]) )
        vy = np.dot( V_sphe, np.array([np.sin(theta) * np.sin(phi), np.cos(theta) * np.sin(phi), np.cos(phi)]) )
        vz = np.dot( V_sphe, np.array([np.cos(theta) , -np.sin(theta), 0.]))

    V_cart=np.array([vx,vy,vz])
    print ('%s is done!'%inspect.stack()[0][3])

    return V_cart

#-------------------------------------
#STREAM LINES FOR DENSITY AND VELOCITY
#-------------------------------------

def streamline(Rd, GRID):

#Cartesian-grid to work in. [xList,yList,zList]

#Values for cos(theta0), see Mendoza 2004

    rRTP = GRID.rRTP
    #------------
    #LISTS TO USE
    #------------
    rList, RList, thetaList = rRTP[:-1] #Phi is not used for streamlines
    rNorm = rList / Rd
    
    #----------------------------
    #MODEL. Mendoza et al. (2004)
    #----------------------------
    costheta = np.cos(thetaList)

    #Case r == Rd     
    costheta0 = np.where( rNorm == 1., costheta**(1./3), 0. )

    #Case r > Rd 
    term1 = abs(rNorm - 1) / 3.  
    term1sqrt = 2 * term1**0.5
    term1_15 = 2 * term1**1.5 

    good_ind = np.where( rNorm > 1. )
    term2 = np.zeros(len(term1))
    term2[good_ind] = np.sinh( 1./3 * np.arcsinh( rNorm[good_ind] * costheta[good_ind] / term1_15[good_ind] ) )

    costheta0 = np.where( rNorm > 1., term1sqrt * term2, costheta0 )

    #Case r < Rd and cond > 0
    cond = ( rNorm * 0.5 * costheta )**2 - term1**3

    good_ind = np.where( (rNorm < 1.) & (cond > 0.) )
    term2[good_ind] = np.cosh( 1./3 * np.arccosh( rNorm[good_ind] * costheta[good_ind] / term1_15[good_ind] ) )

    costheta0 = np.where( (rNorm < 1.) & (cond > 0.) , term1sqrt * term2, costheta0 )

    #Case r < Rd and cond < 0
    good_ind = np.where( (rNorm < 1.) & (cond < 0.) )
    term2[good_ind] = np.cos( 1./3 * np.arccos( rNorm[good_ind] * costheta[good_ind] / term1_15[good_ind] ) )
    costheta0 = np.where( (rNorm < 1.) & (cond < 0.) , term1sqrt * term2, costheta0 )

    #Case costheta0 > 1.0 due to computational waste of the order of 1e-16 above 1.0
    costheta0 = np.where ( costheta0 > 1.0, 1.0, costheta0)

    print ('%s is done!'%inspect.stack()[0][3])

    return costheta0 

#-------------------------------------
#-------------------------------------

#-------------------
#DENSITY FUNCTION
#-------------------

def density_Ulrich(RStar, Rd, rhoE0, Arho, GRID, discFlag=True, envFlag=False, rdisc_max = False, renv_max = False, ang_cavity = False):

#RStar: Star radius
#Rd: Centrifugal radius
#rhoE0: density at Rd and theta=pi/2
#Arho: Factor between envelope and disk densities
#GRID: Cartesian-grid to work in. [xList,yList,zList] 

    XYZgrid, XYZ, rRTP = GRID.XYZgrid, GRID.XYZ, GRID.rRTP
    NPoints = GRID.NPoints
    #------------
    #LISTS TO USE
    #------------
    zList = XYZ[2]
    rList, RList, thetaList = rRTP[:-1] #Phi is not used for density

    #----------------------------------------------
    #MODEL. Ulrich (1976) - Keto,E & Zhang,Q (2010)
    #----------------------------------------------

    #------------
    #DISC PROFILE
    #------------
    if discFlag:
        print ('Calculating Keplerian thin-disc density...')
        
        if not rdisc_max: rdisc_max = Rd #Keto's value for the disc to stop
        rhoD0 = Arho * rhoE0 #Normalization factor based on the envelope
        H0 = 0.01 * RStar #Scaleheight at RStar
        H = H0 * (RList / RStar)**1.25  #Scaleheight
        rhoDISC = np.where( RList <= rdisc_max, rhoD0 * (Rd / RList)**2.25 * np.exp(-0.5 * zList**2 / H**2), 1.0)#1.0e9 / 2)
        rhoDISC = np.where( rhoDISC < 1.0, 1.0, rhoDISC)
    else: 
        print ('No Keplerian thin-disc was invoked!')
        rhoDISC = np.zeros(NPoints)
    
    #----------------
    #ENVELOPE PROFILE
    #----------------
    if ang_cavity: print ('Set cavity for density with aperture %.1f deg'%(ang_cavity * 180 / np.pi))
    
    if envFlag:
        print ('Calculating stream lines for Ulrich envelope...')
        if not renv_max: renv_max = np.max(XYZgrid[0]) #A sphere inscribed in the coordinate X.

        costheta = np.cos(thetaList)
        costheta0 = streamline(Rd,GRID)
        print ('Calculating Envelope density...')
        rhoENV = np.where( (rList <= renv_max) & (thetaList >= ang_cavity), ( (rhoE0 * (rList / Rd)**-1.5) *
                                                                               ( (1 + (costheta / costheta0))**-0.5 ) *
                                                                               (1 + (Rd / rList) * (3. * costheta0**2 - 1))**-1 ), 1.0e9 )
        rhoENV = np.where( rhoENV < 1.0, 1.0, rhoENV)
    else:
        print ('No Envelope was invoked!')
        costheta0 = False
        rhoENV = np.zeros(NPoints)
    
    #----------------------------------------------
    #----------------------------------------------
    
    RHO = rhoDISC + rhoENV 

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'total': RHO, 'disc': rhoDISC, 'env': rhoENV, 'discFlag': discFlag, 'envFlag': envFlag, 
                      'r_disc': rdisc_max, 'r_env': renv_max, 'streamline': costheta0} )


#-------------------
#-------------------

#------------------------------
#DENSITY (Hamburguers) FUNCTION
#------------------------------

def density_Hamburgers(RStar, shFactor, Rd, rhoE0, Arho, GRID, discFlag=True, rdisc_max = False):

#RStar: Star radius
#shFactor: H0 = shFactor * RStar
#Rd: Centrifugal radius
#rhoE0: density at Rd and theta=pi/2
#Arho: Factor between envelope and disk densities
#GRID

    XYZ, rRTP = GRID.XYZ, GRID.rRTP

    #------------
    #LISTS TO USE
    #------------
    zList = XYZ[2]
    rList, RList = rRTP[:-2] 

    #----------------------------------------------
    #MODEL. Galvan-Madrid, A. Izquierdo (2017)
    #----------------------------------------------

    #------------
    #DISC PROFILE
    #------------
    if discFlag:
        print ('Calculating Burger-disc density...')
        if not rdisc_max: rdisc_max = Rd
        rhoD0 = Arho * rhoE0 
        H0 = shFactor * RStar
        H = H0 * (RList / RStar)**1.25  #Scaleheight
        rhoDISC = np.where( RList <= rdisc_max, rhoD0 * (Rd / RList)**2.25 * np.exp(-0.5 * zList**2 / H**2), 1.0)
        rhoDISC = np.where( rhoDISC < 1.0, 1.0, rhoDISC)
    else: 
        print ('No Disc was invoked!')
        rhoDISC = np.zeros(GRID.NPoints)
        
    #----------------------------------------------
    #----------------------------------------------

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'total': rhoDISC, 'disc': rhoDISC, 'env': 0., 'discFlag': True, 'envFlag': False, 'r_disc': rdisc_max, 'r_env': False} ) 

#------------------------------
#------------------------------

#---------------------------
#DENSITY (PowerLaw) FUNCTION
#---------------------------

def density_Powerlaw(r_max, rho_mean, q, GRID):

#rho_mean: Mean density of the Envelope 
#r_max: Maximum radius of the envelope 
#q: power-law for density
#GRID: Cartesian-grid to work in. [xList,yList,zList]     
#Env: On/Off the Envelope

    #------------
    #LISTS TO USE
    #------------
    rList, NPoints = GRID.rRTP[0], GRID.NPoints #Due to spherical symmetry only r is needed

    #------------------------
    #MODEL. Envelope powerlaw
    #------------------------
    print ('Calculating Envelope density with power-law...')
    rqList = np.where(rList <= r_max , rList**q, 0.)

    #As rho_mean = 1/NTotal * np.sum(rho0 * r**q), the normalization rho0 is calculated as follows:  
    rho0 = NPoints * rho_mean / np.sum(rqList)
    rhoENV = np.where( rho0 * rqList < 1.0, 1.0e9, rho0 * rqList )

    #------------------------
    #------------------------

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'total': rhoENV, 'disc': np.zeros(NPoints), 'env': rhoENV, 'discFlag': False, 'envFlag': True, 'r_disc': False, 'r_env': r_max} ) 

#---------------------------
#---------------------------

#---------------------------
#DENSITY (Constant) FUNCTION
#---------------------------

def density_Constant(Rd, GRID, discDens = 0, rdisc_max = False, envDens = 0, renv_max = False):

    rRTP = GRID.rRTP
    NPoints = GRID.NPoints
    #------------
    #LISTS TO USE
    #------------
    rList, RList = rRTP[:-2] 

    #------------
    #DISC PROFILE
    #------------
    if discDens:
        print ('Setting constant Disc density...')
        if not rdisc_max: rdisc_max = Rd
        rhoDISC = np.where( RList <= rdisc_max, discDens, 1.0)
    else: 
        print ('No Disc was invoked!')
        rhoDISC = np.zeros(NPoints)

    #----------------
    #ENVELOPE PROFILE
    #----------------
    if envDens:
        print ('Setting constant Envelope density...')
        if not renv_max: renv_max = np.max(GRID.XYZgrid[0])
        rhoENV = np.where( rList <= renv_max, envDens, 1.0)
    else: 
        print ('No Envelope was invoked!')
        rhoENV = np.zeros(NPoints)
    
    #----------------------------------------------
    #----------------------------------------------
    RHO = rhoDISC + rhoENV
        
    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'total': RHO, 'disc': rhoDISC, 'env': rhoENV, 'discFlag': bool(discDens), 'envFlag': bool(envDens), 'r_disc': rdisc_max, 'r_env': renv_max} ) 

#---------------------------
#---------------------------

#------------------
#ABUNDANCE FUNCTION
#------------------

def abundance(val, NPoints):

    abundList = np.ones(NPoints) * val
    
    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return abundList 

#------------------
#------------------

#-----------------
#GAS-TO-DUST RATIO
#-----------------

def gastodust(val, NPoints):
    
    gtdratioList = np.ones(NPoints) * val
    
    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return gtdratioList 

#-----------------
#-----------------

#----------------------
#TEMPERATURE FUNCTION
#----------------------

def temperature(TStar, Rd, T10Env, RStar, MStar, MRate, BT, p, density, GRID, ang_cavity = False):

#TStar: Star temperature
#T10Env: Envelope temperature at 10AU
#RStar: Star radius
#MStar: Star mass
#MRate: Mass accretion rate
#BT: Disc temperature factor
#p: Temperature power law exponent 
#GRID: [xList,yList,zList]

    rRTP = GRID.rRTP
    #------------
    #LISTS TO USE
    #------------
    rList, RList, thetaList = rRTP[:-1] 
    rhoDISC, rhoENV = density.disc, density.env

    #------------------------------
    #MODEL. Keto,E & Zhang,Q (2010)
    #------------------------------

    #------------
    #DISC Profile
    #------------
    if density.discFlag:
        print ('Calculating Keplerian thin-disc temperature...')
        rdisc = density.r_disc
        if density.envFlag:
            renv = density.r_env
            tempDISC = np.where( (RList <= rdisc) & (rList <= renv) , BT * (3*G * MStar * MRate / (4*np.pi * sigma * RList**3) * (1 - (RStar / RList)**0.5))**0.25, 30.0)
        else: tempDISC = np.where( RList <= rdisc , BT * (3*G * MStar * MRate / (4*np.pi * sigma * RList**3) * (1 - (RStar / RList)**0.5))**0.25, 30.0)
    else: tempDISC = 1.

    #----------------
    #ENVELOPE PROFILE
    #----------------
    if ang_cavity: print ('Set cavity for temperature with aperture %.1f deg'%(ang_cavity * 180 / np.pi))

    if density.envFlag:
        print ('Calculating Envelope temperature...')
        renv = density.r_env
        tempENV = np.where( (rList <= renv) & (thetaList >= ang_cavity), T10Env * 10**0.33 * (rList / AU)**-0.33, 30.0)
        #tempENV = TStar * (RStar / (2.*rList))**(2. / (4+p))
    else: tempENV = 1.
    
    #----------------------------------------------
    #----------------------------------------------

    #Weighted temperature with density 
    TEMP = (tempDISC * rhoDISC + tempENV * rhoENV) / density.total

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'total': TEMP, 'disc': tempDISC, 'env': tempENV, 'discFlag': density.discFlag, 'envFlag': density.envFlag} )

#----------------------
#----------------------

#---------------------------------
#TEMPERATURE (Hamburgers) FUNCTION
#---------------------------------

def temperature_Hamburgers(TStar, RStar, MStar, MRate, Rd, T10Env, shFactor, T_min, BT, p, density, GRID, inverted = False):

#TStar: Star temperature
#T10Env: Envelope temperature at 10AU
#RStar: Star radius
#MStar: Star mass
#MRate: Mass accretion rate
#BT: Disc temperature factor
#p: Temperature power law exponent 
#GRID: [xList,yList,zList]

    XYZ, rRTP = GRID.XYZ, GRID.rRTP
    #------------
    #LISTS TO USE
    #------------
    rList, RList = rRTP[:-2] 
    rhoDISC, rhoENV = density.disc, density.env

    #-----------------------------------------
    #MODEL. Galvan-Madrid, A. Izquierdo (2017)
    #-----------------------------------------

    #------------
    #DISC Profile
    #------------
    if density.discFlag:
        print ('Calculating Burger-disc temperature...')
        zList = XYZ[2]
        rdisc = density.r_disc
        H0 = shFactor * RStar 
        H = H0 * (RList / RStar)**1.25  #Whitney et al. (2003)
        
        if inverted: 
            print ('Set inverted temperature for Burger-disc...')
            tempDISC = np.where( RList <= rdisc, ( BT * (3*G * MStar * MRate / (4*np.pi * sigma * RList**3) * (1 - (RStar / RList)**0.5))**0.25
                                                   * np.exp(- 0.5 * zList**2 / H**2) ), 1.0) #Maximum in z = 0
        else: 
            print ('Set not inverted temperature for Burger-disc...')
            tempDISC = np.where( RList <= rdisc, ( BT * (3*G * MStar * MRate / (4*np.pi * sigma * RList**3) * (1 - (RStar / RList)**0.5))**0.25 
                                                   * np.exp(- 0.5  * (abs(zList) - H)**2 / H**2) ), 1.0) #Maximum in z = H

        tempDISC = np.where( (RList <= rdisc) & (tempDISC <= T_min), T_min, tempDISC)

    else: tempDISC = 1.

    #----------------
    #ENVELOPE PROFILE
    #----------------
    if density.envFlag:
        print ('Calculating Envelope temperature...')
        renv = density.r_env
        tempENV = np.where( rList <= renv, T10Env * 10**0.33 * (rList / AU)**-0.33, 1.0)
        #tempENV = TStar * (RStar / (2.*rList))**(2. / (4+p))
    else: tempENV = 1.
    
    #----------------------------------------------
    #----------------------------------------------

    #Weighted temperature with density 
    TEMP = (tempDISC * rhoDISC + tempENV * rhoENV) / density.total

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'total': TEMP, 'disc': tempDISC, 'env': tempENV, 'discFlag': density.discFlag, 'envFlag': density.envFlag} )

#---------------------------------
#---------------------------------

#-------------------------------
#TEMPERATURE (Constant) FUNCTION
#-------------------------------

def temperature_Constant(density, GRID, discTemp = 0, envTemp = 0):

    rRTP = GRID.rRTP
    NPoints = GRID.NPoints
    #------------
    #LISTS TO USE
    #------------
    rList, RList = rRTP[:-2]
    rhoDISC, rhoENV = density.disc, density.env

    #------------
    #DISC PROFILE
    #------------
    if discTemp:
        if density.discFlag:
            print ('Setting constant Disc temperature...')
            tempDISC = np.where( RList <= density.r_disc, discTemp, 30.0)
        else: sys.exit('ERROR: The disc calculation was turned ON but there is no density distribution for disc!')
    else: tempDISC = 0.

    #----------------
    #ENVELOPE PROFILE
    #----------------
    if envTemp:
        if density.envFlag:
            print ('Setting constant Envelope density...')
            tempENV = np.where( rList <= density.r_env, envTemp, 30.0)
        else: sys.exit('ERROR: The envelope calculation was turned ON but there is no density distribution for envelope!')
    else: tempENV = 0.
        
    #----------------------------------------------
    #----------------------------------------------
    
    #Weighted temperature with density 
    TEMP = (tempDISC * rhoDISC + tempENV * rhoENV) / density.total
        
    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'total': TEMP, 'disc': tempDISC, 'env': tempENV, 'discFlag': bool(discTemp), 'envFlag': bool(envTemp)} )

#-------------------------------
#-------------------------------

#----------------------
#VELOCITY FUNCTION
#----------------------

def velocity(RStar,MStar,Rd,density,GRID):

#MStar: Star mass
#Rd: Centrifugal radius
#GRID: [xList,yList,zList]

    XYZgrid, XYZ, rRTP = GRID.XYZgrid, GRID.XYZ, GRID.rRTP
    NPoints = GRID.NPoints
    #------------
    #LISTS TO USE
    #------------
    rList, RList, thetaList, phiList = rRTP
    rhoDISC, rhoENV = density.disc, density.env
    theta4vel = GRID.theta4vel

    #------------------------------
    #MODEL. Keto,E & Zhang,Q (2010)
    #------------------------------

    #------------
    #DISC Profile
    #------------
    if density.discFlag:
        print ('Calculating Disc velocity...')
        rdisc = density.r_disc
        #Pure azimuthal component. It's assumed that the radial velocity in the rotationally supported disc is comparatively small (Keto 2010).
        vdisc = np.where( RList <= rdisc, (G * MStar / RList)**0.5, 0.)  
    else: vdisc = 0.

    #----------------
    #ENVELOPE PROFILE
    #----------------
    if density.envFlag:
        print ('Calculating Envelope velocity...')
    
        #-------------------------
        #Useful THETA Calculations
        #-------------------------
        costheta0 = density.streamline
        costheta, sintheta, theta0 = np.cos(thetaList), np.sin(thetaList), np.arccos(costheta0)
        sintheta0 = np.sin(theta0)

        #-------------------------
        #-------------------------
        renv = density.r_env
        vr = np.where( rList <= renv, - (G * MStar / rList)**0.5 * (1 + costheta / costheta0)**0.5, 0.)
        
        signo = np.where( theta4vel> np.pi/2, -1, 1) #To respect symmetry. (Using the thetaList from 0 to pi)
        good_ind = np.where( (thetaList != 0.) & (rList <= renv) ) #To avoid the polar angle for theta and phi. In the polar angle all the velocity is radial
        vtheta, vphi = np.zeros(NPoints), np.zeros(NPoints)

        vtheta[good_ind] = signo[good_ind] * ( (G * MStar / rList[good_ind])**0.5 * (costheta0[good_ind] - costheta[good_ind]) / sintheta[good_ind] * 
                                               (1 + costheta[good_ind] / costheta0[good_ind])**0.5 )
        vphi[good_ind] = (G * MStar / rList[good_ind])**0.5 * (sintheta0[good_ind] / sintheta[good_ind]) * (1 - costheta[good_ind] / costheta0[good_ind])**0.5


    else: vr, vtheta, vphi = 0.,0.,0.
        
    #Weighted velocity with density. Vectorial sum.
    vr , vtheta, vphi = map( lambda x,y : (rhoDISC * x + rhoENV * y) / density.total, 
                             [ 0, 0, vdisc ], [vr, vtheta, vphi] ) 

    print ('Converting to cartesian coordinates...') 
    vx, vy, vz = sphe_cart( list( zip(vr, vtheta, vphi) ), theta4vel, phiList)
         
    #----------------------------------------------
    #----------------------------------------------
    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'x': vx, 'y': vy, 'z': vz} )


#----------------------
#----------------------

#--------------------------
#VELOCITY (Random) FUNCTION
#--------------------------

def velocity_random(v_disp,NPoints):

    print ('Calculating random (uniform) velocities...')
    v_disp = v_disp/np.sqrt(3)    

    v_x = v_disp * (2 * np.random.random(NPoints) - 1)  
    v_y = v_disp * (2 * np.random.random(NPoints) - 1)  
    v_z = v_disp * (2 * np.random.random(NPoints) - 1)  

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'x': v_x, 'y': v_y, 'z': v_z} )

#--------------------------
#--------------------------

#--------------------------------
#CENTRAL HOLE FUNCTION (Optional)
#--------------------------------

def MakeHole(T_min,T_max,dens_val,temp_val,abund_val,densList,tempList,abundList):

#The hole is between T_min and T_max
#T_min: Minimum temperature of the hole 
#T_max: Maximum temperature of the hole
#dens_val: Density value for the hole
#temp_val: Temperature value for the hole
#densList: Density list to modify 
#tempList: Temperature list to build the hole
    
    print ('Calculating profiles for a hole...')
    
    densNew = np.where( (tempList >= T_min) & (tempList <= T_max), dens_val, densList)
    tempNew = np.where( (tempList >= T_min) & (tempList <= T_max), temp_val, tempList)
    abundNew = np.where( (tempList >= T_min) & (tempList <= T_max), abund_val, abundList)

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')
    
    return Struct( **{'dens': densNew, 'temp': tempNew, 'abund': abundNew})


#--------------------------------
#--------------------------------    

def PrintProperties(density, temperature, GRID): 

    dx = GRID.XYZgrid[0][1] - GRID.XYZgrid[0][0]
    inddisc = np.where(temperature.disc > 2.)
    indtotal = np.where(temperature.total > 2.)
    Mu_MSun = 2 * Mu/MSun
    
    print ('Total mass (MSun):', np.sum(density.total) * dx**3 * Mu_MSun)
    print ('Mean Total Temperature (Kelvin), weighted by density:', 
           (np.sum(temperature.total[ indtotal ] * density.total[ indtotal ]) 
            / np.sum(density.total[ indtotal ]) ))
    if density.discFlag:
        print ('Total Disc mass (MSun):', np.sum(density.disc) * dx**3 * Mu_MSun)
        print ('Total Envelope mass (MSun):', np.sum(density.env) * dx**3 * Mu_MSun)
        print ('Mean Disc Temperature (Kelvin), weighted by density:', 
               (np.sum(temperature.disc[ inddisc ] * density.disc[ inddisc ])
                / np.sum(density.disc[ inddisc ]) ))
    
    print ('-------------------------------------------------\n-------------------------------------------------')

    
#---------------
#GEOMETRY STUFF 
#---------------   
                           
def Rotation_Matrix(angle_dicts):
    
    Rot = []

    for item in angle_dicts:

        #Rotation along X
        if item['axis'] == 'x':
            
            thX = item['angle'] 
            Rot.append( np.array([[1,0,0],
                                  [0,np.cos(thX),-np.sin(thX)],
                                  [0,np.sin(thX),np.cos(thX)]]) )
            print ("Added: rotation matrix along '%s' axis, %.1f deg"%('X',thX * 180/np.pi) )

        #Rotation along Y
        if item['axis'] == 'y':
            thY = item['angle'] 
            Rot.append( np.array([[np.cos(thY),0,-np.sin(thY)],
                                  [0,1,0],
                                  [np.sin(thY),0,np.cos(thY)]]) )
            print ("Added: rotation matrix along '%s' axis, %.1f deg"%('Y',thY * 180/np.pi) )
    
        #Rotation along Z
        if item['axis'] == 'z':
            thZ = item['angle'] 
            Rot.append( np.array([[np.cos(thZ),-np.sin(thZ),0],
                                  [np.sin(thZ),np.cos(thZ),0],
                                  [0,0,1]]) )
            print ("Added: rotation matrix along '%s' axis, %.1f deg"%('Z',thZ * 180/np.pi) )

    tmp = Rot[0]
    Rot_iter = iter(Rot[1:]) #Iterator for Rot_list from 2nd value: (matriz for matriz in Rot_list[1:])

    for i in range( len(Rot[1:]) ): 
        tmp = np.dot( next(Rot_iter) , tmp )
        
    Rot_total = tmp

    print ('%s is done!'%inspect.stack()[0][3])
   
    return Rot_total



def ChangeGeometry(GRID, center = False ,rot_dict = False, vel = False, vsys = False):

#order: indicates the order of the axes to make the rotation  
 
    rotPos, rotVel, modCenter, modVsys = [False]*4

    POS_vec = np.array(GRID.XYZ).T
    if vel == False and vsys: 
        sys.exit('ERROR: A velocity distribution is needed for the systemic velocity to be added!')
    if vel:
        VEL_vec = np.array([vel.x, vel.y, vel.z]).T

    #---------------
    #ROTATION Matrix
    #---------------

    if rot_dict:
        print ('Calculating Rotation matrix...')
        rot_angles , axis_order = rot_dict['angles'], rot_dict['axis']

        if len(rot_angles) == len(axis_order):
            angle_dicts = ( {'axis': ax, 'angle': ang} for ax,ang in zip(axis_order, rot_angles) )
        else:
            sys.exit('ERROR: rot_angles and axis_order lists must have the same size!') 
            
        NPoints = GRID.NPoints
        Rot_total = Rotation_Matrix(angle_dicts)     
        print ('Rotating Position vectors...')
        XYZ_it = iter(POS_vec)
        POS_vec = np.array([ np.dot( Rot_total, next(XYZ_it) ) for i in range(NPoints) ])
        rotPos = True
        
        if vel:
            print ('Rotating Velocity vectors...')
            VEL_it = iter(VEL_vec)
            VEL_vec = np.array([ np.dot( Rot_total, next(VEL_it) ) for i in range(NPoints) ])
            rotVel = True
        else: 
            print ('#################################################') 
            print ('WARNING: The VELOCITY distribution to be rotated was not provided! If you are not interested in the velocity structure ignore this warning...')
            print ('#################################################')

    if center is not False:
        print ('Moving the object to the new center...')
        center = np.array(center)
        POS_vec = POS_vec + center
        modCenter = True

    if vsys: 
        print ('Adding a systemic velocity along z...')
        vsys_vec = np.array([0,0,vsys])
        VEL_vec = VEL_vec + vsys_vec 
        modVsys = True

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    if (rotPos or modCenter) and (rotVel or modVsys): return Struct( **{ 'newXYZ': POS_vec.T , 'newVEL': VEL_vec.T  } )
    elif rotPos or modCenter: return Struct( **{ 'newXYZ': POS_vec.T } )
    elif rotVel or modVsys: return Struct( **{ 'newVEL': VEL_vec.T })
    else: return 0

#---------------
#---------------

#--------------
#WRITING DATA
#--------------

def DataTab_LIME(dens,temp,vel,abund,gtd,GRID, is_submodel = False, tag = False):
    
    import pandas

    if is_submodel:
        os.system('mkdir Subgrids')
        file0 = './Subgrids/datatab%s.dat'%tag
        file = open(file0,'w')
        x,y,z = GRID.XYZ
        print ('Writing Submodel data on %s'%file0)
        tmp = []
        for i in range(GRID.NPoints): 
            #file.write("%d %e %e %e %e %e %e %e %e %e %e\n"%
             #          (i,x[i],y[i],z[i],dens[i],temp[i],vel['x'][i],vel['y'][i],vel['z'][i],abund[i],gtd[i]))
            tmp.append( "%d %e %e %e %e %e %e %e %e %e %e\n"% (i,x[i],y[i],z[i],dens[i],temp[i],vel.x[i],vel.y[i],vel.z[i],abund[i],gtd[i]))
        file.writelines(tmp)
        
    else:
        files=['datatab.dat','x.dat','y.dat','z.dat']
        sizefile='./npoints.dat'
        print ('Writing grid size on %s'%sizefile)
        sfile = open(sizefile,'w') 
        Ns = GRID.Nodes
        sfile.write("%d %d %d %d"%(Ns[0],Ns[1],Ns[2],GRID.NPoints))
        print ('Writing data on %s'%files[0])
        file = open(files[0],'w')

        for i in range(GRID.NPoints): 
            file.write("%d %e %e %e %e %e %e %e\n"%
                       (i,dens[i],temp[i],vel.x[i],vel.y[i],vel.z[i],abund[i],gtd[i]))

        df = [pandas.DataFrame(GRID.XYZgrid[i]) for i in range(3)]
    
        for i in range(1,4):
            print ('Writing data on %s'%files[i])
            df[i-1].to_csv(files[i],index=False,header=False,float_format='%e') 
        
        sfile.close()
        
    file.close()
    
    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

#--------------
#--------------
