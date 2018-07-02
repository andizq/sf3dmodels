from __future__ import print_function
from .Utils import *

import numpy as np
import random
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
 
def grid(XYZmax, NP, artist = False, radmc3d = False):

    """
    XYZmax: Maximum physical domain of the grid
    Please provide XYZmax as a 3D-list. XYZmax=[Xmax,Ymax,Zmax]
    NP: Number of grid points for each coordinate
    Please provide NP as a 3D-list. NP=[Nx,Ny,Nz]
    """

    print ('-------------------------------------------------\n-------------------------------------------------')
    
    XYZmax = np.array(XYZmax)

    #----------------
    #NUMBER OF POINTS
    #----------------
    #Each NP becomes odd if it's not. This is done to force the grid to contain the midplane too.
    NP_dum = [ NP[i] + 1 if NP[i]%2 == 0 else NP[i] for i in xrange(3) ] 
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
    
    #In case of radmc3d or artist: 
    # a dummy point is created at the end of each coordinate. Artist  won't read them but they are necessary for it to work well! 
    
    if radmc3d or artist: 
        step = 2. * XYZmax / (NP - np.ones(3))
        XYZgrid = [np.linspace(-XYZmax[i], XYZmax[i], NP[i] + 1) for i in xrange(3)]
        #XYZgrid = [np.append( np.linspace(-XYZmax[i], XYZmax[i], NP[i]), (XYZmax[i] + step[i]) ) for i in xrange(3)]
        X, Y ,Z = XYZgrid #The grid must contain an extra node but...
        X = 0.5 * ( X[0:NP[0]] + X[1:NP[0]+1] ) #Moving the node from the corner to the center of the boxel
        Y = 0.5 * ( Y[0:NP[1]] + Y[1:NP[1]+1] )  
        Z = 0.5 * ( Z[0:NP[2]] + Z[1:NP[2]+1] )
  
        #X = X[:-1]; Y = Y[:-1]; Z = Z[:-1] #...the calculations must be done w/o that node 
    else: #lime
        XYZgrid = [np.linspace(-XYZmax[i], XYZmax[i], NP[i]) for i in xrange(3)]
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

    dx = abs(X[1] - X[0]) #Linear step in X
    XYZ = [xList,yList,zList]
    rRTP = [rList,RList,thetaList,phiList]
    

    print ('Number of grid nodes for x,y,z:', NP)
    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

        
    return Struct( **{'XYZgrid': XYZgrid, 'XYZ': XYZ, 'rRTP': rRTP, 'theta4vel': theta4vel, 'NPoints': len(rList), 'Nodes': NP, 'step': dx})

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

def density_Env_Disc(RStar, Rd, rhoE0, Arho, GRID, discFlag=True, envFlag=False, rdisc_max = False, renv_max = False, ang_cavity = False):

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
        print ('Calculating Keplerian flared-disc density...')
        
        if not rdisc_max: rdisc_max = Rd #Keto's value for the disc to stop
        rhoD0 = Arho * rhoE0 #Normalization factor based on the envelope
        H0 = 0.01 * RStar #Scaleheight at RStar
        H = H0 * (RList / RStar)**1.25  #Scaleheight
        rhoDISC = np.where( RList <= rdisc_max, rhoD0 * (Rd / RList)**2.25 * np.exp(-0.5 * zList**2 / H**2), 1.0)#1.0e9 / 2)
        rhoDISC = np.where( rhoDISC < 1.0, 1.0, rhoDISC)
    else: 
        print ('No Keplerian flared-disc was invoked!')
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

def density_Hamburgers(RStar, shFactor, Ro, rhoE0, Arho, GRID, 
                       p = 2.25, q = 0.5, 
                       Rt = False, discFlag=True, rdisc_max = False):

#RStar: Star radius
#shFactor: Scaleheight normalization constant: H0 = shFactor * RStar
#Ro: Outer radius of the disk
#rhoE0: density at Rd and theta=pi/2
#Arho: Density factor 
#GRID

    XYZ, rRTP = GRID.XYZ, GRID.rRTP

    #------------
    #LISTS TO USE
    #------------
    zList = XYZ[2]
    rList, RList = rRTP[:-2] 

    #-----------------------------------------------------
    #MODEL. Chin-Fei Lee, Zhi-Yun Li, Paul Ho, et al. 2017
    #-----------------------------------------------------

    #------------
    #DISC PROFILE
    #------------
    if discFlag:
        print ('Calculating Burger-disc density...')
        if not rdisc_max: rdisc_max = Ro
        rhoD0 = Arho * rhoE0 
        H0 = shFactor * RStar
        print ('Scaleheight normalization constant:', H0 / AU * 1 / ((RStar/AU)**(1 + 0.5*(1-q))))
        H = H0 * (RList / RStar)**(1 + 0.5*(1-q)) #Scaleheight, with no tapering 
        if Rt: H = H * np.exp(-((RList-Rt) / (Ro-Rt))**2)  #Scaleheight, with tapering 
        rhoDISC = np.where( RList <= rdisc_max, rhoD0 * (RList / Ro)**-p * np.exp(-0.5 * zList**2 / H**2), 1.0)
        rhoDISC = np.where( rhoDISC < 1.0, 1.0, rhoDISC)
    else: 
        print ('No Disc was invoked!')
        rhoDISC = np.zeros(GRID.NPoints)
        
    #----------------------------------------------
    #----------------------------------------------

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'total': rhoDISC, 'disc': rhoDISC, 'env': 0., 'H': H,  
                      'discFlag': discFlag, 'envFlag': False, 'Rt': Rt, 'r_disc': rdisc_max, 'r_env': False} ) 

#------------------------------
#------------------------------

#---------------------------
#DENSITY (PowerLaw) FUNCTION
#---------------------------

def density_Powerlaw(r_max, rho_mean, q, GRID):

#r_max: Maximum radius of the envelope 
#rho_mean: Mean density of the Envelope 
#q: power-law for density
#GRID: Grid to work in

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

#-----------------------------------
#DENSITY (Keto2003, HCH_II) FUNCTION
#-----------------------------------

def density_Keto_HII(MStar, r_min, r_max, rho_s, T, GRID, q = 1.5):

#r_min: Minimum radius of the envelope (must be > 0)
#r_max: Maximum radius of the envelope 
#r_s: Reference radius
#rho_s: Density at r_s
#q: power-law for density
#GRID: Grid to work in

    #----------------
    #LOCAL PARAMETERS
    #----------------
    gamma = 5./3 #Heat capacity ratio for a monoatomic gas 
    cs = (gamma * T * kb / Mu)**0.5 #Speed of sound for H_II at T
    rs = 0.5 * G * MStar / cs**2 #Sonic point (v_escape = v_sound)
    if r_min > rs: sys.exit('ERROR: Please define r_min <= rs')

    #------------
    #LISTS TO USE
    #------------
    rList, NPoints = GRID.rRTP[0], GRID.NPoints #Due to spherical symmetry only r is needed

    #-------------------------------------------------
    #MODEL. HCH_II region, Keto 2003 (double gradient)
    #-------------------------------------------------
    print ('Calculating H_II Envelope density with Keto2003 (double gradient)...')
    print ('Speed of sound (cs):', cs/1e3, 'km/s')
    print ('Sonic point (rs):', rs/AU, 'au')
    rhoENV = np.ones(GRID.NPoints)
    ind_gravity = np.where((rList >= r_min) & (rList <= rs))
    ind_pressure = np.where((rList > rs) & (rList <= r_max))
    rhoENV[ind_gravity] = rho_s * (rs / rList[ind_gravity])**q
    rhoENV[ind_pressure] = rho_s * np.exp(-2 * (1 - rs / rList[ind_pressure]))
    #------------------------
    #------------------------

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'total': rhoENV, 'disc': np.zeros(NPoints), 'env': rhoENV, 'discFlag': False, 'envFlag': True, 'r_disc': False, 'r_min': r_min, 'r_env': r_max, 'rs': rs} ) 

#---------------------------
#---------------------------

#-----------------------------------
#DENSITY (PowerLaw, HCH_II) FUNCTION
#-----------------------------------

def density_Powerlaw_HII(r_min, r_max, r_s, rho_s, q, GRID):

#r_min: Minimum radius of the envelope (must be > 0)
#r_max: Maximum radius of the envelope 
#r_s: Reference radius
#rho_s: Density at r_s
#q: power-law for density
#GRID: Grid to work in

    #------------
    #LISTS TO USE
    #------------
    rList, NPoints = GRID.rRTP[0], GRID.NPoints #Due to spherical symmetry only r is needed

    #------------------------------
    #MODEL. HCH_II region, Powerlaw 
    #------------------------------
    print ('Calculating H_II Envelope density with power-law...')
    rhoENV = np.where((rList >= r_min) & (rList <= r_max), rho_s * (r_s / rList)**q, 1.)
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
        if not renv_max: renv_max = Rd #np.max(GRID.XYZgrid[0])
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
#GRID: Grid to work in

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

def temperature_Hamburgers(TStar, RStar, MStar, MRate, Rd, T10Env, T_min, BT, density, GRID, 
                           p = 0.33, inverted = False):

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

    #-----------------------------------------------------------
    #MODEL. Galvan-Madrid et al. 2018 + Chin-Fei Lee et al. 2017
    #-----------------------------------------------------------

    #------------
    #DISC Profile
    #------------
    if density.discFlag:
        print ('Calculating Burger-disc temperature...')
        zList = XYZ[2]
        Rdisc = density.r_disc
        H = density.H
        T_R = BT * (3*G * MStar * MRate / (4*np.pi * sigma * RList**3) * (1 - (RStar / RList)**0.5))**0.25

        if inverted: 
            print ('Set inverted temperature for Burger-disc...')
            if density.Rt: 
                tempDISC = np.where( RList < density.Rt, T_R * np.exp(- 0.5 * zList**2 / H**2), 1.0) 
                tempDISC = np.where( (RList >= density.Rt) & (RList <= Rdisc), T_R, tempDISC)
            else: tempDISC = np.where( RList <= Rdisc, T_R * np.exp(- 0.5 * zList**2 / H**2), 1.0) #Maximum in z = 0
        else: 
            print ('Set not inverted temperature for Burger-disc...')
            if density.Rt: 
                tempDISC = np.where( RList < density.Rt, T_R * np.exp(- 0.5 * (abs(zList) - H)**2 / H**2), 1.0) 
                tempDISC = np.where( (RList >= density.Rt) & (RList <= Rdisc), T_R, tempDISC)
            else: tempDISC = np.where( RList <= Rdisc, T_R * np.exp(- 0.5 * (abs(zList) - H)**2 / H**2), 1.0) #Maximum in z = H

        tempDISC = np.where( (RList <= Rdisc) & (tempDISC <= T_min), T_min, tempDISC)

    else: tempDISC = 1.

    #----------------
    #ENVELOPE PROFILE
    #----------------
    if density.envFlag:
        print ('Calculating Envelope temperature...')
        renv = density.r_env
        tempENV = np.where( rList <= renv, T10Env * 10**p * (rList / AU)**-p, 1.0)
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

def temperature_Constant(density, GRID, discTemp = 0, envTemp = 0, backTemp = 30.0):

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
            tempDISC = np.where( RList <= density.r_disc, discTemp, backTemp)
        else: sys.exit('ERROR: The disc calculation was turned ON but there is no density distribution for disc!')
    else: tempDISC = 0.

    #----------------
    #ENVELOPE PROFILE
    #----------------
    if envTemp:
        if density.envFlag:
            print ('Setting constant Envelope temperature...')
            tempENV = np.where( rList <= density.r_env, envTemp, backTemp)
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

    dx = GRID.step
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

    for i in xrange( len(Rot[1:]) ): 
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
        POS_vec = np.array([ np.dot( Rot_total, next(XYZ_it) ) for i in xrange(NPoints) ])
        rotPos = True
        
        if vel:
            print ('Rotating Velocity vectors...')
            VEL_it = iter(VEL_vec)
            VEL_vec = np.array([ np.dot( Rot_total, next(VEL_it) ) for i in xrange(NPoints) ])
            rotVel = True
        else: 
            print ('=================================================') 
            print ('WARNING: No VELOCITY distribution was provided to rotate!')
            print ('Please provide it if you are going to calculate line emission.')
            print ('=================================================')

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

#------------------------
#WRITING DATA (LIME v1.6)
#------------------------

class Make_Datatab(object):
   
    def __init__(self, prop, GRID):

        if len(np.shape(prop)) == 1: 
            self.prop = [prop]
            self.n = 1
        else: 
            self.prop = prop
            self.n = len(prop)
        self.GRID = GRID
        self.id = np.arange(GRID.NPoints)
        super(Make_Datatab, self).__init__()
    
    def formatter(self, format, tmp = '%d'):

        fmt, type_fmt = format, type(format) 
        nvec = xrange(self.n)

        if type_fmt == str: #If a single format is provided
            print ("Using format '%s'"%fmt) 
            for i in nvec: tmp += ' '+fmt #The same format for all properties
        elif type_fmt == list or type_fmt == np.ndarray: #If a list of formats
            if len(fmt) != self.n: sys.exit('ERROR: The number of formats provided (%d) is not equal to the number of properties to be written (%d)'%(len(fmt),self.n))
            print ('Using format list:', fmt) 
            for f in fmt: tmp += ' '+f
        elif not fmt: #If False
            print ("Using default format '%e'")
            for i in nvec: tmp += ' %e' #Default format for all properties
        else: sys.exit("ERROR: Wrong type: %s. \nPlease provide a valid 'format_list' object (str, list or np.ndarray)"%type_fmt)
        tmp += '\n'
        return tmp

    def submodel(self, tag, format = False, folder = 'Subgrids'):        

        os.system('mkdir %s'%folder)
        file_path = './%s/datatab_%s.dat'%(folder,tag)
        x,y,z = self.GRID.XYZ
        
        tmp = self.formatter(format, tmp = '%d %.8e %.8e %.8e')
        tmp_write = []
        if type(self.prop) == np.ndarray: self.prop = self.prop.tolist()
        list2write = iter(np.array([self.id,x,y,z] + self.prop).T)
        print ('Writing Submodel data in %s'%file_path)        
        for i in self.id: tmp_write.append( tmp % tuple(next(list2write)) )
        file_data = open(file_path, 'w')        
        file_data.writelines(tmp_write)

        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------\n-------------------------------------------------')
                

class Lime(Make_Datatab):
    def __init__(self, prop, GRID):
        print ('Set LIME format')
        super(Lime, self).__init__(prop, GRID)

    def globalgrid(self, format = False, folder = './'):
        
        tmp = self.formatter(format)
        tmp_write = []
        if type(self.prop) == np.ndarray: self.prop = self.prop.tolist()
        list2write = iter(np.array([self.id] + self.prop).T)
        files = [folder + tag for tag in ['datatab.dat','x.dat','y.dat','z.dat']]
        print ('Writing Global grid data in %s'%files[0])        
        for i in self.id: tmp_write.append( tmp % tuple(next(list2write)) )
        file_data = open(files[0],'w')
        file_data.writelines(tmp_write)

        Ns = self.GRID.Nodes
        size_file = folder+'npoints.dat'
        sfile = open(size_file,'w') 
        print ('Writing grid size in %s'%size_file)
        sfile.write("%d %d %d %d"%(Ns[0], Ns[1], Ns[2], self.GRID.NPoints))
        sfile.close()

        print ('Writing space domain in:')
        for i in xrange(1,4):
            print ('%s'%files[i])
            np.savetxt(files[i], self.GRID.XYZgrid[i-1], fmt = '%.8e')

        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------\n-------------------------------------------------')

    
class Radmc3d(object): #RADMC-3D uses the cgs units system
    """
    def __init__(self, prop, GRID):
        print ('Set RADMC-3D format')
        super(Radmc3d, self).__init__(prop, GRID)
    """
    
    def __init__(self, prop_dict, GRID, amr_grid = False, nphot = 1000000):
        print ('Set RADMC-3D format')
        self.prop = prop_dict
        self.GRID = GRID
        self.amr_grid = amr_grid
        self.nphot = nphot
        self.nx, self.ny, self.nz = GRID.Nodes
        self.nn = self.nx * self.ny * self.nz
        self.xi, self.yi, self.zi = np.array(GRID.XYZgrid) * cm #from m to cm
        
    def write_amr_grid(self, iformat = 1, 
                       grid_style = 0, 
                       coord_system = 0,
                       grid_info = 0, 
                       include_dim = [1,1,1]):
        
        #------------------------
        #Write the grid-info file
        #------------------------
        if not self.amr_grid:
            with open('amr_grid.inp','w+') as f:
                f.write('%s\n'%iformat)                             # iformat
                f.write('%s\n'%grid_style)                          # AMR grid style  (0=regular grid, no AMR)
                f.write('%s\n'%coord_system)                        # Coordinate system
                f.write('%s\n'%grid_info)                           # grid_info
                f.write('%s %s %s\n'%tuple(include_dim))            # Include x,y,z coordinate
                f.write('%d %d %d\n'%(self.nx,self.ny,self.nz))     # Size of grid
                for value in self.xi: f.write('%13.6e\n'%(value))   # X coordinates (cell walls)
                for value in self.yi: f.write('%13.6e\n'%(value))   # Y coordinates (cell walls)
                for value in self.zi: f.write('%13.6e\n'%(value))   # Z coordinates (cell walls)
                f.close()

        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------')
        
    def write_electron_numdens(self, dens_elect, format = '%13.6e'):

        #---------------------------------
        #Write the electronic density file
        #---------------------------------
        with open('electron_numdens.inp','w+') as f:
            f.write('1\n')                                          # Format number
            f.write('%d\n'%self.nn)                                 # Nr of cells
            #data = dens_elect.ravel(order='F') # Create a 1-D view, fortran-style indexing
            dens_elect.tofile(f, sep='\n', format=format)
            f.write('\n')

        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------')

    def write_ion_numdens(self, dens_ion, format = '%13.6e'):

        #--------------------------
        #Write the ion density file
        #--------------------------
        with open('ion_numdens.inp','w+') as f:
            f.write('1\n')                                          # Format number
            f.write('%d\n'%self.nn)                                 # Nr of cells
            #data = dens_ion.ravel(order='F') # Create a 1-D view, fortran-style indexing
            dens_ion.tofile(f, sep='\n', format=format)
            f.write('\n')

        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------')

    def write_gas_temperature(self, tgas, format = '%13.6e'):
        
        #-------------------------
        #Write the gas temperature
        #-------------------------
        with open('gas_temperature.inp','w+') as f:
            f.write('1\n')                                          # Format number
            f.write('%d\n'%self.nn)                                 # Nr of cells
            #data = tgas.ravel(order='F') # Create a 1-D view, fortran-style indexing
            tgas.tofile(f, sep='\n', format=format)
            f.write('\n')

        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------')

    def write_radmc3d_control(self, scattering_mode_max = 1, 
                              incl_freefree = 1,
                              incl_dust = 1,
                              tgas_eq_tdust = 1):

        #----------------------------------
        #Write the radmc3d.inp control file
        #----------------------------------
        with open('radmc3d.inp','w+') as f:
            f.write('nphot = %d\n'%(self.nphot))
            f.write('scattering_mode_max = %s\n'%scattering_mode_max)   # Put this to 1 for isotropic scattering
            f.write('incl_freefree = %s\n'%incl_freefree)
            f.write('incl_dust = %s\n'%incl_dust)
            f.write('tgas_eq_tdust = %s'%tgas_eq_tdust)

        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------')

    def write_wavelength_micron(self, lam = [5e2,2e4,4e4,3e5], nxx = [50,50,50], format = '%13.6e'):

        #------------------------------------
        #Write the wavelength_micron.inp file
        #------------------------------------
        len_lam = len(lam)
        if len_lam - 1 == len(nxx):
            lam_list = [np.logspace(np.log10(lam[i]),
                                    np.log10(lam[i+1]),
                                    nxx[i], endpoint=False) 
                        for i in xrange(len_lam-2)]
            lam_list.append(np.logspace(np.log10(lam[-2]),
                                        np.log10(lam[-1]),
                                        nxx[-1], endpoint=True))
            lam_list = np.concatenate(lam_list)
            nlam = lam_list.size
            with open('wavelength_micron.inp','w+') as f:
                f.write('%d\n'%(nlam))
                tmp = format+'\n'
                for value in lam_list: f.write(tmp%(value))

        else: sys.exit("ERROR: Wrong length(s) for input list(s): len(lam)-1 must be equal to len(nxx)")

        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------')
        
    def freefree(self, format = '%13.6e', folder = './'):

        prop = {} #Created a new dict bcause dont want to modify the self.prop variable at the minute
        prop['dens_elect'] = self.prop['dens_elect'] * cm**-3 
        prop['dens_ion'] = self.prop['dens_ion'] * cm**-3
        prop['tgas'] = self.prop['tgas']
        
        self.write_amr_grid()
        self.write_electron_numdens(prop['dens_elect'], format=format)
        self.write_ion_numdens(prop['dens_ion'], format=format)
        self.write_gas_temperature(prop['tgas'])
        self.write_radmc3d_control(incl_dust = 0, tgas_eq_tdust = 0)
        self.write_wavelength_micron()
            
        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------\n-------------------------------------------------')

#-----------------------------
#WRITING DATA (RADMC-3D v0.41)
#-----------------------------

def Datatab_RADMC3D_FreeFree(dens,temp,GRID):

    #dens = 1e-6 * np.where(dens > 10.0, dens, 0) #to cm^-3
    dens = dens / 1e6
    nx,ny,nz = GRID.Nodes
    xi, yi, zi = np.array(GRID.XYZgrid) * 100 #to cm
    nphot = 1000000

#
# Write the grid file
#
    with open('amr_grid.inp','w+') as f:
        f.write('1\n')                       # iformat
        f.write('0\n')                       # AMR grid style  (0=regular grid, no AMR)
        f.write('0\n')                       # Coordinate system
        f.write('0\n')                       # gridinfo
        f.write('1 1 1\n')                   # Include x,y,z coordinate
        f.write('%d %d %d\n'%(nx,ny,nz))     # Size of grid
        for value in xi:
            f.write('%13.6e\n'%(value))      # X coordinates (cell walls)
        for value in yi:
            f.write('%13.6e\n'%(value))      # Y coordinates (cell walls)
        for value in zi:
            f.write('%13.6e\n'%(value))      # Z coordinates (cell walls)
#
# Write the electronic density file.
#
    with open('electron_numdens.inp','w+') as f:
        f.write('1\n')                       # Format number
        f.write('%d\n'%(nx*ny*nz))           # Nr of cells
        #data = rhoelect.ravel(order='F')         # Create a 1-D view, fortran-style indexing
        dens.tofile(f, sep='\n', format="%13.6e")
        f.write('\n')

#
# Write the ion density file.
#
    with open('ion_numdens.inp','w+') as f:
        f.write('1\n')                       # Format number
        f.write('%d\n'%(nx*ny*nz))           # Nr of cells
        #data = rhoelect.ravel(order='F')         # Create a 1-D view, fortran-style indexing
        dens.tofile(f, sep='\n', format="%13.6e")
        f.write('\n')
    
#
# Write the gas temperature
#

    with open('gas_temperature.inp','w+') as f:
        f.write('1\n')                       # Format number
        f.write('%d\n'%(nx*ny*nz))           # Nr of cells
        #data = tgas.ravel(order='F')          # Create a 1-D view, fortran-style indexing
        temp.tofile(f, sep='\n', format="%13.6e")
        f.write('\n')

#
# Write the wavelength_micron.inp file
#

    lam1 = 0.5e3
    lam2 = 2.e4
    lam3 = 4.e4
    lam4 = 3.e5#6.e4
    
    n12      = 50
    n23      = 50
    n34      = 50
    lam12    = np.logspace(np.log10(lam1),np.log10(lam2),n12,endpoint=False)
    lam23    = np.logspace(np.log10(lam2),np.log10(lam3),n23,endpoint=False)
    lam34    = np.logspace(np.log10(lam3),np.log10(lam4),n34,endpoint=True)
    lam      = np.concatenate([lam12,lam23,lam34])
    nlam     = lam.size

#
# Write the wavelength file
#
    with open('wavelength_micron.inp','w+') as f:
        f.write('%d\n'%(nlam))
        for value in lam:
            f.write('%13.6e\n'%(value))

#
# Write the radmc3d.inp control file
#
    with open('radmc3d.inp','w+') as f:
        f.write('nphot = %d\n'%(nphot))
        f.write('scattering_mode_max = 1\n')   # Put this to 1 for isotropic scattering
        f.write('incl_freefree = 1\n')
        f.write('incl_dust = 0\n')
        #f.write('tgas_eq_tdust = 1')

    
    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

#-------------
#-------------

def DataTab_LIME(dens,temp,vel,abund,gtd,GRID, is_submodel = False, tag = False):
    
    import pandas

    if is_submodel:
        os.system('mkdir Subgrids')
        file0 = './Subgrids/datatab%s.dat'%tag
        file = open(file0,'w')
        x,y,z = GRID.XYZ
        print ('Writing Submodel data on %s'%file0)
        tmp = []
        for i in xrange(GRID.NPoints): 
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

        for i in xrange(GRID.NPoints): 
            file.write("%d %e %e %e %e %e %e %e\n"%
                       (i,dens[i],temp[i],vel.x[i],vel.y[i],vel.z[i],abund[i],gtd[i]))

        df = [pandas.DataFrame(GRID.XYZgrid[i]) for i in range(3)]
    
        for i in xrange(1,4):
            print ('Writing data on %s'%files[i])
            df[i-1].to_csv(files[i],index=False,header=False,float_format='%e') 
        
        sfile.close()
        
    file.close()
    
    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

#--------------
#--------------    

def Make_Datatab1(prop_list, GRID, format_list = False, 
                 submodel_tag = False, submodel_folder = 'Subgrids', 
                 lime = True, radmc3d = False):
    
    import pandas
    n = len(prop_list)

    if submodel_tag:
        
        fold, tag, fmt, type_fmt = submodel_folder, submodel_tag, format_list, type(format_list)
        os.system('mkdir %s'%fold)
        file_path = './%s/datatab_%s.dat'%(fold,tag)
        x,y,z = GRID.XYZ
        tmp = '%d %e %e %e'

        if type_fmt == str: #If a single format is provided
            print ("Using format '%s'"%fmt) 
            for i in range(n): tmp += ' '+fmt #The same format for all properties
        elif type_fmt == list or type_fmt == np.ndarray: #If a list of formats
            if len(fmt) != n: sys.exit('ERROR: The number of formats provided (%d) is not equal to the number of properties to be written (%d)'%(len(fmt),n))
            print ('Using format list:', fmt) 
            for f in fmt: tmp += ' '+f
        elif not fmt: #If False
            print ("Using default format '%e'")
            for i in range(n): tmp += ' %e' #Default format for all properties
        else: sys.exit("ERROR: Wrong type: %s. \nPlease provide a valid 'format_list' object (str, list or np.ndarray)"%type_fmt)

        tmp += '\n'
        tmp_write = []
        if type(prop_list) == np.ndarray: prop_list = prop_list.tolist()
        id = np.arange(GRID.NPoints)
        list2write = iter(np.array([id,x,y,z] + prop_list).T)
        file = open(file_path, 'w')
        print ('Writing Submodel data on %s'%file_path)        
        for i in id:    
            tmp_write.append( tmp % tuple(next(list2write)) )
        
        file.writelines(tmp_write)
            
    else:
        if lime:
            pass
        if radmc3d:
            pass

    file.close()
