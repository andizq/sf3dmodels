from __future__ import print_function
import numpy as np
import random
import inspect
import time
import sys
import os

from .Utils import *
from .utils.constants import mH

class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)

#------------------------
#SPATIAL (Spherical-)GRID
#------------------------
 
def grid(XYZmax, NP, rt_code = 'lime', include_zero = True, indexing='ij'):
    """
    Computes the hosting grid for the model(s).
    
    Parameters
    ----------
    
    Returns
    -------    
    
    """
    """
    XYZmax: Maximum physical domain of the grid
    Please provide XYZmax as a 3D-list. XYZmax=[Xmax,Ymax,Zmax]
    NP: Number of grid points for each coordinate
    Please provide NP as a 3D-list. NP=[Nx,Ny,Nz]
    """

    print ('-------------------------------------------------\n-------------------------------------------------')
    
    XYZmax = np.asarray(XYZmax)
    NP = np.asarray(NP).astype('int32')
    #----------------
    #NUMBER OF POINTS
    #----------------
    #if include_zero NP becomes odd where it's not in order to force the grid to contain the midplane(s).
    
    if include_zero:
        NP_dum = np.array([ NP[i] + 1 if NP[i]%2 == 0 else NP[i] for i in range(3) ])
        print('Setting the number of grid points to be {},'.format(NP_dum) +
              '\n so that the planes x=0, y=0, z=0 are all present...' +
              '\nTurn off this feature by setting the flag include_zero=False'
              if (NP != NP_dum).any()
              else '... ... ...'
              )
        NP = NP_dum
    else: 
        for i in range(3): print('Coordinate zero NOT included for axis %d'%i
                                  if NP[i]%2 == 0
                                  else 'Coordinate zero included for axis %d'%i)
    #--------
    #XYZ GRID
    #--------
    print ('Computing Grid...')
    
    Xmax,Ymax,Zmax = XYZmax    
    #epsilon = [RSun / 1.] * 3
    step = [2. * XYZmax[i] / (NP[i]-1) if NP[i]>1 else 0 for i in range(3)]

    if rt_code == 'radmc3d': 
        step = 2. * XYZmax / (NP)
        XYZgrid = [np.linspace(-XYZmax[i], XYZmax[i], NP[i] + 1) for i in range(3)]
        #XYZgrid = [np.append( np.linspace(-XYZmax[i], XYZmax[i], NP[i]), (XYZmax[i] + step[i]) ) for i in range(3)]
        X, Y ,Z = XYZgrid #The grid must contain the cell corners, but the properties are computed at the cell centres. 
        X = 0.5 * ( X[0:NP[0]] + X[1:NP[0]+1] ) #Moving the node from the corner to the center of the boxel. length = lengthofcorners - 1
        Y = 0.5 * ( Y[0:NP[1]] + Y[1:NP[1]+1] )  
        Z = 0.5 * ( Z[0:NP[2]] + Z[1:NP[2]+1] )
        XYZcentres = [X,Y,Z]

        ZYX = np.meshgrid(Z,Y,X, indexing=indexing)
        zList, yList, xList = [zyx.flatten() for zyx in ZYX]
        #X = X[:-1]; Y = Y[:-1]; Z = Z[:-1] #...the calculations must be done w/o that node 

    elif rt_code == 'lime': #lime
        XYZgrid = [np.linspace(-XYZmax[i], XYZmax[i], NP[i]) for i in range(3)]
        print (step)
        X, Y, Z = XYZgrid
        XYZcentres = XYZgrid

        XYZ = np.meshgrid(X,Y,Z, indexing=indexing)
        xList, yList, zList = [xyz.flatten() for xyz in XYZ]

    #--------------------------------------
    #Extended Lists of distance coordinates
    #--------------------------------------
    """    
    if rt_code == 'radmc3d': rRxyzList = np.array([ ((x**2 + y**2 + z**2)**0.5, (x**2 + y**2)**0.5, x,y,z) for z in Z for y in Y for x in X])
    elif rt_code == 'lime': rRxyzList = np.array([ ((x**2 + y**2 + z**2)**0.5, (x**2 + y**2)**0.5, x,y,z) for x in X for y in Y for z in Z])
    rList = rRxyzList[:,0] ; RList = rRxyzList[:,1]; xList = rRxyzList[:,2]; yList = rRxyzList[:,3]; zList = rRxyzList[:,4]
    """
    
    rList = np.linalg.norm([xList,yList,zList], axis = 0)
    RList = np.linalg.norm([xList,yList], axis = 0)
    
    r_ind_zero, = np.where(rList < 1.)
    R_ind_zero, = np.where(RList < 1.)
    if len(r_ind_zero)>0: rList[r_ind_zero] = 0.5*np.unique(rList)[1]
    if len(R_ind_zero)>0: RList[R_ind_zero] = 0.5*np.unique(RList)[1]
  
    """
    rList = np.where(rList < 1., 0.5*rList_unique[1], rList) # If r == 0: use half of the second minimum value of r
    RList = np.where(RList < 1., 0.5*RList_unique[1], RList)
    """
    #-----
    #THETA
    #-----
    thetaList = np.arccos(zList / rList)
    halfPi = 0.5 * np.pi
    theta4vel = thetaList

    #If theta is under the xy plane, take its image above the plane.
    # This is done to avoid negative cosines in func streamline(), 
    #  taking advantage of the symmetry of discs and envelopes along xy.         
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
    
    return Struct( **{'XYZgrid': XYZgrid, 'XYZcentres': XYZcentres, 
                      'XYZ': XYZ, 'rRTP': rRTP, 'theta4vel': theta4vel, 
                      'NPoints': len(rList), 'Nodes': NP, 'step': step,
                      'r_ind_zero': r_ind_zero, 'R_ind_zero': R_ind_zero})

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

#--------------------------------------------
#Vector in Spherical base --> Cartesian base
#--------------------------------------------

def sphe_cart(V_sphe, theta, phi):
    
    #--------------------------------------------------------------
    #Wikipedia: Coordenadas esfericas, subseccion 'Base coordenada'
    #--------------------------------------------------------------
   
    if hasattr(theta, '__iter__') and hasattr(phi, '__iter__'): #A list of vectors
    #np.multiply and operator * are esentially the same even in terms of time
        vx = np.multiply( V_sphe, np.array([np.sin(theta) * np.cos(phi), np.cos(theta) * np.cos(phi), -np.sin(phi)]).T ).sum(1) #along axis 1
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

def density_Env_Disc(RStar, Rd, rhoE0, Arho, GRID, 
                     discFlag = True, envFlag = False, 
                     rdisc_max = False, renv_max = False, 
                     ang_cavity = False, rho_min_env = 1.0e9):

#RStar: stellar radius
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
        print ('Computing Keplerian flared-disc density...')
        
        if not rdisc_max: rdisc_max = Rd #Keto's value for the disc to stop
        rhoD0 = Arho * rhoE0 #Normalization factor based on the envelope
        H0 = 0.01 * RStar #Scaleheight at RStar
        H = H0 * (RList / RStar)**1.25  #Scaleheight
        rhoDISC = np.where( RList <= rdisc_max, rhoD0 * (Rd / RList)**2.25 * np.exp(-0.5 * zList**2 / H**2), 1.0)#1.0e9 / 2)
        rhoDISC = np.where( rhoDISC < 1.0, 1.0, rhoDISC)
    else: 
        print ('No Keplerian flared-disc was invoked!')
        rhoDISC = np.zeros(NPoints)
        H = None
        
    #----------------
    #ENVELOPE PROFILE
    #----------------
    if ang_cavity: print ('Set cavity for density with half-aperture %.1f deg'%(ang_cavity * 180 / np.pi))
    
    if envFlag:
        print ('Computing stream lines for Ulrich envelope...')
        if not renv_max: renv_max = np.max(XYZgrid[0]) #A sphere inscribed in the coordinate X. ##CHANGE this by the smallest of the 3 maximum (1 for each axis)

        costheta = np.cos(thetaList)
        costheta0 = streamline(Rd,GRID)
        print ('Computing Envelope density...')
        rhoENV = np.where( (rList <= renv_max) & (thetaList >= ang_cavity), 
                           ((rhoE0 * (rList / Rd)**-1.5) *
                            ((1 + (costheta / costheta0))**-0.5) *
                            (1 + (Rd / rList) * (3. * costheta0**2 - 1))**-1), 
                           rho_min_env )
        rhoENV = np.where( rhoENV < 1.0, 1.0, rhoENV)
    else:
        print ('No Envelope was invoked!')
        costheta0 = False
        rhoENV = np.zeros(NPoints)
    
    #----------------------------------------------
    #----------------------------------------------
    
    RHO = rhoDISC + rhoENV 

    nonzero_ids = np.where(RHO != 0.0)
    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'total': RHO, 'disc': rhoDISC, 'env': rhoENV, 'H': H,
                      'discFlag': discFlag, 'envFlag': envFlag, 
                      'r_disc': rdisc_max, 'r_env': renv_max, 'streamline': costheta0,
                      'nonzero_ids': nonzero_ids} )

#-------------------------------------
#DENSITY (viscously ev. disc) FUNCTION
#-------------------------------------

def density_lyndenbell_disc(GRID, Rc=100*AU, Ec=30.0, gamma=1.0, H0=6.5*AU, psi=1.25, Ro=500*AU,
                            rho_thres = 10.0*mH, rho_min = 2*mH,
                            discFlag=True, rdisc_max = False):

#Rc: characteristic radius
#Ec: mass surface density normalisation at Rc
#gamma: viscous power-law exponent
#H0: scale height normalisation at Rc
#psi: flaring index
#Ro: outer radius of the disk
#rho_thres: densities below this value will acquire a density val of rho_min 
    
    XYZ, rRTP = GRID.XYZ, GRID.rRTP

    #------------
    #LISTS TO USE
    #------------
    zList = XYZ[2]
    rList, RList = rRTP[:-2] 

    #-------------------------------------------------------
    #MODEL. Hartmann et al. 1998, Izquierdo et al. 2021(A&A)
    #-------------------------------------------------------

    #------------
    #DISC PROFILE
    #------------
    if discFlag:
        print ('Computing Lynden-Bell disc density...')
        if not rdisc_max: rdisc_max = Ro
        print ('Scale-height at Rc (au):', H0/AU)
        H = H0*(RList/Rc)**psi #Scaleheight powerlaw
        rhoDISC = np.where( RList <= rdisc_max,
                            Ec*(RList/Rc)**-gamma*np.exp(-1*(RList/Rc)**(2-gamma)) * np.exp(-0.5*(zList/H)**2)/(np.sqrt(2*np.pi)*H), rho_min)
        rhoDISC = np.where( rhoDISC < rho_thres, rho_min, rhoDISC)
    else: 
        print ('No Disc was invoked!')
        rhoDISC = np.zeros(GRID.NPoints)
        
    #----------------------------------------------
    #----------------------------------------------

    nonzero_ids = np.where(rhoDISC != 0.0)
    #rhoDISC[GRID.R_ind_zero] = 0 #*np.mean(rhoDISC)
    
    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'total': rhoDISC, 'disc': rhoDISC, 'env': 0., 'H': H,  
                      'discFlag': discFlag, 'envFlag': False, 'Rc': Rc, 
                      'r_disc': rdisc_max, 'r_env': False,
                      'nonzero_ids': nonzero_ids} ) 

#------------------------------
#DENSITY (Hamburguers) FUNCTION
#------------------------------

def density_Hamburgers(RStar, shFactor, Ro, rhoE0, Arho, GRID, 
                       p = 2.25, q = 0.5, rho_thres = 10.0, rho_min = 1.0, 
                       Rt = False, discFlag=True, rdisc_max = False):

#RStar: stellar radius
#shFactor: scaleheight normalization constant --> H0 = shFactor * RStar
#Ro: outer radius of the disk
#Rt: radius where tapering commences
#rhoE0: density at Rd and theta=pi/2
#q: scaleheight exponent
#p: density scaling exponent 
#Arho: density scaling factor 
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
        print ('Computing Burger-disc density...')
        if not rdisc_max: rdisc_max = Ro
        rhoD0 = Arho * rhoE0 
        H0 = shFactor * RStar
        print ('Scale-height at RStar (au):', H0 / AU * 1 / ((RStar/AU)**(1 + 0.5*(1-q))))
        H = H0 * (RList / RStar)**(1 + 0.5*(1-q)) #Scaleheight, with no tapering 
        if Rt: H = H * np.exp(-((RList-Rt) / (Ro-Rt))**2)  #Scaleheight, with tapering 
        rhoDISC = np.where( RList <= rdisc_max, rhoD0 * (RList / Ro)**-p * np.exp(-0.5 * zList**2 / H**2), rho_min)
        rhoDISC = np.where( rhoDISC < rho_thres, rho_min, rhoDISC)
    else: 
        print ('No Disc was invoked!')
        rhoDISC = np.zeros(GRID.NPoints)
        
    #----------------------------------------------
    #----------------------------------------------

    nonzero_ids = np.where(rhoDISC != 0.0)
    #rhoDISC[GRID.R_ind_zero] = 0 #*np.mean(rhoDISC)
    
    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'total': rhoDISC, 'disc': rhoDISC, 'env': 0., 'H': H,  
                      'discFlag': discFlag, 'envFlag': False, 'Rt': Rt, 
                      'r_disc': rdisc_max, 'r_env': False,
                      'nonzero_ids': nonzero_ids} ) 

#------------------------------------------
#DENSITY (Hamburguers - piecewise) FUNCTION
#------------------------------------------

def density_Hamburgers_piecewise(RStar, H0, R_list, p_list, rho0, GRID, RH_list = None,
                                 q_list = [0.5], rho_thres = 10.0, rho_min = 0.0, 
                                 Rt = False):

#RStar: stellar radius
#H0: scaleheight normalization constant --> usually H0 = shFactor * R0
#R_list: List of polar limits, length (n,)
#p_list: List of powerlaws in R_list intervals, length (n-1,)
#rho0: density at R_list[0]
#Optionals:
#RH_list: List of limits for piecewise scaleheight, length (n,)
#q_list: List of powerlaws in RH_list intervals, length (n-1,)
#rho_thres: minimum reachable density by the model
#rho_min: background density
#Rt: radius where the disc tapering starts

    XYZ, rRTP, NPoints = GRID.XYZ, GRID.rRTP, GRID.NPoints
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
    print ('Computing Hamburger-disc density profile using power-laws:', p_list)
    Rd = R_list[-1]
    rhoDISC = np.zeros(NPoints)
    rho0_coeff = [rho0]
    for i,R in enumerate(R_list[1:-1],1):
        R_tmp = np.max(RList[RList<=R])
        rho0_coeff.append(rho0_coeff[i-1]*(R_tmp/R_list[i-1])**p_list[i-1])
    for i,p in enumerate(p_list):
        ind, = np.where((RList >= R_list[i]) & (RList <= R_list[i+1]))
        rhoDISC[ind] += rho0_coeff[i]*(RList[ind]/R_list[i])**p  

    if RH_list is None:
        q = q_list[0]
        H = H0 * (RList / RStar)**(1 + 0.5*(1-q)) #Scaleheight, without tapering 
            
    else: 
        H = np.ones(NPoints) #ones instead of zeroes to avoid dividing by zero at empty spaces.
        H_coeff = [H0]
        for i,R in enumerate(RH_list[1:-1],1):
            RH_tmp = np.max(RList[RList<=R])
            H_coeff.append(H_coeff[i-1]*(RH_tmp/RH_list[i-1])**(1 + 0.5 - 0.5*q_list[i-1]))
        for i,q in enumerate(q_list):
            ind, = np.where((RList >= RH_list[i]) & (RList <= RH_list[i+1]))
            H[ind] += H_coeff[i]*(RList[ind]/RH_list[i])**(1 + 0.5 - 0.5*q)  
        
    if Rt: H = H * np.exp(-((RList-Rt) / (Rd-Rt))**2)  #Scaleheight, with tapering 
            
    rhoDISC = np.where( RList <= Rd, rhoDISC * np.exp(-0.5 * zList**2 / H**2), rho_min)
    rhoDISC = np.where( rhoDISC < rho_thres, rho_min, rhoDISC)
            
    #----------------------------------------------
    #----------------------------------------------

    nonzero_ids = np.where(rhoDISC != 0.0)
    #rhoDISC[GRID.R_ind_zero] = 0 #*np.mean(rhoDISC)
    
    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'total': rhoDISC, 'disc': rhoDISC, 'env': 0., 'H': H,  
                      'discFlag': True, 'envFlag': False, 'Rt': Rt, 
                      'r_disc': Rd, 'r_env': False,
                      'nonzero_ids': nonzero_ids} ) 

#-----------------------------------
#DENSITY (PowerLaw-mean_rho) FUNCTION
#-----------------------------------

def density_Powerlaw(r_max, rho_mean, q, GRID, rho_min = 1.0e3):

#r_max: Maximum radius of the envelope 
#rho_mean: Mean density of the Envelope 
#q: power-law for density
#GRID: Grid to work in
#rho_min: Minimum density

    #------------
    #LISTS TO USE
    #------------
    rList, NPoints = GRID.rRTP[0], GRID.NPoints #Due to spherical symmetry only r is needed

    #------------------------
    #MODEL. Envelope powerlaw
    #------------------------
    print ('Computing Envelope density using power-law...')
    rqList = np.where(rList <= r_max , rList**q, 0.)

    #As rho_mean = 1/NTotal * np.sum(rho0 * r**q), the normalization rho0 is calculated as follows:  
    rho0 = NPoints * rho_mean / np.sum(rqList)
    rhoENV = rho0 * rqList
    rhoENV = np.where(rhoENV < rho_min, rho_min, rhoENV)

    #rhoENV = np.where( rho0 * rqList < 1.0, rho_min, rho0 * rqList )
    
    #------------------------
    #------------------------
    
    nonzero_ids = np.where(rhoENV != 0.0)

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'total': rhoENV, 'disc': np.zeros(NPoints), 'env': rhoENV, 
                      'discFlag': False, 'envFlag': True, 'r_disc': False, 'r_env': r_max,
                      'nonzero_ids': nonzero_ids} ) 

#---------------------------
#---------------------------

#------------------------------------
#DENSITY (PowerLaw-standard) FUNCTION
#------------------------------------

def density_Powerlaw2(r_max, r_min, rho0, q, GRID, rho_min = 1.0e3):

#r_max: Maximum radius of the envelope 
#r_min: Minimum radius of the envelope 
#rho0: Density at r_min
#q: power-law for density
#GRID: Grid to work in
#rho_min: Minimum density

    #------------
    #LISTS TO USE
    #------------
    rList, NPoints = GRID.rRTP[0], GRID.NPoints #Due to spherical symmetry only r is needed

    #------------------------
    #MODEL. Envelope powerlaw
    #------------------------
    print ('Computing Envelope density using power-law...')
    rqList = np.where((rList >= r_min) & (rList <= r_max) , rList**q, 0.)
    rhoENV = rho0 * rqList #r_min**-q * rqList
    rhoENV = np.where(rhoENV < rho_min, rho_min, rhoENV)

    #------------------------
    #------------------------
    
    nonzero_ids = np.where(rhoENV != 0.0)

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'total': rhoENV, 'disc': np.zeros(NPoints), 'env': rhoENV, 
                      'discFlag': False, 'envFlag': True, 'r_disc': False, 'r_env': r_max,
                      'nonzero_ids': nonzero_ids} ) 

#---------------------------
#---------------------------

#------------------------------------
#DENSITY (PowerLaw-Shells) FUNCTION
#------------------------------------

def density_PowerlawShells(r_list, p_list, rho0, GRID, rho_min = 1.0e3):

#r_list: List of shells' limits, length (n,)
#p_list: List of powerlaws in r_list intervals, length (n-1,)
#rho0: Density at r_list[0]
#GRID: Grid to work in
#rho_min: Minimum density

    #------------
    #LISTS TO USE
    #------------
    rList, NPoints = GRID.rRTP[0], GRID.NPoints #Due to spherical symmetry only r is needed

    #-------------------------
    #MODEL. Envelope powerlaws
    #-------------------------
    print ('Computing Envelope density using power-laws:', p_list)
    
    rhoENV = np.zeros(NPoints)
    rho0_coeff = [rho0]
    for i,r in enumerate(r_list[1:-1],1):
        r_tmp = np.max(rList[rList<=r])
        rho0_coeff.append(rho0_coeff[i-1]*(r_tmp/r_list[i-1])**p_list[i-1])
    
    for i,p in enumerate(p_list):
        ind, = np.where((rList >= r_list[i]) & (rList <= r_list[i+1]))
        rhoENV[ind] += rho0_coeff[i]*(rList[ind]/r_list[i])**p

    rhoENV = np.where(rhoENV < rho_min, rho_min, rhoENV)

    #------------------------
    #------------------------

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'total': rhoENV, 'disc': np.zeros(NPoints), 'env': rhoENV, 
                      'discFlag': False, 'envFlag': True, 'r_disc': False, 'r_env': r_list[-1]
                      } ) 

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
    print ('Computing H_II Envelope density with Keto2003 (double gradient)...')
    print ('Speed of sound (cs):', cs/1e3, 'km/s')
    print ('Sonic point (rs):', rs/AU, 'au')
    rhoENV = np.ones(GRID.NPoints)
    ind_gravity = np.where((rList >= r_min) & (rList <= rs))
    ind_pressure = np.where((rList > rs) & (rList <= r_max))
    rhoENV[ind_gravity] = rho_s * (rs / rList[ind_gravity])**q
    rhoENV[ind_pressure] = rho_s * np.exp(-2 * (1 - rs / rList[ind_pressure]))
    #------------------------
    #------------------------

    nonzero_ids = np.where(rhoENV != 0.0)

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'total': rhoENV, 'disc': np.zeros(NPoints), 'env': rhoENV, 
                      'discFlag': False, 'envFlag': True, 'r_disc': False, 
                      'r_min': r_min, 'r_env': r_max, 'rs': rs,
                      'nonzero_ids': nonzero_ids} ) 

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
    print ('Computing H_II Envelope density with power-law...')
    rhoENV = np.where((rList >= r_min) & (rList <= r_max), rho_s * (r_s / rList)**q, 1.)
    #------------------------
    #------------------------

    nonzero_ids = np.where(rhoENV != 0.0)

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'total': rhoENV, 'disc': np.zeros(NPoints), 'env': rhoENV, 
                      'discFlag': False, 'envFlag': True, 'r_disc': False, 'r_env': r_max,
                      'nonzero_ids': nonzero_ids} ) 

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
        
    nonzero_ids = np.where(RHO != 0.0)

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'total': RHO, 'disc': rhoDISC, 'env': rhoENV, 
                      'discFlag': bool(discDens), 'envFlag': bool(envDens), 
                      'r_disc': rdisc_max, 'r_env': renv_max,
                      'nonzero_ids': nonzero_ids} ) 

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

#------------------------------------
#DENSITY (PowerLaw-standard) FUNCTION
#------------------------------------

def abundance_Powerlaw(r_max, r_min, ab0, q, GRID, ab_min = 1e-10):

#r_max: Maximum radius of the envelope 
#r_min: Minimum radius of the envelope 
#ab0: scaling coefficient
#q: power-law for density
#GRID: Grid to work in
#rho_min: Minimum density

    #------------
    #LISTS TO USE
    #------------
    rList, NPoints = GRID.rRTP[0], GRID.NPoints #Due to spherical symmetry only r is needed

    #------------------------
    #MODEL. Envelope powerlaw
    #------------------------
    print ('Computing Envelope density using power-law...')
    rqList = np.where((rList >= r_min) & (rList <= r_max) , rList**q, 0.)
    abundList = ab0 * rqList #r_min**-q * rqList
    abundList = np.where(abundList < ab_min, ab_min, abundList)

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return abundList 
    
#------------------------------------
#DENSITY (PowerLaw-Shells) FUNCTION
#------------------------------------

def abundance_PowerlawShells(r_list, p_list, abund0, GRID, abundance_min = 1e-15):

#r_list: List of shells' limits, length (n,)
#p_list: List of powerlaws in r_list intervals, length (n-1,)
#abund0: Density at r_list[0]
#GRID: Grid to work in
#abundance_min: Minimum density

    #------------
    #LISTS TO USE
    #------------
    rList, NPoints = GRID.rRTP[0], GRID.NPoints #Due to spherical symmetry only r is needed

    #-------------------------
    #MODEL. Envelope powerlaws
    #-------------------------
    print ('Computing Envelope density using power-laws:', p_list)
    
    abund = np.zeros(NPoints)
    abund0_coeff = [abund0]
    for i,r in enumerate(r_list[1:-1],1):
        r_tmp = np.max(rList[rList<=r])
        abund0_coeff.append(abund0_coeff[i-1]*(r_tmp/r_list[i-1])**p_list[i-1])
    
    for i,p in enumerate(p_list):
        ind, = np.where((rList >= r_list[i]) & (rList <= r_list[i+1]))
        abund[ind] += abund0_coeff[i]*(rList[ind]/r_list[i])**p

    abund = np.where(abund < abundance_min, abundance_min, abund)

    #------------------------
    #------------------------

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return abund

#---------------------------
#---------------------------

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

def temperature(TStar, Rd, T10Env, RStar, MStar, MRate, BT, density, GRID, 
                Tmin_disc = 30., Tmin_env = 30., p = 0.33, ang_cavity = False):

#TStar: stellar temperature
#T10Env: Envelope temperature at 10AU
#RStar: stellar radius
#MStar: stellar mass
#MRate: Mass accretion rate
#BT: Disc temperature factor
#p: (Envelope) Temperature power law exponent 
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
        print ('Computing Keplerian flared-disc temperature...')
        rdisc = density.r_disc
        #* np.exp(-0.5 * zList**2 / H**2)
        if density.envFlag:
            renv = density.r_env
            tempDISC = np.where( (RList <= rdisc) & (rList <= renv) , 
                                 BT * (3*G * MStar * MRate / (4*np.pi * sigma * RList**3) * (1 - (RStar / RList)**0.5))**0.25, 
                                 Tmin_disc)
        else: tempDISC = np.where(RList <= rdisc , 
                                  BT * (3*G * MStar * MRate / (4*np.pi * sigma * RList**3) * (1 - (RStar / RList)**0.5))**0.25, 
                                  Tmin_disc)
    else: tempDISC = 1.

    #----------------
    #ENVELOPE PROFILE
    #----------------
    if ang_cavity: print ('Set cavity for temperature with half-aperture %.1f deg'%(ang_cavity * 180 / np.pi))

    if density.envFlag:
        print ('Computing Envelope temperature...')
        renv = density.r_env
        tempENV = np.where( (rList <= renv) & (thetaList >= ang_cavity), 
                            T10Env * 10**p * (rList / AU)**-p, 
                            Tmin_env)
        tempENV = np.where(tempENV < Tmin_env, Tmin_env, tempENV)
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

def temperature_Hamburgers(TStar, RStar, MStar, MRate, Rd, T10Env, BT, density, GRID, 
                           p = 0.33, Tmin_disc = 30., Tmin_env = 30., inverted = False):

#TStar: stellar temperature
#T10Env: Envelope temperature at 10AU
#RStar: stellar radius
#MStar: stellar mass
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
        print ('Computing Burger-disc temperature...')
        zList = XYZ[2]
        Rdisc = density.r_disc
        H = density.H
        T_R = BT * (3*G * MStar * MRate / (4*np.pi * sigma * RList**3) * (1 - (RStar / RList)**0.5))**0.25

        if inverted: 
            print ('Set inverted temperature for Burger-disc...')
            T_z = np.exp(- 0.5 * zList**2 / H**2)
            tempDISC = np.where( RList <= Rdisc, T_R * T_z, 1.0) #Maximum in z = 0
        else: 
            print ('Set not inverted temperature for Burger-disc...')
            T_z = np.exp(- 0.5 * (abs(zList) - H)**2 / H**2)
            tempDISC = np.where( RList <= Rdisc, T_R * T_z, 1.0) #Maximum in z = H
            """
            if density.Rt:
               tempDISC = np.where( RList < density.Rt, T_R * T_z, 1.0) 
               tempDISC = np.where( (RList >= density.Rt) & (RList <= Rdisc), T_R, tempDISC)
            else: tempDISC = np.where( RList <= Rdisc, T_R * T_z, 1.0) #Maximum in z = H
            """

        tempDISC = np.where( (RList <= Rdisc) & (tempDISC <= Tmin_disc), Tmin_disc, tempDISC)

    else: tempDISC = 1.

    #----------------
    #ENVELOPE PROFILE
    #----------------
    if density.envFlag:
        print ('Computing Envelope temperature...')
        renv = density.r_env
        tempENV = np.where( rList <= renv, T10Env * 10**p * (rList / AU)**-p, Tmin_env)
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
            tempDISC = np.where( (RList <= density.r_disc), discTemp, backTemp) # & (abs(GRID.XYZ[2]) < density.H)
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

    #zerodens_mask = np.equal(density.total, 0.0)
    if envTemp and discTemp: TEMP = np.where(density.total != 0.0, (tempDISC * rhoDISC + tempENV * rhoENV) / density.total, backTemp)
    elif envTemp: TEMP = np.where(density.total != 0.0, tempENV, backTemp)
    elif discTemp: TEMP = np.where(density.total != 0.0, tempDISC, backTemp)
    
    #TEMP = np.choose(zerodens_mask, ((tempDISC * rhoDISC + tempENV * rhoENV) / density.total, backTemp))
        
    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'total': TEMP, 'disc': tempDISC, 'env': tempENV, 'discFlag': bool(discTemp), 'envFlag': bool(envTemp)} )

#-------------------------------
#-------------------------------

#-----------------------------------
#TEMPERATURE (PowerLaw-mean_rho) FUNCTION
#-----------------------------------

def temperature_Powerlaw(r_max, T_mean, q, GRID, T_min = 2.725):

#r_max: Maximum radius of the envelope 
#T_mean: Mean temperature of the Envelope 
#q: power-law for temperature
#GRID: Grid to work in
#T_min: Minimum temperature

    #------------
    #LISTS TO USE
    #------------
    rList, NPoints = GRID.rRTP[0], GRID.NPoints #Due to spherical symmetry only r is needed

    #------------------------
    #MODEL. Envelope powerlaw
    #------------------------
    print ('Computing Envelope temperature using power-law...')
    rqList = np.where(rList <= r_max , rList**q, 0.)

    #As T_mean = 1/NTotal * np.sum(T0 * r**q), the normalization T0 is calculated as follows:  
    T0 = NPoints * T_mean / np.sum(rqList)
    TENV = T0 * rqList
    TENV = np.where(TENV < T_min, T_min, TENV)

    #TENV = np.where( T0 * rqList < 1.0, T_min, T0 * rqList )
    
    #------------------------
    #------------------------

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'total': TENV, 'disc': np.zeros(NPoints), 'env': TENV, 
                      'discFlag': False, 'envFlag': True
                      } ) 

#---------------------------
#---------------------------

#----------------------------------------
#TEMPERATURE (PowerLaw-standard) FUNCTION
#-----------------------------------------

def temperature_Powerlaw2(r_max, r_min, T0, q, GRID, T_min = 2.725):

#r_max: Maximum radius of the envelope 
#r_min: Minimum radius of the envelope 
#T0: Temperature at r_min
#q: power-law for temperature
#GRID: Grid to work in
#T_min: Minimum temperature

    #------------
    #LISTS TO USE
    #------------
    rList, NPoints = GRID.rRTP[0], GRID.NPoints #Due to spherical symmetry only r is needed

    #------------------------
    #MODEL. Envelope powerlaw
    #------------------------
    print ('Computing Envelope temperature using power-law...')
    rqList = np.where((rList >= r_min) & (rList <= r_max) , rList**q, 0.)
    TENV = T0 * rqList #r_min**-q * rqList
    TENV = np.where(TENV < T_min, T_min, TENV)

    #------------------------
    #------------------------

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'total': TENV, 'disc': np.zeros(NPoints), 'env': TENV, 
                      'discFlag': False, 'envFlag': True
                      } ) 

#---------------------------
#---------------------------

#------------------------------------
#TEMPERATURE (PowerLaw-Shells) FUNCTION
#------------------------------------

def temperature_PowerlawShells(r_list, p_list, T0, GRID, T_min = 1.0e3):

#r_list: List of shells' limits, length (n,)
#p_list: List of powerlaws in r_list intervals, length (n-1,)
#T0: Temperature at r_list[0]
#GRID: Grid to work in
#T_min: Minimum temperature

    #------------
    #LISTS TO USE
    #------------
    rList, NPoints = GRID.rRTP[0], GRID.NPoints #Due to spherical symmetry only r is needed

    #-------------------------
    #MODEL. Envelope powerlaws
    #-------------------------
    print ('Computing Envelope temperature using power-laws:', p_list)

    TENV = np.zeros(NPoints)
    T0_list = [T0]
    for i,r in enumerate(r_list[1:-1],1):
        r_tmp = np.max(rList[rList<=r])
        T0_list.append(T0_list[i-1]*(r_tmp/r_list[i-1])**p_list[i-1])

    for i,p in enumerate(p_list):
        ind, = np.where((rList >= r_list[i]) & (rList <= r_list[i+1]))
        TENV[ind] += T0_list[i]*(rList[ind]/r_list[i])**p  

    TENV = np.where(TENV < T_min, T_min, TENV)

    #------------------------
    #------------------------

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'total': TENV, 'disc': np.zeros(NPoints), 'env': TENV, 
                      'discFlag': False, 'envFlag': True
                      } ) 

#---------------------------
#---------------------------

#----------------------
#VELOCITY FUNCTION
#----------------------

def velocity(RStar,MStar,Rd,density,GRID):

#MStar: Mass of the central source
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

    vr = np.zeros(NPoints)
    vtheta = np.zeros(NPoints)
    vphi = np.zeros(NPoints)
    #------------------------------
    #MODEL. Keto,E & Zhang,Q (2010)
    #------------------------------

    #------------
    #DISC Profile
    #------------
    if density.discFlag:
        print ('Computing Disc velocity...')
        rdisc = density.r_disc
        #Pure azimuthal component. It's assumed that the radial velocity in the rotationally supported disc is comparatively small (Keto 2010).
        vdisc = np.where( RList <= rdisc, (G * MStar / RList)**0.5, 0*-3e8)
        vphi = vdisc
    else: vdisc = 0.

    #----------------
    #ENVELOPE PROFILE
    #----------------
    if density.envFlag:
        print ('Computing Envelope velocity...')
    
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
        good_ind = np.where( (thetaList != 0.) & (rList <= renv) ) #To avoid the polar axis for theta and phi. Along the polar axis all the velocities are all radial.
        vtheta, vphi = np.zeros(NPoints), np.zeros(NPoints)

        vtheta[good_ind] = signo[good_ind] * ( (G * MStar / rList[good_ind])**0.5 * (costheta0[good_ind] - costheta[good_ind]) / sintheta[good_ind] * 
                                               (1 + costheta[good_ind] / costheta0[good_ind])**0.5 )
        vphi[good_ind] = (G * MStar / rList[good_ind])**0.5 * (sintheta0[good_ind] / sintheta[good_ind]) * (1 - costheta[good_ind] / costheta0[good_ind])**0.5


    #Weighted velocity by density. Vectorial sum.
    if density.envFlag and density.discFlag: 
        vr , vtheta, vphi = map( lambda x,y : (rhoDISC * x + rhoENV * y) / density.total, 
                                 [ 0, 0, vdisc ], [vr, vtheta, vphi] ) 

    print ('Converting to cartesian coordinates...') 
    vx, vy, vz = sphe_cart( list( zip(vr, vtheta, vphi) ), theta4vel, phiList)

    for vi in [vx,vy,vz]: vi[GRID.r_ind_zero] = 0.0
    #for vi in [vx,vy]: vi[GRID.R_ind_zero] = 0.0

    #----------------------------------------------
    #----------------------------------------------
    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'x': vx, 'y': vy, 'z': vz} )


#----------------------
#----------------------

def velocity_piecewise(density, GRID, 
                       R_list=None, pR_list=None, v0R=None, #Polar piecewise
                       r_list=None, pr_list=None, v0r=None  #Radial piecewise
                       ):

#density grid
#GRID: grid to work with
#R_list: List of polar piecewise limits, length (n,)
#pR_list: List of powerlaws in R_list intervals, length (n-1,)
#v0R: velocity vector at R_list[0] --> polar rotation, in sph. coords: v = (0,0,vphi(R))

#r_list: List of radial piecewise limits, length (n,)
#pr_list: List of powerlaws in r_list intervals, length (n-1,)
#v0r: velocity vector at r_list[0], can be negative or positive --> radial infall/outflow, defaults to positive (outflow). In sph. coords: v = (vr(r),0,0)

    #------------
    #LISTS TO USE
    #------------
    rList, RList, phiList = GRID.rRTP[0], GRID.rRTP[1], GRID.rRTP[3]
    theta4vel = GRID.theta4vel
    NPoints = GRID.NPoints 

    #-------------------------------------------
    #MODEL. polar and radial piecewise powerlaws
    #-------------------------------------------
    def v_kepler(r, p):
        return (G*MStar/r)**p 

    vr, vtheta, vphi = np.zeros((3,NPoints))        
    if R_list is not None:
        print ('Computing phi velocity with powerlaws', pR_list)
        vphi_coeff = [v0R[2]]
        #Pure azimuthal component. It's assumed that the radial velocity in the rotationally supported disc is comparatively small (Keto 2010).
        for i,R in enumerate(R_list[1:-1],1):
            R_tmp = np.max(RList[RList<=R])
            vphi_coeff.append(vphi_coeff[i-1]*(R_tmp/R_list[i-1])**pR_list[i-1])
    
        for i,p in enumerate(pR_list):
            ind, = np.where((RList >= R_list[i]) & (RList <= R_list[i+1]))
            vphi[ind] += vphi_coeff[i]*(RList[ind]/R_list[i])**p  
        
    if r_list is not None:
        print ('Computing infall velocity with powerlaws', pr_list)
        vr_coeff = [v0r[0]]
        #Infall velocity. Pointing radially inwards.
        for i,r in enumerate(r_list[1:-1],1):
            r_tmp = np.max(rList[rList<=r])
            vr_coeff.append(vr_coeff[i-1]*(r_tmp/r_list[i-1])**pr_list[i-1])
    
        for i,p in enumerate(pr_list):
            ind, = np.where((rList >= r_list[i]) & (rList <= r_list[i+1]))
            vr[ind] += vr_coeff[i]*(rList[ind]/r_list[i])**p  

    print ('Converting to cartesian coordinates...') 
    
    if density.discFlag and density.envFlag:
        vphi[RList>density.r_disc] = 0.0
        vr[rList>density.r_env] = 0.0
    elif density.discFlag: 
        vphi[RList>density.r_disc] = 0.0
        vr[RList>density.r_disc] = 0.0
    elif density.envFlag: 
        vphi[rList>density.r_env] = 0.0
        vr[rList>density.r_env] = 0.0

    vx, vy, vz = sphe_cart( list( zip(vr, vtheta, vphi) ), theta4vel, phiList)

    
    for vi in [vx,vy,vz]: vi[GRID.r_ind_zero] = 0.0        
    #------------------------
    #------------------------

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'x': vx, 'y': vy, 'z': vz} )



#--------------------------
#VELOCITY (Random) FUNCTION
#--------------------------

def velocity_random(v_disp,NPoints):

    print ('Computing random (uniform) velocities...')
    v_disp = v_disp/np.sqrt(3)    

    v_x = v_disp * (2 * np.random.random(NPoints) - 1)  
    v_y = v_disp * (2 * np.random.random(NPoints) - 1)  
    v_z = v_disp * (2 * np.random.random(NPoints) - 1)  

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'x': v_x, 'y': v_y, 'z': v_z} )

#-------------------------------
#VELOCITY-INFALL D.W.Murray+2017
#http://adsabs.harvard.edu/abs/2017MNRAS.465.1316M
#See eq. (4)
#-------------------------------

def velocity_infall(dens_dict, ff_factor, MStar, r_stellar, GRID, v0 = [0.,0.,0.]):

#dens_dict: dictionary containing the density distributions to compute the enclosed mass from
#ff_factor: free-fall normalization factor
#MStar: mass of the star
#r_stellar: stellar sphere of influence
#GRID
#v0: systemic velocity, optional

    from .utils.prop import propTags

    rList = GRID.rRTP[0]
    dx,dy,dz = [GRID.XYZgrid[i][1] - GRID.XYZgrid[i][0] for i in range(3)]

    print ('Computing infall velocities (D.W.Murray+2017)...')
    rUnique = np.unique(rList) #sorted unique values of r

    mass_unit = np.array([propTags.get_dens_mass(dens_name) for dens_name in dens_dict])
    mass = np.sum([dens * mass_unit[i] for i,dens in enumerate(dens_dict.values())], axis=0) * dx*dy*dz
    
    speed = np.zeros(GRID.NPoints)

    r_at_stellar = rUnique[rUnique < r_stellar][-1] #Closest r to r_stellar from the left
    speed_r_stellar = np.sqrt(G*MStar/r_at_stellar)
    speed_menc_stellar = np.sqrt(G*np.sum(mass[rList < r_stellar])/r_at_stellar)
    norm_factor = speed_r_stellar / speed_menc_stellar #Normalization factor to connect both profiles
    
    for r in rUnique: 
        ind_enc = rList <= r
        ind_r = rList == r
        if r < r_stellar: 
            speed_r = np.sqrt(G*MStar/r)
        else:
            mass_enc = np.sum(mass[ind_enc])
            speed_r = norm_factor * np.sqrt(G*mass_enc/r)
        speed[ind_r] = speed_r

    foo = -ff_factor*speed / GRID.rRTP[0]
    v_x = foo*GRID.XYZ[0]
    v_y = foo*GRID.XYZ[1]
    v_z = foo*GRID.XYZ[2]

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

    return Struct( **{'x': v_x, 'y': v_y, 'z': v_z} )


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
    
    print ('Computing profiles for a hole...')
    
    densNew = np.where( (tempList >= T_min) & (tempList <= T_max), dens_val, densList)
    tempNew = np.where( (tempList >= T_min) & (tempList <= T_max), temp_val, tempList)
    abundNew = np.where( (tempList >= T_min) & (tempList <= T_max), abund_val, abundList)

    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')
    
    return Struct( **{'dens': densNew, 'temp': tempNew, 'abund': abundNew})


#--------------------------------
#--------------------------------    

def PrintProperties(density, temperature, GRID, species='dens_H2'): 

    from .utils.prop import propTags

    dv = GRID.step[0]*GRID.step[1]*GRID.step[2]
    inddisc = np.where(temperature.disc > 2.)
    indtotal = np.where(temperature.total > 2.)
    Mu_MSun = propTags.dens_mass[species] / MSun #2 * Mu/MSun
    
    print ('Mass from '+species)
    print ('Total mass (MSun):', np.sum(density.total) * dv * Mu_MSun)
    print ('Mean Total Temperature (Kelvin), weighted by density:', 
           (np.sum(temperature.total[ indtotal ] * density.total[ indtotal ]) 
            / np.sum(density.total[ indtotal ]) ))
    if density.discFlag:
        print ('Total Disc mass (MSun):', np.sum(density.disc) * dv * Mu_MSun)
        print ('Total Envelope mass (MSun):', np.sum(density.env) * dv * Mu_MSun)
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
    Rot_iter = iter(Rot[1:]) #Iterator for Rot_list from 2nd value: (matrix for matrix in Rot_list[1:])

    for i in range( len(Rot[1:]) ): 
        tmp = np.dot( next(Rot_iter) , tmp )
        
    Rot_total = tmp

    print ('%s is done!'%inspect.stack()[0][3])
   
    return Rot_total


def ChangeGeometry(GRID, center = False ,rot_dict = False, vel = False, vsys = False):

#order: indicates the order of the axes along which the rotations will be performed
 
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
        print ('Computing Rotation matrix...')
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
            print ('==========================================================') 
            print ('WARNING: No VELOCITY distribution was provided to be rotated!')
            print ('You should provide it if interested in computing line emission.')
            print ('==========================================================')

    if center is not False:
        print ('Moving the grid to the new center...')
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

#-----------------------------
#WRITING DATA (RADMC-3D v0.41) [OLD Function]
#-----------------------------

def DataTab_RADMC3D_FreeFree(dens,temp,GRID):

    #dens = 1e-6 * np.where(dens > 10.0, dens, 0) #to cm^-3
    dens = dens / 1e6
    nx,ny,nz = GRID.Nodes
    xi, yi, zi = np.array(GRID.XYZgrid) * 100 #to cm
    nphot = 1000000
#
# Writing the grid file
#
    with open('amr_grid.inp','w+') as f:
        f.write('1\n')                       # iformat
        f.write('0\n')                       # AMR grid style  (0=regular grid, no AMR)
        f.write('0\n')                       # Coordinate system
        f.write('0\n')                       # gridinfo
        f.write('1 1 1\n')                   # Include x,y,z coordinate
        f.write('%d %d %d\n'%(nx,ny,nz))     # Size of grid
        tmp = ['%13.6e '*(n+1) for n in [nx,ny,nz]]
        f.write((tmp[0]+'\n')%tuple(xi))
        f.write((tmp[1]+'\n')%tuple(yi))
        f.write((tmp[2]+'\n')%tuple(zi))
        
#
# Writing the electronic density file.
#
    with open('electron_numdens.inp','w+') as f:
        f.write('1\n')                       # Format number
        f.write('%d\n'%(nx*ny*nz))           # Nr of cells
        #data = rhoelect.ravel(order='F')         # Create a 1-D view, fortran-style indexing
        dens.tofile(f, sep='\n', format="%13.6e")
        f.write('\n')

#
# Writing the ion density file.
#
    with open('ion_numdens.inp','w+') as f:
        f.write('1\n')                       # Format number
        f.write('%d\n'%(nx*ny*nz))           # Nr of cells
        #data = rhoelect.ravel(order='F')         # Create a 1-D view, fortran-style indexing
        dens.tofile(f, sep='\n', format="%13.6e")
        f.write('\n')
    
#
# Writing the gas temperature
#

    with open('gas_temperature.inp','w+') as f:
        f.write('1\n')                       # Format number
        f.write('%d\n'%(nx*ny*nz))           # Nr of cells
        #data = temp_gas.ravel(order='F')          # Create a 1-D view, fortran-style indexing
        temp.tofile(f, sep='\n', format="%13.6e")
        f.write('\n')

#
# Writing the wavelength_micron.inp file
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
# Writing the wavelength file
#
    with open('wavelength_micron.inp','w+') as f:
        f.write('%d\n'%(nlam))
        for value in lam:
            f.write('%13.6e\n'%(value))

#
# Writing the radmc3d.inp control file
#
    with open('radmc3d.inp','w+') as f:
        f.write('nphot = %d\n'%(nphot))
        f.write('scattering_mode_max = 1\n')   #Set it to 1 for isotropic scattering
        f.write('incl_freefree = 1\n')
        f.write('incl_dust = 0\n')
        f.write('setthreads = 4\n')
        #f.write('tgas_eq_tdust = 1')
        
    
    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

#-------------
#-------------
def col_ids_LIME(props):
    
    CP_H2 = 1
    CP_p_H2 = 2
    CP_o_H2 = 3
    CP_e = 4
    CP_H = 5
    CP_He = 6
    CP_Hplus = 7
    base = dict(SF3D_id =            0,
                SF3D_x=              1,
                SF3D_y =             2,
                SF3D_z =             3,
                SF3D_dens_H2 =       CP_H2 + 3,     
                SF3D_dens_p_H2 =     CP_p_H2 + 3,
                SF3D_dens_o_H2 =     CP_o_H2 + 3,
                SF3D_dens_e =        CP_e + 3,
                SF3D_dens_H =        CP_H + 3,
                SF3D_dens_He =       CP_He + 3,
                SF3D_dens_Hplus =    CP_Hplus + 3,
                SF3D_temp_gas =      11,
                SF3D_temp_dust =     12,
                SF3D_vel_x =         13,
                SF3D_vel_y =         14,
                SF3D_vel_z =         15,
                SF3D_abundance =     16,
                SF3D_gtdratio =      17,
                SF3D_doppler =       18,
                SF3D_max_cols =      19)
    written_props = [] 
    for col in props:
        if col in base.keys(): written_props.append(base[col])
    written_props.append(4242)
    return written_props


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
        sfile.write("%d %d %d %d %d"%(8,Ns[0],Ns[1],Ns[2],GRID.NPoints))
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
        
        colsfile = './header.dat'
        props = ['SF3D_id', 'SF3D_dens_H2', 'SF3D_temp_gas', 
                 'SF3D_vel_x', 'SF3D_vel_y', 'SF3D_vel_z',
                 'SF3D_abundance', 'SF3D_gtdratio']
        cols2write = np.array([col_ids_LIME(props)]).T
        print ('Writing column ids on %s'%colsfile)
        np.savetxt(colsfile, cols2write, fmt='%d')
        
    file.close()
    
    print ('%s is done!'%inspect.stack()[0][3])
    print ('-------------------------------------------------\n-------------------------------------------------')

#--------------
#--------------    

def DataTab_LIME2(dens_H2,dens_H,dens_Hp,temp,vel,abund,gtd,GRID,tdust = None,doppler=None, is_submodel = False, tag = False, fixed_grid = False):
    
    import pandas

    if is_submodel:
        os.system('mkdir Subgrids')
        file0 = './Subgrids/datatab%s.dat'%tag
        file = open(file0,'w')
        x,y,z = GRID.XYZ 
    
        print ('Writing Submodel data on %s'%file0)
        tmp = []
        colsfile = './Subgrids/header.dat'
        props = 0
        if fixed_grid:
            
            """
            Sorting from least to greatest leads to confusions in the reorderGrid function of Lime (at grid.c and grid_aux.c)
            because the last pIntensity points (the grid points written here) might likely be flagged as sinkPoints since they lie
            over (or close) the domain's radius and therefore their ids (close to Npoints from the left) lie very close to those of the 
            sinkPoints defined first (close to Npoints from the right). This fact, causes almost always an underestimation of the 
            new sinkPoints found by Lime (nExtraSinks), and required to recalculate the parameters par->pIntensity and par->sinkPoints. 
            
            Thus, sorting from greatest to least solves this issue. 
            """
            rr = np.linalg.norm(GRID.XYZ, axis = 0)
            ind_rr = iter(np.argsort(rr)[::-1])

            id = 0       
            
            if isinstance(tdust,list) or isinstance(tdust,np.ndarray):
                for i in ind_rr: 
                    tmp.append( "%d %e %e %e %e %e %e %e %e %e %e %e %e %e\n"% 
                                (id,x[i],y[i],z[i],dens_H2[i],dens_H[i],dens_Hp[i],temp[i],
                                 tdust[i],vel.x[i],vel.y[i],vel.z[i],abund[i],gtd[i]))
                    id+=1
        
                props = ['SF3D_id', 'SF3D_x', 'SF3D_y', 'SF3D_z',
                         'SF3D_dens_H2', 'SF3D_dens_H', 'SF3D_dens_Hplus',
                         'SF3D_temp_gas', 'SF3D_temp_dust',
                         'SF3D_vel_x', 'SF3D_vel_y', 'SF3D_vel_z',
                         'SF3D_abundance', 'SF3D_gtdratio']
        
            else: 
                for i in ind_rr: 
                    tmp.append( "%d %e %e %e %e %e %e %e %e %e %e %e %e\n"% 
                                (id,x[i],y[i],z[i],dens_H2[i],dens_H[i],dens_Hp[i],temp[i],
                                 vel.x[i],vel.y[i],vel.z[i],abund[i],gtd[i]))
                    id+=1
        else:

            if ((isinstance(doppler,list) or isinstance(doppler,np.ndarray)) and 
                (isinstance(tdust,list) or isinstance(tdust,np.ndarray))):
                for i in range(GRID.NPoints):
                    tmp.append( "%d %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n"% 
                                (i,x[i],y[i],z[i],dens_H2[i],dens_H[i],dens_Hp[i],temp[i],
                                 tdust[i],vel.x[i],vel.y[i],vel.z[i],abund[i],gtd[i],doppler[i]))
                
                props = ['SF3D_id', 'SF3D_x', 'SF3D_y', 'SF3D_z',
                         'SF3D_dens_H2', 'SF3D_dens_H', 'SF3D_dens_Hplus',
                         'SF3D_temp_gas', 'SF3D_temp_dust',
                         'SF3D_vel_x', 'SF3D_vel_y', 'SF3D_vel_z',
                         'SF3D_abundance', 'SF3D_gtdratio', 'SF3D_doppler']
           
            elif isinstance(tdust,list) or isinstance(tdust,np.ndarray):
                for i in range(GRID.NPoints): 
                    tmp.append( "%d %e %e %e %e %e %e %e %e %e %e %e %e %e\n"% 
                                (i,x[i],y[i],z[i],dens_H2[i],dens_H[i],dens_Hp[i],temp[i],
                                 tdust[i],vel.x[i],vel.y[i],vel.z[i],abund[i],gtd[i]))
                props = ['SF3D_id', 'SF3D_x', 'SF3D_y', 'SF3D_z',
                         'SF3D_dens_H2', 'SF3D_dens_H', 'SF3D_dens_Hplus',
                         'SF3D_temp_gas', 'SF3D_temp_dust',
                         'SF3D_vel_x', 'SF3D_vel_y', 'SF3D_vel_z',
                         'SF3D_abundance', 'SF3D_gtdratio']

            else: 
                for i in range(GRID.NPoints): tmp.append( "%d %e %e %e %e %e %e %e %e %e %e %e %e\n"% 
                                                           (i,x[i],y[i],z[i],dens_H2[i],dens_H[i],dens_Hp[i],temp[i],
                                                            vel.x[i],vel.y[i],vel.z[i],abund[i],gtd[i]))
        
        file.writelines(tmp)
        print (props)
        cols2write = np.array([col_ids_LIME(props)]).T
        print ('Writing column ids on %s'%colsfile)
        np.savetxt(colsfile, cols2write, fmt='%d')

        sizefile='./Subgrids/npoints.dat'
        print ('Writing grid size on %s'%sizefile)
        sfile = open(sizefile,'w') 
        Ns = [int(GRID.NPoints**(1/3.))]*3
        sfile.write("%d %d %d %d %d"%(len(props),Ns[0],Ns[1],Ns[2],GRID.NPoints))

        
    else:
        files=['datatab.dat','x.dat','y.dat','z.dat']
        sizefile='./npoints.dat'
        print ('Writing grid size on %s'%sizefile)
        sfile = open(sizefile,'w') 
        Ns = GRID.Nodes
        sfile.write("%d %d %d %d %d"%(11,Ns[0],Ns[1],Ns[2],GRID.NPoints))
        print ('Writing data on %s'%files[0])
        file = open(files[0],'w')

        if isinstance(tdust,list) or isinstance(tdust,np.ndarray):
            for i in range(GRID.NPoints): 
                file.write("%d %e %e %e %e %e %e %e %e %e %e\n"%
                           (i,dens_H2[i],dens_H[i],dens_Hp[i],temp[i],tdust[i],vel.x[i],vel.y[i],vel.z[i],abund[i],gtd[i]))
        else:
            for i in range(GRID.NPoints): 
                file.write("%d %e %e %e %e %e %e %e %e %e\n"%
                           (i,dens_H2[i],dens_H[i],dens_Hp[i],temp[i],vel.x[i],vel.y[i],vel.z[i],abund[i],gtd[i]))

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
