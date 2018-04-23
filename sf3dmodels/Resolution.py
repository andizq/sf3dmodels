from __future__ import print_function
import Utils as U
import numpy as np

def FieldView(xmax, radius, NP, d, scale=False):

    """ 
    Parameters
    ----------
    -> xmax (AU): Maximum limit of the model domain 
    -> radius (AU): Maximum limit of the Lime domain
    -> NP : Number of boxels of the model (Number of nodes minus 1)
    -> d (PC): distance to the source
    -> scale (AU): [False by default] Usage: Calc(0,0,0,d,scale = value)
                minimum scale of the grid. In case of irregular grids 
                this must be provided and xmax and NP won't have relevance.  

    Returns: (Scale (AU), Resolution (arcsecs), Number of pixels, Image size (arcsecs)). [type: tuple]


    NOTE: The resultant values are not mandatory to use but a guess to get a relatively good image 
          in a standard computing time. For instance, it is always possible to increase the Number 
          of pixels of the image (being consistent in changing also the Resolution value) to obtain 
          a better SMOOTHING in the result but sacrificing computing time.
    
    exit: press q
    """
        
#xmax: Maximum limit of the model domain (AU)
#NP: Number of divisions (Nx-1)
#radius: Maximum limit of the Lime domain (AU)
#d: distance to the source (PC)
#scale: (AU) minimum scale of the grid. In case of irregular
# grids this must be provided and xmax and NP won't have
#  relevance. 

    if scale==False:
        scale = 2*xmax/NP
            
    d = d*U.PC
    
    resolution = np.arctan(scale/d) # pixel resolution in radians
    phimax = 2*np.arctan(radius/d) #length of the Lime domain in radians  
    pxls = phimax/resolution

    resarcs = resolution * U.arcsecs
    size = resarcs * pxls
    
    print ("\nBoxel size in the model:", scale, "AU")
    print ("Resolution for the Lime image:", resarcs, "arcsecs")
    print ("Number of pixels for the Lime image:", int(round(pxls)))
    print ("Image size:","%.2f x %.2f"%((size,)*2), "arcsecs")

    #round(size,3)
    return scale,resarcs,pxls,size
    

def Rho0(Mr,Rd,MStar):
    
    #print "\n->Mr: Mass accretion rate in MSun/yr"
    #print "->Rd: Centrifugal radius in AU"
    #print "->MStar: Central star mass in MSun"
   
    rho0 = Mr / ( (2*U.Mu) * 4*np.pi * Rd**2 * (U.G*MStar/Rd)**0.5 ) 
    print ('Base density rho0: %e parts/m3, %e parts/cm3'% (rho0, rho0/100**3))
    return rho0
    
#Calculating the image size in pixels (1D)  

def Npxls(resolution,radius,d):    
    
    print ('(Lime domain radius in AU, distance in parsecs, resolution in arcsecs/pixel)')
    
    d=d*U.PC
    phimax = (2*np.arctan(radius/d))*3600*180/np.pi
        #length of the Lime domain in arcsecs  
    pxls = phimax/resolution
    return pxls

#Calculating the radius of the Lime model

def RadiusLime(ax,ay,d):

    print ('(image x-size in arcsecs, image y-size in arcsecs, distance in parsecs )')

    phimax=0
    d=d*U.PC
    
    if ax<=ay: phimax = ax
    else: phimax = ay
    #because tan(phimax/2) = radius/d
    radius = d*np.tan(phimax/2.*1./3600*np.pi/180)
    return radius #in AU

#Converting arsecs to AUs

def Arcsec_AU(Ang,d):
    
    print ('(Angle in arcsecs, distance in parsecs), resultant distance in AU')

    #Exact way:
    
    Ang_torad = Ang * (1./3600) * (np.pi/180)
    r_inpc = d * np.tan(Ang_torad)
    r_inAU = r_inpc * U.PC  #parsecs to AU
        
    #A very good aproximation (small angles):
    """
    r_inAU = Ang * d
    #That's because (1 PC in AU)^-1 = (1./3600) * (np.pi/180)
    """

    return r_inAU


#Calculating the central mass with Newton's law
# given the speed at a certain distance

def Mass(v,r):
    
    # a_c = v**2/r
    # And a_c = GM/r**2
    # Therefore, M = (r*v**2)/G
    
    print ('(velocity in m/s, distance from central mass in AU)')

    r_toSI = r * U.AU
    M = r_toSI * v**2 / U.G
    M_inMSun = M / U.MSun
    return M_inMSun
    
    
    
