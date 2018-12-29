from __future__ import print_function
from .Utils import *
import numpy as np
import matplotlib.pyplot as plt

def RA_to_deg(strval):
    strval = np.array(strval.split(':'),dtype='float')
    strval = strval[0]*15+strval[1]*15/60.+strval[2]*15/3600.
    return strval

def DEC_to_deg(strval):
    strval = np.array(strval.split(':'),dtype='float')
    strval = strval[0]+np.sign(strval[0])*(strval[1]*1/60.+strval[2]*1/3600.)
    return strval

center = ['18:14:39.51','-17:52:00.00']
center = [RA_to_deg(center[0]),DEC_to_deg(center[1])]

def Real_distance(obj1,obj2,d):
    #RA and DEC in Degrees,
    # d: distance to the source in parsecs
    ang_sep = np.sqrt(((obj1[0] - obj2[0])*np.cos(np.pi * center[1] / 180))**2 + (obj1[1] - obj2[1])**2) 
    #cosine in the right ascension to scale the angle as if it were in the celestial equator 
    distance = ang_sep * 3600 * d 
    #ang_sep in arcsecs times the distance in parsecs results in the separation in AU between a pair of points   
    return distance 

def XY_to_RADEC(r,d):
    
    deltaRA = r[0]/d / np.cos(np.pi * center[1] / 180) #in arcsecs (tan(alpha) = r/d  /cos(declination)) 
    deltaDEC = r[1]/d #in arcsecs
    RA = center[0] + (-np.sign(deltaRA)) * (abs(deltaRA) / 3600)
    DEC = center[1] + deltaDEC / 3600
    
    return RA,DEC

def RADEC_to_XY(deg,d):

    RA = deg[0]
    DEC = deg[1]
    
    deltaRA = (RA - center[0]) * 3600 #to arcsecs
    deltaDEC = (DEC - center[1]) * 3600 #to arcsecs
    
    x = np.cos(np.pi * center[1] / 180) * d * - np.sign(deltaRA) * abs(deltaRA)
    y = d * np.sign(deltaDEC) * abs(deltaDEC)
    
    return x,y

"""
ridge = ['18:14:39.47776', '-17:51:59.4556']
ridge = [RA_to_deg(ridge[0]),DEC_to_deg(ridge[1])]
ridge = RADEC_to_XY(ridge,2400)
print ridge
"""
"""
MNE_Ulrich = ['18:14:39.52439', '-17:52:00.0252']
MNE_Ulrich = [RA_to_deg(MNE_Ulrich[0]),DEC_to_deg(MNE_Ulrich[1])]
MNE_Ulrich = RADEC_to_XY(MNE_Ulrich,2400)
print MNE_Ulrich
"""
"""
E_envelope = ['18:14:39.56600', '-17:52:00.0100']
E_envelope = [RA_to_deg(E_envelope[0]),DEC_to_deg(E_envelope[1])]
E_envelope = RADEC_to_XY(E_envelope,2400)
print ("E coords.:", E_envelope) 

MNE_new = ['18:14:39.51794', '-17:51:59.9715']
MNE_new = [RA_to_deg(MNE_new[0]),DEC_to_deg(MNE_new[1])]
MNE_new = RADEC_to_XY(MNE_new,2400)
print ("MNE coords.:", MNE_new)

MSW = ['18:14:39.50119', '-17:52:00.3443']
MSW = [RA_to_deg(MSW[0]),DEC_to_deg(MSW[1])]
MSW = RADEC_to_XY(MSW,2400)
print ("SNW coords.:", MSW)
"""
def Rotation_Matrix(angles):

    thX,thY,thZ = [angle for angle in angles]

    Rot = np.zeros((3,3,3))

    #Rotation along X                                 
    Rot[0] = np.array([[1,0,0],
                       [0,np.cos(thX),-np.sin(thX)],
                       [0,np.sin(thX),np.cos(thX)]])

    #Rotation along Y                                 
    Rot[1] = np.array([[np.cos(thY),0,-np.sin(thY)],
                       [0,1,0],
                       [np.sin(thY),0,np.cos(thY)]])

    #Rotation along Z                                 
    Rot[2] = np.array([[np.cos(thZ),-np.sin(thZ),0],
                       [np.sin(thZ),np.cos(thZ),0],
                       [0,0,1]])

    return Rot

def length_parabola(xlims, p):
    x0 = xlims[0] / AU
    x1 = xlims[1] / AU
    dx = 1. #AU
    p = p / AU
    x = np.arange(x0, x1, dx)
    dydx = x / (2*p)
    integral = dx * np.sum( np.sqrt(1 + dydx**2) )
    return integral

def make_parabola(p,Mass,xlims,width,drBIGGRID,angles,order,traslation,function_R_rho,function_R_T,vsys=0,vradial=0,abund=5e-8,gtd=100.,name = 'parabola_ridge.dat'):
    
    #-------------------------------------
    #ROTATION STUFF
    #-------------------------------------

    Rot = Rotation_Matrix(angles)
    Rot_list = []

    for axis in order:
        if axis == 'x': Rot_list.append(0)
        if axis == 'y': Rot_list.append(1)
        if axis == 'z': Rot_list.append(2)
    

    l,m,n = Rot_list[::-1]
    
    Rot_total = np.dot(Rot[l],np.dot(Rot[m],Rot[n]))

    #print Rot_total

    #-------------------------------------
    #-------------------------------------


    #-------------------------------------
    #POPULATING THE PARABOLA
    #-------------------------------------

    drmax = drBIGGRID
    dr = drmax/4.
    

    midpoint = np.sum(xlims) / 2 #Finding a third point in the parabola 
    xflags = np.sort(np.append(xlims, midpoint)) #Adding that point to a list with the referent points of the parabola
    yflags = xflags ** 2 / (4*p) #Calculating the y value for that list of points
    vecflags = [np.array([xflags[i],yflags[i]]) for i in range(len(xflags))] #Joining each pair of x,y in a vector
    
    #print xflags, vecflags
    
    refvecs = [vecflags[1] - vecflags[0], vecflags[1] - vecflags[2]] #Calculating the vectors between the xlims,ylims points and the midpoint
    r_seg_mag = np.sum( [np.linalg.norm(vec) for vec in refvecs] ) 
    #Calculating the magnitude of those vectors and adding it together to obtain an rough estimate ('guesstimate') of the length of the parabola
    
    
    #print (np.linalg.norm(vecflags[2]-vecflags[0])) # print the distance between the boundary points 
    
    #print (refvecs)
    #print (r_seg_mag)
    
    Npoints = int(r_seg_mag/dr * width/dr * width/dr)
    print ("Npoints:", Npoints)
    
    file = open(name,'w')
    x,y,z = (0,0,0)
    vx,vy,vz = (0,0,0)
    vel = 0
    rho = 0
    T = 0

    xList = []
    yList = []
    vmean = []
    
    for i in range(Npoints):
        #break
        xrand = np.random.uniform(xlims[0] , xlims[1]) #Random x points between the given limits.
        yrand = xrand ** 2 / (4*p)
        curve_vec = np.array([xrand,yrand,0])

        R = np.random.uniform(0,width) #Random R from the segment 
       
        theta = np.arctan2(1./(2*p),1./xrand)
        tangent_vec = np.array([np.cos(theta),np.sin(theta),0]) #Unitary tangent vector to the parabola 
        
        rand_vec = np.random.uniform(-1,1,3) #Random vector to generate the perpendicular vector to the segment
        rand_vec = rand_vec/np.linalg.norm(rand_vec)
        
        rand_vec_plane =  R * np.cross(tangent_vec,rand_vec) #Perpendicular (random) vector to the segment (normal)
  
        r_real = curve_vec + rand_vec_plane #r vector of the point
    
        dist = np.linalg.norm(np.array([0,p,0]) - r_real) #Distance between the center of mass and the point 
        vel = np.sqrt(2*G*Mass/dist)

        #print dist/AU, vel
        
        vx,vy,vz = vel * - np.sign(tangent_vec[0]) * tangent_vec + vradial * - rand_vec_plane / R #Velocity of the given r_real 
        
        rho = function_R_rho(R) #Density of the given R
        T = function_R_T(R) #Temperature of the given R

        r_real = np.dot(Rot_total,r_real) #Rotated r vector
        vx,vy,vz = np.dot(Rot_total,np.array([vx,vy,vz])) #Rotated velocity vector
        vz = vz + vsys
       
        #if abs(vz) < 100: print vz
        x,y,z = r_real
        x,y,z = np.array([x,y,z]) + traslation
        
        file.write('%d %e %e %e %e %e %e %e %e %e %e\n'%(i,x,y,z,rho,T,vx,vy,vz,abund,gtd))

        xList.append(x)
        yList.append(y)
        vmean.append(vel)
        #break

    file.close()
    vmean = np.sum( np.array(vmean) ) / len(vmean)
    print ("Mean tangential velocity:", vmean, "m/s")
    cross_sec = np.pi * width**2
    inflow_rate = vmean * (rho * 2*Mu) * cross_sec / MSun * yr
    print ("Mass inflow rate:", inflow_rate, "MSun/yr")

    return np.array(xList),np.array(yList)


def make_paraboloid(pos_c, pos_f, r_min, drBIGGRID, a, b, dens, temp, width = 0, 
                    r_max = 0, vsys=0, name = 'paraboloid0.dat'):
    
    """
    pos_c: Position of the paraboloid's center
    pos_f: Position of the farthest point from the center along the paraboloid's axis   
    drBIGGRID: Maximum separation between nodes in the Global Grid

    Note: If r_max is set, pos_f will just store the vectorial direction of the paraboloid's axis: pos_f = pos_i + dir. 
    For example: pos_i = np.array([100*AU, 0, 0]), r_max = 2000 * AU, then an outflow pointing towards the Y direction should have
    pos_f = pos_i + np.array([0, 1, 0]). The equivalent solution (without setting r_max) would be pos_f = np.array([100*AU, 2000*AU, 0]).
    
    """
    #------------------------------------------
    #Paraboloid Model (Izquierdo-Ginsburg 2018)
    #------------------------------------------
    r_min = float(r_min)
    r0 = r_min

    """
    cw = w[0] * r0**-w[1]
    def width(r): #Jet half-width 
        return cw * r**w[1]
    
    cd = dens[0] * r0**-dens[1]
    def density(r): #Jet density
        return cd * r**dens[1]

    ct = temp[0] * r0**-temp[1]
    def temperature(r): #Jet temperature
        return ct * r**temp[1]

    qv = -2 * w[1] - dens[1]
    print(qv)
    cv = v0 * r0**-qv
    def velocity(r): #Jet velocity
        return cv * r**qv

    ci = ionfrac[0] * r0**-ionfrac[1]
    def ionfraction(r): #Jet ionization fraction
        return ci * r**ionfrac[1]
    """

    def ellipse_func_y(z,x): #Returns the positive Y for a given Z,X
        y = np.sqrt(b**2 * (z - (x / a)**2))
        return y    

    cd = dens[0] * r0**-dens[1]
    def density(r,y):
        return cd * r**dens[1] #* (y/(300*AU))**1.2

    a2 = a**2; b2 = b**2
    def normal_vec(x,y,z):
        vec = np.array([2*x/a2, 2*y/b2, 1])
        nor = np.linalg.norm(vec)
        dir = vec/nor
        return dir

    pos_c = np.array(pos_c)
    pos_f = np.array(pos_f)
    
    r_seg = pos_f - pos_c
    r_seg_mag = np.linalg.norm(r_seg) #Magnitude of the segment between low mass and high mass
    r_seg_dir = r_seg/r_seg_mag #Direction of the segment

    r_seg = r_seg - r_min * r_seg_dir
    
    if r_max: r_seg = (r_max - r_min) * r_seg_dir  

    r_seg_mag = np.linalg.norm(r_seg) #Magnitude of the segment between low mass and high mass
    r_seg_dir = r_seg/r_seg_mag #Direction of the segment

    #Guess of the number of points to generate: (approximation to a rectangular region)
    # Number of divisions along the main segment, 
    #  times the number of divisions perpendicular to the segment,
    #   times the number of divisions in the perpendicular to both segments above
    
    drmax = drBIGGRID
    dr = drmax/4.

    mean_w = 0.5 * r_seg_mag
    Npoints = int(r_seg_mag/dr * (mean_w/dr)**2 ) 
    print ('Number of points generated:', Npoints)
    
    flag = True
        
    file = open(name,'w')

    hwidth = 0.5*width
    x,y,z = (0,0,0)
    vx,vy,vz = (0,0,0)
    speed = 0
    vmean, nemean, Tmean = [], [], []
    rho = 0
    T = 0
    k = 0
    
    for i in range(Npoints):
        
        r = np.random.uniform(r_min, r_seg_mag) #Random z along the paraboloid's axis

        if width: 
            X = np.random.uniform(-r**0.5*a, r**0.5*a) #Random X. Must be inside the range of the projected ellipse 
            Y = ellipse_func_y(r,X) #Resultant +Y
            Y = Y * np.random.choice([-1,1]) #Random choice between + and -
            normal = normal_vec(X,Y,r) #Normal vector to the surface at x,y,z
            R = np.random.uniform(-hwidth, hwidth) #Random R from the surface
            surf_vec = R*normal 
            res_vec = np.array([X,Y,r]) + surf_vec
            x,y,z = res_vec
            
        else:
            X = np.random.uniform(-r**0.5*a, r**0.5*a) #Random X from the segment 
            Y = ellipse_func_y(r,X) #Resultant (maximum) +Y
            Y = np.random.uniform(0,Y)
            Y = Y * np.random.choice([-1,1]) #Random choice between + and -
            x,y,z = X,Y,r
            
     #   xn,yn,zn = r_real_n
            
        rr = np.sqrt(x**2 + y**2 + z**2)
        speed = 0#velocity(r)
        (vx,vy,vz) = speed * r_seg_dir #Velocity of the given r_real
        vz = vz + vsys
        
        ne = density(rr, abs(y)) #Electronic density of the given r
        T = temp #temperature(r) #Temperature of the given r

        vmean.append(speed)
        nemean.append(ne)
        Tmean.append(T)

        file.write('%d %e %e %e %e %e %e %e %e %e %e\n'%(i,x,y,z,ne,T,vx,vy,vz,0,0))
        """
        xp,yp,zp = [x,xn],[y,yn],[z,zn]
        for j in range(2): 
            k+=1
            file.write('%d %e %e %e %e %e %e %e %e %e %e\n'%(k,xp[j],yp[j],zp[j],ne,T,vx,vy,vz,0,0))
        """
        """
        while flag:
            x,y,z = np.random.uniform(0,500,3)
            #Equation of a plane with the normal unitary vector (a,b,c) in (x0,y0,z0): 
            # f = a*(x-x0) + b*(y-y0) + c*(z-z0) = 0 
        """
    file.close()
    vmean0 = np.sum( np.array(vmean) ) / len(vmean)
    netot = np.sum( np.array(nemean) )
    nemean0 = netot / len(nemean)
    Tmean0 = np.sum(np.array(nemean) * np.array(Tmean)) / netot  
    print ("Mean tangential velocity: %.2f km/s"%(vmean0 * 1e-3) )
    print ("Mean density: %.2e e-/cm^3"%(nemean0 * 1e-6) )
    print ("Mean temperature: %.1f K"%(Tmean0) )
    #cross_sec = np.pi * width(0.5 * r_seg_mag)**2
    #inflow_rate = vmean0 * (rho * 2*Mu) * cross_sec / MSun_yr
    #print ("Mass inflow rate:", inflow_rate, "MSun/yr")
    
    return r_seg_mag
