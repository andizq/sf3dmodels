from __future__ import print_function
from Utils import *
import numpy as np

def velocity(M,R_l,r_seg,r): #(Big mass,Low mass position,Particle position)                        
    
    R_h = R_l + r_seg #High mass position = Low mass position + The segment between both masses 
    r_rel = np.array(R_h) - np.array(r)

    #Virial theorem: (1/2)*<T> = -(3/5)*<V_total> in t -> inf.             
    # in wiki english: 1/2 and 3/5 can be neglected                        

    d = np.linalg.norm(r_rel)
    speed = np.sqrt(G*M/d)
    dir = r_rel / d
    vel = speed * dir

    return speed, vel #returns speed, velocity


def segment_between_stars(M_pair,R_pair):
    
    M_pair = np.array(M_pair)
    R_pair = np.array(R_pair)
    
    ind = [1,0]
    if M_pair[0] > M_pair[1]: ind = [0,1]
    
    # Forcing the vector of the segment to point towards the bigger mass  
    r_seg = R_pair[ind[0]] - R_pair[ind[1]] 
    
    # Returns the position of the lower mass and the vector of the segment
    return R_pair[ind[1]] , r_seg

    
def make_cylinder(M_pair,DiskSizes,R_l,r_seg,drBIGGRID,width, function_R_rho, function_R_T, vsys=0, abund=5e-8, gtd=100., name = 'cylinder0.dat'):
    
    #width_x = width[0]
    #width_y = width[1]
    
    M = 0

    if M_pair[0] > M_pair[1]: M = M_pair[0]
    else: M = M_pair[1]
        

    Rdisk_l = DiskSizes[0]
    Rdisk_H = DiskSizes[1]
    
    drmax = drBIGGRID
    dr = drmax/4.

    r_seg_mag = np.linalg.norm(r_seg) #Magnitude of the segment between low mass and high mass
    r_seg_dir = r_seg/r_seg_mag #Direction of the segment

    #Guess of the number of points to generate: (approximation to a rectangular region)
    # Number of divisions along the main segment, 
    #  times the number of divisions perpendicular to the segment,
    #   times the number of divisions in the perpendicular to both segments above
    
    Npoints = int(r_seg_mag/dr * width/dr * width/dr) 
    print Npoints
    
    flag = True
    
    file = open(name,'w')

    x,y,z = (0,0,0)
    vx,vy,vz = (0,0,0)
    speed = 0
    vmean = []
    rho = 0
    T = 0

    for i in range(Npoints):
        
        r = np.random.uniform(Rdisk_l, r_seg_mag - Rdisk_H) #Random r from low_mass disk border until high_mass disk border
        R = np.random.uniform(0, width) #Random R from the segment 

        r_vec = r * r_seg_dir #Random vector along the segment
        
        #a,b,c = r_vec
        
        rand_vec = np.random.uniform(-1,1,3) #Random vector to generate the perpendicular vector to the segment
        rand_vec = rand_vec/np.linalg.norm(rand_vec)
        
        rand_vec_plane =  R * np.cross(r_seg_dir,rand_vec) #Perpendicular (random) vector to the segment
        
        Rl_r = r_vec + rand_vec_plane #Vector from low mass to the generated point in the cylinder 
        r_real = R_l + Rl_r  #Real position from the origin of coordinates 
        
        x,y,z = r_real
        speed, (vx,vy,vz) = velocity(M,R_l,r_seg,r_real) #Velocity of the given r_real
        vmean.append(speed)
        
        vz = vz + vsys

        rho = function_R_rho(R) #Density of the given R
        T = function_R_T(R) #Temperature of the given R

        file.write('%d %e %e %e %e %e %e %e %e %e %e\n'%(i,x,y,z,rho,T,vx,vy,vz,abund,gtd))
        
        """
        while flag:
            x,y,z = np.random.uniform(0,500,3)
            #Equation of a plane with the normal unitary vector (a,b,c) in (x0,y0,z0): 
            # f = a*(x-x0) + b*(y-y0) + c*(z-z0) = 0 
        """
    file.close()
    vmean = np.sum( np.array(vmean) ) / len(vmean)
    print ("Mean tangential velocity:", vmean, "m/s")
    cross_sec = np.pi * width**2
    inflow_rate = vmean * (rho * 2*Mu) * cross_sec / MSun * yr
    print ("Mass inflow rate:", inflow_rate, "MSun/yr")
    
    return r_seg_mag
    
    












