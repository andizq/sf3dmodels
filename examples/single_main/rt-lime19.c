/*
 *  model.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015 The LIME development team
 *
 */

#include "lime.h"

/******************************************************************************/

void
input(inputPars *par, image *img){
  int i;
  /*
   * Basic parameters. See cheat sheet for details.
   */
  par->radius                   = 500*AU;
  par->minScale                 = 2.5*AU; // 2 * sizex / Nx / 2
  par->pIntensity               = 20000; 
  par->sinkPoints               = 5000; 
  par->dust                     = "opacities_k05_230GHz_B_1_7.tab";
  par->moldatfile[0]            = "ch3cn.dat";
  //par->antialias                = 1;
  par->sampling                 = 1; // if 0: log distr. for radius, directions distr. uniformly on a sphere.

  par->lte_only = 1;
  //  par->outputfile               = "populations.pop";
  //par->binoutputfile            = "restart.pop";
  par->gridfile                 = "grid.vtk";

  /* If required, in the next line you can set manually the normalization density (density^0.2) 
     for Lime to accept or reject pintensity points. 0.2 is the default normalization exponent.
     It can be modified in the pars file lime.h */
  
  //par->gridDensMaxValues[0] = 1400.; 

  /*
  par->collPartIds[0]           = CP_H2;
  par->nMolWeights[0]           = 1.0;
  par->dustWeights[0]           = 1.0;
  */

  /*
   * Definitions for image #0. Add blocks with successive values of i for additional images.
   */
  
  //Line image
  i=0;
  img[i].nchan                  = 62*1;             // Number of channels
  img[i].velres                 = 420.;           // Channel resolution in m/s
  img[i].trans                  = 127;  //K4            // zero-indexed J quantum number
  img[i].pxls                   = 50*4;               // Pixels per dimension
  img[i].imgres                 = 0.00833/4;           // Resolution in arc seconds
  img[i].theta                  = M_PI/4;            // 0: face-on, pi/2: edge-on
  img[i].phi                    = 0.;            // Azimuthal angle
  img[i].distance               = 2400*PC;         // source distance in m
  img[i].source_vel             = 0;              // source velocity in m/s
  img[i].unit                   = 1;              // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[i].filename               = "img_ch3cnK4_main19.fits"; //Output file


  //1st continuum image
  i=1;
  img[i].freq                   = 349.345849e9; //K4         //Continuum central frequency                                     
  img[i].pxls                   = 50*4;               // Pixels per dimension                                            
  img[i].imgres                 = 0.00833/4;           // Resolution in arc seconds                                        
  img[i].theta                  = M_PI/4;            // 0: face-on, pi/2: edge-on                                      
  img[i].phi                    = 0.;            // Azimuthal angle                                                
  img[i].distance               = 2400*PC;         // source distance in m                                                
  img[i].source_vel             = 0;              // source velocity in m/s                                        
  img[i].unit                   = 1;              // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau                      
  img[i].filename               = "img_contK4_main19.fits";

  //2nd continuum image
  i=2;
  img[i].freq                   = 342.76e9;          //Continuum central frequency                                     
  img[i].pxls                   = 50*4;               // Pixels per dimension                                            
  img[i].imgres                 = 0.00833/4;           // Resolution in arc seconds                                        
  img[i].theta                  = M_PI/4;            // 0: face-on, pi/2: edge-on                                      
  img[i].phi                    = 0.;            // Azimuthal angle                                                
  img[i].distance               = 2400*PC;         // source distance in m                                                
  img[i].source_vel             = 0;              // source velocity in m/s                                        
  img[i].unit                   = 1;              // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau                      
  img[i].filename               = "img_cont342ghz_main19.fits";

  //Line image (Tau map)
  i=3;
  img[i].nchan                  = 62*1;             // Number of channels
  img[i].velres                 = 420.;           // Channel resolution in m/s
  img[i].trans                  = 127;  //K4            // zero-indexed J quantum number
  img[i].pxls                   = 50*4;               // Pixels per dimension
  img[i].imgres                 = 0.00833/4;           // Resolution in arc seconds
  img[i].theta                  = M_PI/4;            // 0: face-on, pi/2: edge-on
  img[i].phi                    = 0.;            // Azimuthal angle
  img[i].distance               = 2400*PC;         // source distance in m
  img[i].source_vel             = 0;              // source velocity in m/s
  img[i].unit                   = 4;              // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[i].filename               = "img_tau_main19.fits"; //Output file


}

/******************************************************************************/

int Ndata, Nx, Ny, Nz, *ID;
double *xm, *ym, *zm;
double *DENS, *TEMP, *VEL_x, *VEL_y, *VEL_z, *ABUND, *GTD;


double mindistance(double x,double *xma,int Nx){

  double mindist,distx;

  mindist=100000*AU;
  int i,j;
  for( i = 0; i < Nx; i++){
    distx = fabs(x-xma[i]);

    if (distx<mindist){

      mindist=distx;
      j=i;

    }

  }

  return j;
}

void
density(double x, double y, double z, double *density){
  
  int i,j,k,Num;

  // printf("%e %e %e\n",x,y,z);
  i=mindistance(x,xm,Nx);
  j=mindistance(y,ym,Ny);
  k=mindistance(z,zm,Nz);

  Num=i*(Ny)*(Nz)+j*(Nz)+k;

  density[0] = DENS[Num]*1.0; 
}

/******************************************************************************/

void
temperature(double x, double y, double z, double *temperature){

  int i,j,k,Num;

  i=mindistance(x,xm,Nx);
  j=mindistance(y,ym,Ny);
  k=mindistance(z,zm,Nz);

  Num=i*(Ny)*(Nz)+j*(Nz)+k;

  temperature[0] = TEMP[Num];
}

/******************************************************************************/

void
abundance(double x, double y, double z, double *abundance){

  int i,j,k,Num;

  i=mindistance(x,xm,Nx);
  j=mindistance(y,ym,Ny);
  k=mindistance(z,zm,Nz);

  Num=i*(Ny)*(Nz)+j*(Nz)+k;
  abundance[0] = ABUND[Num];
}

/******************************************************************************/

void
doppler(double x, double y, double z, double *doppler){
  /*
   * 200 m/s as the doppler b-parameter. This
   * can be a function of (x,y,z) as well.
   * Note that *doppler is a pointer, not an array.
   * Remember the * in front of doppler.
   */
  *doppler = 200.;
}

/******************************************************************************/

void
velocity(double x, double y, double z, double *vel){
  
  int i,j,k,Num;

  i=mindistance(x,xm,Nx);
  j=mindistance(y,ym,Ny);
  k=mindistance(z,zm,Nz);
  

  Num=i*(Ny)*(Nz)+j*(Nz)+k;
  //  printf("%d %d %d %lf %lf %lf %lf %lf %lf %d\n", i,j,k,x,y,z,Num,VEL_x[Num],VEL_y[Num],VEL_z[Num]);

  vel[0] = VEL_x[Num];
  vel[1] = VEL_y[Num];
  vel[2] = VEL_z[Num];
}

/******************************************************************************/


