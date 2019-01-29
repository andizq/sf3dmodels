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
  par->pIntensity               = 10000; 
  par->sinkPoints               = 2000; 
  par->dust                     = "opacities_k05_230GHz_B_1_7.tab";
  par->moldatfile[0]            = "ch3cn.dat";
  //par->antialias                = 1;
  par->sampling                 = 1; // if 0: log distr. for radius, directions distr. uniformly on a sphere.
  // par->samplingAlgorithm                 = 1; //This option allows the user to build their own gridDensity (fracdensity) function to reject or accept points.

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


void
density(double dummy0, double dummy1, double id, double *density){
  int id_int=ceil(id);
  density[0] = sf3d->dens_H2[id_int]; 
}

/******************************************************************************/

void
temperature(double dummy0, double dummy1, double id, double *temperature){
  int id_int=ceil(id);
  temperature[0] = sf3d->temp_gas[id_int];
}

/******************************************************************************/

void
abundance(double dummy0, double dummy1, double id, double *abundance){
  int id_int=ceil(id);
  abundance[0] = sf3d->abundance[0][id_int];
}

/******************************************************************************/

void
doppler(double dummy0, double dummy1, double id, double *doppler){
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
velocity(double dummy0, double dummy1, double id, double *vel){
  int id_int=ceil(id);
  vel[0] = sf3d->vel_x[id_int];
  vel[1] = sf3d->vel_y[id_int];
  vel[2] = sf3d->vel_z[id_int]; 
}

/******************************************************************************/
