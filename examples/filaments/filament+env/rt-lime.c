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
  par->radius                   = 0.3*PC; //130*PC; //200*PC; 
  par->minScale                 = 0.01*PC; // 2 * sizex / Nx / 2
  par->pIntensity               = 40000; //5*10000; //3419520; //0.25*20000; 
  par->sinkPoints               = 10000; //5*2000; //1000000; //0.25*5000; 
  par->dust                     = "../../opacities_k05_230GHz_B_1_7.tab";
  par->moldatfile[0]            = "co.dat";
  par->sampling                 = 1; // if 0: log distr. for radius, directions distr. uniformly on a sphere.

  par->lte_only = 1;
  par->gridfile                 = "grid.vtk";

  /* If required, in the next line you can set manually the normalization density (density^0.2) 
   * for Lime to accept or reject pintensity points. 0.2 is the default normalization exponent.
   * It can be modified in the pars file lime_config.h 
   */
  
  par->gridDensMaxValues[0] = 1400.;

  par->collPartIds[0]           = CP_H;  
  par->nMolWeights[0]           = 1.0; 

  /*
   * Definitions for image #0. Add blocks with successive values of i for additional images.
   */
    
  //Face-on line image
  i=0;
  img[i].nchan                  = 60;             // Number of channels
  img[i].velres                 = 60.;           // Channel resolution in m/s
  img[i].trans                  = 0;              // zero-indexed J quantum number
  img[i].pxls                   = 350;               // Pixels per dimension
  img[i].imgres                 = 0.15;           // Resolution in arc seconds
  img[i].theta                  = 0.;            // 0: face-on, pi/2: edge-on
  img[i].phi                    = 0.;            // Azimuthal angle
  img[i].distance               = 2400*PC;         // source distance in m
  img[i].source_vel             = 0;              // source velocity in m/s
  img[i].unit                   = 1;              // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[i].filename               = "img_fil+env_CO_J1-0_LTE_jypxl_4e4pInt.fits"; //Output file

  //Face-on dust continuum image
  i=1;
  img[i].freq                   = 115.27e9;          //Continuum central frequency                                     
  img[i].pxls                   = 350;               // Pixels per dimension                                            
  img[i].imgres                 = 0.15;           // Resolution in arc seconds                                        
  img[i].theta                  = 0.;            // 0: face-on, pi/2: edge-on                                      
  img[i].phi                    = 0.;            // Azimuthal angle                                                
  img[i].distance               = 2400*PC;         // source distance in m                                                
  img[i].source_vel             = 0;              // source velocity in m/s                                        
  img[i].unit                   = 1;              // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau                      
  img[i].filename               = "img_fil+env_cont115ghz_jypxl_4e4pInt.fits";
  
  //Edge-on line image
  i=2;
  img[i].nchan                  = 60;             // Number of channels
  img[i].velres                 = 60.;           // Channel resolution in m/s
  img[i].trans                  = 0;              // zero-indexed J quantum number
  img[i].pxls                   = 350;               // Pixels per dimension
  img[i].imgres                 = 0.15;           // Resolution in arc seconds
  img[i].theta                  = 0.;            // 0: face-on, pi/2: edge-on
  img[i].phi                    = M_PI/2;            // Azimuthal angle
  img[i].distance               = 2400*PC;         // source distance in m
  img[i].source_vel             = 0;              // source velocity in m/s
  img[i].unit                   = 1;              // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[i].filename               = "img_fil+env_CO_J1-0_LTE_jypxl_4e4pInt_edgeon.fits"; //Output file

  //Edge-on dust continuum image  
  i=3;
  img[i].freq                   = 115.27e9;          //Continuum central frequency                                     
  img[i].pxls                   = 350;               // Pixels per dimension                                            
  img[i].imgres                 = 0.15;           // Resolution in arc seconds                                        
  img[i].theta                  = 0.;            // 0: face-on, pi/2: edge-on                                      
  img[i].phi                    = M_PI/2;            // Azimuthal angle                                                
  img[i].distance               = 2400*PC;         // source distance in m                                                
  img[i].source_vel             = 0;              // source velocity in m/s                                        
  img[i].unit                   = 1;              // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau                      
  img[i].filename               = "img_fil+env_cont115ghz_jypxl_4e4pInt_edgeon.fits";

}

/******************************************************************************/

void
density(double dummy0, double dummy1, double id, double *density){

  int id_int;
  id_int=ceil(id);
  density[0] = sf3d->dens_H[id_int];                                           //H
}

/******************************************************************************/

void
temperature(double dummy0, double dummy1, double id, double *temperature){

  int id_int;
  id_int=ceil(id);
  temperature[0] = sf3d->temp_gas[id_int];
}

/******************************************************************************/

void
abundance(double dummy0, double dummy1, double id, double *abundance){

  int id_int;
  id_int=ceil(id);
  abundance[0] = sf3d->abundance[0][id_int];
}

/******************************************************************************/

void
gasIIdust(double dummy0, double dummy1, double id, double *gtd){

  int id_int;
  id_int=ceil(id);
  *gtd = sf3d->gtdratio[id_int];
}

/******************************************************************************/

//Required for Lime
void
doppler(double dummy0, double dummy1, double id, double *doppler){

  *doppler = 0.0;
}

/******************************************************************************/

void
velocity(double dummy0, double dummy1, double id, double *vel){

  int id_int;
  id_int=ceil(id);
  vel[0] = sf3d->vel_x[id_int];
  vel[1] = sf3d->vel_y[id_int];
  vel[2] = sf3d->vel_z[id_int];
}

/******************************************************************************/


