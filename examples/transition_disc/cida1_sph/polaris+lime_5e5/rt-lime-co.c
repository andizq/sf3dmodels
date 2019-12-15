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
  par->radius                   = 100*AU;
  par->minScale                 = 0.1*AU; // 2 * sizex / Nx / 2
  par->pIntensity               = 239098;//273775;
  par->sinkPoints               = 60000; 
  par->dust                     = "../../../opacities_k05_230GHz_B_1_7.tab";
  par->moldatfile[0]            = "co.dat";
  par->moldatfile[1]            = "ch3cn.dat";
  
  par->antialias                = 1;
  par->sampling                 = 1; // if 0: log distr. for radius, directions distr. uniformly on a sphere.

  par->lte_only                 = 1;
  //par->nSolveIters              = 10;  
  par->gridfile                 = "grid.vtk";

  par->collPartIds[0]           = CP_H2; //Coll part

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

  //Dust continuum image
  i=0;
  img[i].pxls                   = 500;                // Pixels per dimension                                            
  img[i].imgres                 = 0.0029;             // Resolution in arc seconds                                        
  img[i].azimuth                = 0;                  // Azimuth                                      
  img[i].incl                   = 0.651;              // Inclination     
  img[i].posang                 = 0.21;               // Position angle                                                
  img[i].distance               = 140*PC;             // source distance in m                                                
  img[i].source_vel             = 0;                  // source velocity in m/s                                        
  img[i].unit                   = 1;                  // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau                      
  img[i].filename               = "img_cont_nodummies_band7.fits";
  img[i].freq                   = 345.7959899e9;      //Continuum central frequency 


  //Dust continuum image
  i=1;
  img[i].pxls                   = 500;                // Pixels per dimension                                            
  img[i].imgres                 = 0.0029;             // Resolution in arc seconds                                        
  img[i].azimuth                = 0;                  // Azimuth                                      
  img[i].incl                   = 0.651;              // Inclination     
  img[i].posang                 = 0.21;               // Position angle                                                
  img[i].distance               = 140*PC;             // source distance in m                                                 
  img[i].source_vel             = 0;                  // source velocity in m/s                                        
  img[i].unit                   = 4;                  // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau                      
  img[i].filename               = "img_cont_nodummies_band7_tau.fits";
  img[i].freq                   = 345.7959899e9;      //Continuum central frequency                                    


  //Dust continuum image
  i=2;
  img[i].pxls                   = 500;                // Pixels per dimension                                            
  img[i].imgres                 = 0.0029;             // Resolution in arc seconds                                        
  img[i].azimuth                = 0;                  // Azimuth                                      
  img[i].incl                   = 0.651;              // Inclination     
  img[i].posang                 = 0.21;               // Position angle                                                
  img[i].distance               = 140*PC;             // source distance in m                                                
  img[i].source_vel             = 0;                  // source velocity in m/s                                        
  img[i].unit                   = 1;                  // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau                      
  img[i].filename               = "img_cont_nodummies_band6.fits";
  img[i].freq                   = 230.0e9;            //Continuum central frequency                                                       

                                       
  //Line image
  i=3;
  img[i].trans                  = 127;                //Line transition
  img[i].molI                   = 1;                  //mol ID
  img[i].pxls                   = 500;                // Pixels per dimension                                            
  img[i].nchan                  = 81;
  img[i].velres                 = 100.;
  img[i].imgres                 = 0.0029;             // Resolution in arc seconds                                        
  img[i].azimuth                = 0;                  // Azimuth                                      
  img[i].incl                   = 0.651;              // Inclination     
  img[i].posang                 = 0.21;               // Position angle                                                
  img[i].distance               = 140*PC;             // source distance in m                                               
  img[i].source_vel             = 0;                  // source velocity in m/s                                        
  img[i].unit                   = 1;                  // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau                      
  img[i].filename               = "img_CH3CN_J19-18K4_LTE_kelvin_nodummies.fits";

  
  //Line image
  i=4;
  img[i].trans                  = 2;                  //Line transition
  img[i].molI                   = 0;                  //mol ID
  img[i].pxls                   = 500;                // Pixels per dimension                                            
  img[i].nchan                  = 81;
  img[i].velres                 = 100.;
  img[i].imgres                 = 0.0029;             // Resolution in arc seconds                                        
  img[i].azimuth                = 0;                  // Azimuth                                      
  img[i].incl                   = 0.651;              // Inclination     
  img[i].posang                 = 0.21;               // Position angle                                                
  img[i].distance               = 140*PC;             // source distance in m                                                
  img[i].source_vel             = 0;                  // source velocity in m/s                                        
  img[i].unit                   = 1;                  // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau                      
  img[i].filename               = "img_CO_J3-2_LTE_kelvin_nodummies.fits";
  
  //Line image
  i=5;
  img[i].trans                  = 2;                  //Line transition
  img[i].molI                   = 0;                  //mol ID
  img[i].pxls                   = 500;                // Pixels per dimension                                            
  img[i].nchan                  = 81;
  img[i].velres                 = 100.;
  img[i].imgres                 = 0.0029;             // Resolution in arc seconds                                        
  img[i].azimuth                = 0;                  // Azimuth                                      
  img[i].incl                   = 0.651;              // Inclination     
  img[i].posang                 = 0.21;               // Position angle                                                
  img[i].distance               = 140*PC;             // source distance in m                                                
  img[i].source_vel             = 0;                  // source velocity in m/s                                        
  img[i].unit                   = 4;                  // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau                      
  img[i].filename               = "img_CO_J3-2_LTE_tau_nodummies.fits";
  
}

/******************************************************************************/


void
density(double dummy0, double dummy1, double id, double *density){
  int id_int=round(id);
  density[0] = sf3d->dens_H2[id_int]; 
}

/******************************************************************************/

void
temperature(double dummy0, double dummy1, double id, double *temperature){
  int id_int=round(id);
  temperature[0] = sf3d->temp_gas[id_int];
  temperature[1] = sf3d->temp_dust[id_int];
}

/******************************************************************************/

void
abundance(double dummy0, double dummy1, double id, double *abundance){
  int id_int=round(id);
  abundance[0] = sf3d->abundance[0][id_int];
  abundance[1] = 10*sf3d->abundance[1][id_int];
}

/******************************************************************************/

void
gasIIdust(double dummy0, double dummy1, double id, double *gtd){
  int id_int=round(id);
  *gtd = 100;
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
  *doppler = 10.;
}

/******************************************************************************/

void
velocity(double dummy0, double dummy1, double id, double *vel){
  int id_int=round(id);
  vel[0] = sf3d->vel_x[id_int];
  vel[1] = sf3d->vel_y[id_int];
  vel[2] = sf3d->vel_z[id_int]; 
}

/******************************************************************************/
