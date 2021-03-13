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
   * Basic parameters. See LIME cheat sheet for details.
   */
  par->radius                   = 110*PC; 
  par->minScale                 = 0.5*PC; 
  par->pIntensity               = 145792; 
  par->sinkPoints               = 35000; 
  //par->dust                     = "opacities_k05_230GHz_B_1_7.tab";
  par->moldatfile[0]            = "co.dat";
  //par->sampling                 = 1; // if 0: log distr. for radius, directions distr. uniformly on a sphere.

  par->lte_only = 1; // make this 1 for LTE
  par->nSolveIters              = 50;
  par->gridfile                 = "grid.vtk";

  /* Collisional partners and H-species (and e-,He) whose density will be inputted */

  par->collPartIds[0]           = CP_H;
  par->collPartIds[1]           = CP_p_H2; //Coll part
  par->collPartIds[2]           = CP_o_H2; //Coll part
  par->collPartIds[3]           = CP_Hplus;
  par->collPartIds[4]           = CP_He;
    
  /* Species weights to compute final molecular densities. 
   */
  /* --> nCO = xCO * sum(w0*nH + w1*npH2 + w2*noH2) 
   */
  
  par->nMolWeights[0]           = 1.0; 
  par->nMolWeights[1]           = 2.0;
  par->nMolWeights[2]           = 2.0;
  par->nMolWeights[3]           = 1.0;
  par->nMolWeights[4]           = 0; //0 for helium. The CO abundance only depends on H nucleons
  
  /*
   * Definitions for image #0. Add blocks with successive values of i for additional images.
   */

  //Face-on line image
  i=0;
  img[i].nchan                  = 101;             // Number of channels
  img[i].velres                 = 160.;           // Channel resolution in m/s
  img[i].trans                  = 0;              // zero-indexed J quantum number
  img[i].imgres                 = 26;           // Resolution in arc seconds
  img[i].distance               = 2400*PC;         // source distance in m
  img[i].pxls                   = ceil((2*par->radius/AU)/(img[i].distance/PC)/img[i].imgres);
  img[i].theta                  = 0.;            // 0: face-on, pi/2: edge-on
  img[i].phi                    = 0.;            // Azimuthal angle
  img[i].source_vel             = 0;              // source velocity in m/s
  img[i].unit                   = 1;              // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[i].filename               = "img_CO_J1-0_nonLTE_jypxl_faceon.fits"; //Output file
  
  //Edge-on line image
  i=1;
  img[i].nchan                  = 101;             // Number of channels
  img[i].velres                 = 120.;           // Channel resolution in m/s
  img[i].trans                  = 0;              // zero-indexed J quantum number
  img[i].imgres                 = 26;           // Resolution in arc seconds
  img[i].distance               = 2400*PC;         // source distance in m
  img[i].pxls                   = ceil((2*par->radius/AU)/(img[i].distance/PC)/img[i].imgres);
  img[i].theta                  = -M_PI/2.;            // 0: face-on, pi/2: edge-on
  img[i].phi                    = 0.;            // Azimuthal angle
  img[i].source_vel             = 0;              // source velocity in m/s
  img[i].unit                   = 1;              // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[i].filename               = "img_CO_J1-0_nonLTE_jypxl_edgeon.fits"; //Output file

  //Edge-on phi=90 line image
  i=2;
  img[i].nchan                  = 101;             // Number of channels
  img[i].velres                 = 140.;           // Channel resolution in m/s
  img[i].trans                  = 0;              // zero-indexed J quantum number
  img[i].imgres                 = 26;           // Resolution in arc seconds
  img[i].distance               = 2400*PC;         // source distance in m
  img[i].pxls                   = ceil((2*par->radius/AU)/(img[i].distance/PC)/img[i].imgres);
  img[i].theta                  = -M_PI/2.;            // 0: face-on, pi/2: edge-on
  img[i].phi                    = -M_PI/2.;            // Azimuthal angle
  img[i].source_vel             = 0;              // source velocity in m/s
  img[i].unit                   = 1;              // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[i].filename               = "img_CO_J1-0_nonLTE_jypxl_edgeon_phi90.fits"; //Output file

  //Edge-on tau image
  i=3;
  img[i].nchan                  = 101;             // Number of channels
  img[i].velres                 = 140.;           // Channel resolution in m/s
  img[i].trans                  = 0;              // zero-indexed J quantum number
  img[i].imgres                 = 26;           // Resolution in arc seconds
  img[i].distance               = 2400*PC;         // source distance in m
  img[i].pxls                   = ceil((2*par->radius/AU)/(img[i].distance/PC)/img[i].imgres);
  img[i].theta                  = -M_PI/2.;            // 0: face-on, pi/2: edge-on
  img[i].phi                    = -M_PI/2.;            // Azimuthal angle
  img[i].source_vel             = 0;              // source velocity in m/s
  img[i].unit                   = 4;              // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[i].filename               = "img_CO_J1-0_nonLTE_tau_edgeon_phi90.fits"; //Output file
}

/******************************************************************************/

void
density(double dummy0, double dummy1, double id, double *density){

  int id_int;
  id_int=round(id);
  density[0] = sf3d->dens_H[id_int];                                           //H
  density[1] = 0.5 * sf3d->dens_H2[id_int];                                    //parallel_H2
  density[2] = 0.5 * sf3d->dens_H2[id_int];                                    //orthogonal_H2 
  density[3] = sf3d->dens_Hplus[id_int];                                       //H+ or HII
  density[4] = 0.1 * (density[0] + 2*(density[1] + density[2]) + density[3]);  //Helium 
}

/******************************************************************************/

void
temperature(double dummy0, double dummy1, double id, double *temperature){

  int id_int;
  id_int=round(id);
  temperature[0] = sf3d->temp_gas[id_int];
  temperature[1] = sf3d->temp_dust[id_int]; //Not using it
}

/******************************************************************************/

void
abundance(double dummy0, double dummy1, double id, double *abundance){

  int id_int;
  id_int=round(id);
  abundance[0] = sf3d->abundance[0][id_int];
}

/******************************************************************************/

void
doppler(double dummy0, double dummy1, double id, double *doppler){

  int id_int;
  id_int=round(id);
  *doppler = sf3d->doppler[id_int];
}

/******************************************************************************/

void
velocity(double dummy0, double dummy1, double id, double *vel){

  int id_int;
  id_int=round(id);
  vel[0] = sf3d->vel_x[id_int];
  vel[1] = sf3d->vel_y[id_int];
  vel[2] = sf3d->vel_z[id_int];
}

/******************************************************************************/


