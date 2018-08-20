#ifndef _readdata_h
#define _readdata_h

#include <stdio.h>
#include <stdlib.h>

extern int Ndata, Nx, Ny, Nz, *ID;
extern double *xm, *ym, *zm;
extern double *DENS, *TEMP, *VEL_x, *VEL_y, *VEL_z, *ABUND, *GTD;
extern void readDatatab();

#endif

/*
int Ndata=421875;
int Nx=75;
int Ny=75;
int Nz=75;
*/

int Ndata, Nx, Ny, Nz, *ID;
double *xm, *ym, *zm;
double *DENS, *TEMP, *VEL_x, *VEL_y, *VEL_z, *ABUND, *GTD;

void readDatatab(){

  int noo;

  FILE *gridsize = fopen("npoints.dat", "r");
  fscanf(gridsize,"%d %d %d %d",&Nx,&Ny,&Nz,&Ndata);
  printf("Nx,Ny,Nz,N: %d %d %d %d\n",Nx,Ny,Nz,Ndata);
  xm = malloc (sizeof(double) * Nx);
  FILE *fx  = fopen("x.dat", "r");
  for( noo = 0; noo < Nx; noo++ ){

    fscanf(fx,"%lf",&xm[noo]);
   
  }
 
  /////////////////////////////////////////

  ym = malloc (sizeof(double) * Ny);
  FILE *fy  = fopen("y.dat", "r");
  for( noo = 0; noo < Ny; noo++ ){

    fscanf(fy,"%lf",&ym[noo]);
   
  }

  /////////////////////////////////////////

  zm = malloc (sizeof(double) * Nz);
  FILE *fz  = fopen("z.dat", "r");
  for( noo = 0; noo < Nz; noo++ ){

    fscanf(fz,"%lf",&zm[noo]);
    
  }

  /////////////////////////////////////////

  ID = malloc (sizeof(int) * Ndata);
  DENS = malloc (sizeof(double) * Ndata);
  TEMP = malloc (sizeof(double) * Ndata);
  VEL_x = malloc (sizeof(double) * Ndata);
  VEL_y = malloc (sizeof(double) * Ndata);
  VEL_z = malloc (sizeof(double) * Ndata);
  ABUND = malloc (sizeof(double) * Ndata);
  GTD = malloc (sizeof(double) * Ndata);

  FILE *fp  = fopen("datatab.dat", "r");

  for( noo = 0; noo < Ndata; noo++ ){

    fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf",&ID[noo],&DENS[noo],&TEMP[noo],&VEL_x[noo],&VEL_y[noo],&VEL_z[noo],&ABUND[noo],&GTD[noo]);
    
  }
  
  fclose(fx);
  fclose(fy);
  fclose(fz);
  fclose(fp);
  
}

