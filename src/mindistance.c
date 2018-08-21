#ifndef _mindistance_h
#define _mindistance_h

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
