#ifndef READDATA_H
#define READDATA_H

#include <stdio.h>
#include <stdlib.h>

int Ndata, Nx, Ny, Nz, *ID;
double *xm, *ym, *zm;
double *DENS, *TEMP, *VEL_x, *VEL_y, *VEL_z, *ABUND, *GTD;
void readDatatab();

#endif
