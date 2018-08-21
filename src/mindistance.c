#include "lime.h"
#include "mindistance.h"

double mindistance(double x,double *xma,int Nx){

  double mindist,distx;

  mindist=1000000*AU;
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
