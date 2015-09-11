/*
 *  raytrace.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015 The LIME development team
 *
 */

#include "lime.h"


void
velocityspline2(double x[3], double dx[3], double ds, double binv, double deltav, double *vfac){
  int i,steps=10;
  double v,d,val,vel[3];

  *vfac=0.;
  for(i=0;i<steps;i++){
    d=i*ds/steps;
    velocity(x[0]+(dx[0]*d),x[1]+(dx[1]*d),x[2]+(dx[2]*d),vel);
    v=deltav+veloproject(dx,vel);
    val=fabs(v)*binv;
    if(val <=  2500.){
#ifdef FASTEXP
      *vfac+= FastExp(val*val);
#else
      *vfac+=   exp(-(val*val));
#endif
    }
  }
  *vfac=*vfac/steps;
  return;
}


void
line_plane_intersect(struct grid *g, double *ds, int posn, int *nposn, double *dx, double *x, double cutoff){
  double newdist, numerator, denominator ;
  int i;

  for(i=0;i<g[posn].numNeigh;i++) {
    /* Find the shortest distance between (x,y,z) and any of the posn Voronoi faces */
    /* ds=(p0-l0) dot n / l dot n */

    numerator=((g[posn].x[0]+g[posn].dir[i].x[0]/2. - x[0]) * g[posn].dir[i].x[0]+
               (g[posn].x[1]+g[posn].dir[i].x[1]/2. - x[1]) * g[posn].dir[i].x[1]+
               (g[posn].x[2]+g[posn].dir[i].x[2]/2. - x[2]) * g[posn].dir[i].x[2]);

    denominator=(dx[0]*g[posn].dir[i].x[0]+dx[1]*g[posn].dir[i].x[1]+dx[2]*g[posn].dir[i].x[2]);
    
    if(fabs(denominator) > 0){
      newdist=numerator/denominator;
      if(newdist<*ds && newdist > cutoff){
        *ds=newdist;
        *nposn=g[posn].neigh[i]->id;
      }
    }
  }
  if(*nposn==-1) *nposn=posn;
}


void
traceray(rayData ray, int tmptrans, int im, inputPars *par, struct grid *g, molData *m, image *img, int nlinetot, int *counta, int *countb, double cutoff){
  int ichan,posn,nposn,i,iline;
  double vfac=0.,x[3],dx[3],vThisChan;
  double deltav,ds,dist2,ndist2,xp,yp,zp,col,shift,jnu,alpha,remnantSnu,dtau,expDTau,snu_pol[3];

  for(ichan=0;ichan<img[im].nchan;ichan++){
    ray.tau[ichan]=0.0;
    ray.intensity[ichan]=0.0;
  }

  xp=ray.x;
  yp=ray.y;

  if((xp*xp+yp*yp)/par->radiusSqu <= 1 ) {
    zp=sqrt(par->radiusSqu-(xp*xp+yp*yp));

    for(i=0;i<3;i++){
      x[i]=xp*img[im].rotMat[i][0] + yp*img[im].rotMat[i][1] + zp*img[im].rotMat[i][2];
      dx[i]= -img[im].rotMat[i][2];
    }

    i=0;
    dist2=(x[0]-g[i].x[0])*(x[0]-g[i].x[0]) + (x[1]-g[i].x[1])*(x[1]-g[i].x[1]) + (x[2]-g[i].x[2])*(x[2]-g[i].x[2]);
    posn=i;
    for(i=1;i<par->ncell;i++){
      ndist2=(x[0]-g[i].x[0])*(x[0]-g[i].x[0]) + (x[1]-g[i].x[1])*(x[1]-g[i].x[1]) + (x[2]-g[i].x[2])*(x[2]-g[i].x[2]);
      if(ndist2<dist2){
        posn=i;
        dist2=ndist2;
      }
    }

    col=0;
    do{
      ds=2.*zp-col;
      nposn=-1;
      line_plane_intersect(g,&ds,posn,&nposn,dx,x,cutoff);
      if(par->polarization){
        for(ichan=0;ichan<img[im].nchan;ichan++){
          sourceFunc_pol(snu_pol,&dtau,ds,m,vfac,g,posn,0,0,img[im].theta);
#ifdef FASTEXP
          ray.intensity[ichan]+=FastExp(ray.tau[ichan])*(1.-exp(-dtau))*snu_pol[ichan];
#else
          ray.intensity[ichan]+=   exp(-ray.tau[ichan])*(1.-exp(-dtau))*snu_pol[ichan];
#endif
          ray.tau[ichan]+=dtau;
        }
      } else {
        for(ichan=0;ichan<img[im].nchan;ichan++){
          jnu=.0;
          alpha=0.;

          for(iline=0;iline<nlinetot;iline++){
            if(img[im].doline && m[counta[iline]].freq[countb[iline]] > img[im].freq-img[im].bandwidth/2.
            && m[counta[iline]].freq[countb[iline]] < img[im].freq+img[im].bandwidth/2.){
              if(img[im].trans > -1){
                shift=(m[counta[iline]].freq[countb[iline]]-m[counta[iline]].freq[img[im].trans])/m[counta[iline]].freq[img[im].trans]*CLIGHT;
              } else {
                shift=(m[counta[iline]].freq[countb[iline]]-img[im].freq)/img[im].freq*CLIGHT;
              }
              vThisChan=(ichan-(img[im].nchan-1)/2.)*img[im].velres; // Consistent with the WCS definition in writefits().
              deltav = vThisChan - img[im].source_vel + shift;

              if(!par->pregrid) velocityspline2(x,dx,ds,g[posn].mol[counta[iline]].binv,deltav,&vfac);
              else vfac=gaussline(deltav+veloproject(dx,g[posn].vel),g[posn].mol[counta[iline]].binv);

              sourceFunc_line(&jnu,&alpha,m,vfac,g,posn,counta[iline],countb[iline]);
            }
          }

          if(img[im].doline && img[im].trans > -1) sourceFunc_cont(&jnu,&alpha,g,posn,0,img[im].trans);
          else if(img[im].doline && img[im].trans == -1) sourceFunc_cont(&jnu,&alpha,g,posn,0,tmptrans);
          else sourceFunc_cont(&jnu,&alpha,g,posn,0,0);
          dtau=alpha*ds;
//???          if(dtau < -30) dtau = -30; // as in photon()?
          calcSourceFn(dtau, par, &remnantSnu, &expDTau);
          remnantSnu *= jnu*m[0].norminv*ds;
#ifdef FASTEXP
          ray.intensity[ichan]+=FastExp(ray.tau[ichan])*remnantSnu;
#else
          ray.intensity[ichan]+=   exp(-ray.tau[ichan])*remnantSnu;
#endif
          ray.tau[ichan]+=dtau;
        }
      }

      /* new coordinates */
      for(i=0;i<3;i++) x[i]+=ds*dx[i];
      col+=ds;
      posn=nposn;
    } while(col < 2*zp);

    /* add or subtract cmb */
#ifdef FASTEXP
    for(ichan=0;ichan<img[im].nchan;ichan++){
      ray.intensity[ichan]+=(FastExp(ray.tau[ichan])-1.)*m[0].local_cmb[tmptrans];
    }
#else
    for(ichan=0;ichan<img[im].nchan;ichan++){
      ray.intensity[ichan]+=(exp(-ray.tau[ichan])-1.)*m[0].local_cmb[tmptrans];
    }
#endif
  }
}


void
raytrace(int im, inputPars *par, struct grid *g, molData *m, image *img){
  int *counta,*countb,nlinetot,aa,numRaysThisPixel;
  int ichan,gi,iline,tmptrans,i,threadI,nRaysDone;
  double size,minfreq,absDeltaFreq;
  double cutoff;
  unsigned int totalNumImagePixels,ppi;
  double imgXiCentre,imgYiCentre,x[3];
  int xi,yi,ri,j;
  const gsl_rng_type *ranNumGenType = gsl_rng_ranlxs2;

  gsl_rng *randGen = gsl_rng_alloc(ranNumGenType);	/* Random number generator */
#ifdef TEST
  gsl_rng_set(randGen,178490);
#else
  gsl_rng_set(randGen,time(0));
#endif

  gsl_rng **threadRans;
  threadRans = malloc(sizeof(gsl_rng *)*par->nThreads);

  for (i=0;i<par->nThreads;i++){
    threadRans[i] = gsl_rng_alloc(ranNumGenType);
    gsl_rng_set(threadRans[i],(int)gsl_rng_uniform(randGen)*1e6);
  }

  size=img[im].distance*img[im].imgres;
  totalNumImagePixels = img[im].pxls*img[im].pxls;
  imgXiCentre = img[im].pxls/2.0;
  imgYiCentre = img[im].pxls/2.0;

  /* Determine whether there are blended lines or not */
  lineCount(par->nSpecies, m, &counta, &countb, &nlinetot);
  if(img[im].doline==0) nlinetot=1;

  /* Fix the image parameters */
  if(img[im].freq < 0) img[im].freq=m[0].freq[img[im].trans];
  if(img[im].nchan == 0 && img[im].bandwidth>0){
    img[im].nchan=(int) (img[im].bandwidth/(img[im].velres/CLIGHT*img[im].freq));
  } else if (img[im].velres<0 && img[im].bandwidth>0){
    img[im].velres = img[im].bandwidth*CLIGHT/img[im].freq/img[im].nchan;
  } else img[im].bandwidth = img[im].nchan*img[im].velres/CLIGHT * img[im].freq;

  if(img[im].trans<0){
    iline=0;
    minfreq=fabs(img[im].freq-m[0].freq[iline]);
    tmptrans=iline;
    for(iline=1;iline<m[0].nline;iline++){
      absDeltaFreq=fabs(img[im].freq-m[0].freq[iline]);
      if(absDeltaFreq<minfreq){
        minfreq=absDeltaFreq;
        tmptrans=iline;
      }
    }
  } else tmptrans=img[im].trans;

  for(ppi=0;ppi<totalNumImagePixels;ppi++){
    for(ichan=0;ichan<img[im].nchan;ichan++){
      img[im].pixel[ppi].intense[ichan]=0.0;
      img[im].pixel[ppi].tau[ichan]=0.0;
    }
  }

  if(par->raytraceAlgorithm==0){
    for(ppi=0;ppi<totalNumImagePixels;ppi++){
      img[im].pixel[ppi].numRays=par->antialias;
    }

  } else if(par->raytraceAlgorithm==1){
    /*
In this algorithm we attempt to deal better with image pixels which cover areas of the model where grid points are packed much more densely than the pixel spacing. In algorithm 0, such regions are sampled by only at par->antialias rays per pixel, which risks missing extreme values. The approach here is to generate approximately the same number of rays per pixel as the number of projected grid points in that pixel (with par->antialias providing a minimum floor). The set of values per pixel are then averaged to obtain the image value for that pixel.
    */

    for(ppi=0;ppi<totalNumImagePixels;ppi++){
      img[im].pixel[ppi].numRays=0;
    }

    for(gi=0;gi<par->pIntensity;gi++){
      /* Apply the inverse (i.e. transpose) rotation matrix. (We use the inverse matrix here because here we want to rotate grid coordinates to the observer frame, whereas inside traceray() we rotate observer coordinates to the grid frame.)
      */
      for(i=0;i<2;i++){
        x[i]=0.0;
        for(j=0;j<DIM;j++){
          x[i] += g[gi].x[j]*img[im].rotMat[j][i];
        }
      }

      /* Calculate which pixel the projected position (x[0],x[1]) falls within. */
      xi = floor(x[0]/size + imgXiCentre);
      yi = floor(x[1]/size + imgYiCentre);
      ppi = yi*img[im].pxls + xi;
      if(ppi>=0 && ppi<totalNumImagePixels)
        img[im].pixel[ppi].numRays++;
    }

  } else {
    if(!silent) bail_out("Unrecognized value of par.raytraceAlgorithm");
    exit(1);
  }

  cutoff = par->minScale*1.0e-7;

  nRaysDone=0;
  omp_set_dynamic(0);
  #pragma omp parallel private(ppi,aa,threadI) num_threads(par->nThreads)
  {
    threadI = omp_get_thread_num();

    /* Declaration of thread-private pointers */
    rayData ray;
    ray.intensity= malloc(sizeof(double)*img[im].nchan);
    ray.tau      = malloc(sizeof(double)*img[im].nchan);

    #pragma omp for
    /* Main loop through pixel grid. ***NOTE*** though that the division of labour between threads can be uneven if par->raytraceAlgorithm==1. Probably should rearrange the code, maybe so there is only 1 'for' loop over all the rays?
   */
    for(ppi=0;ppi<totalNumImagePixels;ppi++){
      xi = ppi%img[im].pxls;
      yi = ppi/img[im].pxls;

      if(img[im].pixel[ppi].numRays>par->antialias)
        numRaysThisPixel = img[im].pixel[ppi].numRays;
      else
        numRaysThisPixel = par->antialias;

      #pragma omp atomic
      ++nRaysDone;

      for(aa=0;aa<numRaysThisPixel;aa++){
        ray.x = size*(gsl_rng_uniform(threadRans[threadI]) + xi - imgXiCentre);
        ray.y = size*(gsl_rng_uniform(threadRans[threadI]) + yi - imgYiCentre);

        traceray(ray, tmptrans, im, par, g, m, img, nlinetot, counta, countb, cutoff);

        #pragma omp critical
        {
          for(ichan=0;ichan<img[im].nchan;ichan++){
            img[im].pixel[ppi].intense[ichan] += ray.intensity[ichan]/(double)numRaysThisPixel;
            img[im].pixel[ppi].tau[    ichan] += ray.tau[      ichan]/(double)numRaysThisPixel;
          }
        }
      }
      if (threadI == 0){ // i.e., is master thread
        if(!silent) progressbar((double)(nRaysDone)/(double)(totalNumImagePixels-1), 13);
      }
    }

    free(ray.tau);
    free(ray.intensity);
  } // end of parallel block.

  img[im].trans=tmptrans;

  free(counta);
  free(countb);
  for (i=0;i<par->nThreads;i++){
    gsl_rng_free(threadRans[i]);
  }
  free(threadRans);
  gsl_rng_free(randGen);
}

