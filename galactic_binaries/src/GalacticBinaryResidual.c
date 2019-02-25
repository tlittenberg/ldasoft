//
//  GalacticBinaryResidual.c
//
//
//  Created by Littenberg, Tyson B. (MSFC-ST12) on 8/24/18.
//
//

/***************************  REQUIRED LIBRARIES  ***************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/*************  PROTOTYPE DECLARATIONS FOR INTERNAL FUNCTIONS  **************/

#include "LISA.h"
#include "Constants.h"
#include "GalacticBinary.h"
#include "GalacticBinaryIO.h"
#include "GalacticBinaryData.h"
#include "GalacticBinaryMath.h"
#include "GalacticBinaryModel.h"
#include "GalacticBinaryWaveform.h"


int main(int argc, char *argv[])
{
  int NMAX   = 1;   //max number of waveforms & time segments
  
  time_t start, stop;
  start = time(NULL);
  
  
  /* Allocate data structures */
  struct Flags *flags       = malloc(sizeof(struct Flags));
  struct Orbit *orbit       = malloc(sizeof(struct Orbit));
  struct Data  **data_ptr  = malloc(sizeof(struct Data**)*NMAX);
  struct Chain *chain       = malloc(sizeof(struct Chain));
  
  /* Parse command line and set defaults/flags */
  for(int i=0; i<NMAX; i++)
  {
    data_ptr[i] = malloc(sizeof(struct Data));
  }
  parse(argc,argv,data_ptr,orbit,flags,chain,NMAX,NMAX);
  
  /* Load spacecraft ephemerides */
  switch(flags->orbit)
  {
    case 0:
      initialize_analytic_orbit(orbit);
      break;
    case 1:
      initialize_numeric_orbit(orbit);
      break;
    default:
      fprintf(stderr,"unsupported orbit type\n");
      return(1);
      break;
  }
  
  /* Initialize data structures */
  alloc_data(data_ptr, flags);
  
  struct Data *data = data_ptr[0];
  
  fprintf(stdout,"\n==== GalacticBinarySubtractDetectedSources ====\n");
  
  /* Get injection parameters */
  double f0,dfdt,theta,phi,amp,iota,psi,phi0; //read from injection file
  
  FILE *catalogFile = fopen(flags->injFile[0],"r");
  if(!catalogFile)
    fprintf(stderr,"Missing catalog file %s\n",flags->injFile[0]);
  else
    fprintf(stdout,"Subtracting binary catalog %s\n",flags->injFile[0]);
  
  /* Get Galaxy File */
  GalacticBinaryReadData(data_ptr,orbit,flags);
  
  //count sources in file
  int N=0;
  fprintf(stdout,"Counting sources in catalog file %s:\n",flags->injFile[0]);
  while(!feof(catalogFile))
  {
    fscanf(catalogFile,"%lg %lg %lg %lg %lg %lg %lg %lg",&f0,&dfdt,&amp,&phi,&theta,&iota,&psi,&phi0);
    N++;
  }
  rewind(catalogFile);
  N--;
  
  fprintf(stdout,"Found %i sources in %s\n",N,flags->injFile[0]);
  

  fprintf(stdout,"Generating waveforms:\n");
  
  
  
  
  
  
  
  
  
  
  
  for(int i=0; i<data->N; i++)
  {
    data->tdi[0]->A[2*i]   = 0.0;// inj->tdi->A[2*n];
    data->tdi[0]->E[2*i]   = 0.0;//inj->tdi->E[2*n];
    data->tdi[0]->A[2*i+1] = 0.0;//inj->tdi->A[2*n+1];
    data->tdi[0]->E[2*i+1] = 0.0;//inj->tdi->E[2*n+1];
  }

  
  
  
  
  
  
  
  
  
  
  
  for(int nn=0; nn<N; nn++)
  {
    if(nn%(N/100)==0)printProgress( (double)nn / (double)N );
    fscanf(catalogFile,"%lg %lg %lg %lg %lg %lg %lg %lg",&f0,&dfdt,&amp,&phi,&theta,&iota,&psi,&phi0);
    
    //set bandwidth of data segment centered on injection
    data->fmin = f0 - (data->N/2)/data->T;
    data->fmax = f0 + (data->N/2)/data->T;
    data->qmin = (int)(data->fmin*data->T);
    data->qmax = data->qmin+data->N;
    
    struct Source *inj = data->inj;
    
    for(int n=0; n<2*data->N; n++)
    {
      inj->tdi->A[n] = 0.0;
      inj->tdi->E[n] = 0.0;
      inj->tdi->X[n] = 0.0;
    }
    
    //map parameters to vector
    inj->f0       = f0;
    inj->dfdt     = dfdt;
    inj->costheta = theta;
    inj->phi      = phi;
    inj->amp      = amp;
    inj->cosi     = iota;
    inj->phi0     = phi0;
    inj->psi      = psi;
    
    map_params_to_array(inj, inj->params, data->T);
    
    //Book-keeping of injection time-frequency volume
    galactic_binary_alignment(orbit, data, inj);
    
    //Simulate gravitational wave signal
    double t0 = data->t0[0];
    galactic_binary(orbit, data->format, data->T, t0, inj->params, 8, inj->tdi->X, inj->tdi->A, inj->tdi->E, inj->BW, 2);
    
    
    if(inj->BW > data->N) printf("WARNING:  Bandwidth %i wider than N %i at f=%.2e\n",inj->BW,data->N,data->fmin);

    for(int n=0; n<inj->BW; n++)
    {
      int i = inj->qmin+n;
      if(i>0 && i<data->N)
      {
//        data->tdi[0]->A[2*i]   -= inj->tdi->A[2*n];
//        data->tdi[0]->E[2*i]   -= inj->tdi->E[2*n];
//        data->tdi[0]->A[2*i+1] -= inj->tdi->A[2*n+1];
//        data->tdi[0]->E[2*i+1] -= inj->tdi->E[2*n+1];
        data->tdi[0]->A[2*i]   += inj->tdi->A[2*n];
        data->tdi[0]->E[2*i]   += inj->tdi->E[2*n];
        data->tdi[0]->A[2*i+1] += inj->tdi->A[2*n+1];
        data->tdi[0]->E[2*i+1] += inj->tdi->E[2*n+1];
      }
    }
  }
  fprintf(stdout,"\nFinished subtracting galaxy catalog\n");
    
  //print full galaxy and density file
  FILE *galaxyFile  = fopen("galaxy_data_AE_clean.dat","w");

  for(int n=0; n<data->N; n++)
  {
    double f = (double)n/data->T;
    fprintf(galaxyFile,"%.12g ",f);
    fprintf(galaxyFile,"%.12g %.12g ",data->tdi[0]->A[2*n],data->tdi[0]->A[2*n+1]);
    fprintf(galaxyFile,"%.12g %.12g\n",data->tdi[0]->E[2*n],data->tdi[0]->E[2*n+1]);
  }
  fclose(galaxyFile);

    
  fprintf(stdout,"================================================\n\n");
  
  //print total run time
  stop = time(NULL);
  
  printf(" ELAPSED TIME = %g second\n",(double)(stop-start));
  
  return 0;
}
