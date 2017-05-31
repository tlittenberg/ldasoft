//
//  GalacticBinaryFisher.c
//  
//
//  Created by Littenberg, Tyson B. (MSFC-ZP12) on 4/4/17.
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
  parse(argc,argv,data_ptr,orbit,flags,chain,NMAX);
  
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
  
  fprintf(stdout,"\n==== GalacticBinaryInjectSimulatedSource ====\n");
  
  /* Get injection parameters */
  double f0,dfdt,theta,phi,amp,iota,psi,phi0; //read from injection file
  
  FILE *injectionFile = fopen(flags->injFile[0],"r");
  if(!injectionFile)
    fprintf(stderr,"Missing catalog file %s\n",flags->injFile[0]);
  else
    fprintf(stdout,"Simulateing binary catalog %s\n",flags->injFile[0]);
  
  //count sources in file
  int N=0;
  while(!feof(injectionFile))
  {
    fscanf(injectionFile,"%lg %lg %lg %lg %lg %lg %lg %lg",&f0,&dfdt,&theta,&phi,&amp,&iota,&psi,&phi0);
    N++;
  }
  rewind(injectionFile);
  N--;
  
  fprintf(stdout,"Found %i sources in %s\n",N,flags->injFile[0]);
  
  FILE *outfile = fopen("snr.dat","w");
  int goodPE = 0;
  int goodCal = 0;

  FILE *peFile  = fopen("PrecisionBinaries.txt","w");
  FILE *calFile = fopen("CalibrationBinaries.txt","w");
  FILE *pFile;
  char filename[128];
  for(int n=0; n<N; n++)
  {
    
    fscanf(injectionFile,"%lg %lg %lg %lg %lg %lg %lg %lg",&f0,&dfdt,&theta,&phi,&amp,&iota,&psi,&phi0);
    
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
    inj->costheta = cos(M_PI/2. - theta);
    inj->phi      = phi;
    inj->amp      = amp;
    inj->cosi     = cos(iota);
    inj->phi0     = phi0;
    inj->psi      = psi;
    
    map_params_to_array(inj, inj->params, data->T);
    
    //Book-keeping of injection time-frequency volume
    galactic_binary_alignment(orbit, data, inj);
    
    //Simulate gravitational wave signal
    double t0 = data->t0[0];
    galactic_binary(orbit, data->T, t0, inj->params, 8, inj->tdi->X, inj->tdi->A, inj->tdi->E, inj->BW, 2);
    
    //Get noise spectrum for data segment
    for(int n=0; n<data->N; n++)
    {
      double f = data->fmin + (double)(n)/data->T;
      data->noise[0]->SnA[n] = AEnoise(orbit->L, orbit->fstar, f);
      data->noise[0]->SnE[n] = AEnoise(orbit->L, orbit->fstar, f);
    }
    
    //Get injected SNR
    double SNR = snr(inj, data->noise[0]);
    double Mc  = galactic_binary_Mc(f0, dfdt, data->T);
    double dL  = galactic_binary_dL(f0, dfdt, amp, data->T);
    
    fprintf(outfile,"%g %g %g %g %g %g %g %g %g\n",f0,dfdt,amp,cos(iota),Mc,dL,cos(M_PI/2 - theta),phi,SNR);
    
    //Check if sources meet the "good PE candidate" requirement
    double eclipse = fabs(cos(iota)/pow((f0/3.5e-3),(2./3.)));
    double fddot = (11.0/3.0)*dfdt*dfdt/f0;
    if(eclipse<0.3 && fddot*pow(5*YEAR,3.)>0.1)
    {
      fprintf(stdout,"PE: %g %g %g %g %g %g %g %g %g\n",f0,dfdt,amp,cos(iota),Mc,dL,cos(M_PI/2 - theta),phi,SNR);
      fprintf(peFile,"%lg %lg %lg %lg %lg %lg %lg %lg",f0,dfdt,theta,phi,amp,iota,psi,phi0);

      sprintf(filename,"PrecisionSource_%i.txt",goodPE);
      pFile = fopen(filename,"w");
      fprintf(pFile,"%lg %lg %lg %lg %lg %lg %lg %lg\n",f0,dfdt,theta,phi,amp,iota,psi,phi0);
      fclose(pFile);
      
      goodPE++;
    }
    
    //Check if sources meet the "good standard sirens" requirement
    if(SNR>1000 && dfdt>0)
    {
      fprintf(stdout,"Calib: %g %g %g %g %g %g %g %g %g\n",f0,dfdt,amp,cos(iota),Mc,dL,cos(M_PI/2 - theta),phi,SNR);
      fprintf(calFile,"%lg %lg %lg %lg %lg %lg %lg %lg\n",f0,dfdt,theta,phi,amp,iota,psi,phi0);
      
      sprintf(filename,"CalibrationSource_%i.txt",goodCal);
      pFile = fopen(filename,"w");
      fprintf(pFile,"%lg %lg %lg %lg %lg %lg %lg %lg\n",f0,dfdt,theta,phi,amp,iota,psi,phi0);
      fclose(pFile);
      goodCal++;
    }
  }
  
  printf("Good PE sources:  %i\n",goodPE);
  printf("Good Calibration soruces:  %i\n",goodCal);
  
  fclose(injectionFile);
  
  fprintf(stdout,"================================================\n\n");
  
  //print total run time
  stop = time(NULL);
  
  printf(" ELAPSED TIME = %g second\n",(double)(stop-start));

  return 0;
}
