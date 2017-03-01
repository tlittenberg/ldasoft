//
//  GalacticBinaryData.c
//  
//
//  Created by Littenberg, Tyson B. (MSFC-ZP12) on 2/3/17.
//
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "LISA.h"
#include "GalacticBinary.h"
#include "GalacticBinaryData.h"
#include "GalacticBinaryMath.h"
#include "GalacticBinaryModel.h"
#include "GalacticBinaryWaveform.h"

void GalacticBinaryReadData(struct Data *data)
{
  
}
void GalacticBinarySimulateData(struct Data *data)
{
  
}

void GalacticBinaryInjectVerificationSource(struct Data **data_vec, struct Orbit *orbit, struct Flags *flags)
{
  //TODO: support Michelson-only injection
  fprintf(stdout,"\n==== GalacticBinaryInjectVerificationSource ====\n");
  
  FILE *fptr;
  
  /* Get injection parameters */
  double f0,dfdt,costheta,phi,m1,m2,D; //read from injection file
  double cosi,phi0,psi;                //drawn from prior
  double Mc,amp;                       //calculated
  
  FILE *injectionFile;
  char filename[1024];
  
  for(int ii = 0; ii<flags->injection; ii++)
  {
    struct Data *data  = data_vec[ii];
    struct Source *inj = data->inj;
  
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc(T);
    gsl_rng_env_setup();
    gsl_rng_set (r, data->iseed);
  
    injectionFile = fopen(flags->injFile[ii],"r");
    if(!injectionFile)
      fprintf(stderr,"Missing injection file %s\n",flags->injFile[ii]);
    else
      fprintf(stdout,"Injecting verification binary %s  (%i/%i)\n",flags->injFile[ii],ii+1, flags->injection);
    
    fscanf(injectionFile,"%lg %lg %lg %lg %lg %lg %lg",&f0,&dfdt,&costheta,&phi,&m1,&m2,&D);
    
    //set bandwidth of data segment centered on injection
    data->fmin = f0 - (data->N/2)/data->T;
    data->fmax = f0 + (data->N/2)/data->T;
    data->qmin = (int)(data->fmin*data->T);
    data->qmax = data->qmin+data->N;
    
    //draw extrinsic parameters
    //TODO: support for verification binary priors
    cosi = -1.0 + gsl_rng_uniform(r)*2.0;
    phi0 = gsl_rng_uniform(r)*M_PI*2.;
    psi  = gsl_rng_uniform(r)*M_PI/4.;
    
    //compute derived parameters
    Mc  = chirpmass(m1,m2);
    amp = galactic_binary_Amp(Mc, f0, D, data->T);
    
    //map parameters to vector
    inj->f0       = f0;
    inj->dfdt     = dfdt;
    inj->costheta = costheta;
    inj->phi      = phi;
    inj->amp      = amp;
    inj->cosi     = cosi;
    inj->phi0     = phi0;
    inj->psi      = psi;
    map_params_to_array(inj, inj->params, data->T);
    
    //Book-keeping of injection time-frequency volume
    galactic_binary_alignment(orbit, data, inj);
    
    //Simulate gravitational wave signal
    galactic_binary(orbit, data->T, data->t0, inj->params, inj->tdi->X, inj->tdi->A, inj->tdi->E, inj->BW, 2);
    
    //Add waveform to data TDI channels
    for(int n=0; n<inj->BW; n++)
    {
      int i = n+inj->imin;
      
      data->tdi->X[2*i]   = inj->tdi->X[2*n];
      data->tdi->X[2*i+1] = inj->tdi->X[2*n+1];
      
      data->tdi->A[2*i]   = inj->tdi->A[2*n];
      data->tdi->A[2*i+1] = inj->tdi->A[2*n+1];
      
      data->tdi->E[2*i]   = inj->tdi->E[2*n];
      data->tdi->E[2*i+1] = inj->tdi->E[2*n+1];
    }
    
    sprintf(filename,"power_injection_%i.dat",ii);
    fptr=fopen(filename,"w");
    for(int i=0; i<data->N; i++)
    {
      double f = (double)(i+data->qmin)/data->T;
      fprintf(fptr,"%lg %lg %lg ",
              f,
              data->tdi->A[2*i]*data->tdi->A[2*i]+data->tdi->A[2*i+1]*data->tdi->A[2*i+1],
              data->tdi->E[2*i]*data->tdi->E[2*i]+data->tdi->E[2*i+1]*data->tdi->E[2*i+1]);
      fprintf(fptr,"\n");
    }
    fclose(fptr);
    
    //Get noise spectrum for data segment
    for(int n=0; n<data->N; n++)
    {
      double f = data->fmin + (double)(n)/data->T;
      data->noise->SnA[n] = AEnoise(orbit->L, orbit->fstar, f);
      data->noise->SnE[n] = AEnoise(orbit->L, orbit->fstar, f);
    }
    
    //Get injected SNR
    fprintf(stdout,"   ...injected SNR=%g\n",snr(inj, data->noise));
    
    //Add Gaussian noise to injection
    gsl_rng_set (r, data->nseed);
    
    if(!flags->zeroNoise)
    {
      printf("   ...adding Gaussian noise realization\n");
      
      for(int n=0; n<data->N; n++)
      {
        data->tdi->A[2*n]   += gsl_ran_gaussian (r, 1)*sqrt(data->noise->SnA[n])/2.;
        data->tdi->A[2*n+1] += gsl_ran_gaussian (r, 1)*sqrt(data->noise->SnA[n])/2.;
        
        data->tdi->E[2*n]   += gsl_ran_gaussian (r, sqrt(data->noise->SnE[n])/2.);
        data->tdi->E[2*n+1] += gsl_ran_gaussian (r, sqrt(data->noise->SnE[n])/2.);
      }
    }
    
    //Compute fisher information matrix of injection
    printf("   ...computing Fisher Information Matrix of injection\n");
    
    galactic_binary_fisher(orbit, data, inj, data->noise);
    
    /*
     printf("\n Fisher Matrix:\n");
     for(int i=0; i<8; i++)
     {
     fprintf(stdout," ");
     for(int j=0; j<8; j++)
     {
     if(inj->fisher_matrix[i][j]<0)fprintf(stdout,"%.2e ", inj->fisher_matrix[i][j]);
     else                          fprintf(stdout,"+%.2e ",inj->fisher_matrix[i][j]);
     }
     fprintf(stdout,"\n");
     }
     
     printf("\n Fisher std. errors:\n");
     for(int j=0; j<8; j++)  fprintf(stdout," %.4e\n", 1./sqrt(inj->fisher_evalue[j]));
     */
    
    
    sprintf(filename,"power_data_%i.dat",ii);
    fptr=fopen(filename,"w");

    for(int i=0; i<data->N; i++)
    {
      double f = (double)(i+data->qmin)/data->T;
      fprintf(fptr,"%lg %lg %lg ",
              f,
              data->tdi->A[2*i]*data->tdi->A[2*i]+data->tdi->A[2*i+1]*data->tdi->A[2*i+1],
              data->tdi->E[2*i]*data->tdi->E[2*i]+data->tdi->E[2*i+1]*data->tdi->E[2*i+1]);
      fprintf(fptr,"\n");
    }
    fclose(fptr);
    fclose(injectionFile);
    gsl_rng_free(r);
  }
  
  fprintf(stdout,"================================================\n\n");
}

void GalacticBinaryInjectSimulatedSource(struct Data *data)
{
  
}
