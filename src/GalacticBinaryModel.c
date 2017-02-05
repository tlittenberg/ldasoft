//
//  GalacticBinaryModel.c
//  
//
//  Created by Littenberg, Tyson B. (MSFC-ZP12) on 1/15/17.
//
//
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "LISA.h"
#include "Constants.h"
#include "GalacticBinary.h"
#include "GalacticBinaryModel.h"

void map_array_to_params(struct Source *source, double *params)
{
  source->f0       = params[0];
  source->costheta = params[1];
  source->phi      = params[2];
  source->amp      = params[3];
  source->cosi     = params[4];
  source->psi      = params[5];
  source->phi0     = params[6];
  source->dfdt     = params[7];
}

void map_params_to_array(struct Source *source, double *params)
{
  params[0] = source->f0;
  params[1] = source->costheta;
  params[2] = source->phi;
  params[3] = source->amp;
  params[4] = source->cosi;
  params[5] = source->psi;
  params[6] = source->phi0;
  params[7] = source->dfdt;
}

void initialize_chain(struct Chain *chain, long seed, int i)
{
  chain->index = i;
  chain->T = gsl_rng_default;
  chain->r = gsl_rng_alloc(chain->T);
  gsl_rng_env_setup();
  gsl_rng_set (chain->r, seed);
}

void initialize_model(struct Model *model, int Nmax, int NFFT, int Nchannel)
{
  model->Nmax = Nmax;
  model->source = malloc(model->Nmax*sizeof(struct Source *));
  model->noise  = malloc(sizeof(struct Noise));
  model->tdi    = malloc(sizeof(struct TDI));
  
  int n;

  for(n=0; n<model->Nmax; n++)
  {
    model->source[n] = malloc(sizeof(struct Source));
    initialize_source(model->source[n],NFFT,Nchannel);
  }
  
  initialize_tdi(model->tdi,NFFT, Nchannel);

  initialize_noise(model->noise,NFFT);
}

void free_model(struct Model *model)
{
  int n;
  for(n=0; n<model->Nmax; n++)
  {
    free(model->source[n]);
  }
  free(model->source);
  
  free(model->tdi);
  free(model->noise);
  free(model);
}

void initialize_tdi(struct TDI *tdi, int NFFT, int Nchannel)
{
  //Number of frequency bins (2*N samples)
  tdi->N = NFFT;

  //Michelson
  tdi->X = malloc(2*tdi->N*sizeof(double));
  tdi->Y = malloc(2*tdi->N*sizeof(double));
  tdi->Z = malloc(2*tdi->N*sizeof(double));
  
  //Noise-orthogonal
  tdi->A = malloc(2*tdi->N*sizeof(double));
  tdi->E = malloc(2*tdi->N*sizeof(double));
  tdi->T = malloc(2*tdi->N*sizeof(double));
  
  int n;
  for(n=0; n<2*tdi->N; n++)
  {
    tdi->X[n] = 0.0;
    tdi->Y[n] = 0.0;
    tdi->Z[n] = 0.0;
    tdi->A[n] = 0.0;
    tdi->E[n] = 0.0;
    tdi->T[n] = 0.0;
  }
  
  //Number of TDI channels (X or A&E or maybe one day A,E,&T)
  tdi->Nchannel = Nchannel;
}

void free_tdi(struct TDI *tdi)
{
  free(tdi->X);
  free(tdi->Y);
  free(tdi->Z);
  free(tdi->A);
  free(tdi->E);
  free(tdi->T);
  free(tdi);
}

void initialize_noise(struct Noise *noise, int NFFT)
{
  noise->etaA = 1.0;
  noise->etaE = 1.0;
  noise->etaX = 1.0;
  
  noise->SnA = malloc(NFFT*sizeof(double));
  noise->SnE = malloc(NFFT*sizeof(double));
  noise->SnX = malloc(NFFT*sizeof(double));
  
  int n;
  for(n=0; n<NFFT; n++)
  {
    noise->SnA[n]=1.0;
    noise->SnE[n]=1.0;
    noise->SnX[n]=1.0;
  }
}

void free_noise(struct Noise *noise)
{
  free(noise->SnA);
  free(noise->SnE);
  free(noise->SnX);
  free(noise);
}

void initialize_source(struct Source *source, int NFFT, int Nchannel)
{
  //Intrinsic
  source->m1=1.;
  source->m2=1.;
  source->f0=0.;
  
  //Extrinisic
  source->psi=0.0;
  source->cosi=0.0;
  source->phi0=0.0;
  
  source->D=1.0;
  source->phi=0.0;
  source->costheta=0.0;
  
  //Derived
  source->amp=1.;
  source->Mc=1.;
  source->dfdt=0.;
  
  //Package parameters for waveform generator
  source->params=malloc(8*sizeof(double));
  
  //Response
  source->tdi = malloc(sizeof(struct TDI));
  initialize_tdi(source->tdi,NFFT, Nchannel);
  
};

void free_source(struct Source *source)
{
  free(source->tdi);
  free(source);
}


void simualte_data(struct Data *data, struct Flags *flags, struct Source **injections, int Ninj)
{
  int i;
  
  double f,fmin,fmax;
  
  //find minimum and frequency of injections
  fmin = 1.0e6;
  fmax =-1.0e6;
  
  for(i=0; i<Ninj; i++)
  {
    f = injections[i]->f0;
    if(f<fmin) fmin=f;
    if(f>fmax) fmax=f;
  }
  
  //get boundaries of data segment
  data->fmin = fmin;
  data->fmin = fmax;
  data->qmin = (int)floor(data->fmin*data->T) - data->N/2;
  data->qmax = data->qmin + data->N;
  
  //check we have enough data for all of the injections
  if( data->qmax < (int)floor(data->fmax*data->T) )
  {
    fprintf(stdout,"Error:  Data bandwdith does not contain all of the injections\n");
    fprintf(stdout,"        Injections covere frequences [%g,%g]\n",data->fmin,data->fmax);
    fprintf(stdout,"        Data covers frequencies      [%g,%g]\n",data->qmin/data->T,data->qmax/data->T);
    fprintf(stdout,"        Use large N or less injections\n");
    exit(255);
  }
  
}
