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
#include "GalacticBinaryMath.h"
#include "GalacticBinaryModel.h"
#include "GalacticBinaryWaveform.h"

void map_array_to_params(struct Source *source, double *params, double T)
{
  source->f0       = params[0]/T;
  source->costheta = params[1];
  source->phi      = params[2];
  source->amp      = exp(params[3]);
  source->cosi     = params[4];
  source->psi      = params[5];
  source->phi0     = params[6];
  source->dfdt     = params[7]/(T*T);
}

void map_params_to_array(struct Source *source, double *params, double T)
{
  params[0] = source->f0*T;
  params[1] = source->costheta;
  params[2] = source->phi;
  params[3] = log(source->amp);
  params[4] = source->cosi;
  params[5] = source->psi;
  params[6] = source->phi0;
  params[7] = source->dfdt*T*T;
}

void initialize_chain(struct Chain *chain, long *seed, int NC)
{
  int ic;
  chain->NC = NC;
  chain->index = malloc(NC*sizeof(int));
  chain->temperature = malloc(NC*sizeof(double));
  for(ic=0; ic<NC; ic++)
  {
    chain->index[ic]=ic;
    chain->temperature[ic] = pow(1.2,(double)ic);
  }
  
  
  chain->r = malloc(NC*sizeof(gsl_rng *));
  chain->T = malloc(NC*sizeof(const gsl_rng_type *));
  
  for(ic=0; ic<NC; ic++)
  {
    chain->T[ic] = gsl_rng_default;
    chain->r[ic] = gsl_rng_alloc(chain->T[ic]);
    gsl_rng_env_setup();
    gsl_rng_set (chain->r[ic], *seed);
    *seed = (long)gsl_rng_get(chain->r[ic]);
  }
}

void free_chain(struct Chain *chain)
{
  free(chain->index);
  
  for(int ic=0; ic<chain->NC; ic++) gsl_rng_free(chain->r[ic]);
  free(chain->r);
  free(chain);
}

void alloc_model(struct Model *model, int Nmax, int NFFT, int Nchannel)
{
  model->Nlive  = 1;
  model->Nmax   = Nmax;
  model->source = malloc(model->Nmax*sizeof(struct Source *));
  model->noise  = malloc(sizeof(struct Noise));
  model->tdi    = malloc(sizeof(struct TDI));
  
  int n;

  for(n=0; n<model->Nmax; n++)
  {
    model->source[n] = malloc(sizeof(struct Source));
    alloc_source(model->source[n],NFFT,Nchannel);
  }
  
  alloc_tdi(model->tdi,NFFT, Nchannel);

  alloc_noise(model->noise,NFFT);
  
  model->prior = malloc(8*sizeof(double *));
  for(n=0; n<8; n++) model->prior[n] = malloc(2*sizeof(double));
}

void copy_model(struct Model *origin, struct Model *copy)
{
  copy->logL           = origin->logL;
  copy->Nlive          = origin->Nlive;
  copy->Nmax           = origin->Nmax;
  copy->logPriorVolume = origin->logPriorVolume;

  for(int n=0; n<origin->Nmax; n++) copy_source(origin->source[n],copy->source[n]);
  
  //copy_tdi(origin->tdi,copy->tdi);
  copy_noise(origin->noise,copy->noise);
  
  for(int n=0; n<8; n++) for(int j=0; j<2; j++) copy->prior[n][j] = origin->prior[n][j];
}

void free_model(struct Model *model)
{
  int n;
  for(n=0; n<model->Nmax; n++)
  {
    free(model->source[n]);
  }
  free(model->source);
  
  for(n=0; n<8; n++) free(model->prior[n]);
  free(model->prior);
  
  free(model->tdi);
  free(model->noise);
  free(model);
}

void alloc_tdi(struct TDI *tdi, int NFFT, int Nchannel)
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

void copy_tdi(struct TDI *origin, struct TDI *copy)
{
  copy->N        = origin->N;
  copy->Nchannel = origin->Nchannel;
  
  for(int n=0; n<2*origin->N; n++)
  {
    copy->X[n] = origin->X[n];
    copy->Y[n] = origin->Y[n];
    copy->Z[n] = origin->Z[n];
    copy->A[n] = origin->A[n];
    copy->E[n] = origin->E[n];
    copy->T[n] = origin->T[n];
  }
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

void alloc_noise(struct Noise *noise, int NFFT)
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

void copy_noise(struct Noise *origin, struct Noise *copy)
{
  copy->etaA = origin->etaA;
  copy->etaE = origin->etaE;
  copy->etaX = origin->etaX;
}

void free_noise(struct Noise *noise)
{
  free(noise->SnA);
  free(noise->SnE);
  free(noise->SnX);
  free(noise);
}

void alloc_source(struct Source *source, int NFFT, int Nchannel)
{
  int NP = 8;
  
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
  
  //Book-keeping
  source->BW   = NFFT;
  source->qmin = 0;
  source->qmax = NFFT;
  source->imin = 0;
  source->imax = NFFT;
  
  source->t0 = 0.0;

  
  //Package parameters for waveform generator
  source->params=malloc(NP*sizeof(double));
  
  //Response
  source->tdi = malloc(sizeof(struct TDI));
  alloc_tdi(source->tdi,NFFT, Nchannel);
  
  //FIsher
  source->fisher_matrix = malloc(NP*sizeof(double *));
  source->fisher_evectr = malloc(NP*sizeof(double *));
  source->fisher_evalue = malloc(NP*sizeof(double));
  for(int i=0; i<NP; i++)
  {
    source->fisher_matrix[i] = malloc(NP*sizeof(double));
    source->fisher_evectr[i] = malloc(NP*sizeof(double));
  }
};

void copy_source(struct Source *origin, struct Source *copy)
{
  int NP = 8;
  
  //Intrinsic
  copy->m1 = origin->m1;
  copy->m2 = origin->m2;
  copy->f0 = origin->f0;
  
  //Extrinisic
  copy->psi  = origin->psi;
  copy->cosi = origin->cosi;
  copy->phi0 = origin->phi0;
  
  copy->D        = origin->D;
  copy->phi      = origin->phi;
  copy->costheta = origin->costheta;
  
  //Derived
  copy->amp  = origin->amp;
  copy->Mc   = origin->Mc;
  copy->dfdt = origin->dfdt;
  
  //Book-keeping
  copy->BW   = origin->BW;
  copy->qmin = origin->qmin;
  copy->qmax = origin->qmax;
  copy->imin = origin->imin;
  copy->imax = origin->imax;
  
  copy->t0 = origin->t0;
  
  
  //Response
  //copy_tdi(origin->tdi,copy->tdi);
  
  //FIsher
  for(int i=0; i<NP; i++)
  {
    for(int j=0; j<NP; j++)
    {
      copy->fisher_matrix[i][j] = origin->fisher_evectr[i][j];
      copy->fisher_evectr[i][j] = origin->fisher_evectr[i][j];
    }
    copy->fisher_evalue[i] = origin->fisher_evalue[i];
    copy->params[i]        = origin->params[i];
  }
}

void free_source(struct Source *source)
{
  free(source->tdi);
  free(source);
  
  int NP=8;
  for(int i=0; i<NP; i++)
  {
    free(source->fisher_matrix[i]);
    free(source->fisher_evectr[i]);
  }
  free(source->fisher_matrix);
  free(source->fisher_evectr);
  free(source->fisher_evalue);

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

double gaussian_log_likelihood(struct Orbit *orbit, struct Data *data, struct Model *model)
{
  
  int n;
  int i,j;
  int N2=data->N*2;
  
  /**************************/
  /*                        */
  /*  Form master template  */
  /*                        */
  /**************************/

  for(n=0; n<N2; n++)
  {
    model->tdi->X[n]=0.0;
    model->tdi->A[n]=0.0;
    model->tdi->E[n]=0.0;
  }

  struct Source *source;

  //Loop over signals in model
  for(n=0; n<model->Nlive; n++)
  {
    source = model->source[n];
    
    //Book-keeping of injection time-frequency volume
    galactic_binary_alignment(orbit, data, source);
    
    //Simulate gravitational wave signal
    galactic_binary(orbit, data->T, source->t0, source->params, source->tdi->X, source->tdi->A, source->tdi->E, source->BW, source->tdi->Nchannel);
    
    //Add waveform to model TDI channels
    for(i=0; i<source->BW; i++)
    {
      j = i+source->imin;
            
      if(j>-1 && j<data->N)
      {
        int i_re = 2*i;
        int i_im = i_re+1;
        int j_re = 2*j;
        int j_im = j_re+1;
        
        model->tdi->X[j_re] += source->tdi->X[i_re];
        model->tdi->X[j_im] += source->tdi->X[i_im];
        
        model->tdi->A[j_re] += source->tdi->A[i_re];
        model->tdi->A[j_im] += source->tdi->A[i_im];
        
        model->tdi->E[j_re] += source->tdi->E[i_re];
        model->tdi->E[j_im] += source->tdi->E[i_im];
      }//check that index is in range
    }//loop over waveform bins
  }//loop over sources

  /**************************/
  /*                        */
  /* Form residual and sum  */
  /*                        */
  /**************************/
  
  struct TDI *residual = malloc(sizeof(struct TDI));
  alloc_tdi(residual, data->N, data->Nchannel);
  
  for(i=0; i<N2; i++)
  {
    residual->X[i] = data->tdi->X[i] - model->tdi->X[i];
    residual->A[i] = data->tdi->A[i] - model->tdi->A[i];
    residual->E[i] = data->tdi->E[i] - model->tdi->E[i];
  }

  double logL = 0.0;
  switch(data->Nchannel)
  {
    case 1:
      logL += -0.5*fourier_nwip(residual->X, residual->X, data->noise->SnX, data->N)/model->noise->etaX;
      logL +=  2.0*(double)data->N*log(model->noise->etaX);
      break;
    case 2:
      logL += -0.5*fourier_nwip(residual->A, residual->A, data->noise->SnA, data->N)/model->noise->etaA;
      logL += -0.5*fourier_nwip(residual->E, residual->E, data->noise->SnE, data->N)/model->noise->etaE;
      logL +=  2.0*(double)data->N*log(model->noise->etaA);
      logL +=  2.0*(double)data->N*log(model->noise->etaE);
      break;
    default:
      fprintf(stderr,"Unsupported number of channels in gaussian_log_likelihood()\n");
      exit(1);
  }

  free_tdi(residual);
  return logL;
}

