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

#define FIXME 0

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

void alloc_data(struct Data **data_vec, int NMCMC, int NMAX)
{
  for(int n=0; n<NMAX; n++)
  {
    struct Data *data = data_vec[n];
    
    data->tdi   = malloc(sizeof(struct TDI*)    * data->Nsegment);
    data->noise = malloc(sizeof(struct Noise*)  * data->Nsegment);
    data->inj   = malloc(sizeof(struct Source*) * data->Nsegment);

    for(int s=0; s<data->Nsegment; s++)
    {
      data->tdi[s]   = malloc(sizeof(struct TDI));
      data->noise[s] = malloc(sizeof(struct Noise));
      data->inj[s]   = malloc(sizeof(struct Source));

      alloc_tdi(data->tdi[s], data->N, data->Nchannel);
      alloc_noise(data->noise[s], data->N);
      alloc_source(data->inj[s],data->N,data->Nchannel);
    }
    
    //reconstructed signal model
    int n_re,n_im;
    data->h_rec = malloc(data->N*2*sizeof(double **));
    data->h_res = malloc(data->N*sizeof(double **));
    data->h_pow = malloc(data->N*sizeof(double **));
    data->S_pow = malloc(data->N*sizeof(double **));
    
    //number of waveform samples to save
    data->Nwave=100;
    
    //downsampling rate of post-burn-in samples
    data->downsample = NMCMC/data->Nwave;
    
    for(int n=0; n<data->N; n++)
    {
      n_re = 2*n;
      n_im = n_re+1;
      
      data->S_pow[n]    = malloc(data->Nchannel*sizeof(double *));
      data->h_pow[n]    = malloc(data->Nchannel*sizeof(double *));
      data->h_res[n]    = malloc(data->Nchannel*sizeof(double *));
      data->h_rec[n_re] = malloc(data->Nchannel*sizeof(double *));
      data->h_rec[n_im] = malloc(data->Nchannel*sizeof(double *));
      for(int l=0; l<data->Nchannel; l++)
      {
        data->S_pow[n][l]    = malloc(data->Nwave*sizeof(double));
        data->h_pow[n][l]    = malloc(data->Nwave*sizeof(double));
        data->h_res[n][l]    = malloc(data->Nwave*sizeof(double));
        data->h_rec[n_re][l] = malloc(data->Nwave*sizeof(double));
        data->h_rec[n_im][l] = malloc(data->Nwave*sizeof(double));
        for(int m=0; m<data->Nwave; m++)
        {
          data->h_rec[n_re][l][m] = 0.0;
          data->h_rec[n_im][l][m] = 0.0;
          data->h_res[n][l][m]    = 0.0;
          data->h_pow[n][l][m]    = 0.0;
          data->S_pow[n][l][m]    = 0.0;
        }
      }
    }
  }
}

void initialize_chain(struct Chain *chain, struct Flags *flags, long *seed, int NC)
{
  int ic;
  chain->NC = NC;
  chain->index = malloc(NC*sizeof(int));
  chain->acceptance = malloc(NC*sizeof(double));
  chain->temperature = malloc(NC*sizeof(double));
  chain->avgLogL     = malloc(NC*sizeof(double));
  for(ic=0; ic<NC; ic++)
  {
    chain->index[ic]=ic;
    chain->acceptance[ic] = 1.0;
    chain->temperature[ic] = pow(1.2,(double)ic);
    chain->avgLogL[ic] = 0.0;
  }
  //set hottest chain to ~infinite temperature
  chain->temperature[NC-1] = 1e6;
  chain->logLmax = 0.0;
  
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

  chain->likelihoodFile = fopen("log_likelihood_chain.dat","w");
  
  chain->temperatureFile = fopen("temperature_chain.dat","w");

  chain->parameterFile = malloc(NC*sizeof(FILE *));
  chain->parameterFile[0] = fopen("parameter_chain.dat.0","w");

  chain->noiseFile = malloc(NC*sizeof(FILE *));
  chain->noiseFile[0] = fopen("noise_chain.dat.0","w");

  if(flags->verbose)
  {
    char filename[1024];
    for(ic=1; ic<NC; ic++)
    {
      sprintf(filename,"parameter_chain.dat.%i",ic);
      chain->parameterFile[ic] = fopen(filename,"w");

      sprintf(filename,"noise_chain.dat.%i",ic);
      chain->noiseFile[ic] = fopen(filename,"w");
    }
  }
}

void free_chain(struct Chain *chain, struct Flags *flags)
{
  free(chain->index);
  free(chain->temperature);
  free(chain->acceptance);
  free(chain->avgLogL);
  
  for(int ic=0; ic<chain->NC; ic++) gsl_rng_free(chain->r[ic]);
  free(chain->r);
  free(chain->T);
  
  fclose(chain->likelihoodFile);
  fclose(chain->parameterFile[0]);
  if(flags->verbose)
  {
    for(int ic=1; ic<chain->NC; ic++) fclose(chain->parameterFile[ic]);
  }
  free(chain->parameterFile);
  
  free(chain);
}

void alloc_model(struct Model *model, int Nmax, int NFFT, int Nchannel)
{
  model->Nlive  = 1;
  model->Nmax   = Nmax;
  model->source = malloc(model->Nmax*sizeof(struct Source *));
  model->noise  = malloc(sizeof(struct Noise));
  model->tdi    = malloc(sizeof(struct TDI));
  model->t0     = 0.0;

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
  //Source parameters
  copy->Nmax           = origin->Nmax;
  copy->Nlive          = origin->Nlive;
  for(int n=0; n<origin->Nmax; n++) copy_source(origin->source[n],copy->source[n]);
  
  //Noise parameters
  copy_noise(origin->noise,copy->noise);
  
  //TDI
  copy_tdi(origin->tdi,copy->tdi);

  //Start time for segment for model
  copy->t0 = origin->t0;
  copy->t0_min = origin->t0_min;
  copy->t0_max = origin->t0_max;

  //Source parameter priors
  for(int n=0; n<8; n++) for(int j=0; j<2; j++) copy->prior[n][j] = origin->prior[n][j];
  copy->logPriorVolume = origin->logPriorVolume;

  //Model likelihood
  copy->logL           = origin->logL;
  copy->logLnorm       = origin->logLnorm;
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
  noise->N    = NFFT;
  
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
  
  for(int n=0; n<origin->N; n++)
  {
    copy->SnX[n] = origin->SnX[n];
    copy->SnA[n] = origin->SnA[n];
    copy->SnE[n] = origin->SnE[n];
  }
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
  
  
  //Response
  copy_tdi(origin->tdi,copy->tdi);
  
  //FIsher
  for(int i=0; i<NP; i++)
  {
    for(int j=0; j<NP; j++)
    {
      copy->fisher_matrix[i][j] = origin->fisher_matrix[i][j];
      copy->fisher_evectr[i][j] = origin->fisher_evectr[i][j];
    }
    copy->fisher_evalue[i] = origin->fisher_evalue[i];
    copy->params[i]        = origin->params[i];
  }
}

void free_source(struct Source *source)
{
  int NP=8;
  for(int i=0; i<NP; i++)
  {
    free(source->fisher_matrix[i]);
    free(source->fisher_evectr[i]);
  }
  free(source->fisher_matrix);
  free(source->fisher_evectr);
  free(source->fisher_evalue);
  free(source->params);
  
  free_tdi(source->tdi);
  
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

void generate_signal_model(struct Orbit *orbit, struct Data *data, struct Model *model)
{
  int i,j,n;
  int N2=data->N*2;

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
    galactic_binary(orbit, data->T, model->t0, source->params, source->tdi->X, source->tdi->A, source->tdi->E, source->BW, source->tdi->Nchannel);
    
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
}

void generate_noise_model(struct Data *data, struct Model *model)
{
  switch(data->Nchannel)
  {
    case 1:
      for(int n=0; n<data->N; n++) model->noise->SnX[n] = data->noise[FIXME]->SnX[n]*model->noise->etaX;
      break;
    case 2:
      for(int n=0; n<data->N; n++)
      {
        model->noise->SnA[n] = data->noise[FIXME]->SnA[n]*model->noise->etaA;
        model->noise->SnE[n] = data->noise[FIXME]->SnE[n]*model->noise->etaE;
      }
      break;
    default:
      break;
  }
}

double gaussian_log_likelihood(struct Orbit *orbit, struct Data *data, struct Model *model, int n)
{
  
  int N2=data->N*2;
  

  /**************************/
  /*                        */
  /* Form residual and sum  */
  /*                        */
  /**************************/
  
  struct TDI *residual = malloc(sizeof(struct TDI));
  alloc_tdi(residual, data->N, data->Nchannel);
  
  for(int i=0; i<N2; i++)
  {
    residual->X[i] = data->tdi[n]->X[i] - model->tdi->X[i];
    residual->A[i] = data->tdi[n]->A[i] - model->tdi->A[i];
    residual->E[i] = data->tdi[n]->E[i] - model->tdi->E[i];
  }

  double logL = 0.0;
  switch(data->Nchannel)
  {
    case 1:
      logL += -0.5*fourier_nwip(residual->X, residual->X, model->noise->SnX, data->N);
      break;
    case 2:
      logL += -0.5*fourier_nwip(residual->A, residual->A, model->noise->SnA, data->N);
      logL += -0.5*fourier_nwip(residual->E, residual->E, model->noise->SnE, data->N);
      break;
    default:
      fprintf(stderr,"Unsupported number of channels in gaussian_log_likelihood()\n");
      exit(1);
  }

  free_tdi(residual);

  return logL;
}

double gaussian_log_likelihood_constant_norm(struct Data *data, struct Model *model)
{
  
  double logLnorm = 0.0;
  
  switch(data->Nchannel)
  {
    case 1:
      logLnorm -= (double)data->N*log(model->noise->etaX);
      break;
    case 2:
      logLnorm -= (double)data->N*log(model->noise->etaA);
      logLnorm -= (double)data->N*log(model->noise->etaE);
      break;
    default:
      fprintf(stderr,"Unsupported number of channels in gaussian_log_likelihood()\n");
      exit(1);
  }
  
  return logLnorm;
}

double gaussian_log_likelihood_model_norm(struct Data *data, struct Model *model)
{
  
  double logLnorm = 0.0;
  
  switch(data->Nchannel)
  {
    case 1:
      for(int n=0; n<data->N; n++) logLnorm -= log(model->noise->SnX[n]);
      break;
    case 2:
      for(int n=0; n<data->N; n++)
      {
        logLnorm -= log(model->noise->SnA[n]);
        logLnorm -= log(model->noise->SnE[n]);
      }
      break;
    default:
      fprintf(stderr,"Unsupported number of channels in gaussian_log_likelihood()\n");
      exit(1);
  }
  
  return logLnorm;
}

void update_max_log_likelihood(struct Model ****model, struct Chain *chain, struct Flags *flags)
{
  int n = chain->index[0];
  int N = flags->injection;
  
  double logL = 0.0;
  double dlogL= 0.0;
  
  // get full likelihood
  for(int i=0; i<flags->injection; i++) for(int j=0; j<flags->segment; j++) logL += model[chain->index[0]][i][j]->logL + model[chain->index[0]][i][j]->logLnorm;
  
  // update max
  if(logL > chain->logLmax)
  {
    dlogL = logL - chain->logLmax;
    chain->logLmax = logL;
    
    //clone chains if new mode is found (dlogL > D/2)
    if( dlogL > (double)(8*N/2) )
    {
      for(int ic=1; ic<chain->NC; ic++)
      {
        for(int i=0; i<N; i++) for(int j=0; j<flags->segment; j++) copy_model(model[n][i][j],model[chain->index[ic]][i][j]);
      }
    }
  }
}

