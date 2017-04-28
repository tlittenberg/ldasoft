//
//  GalacticBinaryProposal.c
//  
//
//  Created by Littenberg, Tyson B. (MSFC-ZP12) on 2/6/17.
//
//

#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "LISA.h"
#include "Constants.h"
#include "GalacticBinary.h"
#include "GalacticBinaryPrior.h"
#include "GalacticBinaryWaveform.h"
#include "GalacticBinaryProposal.h"


void setup_frequency_proposal(struct Data *data)
{
  int BW = 20;
  double *power = data->p;
  double total  = 0.0;
  FILE *temp = fopen("temp2.dat","w");
  for(int i=0; i<data->N-BW; i++)
  {
    power[i]=0.0;
    for(int n=i; n<i+BW; n++)
    {
      double SnA = data->noise->SnA[n];
      double SnE = data->noise->SnE[n];
      
      double AA = data->tdi->A[2*n]*data->tdi->A[2*n]+data->tdi->A[2*n+1]*data->tdi->A[2*n+1];
      double EE = data->tdi->E[2*n]*data->tdi->E[2*n]+data->tdi->E[2*n+1]*data->tdi->E[2*n+1];
      
      power[i] += AA/SnA + EE/SnE;
      total += power[i];
    }
  }
  for(int i=data->N-BW; i<data->N; i++)
  {
    power[i] = power[data->N-BW-1];
    total += power[i];
  }
  
  data->pmax = 0.0;
  for(int i=0; i<data->N; i++)
  {
    fprintf(temp,"%i %lg\n",i,power[i]);
    if(power[i]>data->pmax) data->pmax = power[i];
  }
  fclose(temp);
  
  
  //also get SNR^2 of data
  total = 0.0;
  for(int n=0; n<data->N; n++)
  {
      double SnA = data->noise->SnA[n];
      double SnE = data->noise->SnE[n];
      
      double AA = data->tdi->A[2*n]*data->tdi->A[2*n]+data->tdi->A[2*n+1]*data->tdi->A[2*n+1];
      double EE = data->tdi->E[2*n]*data->tdi->E[2*n]+data->tdi->E[2*n+1]*data->tdi->E[2*n+1];
      
      total += AA/SnA + EE/SnE;
  }

  data->SNR2 = total - data->N;
  printf("total=%g, N=%i\n",total, data->N);
  if(data->SNR2<0.0)data->SNR2=0.0;
  data->SNR2*=4.0;//why the factor of 4?
  printf("data-based SNR^2:  %g (%g)\n", data->SNR2, sqrt(data->SNR2));

}

double draw_from_spectrum(struct Data *data, struct Model *model, struct Source *source, double *params, gsl_rng *seed)
{
  //TODO: Work in amplitude
  
  //rejections ample for f
  int check = 1;
  double alpha;
  int q;
  int count=0;
  while(check)
  {
    params[0] = model->prior[0][0] + gsl_rng_uniform(seed)*(model->prior[0][1]-model->prior[0][0]);
    alpha     = gsl_rng_uniform(seed)*data->pmax;
    q = (int)(params[0]-data->qmin);
    if(alpha<data->p[q]) check = 0;
    count++;
  }
  
  //random draws for other parameters
  for(int n=1; n<8; n++) params[n] = model->prior[n][0] + gsl_rng_uniform(seed)*(model->prior[n][1]-model->prior[n][0]);

  return 0;
}

double draw_from_prior(UNUSED struct Data *data, struct Model *model, UNUSED struct Source *source, double *params, gsl_rng *seed)
{
  for(int n=0; n<source->NP; n++) params[n] = model->prior[n][0] + gsl_rng_uniform(seed)*(model->prior[n][1]-model->prior[n][0]);
  
  for(int j=0; j<source->NP; j++)
  {
    if(params[j]!=params[j]) fprintf(stderr,"draw_from_prior: params[%i]=%g, U[%g,%g]\n",j,params[j],model->prior[j][0],model->prior[j][1]);
  }

  return model->logPriorVolume;
}

double draw_from_extrinsic_prior(UNUSED struct Data *data, struct Model *model, UNUSED struct Source *source, double *params, gsl_rng *seed)
{
  for(int n=1; n<3; n++) params[n] = model->prior[n][0] + gsl_rng_uniform(seed)*(model->prior[n][1]-model->prior[n][0]);
  for(int n=4; n<source->NP; n++) params[n] = model->prior[n][0] + gsl_rng_uniform(seed)*(model->prior[n][1]-model->prior[n][0]);
  
  for(int j=0; j<8; j++)
  {
    if(params[j]!=params[j]) fprintf(stderr,"draw_from_prior: params[%i]=%g, U[%g,%g]\n",j,params[j],model->prior[j][0],model->prior[j][1]);
  }
  
  return model->logPriorVolume;
}

double draw_from_fisher(UNUSED struct Data *data, struct Model *model, struct Source *source, double *params, gsl_rng *seed)
{
  int i,j;
  int NP=source->NP;
  double sqNP = 2.82842712474619; //sqrt(8)
  double Amps[NP];
  double jump[NP];
  
  //draw the eigen-jump amplitudes from N[0,1] scaled by evalue & dimension
  for(i=0; i<NP; i++)
  {
    //Amps[i] = gsl_ran_gaussian(seed,1)/sqrt(source->fisher_evalue[i])/(double)NP;
    Amps[i] = gsl_ran_gaussian(seed,1)/sqrt(source->fisher_evalue[i])/sqNP;
    jump[i] = 0.0;
  }
  
  //decompose eigenjumps into paramter directions
  /*
  for(i=0; i<NP; i++) for (j=0; j<NP; j++)
  {
    jump[j] += Amps[i]*source->fisher_evectr[j][i];
    if(jump[j]!=jump[j])jump[j]=0.0;
  }*/
  
  //choose one eigenvector to jump along
  i = (int)(gsl_rng_uniform(seed)*(double)NP);
  for (j=0; j<NP; j++) jump[j] += Amps[i]*source->fisher_evectr[j][i];
  
  //jump from current position
  for(i=0; i<NP; i++) params[i] = source->params[i] + jump[i];
  
  for(int j=0; j<NP; j++)
  {
    if(params[j]!=params[j]) fprintf(stderr,"draw_from_fisher: params[%i]=%g, N[%g,%g]\n",j,params[j],source->params[j],jump[j]);
  }

  //not updating Fisher between moves, proposal is symmetric
  return 0.0;
}

double fm_shift(struct Data *data, struct Model *model, struct Source *source, double *params, gsl_rng *seed)
{
  //doppler modulation frequency (in bins)
  double fm = data->T/YEAR;
  
  //update all parameters with a draw from the fisher
  //draw_from_fisher(data, model, source, params, seed);
  for(int n=1; n<source->NP; n++) params[n] = model->prior[n][0] + gsl_rng_uniform(seed)*(model->prior[n][1]-model->prior[n][0]);

  //perturb frequency by 1 fm
  double scale = floor(6*gsl_ran_gaussian(seed,1));
  
  params[0] += scale*fm;
  //params[7] += scale*fm*fm;
  
  //fm shift is symmetric
  return 0.0;
}

double t0_shift(UNUSED struct Data *data, struct Model *model, UNUSED struct Source *source, UNUSED double *params, gsl_rng *seed)
{
  //uniform draw
  if(gsl_rng_uniform(seed) < 0.5 )
    model->t0 = model->t0_min + gsl_rng_uniform(seed)*(model->t0_max - model->t0_min);
  
  //gaussian draw
  else
    model->t0 += 3.0*gsl_ran_gaussian(seed,1);
  
  
  //t0 shift is symmetric
  if(model->t0 < model->t0_min || model->t0 >= model->t0_max) return -INFINITY;
  else return 0.0;
}
