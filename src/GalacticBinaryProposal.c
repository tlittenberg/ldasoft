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


double draw_from_prior(struct Model *model, UNUSED struct Source *source, double *params, gsl_rng *seed)
{
  for(int n=0; n<8; n++) params[n] = model->prior[n][0] + gsl_rng_uniform(seed)*(model->prior[n][1]-model->prior[n][0]);
  
  for(int j=0; j<8; j++)
  {
    if(params[j]!=params[j]) fprintf(stderr,"draw_from_prior: params[%i]=%g, U[%g,%g]\n",j,params[j],model->prior[j][0],model->prior[j][1]);
  }

  return model->logPriorVolume;
}

double draw_from_fisher(struct Model *model, struct Source *source, double *params, gsl_rng *seed)
{
  int i,j;
  int NP=8;
  double Amps[NP];
  double jump[NP];
  
  //draw the eigen-jump amplitudes from N[0,1] scaled by evalue & dimension
  for(i=0; i<NP; i++)
  {
    Amps[i] = gsl_ran_gaussian(seed,1)/sqrt(source->fisher_evalue[i]);
    jump[i] = 0.0;
  }
  
  //decompose eigenjumps into paramter directions
  for(i=0; i<NP; i++) for (j=0; j<NP; j++)
  {
    jump[j] += Amps[i]*source->fisher_evectr[j][i];
    if(jump[j]!=jump[j])jump[j]=0.0;
  }
  
  //choose one eigenvector to jump along
  //i = (int)(gsl_rng_uniform(seed)*(double)NF);
  //for (j=0; j<NF; j++) jump[j] += Amps[i]*fisher->evector[b][j][i];
  
  //jump from current position
  for(i=0; i<NP; i++) params[i] = source->params[i] + jump[i];
  
  for(int j=0; j<8; j++)
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
  draw_from_fisher(model, source, params, seed);
  
  //perturb frequency by 1 fm
  double sign = -1. + 2*gsl_rng_uniform(seed);
  if(sign<0) sign = -1.0;
  else       sign =  1.0;
  
  double scale = 1.0 + floor(4.*gsl_rng_uniform(seed));
  
  params[0] += sign*scale*fm;
  
  //fm shift is symmetric
  return 0.0;
}
