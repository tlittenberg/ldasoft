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
#include "GalacticBinary.h"
#include "GalacticBinaryPrior.h"
#include "GalacticBinaryWaveform.h"
#include "GalacticBinaryProposal.h"

double draw_from_prior(struct Model *model, double *params, gsl_rng *seed)
{
  for(int n=0; n<8; n++) params[n] = model->prior[n][0] + gsl_rng_uniform(seed)*(model->prior[n][1]-model->prior[n][0]);
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
  for(i=0; i<NP; i++) for (j=0; j<NP; j++) jump[j] += Amps[i]*source->fisher_evectr[j][i];
  
  //choose one eigenvector to jump along
  //i = (int)(gsl_rng_uniform(seed)*(double)NF);
  //for (j=0; j<NF; j++) jump[j] += Amps[i]*fisher->evector[b][j][i];
  
  //jump from current position
  for(i=0; i<NP; i++) params[i] = source->params[i] + jump[i];
  
  //ensure all new parameters are in range
  if(evaluate_uniform_prior(model, params)>-INFINITY) return 0.0;
  else return -INFINITY;
}
