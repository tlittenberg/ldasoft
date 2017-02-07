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

#include "GalacticBinary.h"
#include "GalacticBinaryPrior.h"
#include "GalacticBinaryProposal.h"

double draw_from_prior(struct Model *model, double *params, gsl_rng *seed)
{
  for(int n=0; n<8; n++) params[n] = model->prior[n][0] + gsl_rng_uniform(seed)*(model->prior[n][1]-model->prior[n][0]);
  return model->logPriorVolume;
}
