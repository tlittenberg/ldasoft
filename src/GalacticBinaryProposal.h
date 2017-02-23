//
//  GalacticBinaryProposal.h
//  
//
//  Created by Littenberg, Tyson B. (MSFC-ZP12) on 2/6/17.
//
//

#ifndef GalacticBinaryProposal_h
#define GalacticBinaryProposal_h

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include <stdio.h>

double draw_from_prior(struct Model *model, UNUSED struct Source *source, double *params, gsl_rng *seed);
double draw_from_fisher(struct Model *model, struct Source *source, double *params, gsl_rng *seed);
double fm_shift(struct Data *data, struct Model *model, struct Source *source, double *params, gsl_rng *seed);

#endif /* GalacticBinaryProposal_h */
