//
//  GalacticBinaryProposal.h
//  
//
//  Created by Littenberg, Tyson B. (MSFC-ZP12) on 2/6/17.
//
//

#ifndef GalacticBinaryProposal_h
#define GalacticBinaryProposal_h

#include <stdio.h>

double draw_from_prior(struct Model *model, double *params, gsl_rng *seed);

#endif /* GalacticBinaryProposal_h */
