//
//  GalacticBinaryPrior.h
//  
//
//  Created by Littenberg, Tyson B. (MSFC-ZP12) on 2/6/17.
//
//

#ifndef GalacticBinaryPrior_h
#define GalacticBinaryPrior_h

#include <stdio.h>

void set_uniform_prior(struct Model *model, struct Data *data);
double evaluate_uniform_prior(struct Model *model, double *params);

void setup_galaxy_prior(struct Flags *flags, double *skyhist, int Nth, int Nph);

#endif /* GalacticBinaryPrior_h */
