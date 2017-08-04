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

/*
 proposal prototype
 (func*)(struct, struct, struct, double, gsl_rng)
 */
struct Proposal
{
  double (*function)(struct Data*,struct Model*,struct Source*,struct Proposal*,double*,gsl_rng*);
  int *trial;
  int *accept;
  char name[128];
  double weight; /* between 0 and 1 */

  int size;
  double *vector;
  double **matrix;
  double ***tensor;
};

void setup_frequency_proposal(struct Data *data);
void print_acceptance_rates(struct Proposal **proposal, int NP, int ic, FILE *fptr);
double draw_from_spectrum(struct Data *data, struct Model *model, struct Source *source, UNUSED struct Proposal *proposal, double *params, gsl_rng *seed);
double draw_from_prior(UNUSED struct Data *data, struct Model *model, UNUSED struct Source *source, struct Proposal *proposal, double *params, gsl_rng *seed);
double draw_from_extrinsic_prior(UNUSED struct Data *data, struct Model *model, UNUSED struct Source *source, struct Proposal *proposal, double *params, gsl_rng *seed);
double draw_from_fisher(UNUSED struct Data *data, struct Model *model, struct Source *source, struct Proposal *proposal, double *params, gsl_rng *seed);
double draw_from_cdf(UNUSED struct Data *data, struct Model *model, struct Source *source, struct Proposal *proposal, double *params, gsl_rng *seed);
double fm_shift(struct Data *data, struct Model *model, struct Source *source, struct Proposal *proposal, double *params, gsl_rng *seed);
double t0_shift(UNUSED struct Data *data, struct Model *model, UNUSED struct Source *source, UNUSED struct Proposal *proposal, UNUSED double *params, gsl_rng *seed);

double cdf_density(struct Model *model, struct Source *source, struct Proposal *proposal);

void initialize_proposal(struct Orbit *orbit, struct Data *data, struct Chain *chain, struct Flags *flags, struct Proposal **proposal, int NMAX);

void setup_fstatistic_proposal(struct Orbit *orbit, struct Data *data, struct Flags *flags, struct Proposal *proposal);

#endif /* GalacticBinaryProposal_h */
