//
//  Noise.h
//
//
//  Created by Tyson Littenberg on 4/05/21.
//


/**
 @file Noise.h
 \brief Sampling routines for Noise module
 
 Including
 - BayesLine spline-only
 */

#ifndef Noise_h
#define Noise_h

struct SplineModel
{
    int Nmin; ///<! Minimum number of spline control points
    int Nmax; ///<! Maximum number of spline control points
    double logL; ///<! Log likelhood of spline model
    struct Noise *spline; ///<! Spline control points
    struct Noise *psd; ///<! Reconstructed noise model
};

void alloc_spline_model(struct SplineModel *model, int Ndata, int Nspline);
void free_spline_model(struct SplineModel *model);
void copy_spline_model(struct SplineModel *origin, struct SplineModel *copy);
void print_spline_state(struct SplineModel *model, FILE *fptr, int step);

void CubicSplineGSL(int N, double *x, double *y, int Nint, double *xint, double *yint);
void generate_spline_noise_model(struct SplineModel *model);

double noise_log_likelihood(struct Data *data, struct SplineModel *model);

void spline_ptmcmc(struct SplineModel **model, struct Chain *chain, struct Flags *flags);
void noise_spline_model_mcmc(struct Orbit *orbit, struct Data *data, struct SplineModel *model, struct Chain *chain, struct Flags *flags, int ic);
void noise_spline_model_rjmcmc(struct Orbit *orbit, struct Data *data, struct SplineModel *model, struct Chain *chain, struct Flags *flags, int ic);

void initialize_spline_model(struct Orbit *orbit, struct Data *data, struct SplineModel *model, int Nspline);
void print_noise_model(struct Noise *noise, char filename[]);
#endif /* NoiseMCMC_h */
