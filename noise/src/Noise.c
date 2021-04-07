//
//  Noise.c
//
//
//  Created by Tyson Littenberg on 4/05/21.
//


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <omp.h>

#include <LISA.h>
#include <GalacticBinary.h>
#include <GalacticBinaryIO.h>
#include <GalacticBinaryMath.h>
#include <GalacticBinaryData.h>
#include <GalacticBinaryModel.h>

#include "Noise.h"


void alloc_spline_model(struct SplineModel *model, int Ndata, int Nspline)
{
    model->psd = malloc(sizeof(struct Noise));
    model->spline=malloc(sizeof(struct Noise));
    
    alloc_noise(model->psd, Ndata);
    alloc_noise(model->spline, Nspline);
}

void free_spline_model(struct SplineModel *model)
{
    free_noise(model->psd);
    free_noise(model->spline);
    free(model);
}

void copy_spline_model(struct SplineModel *origin, struct SplineModel *copy)
{
    //Spline parameters
    copy_noise(origin->spline,copy->spline);
    
    //Noise model parameters
    copy_noise(origin->psd,copy->psd);
    
    //Model likelihood
    copy->logL = origin->logL;
}

void CubicSplineGSL(int N, double *x, double *y, int Nint, double *xint, double *yint)
{
    int n;
    
    /* set up GSL cubic spline */
    gsl_spline       *cspline = gsl_spline_alloc(gsl_interp_cspline, N);
    gsl_interp_accel *acc    = gsl_interp_accel_alloc();
    
    /* get derivatives */
    gsl_spline_init(cspline,x,y,N);
    
    for(n=0; n<Nint; n++) yint[n]=gsl_spline_eval(cspline,xint[n],acc);
    
    
    gsl_spline_free (cspline);
    gsl_interp_accel_free (acc);
    
}

void generate_spline_noise_model(struct SplineModel *model)
{
    struct Noise *psd = model->psd;
    struct Noise *spline = model->spline;
    
    CubicSplineGSL(spline->N, spline->f, spline->SnA, psd->N, psd->f, psd->SnA);
    CubicSplineGSL(spline->N, spline->f, spline->SnE, psd->N, psd->f, psd->SnE);
}

double noise_log_likelihood(struct Data *data, struct SplineModel *model)
{
    double logL = 0.0;
    
    struct TDI *tdi = data->tdi[0];
    struct Noise *psd = model->psd;
    
    int N = data->N;
    
    logL += -0.5*fourier_nwip(tdi->A, tdi->A, psd->SnA, N);
    logL += -0.5*fourier_nwip(tdi->E, tdi->E, psd->SnE, N);
    
    for(int n=0; n<N; n++)
    {
        logL -= log(psd->SnA[n]);
        logL -= log(psd->SnE[n]);
    }
    
    return logL;
}

void noise_spline_model_mcmc(struct Orbit *orbit, struct Data *data, struct SplineModel *model, struct SplineModel *trial, struct Chain *chain, struct Flags *flags, int ic)
{
    double logH  = 0.0; //(log) Hastings ratio
    double loga  = 1.0; //(log) transition probability
    
    double logPx  = 0.0; //(log) prior density for model x (current state)
    double logPy  = 0.0; //(log) prior density for model y (proposed state)
    
    //shorthand pointers
    struct SplineModel *model_x = model;
    struct SplineModel *model_y = trial;
    
    //copy current state into trial
    copy_spline_model(model_x, model_y);
    
    //pick a point any point
    int k = (int)floor(gsl_rng_uniform(chain->r[0])*(double)model_y->spline->N);
    
    double Sn = AEnoise_FF(orbit->L, orbit->fstar, model_y->spline->f[k]);
    model_y->spline->SnA[k] += 0.2*Sn*gsl_ran_gaussian(chain->r[0],1);
    model_y->spline->SnE[k] += 0.2*Sn*gsl_ran_gaussian(chain->r[0],1);
    
    //compute spline
    if(!flags->prior)
    {
        //generate_noise_model(data, model_y);
        generate_spline_noise_model(model_y);
        
        //compute likelihood
        model_y->logL = noise_log_likelihood(data, model_y);
        
        /*
         H = [p(d|y)/p(d|x)]/T x p(y)/p(x) x q(x|y)/q(y|x)
         */
        logH += (model_y->logL - model_x->logL)/chain->temperature[ic]; //delta logL
    }
    logH += logPy - logPx; //priors
    
    loga = log(gsl_rng_uniform(chain->r[ic]));
    if(logH > loga)
    {
        copy_spline_model(model_y, model_x);
    }
}

void initialize_spline_model(struct Orbit *orbit, struct Data *data, struct SplineModel *model, int Nspline)
{
    
    /* Initialize data models */
    alloc_spline_model(model, data->N, Nspline);
    
    //set up psd frequency grid
    for(int n=0; n<model->psd->N; n++)
    {
        double f = data->fmin + (double)n/data->T;
        model->psd->f[n] = f;
    }
    
    //divide into Nspline control points
    double df = (data->fmax - data->fmin)/(Nspline-1);
    
    for(int i=0; i<Nspline; i++)
    {
        double f = data->fmin + (double)i*df;
        model->spline->f[i] = f;
        
        model->spline->SnA[i] = AEnoise_FF(orbit->L, orbit->fstar, f);
        model->spline->SnE[i] = AEnoise_FF(orbit->L, orbit->fstar, f);
    }
    //shift first spline control point by half a bin to avoid rounding problems
    model->spline->f[0] -= 0.5/data->T;
    
    generate_spline_noise_model(model);
    model->logL = noise_log_likelihood(data, model);
}

void print_noise_model(struct Noise *noise, char filename[])
{
    FILE *fptr = fopen(filename,"w");
    for(int i=0; i<noise->N; i++)
    {
        fprintf(fptr,"%lg ",noise->f[i]);
        fprintf(fptr,"%lg ",noise->SnA[i]);
        fprintf(fptr,"%lg ",noise->SnE[i]);
        fprintf(fptr,"\n");
    }
    fclose(fptr);
}

