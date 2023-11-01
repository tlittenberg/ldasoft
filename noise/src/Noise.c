/*
 *  Copyright (C) 2021 Tyson B. Littenberg (MSFC-ST12)
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

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

#define MIN_SPLINE_STENCIL 5
#define MIN_SPLINE_SPACING 256.

void map_array_to_noise_params(struct InstrumentModel *model)
{
    model->sacc12 = model->sacc[0];
    model->sacc13 = model->sacc[1];
    model->sacc21 = model->sacc[2];
    model->sacc23 = model->sacc[3];
    model->sacc31 = model->sacc[4];
    model->sacc32 = model->sacc[5];
    
    model->soms12 = model->soms[0];
    model->soms13 = model->soms[1];
    model->soms21 = model->soms[2];
    model->soms23 = model->soms[3];
    model->soms31 = model->soms[4];
    model->soms32 = model->soms[5];
}

void map_noise_params_to_array(struct InstrumentModel *model)
{
    model->sacc[0] = model->sacc12;
    model->sacc[1] = model->sacc13;
    model->sacc[2] = model->sacc21;
    model->sacc[3] = model->sacc23;
    model->sacc[4] = model->sacc31;
    model->sacc[5] = model->sacc32;
    
    model->soms[0] = model->soms12;
    model->soms[1] = model->soms13;
    model->soms[2] = model->soms21;
    model->soms[3] = model->soms23;
    model->soms[4] = model->soms31;
    model->soms[5] = model->soms32;
}

void map_array_to_foreground_params(struct ForegroundModel *model)
{
    model->Amp   = exp(model->sgal[0]);
    model->f1    = exp(model->sgal[1]);
    model->alpha = model->sgal[2];
    model->fk    = exp(model->sgal[3]);
    model->f2    = exp(model->sgal[4]);
}

void map_foreground_params_to_array(struct ForegroundModel *model)
{
    model->sgal[0] = log(model->Amp);
    model->sgal[1] = log(model->f1);
    model->sgal[2] = model->alpha;
    model->sgal[3] = log(model->fk);
    model->sgal[4] = log(model->f2);
}

void alloc_spline_model(struct SplineModel *model, int Ndata, int Nchannel, int Nspline)
{
    model->psd = malloc(sizeof(struct Noise));
    model->spline=malloc(sizeof(struct Noise));
    
    alloc_noise(model->psd, Nchannel, Ndata);
    alloc_noise(model->spline, Nchannel, Nspline);
}

void alloc_instrument_model(struct InstrumentModel *model, int Ndata, int Nchannel)
{
    model->Nlink=6;
    model->soms = malloc(model->Nlink*sizeof(double));
    model->sacc = malloc(model->Nlink*sizeof(double));
    model->psd = malloc(sizeof(struct Noise));
    alloc_noise(model->psd, Ndata, Nchannel);
}

void alloc_foreground_model(struct ForegroundModel *model, int Ndata, int Nchannel)
{
    model->Nparams=5;
    model->sgal = malloc(model->Nparams*sizeof(double));
    model->psd = malloc(sizeof(struct Noise));
    alloc_noise(model->psd, Ndata, Nchannel);
}

void free_spline_model(struct SplineModel *model)
{
    free_noise(model->psd);
    free_noise(model->spline);
    free(model);
}

void free_instrument_model(struct InstrumentModel *model)
{
    free(model->soms);
    free(model->sacc);
    free_noise(model->psd);
    free(model);
}

void free_foreground_model(struct ForegroundModel *model)
{
    free_noise(model->psd);
    free(model->sgal);
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
    
    //Priors
    copy->Nmin = origin->Nmin;
    copy->Nmax = origin->Nmax;
}

void copy_instrument_model(struct InstrumentModel *origin, struct InstrumentModel *copy)
{
    //Noise model parameters
    copy_noise(origin->psd,copy->psd);

    //Model likelihood
    copy->logL = origin->logL;

    //Instrument model Parameters
    copy->Nlink = origin->Nlink;
    memcpy(copy->soms, origin->soms, origin->Nlink*sizeof(double));
    memcpy(copy->sacc, origin->sacc, origin->Nlink*sizeof(double));
}

void copy_foreground_model(struct ForegroundModel *origin, struct ForegroundModel *copy)
{
    //Noise model parameters
    copy_noise(origin->psd,copy->psd);
    
    //Model likelihood
    copy->logL = origin->logL;

    //Foreground model Parameters
    copy->Tobs = origin->Tobs;
    copy->Nparams = origin->Nparams;
    memcpy(copy->sgal, origin->sgal, origin->Nparams*sizeof(double));
    
}

void print_spline_state(struct SplineModel *model, FILE *fptr, int step)
{
    fprintf(fptr,"%i %.12g %i\n",step,model->logL,model->spline->N);
}

void print_instrument_state(struct InstrumentModel *model, FILE *fptr, int step)
{
    fprintf(fptr,"%i %.12g ",step,model->logL);
    for(int i=0; i<model->Nlink; i++)fprintf(fptr,"%.12g ", model->sacc[i]);
    for(int i=0; i<model->Nlink; i++)fprintf(fptr,"%.12g ", model->soms[i]);
    fprintf(fptr,"\n");
}

void print_foreground_state(struct ForegroundModel *model, FILE *fptr, int step)
{
    fprintf(fptr,"%i %.12g ",step,model->logL);
    for(int i=0; i<model->Nparams; i++)fprintf(fptr,"%.12g ", model->sgal[i]);
    fprintf(fptr,"\n");
}


void update_spline_noise_model(struct SplineModel *model, int new_knot, int min_knot, int max_knot)
{
    struct Noise *psd = model->psd;
    struct Noise *spline = model->spline;
    
    /* find location in data vector of knot */
    double T = 1./(psd->f[1] - psd->f[0]);

    gsl_spline *cspline_A = gsl_spline_alloc(gsl_interp_akima, spline->N);
    gsl_spline *cspline_E = gsl_spline_alloc(gsl_interp_akima, spline->N);
    gsl_interp_accel *acc_A = gsl_interp_accel_alloc();
    gsl_interp_accel *acc_E = gsl_interp_accel_alloc();

    /* have to recompute the spline everywhere (derivatives on boundary) */
    gsl_spline_init(cspline_A,spline->f,spline->C[0][0],spline->N);
    gsl_spline_init(cspline_E,spline->f,spline->C[1][1],spline->N);

    
    int imin = (int)((spline->f[min_knot]-psd->f[0])*T);
    int imax = (int)((spline->f[max_knot]-psd->f[0])*T);
    
    
    for(int n=imin; n<imax; n++)
    {
        psd->C[0][0][n]=gsl_spline_eval(cspline_A,psd->f[n],acc_A);
        psd->C[1][1][n]=gsl_spline_eval(cspline_E,psd->f[n],acc_E);
        
        /*
         apply transfer function
         -this catches the sharp features in the spectrum from f/fstar
         -without needing to interpolate
         */
        
        psd->C[0][0][n]+=psd->transfer[n];
        psd->C[1][1][n]+=psd->transfer[n];
        
    }
    invert_noise_covariance_matrix(psd);

    gsl_spline_free (cspline_A);
    gsl_spline_free (cspline_E);
    gsl_interp_accel_free (acc_A);
    gsl_interp_accel_free (acc_E);

}


void generate_spline_noise_model(struct SplineModel *model)
{
    struct Noise *psd = model->psd;
    struct Noise *spline = model->spline;
    
    CubicSplineGSL(spline->N, spline->f, spline->C[0][0], psd->N, psd->f, psd->C[0][0]);
    CubicSplineGSL(spline->N, spline->f, spline->C[1][1], psd->N, psd->f, psd->C[1][1]);

    //set a floor on Sn so likelihood doesn't go crazy
    for(int n=0; n<psd->N; n++)
    {
        /*
         apply transfer function
         -this catches the sharp features in the spectrum from f/fstar
         -without needing to interpolate
         */
        psd->C[0][0][n]+=psd->transfer[n];
        psd->C[1][1][n]+=psd->transfer[n];
        
    }
    invert_noise_covariance_matrix(psd);
}

void generate_instrument_noise_model(struct Data *data, struct Orbit *orbit, struct InstrumentModel *model)
{
    double f;
    double x;
    double cosx;
    double oms_transfer_function;
    double acc_transfer_function;
    double tdi_transfer_function;
    double Sacc;
    double Soms;
    double Sacc12,Sacc21,Sacc13,Sacc31,Sacc23,Sacc32;
    double Soms12,Soms21,Soms13,Soms31,Soms23,Soms32;

    map_array_to_noise_params(model);

    for(int n=0; n<data->N; n++)
    {
        f = model->psd->f[n];
        x = f/orbit->fstar;
        cosx = cos(x);
        acc_transfer_function = 1./(PI2*f*CLIGHT)/(PI2*f*CLIGHT) * (1.0 + pow(0.4e-3/f,2)) * (1.0 + pow(f/8.0e-3,4));
        oms_transfer_function = (PI2*f/CLIGHT)*(PI2*f/CLIGHT) * (1.0 + pow(2.0e-3/f,4));
        tdi_transfer_function = noise_transfer_function(x);
        
        switch(data->Nchannel)
        {
            case 1:
                Sacc = model->sacc12 * acc_transfer_function;
                Soms = model->soms12 * oms_transfer_function;
                model->psd->C[0][0][n] = XYZnoise_FF(orbit->L, orbit->fstar, f, Sacc, Soms);
                break;
            case 2:
                Sacc = model->sacc12 * acc_transfer_function;
                Soms = model->soms12 * oms_transfer_function;
                model->psd->C[0][0][n] = AEnoise_FF(orbit->L, orbit->fstar, f, Sacc, Soms);
                model->psd->C[1][1][n] = AEnoise_FF(orbit->L, orbit->fstar, f, Sacc, Soms);
                model->psd->C[0][1][n] = 0.0;
                break;
            case 3:
                
                Sacc12 = model->sacc12 * acc_transfer_function;
                Sacc21 = model->sacc21 * acc_transfer_function;
                Sacc23 = model->sacc23 * acc_transfer_function;
                Sacc32 = model->sacc32 * acc_transfer_function;
                Sacc13 = model->sacc13 * acc_transfer_function;
                Sacc31 = model->sacc31 * acc_transfer_function;
                
                Soms12 = model->soms12 * oms_transfer_function;
                Soms21 = model->soms21 * oms_transfer_function;
                Soms23 = model->soms23 * oms_transfer_function;
                Soms32 = model->soms32 * oms_transfer_function;
                Soms13 = model->soms13 * oms_transfer_function;
                Soms31 = model->soms31 * oms_transfer_function;
                
                
                //OMS noise
                model->psd->C[0][0][n] = (Soms12+Soms21+Soms13+Soms31)*.25;
                model->psd->C[1][1][n] = (Soms23+Soms32+Soms21+Soms12)*.25;
                model->psd->C[2][2][n] = (Soms31+Soms13+Soms32+Soms23)*.25;

                model->psd->C[0][1][n] = (Soms12 + Soms21)*.5;
                model->psd->C[0][2][n] = (Soms13 + Soms31)*.5;
                model->psd->C[1][2][n] = (Soms32 + Soms23)*.5;
                
                //Acceleration noise
                model->psd->C[0][0][n] += 2. * ( (Sacc12 + Sacc13)*.5 + ((Sacc21 + Sacc31)*.5)*cosx*cosx );
                model->psd->C[1][1][n] += 2. * ( (Sacc23 + Sacc21)*.5 + ((Sacc32 + Sacc12)*.5)*cosx*cosx );
                model->psd->C[2][2][n] += 2. * ( (Sacc31 + Sacc32)*.5 + ((Sacc13 + Sacc23)*.5)*cosx*cosx );

                model->psd->C[0][1][n] += 4. * ( (Sacc12 + Sacc21)*.5 );
                model->psd->C[0][2][n] += 4. * ( (Sacc13 + Sacc31)*.5 );
                model->psd->C[1][2][n] += 4. * ( (Sacc32 + Sacc23)*.5 );

                //TDI transfer functions
                model->psd->C[0][0][n] *= 16. * tdi_transfer_function;
                model->psd->C[1][1][n] *= 16. * tdi_transfer_function;
                model->psd->C[2][2][n] *= 16. * tdi_transfer_function;

                model->psd->C[0][1][n] *= -8. * tdi_transfer_function * cosx ;
                model->psd->C[0][2][n] *= -8. * tdi_transfer_function * cosx ;
                model->psd->C[1][2][n] *= -8. * tdi_transfer_function * cosx ;
                
                //Symmetry
                model->psd->C[1][0][n] = model->psd->C[0][1][n];
                model->psd->C[2][0][n] = model->psd->C[0][2][n];
                model->psd->C[2][1][n] = model->psd->C[1][2][n];
                
                break;
        }
        
        //Normalization
        for(int i=0; i<data->Nchannel; i++)
            for(int j=0; j<data->Nchannel; j++)
                model->psd->C[i][j][n] *= 2.0;
    }
}

void generate_galactic_foreground_model(struct Data *data, struct Orbit *orbit, struct ForegroundModel *model)
{
    double f;
    double fk,f1;
    double Sgal;
    
    map_array_to_foreground_params(model);
    
    for(int n=0; n<data->N; n++)
    {
        f = model->psd->f[n];
        
        //Sgal = model->Amp*pow(f,5./3.) * exp(-pow(f/model->f1,model->alpha)) * ( 1. + tanh( (model->fk - f)/model->f2 ) );
        Sgal = model->Amp*pow(f,5./3.) * ( 1. + tanh( (model->fk - f)/model->f2 ) );

        switch(data->Nchannel)
        {
            case 1:
                model->psd->C[0][0][n] = Sgal;
                break;
            case 2:
                model->psd->C[0][0][n] = model->psd->C[1][1][n] = 1.5*Sgal;
                model->psd->C[0][1][n] = model->psd->C[1][0][n] = 0;
                break;
            case 3:
                model->psd->C[0][0][n] = model->psd->C[1][1][n] = model->psd->C[2][2][n] = Sgal;
                model->psd->C[0][1][n] = model->psd->C[0][2][n] = model->psd->C[1][2][n] = -0.5*Sgal;
                model->psd->C[1][0][n] = model->psd->C[2][0][n] = model->psd->C[2][1][n] = -0.5*Sgal;
                break;
        }
        
        //Normalization
        for(int i=0; i<data->Nchannel; i++)
            for(int j=0; j<data->Nchannel; j++)
                model->psd->C[i][j][n] *= 2.0;

    }
}

void generate_full_covariance_matrix(struct Noise *full, struct Noise *component, int Nchannel)
{
    for(int n=0; n<full->N; n++)
    {
        for(int i=0; i<Nchannel; i++)
        {
            for(int j=0; j<Nchannel; j++)
            {
                full->C[i][j][n] += component->C[i][j][n];
            }
        }
    }
}

double noise_log_likelihood(struct Data *data, struct Noise *noise)
{
    double logL = 0.0;
    
    struct TDI *tdi = data->tdi[0];
    
    int N = data->N;
    
    switch(data->Nchannel)
    {
        case 1:
            logL += -0.5*fourier_nwip(tdi->X, tdi->X, noise->invC[0][0], N);
            break;
        case 2:
            logL += -0.5*fourier_nwip(tdi->A, tdi->A, noise->invC[0][0], N);
            logL += -0.5*fourier_nwip(tdi->E, tdi->E, noise->invC[1][1], N);
            break;
        case 3:
            logL += -0.5*fourier_nwip(tdi->X, tdi->X, noise->invC[0][0], N);
            logL += -0.5*fourier_nwip(tdi->Y, tdi->Y, noise->invC[1][1], N);
            logL += -0.5*fourier_nwip(tdi->Z, tdi->Z, noise->invC[2][2], N);
            logL += -fourier_nwip(tdi->X, tdi->Y, noise->invC[0][1], N);
            logL += -fourier_nwip(tdi->X, tdi->Z, noise->invC[0][2], N);
            logL += -fourier_nwip(tdi->Y, tdi->Z, noise->invC[1][2], N);
            break;
    }
    for(int n=0; n<N; n++)
        logL -= 0.5*log(noise->detC[n]);
    
    return logL;
}

double noise_delta_log_likelihood(struct Data *data, struct SplineModel *model_x, struct SplineModel *model_y, double fmin, double fmax,int ic)
{
    double dlogL = 0.0;
    
    struct TDI *tdi = data->tdi[0];
    struct Noise *psd_x = model_x->psd;
    struct Noise *psd_y = model_y->psd;

    int N = (int)floor((fmax-fmin)*data->T);
    int imin = (int)floor((fmin-psd_x->f[0])*data->T);
    if(imin<0)imin=0;
    
    /* remove contribution for current state x */
    dlogL -= -0.5*fourier_nwip(tdi->A+2*imin, tdi->A+2*imin, psd_x->invC[0][0]+imin, N);
    dlogL -= -0.5*fourier_nwip(tdi->E+2*imin, tdi->E+2*imin, psd_x->invC[1][1]+imin, N);
    for(int n=imin; n<imin+N; n++)
        dlogL += log(psd_x->detC[n]);

    /* add contribution for proposed state y */
    dlogL += -0.5*fourier_nwip(tdi->A+2*imin, tdi->A+2*imin, psd_y->invC[0][0]+imin, N);
    dlogL += -0.5*fourier_nwip(tdi->E+2*imin, tdi->E+2*imin, psd_y->invC[1][1]+imin, N);
    for(int n=imin; n<imin+N; n++)
        dlogL -= log(psd_y->detC[n]);
    
    return dlogL;
}


void spline_ptmcmc(struct SplineModel **model, struct Chain *chain, struct Flags *flags)
{
    int a, b;
    int olda, oldb;
    
    double heat1, heat2;
    double logL1, logL2;
    double dlogL;
    double H;
    double alpha;
    double beta;
    
    int NC = chain->NC;
    
    for(b=NC-1; b>0; b--)
    {
        a = b - 1;
        chain->acceptance[a]=0;
        
        olda = chain->index[a];
        oldb = chain->index[b];
        
        heat1 = chain->temperature[a];
        heat2 = chain->temperature[b];
        
        logL1 = model[olda]->logL;
        logL2 = model[oldb]->logL;
        
        //Hot chains jump more rarely
        if(gsl_rng_uniform(chain->r[a])<1.0)
        {
            dlogL = logL2 - logL1;
            H  = (heat2 - heat1)/(heat2*heat1);
            
            alpha = exp(dlogL*H);
            beta  = gsl_rng_uniform(chain->r[a]);
            
            if(alpha >= beta)
            {
                chain->index[a] = oldb;
                chain->index[b] = olda;
                chain->acceptance[a]=1;
            }
        }
    }
}

void noise_ptmcmc(struct InstrumentModel **model, struct Chain *chain, struct Flags *flags)
{
    int a, b;
    int olda, oldb;
    
    double heat1, heat2;
    double logL1, logL2;
    double dlogL;
    double H;
    double alpha;
    double beta;
    
    int NC = chain->NC;
    
    for(b=NC-1; b>0; b--)
    {
        a = b - 1;
        chain->acceptance[a]=0;
        
        olda = chain->index[a];
        oldb = chain->index[b];
        
        heat1 = chain->temperature[a];
        heat2 = chain->temperature[b];
        
        logL1 = model[olda]->logL;
        logL2 = model[oldb]->logL;
        
        //Hot chains jump more rarely
        if(gsl_rng_uniform(chain->r[a])<1.0)
        {
            dlogL = logL2 - logL1;
            H  = (heat2 - heat1)/(heat2*heat1);
            
            alpha = exp(dlogL*H);
            beta  = gsl_rng_uniform(chain->r[a]);
            
            if(alpha >= beta)
            {
                chain->index[a] = oldb;
                chain->index[b] = olda;
                chain->acceptance[a]=1;
            }
        }
    }
}

static double uniform_frequency_draw(double fmin, double fmax, gsl_rng *r)
{
    return exp(log(fmin) + (log(fmax) - log(fmin))*gsl_rng_uniform(r));
}

static int check_frequency_spacing(double *f, int k, double T)
{
    if (
        (f[k]-f[k-1])*T < MIN_SPLINE_SPACING ||
        (f[k+1]-f[k])*T < MIN_SPLINE_SPACING
        ) return 1;
    else return 0;
}

void noise_spline_model_mcmc(struct Orbit *orbit, struct Data *data, struct SplineModel *model, struct Chain *chain, struct Flags *flags, int ic)
{
    double logH  = 0.0; //(log) Hastings ratio
    double loga  = 1.0; //(log) transition probability
    
    double logPx  = 0.0; //(log) prior density for model x (current state)
    double logPy  = 0.0; //(log) prior density for model y (proposed state)
    
    //shorthand pointers
    struct SplineModel *model_x = model;
    struct SplineModel *model_y = malloc(sizeof(struct SplineModel));
    alloc_spline_model(model_y, model_x->psd->N, data->Nchannel, model_x->spline->N);

    //alisases for pointers to frequency vectors
    double *fx = model_x->spline->f;
    double *fy = model_y->spline->f;
    
    //copy current state into trial
    copy_spline_model(model_x, model_y);
    
    //pick a point any point
    int k = (int)floor(gsl_rng_uniform(chain->r[ic])*(double)model_y->spline->N);

    //find the minimum and maximum indecies for stencil
    int half_stencil = (MIN_SPLINE_STENCIL-1)/2;
    int kmin = (k-half_stencil < 0) ? 0 : k-half_stencil;
    int kmax = (k+half_stencil > model_x->spline->N-1) ? model_x->spline->N-1 : k+half_stencil;
    
    //update frequency
    model_y->spline->f[k] = model_x->spline->f[k];
    if(k>0 && k<model_y->spline->N-1)
    {
        /* normal draw */
        if(gsl_rng_uniform(chain->r[ic])<0.8)
        {
            
            //get shortest distance between neighboring points
            double df_left = log(fx[k]) - log(fx[k-1]);
            double df_right= log(fx[k+1]) - log(fx[k]);
            double df = (df_left < df_right) ? df_left : df_right;
            
            //set shortest distance as 3-sigma jump for Gaussian draw
            double sigma = df/3.;
            
            //draw new frequency (illegally) checking that it stays between existing points
            do fy[k] = exp(log(fx[k]) + sigma*gsl_ran_gaussian(chain->r[ic],1));
            while( fy[k]<fy[k-1] || fy[k]>fy[k+1]);
        }
        else
        {
            /* uniform draw */
            fy[k] = uniform_frequency_draw(fx[k-1], fx[k+1], chain->r[ic]);
        }

        //check frequency prior
        if(check_frequency_spacing(fy, k, data->T)) logPy = -INFINITY;
        
    }

    //update amplitude
    double Sop, Spm;
    get_noise_levels("sangria", model_y->spline->f[k], &Spm, &Sop);
    double Sn = AEnoise_FF(orbit->L, orbit->fstar, model_y->spline->f[k], Spm, Sop);
    double scale = pow(10., -2.0 + 2.0*gsl_rng_uniform(chain->r[ic]));
    model_y->spline->C[0][0][k] += scale*Sn*gsl_ran_gaussian(chain->r[ic],1);
    model_y->spline->C[1][1][k] += scale*Sn*gsl_ran_gaussian(chain->r[ic],1);
    

    
    //compute spline
    if(!flags->prior)
    {
        /* compute spline model */
        //generate_spline_noise_model(model_y); //full interpolation
        update_spline_noise_model(model_y, k, kmin, kmax); //interpolation over stencil
        
        /* get spline model likelihood */
        model_y->logL = noise_log_likelihood(data, model_y->psd);
        //model_y->logL = model_x->logL + noise_delta_log_likelihood(data, model_x, model_y, model_x->spline->f[kmin] , model_x->spline->f[kmax],ic);

        
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
    
    free_spline_model(model_y);
}

void noise_spline_model_rjmcmc(struct Orbit *orbit, struct Data *data, struct SplineModel *model, struct Chain *chain, struct Flags *flags, int ic)
{
    double logH  = 0.0; //(log) Hastings ratio
    double loga  = 1.0; //(log) transition probability
    
    double logPx  = 0.0; //(log) prior density for model x (current state)
    double logPy  = 0.0; //(log) prior density for model y (proposed state)
    
    //shorthand pointers
    struct SplineModel *model_x = model;
    struct SplineModel *model_y = malloc(sizeof(struct SplineModel));
    
    //decide if doing a birth or death move
    int Nspline;
    char move;
    if(gsl_rng_uniform(chain->r[ic])>0.5)
    {
        Nspline = model_x->spline->N + 1; //birth move
        move = 'B';
    }
    else
    {
        Nspline = model_x->spline->N - 1; //death move
        move = 'D';
    }
    alloc_spline_model(model_y, model_x->psd->N, data->Nchannel, Nspline);
    model_y->Nmin = model_x->Nmin;
    model_y->Nmax = model_x->Nmax;

    //alisases for pointers to frequency vectors
    double *fx = model_x->spline->f;
    double *fy = model_y->spline->f;
    
    copy_noise(model_x->psd,model_y->psd);

    //check range
    if(Nspline < model_x->Nmin || Nspline >= model_x->Nmax)
        move = 'R'; //reject move
    
    int kmin,kmax;
    switch(move)
    {
        case 'B':
            /*
             the move is to slot a new spline point between two existing
             points [kmin,kmax]
             */
            
            // first pick the left most point (it can't be the last point)
            kmin = (int)floor(gsl_rng_uniform(chain->r[0])*(double)(model_x->spline->N-1));
            kmax = kmin+1;
            
            //copy current state into trial
            for(int k=0; k<=kmin; k++)
            {
                model_y->spline->f[k] = model_x->spline->f[k];
                model_y->spline->C[0][0][k] = model_x->spline->C[0][0][k];
                model_y->spline->C[1][1][k] = model_x->spline->C[1][1][k];
                model_y->spline->invC[0][0][k] = model_x->spline->invC[0][0][k];
                model_y->spline->invC[1][1][k] = model_x->spline->invC[1][1][k];
                model_y->spline->detC[k] = model_x->spline->detC[k];
            }
            
            //get grid place for new point
            int birth = kmin+1;
            fy[birth] = uniform_frequency_draw(fx[kmin], fx[kmax], chain->r[ic]);
            
            //check frequency prior
            if(check_frequency_spacing(fy, birth, data->T)) logPy = -INFINITY;

            double Spm,Sop;
            get_noise_levels("sangria", model_y->spline->f[birth], &Spm, &Sop);

            double Sn = AEnoise_FF(orbit->L, orbit->fstar, model_y->spline->f[birth], Spm, Sop);//noise_transfer_function(model_y->spline->f[birth]/orbit->fstar);
            double Snmin = -Sn*100.;
            double Snmax =  Sn*100.;
            model_y->spline->C[0][0][birth] = Snmin + (Snmax - Snmin)*gsl_rng_uniform(chain->r[ic]);
            model_y->spline->C[1][1][birth] = Snmin + (Snmax - Snmin)*gsl_rng_uniform(chain->r[ic]);
            
            // now fill in all higher points over k index
            for(int k=kmax; k<model_x->spline->N; k++)
            {
                model_y->spline->f[k+1] = model_x->spline->f[k];
                model_y->spline->C[0][0][k+1] = model_x->spline->C[0][0][k];
                model_y->spline->C[1][1][k+1] = model_x->spline->C[1][1][k];
            }
            break;

        case 'D':
            /*
             the move is to remove any existing point between the end points (0,Nspline)
             */
            
            // first pick the left most point (it can't be the last point)
            kmin = 1;
            kmax = model_x->spline->N - 1;
            int kill = kmin + (int)floor( (double)(kmax-kmin)*gsl_rng_uniform(chain->r[ic]) );

            //copy current state into trial
            for(int k=0; k<kill; k++)
            {
                model_y->spline->f[k] = model_x->spline->f[k];
                model_y->spline->C[0][0][k] = model_x->spline->C[0][0][k];
                model_y->spline->C[1][1][k] = model_x->spline->C[1][1][k];
            }
            
            // now fill in all higher points over k index
            for(int k=kill; k<model_y->spline->N; k++)
            {
                model_y->spline->f[k] = model_x->spline->f[k+1];
                model_y->spline->C[0][0][k] = model_x->spline->C[0][0][k+1];
                model_y->spline->C[1][1][k] = model_x->spline->C[1][1][k+1];
            }
            
            break;
            
        case 'R':
            logPy = -INFINITY;
            break;
            
        default:
            printf("Invalid case %c\n",move);
    }
    
    fflush(stdout);
    
    //compute Hasting's ratio
    if(logPy > -INFINITY && !flags->prior)
    {

        generate_spline_noise_model(model_y);
        
        //compute likelihood
        model_y->logL = noise_log_likelihood(data, model_y->psd);
        
        /*
         H = [p(d|y)/p(d|x)]/T x p(y)/p(x) x q(x|y)/q(y|x)
         */
        logH += (model_y->logL - model_x->logL)/chain->temperature[ic]; //delta logL
    }
    logH += logPy - logPx; //priors
    
    loga = log(gsl_rng_uniform(chain->r[ic]));
    if(logH > loga)
    {
        //reallocate noise structure to size of new spline model
        //realloc_noise(model_x->spline, model_y->spline->N);
        
        free_noise(model_x->spline);
        model_x->spline = malloc(sizeof(struct Noise));
        alloc_noise(model_x->spline,data->Nchannel, model_y->spline->N);
         
        
        copy_spline_model(model_y, model_x);
    }
    
    free_spline_model(model_y);
}

void noise_instrument_model_mcmc(struct Orbit *orbit, struct Data *data, struct InstrumentModel *model, struct ForegroundModel *galaxy, struct Chain *chain, struct Flags *flags, int ic)
{
    double logH  = 0.0; //(log) Hastings ratio
    double loga  = 1.0; //(log) transition probability
    
    double logPx  = 0.0; //(log) prior density for model x (current state)
    double logPy  = 0.0; //(log) prior density for model y (proposed state)
    
    //shorthand pointers
    struct InstrumentModel *model_x = model;
    struct InstrumentModel *model_y = malloc(sizeof(struct InstrumentModel));
    alloc_instrument_model(model_y, data->N, data->Nchannel);
    copy_instrument_model(model_x,model_y);
    
    //structure for full noise covariance matrix
    struct Noise *psd =  malloc(sizeof(struct Noise));
    alloc_noise(psd, data->N, data->Nchannel);
    
    //set priors
    double Sacc = 9.00e-30;
    double Soms = 2.25e-22;
    double Sacc_min = 9.00e-30/10;
    double Sacc_max = 9.00e-30*10;
    double Soms_min = 2.25e-22/10;
    double Soms_max = 2.25e-22*10;

    //get jump sizes
    double scale;
    if(gsl_rng_uniform(chain->r[ic])>0.75)
        scale = 0.1;
    else if(gsl_rng_uniform(chain->r[ic])>0.5)
        scale = 0.01;
    else if(gsl_rng_uniform(chain->r[ic])>0.25)
        scale = 0.001;
    else
        scale = 0.001;
    
    /* get proposed noise parameters */
    
    //update one link at a time
    int i = (int)(gsl_rng_uniform(chain->r[ic])* (double)model_x->Nlink);
    model_y->sacc[i] = model_x->sacc[i] + scale * Sacc * gsl_ran_gaussian(chain->r[ic],1);
    model_y->soms[i] = model_x->soms[i] + scale * Soms * gsl_ran_gaussian(chain->r[ic],1);
    
    //check priors
    if(model_y->sacc[i] < Sacc_min || model_y->sacc[i] > Sacc_max) logPy = -INFINITY;
    if(model_y->soms[i] < Soms_min || model_y->soms[i] > Soms_max) logPy = -INFINITY;

    //OMS noise is degenerate on a link
    model_y->soms[1] = model_y->soms[0]; //Soms12 and Soms21
    model_y->soms[3] = model_y->soms[2]; //Soms23 and Soms32
    model_y->soms[5] = model_y->soms[4]; //Soms13 and Soms31

    //get noise covariance matrix for initial parameters
    if(logPy > -INFINITY && !flags->prior)
    {
        generate_instrument_noise_model(data,orbit,model_y);
        copy_noise(model_y->psd,psd);

        //add foreground noise contribution
        if(flags->confNoise) 
            generate_full_covariance_matrix(psd,galaxy->psd, data->Nchannel);
        
        invert_noise_covariance_matrix(psd);

        model_y->logL = noise_log_likelihood(data, psd);

        logH += (model_y->logL - model_x->logL)/chain->temperature[ic]; //delta logL
    }
    logH += logPy - logPx; //priors

    loga = log(gsl_rng_uniform(chain->r[ic]));
    if(logH > loga)
    {
        copy_instrument_model(model_y, model_x);
        if(flags->confNoise) galaxy->logL = model_x->logL;
    }
    
    free_noise(psd);
    free_instrument_model(model_y);
}

void noise_foreground_model_mcmc(struct Orbit *orbit, struct Data *data, struct InstrumentModel *noise, struct ForegroundModel *model, struct Chain *chain, struct Flags *flags, int ic)
{
    double logH  = 0.0; //(log) Hastings ratio
    double loga  = 1.0; //(log) transition probability
    
    double logPx  = 0.0; //(log) prior density for model x (current state)
    double logPy  = 0.0; //(log) prior density for model y (proposed state)
    
    //shorthand pointers
    struct ForegroundModel *model_x = model;
    struct ForegroundModel *model_y = malloc(sizeof(struct ForegroundModel));
    alloc_foreground_model(model_y, data->N, data->Nchannel);
    copy_foreground_model(model_x,model_y);
    
    //structure for full noise covariance matrix
    struct Noise *psd =  malloc(sizeof(struct Noise));
    alloc_noise(psd, data->N, data->Nchannel);

    /* set priors */
    double **prior = malloc(sizeof(double *)*model_x->Nparams);
    for(int n=0; n<model_x->Nparams; n++)
        prior[n] = malloc(sizeof(double)*2);
    
    //log(A)
    prior[0][0] = -86.0;
    prior[0][1] = -82.0;
    
    //f1
    prior[1][0] = log(0.0001);
    prior[1][1] = log(0.01);
    
    //alpha
    prior[2][0] = 0.0;
    prior[2][1] = 3.0;
    
    //fk
    prior[3][0] = log(0.0001);
    prior[3][1] = log(0.01);
    
    //f2
    prior[4][0] = log(0.0001);
    prior[4][1] = log(0.01);
        
    /* get proposed noise parameters */

    //get jump sizes
    double scale;
    if(gsl_rng_uniform(chain->r[ic])>0.75)
        scale = 1;
    else if(gsl_rng_uniform(chain->r[ic])>0.5)
        scale = 0.1;
    else if(gsl_rng_uniform(chain->r[ic])>0.25)
        scale = 0.01;
    else
        scale = 0.01;

    //pick which parameter to update
    int i = (int)(gsl_rng_uniform(chain->r[ic])* (double)model_x->Nparams);
    model_y->sgal[i] = model_x->sgal[i] + scale * 0.5*(prior[i][1]-prior[i][0]) * gsl_ran_gaussian(chain->r[ic],1);

    //check priors
    for(int n=0; n<model_y->Nparams; n++)
        if(model_y->sgal[n] < prior[n][0] || model_y->sgal[n] > prior[n][1]) 
            logPy = -INFINITY;

    //get noise covariance matrix for initial parameters
    if(logPy > -INFINITY && !flags->prior)
    {
        generate_galactic_foreground_model(data,orbit,model_y);
        copy_noise(model_y->psd,psd);

        //add instrument noise contribution
        generate_full_covariance_matrix(psd, noise->psd, data->Nchannel);
        
        invert_noise_covariance_matrix(psd);

        model_y->logL = noise_log_likelihood(data, psd);

        logH += (model_y->logL - model_x->logL)/chain->temperature[ic]; //delta logL
    }
    logH += logPy - logPx; //priors

    loga = log(gsl_rng_uniform(chain->r[ic]));
    if(logH > loga)
    {
        copy_foreground_model(model_y, model_x);
        noise->logL = model_x->logL;
    }
    
    free_noise(psd);
    free_foreground_model(model_y);
    for(int n=0; n<model_x->Nparams; n++) free(prior[n]);
    free(prior);

}


void initialize_spline_model(struct Orbit *orbit, struct Data *data, struct SplineModel *model, int Nspline)
{
    
    // Initialize data models
    alloc_spline_model(model, data->N, data->Nchannel, Nspline);
    
    //set max and min spline points
    model->Nmin = MIN_SPLINE_STENCIL;
    model->Nmax = Nspline;
    
    //set up psd frequency grid
    for(int n=0; n<model->psd->N; n++)
    {
        double f = data->fmin + (double)n/data->T;
        double Spm, Sop;
        get_noise_levels("sangria", f, &Spm, &Sop);
        model->psd->f[n] = f;
        model->psd->transfer[n] = AEnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);//noise_transfer_function(f/orbit->fstar);
    }
    
    //divide into Nspline control points
    double logdf = (log(data->fmax) - log(data->fmin))/(Nspline-1);
    for(int i=0; i<Nspline; i++)
    {
        double f = exp(log(data->fmin) + (double)i*logdf);
        model->spline->f[i] = f;
        
        /* initialize model to theoretical level without transfer function applied */
        model->spline->C[0][0][i] = 0.0;//AEnoise_FF(orbit->L, orbit->fstar, f)/noise_transfer_function(f/orbit->fstar);
        model->spline->C[1][1][i] = 0.0;//AEnoise_FF(orbit->L, orbit->fstar, f)/noise_transfer_function(f/orbit->fstar);
    }
    //shift first spline control point by half a bin to avoid rounding problems
    model->spline->f[0] -= 0.5/data->T;
    
    generate_spline_noise_model(model);
    model->logL = noise_log_likelihood(data, model->psd);
}

void initialize_instrument_model(struct Orbit *orbit, struct Data *data, struct InstrumentModel *model)
{
    // initialize data models
    alloc_instrument_model(model, data->N, data->Nchannel);
    
    // set up psd frequency grid
    for(int n=0; n<model->psd->N; n++)
        model->psd->f[n] = data->fmin + (double)n/data->T;

    // initialize noise levels
    for(int i=0; i<model->Nlink; i++)
    {
        model->soms[i] = 2.25e-22;
        model->sacc[i] = 9.00e-30;
    }
    
    // get noise covariance matrix for initial parameters
    generate_instrument_noise_model(data,orbit,model);
}

void initialize_foreground_model(struct Orbit *orbit, struct Data *data, struct ForegroundModel *model)
{
    // initialize data models
    alloc_foreground_model(model, data->N, data->Nchannel);
    
    // set up psd frequency grid
    for(int n=0; n<model->psd->N; n++)
        model->psd->f[n] = data->fmin + (double)n/data->T;

    // initialize foreground parameters levels
    model->Tobs  =  data->T;
    model->Amp   =  1e-36;
    model->f1    =  0.001;
    model->alpha =  1.58;
    model->fk    =  0.001;
    model->f2    =  0.001;
    map_foreground_params_to_array(model);
    
    // get noise covariance matrix for initial parameters
    generate_galactic_foreground_model(data,orbit,model);
}

void print_noise_model(struct Noise *noise, char filename[])
{
    FILE *fptr = fopen(filename,"w");
    for(int i=0; i<noise->N; i++)
    {
        fprintf(fptr,"%lg ",noise->f[i]);
        for(int j=0; j<noise->Nchannel; j++)
            fprintf(fptr,"%lg ",noise->C[j][j][i]);
        fprintf(fptr,"%lg ",noise->C[0][1][i]);
        fprintf(fptr,"%lg ",noise->C[0][2][i]);
        fprintf(fptr,"%lg ",noise->C[1][2][i]);
        fprintf(fptr,"\n");
    }
    fclose(fptr);
    
}

void print_whitened_data(struct Data *data, struct Noise *noise, char filename[])
{
    FILE *fptr = fopen(filename,"w");
    for(int i=0; i<noise->N; i++)
    {
        fprintf(fptr,"%lg ",noise->f[i]);
        fprintf(fptr,"%lg %lg ",data->tdi[0]->X[2*i]/sqrt(noise->C[0][0][i]),data->tdi[0]->X[2*i+1]/sqrt(noise->C[0][0][i]));
        fprintf(fptr,"%lg %lg ",data->tdi[0]->Y[2*i]/sqrt(noise->C[1][1][i]),data->tdi[0]->Y[2*i+1]/sqrt(noise->C[1][1][i]));
        fprintf(fptr,"%lg %lg ",data->tdi[0]->Z[2*i]/sqrt(noise->C[2][2][i]),data->tdi[0]->Z[2*i+1]/sqrt(noise->C[2][2][i]));
        fprintf(fptr,"\n");
    }
        
    fclose(fptr);
}
