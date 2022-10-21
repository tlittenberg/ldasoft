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
    
    //Priors
    copy->Nmin = origin->Nmin;
    copy->Nmax = origin->Nmax;
}

void print_spline_state(struct SplineModel *model, FILE *fptr, int step)
{
    fprintf(fptr,"%i %.12g %i\n",step,model->logL,model->spline->N);
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
    gsl_spline_init(cspline_A,spline->f,spline->SnA,spline->N);
    gsl_spline_init(cspline_E,spline->f,spline->SnE,spline->N);

    
    int imin = (int)((spline->f[min_knot]-psd->f[0])*T);
    int imax = (int)((spline->f[max_knot]-psd->f[0])*T);
    
    
    for(int n=imin; n<imax; n++)
    {
        psd->SnA[n]=gsl_spline_eval(cspline_A,psd->f[n],acc_A);
        psd->SnE[n]=gsl_spline_eval(cspline_E,psd->f[n],acc_E);
        
        /*
         apply transfer function
         -this catches the sharp features in the spectrum from f/fstar
         -without needing to interpolate
         */
        
        psd->SnA[n]+=psd->transfer[n];
        psd->SnE[n]+=psd->transfer[n];
    }
    
    gsl_spline_free (cspline_A);
    gsl_spline_free (cspline_E);
    gsl_interp_accel_free (acc_A);
    gsl_interp_accel_free (acc_E);

}


void generate_spline_noise_model(struct SplineModel *model)
{
    struct Noise *psd = model->psd;
    struct Noise *spline = model->spline;
    
    CubicSplineGSL(spline->N, spline->f, spline->SnA, psd->N, psd->f, psd->SnA);
    CubicSplineGSL(spline->N, spline->f, spline->SnE, psd->N, psd->f, psd->SnE);

    //set a floor on Sn so likelihood doesn't go crazy
    for(int n=0; n<psd->N; n++)
    {
        /*
         apply transfer function
         -this catches the sharp features in the spectrum from f/fstar
         -without needing to interpolate
         */
        psd->SnA[n]+=psd->transfer[n];
        psd->SnE[n]+=psd->transfer[n];
    }
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
    dlogL -= -0.5*fourier_nwip(tdi->A+2*imin, tdi->A+2*imin, psd_x->SnA+imin, N);
    dlogL -= -0.5*fourier_nwip(tdi->E+2*imin, tdi->E+2*imin, psd_x->SnE+imin, N);
    for(int n=imin; n<imin+N; n++)
    {
        dlogL += log(psd_x->SnA[n]);
        dlogL += log(psd_x->SnE[n]);
    }

    /* add contribution for proposed state y */
    dlogL += -0.5*fourier_nwip(tdi->A+2*imin, tdi->A+2*imin, psd_y->SnA+imin, N);
    dlogL += -0.5*fourier_nwip(tdi->E+2*imin, tdi->E+2*imin, psd_y->SnE+imin, N);
    for(int n=imin; n<imin+N; n++)
    {
        dlogL -= log(psd_y->SnA[n]);
        dlogL -= log(psd_y->SnE[n]);
    }
    
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
    alloc_spline_model(model_y, model_x->psd->N, model_x->spline->N);

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
    double Sn = AEnoise_FF(orbit->L, orbit->fstar, model_y->spline->f[k]);
    double scale = pow(10., -2.0 + 2.0*gsl_rng_uniform(chain->r[ic]));
    model_y->spline->SnA[k] += scale*Sn*gsl_ran_gaussian(chain->r[ic],1);
    model_y->spline->SnE[k] += scale*Sn*gsl_ran_gaussian(chain->r[ic],1);
    

    
    //compute spline
    if(!flags->prior)
    {
        /* compute spline model */
        //generate_spline_noise_model(model_y); //full interpolation
        update_spline_noise_model(model_y, k, kmin, kmax); //interpolation over stencil
        
        /* get spline model likelihood */
        model_y->logL = noise_log_likelihood(data, model_y);
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
    alloc_spline_model(model_y, model_x->psd->N, Nspline);
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
                model_y->spline->SnA[k] = model_x->spline->SnA[k];
                model_y->spline->SnE[k] = model_x->spline->SnE[k];
            }
            
            //get grid place for new point
            int birth = kmin+1;
            fy[birth] = uniform_frequency_draw(fx[kmin], fx[kmax], chain->r[ic]);
            
            //check frequency prior
            if(check_frequency_spacing(fy, birth, data->T)) logPy = -INFINITY;

            
            double Sn = AEnoise_FF(orbit->L, orbit->fstar, model_y->spline->f[birth]);//noise_transfer_function(model_y->spline->f[birth]/orbit->fstar);
            double Snmin = -Sn*100.;
            double Snmax =  Sn*100.;
            model_y->spline->SnA[birth] = Snmin + (Snmax - Snmin)*gsl_rng_uniform(chain->r[ic]);
            model_y->spline->SnE[birth] = Snmin + (Snmax - Snmin)*gsl_rng_uniform(chain->r[ic]);
            
            // now fill in all higher points over k index
            for(int k=kmax; k<model_x->spline->N; k++)
            {
                model_y->spline->f[k+1] = model_x->spline->f[k];
                model_y->spline->SnA[k+1] = model_x->spline->SnA[k];
                model_y->spline->SnE[k+1] = model_x->spline->SnE[k];
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
                model_y->spline->SnA[k] = model_x->spline->SnA[k];
                model_y->spline->SnE[k] = model_x->spline->SnE[k];
            }
            
            // now fill in all higher points over k index
            for(int k=kill; k<model_y->spline->N; k++)
            {
                model_y->spline->f[k] = model_x->spline->f[k+1];
                model_y->spline->SnA[k] = model_x->spline->SnA[k+1];
                model_y->spline->SnE[k] = model_x->spline->SnE[k+1];
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
        //reallocate noise structure to size of new spline model
        //realloc_noise(model_x->spline, model_y->spline->N);
        
        free_noise(model_x->spline);
        model_x->spline = malloc(sizeof(struct Noise));
        alloc_noise(model_x->spline,model_y->spline->N);
         
        
        copy_spline_model(model_y, model_x);
    }
    
    free_spline_model(model_y);
}

void initialize_spline_model(struct Orbit *orbit, struct Data *data, struct SplineModel *model, int Nspline)
{
    
    // Initialize data models
    alloc_spline_model(model, data->N, Nspline);
    
    //set max and min spline points
    model->Nmin = MIN_SPLINE_STENCIL;
    model->Nmax = Nspline;
    
    //set up psd frequency grid
    for(int n=0; n<model->psd->N; n++)
    {
        double f = data->fmin + (double)n/data->T;
        model->psd->f[n] = f;
        model->psd->transfer[n] = AEnoise_FF(orbit->L, orbit->fstar, f);//noise_transfer_function(f/orbit->fstar);
    }
    
    //divide into Nspline control points
    double logdf = (log(data->fmax) - log(data->fmin))/(Nspline-1);
    for(int i=0; i<Nspline; i++)
    {
        double f = exp(log(data->fmin) + (double)i*logdf);
        model->spline->f[i] = f;
        
        /* initialize model to theoretical level without transfer function applied */
        model->spline->SnA[i] = 0.0;//AEnoise_FF(orbit->L, orbit->fstar, f)/noise_transfer_function(f/orbit->fstar);
        model->spline->SnE[i] = 0.0;//AEnoise_FF(orbit->L, orbit->fstar, f)/noise_transfer_function(f/orbit->fstar);
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

