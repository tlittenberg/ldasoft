/*
 *  Copyright (C) 2023 Tyson B. Littenberg (MSFC-ST12)
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

#include <glass_utils.h>

#include "glass_noise.h"


void map_array_to_noise_params(struct InstrumentModel *model)
{
    model->sacc12 = model->sacc[0];
    model->sacc21 = model->sacc[1];

    model->sacc23 = model->sacc[2];
    model->sacc32 = model->sacc[3];

    model->sacc13 = model->sacc[4];
    model->sacc31 = model->sacc[5];

    model->soms12 = model->soms[0];
    model->soms21 = model->soms[1];
        
    model->soms23 = model->soms[2];
    model->soms32 = model->soms[3];
    
    model->soms13 = model->soms[4];
    model->soms31 = model->soms[5];

}

void map_noise_params_to_array(struct InstrumentModel *model)
{
    model->sacc[0] = model->sacc12;
    model->sacc[1] = model->sacc21;
        
    model->sacc[2] = model->sacc23;
    model->sacc[3] = model->sacc32;
    
    model->sacc[4] = model->sacc13;
    model->sacc[5] = model->sacc31;

    model->soms[0] = model->soms12;
    model->soms[1] = model->soms21;
    
    model->soms[2] = model->soms23;
    model->soms[3] = model->soms32;

    model->soms[4] = model->soms13;
    model->soms[5] = model->soms31;
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
    model->Nchannel = Nchannel;
    
    model->psd = malloc(sizeof(struct Noise));
    model->spline=malloc(sizeof(struct Noise));
    
    alloc_noise(model->psd, Ndata, Nchannel);
    alloc_noise(model->spline, Nspline, Nchannel);
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
    copy->Nchannel = origin->Nchannel;

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

void update_spline_noise_model(struct SplineModel *model, int new_knot, int min_knot, int max_knot)
{
    struct Noise *psd = model->psd;
    struct Noise *spline = model->spline;
    
    /* find location in data vector of knot */
    double T = 1./(psd->f[1] - psd->f[0]);

    gsl_spline **cspline = malloc(model->Nchannel*sizeof(gsl_spline *));
    gsl_interp_accel **acc = malloc(model->Nchannel*sizeof(gsl_interp_accel *));
    
    /* have to recompute the spline everywhere (derivatives on boundary) */
    for(int n=0; n<model->Nchannel; n++)
    {
        cspline[n] = gsl_spline_alloc(gsl_interp_akima, spline->N);
        acc[n] = gsl_interp_accel_alloc();
        gsl_spline_init(cspline[n],spline->f,spline->C[n][n],spline->N);
    }
    
    int imin = (int)((spline->f[min_knot]-psd->f[0])*T);
    int imax = (int)((spline->f[max_knot]-psd->f[0])*T);
    
    
    for(int i=imin; i<imax; i++)
    {
        for(int n=0; n<model->Nchannel; n++)
        {
            psd->C[n][n][i]=gsl_spline_eval(cspline[n],psd->f[i],acc[n]);
            
            /*
             apply transfer function
             -this catches the sharp features in the spectrum from f/fstar
             -without needing to interpolate
             */
            psd->C[n][n][i]+=psd->transfer[i];
        }
    }
    invert_noise_covariance_matrix(psd);

    for(int n=0; n<model->Nchannel; n++)
    {
        gsl_spline_free(cspline[n]);
        gsl_interp_accel_free(acc[n]);
    }
    free(cspline);
    free(acc);

}


void generate_spline_noise_model(struct SplineModel *model)
{
    struct Noise *psd = model->psd;
    struct Noise *spline = model->spline;
    
    for(int i=0; i<model->Nchannel; i++)
        CubicSplineGSL(spline->N, spline->f, spline->C[i][i], psd->N, psd->f, psd->C[i][i]);

    //set a floor on Sn so likelihood doesn't go crazy
    for(int n=0; n<psd->N; n++)
    {
        /*
         apply transfer function
         -this catches the sharp features in the spectrum from f/fstar
         -without needing to interpolate
         */
        for(int i=0; i<model->Nchannel; i++)
            psd->C[i][i][n]+=psd->transfer[n];
        
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
    double Sgal;
    
    map_array_to_foreground_params(model);
    
    for(int n=0; n<data->N; n++)
    {
        f = model->psd->f[n];
        
        Sgal = model->Amp*pow(f,5./3.) * exp(-pow(f/model->f1,model->alpha)) * 0.5*( 1. + tanh( (model->fk - f)/model->f2 ) );

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
    
    struct TDI *tdi = data->tdi;
    
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
    
    struct TDI *tdi = data->tdi;
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

void initialize_spline_model(struct Orbit *orbit, struct Data *data, struct SplineModel *model, int Nspline)
{
    
    // Initialize data models
    alloc_spline_model(model, data->N, data->Nchannel, Nspline);
    
    //set max and min spline points
    model->Nmin = MIN_SPLINE_STENCIL;
    model->Nmax = Nspline;
    model->Nchannel = data->Nchannel;
    
    //set up psd frequency grid
    for(int n=0; n<model->psd->N; n++)
    {
        double f = data->fmin + (double)n/data->T;
        double Spm, Sop;
        get_noise_levels("sangria", f, &Spm, &Sop);
        model->psd->f[n] = f;
        if(model->Nchannel==2) model->psd->transfer[n] = AEnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);//noise_transfer_function(f/orbit->fstar);
        else model->psd->transfer[n] = XYZnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);
    }
    
    //divide into Nspline control points
    double logdf = (log(data->fmax) - log(data->fmin))/(Nspline-1);
    for(int i=0; i<Nspline; i++)
    {
        double f = exp(log(data->fmin) + (double)i*logdf);
        model->spline->f[i] = f;
        
        for(int n=0; n<model->Nchannel; n++)
        {
            /* initialize model to theoretical level without transfer function applied */
            for(int m=0; m<model->Nchannel; m++)
                model->spline->C[n][m][i] = 0.0;
        }
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
    model->Amp   =  3.86677e-37;
    model->f1    =  0.00344439;
    model->alpha =  1.629667;
    model->fk    =  0.0102644;
    model->f2    =  4.810781e-4;
    map_foreground_params_to_array(model);
    
    // get noise covariance matrix for initial parameters
    generate_galactic_foreground_model(data,orbit,model);
}
