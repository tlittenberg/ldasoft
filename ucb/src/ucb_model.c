/*
 *  Copyright (C) 2019 Tyson B. Littenberg (MSFC-ST12), Neil J. Cornish
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


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lisa.h>
#include <data.h>

#include "ucb_model.h"
#include "ucb_catalog.h"
#include "ucb_io.h"
#include "ucb_prior.h"
#include "ucb_waveform.h"
#include "ucb_fstatistic.h"

#define FIXME 0
#define find_max(x,y) (((x) >= (y)) ? (x) : (y))
#define find_min(x,y) (((x) <= (y)) ? (x) : (y))

void map_array_to_params(struct Source *source, double *params, double T)
{
    source->f0       = params[0]/T;
    source->costheta = params[1];
    source->phi      = params[2];
    source->amp      = exp(params[3]);
    source->cosi     = params[4];
    source->psi      = params[5];
    source->phi0     = params[6];
    if(NP>7)
        source->dfdt   = params[7]/(T*T);
    if(NP>8)
        source->d2fdt2 = params[8]/(T*T*T);
}

void map_params_to_array(struct Source *source, double *params, double T)
{
    params[0] = source->f0*T;
    params[1] = source->costheta;
    params[2] = source->phi;
    params[3] = log(source->amp);
    params[4] = source->cosi;
    params[5] = source->psi;
    params[6] = source->phi0;
    if(NP>7)
        params[7] = source->dfdt*T*T;
    if(NP>8)
        params[8] = source->d2fdt2*T*T*T;
}

void alloc_model(struct Model *model, int Nmax, int NFFT, int Nchannel)
{
    int n;
    
    model->Nlive  = 1;
    model->Nmax   = Nmax;
    model->t0 = 0.0;

    model->source = malloc(model->Nmax*sizeof(struct Source *));
    
    model->calibration = malloc(sizeof(struct Calibration) );
    model->noise       = malloc(sizeof(struct Noise)       );
    model->tdi         = malloc(sizeof(struct TDI)         );
    model->residual    = malloc(sizeof(struct TDI)         );
    
    alloc_noise(model->noise,NFFT, Nchannel);
    alloc_tdi(model->tdi, NFFT, Nchannel);
    alloc_tdi(model->residual, NFFT, Nchannel);
    alloc_calibration(model->calibration);
    
    for(n=0; n<model->Nmax; n++)
    {
        model->source[n] = malloc(sizeof(struct Source));
        alloc_source(model->source[n],NFFT,Nchannel);
    }
    
    model->logPriorVolume = calloc(NP,sizeof(double));
    model->prior = malloc(NP*sizeof(double *));
    for(n=0; n<NP; n++) model->prior[n] = calloc(2,sizeof(double));
}

void copy_model(struct Model *origin, struct Model *copy)
{
    //Source parameters
    copy->Nmax           = origin->Nmax;
    copy->Nlive          = origin->Nlive;
    for(int n=0; n<origin->Nlive; n++)
        copy_source(origin->source[n],copy->source[n]);
    
    //Noise parameters
    copy_noise(origin->noise,copy->noise);
    
    //TDI
    copy_tdi(origin->tdi,copy->tdi);
    
    //Calibration parameters
    copy_calibration(origin->calibration,copy->calibration);
    
    //Residual
    copy_tdi(origin->residual,copy->residual);
    
    //Start time for segment for model
    copy->t0 = origin->t0;
    
    
    //Source parameter priors
    for(int n=0; n<NP; n++)
    {
        for(int j=0; j<2; j++) copy->prior[n][j] = origin->prior[n][j];
        copy->logPriorVolume[n] = origin->logPriorVolume[n];
    }
    //Model likelihood
    copy->logL           = origin->logL;
    copy->logLnorm       = origin->logLnorm;
}

void copy_model_lite(struct Model *origin, struct Model *copy)
{
    copy->Nlive = origin->Nlive;

    //Source parameters
    for(int n=0; n<origin->Nlive; n++)
        copy_source(origin->source[n],copy->source[n]);
    
    //Source waveforms
    copy_tdi(origin->tdi,copy->tdi);
    copy_tdi(origin->residual,copy->residual);
    
    //Model likelihood
    copy->logL = origin->logL;
}

int compare_model(struct Model *a, struct Model *b)
{
    
    //Source parameters
    if(a->Nmax  != b->Nmax)  return 1;  //maximum number of signals in model
    if(a->Nlive != b->Nlive) return 1;  //current number of signals in model
    
    //struct Source **source;
    for(int i=0; i<a->Nlive; i++)
    {
        struct Source *sa = a->source[i];
        struct Source *sb = b->source[i];
        
        //Intrinsic
        if(sa->m1 != sb->m1) return 1;
        if(sa->m2 != sb->m2) return 1;
        if(sa->f0 != sb->f0) return 1;
        
        //Extrinisic
        if(sa->psi  != sb->psi)  return 1;
        if(sa->cosi != sb->cosi) return 1;
        if(sa->phi0 != sb->phi0) return 1;
        
        if(sa->D        != sb->D)        return 1;
        if(sa->phi      != sb->phi)      return 1;
        if(sa->costheta != sb->costheta) return 1;
        
        //Derived
        if(sa->amp    != sb->amp)    return 1;
        if(sa->dfdt   != sb->dfdt)   return 1;
        if(sa->d2fdt2 != sb->d2fdt2) return 1;
        if(sa->Mc     != sb->Mc)     return 1;
        
        //Book-keeping
        if(sa->BW   != sb->BW)   return 1;
        if(sa->qmin != sb->qmin) return 1;
        if(sa->qmax != sb->qmax) return 1;
        if(sa->imin != sb->imin) return 1;
        if(sa->imax != sb->imax) return 1;
        
        //Response
        //TDI
        struct TDI *tsa = sa->tdi;
        struct TDI *tsb = sb->tdi;
        
        if(tsa->N != tsb->N) return 1;
        if(tsa->Nchannel != tsb->Nchannel) return 1;
        
        for(int i=0; i<2*tsa->N; i++)
        {
            
            //Michelson
            if(tsa->X[i] != tsb->X[i]) return 1;
            if(tsa->Y[i] != tsb->Y[i]) return 1;
            if(tsa->Z[i] != tsb->Z[i]) return 1;
            
            //Noise-orthogonal
            if(tsa->A[i] != tsb->A[i]) return 1;
            if(tsa->E[i] != tsb->E[i]) return 1;
            if(tsa->T[i] != tsb->T[i]) return 1;
        }
        
        //Package parameters for waveform generator
        for(int j=0; j<NP; j++) if(sa->params[j] != sb->params[j]) return 1;
        
    }
    
    //Noise parameters
    struct Noise *na = a->noise;
    struct Noise *nb = b->noise;
    
    if(na->N != nb->N) return 1;
    if(na->Nchannel != nb->Nchannel) return 1;

    for(int n=0; n<na->Nchannel; n++)
    {
        if(na->eta[n] != nb->eta[n]) return 1;
        for(int i=0; i<na->N; i++) if(na->detC[i] != nb->detC[i]) return 1;
        for(int m=n; m<na->Nchannel; m++)
        {
            for(int i=0; i<na->N; i++)
            {
                if(na->C[n][m][i] != nb->C[n][m][i]) return 1;
                if(na->invC[n][m][i] != nb->invC[n][m][i]) return 1;
            }
        }
    }
    
    //TDI
    struct TDI *ta = a->tdi;
    struct TDI *tb = b->tdi;
    
    if(ta->N != tb->N) return 1;
    if(ta->Nchannel != tb->Nchannel) return 1;
    
    for(int i=0; i<2*ta->N; i++)
    {
        
        //Michelson
        if(ta->X[i] != tb->X[i]) return 1;
        if(ta->Y[i] != tb->Y[i]) return 1;
        if(ta->Z[i] != tb->Z[i]) return 1;
        
        //Noise-orthogonal
        if(ta->A[i] != tb->A[i]) return 1;
        if(ta->E[i] != tb->E[i]) return 1;
        if(ta->T[i] != tb->T[i]) return 1;
    }
    
    //Start time for segment for model
    if(a->t0 != b->t0)     return 1;
    
    //Source parameter priors
    //double **prior;
    if(a->logPriorVolume != b->logPriorVolume) return 1;
    
    //Model likelihood
    if(a->logL     != b->logL)     return 1;
    if(a->logLnorm != b->logLnorm) return 1;
    
    return 0;
}



void free_model(struct Model *model)
{
    int n;
    for(n=0; n<model->Nmax; n++)
    {
        free_source(model->source[n]);
    }
    free(model->source);

    for(n=0; n<NP; n++) free(model->prior[n]);
    free(model->prior);
    free(model->logPriorVolume);

    free_tdi(model->tdi);
    free_tdi(model->residual);
    free_noise(model->noise);
    free_calibration(model->calibration);

    free(model->noise);
    free(model->tdi);
    free(model->residual);
    free(model->calibration);
    
    free(model);
}

void alloc_source(struct Source *source, int NFFT, int Nchannel)
{
    //Intrinsic
    source->m1=1.;
    source->m2=1.;
    source->f0=0.;
    
    //Extrinisic
    source->psi=0.0;
    source->cosi=0.0;
    source->phi0=0.0;
    
    source->D=1.0;
    source->phi=0.0;
    source->costheta=0.0;
    
    //Derived
    source->amp=1.;
    source->Mc=1.;
    source->dfdt=0.;
    source->d2fdt2=0.;
    
    //Book-keeping
    source->BW   = NFFT;
    source->qmin = 0;
    source->qmax = NFFT;
    source->imin = 0;
    source->imax = NFFT;
    
    
    //Package parameters for waveform generator
    source->params=calloc(NP,sizeof(double));
    
    //Response
    source->tdi = malloc(sizeof(struct TDI));
    alloc_tdi(source->tdi,NFFT, Nchannel);
    
    //FIsher
    source->fisher_matrix = malloc(NP*sizeof(double *));
    source->fisher_evectr = malloc(NP*sizeof(double *));
    source->fisher_evalue = calloc(NP,sizeof(double));
    for(int i=0; i<NP; i++)
    {
        source->fisher_matrix[i] = calloc(NP,sizeof(double));
        source->fisher_evectr[i] = calloc(NP,sizeof(double));
    }
};

void copy_source(struct Source *origin, struct Source *copy)
{
    //Intrinsic
    copy->m1 = origin->m1;
    copy->m2 = origin->m2;
    copy->f0 = origin->f0;
    
    //Extrinisic
    copy->psi  = origin->psi;
    copy->cosi = origin->cosi;
    copy->phi0 = origin->phi0;
    
    copy->D        = origin->D;
    copy->phi      = origin->phi;
    copy->costheta = origin->costheta;
    
    //Derived
    copy->amp    = origin->amp;
    copy->Mc     = origin->Mc;
    copy->dfdt   = origin->dfdt;
    copy->d2fdt2 = origin->d2fdt2;
    
    //Book-keeping
    copy->BW   = origin->BW;
    copy->qmin = origin->qmin;
    copy->qmax = origin->qmax;
    copy->imin = origin->imin;
    copy->imax = origin->imax;
    
    //Response
    copy_tdi(origin->tdi,copy->tdi);

    
    //Fisher
    memcpy(copy->fisher_evalue, origin->fisher_evalue, NP*sizeof(double));
    memcpy(copy->params, origin->params, NP*sizeof(double));
    copy->fisher_update_flag = origin->fisher_update_flag;
    
    for(int i=0; i<NP; i++)
    {
        memcpy(copy->fisher_matrix[i], origin->fisher_matrix[i], NP*sizeof(double));
        memcpy(copy->fisher_evectr[i], origin->fisher_evectr[i], NP*sizeof(double));
    }
    
}

void free_source(struct Source *source)
{
    for(int i=0; i<NP; i++)
    {
        free(source->fisher_matrix[i]);
        free(source->fisher_evectr[i]);
    }
    free(source->fisher_matrix);
    free(source->fisher_evectr);
    free(source->fisher_evalue);
    free(source->params);
    
    free_tdi(source->tdi);
    
    free(source);
}

void generate_signal_model(struct Orbit *orbit, struct Data *data, struct Model *model, int source_id)
{
    int i,j,n;
    int N2=data->N*2;
    
    for(n=0; n<N2; n++)
    {
        model->tdi->X[n]=0.0;
        model->tdi->Y[n]=0.0;
        model->tdi->Z[n]=0.0;
        model->tdi->A[n]=0.0;
        model->tdi->E[n]=0.0;
    }
    
    //Loop over signals in model
    for(n=0; n<model->Nlive; n++)
    {
        struct Source *source = model->source[n];
        
        if(source_id==-1 || source_id==n)
        {
            for(i=0; i<N2; i++)
            {
                //source->tdi->X[i]=0.0;
                source->tdi->A[i]=0.0;
                source->tdi->E[i]=0.0;
            }
            
            map_array_to_params(source, source->params, data->T);
            
            //Book-keeping of injection time-frequency volume
            galactic_binary_alignment(orbit, data, source);
        }
        
        //Simulate gravitational wave signal
        /* the source_id = -1 condition is redundent if the model->tdi structure is up to date...*/
        if(source_id==-1 || source_id==n) galactic_binary(orbit, data->format, data->T, model->t0, source->params, NP, source->tdi->X, source->tdi->Y, source->tdi->Z, source->tdi->A, source->tdi->E, source->BW, source->tdi->Nchannel);
        
        //Add waveform to model TDI channels
        for(i=0; i<source->BW; i++)
        {
            j = i+source->imin;
            
            if(j>-1 && j<data->N)
            {
                int i_re = 2*i;
                int i_im = i_re+1;
                int j_re = 2*j;
                int j_im = j_re+1;
                
                switch(data->Nchannel)
                {
                    case 1:
                         model->tdi->X[j_re] += source->tdi->X[i_re];
                         model->tdi->X[j_im] += source->tdi->X[i_im];
                        break;
                        
                    case 2:
                        model->tdi->A[j_re] += source->tdi->A[i_re];
                        model->tdi->A[j_im] += source->tdi->A[i_im];
                        
                        model->tdi->E[j_re] += source->tdi->E[i_re];
                        model->tdi->E[j_im] += source->tdi->E[i_im];
                        break;
                    case 3:
                        model->tdi->X[j_re] += source->tdi->X[i_re];
                        model->tdi->X[j_im] += source->tdi->X[i_im];
                        
                        model->tdi->Y[j_re] += source->tdi->Y[i_re];
                        model->tdi->Y[j_im] += source->tdi->Y[i_im];
                        
                        model->tdi->Z[j_re] += source->tdi->Z[i_re];
                        model->tdi->Z[j_im] += source->tdi->Z[i_im];
                        break;
                }
            }//check that source_id is in range
        }//loop over waveform bins
}//loop over sources
}

void update_signal_model(struct Orbit *orbit, struct Data *data, struct Model *model_x, struct Model *model_y, int source_id)
{
    int i,j;
    int N2=data->N*2;
    struct Source *source_x = model_x->source[source_id];
    struct Source *source_y = model_y->source[source_id];
    
    //subtract current nth source from  model
    for(i=0; i<source_x->BW; i++)
    {
        j = i+source_x->imin;
        
        if(j>-1 && j<data->N)
        {
            int i_re = 2*i;
            int i_im = i_re+1;
            int j_re = 2*j;
            int j_im = j_re+1;
            
            switch(data->Nchannel)
            {
                case 1:
                    model_y->tdi->X[j_re] -= source_x->tdi->X[i_re];
                    model_y->tdi->X[j_im] -= source_x->tdi->X[i_im];
                    break;
                case 2:
                    model_y->tdi->A[j_re] -= source_x->tdi->A[i_re];
                    model_y->tdi->A[j_im] -= source_x->tdi->A[i_im];
                    
                    model_y->tdi->E[j_re] -= source_x->tdi->E[i_re];
                    model_y->tdi->E[j_im] -= source_x->tdi->E[i_im];
                    break;
                case 3:
                    model_y->tdi->X[j_re] -= source_x->tdi->X[i_re];
                    model_y->tdi->X[j_im] -= source_x->tdi->X[i_im];
                    
                    model_y->tdi->Y[j_re] -= source_x->tdi->Y[i_re];
                    model_y->tdi->Y[j_im] -= source_x->tdi->Y[i_im];

                    model_y->tdi->Z[j_re] -= source_x->tdi->Z[i_re];
                    model_y->tdi->Z[j_im] -= source_x->tdi->Z[i_im];
                    break;
            }
        }//check that source_id is in range
    }//loop over waveform bins

    //generate proposed signal model
    for(i=0; i<N2; i++)
    {
        switch(data->Nchannel)
        {
            case 1:
                source_y->tdi->X[i]=0.0;
                break;
            case 2:
                source_y->tdi->A[i]=0.0;
                source_y->tdi->E[i]=0.0;
                break;
            case 3:
                source_y->tdi->X[i]=0.0;
                source_y->tdi->Y[i]=0.0;
                source_y->tdi->Z[i]=0.0;
                break;
        }
    }

    map_array_to_params(source_y, source_y->params, data->T);
    galactic_binary_alignment(orbit, data, source_y);

    galactic_binary(orbit, data->format, data->T, model_y->t0, source_y->params, NP, source_y->tdi->X,source_y->tdi->Y,source_y->tdi->Z, source_y->tdi->A, source_y->tdi->E, source_y->BW, source_y->tdi->Nchannel);

    //subtract proposed nth source to model
    for(i=0; i<source_y->BW; i++)
    {
        j = i+source_y->imin;
        
        if(j>-1 && j<data->N)
        {
            int i_re = 2*i;
            int i_im = i_re+1;
            int j_re = 2*j;
            int j_im = j_re+1;
            
            switch(data->Nchannel)
            {
                case 1:
                    model_y->tdi->X[j_re] += source_y->tdi->X[i_re];
                    model_y->tdi->X[j_im] += source_y->tdi->X[i_im];
                    break;
                case 2:
                    model_y->tdi->A[j_re] += source_y->tdi->A[i_re];
                    model_y->tdi->A[j_im] += source_y->tdi->A[i_im];
                    
                    model_y->tdi->E[j_re] += source_y->tdi->E[i_re];
                    model_y->tdi->E[j_im] += source_y->tdi->E[i_im];
                    break;
                case 3:
                    model_y->tdi->X[j_re] += source_y->tdi->X[i_re];
                    model_y->tdi->X[j_im] += source_y->tdi->X[i_im];
                    
                    model_y->tdi->Y[j_re] += source_y->tdi->Y[i_re];
                    model_y->tdi->Y[j_im] += source_y->tdi->Y[i_im];
                    
                    model_y->tdi->Z[j_re] += source_y->tdi->Z[i_re];
                    model_y->tdi->Z[j_im] += source_y->tdi->Z[i_im];
                    break;
            }
        }//check that source_id is in range
    }//loop over waveform bins
}

void generate_noise_model(struct Data *data, struct Model *model)
{
    for(int n=0; n<data->N; n++)
    {
        model->noise->f[n] = data->fmin + (double)n/data->T;

        for(int i=0; i<data->Nchannel; i++)
            for(int j=i; j<data->Nchannel; j++)
                model->noise->C[i][j][n] = model->noise->C[j][i][n] = data->noise->C[i][j][n]*sqrt(model->noise->eta[i]*model->noise->eta[j]);

    }
    invert_noise_covariance_matrix(model->noise);
}

void generate_calibration_model(struct Data *data, struct Model *model)
{
    struct Calibration *calibration = model->calibration;
    switch(data->Nchannel)
    {
        case 1:
            calibration->real_dphiX = cos(calibration->dphiX);
            calibration->imag_dphiX = sin(calibration->dphiX);
            break;
        case 2:
            calibration->real_dphiA = cos(calibration->dphiA);
            calibration->imag_dphiA = sin(calibration->dphiA);
            
            calibration->real_dphiE = cos(calibration->dphiE);
            calibration->imag_dphiE = sin(calibration->dphiE);
            break;
        case 3:
            calibration->real_dphiX = cos(calibration->dphiX);
            calibration->imag_dphiX = sin(calibration->dphiX);
            
            calibration->real_dphiY = cos(calibration->dphiY);
            calibration->imag_dphiY = sin(calibration->dphiY);

            calibration->real_dphiZ = cos(calibration->dphiZ);
            calibration->imag_dphiZ = sin(calibration->dphiZ);
            break;
        default:
            break;
    }
}

void apply_calibration_model(struct Data *data, struct Model *model)
{
    double dA;
    double cal_re;
    double cal_im;
    double h_re;
    double h_im;
    int i_re;
    int i_im;
    int i;
    
//apply calibration error to full signal model
    for(i=0; i<data->N; i++)
    {
        i_re = 2*i;
        i_im = i_re+1;
        
        switch(data->Nchannel)
        {
            case 1:
                h_re = model->tdi->X[i_re];
                h_im = model->tdi->X[i_im];
                
                dA     = (1.0 + model->calibration->dampX);
                cal_re = model->calibration->real_dphiX;
                cal_im = model->calibration->imag_dphiX;
                
                model->tdi->X[i_re] = dA*(h_re*cal_re - h_im*cal_im);
                model->tdi->X[i_im] = dA*(h_re*cal_im + h_im*cal_re);
                break;
            case 2:
                h_re = model->tdi->A[i_re];
                h_im = model->tdi->A[i_im];
                
                dA     = (1.0 + model->calibration->dampA);
                cal_re = model->calibration->real_dphiA;
                cal_im = model->calibration->imag_dphiA;
                
                model->tdi->A[i_re] = dA*(h_re*cal_re - h_im*cal_im);
                model->tdi->A[i_im] = dA*(h_re*cal_im + h_im*cal_re);
                
                
                h_re = model->tdi->E[i_re];
                h_im = model->tdi->E[i_im];
                
                dA     = (1.0 + model->calibration->dampE);
                cal_re = model->calibration->real_dphiE;
                cal_im = model->calibration->imag_dphiE;
                
                model->tdi->E[i_re] = dA*(h_re*cal_re - h_im*cal_im);
                model->tdi->E[i_im] = dA*(h_re*cal_im + h_im*cal_re);
                break;
            case 3:
                h_re = model->tdi->X[i_re];
                h_im = model->tdi->X[i_im];
                
                dA     = (1.0 + model->calibration->dampX);
                cal_re = model->calibration->real_dphiX;
                cal_im = model->calibration->imag_dphiX;
                
                model->tdi->X[i_re] = dA*(h_re*cal_re - h_im*cal_im);
                model->tdi->X[i_im] = dA*(h_re*cal_im + h_im*cal_re);
                
                
                h_re = model->tdi->Y[i_re];
                h_im = model->tdi->Y[i_im];
                
                dA     = (1.0 + model->calibration->dampY);
                cal_re = model->calibration->real_dphiY;
                cal_im = model->calibration->imag_dphiY;
                
                model->tdi->Y[i_re] = dA*(h_re*cal_re - h_im*cal_im);
                model->tdi->Y[i_im] = dA*(h_re*cal_im + h_im*cal_re);

                h_re = model->tdi->Z[i_re];
                h_im = model->tdi->Z[i_im];
                
                dA     = (1.0 + model->calibration->dampZ);
                cal_re = model->calibration->real_dphiZ;
                cal_im = model->calibration->imag_dphiZ;
                
                model->tdi->Z[i_re] = dA*(h_re*cal_re - h_im*cal_im);
                model->tdi->Z[i_im] = dA*(h_re*cal_im + h_im*cal_re);
                break;
            default:
                break;
        }//end switch
    }//end loop over data
}

void maximize_signal_model(struct Orbit *orbit, struct Data *data, struct Model *model, int source_id)
{
    if(source_id < model->Nlive)
    {
        double *Fparams = calloc(NP,sizeof(double));
        
        struct Source *source = model->source[source_id];
        
        /* save original data */
        //    struct TDI *data_save = malloc(sizeof(struct TDI));
        //    alloc_tdi(data_save, data->N, data->Nchannel);
        //    copy_tdi(data->tdi[FIXME],data_save);
        //
        //    /* save original noise */
        //    struct Noise *noise_save = malloc(sizeof(struct Noise));
        //    alloc_noise(noise_save, data->N);
        //    copy_noise(data->noise[FIXME], noise_save);
        //
        //    /* put current noise model in data structure for get_Fstat_logL() */
        //    copy_noise(model->noise[FIXME],data->noise[FIXME]);
        
        
        /* create residual of all sources but n for F-statistic */
        //    for(int m=0; m<model->Nlive; m++)
        //    {
        //        if(m!=source_id)
        //        {
        //            for(int i=0; i<data->N*2; i++)
        //            {
        //                data->tdi[FIXME]->A[i] -= model->source[m]->tdi->A[i];
        //                data->tdi[FIXME]->E[i] -= model->source[m]->tdi->E[i];
        //                data->tdi[FIXME]->X[i] -= model->source[m]->tdi->X[i];
        //            }
        //        }
        //    }
        
        
        if(!check_range(source->params, model->prior))
        {
            /* maximize parameters w/ F-statistic */
            get_Fstat_xmax(orbit, data, source->params, Fparams);
            
            /* unpack maximized parameters */
            source->amp  = exp(Fparams[3]);
            source->cosi = Fparams[4];
            source->psi  = Fparams[5];
            source->phi0 = Fparams[6];
            map_params_to_array(source, source->params, data->T);
        }
        
        /* restore original data */
        //    copy_tdi(data_save,data->tdi[FIXME]);
        //    free_tdi(data_save);
        //
        //    /* restore original noise */
        //    copy_noise(noise_save,data->noise[FIXME]);
        //    free_noise(noise_save);
        free(Fparams);
    }
}
double gaussian_log_likelihood(struct Data *data, struct Model *model)
{
    
    /*
    *
    * Form residual and sum
    *
    */
    
    int N2 = data->N*2;
    double chi2 = 0.0; //chi^2, where logL = -chi^2 / 2
    
    //loop over time segments
    struct TDI *residual = model->residual;
    
    for(int i=0; i<N2; i++)
    {
        residual->X[i] = data->tdi->X[i] - model->tdi->X[i];
        residual->Y[i] = data->tdi->Y[i] - model->tdi->Y[i];
        residual->Z[i] = data->tdi->Z[i] - model->tdi->Z[i];
        residual->A[i] = data->tdi->A[i] - model->tdi->A[i];
        residual->E[i] = data->tdi->E[i] - model->tdi->E[i];
    }
            
    switch(data->Nchannel)
    {
        case 1:
            chi2 += fourier_nwip(residual->X, residual->X, model->noise->invC[0][0], data->N);
            break;
        case 2:
            chi2 += fourier_nwip(residual->A, residual->A, model->noise->invC[0][0], data->N);
            chi2 += fourier_nwip(residual->E, residual->E, model->noise->invC[1][1], data->N);
            break;
        case 3:
            chi2 += fourier_nwip(residual->X, residual->X, model->noise->invC[0][0], data->N);
            chi2 += fourier_nwip(residual->Y, residual->Y, model->noise->invC[1][1], data->N);
            chi2 += fourier_nwip(residual->Z, residual->Z, model->noise->invC[2][2], data->N);

            chi2 += 2.0*fourier_nwip(residual->X, residual->Y, model->noise->invC[0][1], data->N);
            chi2 += 2.0*fourier_nwip(residual->X, residual->Z, model->noise->invC[0][2], data->N);
            chi2 += 2.0*fourier_nwip(residual->Y, residual->Z, model->noise->invC[1][2], data->N);

            break;
        default:
            fprintf(stderr,"Unsupported number of channels in gaussian_log_likelihood()\n");
            exit(1);
    }

    return -0.5*chi2;
}

double gaussian_log_likelihood_constant_norm(struct Data *data, struct Model *model)
{
    
    double logLnorm = 0.0;
    
    //loop over time segments
    logLnorm -= (double)data->N*log(model->noise->detC[0]);
    
    return logLnorm;
}

double gaussian_log_likelihood_model_norm(struct Data *data, struct Model *model)
{
    
    double logLnorm = 0.0;
    
    //loop over time segments
    for(int n=0; n<data->N; n++)
        logLnorm -= log(model->noise->detC[n]);

    return logLnorm;
}

double delta_log_likelihood(struct Data *data, struct Model *model_x, struct Model *model_y, int source_id)
{
    /*
    *
    * Update residual and only sum over affected bins
    *
    */

    double deltalogL = 0.0;
    struct Source *source_x = model_x->source[source_id];
    struct Source *source_y = model_y->source[source_id];

    //loop over time segments
    struct TDI *residual_x = model_x->residual;
    struct TDI *residual_y = model_y->residual;
    
    //add current source back into residual
    for(int i=0; i<source_x->BW; i++)
    {
        int j = i+source_x->imin;
        if(j>-1 && j<data->N)
        {
            int i_re = 2*i;
            int i_im = i_re+1;
            int j_re = 2*j;
            int j_im = j_re+1;
            
            residual_y->X[i_re] += source_x->tdi->X[i_re];
            residual_y->X[i_im] += source_x->tdi->X[i_im];
            residual_y->Y[i_re] += source_x->tdi->Y[i_re];
            residual_y->Y[i_im] += source_x->tdi->Y[i_im];
            residual_y->Z[i_re] += source_x->tdi->Z[i_re];
            residual_y->Z[i_im] += source_x->tdi->Z[i_im];
            residual_y->A[j_re] += source_x->tdi->A[i_re];
            residual_y->A[j_im] += source_x->tdi->A[i_im];
            residual_y->E[j_re] += source_x->tdi->E[i_re];
            residual_y->E[j_im] += source_x->tdi->E[i_im];
        }
    }

    //form proposed residual by subtracting new source
    for(int i=0; i<source_y->BW; i++)
    {
        int j = i+source_y->imin;
        if(j>-1 && j<data->N)
        {
            int i_re = 2*i;
            int i_im = i_re+1;
            int j_re = 2*j;
            int j_im = j_re+1;
            
            residual_y->X[i_re] -= source_y->tdi->X[i_re];
            residual_y->X[i_im] -= source_y->tdi->X[i_im];
            residual_y->Y[i_re] -= source_y->tdi->Y[i_re];
            residual_y->Y[i_im] -= source_y->tdi->Y[i_im];
            residual_y->Z[i_re] -= source_y->tdi->Z[i_re];
            residual_y->Z[i_im] -= source_y->tdi->Z[i_im];
            residual_y->A[j_re] -= source_y->tdi->A[i_re];
            residual_y->A[j_im] -= source_y->tdi->A[i_im];
            residual_y->E[j_re] -= source_y->tdi->E[i_re];
            residual_y->E[j_im] -= source_y->tdi->E[i_im];
        }
    }

    //find range of integration
    int imin = find_min(source_x->imin,source_y->imin);
    int imax = find_max(source_y->imin+source_y->BW,source_x->imin+source_x->BW);
    
    //keep it in bounds
    if(imax>data->N)imax=data->N;
    if(imin<0)imin=0;

    //the complex array elements to skip in the sum
    int skip=2*imin;
    
    switch(data->Nchannel)
    {
        case 1:
            deltalogL -= fourier_nwip(residual_x->X+skip, residual_x->X+skip, model_x->noise->invC[0][0]+imin, imax-imin);
            deltalogL += fourier_nwip(residual_y->X+skip, residual_y->X+skip, model_x->noise->invC[0][0]+imin, imax-imin);
            break;
            
        case 2:
            deltalogL -= fourier_nwip(residual_x->A+skip, residual_x->A+skip, model_x->noise->invC[0][0]+imin, imax-imin);
            deltalogL -= fourier_nwip(residual_x->E+skip, residual_x->E+skip, model_x->noise->invC[1][1]+imin, imax-imin);

            deltalogL += fourier_nwip(residual_y->A+skip, residual_y->A+skip, model_y->noise->invC[0][0]+imin, imax-imin);
            deltalogL += fourier_nwip(residual_y->E+skip, residual_y->E+skip, model_y->noise->invC[1][1]+imin, imax-imin);

            break;
        case 3:
            deltalogL -= fourier_nwip(residual_x->X+skip, residual_x->X+skip, model_x->noise->invC[0][0]+imin, imax-imin);
            deltalogL -= fourier_nwip(residual_x->Y+skip, residual_x->Y+skip, model_x->noise->invC[1][1]+imin, imax-imin);
            deltalogL -= fourier_nwip(residual_x->Z+skip, residual_x->Z+skip, model_x->noise->invC[2][2]+imin, imax-imin);
            deltalogL -= 2.0*fourier_nwip(residual_x->X+skip, residual_x->Y+skip, model_x->noise->invC[0][1]+imin, imax-imin);
            deltalogL -= 2.0*fourier_nwip(residual_x->X+skip, residual_x->Z+skip, model_x->noise->invC[0][2]+imin, imax-imin);
            deltalogL -= 2.0*fourier_nwip(residual_x->Y+skip, residual_x->Z+skip, model_x->noise->invC[1][2]+imin, imax-imin);

            deltalogL += fourier_nwip(residual_y->X+skip, residual_y->X+skip, model_y->noise->invC[0][0]+imin, imax-imin);
            deltalogL += fourier_nwip(residual_y->Y+skip, residual_y->Y+skip, model_y->noise->invC[1][1]+imin, imax-imin);
            deltalogL += fourier_nwip(residual_y->Z+skip, residual_y->Z+skip, model_y->noise->invC[2][2]+imin, imax-imin);
            deltalogL += 2.0*fourier_nwip(residual_y->X+skip, residual_y->Y+skip, model_y->noise->invC[0][1]+imin, imax-imin);
            deltalogL += 2.0*fourier_nwip(residual_y->X+skip, residual_y->Z+skip, model_y->noise->invC[0][2]+imin, imax-imin);
            deltalogL += 2.0*fourier_nwip(residual_y->Y+skip, residual_y->Z+skip, model_y->noise->invC[1][2]+imin, imax-imin);

            break;
        default:
            fprintf(stderr,"Unsupported number of channels in delta_log_likelihood()\n");
            exit(1);
    }
        
    return -0.5*deltalogL;

}

int update_max_log_likelihood(struct Model **model, struct Chain *chain, struct Flags *flags)
{
    int n = chain->index[0];
    int N = model[n]->Nlive;
    
    double logL = 0.0;
    double dlogL= 0.0;
    
    // get full likelihood
    logL = model[n]->logL + model[n]->logLnorm;
    
    // update max
    if(logL > chain->logLmax)
    {
        dlogL = logL - chain->logLmax;
        
        //clone chains if new mode is found (dlogL > D/2)
        if( dlogL > (double)(8*N/2) )
        {
            chain->logLmax = logL;
            
            for(int ic=1; ic<chain->NC; ic++)
            {
                int m = chain->index[ic];
                copy_model(model[n],model[m]);
            }
            if(flags->burnin)return 1;
        }
    }
    
    return 0;
}

