/*
 * Copyright 2019 Tyson B. Littenberg & Neil J. Cornish 
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <glass_utils.h>

#include "glass_ucb_model.h"
#include "glass_ucb_catalog.h"
#include "glass_ucb_io.h"
#include "glass_ucb_prior.h"
#include "glass_ucb_waveform.h"
#include "glass_ucb_fstatistic.h"

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
    if(UCB_MODEL_NP>7)
        source->dfdt   = params[7]/(T*T);
    if(UCB_MODEL_NP>8)
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
    if(UCB_MODEL_NP>7)
        params[7] = source->dfdt*T*T;
    if(UCB_MODEL_NP>8)
        params[8] = source->d2fdt2*T*T*T;
}

void alloc_model(struct Data *data, struct Model *model, int Nmax)
{
    int n;
    
    model->Nlive = 1;
    model->Nmax  = Nmax;
    model->Neff  = 2;
    model->t0 = 0.0;

    //Wavelet bookkeeping
    model->Nlist = 0;
    model->list = calloc(data->N,sizeof(int));

    model->source = malloc(model->Nmax*sizeof(struct Source *));
    
    model->calibration = malloc(sizeof(struct Calibration) );
    model->noise       = malloc(sizeof(struct Noise)       );
    model->tdi         = malloc(sizeof(struct TDI)         );
    model->residual    = malloc(sizeof(struct TDI)         );
    
    if(!strcmp(data->basis,"fourier")) alloc_noise(model->noise,data->NFFT, data->Nlayer, data->Nchannel);
    if(!strcmp(data->basis,"wavelet")) alloc_noise(model->noise,data->N, data->Nlayer, data->Nchannel);

    alloc_tdi(model->tdi, data->N, data->Nchannel);
    alloc_tdi(model->residual, data->N, data->Nchannel);
    alloc_calibration(model->calibration);
    
    for(n=0; n<model->Nmax; n++)
    {
        model->source[n] = malloc(sizeof(struct Source));
        alloc_source(model->source[n],data->N,data->Nchannel);
    }
    
    model->logPriorVolume = calloc(UCB_MODEL_NP,sizeof(double));
    model->prior = malloc(UCB_MODEL_NP*sizeof(double *));
    for(n=0; n<UCB_MODEL_NP; n++) model->prior[n] = calloc(2,sizeof(double));

}

void copy_model(struct Model *origin, struct Model *copy)
{
    //Source parameters
    copy->Nmax  = origin->Nmax;
    copy->Neff  = origin->Neff;
    copy->Nlive = origin->Nlive;
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
    for(int n=0; n<UCB_MODEL_NP; n++)
    {
        for(int j=0; j<2; j++) copy->prior[n][j] = origin->prior[n][j];
        copy->logPriorVolume[n] = origin->logPriorVolume[n];
    }
    //Model likelihood
    copy->logL           = origin->logL;
    copy->logLnorm       = origin->logLnorm;

    //Wavelet bookkeeping
    copy->Nlist = origin->Nlist;
    memcpy(copy->list,origin->list,origin->Nlist*sizeof(int));

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

    //Wavelet bookkeeping
    copy->Nlist = origin->Nlist;
    memcpy(copy->list,origin->list,origin->Nlist*sizeof(int));
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
        for(int j=0; j<UCB_MODEL_NP; j++) if(sa->params[j] != sb->params[j]) return 1;
        
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

    for(n=0; n<UCB_MODEL_NP; n++) free(model->prior[n]);
    free(model->prior);
    free(model->logPriorVolume);

    free_tdi(model->tdi);
    free_tdi(model->residual);
    free_noise(model->noise);
    free_calibration(model->calibration);

    free(model->list);
    
    free(model);
}

void alloc_source(struct Source *source, int N, int Nchannel)
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
    source->BW   = N/2;
    source->qmin = 0;
    source->qmax = N/2;
    source->imin = 0;
    source->imax = N/2;
    
    
    //Package parameters for waveform generator
    source->params=calloc(UCB_MODEL_NP,sizeof(double));
    
    //Response
    source->tdi = malloc(sizeof(struct TDI));
    alloc_tdi(source->tdi,N, Nchannel);
    
    //Fisher
    source->fisher_matrix = malloc(UCB_MODEL_NP*sizeof(double *));
    source->fisher_evectr = malloc(UCB_MODEL_NP*sizeof(double *));
    source->fisher_evalue = calloc(UCB_MODEL_NP,sizeof(double));
    for(int i=0; i<UCB_MODEL_NP; i++)
    {
        source->fisher_matrix[i] = calloc(UCB_MODEL_NP,sizeof(double));
        source->fisher_evectr[i] = calloc(UCB_MODEL_NP,sizeof(double));
    }

    //Wavelet bookkeeping
    source->Nlist = 0;
    source->list = calloc(N,sizeof(int));
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
    memcpy(copy->fisher_evalue, origin->fisher_evalue, UCB_MODEL_NP*sizeof(double));
    memcpy(copy->params, origin->params, UCB_MODEL_NP*sizeof(double));
    copy->fisher_update_flag = origin->fisher_update_flag;
    
    for(int i=0; i<UCB_MODEL_NP; i++)
    {
        memcpy(copy->fisher_matrix[i], origin->fisher_matrix[i], UCB_MODEL_NP*sizeof(double));
        memcpy(copy->fisher_evectr[i], origin->fisher_evectr[i], UCB_MODEL_NP*sizeof(double));
    }

    copy->Nlist = origin->Nlist;
    memcpy(copy->list, origin->list, origin->Nlist*sizeof(int));
    
}

void free_source(struct Source *source)
{
    for(int i=0; i<UCB_MODEL_NP; i++)
    {
        free(source->fisher_matrix[i]);
        free(source->fisher_evectr[i]);
    }
    free(source->fisher_matrix);
    free(source->fisher_evectr);
    free(source->fisher_evalue);
    free(source->params);
    
    free_tdi(source->tdi);

    free(source->list);
    
    free(source);
}

void generate_signal_model(struct Orbit *orbit, struct Data *data, struct Model *model, int source_id)
{
    int i,n;

    for(n=0; n<data->N; n++)
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

            for(i=0; i<data->N; i++)
            {
                source->tdi->X[i]=0.0;
                source->tdi->Y[i]=0.0;
                source->tdi->Z[i]=0.0;
                source->tdi->A[i]=0.0;
                source->tdi->E[i]=0.0;
            }

            map_array_to_params(source, source->params, data->T);
            

            //Book-keeping of injection time-frequency volume
            ucb_alignment(orbit, data, source);

        }
        
        //Simulate gravitational wave signal
        /* the source_id = -1 condition is redundent if the model->tdi structure is up to date...*/

        if(source_id==-1 || source_id==n) ucb_waveform(orbit, data->format, data->T, model->t0, source->params, UCB_MODEL_NP, source->tdi->X, source->tdi->Y, source->tdi->Z, source->tdi->A, source->tdi->E, source->BW, source->tdi->Nchannel);
        
        //Add waveform to model TDI channels
        add_signal_model(data,model,source);

        
    }//loop over sources
}

void generate_signal_model_wavelet(struct Orbit *orbit, struct Data *data, struct Model *model, int source_id)
{
    int i,n;
    
    for(n=0; n<data->N; n++)
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
            for(i=0; i<data->N; i++)
            {
                source->tdi->X[i]=0.0;
                source->tdi->Y[i]=0.0;
                source->tdi->Z[i]=0.0;
                source->tdi->A[i]=0.0;
                source->tdi->E[i]=0.0;
                source->list[i]=0;
            }
            
            map_array_to_params(source, source->params, data->T);            
        }
        
        //Simulate gravitational wave signal
        /* the source_id = -1 condition is redundent if the model->tdi structure is up to date...*/
        if(source_id==-1 || source_id==n) ucb_waveform_wavelet(orbit, data->wdm, data->T, model->t0, source->params, source->list, &source->Nlist, source->tdi->X, source->tdi->Y, source->tdi->Z);

        //Add waveform to model TDI channels
        add_signal_model_wavelet(data,model,source);

    }//loop over sources
}

void remove_signal_model(struct Data *data, struct Model *model, struct Source *source)
{
    //subtract current nth source from  model
    for(int i=0; i<source->BW; i++)
    {
        int j = i+source->imin;
        
        if(j>-1 && j<data->NFFT)
        {
            int i_re = 2*i;
            int i_im = i_re+1;
            int j_re = 2*j;
            int j_im = j_re+1;
            
            switch(data->Nchannel)
            {
                case 1:
                    model->tdi->X[j_re] -= source->tdi->X[i_re];
                    model->tdi->X[j_im] -= source->tdi->X[i_im];
                    break;
                case 2:
                    model->tdi->A[j_re] -= source->tdi->A[i_re];
                    model->tdi->A[j_im] -= source->tdi->A[i_im];
                    
                    model->tdi->E[j_re] -= source->tdi->E[i_re];
                    model->tdi->E[j_im] -= source->tdi->E[i_im];
                    break;
                case 3:
                    model->tdi->X[j_re] -= source->tdi->X[i_re];
                    model->tdi->X[j_im] -= source->tdi->X[i_im];
                    
                    model->tdi->Y[j_re] -= source->tdi->Y[i_re];
                    model->tdi->Y[j_im] -= source->tdi->Y[i_im];

                    model->tdi->Z[j_re] -= source->tdi->Z[i_re];
                    model->tdi->Z[j_im] -= source->tdi->Z[i_im];
                    break;
            }
        }//check that source_id is in range
    }//loop over waveform bins
}


void add_signal_model(struct Data *data, struct Model *model, struct Source *source)
{
    //add current nth source to  model
    for(int i=0; i<source->BW; i++)
    {
        int j = i+source->imin;
        
        if(j>-1 && j<data->NFFT)
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
}
void add_signal_model_wavelet(struct Data *data, struct Model *model, struct Source *source)
{
    //insert source into model 
    for(int n=0; n<source->Nlist; n++)
    {
        int k = source->list[n];
        model->tdi->X[k] += source->tdi->X[k];
        model->tdi->Y[k] += source->tdi->Y[k];
        model->tdi->Z[k] += source->tdi->Z[k];
    }

    //get union of list
    list_union(model->list, source->list, model->Nlist, source->Nlist, model->list, &model->Nlist); 
}

void remove_signal_model_wavelet(struct Data *data, struct Model *model, struct Source *source)
{
    //insert source into model 
    for(int n=0; n<source->Nlist; n++)
    {
        int k=source->list[n];
        if(k>=0 && k<data->N)
        {
            model->tdi->X[k] -= source->tdi->X[k];
            model->tdi->Y[k] -= source->tdi->Y[k];
            model->tdi->Z[k] -= source->tdi->Z[k];
        }
    }

    //get union of list
    if(model->Nlive == 0) model->Nlist = 0;
    else
    {
        for(int n=0; n<model->Nlive; n++)
            list_union(model->list, model->source[n]->list, model->Nlist, model->source[n]->Nlist, model->list, &model->Nlist); 
    }
}

void update_signal_model(struct Orbit *orbit, struct Data *data, struct Model *model_x, struct Model *model_y, int source_id)
{
    int i;
    struct Source *source_x = model_x->source[source_id];
    struct Source *source_y = model_y->source[source_id];
    
    //subtract current nth source from  model
    remove_signal_model(data,model_y,source_x);

    //generate proposed signal model
    for(i=0; i<data->N; i++)
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
    ucb_alignment(orbit, data, source_y);

    ucb_waveform(orbit, data->format, data->T, model_y->t0, source_y->params, UCB_MODEL_NP, source_y->tdi->X,source_y->tdi->Y,source_y->tdi->Z, source_y->tdi->A, source_y->tdi->E, source_y->BW, source_y->tdi->Nchannel);

    //add proposed nth source to model
    add_signal_model(data,model_y,source_y);
    
}

void update_signal_model_wavelet(struct Orbit *orbit, struct Data *data, struct Model *model_x, struct Model *model_y, int source_id)
{
    int i;
    int N=data->N;
    struct Source *source_x = model_x->source[source_id];
    struct Source *source_y = model_y->source[source_id];
    
    //subtract current nth source from  model
    remove_signal_model_wavelet(data,model_y,source_x);

    //generate proposed signal model
    for(i=0; i<N; i++)
    {
        source_y->tdi->X[i]=0.0;
        source_y->tdi->Y[i]=0.0;
        source_y->tdi->Z[i]=0.0;
    }

    map_array_to_params(source_y, source_y->params, data->T);
    ucb_waveform_wavelet(orbit, data->wdm, data->T, model_y->t0, source_y->params, source_y->list, &source_y->Nlist, source_y->tdi->X, source_y->tdi->Y, source_y->tdi->Z);

    //add proposed nth source to model
    add_signal_model_wavelet(data,model_y,source_y);
    
}

void generate_noise_model(struct Data *data, struct Model *model)
{
    for(int n=0; n<data->NFFT; n++)
    {
        for(int i=0; i<data->Nchannel; i++)
            for(int j=i; j<data->Nchannel; j++)
                model->noise->C[i][j][n] = model->noise->C[j][i][n] = data->noise->C[i][j][n]*sqrt(model->noise->eta[i]*model->noise->eta[j]);

    }
    invert_noise_covariance_matrix(model->noise);
}

void generate_noise_model_wavelet(struct Data *data, struct Model *model)
{
    int Nlayers = data->Nlayer;         //number of frequency layers
    int Nslices = data->N/data->Nlayer; //number of time slices
    for(int n=0; n<Nlayers; n++)
    {
        for(int m=0; m<Nslices; m++)
        {
            int k = n*Nslices+m;
            for(int i=0; i<data->Nchannel; i++)
            {
                for(int j=i; j<data->Nchannel; j++)
                {
                    model->noise->C[i][j][k] = model->noise->C[j][i][k] = data->noise->C[i][j][k]*sqrt(model->noise->eta[i*Nlayers+n]*model->noise->eta[j*Nlayers+n]);
                }
            }
        }

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
        double *Fparams = calloc(UCB_MODEL_NP,sizeof(double));
        
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
    
    double chi2 = 0.0; //chi^2, where logL = -chi^2 / 2
    
    //loop over time segments
    struct TDI *residual = model->residual;
    
    for(int i=0; i<data->N; i++)
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
            chi2 += fourier_nwip(residual->X, residual->X, model->noise->invC[0][0], data->NFFT);
            break;
        case 2:
            chi2 += fourier_nwip(residual->A, residual->A, model->noise->invC[0][0], data->NFFT);
            chi2 += fourier_nwip(residual->E, residual->E, model->noise->invC[1][1], data->NFFT);
            break;
        case 3:
            chi2 += fourier_nwip(residual->X, residual->X, model->noise->invC[0][0], data->NFFT);
            chi2 += fourier_nwip(residual->Y, residual->Y, model->noise->invC[1][1], data->NFFT);
            chi2 += fourier_nwip(residual->Z, residual->Z, model->noise->invC[2][2], data->NFFT);

            chi2 += 2.0*fourier_nwip(residual->X, residual->Y, model->noise->invC[0][1], data->NFFT);
            chi2 += 2.0*fourier_nwip(residual->X, residual->Z, model->noise->invC[0][2], data->NFFT);
            chi2 += 2.0*fourier_nwip(residual->Y, residual->Z, model->noise->invC[1][2], data->NFFT);

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
    if(!strcmp(data->basis,"fourier")) logLnorm -= (double)data->NFFT*log(model->noise->detC[0]);
    if(!strcmp(data->basis,"wavelet")) logLnorm -= 0.5*(double)data->N*log(model->noise->detC[0]);
    
    return logLnorm;
}

double gaussian_log_likelihood_model_norm(struct Data *data, struct Model *model)
{
    
    double logLnorm = 0.0;
    int N;
    if(!strcmp(data->basis,"fourier")) N=data->NFFT;
    if(!strcmp(data->basis,"wavelet")) N=data->N;

    //loop over time segments
    for(int n=0; n<N; n++)
        logLnorm -= log(model->noise->detC[n]);

    if(!strcmp(data->basis,"wavelet")) logLnorm *= 0.5;  //normalization of 1/2 for wavelet domain (sum over N, not N/2)
    
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
        if(j>-1 && j<data->NFFT)
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
        if(j>-1 && j<data->NFFT)
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
    if(imax>data->NFFT)imax=data->NFFT;
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

double gaussian_log_likelhood_wavelet(struct Data *data, struct Model *model)
{
    /*
    Form residual and sum
    */

    double chi2 = 0.0;

    struct TDI *residual = model->residual;

    
    for(int n=0; n<data->N; n++)
    {
        residual->X[n] = data->tdi->X[n];
        residual->Y[n] = data->tdi->Y[n];
        residual->Z[n] = data->tdi->Z[n];
    }
    
    for(int n=0; n<model->Nlist; n++)
    {
        int k = model->list[n];
        if(k>=0 && k<data->N)
        {
            residual->X[k] -= model->tdi->X[k];
            residual->Y[k] -= model->tdi->Y[k];
            residual->Z[k] -= model->tdi->Z[k];
        }
    }
    
    int *list = int_vector(data->N);
    for(int n=0; n<data->N; n++) list[n]=n;

    chi2 += wavelet_nwip(residual->X, residual->X, model->noise->invC[0][0], list, data->N);
    chi2 += wavelet_nwip(residual->Y, residual->Y, model->noise->invC[1][1], list, data->N);
    chi2 += wavelet_nwip(residual->Z, residual->Z, model->noise->invC[2][2], list, data->N);
    chi2 += wavelet_nwip(residual->X, residual->Y, model->noise->invC[0][1], list, data->N)*2;
    chi2 += wavelet_nwip(residual->X, residual->Z, model->noise->invC[0][2], list, data->N)*2;
    chi2 += wavelet_nwip(residual->Y, residual->Z, model->noise->invC[1][2], list, data->N)*2;

    free_int_vector(list);
    
    return -0.5*chi2;

}
