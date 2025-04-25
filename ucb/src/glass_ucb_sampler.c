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

#include <glass_utils.h>
#include <glass_noise.h>

#include "glass_ucb_model.h"
#include "glass_ucb_catalog.h"
#include "glass_ucb_io.h"
#include "glass_ucb_data.h"
#include "glass_ucb_prior.h"
#include "glass_ucb_proposal.h"
#include "glass_ucb_waveform.h"
#include "glass_ucb_sampler.h"

static double maximization_penalty(int dimension, int bandwidth)
{
    //return -0.5*(double)dimension*log(2.*(double)bandwidth); //bic-inspired = -klog(N)/2
    return -(double)dimension; //likelihood-inspired = -k
    //return 0.0; //no penalty
}

void ptmcmc(struct Model **model, struct Chain *chain, struct Flags *flags)
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
    
    //b = (int)(ran2(seed)*((double)(chain->NC-1)));
    for(b=NC-1; b>0; b--)
    {
        a = b - 1;
        chain->acceptance[a]=0;
        
        olda = chain->index[a];
        oldb = chain->index[b];
        
        heat1 = chain->temperature[a];
        heat2 = chain->temperature[b];
        
        logL1 = model[olda]->logL + model[olda]->logLnorm;
        logL2 = model[oldb]->logL + model[oldb]->logLnorm;
        
        //Hot chains jump more rarely
        if(rand_r_U_0_1(&chain->r[a])<1.0)
        {
            dlogL = logL2 - logL1;
            H  = (heat2 - heat1)/(heat2*heat1);
            
            alpha = exp(dlogL*H);
            beta  = rand_r_U_0_1(&chain->r[a]);
            
            if(alpha >= beta)
            {
                chain->index[a] = oldb;
                chain->index[b] = olda;
                chain->acceptance[a]=1;
            }
        }
    }
}

void adapt_temperature_ladder(struct Chain *chain, int mcmc)
{
    int ic;
    
    int NC = chain->NC;
    
    double S[NC];
    double A[NC][2];
    
    double nu=10;
    double t0=10000.;
    
    for(ic=1; ic<NC-1; ic++)
    {
        S[ic] = log(chain->temperature[ic] - chain->temperature[ic-1]);
        A[ic][0] = chain->acceptance[ic-1];
        A[ic][1] = chain->acceptance[ic];
    }
    
    for(ic=1; ic<NC-1; ic++)
    {
        S[ic] += (A[ic][0] - A[ic][1])*(t0/((double)mcmc+t0))/nu;
        
        chain->temperature[ic] = chain->temperature[ic-1] + exp(S[ic]);
        
        if(chain->temperature[ic]/chain->temperature[ic-1] < 1.1) chain->temperature[ic] = chain->temperature[ic-1]*1.1;
    }//end loop over ic
}//end adapt function

void noise_model_mcmc(struct Orbit *orbit, struct Data *data, struct Model *model, struct Model *trial, struct Chain *chain, struct Flags *flags, int ic)
{
    double logH  = 0.0; //(log) Hastings ratio
    double loga  = 1.0; //(log) transition probability
    
    double logPx  = 0.0; //(log) prior density for model x (current state)
    double logPy  = 0.0; //(log) prior density for model y (proposed state)
    
    //shorthand pointers
    struct Model *model_x = model;
    struct Model *model_y = trial;
    
    copy_model(model_x,model_y);
    
    //choose proposal distribution
    for(int n=0; n<data->Nchannel*data->Nlayer; n++)
    {
        model_y->noise->eta[n] = model_x->noise->eta[n] + 0.05*rand_r_N_0_1(&chain->r[ic]);
        if(model_y->noise->eta[n] < 0.01 || model_y->noise->eta[n]>100) logPy=-INFINITY;
    }
    
    
    if(!flags->prior)
    {
        //  Form master template
        if(!strcmp("fourier",data->basis)) generate_noise_model(data, model_y);
        if(!strcmp("wavelet",data->basis)) generate_noise_model_wavelet(data, model_y);
        
        //get likelihood for y
        if(!strcmp("fourier",data->basis)) model_y->logL = gaussian_log_likelihood(data, model_y);
        if(!strcmp("wavelet",data->basis)) model_y->logL = gaussian_log_likelhood_wavelet(data, model_y);

        //model_y->logLnorm = gaussian_log_likelihood_constant_norm(data, model_y);
        model_y->logLnorm = gaussian_log_likelihood_model_norm(data, model_y);

        /*
         H = [p(d|y)/p(d|x)]/T x p(y)/p(x) x q(x|y)/q(y|x)
         */
        logH += ( (model_y->logL+model_y->logLnorm) - (model_x->logL+model_x->logLnorm) )/chain->temperature[ic]; //delta logL
    }
    logH += logPy  - logPx; //priors
    
    loga = log(rand_r_U_0_1(&chain->r[ic]));
    if(logH > loga) copy_model(model_y,model_x);
    
}

void ucb_mcmc(struct Orbit *orbit, struct Data *data, struct Model *model, struct Model *trial, struct Chain *chain, struct Flags *flags, struct Prior *prior, struct Proposal **proposal, int ic)
{
    double logH  = 0.0; //(log) Hastings ratio
    double loga  = 1.0; //(log) transition probability
    
    double logPx  = 0.0; //(log) prior density for model x (current state)
    double logPy  = 0.0; //(log) prior density for model y (proposed state)
    double logQyx = 0.0; //(log) proposal denstiy from x->y
    double logQxy = 0.0; //(log) proposal density from y->x
    
    //shorthand pointers
    struct Model *model_x = model;
    struct Model *model_y = trial;
    
    copy_model(model_x,model_y);

    //pick a source to update
    int n = (int)(rand_r_U_0_1(&chain->r[ic])*(double)model_x->Nlive);
    
    //more shorthand pointers
    struct Source *source_x = model_x->source[n];
    struct Source *source_y = model_y->source[n];
    
    //choose proposal distribution
    double draw;
    int nprop;
    do
    {
        nprop = (int)floor((UCB_PROPOSAL_NPROP)*rand_r_U_0_1(&chain->r[ic]));
        draw = rand_r_U_0_1(&chain->r[ic]);
    }while(proposal[nprop]->weight <= draw);
    
    proposal[nprop]->trial[ic]++;

    //call proposal function to update source parameters
    (*proposal[nprop]->function)(data, model_x, source_y, proposal[nprop], source_y->params, &chain->r[ic]);

    //hold sky position fixed to injected value
    if(flags->fixSky)
    {
        source_y->params[1] = source_x->params[1];
        source_y->params[2] = source_x->params[2];
    }
    
    //hold frequencies fixed to injected value
    if(flags->fixFreq) source_y->params[0] = source_x->params[0];
    if(flags->fixFdot) source_y->params[7] = source_x->params[7];
    
    //call associated proposal density functions
    logQyx = (*proposal[nprop]->density)(data, model_x, source_y, proposal[nprop], source_y->params);
    logQxy = (*proposal[nprop]->density)(data, model_x, source_x, proposal[nprop], source_x->params);

    map_array_to_params(source_y, source_y->params, data->T);

    //get priors for x and y
    logPx = evaluate_prior(flags, data, model_x, prior, source_x->params);
    logPy = evaluate_prior(flags, data, model_y, prior, source_y->params);

   
    if(logPy > -INFINITY)
    {
        if(!flags->prior)
        {
            //form master template
            if(!strcmp("fourier",data->basis)) generate_signal_model(orbit, data, model_y, n);
            if(!strcmp("wavelet",data->basis)) generate_signal_model_wavelet(orbit, data, model_y, n);  
            //update master template
            //update_signal_model(orbit, data, model_x, model_y, n);
            

            //rejection sample on SNR?
            double sn;
            if(!strcmp("fourier",data->basis)) sn = snr(model_y->source[n], model_y->noise);
            if(!strcmp("wavelet",data->basis)) sn = snr_wavelet(model_y->source[n], model_y->noise);
            if(sn<5.0) logPy = -INFINITY;

            //add calibration error
            if(flags->calibration)
            {
                generate_calibration_model(data, model_y);
                apply_calibration_model(data, model_y);
            }
            

            //get likelihood for y
            if(!strcmp("fourier",data->basis)) model_y->logL = gaussian_log_likelihood(data, model_y);
            if(!strcmp("wavelet",data->basis)) model_y->logL = gaussian_log_likelhood_wavelet(data, model_y);

            //get delta log likelihood
            //model_y->logL = model_x->logL + delta_log_likelihood(data, model_x, model_y, n);
            
            /*
             H = [p(d|y)/p(d|x)]/T x p(y)/p(x) x q(x|y)/q(y|x)
             */
            logH += (model_y->logL - model_x->logL)/chain->temperature[ic]; //delta logL
        }
        logH += logPy  - logPx;  //priors
        logH += logQxy - logQyx; //proposals
        
        loga = log(rand_r_U_0_1(&chain->r[ic]));
        
        if(isfinite(logH) && logH > loga)
        {

            proposal[nprop]->accept[ic]++;
            copy_model_lite(model_y,model_x);

        }
    }
}

static void rj_birth_death(struct Orbit *orbit, struct Data *data, struct Model *model_x, struct Model *model_y, struct Chain *chain, struct Flags *flags, struct Prior *prior, struct Proposal *proposal, int ic, double *logQxy, double *logQyx, double *logPy, double *penalty)
{
    /* pick birth or death move */
    if(rand_r_U_0_1(&chain->r[ic])<0.5)/* birth move */
    {
        //ny=nx+1
        model_y->Nlive++;
        
        //slot new source in at end of live  source array
        int create = model_y->Nlive-1;
        
        if(model_y->Nlive < model_x->Neff)
        {
            //draw new parameters
            //TODO: insert draw from galaxy prior into draw_from_uniform_prior()
            *logQyx = (*proposal->function)(data, model_y, model_y->source[create], proposal, model_y->source[create]->params, &chain->r[ic]);
            *logQxy = 0;
            map_array_to_params(model_y->source[create], model_y->source[create]->params, data->T);
            
            if(flags->maximize)
            {
                maximize_signal_model(orbit, data, model_y, create);
                *penalty = maximization_penalty(4,2*model_y->source[create]->BW);
            }
            
            if(!strcmp("fourier",data->basis)) generate_signal_model(orbit, data, model_y, create);
            if(!strcmp("wavelet",data->basis)) generate_signal_model_wavelet(orbit, data, model_y, create);

            //rejection sample on SNR?
            if(!flags->prior)
            {
                double sn;
                if(!strcmp("fourier",data->basis)) sn = snr(model_y->source[create], model_y->noise);
                if(!strcmp("wavelet",data->basis)) sn = snr_wavelet(model_y->source[create], model_y->noise);
                if(sn<5.0) *logPy = -INFINITY;
            }
                

            model_y->source[create]->fisher_update_flag = 1;
        }
        else *logPy = -INFINITY;
    }
    else /* death move */
    {
        //ny=nx-1
        model_y->Nlive--;
        
        //pick source to kill
        int kill = (int)(rand_r_U_0_1(&chain->r[ic])*(double)model_x->Nlive);
        if(!strcmp("fourier",data->basis)) remove_signal_model(data, model_y,  model_y->source[kill]);
        if(!strcmp("wavelet",data->basis)) remove_signal_model_wavelet(data, model_y, model_y->source[kill]);
        
        if(model_y->Nlive>-1)
        {
            *logQyx = 0;
            *logQxy = (*proposal->density)(data, model_y, model_y->source[kill], proposal, model_y->source[kill]->params);
            
            //consolodiate parameter structure
            for(int j=kill; j<model_y->Nlive; j++)
            {
                copy_source(model_x->source[j+1],model_y->source[j]);
            }
            
            if(flags->maximize)
            {
                *penalty = maximization_penalty(4,2*model_x->source[kill]->BW);
            }
            
            model_y->source[model_y->Nlive]->fisher_update_flag = 1;
        }
        else *logPy = -INFINITY;
    }
}

static void rj_split_merge(struct Orbit *orbit, struct Data *data, struct Model *model_x, struct Model *model_y, struct Chain *chain, struct Flags *flags, struct Prior *prior, struct Proposal *proposal, int ic, double *logQxy, double *logQyx, double *logPy, double *penalty)
{
    /* pick split or merge move */
    int branch[2];
    int trunk;
    
    if(rand_r_U_0_1(&chain->r[ic])<0.5)/* split move */
    {
        //ny=nx+1
        model_y->Nlive++;
        
        //pick source to split
        trunk = (int)(rand_r_U_0_1(&chain->r[ic])*(double)model_x->Nlive);
        
        //slot new sources in place and at end of live source array
        branch[0] = trunk;
        branch[1] = model_y->Nlive-1;
        
        for(int n=0; n<data->N; n++)
        {
            for(int m=0; m<2; m++)
            {
                model_y->source[branch[m]]->tdi->X[n]=0.0;
                model_y->source[branch[m]]->tdi->Y[n]=0.0;
                model_y->source[branch[m]]->tdi->Z[n]=0.0;
                model_y->source[branch[m]]->tdi->A[n]=0.0;
                model_y->source[branch[m]]->tdi->E[n]=0.0;
            }
        }

        
        if(model_y->Nlive<model_x->Neff)
        {
            
            //draw parameters for new sources (branches of trunk)
            for(int n=0; n<2; n++)
            {
                *logQyx += (*proposal->function)(data, model_y, model_y->source[branch[n]], proposal, model_y->source[branch[n]]->params, &chain->r[ic]);
                map_array_to_params(model_y->source[branch[n]], model_y->source[branch[n]]->params, data->T);
                if(flags->maximize)
                {
                    maximize_signal_model(orbit, data, model_y, branch[n]);
                    *penalty += maximization_penalty(4,2*model_y->source[branch[n]]->BW);
                }
                generate_signal_model(orbit, data, model_y, branch[n]);

                //rejection sample on SNR?
                if(snr(model_y->source[branch[n]], model_y->noise) < 5.0) *logPy = -INFINITY;

                model_y->source[branch[n]]->fisher_update_flag = 1;
            }
            
            //get reverse move (merge branches to trunk)
            *logQxy += (*proposal->density)(data, model_x, model_x->source[trunk], proposal, model_x->source[trunk]->params);

            if(flags->maximize)
            {
                *penalty -= maximization_penalty(4,2*model_x->source[trunk]->BW);
            }

        }
        else *logPy = -INFINITY;
    }
    else /* merge move */
    {
        //ny=nx-1
        model_y->Nlive--;
        
        if(model_y->Nlive>0)
        {
            //pick sources to merge
            do
            {
                branch[0] = (int)(rand_r_U_0_1(&chain->r[ic])*(double)model_x->Nlive);
                branch[1] = (int)(rand_r_U_0_1(&chain->r[ic])*(double)model_x->Nlive);
            }while(branch[0]==branch[1]);
            
            
            //pick source to replace with merged sources
            trunk = branch[0];
                        
            //consolodate parameter structure
            for(int j=branch[1]; j<model_y->Nlive; j++)
            {
                copy_source(model_x->source[j+1],model_y->source[j]);
            }
            
            //draw parameters of trunk
            *logQyx += (*proposal->function)(data, model_y, model_y->source[trunk], proposal, model_y->source[trunk]->params, &chain->r[ic]);
            map_array_to_params(model_y->source[trunk], model_y->source[trunk]->params, data->T);
            if(flags->maximize)
            {
                maximize_signal_model(orbit, data, model_y, trunk);
                *penalty -= maximization_penalty(4,2*model_y->source[trunk]->BW);
            }
            generate_signal_model(orbit, data, model_y, trunk);
            
            //rejection sample on SNR?
            if(snr(model_y->source[trunk], model_y->noise) < 5.0) *logPy = -INFINITY;

            model_y->source[trunk]->fisher_update_flag = 1;
            
            //get reverse move (split trunk into branches)
            for(int n=0; n<2; n++)
            {
                *logQxy += (*proposal->density)(data, model_x, model_x->source[branch[n]], proposal, model_x->source[branch[n]]->params);
                
                if(flags->maximize)
                {
                    *penalty += maximization_penalty(4,2*model_x->source[branch[n]]->BW);
                }

            }
        }
        else *logPy = -INFINITY;
        
    }

}

static void rj_cluster_bomb(struct Orbit *orbit, struct Data *data, struct Model *model_x, struct Model *model_y, struct Chain *chain, struct Flags *flags, struct Prior *prior, struct Proposal *proposal, int ic, double *logQxy, double *logQyx, double *logPy, double *penalty)
{
    
    int N = model_x->Nlive;
    double * f = double_vector(N);
    for(int n=0; n<N; n++) f[n] = model_x->source[n]->f0;
    
    //DBSCAN clustering of source frequencies
    int K=0;                //number of clusters
    double eps = model_x->source[0]->BW/data->T; //use typical bandwidth as the max cluster spacing
    //double eps = 2/data->T; //use typical bandwidth as the max cluster spacing
    int min = 2;          //minimum occupation number for cluster
    int C[N];             //cluster assignments
    if(N>=min) dbscan(f,eps,min,C,&K,N);
    
    model_y->Nlive=model_x->Nlive;
    for(int n=0; n<N; n++)
        copy_source(model_x->source[n],model_y->source[n]);

    if(K>0) // if there are clusters to try killing
    {
        //pick a cluster to kill
        int ckill = (int)(rand_r_U_0_1(&chain->r[ic]) * (double)K);
        double qmax = 0.0;
        double qmin = data->T;
        
        //consolodate parameter structure
        model_y->Nlive = 0;
        for(int n=0; n<N; n++)
        {
            if(C[n]!=ckill)
            {
                copy_source(model_x->source[n],model_y->source[model_y->Nlive]);
                model_y->Nlive++;
            }
            else
            {
                if(model_x->source[n]->params[0] < qmin) qmin = model_x->source[n]->params[0];
                if(model_x->source[n]->params[0] > qmax) qmax = model_x->source[n]->params[0];
            }
        }
        //restrict prior to cluster width
        model_y->prior[0][0] = qmin-0.25;
        model_y->prior[0][1] = qmax+0.25;

        
        //draw parameters and generate signal model for replacement
        *logQyx += (*proposal->function)(data, model_y, model_y->source[model_y->Nlive], proposal, model_y->source[model_y->Nlive]->params, &chain->r[ic]);
        

        if(flags->maximize)
        {
            maximize_signal_model(orbit, data, model_y, model_y->Nlive);
            *penalty -= maximization_penalty(4,2*model_y->source[model_y->Nlive]->BW);
        }
        model_y->Nlive++;

        generate_signal_model(orbit, data, model_y, -1);
        
        //rejection sample on SNR?
        if(snr(model_y->source[model_y->Nlive], model_y->noise) < 5.0) *logPy = -INFINITY;
        
        for(int n=0; n<model_y->Nlive; n++)
            model_y->source[n]->fisher_update_flag = 1;

        //get reverse move (propose every source in cluster)
        for(int n=0; n<N; n++)
        {
            if(C[n]==ckill)
            {
                *logQxy += (*proposal->density)(data, model_x, model_x->source[n], proposal, model_x->source[n]->params);
                if(flags->maximize)
                    *penalty += maximization_penalty(4,2*model_x->source[n]->BW);
            }
        }

        //restore full prior range
        model_y->prior[0][0] = model_x->prior[0][0];
        model_y->prior[0][1] = model_x->prior[0][1];

    }
    else *logPy = -INFINITY;
    
    free_double_vector(f);
}

void ucb_rjmcmc(struct Orbit *orbit, struct Data *data, struct Model *model, struct Model *trial, struct Chain *chain, struct Flags *flags, struct Prior *prior, struct Proposal **proposal, int ic)
{
    double logH  = 0.0; //(log) Hastings ratio
    double loga  = 1.0; //(log) transition probability
    
    double logPx  = 0.0; //(log) prior density for model x (current state)
    double logPy  = 0.0; //(log) prior density for model y (proposed state)
    double logQyx = 0.0; //(log) proposal denstiy from x->y
    double logQxy = 0.0; //(log) proposal density from y->x
    
    double dlogL  = 0.0; //delta log likelihood
    
    /*
     * BIC-inspired likelihood penalty -klogN/2
     * k = [A, inc, psi, phi]
     * N is the bandwidth of the signal
     */
    double penalty = 0.0;

    //shorthand pointers
    struct Model *model_x = model;
    struct Model *model_y = trial;
    
    copy_model(model_x,model_y);
    
    int nprop;
    do nprop = (int)floor((UCB_PROPOSAL_NPROP)*rand_r_U_0_1(&chain->r[ic]));
    while(proposal[nprop]->rjweight <= rand_r_U_0_1(&chain->r[ic]));
    
    proposal[nprop]->trial[ic]++;
        
    /* Choose birth/death move, or split/merge move */
    if( rand_r_U_0_1(&chain->r[ic]) < 1.5)/* birth/death move */
        rj_birth_death(orbit, data, model_x, model_y, chain, flags, prior, proposal[nprop], ic, &logQxy, &logQyx, &logPy, &penalty);
    else if( rand_r_U_0_1(&chain->r[ic]) < 1.5) /* birth/death move */
        rj_split_merge(orbit, data, model_x, model_y, chain, flags, prior, proposal[nprop], ic, &logQxy, &logQyx, &logPy, &penalty);
    else
        rj_cluster_bomb(orbit, data, model_x, model_y, chain, flags, prior, proposal[nprop], ic, &logQxy, &logQyx, &logPy, &penalty);


    if(logPy > -INFINITY )
    {
        for(int n=0; n<model_x->Nlive; n++) logPx += evaluate_prior(flags, data, model_x, prior, model_x->source[n]->params);
        for(int n=0; n<model_y->Nlive; n++) logPy += evaluate_prior(flags, data, model_y, prior, model_y->source[n]->params);
    }

    /* Hasting's ratio */
    if(logPy > -INFINITY && !flags->prior)
    {
        //  Form master template
        /*
         generate_signal_model() is called inside rj proposal wrappers
         rj_birth_death() and rj_split_merge()
         */
        
        //calibration error
        if(flags->calibration)
        {
            generate_calibration_model(data, model_y);
            apply_calibration_model(data, model_y);
        }
        
        //get likelihood for y
        if(!strcmp("fourier",data->basis)) model_y->logL = gaussian_log_likelihood(data, model_y);
        if(!strcmp("wavelet",data->basis)) model_y->logL = gaussian_log_likelhood_wavelet(data, model_y);

        //get likelihood difference
        dlogL = model_y->logL - model_x->logL;
        
        //penalize likelihood when using maximized parameters
        if(flags->maximize) dlogL += (model_y->Nlive - model_x->Nlive)*penalty;

        /*
         H = [p(d|y)/p(d|x)]/T x p(y)/p(x) x q(x|y)/q(y|x)
         */
        logH += dlogL/chain->temperature[ic]; //delta logL
    }
    
    
    logH += logPy  - logPx;  //priors
    logH += logQxy - logQyx; //proposals
    
////      if(model_y->Nlive > model_x->Nlive && ic==0 && sm_flag)
////        if(ic==0)
//        if(sm_flag && ic==0)
//        {
////          FILE *fptr = fopen("proposal.dat","a");
//              fprintf(stdout,"%lg %lg %lg %lg %g ",model_y->logL+model_y->logLnorm, model_x->logL+model_x->logLnorm, logH, logPy  - logPx, logQxy - logQyx);
//              fprintf(stdout,"%i -> %i \n",model_x->Nlive, model_y->Nlive);
////          print_source_params(data, model_y->source[model_y->Nlive-1], fptr);
////          fprintf(fptr,"\n");
////          fclose(fptr);
//        }
    
    loga = log(rand_r_U_0_1(&chain->r[ic]));
    if(isfinite(logH) && logH > loga)
    {
        //update FIMs for sources that have changed in the proposed state
        for(int n=0; n<model_y->Nlive; n++)
        {
            //only compute FIM for sources that need it
            if(model_y->source[n]->fisher_update_flag)
            {
                ucb_fisher(orbit, data, model_y->source[n], data->noise);

                //reset flag indicating FIM is up-to-date
                model_y->source[n]->fisher_update_flag = 0;
            }
        }
        proposal[nprop]->accept[ic]++;
        copy_model(model_y,model_x);
    }
}

void initialize_ucb_state(struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Chain *chain, struct Proposal **proposal, struct Model **model, struct Model **trial, struct Source **inj_vec)
{
    int NC = chain->NC;
    int DMAX = flags->DMAX;
    
    //keep track of number of injections
    struct Source *inj;
    int Ninj=0;
    
    for(int ic=0; ic<NC; ic++)
    {
        trial[ic] = malloc(sizeof(struct Model));
        model[ic] = malloc(sizeof(struct Model));
        
        alloc_model(data,trial[ic],DMAX);
        alloc_model(data,model[ic],DMAX);

        if(ic==0)set_uniform_prior(flags, model[ic], data, 1);
        else     set_uniform_prior(flags, model[ic], data, 0);
                
        //override noise model w/ stationary version
        if(!strcmp("wavelet",data->basis) && flags->stationary) GetStationaryNoiseModel(data, orbit, flags, data->noise);

        //set noise model
        copy_noise(data->noise, model[ic]->noise);
        
        //draw signal model
        for(int n=0; n<DMAX; n++)
        {
            if(flags->cheat && inj_vec[n]->f0>0)
            {
                
                if(ic==0)Ninj++;
                inj=inj_vec[n];
                
                //map parameters to vector
                model[ic]->source[n]->f0       = inj->f0;
                model[ic]->source[n]->dfdt     = inj->dfdt;
                model[ic]->source[n]->costheta = inj->costheta;
                model[ic]->source[n]->phi      = inj->phi;
                model[ic]->source[n]->amp      = inj->amp;
                model[ic]->source[n]->cosi     = inj->cosi;
                model[ic]->source[n]->phi0     = inj->phi0;
                model[ic]->source[n]->psi      = inj->psi;
                model[ic]->source[n]->d2fdt2   = inj->d2fdt2;
                map_params_to_array(model[ic]->source[n], model[ic]->source[n]->params, data->T);
                if(!strcmp("wavelet",data->basis) && ic==0 && n==0 && flags->stationary)
                    fprintf(stdout,"   ...stationary SNR for source %i=%g\n",n,snr_wavelet(inj,model[ic]->noise));
            }
            else
            {
                draw_from_uniform_prior(data, model[ic], model[ic]->source[n], proposal[0], model[ic]->source[n]->params, &chain->r[ic]);
            }
            map_array_to_params(model[ic]->source[n], model[ic]->source[n]->params, data->T);

            if(!strcmp("fourier",data->basis)) ucb_fisher(orbit, data, model[ic]->source[n], data->noise);
            if(!strcmp("wavelet",data->basis)) ucb_fisher_wavelet(orbit, data, model[ic]->source[n], data->noise);

            model[ic]->source[n]->fisher_update_flag=0;
        }
        
        //initialize sampler to proper size of model
        if(flags->cheat) 
        {
            model[ic]->Neff = Ninj+1;
            model[ic]->Nlive= Ninj;
        }
        
        // Form master model & compute likelihood of starting position
        if(!strcmp("fourier",data->basis)) 
        {
            generate_noise_model(data, model[ic]);
            generate_signal_model(orbit, data, model[ic], -1);
        }
        if(!strcmp("wavelet",data->basis)) 
        {
            generate_noise_model_wavelet(data, model[ic]);
            generate_signal_model_wavelet(orbit, data, model[ic], -1);
        }

        //calibration error
        if(flags->calibration)
        {
            draw_calibration_parameters(data, model[ic], &chain->r[ic]);
            generate_calibration_model(data, model[ic]);
            apply_calibration_model(data, model[ic]);
        }
        if(!flags->prior)
        {
            if(!strcmp("fourier",data->basis))model[ic]->logL = gaussian_log_likelihood(data, model[ic]);
            if(!strcmp("wavelet",data->basis))model[ic]->logL = gaussian_log_likelhood_wavelet(data, model[ic]);
            model[ic]->logLnorm = gaussian_log_likelihood_model_norm(data,model[ic]);
        
        }
        else model[ic]->logL = model[ic]->logLnorm = 0.0;
        
        if(ic==0) chain->logLmax += model[ic]->logL + model[ic]->logLnorm;
        
    }//end loop over chains

}


