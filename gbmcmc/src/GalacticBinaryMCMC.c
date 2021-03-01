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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <omp.h>

#include <LISA.h>

#include "GalacticBinary.h"
#include "GalacticBinaryIO.h"
#include "GalacticBinaryData.h"
#include "GalacticBinaryPrior.h"
#include "GalacticBinaryModel.h"
#include "GalacticBinaryProposal.h"
#include "GalacticBinaryWaveform.h"
#include "GalacticBinaryMCMC.h"

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
    for(int i=0; i<flags->NT; i++)
    {
        switch(data->Nchannel)
        {
            case 1:
                model_y->noise[i]->etaX = model_x->noise[i]->etaX + 0.1*gsl_ran_gaussian(chain->r[ic],1);
                break;
            case 2:
                model_y->noise[i]->etaA = model_x->noise[i]->etaA + 0.1*gsl_ran_gaussian(chain->r[ic],1);
                model_y->noise[i]->etaE = model_x->noise[i]->etaE + 0.1*gsl_ran_gaussian(chain->r[ic],1);
                break;
        }
        
        //get priors for x and y
        switch(data->Nchannel)
        {
            case 1:
                if(model_y->noise[i]->etaX < 0.01 || model_y->noise[i]->etaX>100) logPy=-INFINITY;
                break;
            case 2:
                if(model_y->noise[i]->etaA < 0.01 || model_y->noise[i]->etaA>100.) logPy=-INFINITY;
                if(model_y->noise[i]->etaE < 0.01 || model_y->noise[i]->etaE>100.) logPy=-INFINITY;
                break;
        }
    }
    
    
    if(!flags->prior)
    {
        //  Form master template
        generate_noise_model(data, model_y);
        
        //get likelihood for y
        model_y->logL     = gaussian_log_likelihood(data, model_y);
        model_y->logLnorm = gaussian_log_likelihood_constant_norm(data, model_y);
        
        /*
         H = [p(d|y)/p(d|x)]/T x p(y)/p(x) x q(x|y)/q(y|x)
         */
        logH += ( (model_y->logL+model_y->logLnorm) - (model_x->logL+model_x->logLnorm) )/chain->temperature[ic]; //delta logL
    }
    logH += logPy  - logPx;                                         //priors
    
    loga = log(gsl_rng_uniform(chain->r[ic]));
    if(logH > loga) copy_model(model_y,model_x);
    
}

void galactic_binary_mcmc(struct Orbit *orbit, struct Data *data, struct Model *model, struct Model *trial, struct Chain *chain, struct Flags *flags, struct Prior *prior, struct Proposal **proposal, int ic)
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
    int n = (int)(gsl_rng_uniform(chain->r[ic])*(double)model_x->Nlive);
    
    //more shorthand pointers
    struct Source *source_x = model_x->source[n];
    struct Source *source_y = model_y->source[n];
    
    
    //choose proposal distribution
    int trial_n;
    double trial_w;
    int nprop=-1;
    
    while(nprop<0)
    {
        trial_n = (int)floor((chain->NP)*gsl_rng_uniform(chain->r[ic]));
        trial_w = gsl_rng_uniform(chain->r[ic]);
        if(trial_w < proposal[trial_n]->weight) nprop = trial_n;
    }
    proposal[nprop]->trial[ic]++;
    
    //call proposal function to update source parameters
    (*proposal[nprop]->function)(data, model_x, source_y, proposal[nprop], source_y->params, chain->r[ic]);
    
    //hold sky position fixed to injected value
    if(flags->fixSky)
    {
        source_y->params[1] = data->inj->costheta;
        source_y->params[2] = data->inj->phi;
    }
    
    //hold frequencies fixed to injected value
    if(flags->fixFreq) source_y->params[0] = data->inj->f0*data->T;
    if(flags->fixFdot) source_y->params[7] = data->inj->dfdt*data->T*data->T;
    
    //call associated proposal density functions
    logQyx = (*proposal[nprop]->density)(data, model_x, source_y, proposal[nprop], source_y->params);
    logQxy = (*proposal[nprop]->density)(data, model_x, source_x, proposal[nprop], source_x->params);
        
    map_array_to_params(source_y, source_y->params, data->T);

    /*
     if(flags->maximize &&
       !check_range(source_y->params, model->prior, model->NP) )
        maximize_signal_model(orbit, data, model_y, n);
     */

    
    //update calibration parameters
    if(flags->calibration) draw_calibration_parameters(data, model_y, chain->r[ic]);
    /*
     no proposal density for calibration parameters
     because we are always drawing from prior...for now
     */
    
    //TODO:copy params for segment 0 into higher segments
    //copy_source(model_y->source[n],model_y->source[n]);
    //map_params_to_array(model_y->source[n], model_y->source[n]->params, data->T);
    
    //get priors for x and y
    logPx = evaluate_prior(flags, data, model_x, prior, source_x->params);
    logPy = evaluate_prior(flags, data, model_y, prior, source_y->params);
    
    //add calibration source parameters
    /*
     no prior density for calibration parameters
     because we are always drawing from prior...for now
     */
    
    if(logPy > -INFINITY)
    {
        if(!flags->prior)
        {
            //  Form master template
            generate_signal_model(orbit, data, model_y, n);
            
            //calibration error
            if(flags->calibration)
            {
                generate_calibration_model(data, model_y);
                apply_calibration_model(data, model_y);
            }
            
            //get likelihood for y
            model_y->logL = gaussian_log_likelihood(data, model_y);
            
            /*
             H = [p(d|y)/p(d|x)]/T x p(y)/p(x) x q(x|y)/q(y|x)
             */
            logH += (model_y->logL - model_x->logL)/chain->temperature[ic]; //delta logL
        }
        logH += logPy  - logPx;  //priors
        logH += logQxy - logQyx; //proposals
        
        loga = log(gsl_rng_uniform(chain->r[ic]));
        
        if(isfinite(logH) && logH > loga)
        {
            proposal[nprop]->accept[ic]++;
            copy_model(model_y,model_x);
        }
    }
}

void galactic_binary_rjmcmc(struct Orbit *orbit, struct Data *data, struct Model *model, struct Model *trial, struct Chain *chain, struct Flags *flags, struct Prior *prior, struct Proposal **proposal, int ic)
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
     */
    double penalty= -2.;//1.*data->logN;

    //shorthand pointers
    struct Model *model_x = model;
    struct Model *model_y = trial;
    
    copy_model(model_x,model_y);
    
    int nprop = -1;
    int trial_n;
    double trial_w;
    while(nprop<0)
    {
        trial_n = (int)floor((chain->NP)*gsl_rng_uniform(chain->r[ic]));
        trial_w = gsl_rng_uniform(chain->r[ic]);
        if(trial_w < proposal[trial_n]->rjweight) nprop = trial_n;
    }
    
    proposal[nprop]->trial[ic]++;
    
    /* pick birth or death move */
    if(gsl_rng_uniform(chain->r[ic])<0.5)/* birth move */
    {
        //ny=nx+1
        model_y->Nlive++;
        
        //slot new source in at end of live  source array
        int create = model_y->Nlive-1;
        
        if(model_y->Nlive<model_x->Nmax)
        {
            //draw new parameters
            //TODO: insert draw from galaxy prior into draw_from_uniform_prior()
            logQyx = (*proposal[nprop]->function)(data, model_y, model_y->source[create], proposal[nprop], model_y->source[create]->params, chain->r[ic]);
            logQxy = 0;
            
            map_array_to_params(model_y->source[create], model_y->source[create]->params, data->T);
            
            if(flags->maximize) maximize_signal_model(orbit, data, model_y, create);

        }
        else logPy = -INFINITY;
    }
    else /* death move */
    {
        //ny=nx-1
        model_y->Nlive--;
        
        //pick source to kill
        int kill = (int)(gsl_rng_uniform(chain->r[ic])*(double)model_x->Nlive);
        
        if(model_y->Nlive>-1)
        {
            logQyx = 0;
            logQxy = (*proposal[nprop]->density)(data, model_y, model_y->source[kill], proposal[nprop], model_y->source[kill]->params);
            
            //consolodiate parameter structure
            for(int j=kill; j<model_x->Nlive; j++)
            {
                copy_source(model_x->source[j+1],model_y->source[j]);
            }
        }
        else logPy = -INFINITY;
    }
    
    for(int n=0; n<model_x->Nlive; n++) logPx +=  evaluate_prior(flags, data, model_x, prior, model_x->source[n]->params);
    for(int n=0; n<model_y->Nlive; n++) logPy +=  evaluate_prior(flags, data, model_y, prior, model_y->source[n]->params);
    
    
    /* Hasting's ratio */
    if(logPy > -INFINITY && !flags->prior)
    {
        //  Form master template
        /*
         generate_signal_model is passed an integer telling it which source to update.
         passing model_x->Nlive is a trick to skip waveform generation for kill move
         and to only calculate new source for create move
         */
        generate_signal_model(orbit, data, model_y, model_x->Nlive);
        
        //calibration error
        if(flags->calibration)
        {
            generate_calibration_model(data, model_y);
            apply_calibration_model(data, model_y);
        }
        
        //get likelihood for y
        model_y->logL = gaussian_log_likelihood(data, model_y);
        
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
    
    //  if(model_y->Nlive > model_x->Nlive && ic==0)
    //    if(ic==0)
    //    {
    //      FILE *fptr = fopen("proposal.dat","a");
    //          fprintf(stdout,"%lg %lg %lg %lg %g ",model_y->logL+model_y->logLnorm, model_x->logL+model_x->logLnorm, logH, logPy  - logPx, logQxy - logQyx);
    //          fprintf(stdout,"%i -> %i \n",model_x->Nlive, model_y->Nlive);
    //      print_source_params(data, model_y->source[model_y->Nlive-1], fptr);
    //      fprintf(fptr,"\n");
    //      fclose(fptr);
    //    }
    
    loga = log(gsl_rng_uniform(chain->r[ic]));
    if(isfinite(logH) && logH > loga)
    {
        proposal[nprop]->accept[ic]++;
        copy_model(model_y,model_x);
    }
    
}

void data_mcmc(struct Orbit *orbit, struct Data *data, struct Model *model, struct Chain *chain, struct Flags *flags, struct Proposal **proposal, int ic)
{
    double logH  = 0.0; //(log) Hastings ratio
    double loga  = 1.0; //(log) transition probability
    double logQ  = 0.0;
    
    struct Model *trial = malloc(sizeof(struct Model));
    
    
    alloc_model(trial,model->Nmax,data->N,data->Nchannel, data->NP, flags->NT);
    
    set_uniform_prior(flags, trial, data, 0);
    
    copy_model(model,trial);
    
    logQ += t0_shift(data, trial, trial->source[0], proposal[0], trial->source[0]->params, chain->r[ic]);
    
    // Form master template
    /*
     passing generate_signal_model -1 results in full recalculation of waveform model
     */
    generate_signal_model(orbit, data, trial, -1);
    
    /*
     H = [p(d|y)/p(d|x)]/T x p(y)/p(x) x q(x|y)/q(y|x)
     */
    if(!flags->prior)
    {
        // get likelihood for y
        trial->logL = gaussian_log_likelihood(data, trial);
        
        logH += (trial->logL - model->logL)/chain->temperature[ic];
    }
    logH += logQ; //delta logL
    
    
    loga = log(gsl_rng_uniform(chain->r[ic]));
    
    if(logH > loga) copy_model(trial,model);
    
    free_model(trial);
    
}

void initialize_gbmcmc_state(struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Chain *chain, struct Proposal **proposal, struct Model **model, struct Model **trial)
{
    int NC = chain->NC;
    int DMAX = flags->DMAX;
    for(int ic=0; ic<NC; ic++)
    {
        
        trial[ic] = malloc(sizeof(struct Model));
        model[ic] = malloc(sizeof(struct Model));
        
        alloc_model(trial[ic],DMAX,data->N,data->Nchannel,data->NP, data->NT);
        alloc_model(model[ic],DMAX,data->N,data->Nchannel, data->NP, flags->NT);
                
        if(ic==0)set_uniform_prior(flags, model[ic], data, 1);
        else     set_uniform_prior(flags, model[ic], data, 0);
        
        //set noise model
        for(int j=0; j<flags->NT; j++) copy_noise(data->noise[j], model[ic]->noise[j]);
        
        //draw signal model
        for(int n=0; n<DMAX; n++)
        {
            if(flags->cheat)
            {
                struct Source *inj = data->inj;
                //map parameters to vector
                model[ic]->source[n]->NP       = inj->NP;
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
                
            }
            else if(flags->updateCov)
            {
                while ( !isfinite(draw_from_cov(data, model[ic], model[ic]->source[n], proposal[8], model[ic]->source[n]->params , chain->r[ic])));
            }
            else if(flags->update)
            {
                draw_from_gmm_prior(data, model[ic], model[ic]->source[n], proposal[7], model[ic]->source[n]->params , chain->r[ic]);
            }
            else
            {
                draw_from_uniform_prior(data, model[ic], model[ic]->source[n], proposal[0], model[ic]->source[n]->params , chain->r[ic]);
            }
            map_array_to_params(model[ic]->source[n], model[ic]->source[n]->params, data->T);
            galactic_binary_fisher(orbit, data, model[ic]->source[n], data->noise[0]);
        }
        
        // Form master model & compute likelihood of starting position
        generate_noise_model(data, model[ic]);
        generate_signal_model(orbit, data, model[ic], -1);
        
        //calibration error
        if(flags->calibration)
        {
            draw_calibration_parameters(data, model[ic], chain->r[ic]);
            generate_calibration_model(data, model[ic]);
            apply_calibration_model(data, model[ic]);
        }
        if(!flags->prior)
        {
            model[ic]->logL     = gaussian_log_likelihood(data, model[ic]);
            model[ic]->logLnorm = gaussian_log_likelihood_constant_norm(data, model[ic]);
        }
        else model[ic]->logL = model[ic]->logLnorm = 0.0;
        
        if(ic==0) chain->logLmax += model[ic]->logL + model[ic]->logLnorm;
        
    }//end loop over chains

}


