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

/**
 @file GalacticBinaryMCMC.c
 \brief Sampling routines and main function
 
 Including
  - Fixed dimension galactic binary MCMC
  - Trans-dimension galactic binary RJMCMC
  - Parallel tempering chain exchanges ptmcmc()
  - Fixed dimension noise model MCMC
  - Fixed dimension data model MCMC
 */

/***************************  REQUIRED LIBRARIES  ***************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//#include <omp.h>

/*************  PROTOTYPE DECLARATIONS FOR INTERNAL FUNCTIONS  **************/

#include "LISA.h"
#include "Constants.h"
#include "BayesLine.h"
#include "GalacticBinary.h"
#include "GalacticBinaryIO.h"
#include "GalacticBinaryData.h"
#include "GalacticBinaryPrior.h"
#include "GalacticBinaryModel.h"
#include "GalacticBinaryProposal.h"
#include "GalacticBinaryWaveform.h"

/** \brief Parallel tempering exchange
 
 Cycles through all chains and proposes swaps between adjacent pairs.
*/
void ptmcmc(struct Model ***model, struct Chain *chain, struct Flags *flags);

/**
 \brief Adaptive temperature spacing
 
 Adjusts temperature spacing between tempered chains according with the goal of having an even acceptance rate between all pairs. The sensitivity of the adjustment asymptotically goes to zero as the sampler approaches the end of the burn-in phase.
*/
void adapt_temperature_ladder(struct Chain *chain, int mcmc);

/**
\brief Fixed dimension galactic binary MCMC
 
 One fixed-dimension MCMC step for galactic binary model.
 The sampler chooses a source at random from the full model to update.
 @param[in] ic index for which chain is being updated
*/
void galactic_binary_mcmc(struct Orbit *orbit, struct Data *data, struct Model *model, struct Model *trial, struct Chain *chain, struct Flags *flags, struct Prior *prior, struct Proposal **proposal, int ic);

/**
\brief Trans-dimension galactic binary RJMCMC
 
 One trans-dimension RJMCMC step for galactic binary model.
 The sampler chooses to either add or remove a galactic binary signal to the model.
 @param[in] ic index for which chain is being updated
*/
void galactic_binary_rjmcmc(struct Orbit *orbit, struct Data *data, struct Model *model, struct Model *trial, struct Chain *chain, struct Flags *flags, struct Prior *prior, struct Proposal **proposal, int ic);

/**
\brief Data model MCMC
 
 One fixed-dimension MCMC step to update the data model.
 Currently this includes calibration parameters but could be extended to TDI, etc.
 @param[in] ic index for which chain is being updated
*/
void data_mcmc(struct Orbit *orbit, struct Data **data, struct Model **model, struct Chain *chain, struct Flags *flags, struct Proposal **proposal, int ic);

/**
\brief Noise model MCMC
 
 One fixed-dimension MCMC step to update the noise model.
 Currently this only the variance of the orthogonal A and E TDI channels. This will be generalized to include more complete covariance matrices, including non-stationary and non-orthogonal noise.
 @param[in] ic index for which chain is being updated
*/
void noise_model_mcmc(struct Orbit *orbit, struct Data *data, struct Model *model, struct Model *trial, struct Chain *chain, struct Flags *flags, int ic);

/* ============================  MAIN PROGRAM  ============================ */

/**
 * This is the main function
 *
*/
int main(int argc, char *argv[])
{
    
    time_t start, stop;
    start = time(NULL);
    
    int NMAX = 10;   //max number of frequency & time segments
    
    /* Allocate data structures */
    struct Flags *flags = malloc(sizeof(struct Flags));
    struct Orbit *orbit = malloc(sizeof(struct Orbit));
    struct Chain *chain = malloc(sizeof(struct Chain));
    struct Data  **data = malloc(sizeof(struct Data*)*NMAX); //data[NF]
    
    
    /* Parse command line and set defaults/flags */
    for(int i=0; i<NMAX; i++)
    {
        data[i] = malloc(sizeof(struct Data));
        data[i]->t0   = calloc( NMAX , sizeof(double) );
        data[i]->tgap = calloc( NMAX , sizeof(double) );
    }
    parse(argc,argv,data,orbit,flags,chain,NMAX);
    int NC = chain->NC;
    int DMAX = flags->DMAX;
    int mcmc_start = -flags->NBURN;
    
    /* Allocate model structures */
    struct Model **trial = malloc(sizeof(struct Model*)*NC);//trial[chain]
    struct Model ***model= malloc(sizeof(struct Model**)*NC); //model[chain][source][segment]
    
    /* Load spacecraft ephemerides */
    switch(flags->orbit)
    {
        case 0:
            initialize_analytic_orbit(orbit);
            break;
        case 1:
            initialize_numeric_orbit(orbit);
            break;
        default:
            fprintf(stderr,"unsupported orbit type\n");
            return(1);
            break;
    }
    
    /* Initialize data structures */
    alloc_data(data, flags);
    
    /* Inject strain data */
    if(flags->strainData)
    {
        GalacticBinaryReadData(data,orbit,flags);
    }
    else
    {
        /* Inject gravitational wave signal */
        if(flags->knownSource)
            GalacticBinaryInjectVerificationSource(data,orbit,flags);
        else
            GalacticBinaryInjectSimulatedSource(data,orbit,flags);
    }
    
    /* set approximate f/fstar for segment */
    for(int i=0; i<flags->NDATA; i++)
        data[i]->sine_f_on_fstar = sin((data[i]->fmin + (data[i]->fmax-data[i]->fmin)/2.)/orbit->fstar);
    
    /* Remove sources outside of window */
    if(flags->catalog)
        GalacticBinaryCleanEdges(data, orbit, flags);
    
    /* Initialize data-dependent proposal */
    setup_frequency_proposal(data[0]);
    
    /* Initialize parallel chain */
    if(flags->resume)
        initialize_chain(chain, flags, &data[0]->cseed, "a");
    else
        initialize_chain(chain, flags, &data[0]->cseed, "w");
    
    /* Initialize priors */
    struct Prior *prior = malloc(sizeof(struct Prior));
    if(flags->galaxyPrior) set_galaxy_prior(flags, prior);
    
    /* Initialize MCMC proposals */
    struct Proposal ***proposal = malloc(NMAX*sizeof(struct Proposal**));
    for(int j=0; j<NMAX; j++)
    {
        proposal[j] = malloc((chain->NP)*sizeof(struct Proposal*));
        for(int i=0; i<chain->NP; i++) proposal[j][i] = malloc(sizeof(struct Proposal));
    }
    for(int j=0; j<flags->NDATA; j++) initialize_proposal(orbit, data[j], prior, chain, flags, proposal[j], DMAX);
    
    /* Test noise model */
    //test_noise_model(orbit);
    
    /* Initialize data models */
    for(int ic=0; ic<NC; ic++)
    {
        //printf("initialize model\n");
        
        trial[ic] = malloc(sizeof(struct Model));
        alloc_model(trial[ic],DMAX,data[0]->N,data[0]->Nchannel,data[0]->NP, data[0]->NT);
        
        model[ic] = malloc(sizeof(struct Model *) * flags->NDATA);
        
        //loop over frequency segments
        for(int i=0; i<flags->NDATA; i++)
        {
            //printf("frequency segment %i\n",i);
            
            model[ic][i] = malloc(sizeof(struct Model));
            
            struct Model *model_ptr = model[ic][i];
            struct Data  *data_ptr  = data[i];
            
            alloc_model(model_ptr,DMAX,data_ptr->N,data_ptr->Nchannel, data_ptr->NP, flags->NT);
            
            if(ic==0)set_uniform_prior(flags, model_ptr, data_ptr, 1);
            else     set_uniform_prior(flags, model_ptr, data_ptr, 0);
            
            //set noise model
            for(int j=0; j<flags->NT; j++) copy_noise(data_ptr->noise[j], model_ptr->noise[j]);
            
            //draw signal model
            for(int n=0; n<DMAX; n++)
            {
                if(flags->cheat)
                {
                    struct Source *inj = data_ptr->inj;
                    //map parameters to vector
                    model_ptr->source[n]->NP       = inj->NP;
                    model_ptr->source[n]->f0       = inj->f0;
                    model_ptr->source[n]->dfdt     = inj->dfdt;
                    model_ptr->source[n]->costheta = inj->costheta;
                    model_ptr->source[n]->phi      = inj->phi;
                    model_ptr->source[n]->amp      = inj->amp;
                    model_ptr->source[n]->cosi     = inj->cosi;
                    model_ptr->source[n]->phi0     = inj->phi0;
                    model_ptr->source[n]->psi      = inj->psi;
                    model_ptr->source[n]->d2fdt2   = inj->d2fdt2;
                    map_params_to_array(model_ptr->source[n], model_ptr->source[n]->params, data_ptr->T);
                    
                }
                else if(flags->updateCov)
                {
                    while ( !isfinite(draw_from_cov(data_ptr, model_ptr, model_ptr->source[n], proposal[i][8], model_ptr->source[n]->params , chain->r[ic])));
                }
                else if(flags->update)
                {
                    draw_from_cdf(data_ptr, model_ptr, model_ptr->source[n], proposal[i][7], model_ptr->source[n]->params , chain->r[ic]);
                }
                else
                {
                    draw_from_prior(data_ptr, model_ptr, model_ptr->source[n], proposal[i][0], model_ptr->source[n]->params , chain->r[ic]);
                }
                map_array_to_params(model_ptr->source[n], model_ptr->source[n]->params, data_ptr->T);
                galactic_binary_fisher(orbit, data_ptr, model_ptr->source[n], data_ptr->noise[0]);
            }
            
            // Form master model & compute likelihood of starting position
            generate_noise_model(data_ptr, model_ptr);
            generate_signal_model(orbit, data_ptr, model_ptr, -1);
            
            //calibration error
            if(flags->calibration)
            {
                draw_calibration_parameters(data_ptr, model_ptr, chain->r[ic]);
                generate_calibration_model(data_ptr, model_ptr);
                apply_calibration_model(data_ptr, model_ptr);
            }
            if(!flags->prior)
            {
                model_ptr->logL     = gaussian_log_likelihood(orbit, data_ptr, model_ptr);
                model_ptr->logLnorm = gaussian_log_likelihood_constant_norm(data_ptr, model_ptr);
            }
            else model_ptr->logL = model_ptr->logLnorm = 0.0;
            
            if(ic==0) chain->logLmax += model_ptr->logL + model_ptr->logLnorm;
            
        }//end loop over frequency segments
    }//end loop over chains
    
    /* Start analysis from saved chain state */
    if(flags->resume)
    {
        fprintf(stdout,"\n=============== Checkpointing ===============\n");

        //check for files needed to resume
        FILE *fptr = NULL;
        char filename[MAXSTRINGSIZE];
        int file_error = 0;
        
        for(int ic=0; ic<chain->NC; ic++)
        {
            sprintf(filename,"checkpoint/chain_state_%i.dat",ic);
            
            if( (fptr = fopen(filename,"r")) == NULL )
            {
                fprintf(stderr,"Warning: Could not checkpoint run state\n");
                fprintf(stderr,"         Parameter file %s does not exist\n",filename);
                file_error++;
                break;
            }
        }
        
        //if all of the files exist resume run from checkpointed state
        if(!file_error)
        {
            fprintf(stdout,"   Checkpoint files found. Resuming chain\n");
            restore_chain_state(orbit, data, model, chain, flags, &mcmc_start);
        }
        fprintf(stdout,"============================================\n\n");
    }
    
    //test covariance proposal
    /*
     FILE *test=fopen("proposal_test.dat","w");
     for(int i=0; i<100000; i++)
     {
     double logP = draw_from_cov(data[0], model[0][0], model[0][0]->source[0], proposal[0][8], model[0][0]->source[0]->params, chain->r[0]);
     print_source_params(data[0], model[0][0]->source[0], test);
     fprintf(test,"%lg\n",logP);
     }
     fclose(test);
     */
    //test covariance proposal
    if(flags->updateCov) test_covariance_proposal(data[0], flags, model[0][0], prior, proposal[0][8], chain->r[0]);

    
    /* Write example gb_catalog bash script in run directory */
    print_gb_catalog_script(flags, data[0], orbit);

    
    /* The MCMC loop */
    for(int mcmc = mcmc_start; mcmc < flags->NMCMC; mcmc++)
    {
        if(mcmc<0) flags->burnin=1;
        else       flags->burnin=0;
        
        //set annealinging tempurature during burnin
        /*
         if(flags->burnin)
         {
         chain->annealing = data[0]->SNR2*pow(data[0]->SNR2,-((double)mcmc+(double)flags->NBURN)/((double)flags->NBURN/(double)10))/40.;
         if(chain->annealing<1.0)chain->annealing=1.0;
         chain->annealing=1.0;
         //printf("annealing=%g\n",chain->annealing);
         }
         */
        chain->annealing=1.0;
        
        // (parallel) loop over chains
        //#pragma omp parallel for private(ic) shared(flags,model,trial,chain,orbit,proposal)
        for(int ic=0; ic<NC; ic++)
        {
            
            //loop over frequency segments
            for(int i=0; i<flags->NDATA; i++)
            {
                struct Model *model_ptr = model[chain->index[ic]][i];
                struct Model *trial_ptr = trial[chain->index[ic]];
                struct Data  *data_ptr  = data[i];
                
                
                for(int steps=0; steps < 100; steps++)
                {
                    //for(int j=0; j<model_ptr->Nlive; j++)
                    galactic_binary_mcmc(orbit, data_ptr, model_ptr, trial_ptr, chain, flags, prior, proposal[i], ic);
                    
                    noise_model_mcmc(orbit, data_ptr, model_ptr, trial_ptr, chain, flags, ic);
                }//loop over MCMC steps
                
                //reverse jump birth/death move
                if(flags->rj)galactic_binary_rjmcmc(orbit, data_ptr, model_ptr, trial_ptr, chain, flags, prior, proposal[i], ic);
                
                //update fisher matrix for each chain
                if(mcmc%100==0)
                {
                    for(int n=0; n<model_ptr->Nlive; n++)
                    {
                        galactic_binary_fisher(orbit, data_ptr, model_ptr->source[n], data_ptr->noise[FIXME]);
                    }
                }
            }//end loop over frequency segments
            
            //update start time for data segments
            if(flags->gap) data_mcmc(orbit, data, model[chain->index[ic]], chain, flags, proposal[0], ic);
            
        }// end (parallel) loop over chains
        
        ptmcmc(model,chain,flags);
        adapt_temperature_ladder(chain, mcmc+flags->NBURN);
        
        print_chain_files(data[FIXME], model, chain, flags, mcmc);
        
        //track maximum log Likelihood
        if(mcmc%100)
        {
            if(update_max_log_likelihood(model, chain, flags)) mcmc = -flags->NBURN;
        }
        
        //store reconstructed waveform
        if(!flags->quiet) print_waveform_draw(data, model[chain->index[0]], flags);
        
        //update run status
        if(mcmc%data[FIXME]->downsample==0)
        {
            
            if(!flags->quiet)
            {
                for(int i=0; i<flags->NDATA; i++)
                {
                    print_chain_state(data[i], chain, model[chain->index[0]][i], flags, stdout, mcmc);
                    fprintf(stdout,"Sources: %i\n",model[chain->index[0]][i]->Nlive);
                    print_acceptance_rates(proposal[i], chain->NP, 0, stdout);
                }
            }
            
            //save chain state to resume sampler
            save_chain_state(data, model, chain, flags, mcmc);
            
        }
        
        //dump waveforms to file, update avgLogL for thermodynamic integration
        if(mcmc>0 && mcmc%data[FIXME]->downsample==0)
        {
            for(int i=0; i<flags->NDATA; i++)save_waveforms(data[i], model[chain->index[0]][i], mcmc/data[i]->downsample);
            for(int ic=0; ic<NC; ic++)
            {
                chain->dimension[ic][model[chain->index[ic]][0]->Nlive]++;
                for(int i=0; i<flags->NDATA; i++)
                    chain->avgLogL[ic] += model[chain->index[ic]][i]->logL + model[chain->index[ic]][i]->logLnorm;
            }
        }
        
    }// end MCMC loop
    
    //print aggregate run files/results
    for(int i=0; i<flags->NDATA; i++)print_waveforms_reconstruction(data[i],i);
    
    FILE *chainFile = fopen("avg_log_likelihood.dat","w");
    for(int ic=0; ic<NC; ic++) fprintf(chainFile,"%lg %lg\n",1./chain->temperature[ic],chain->avgLogL[ic]/(double)(flags->NMCMC/data[FIXME]->downsample));
    fclose(chainFile);
    
    FILE *zFile = fopen("evidence.dat","w");
    for(int i=0; i<DMAX; i++) fprintf(zFile,"%i %i\n",i,chain->dimension[0][i]);
    fclose(zFile);
    
    //print total run time
    stop = time(NULL);
    
    printf(" ELAPSED TIME = %g second\n",(double)(stop-start));
    FILE *runlog = fopen("gb_mcmc.log","a");
    fprintf(runlog," ELAPSED TIME = %g second\n",(double)(stop-start));
    fclose(runlog);
    
    //free memory and exit cleanly
    for(int ic=0; ic<NC; ic++)
    {
        for(int i=0; i<flags->NDATA; i++) free_model(model[ic][i]);
        free_model(trial[ic]);
    }
    if(flags->orbit)free_orbit(orbit);
    //free_noise(data[0]->noise[FIXME]);
    //free_tdi(data[0]->tdi[FIXME]);
    free_chain(chain,flags);
    //free(model[FIXME][FIXME]);
    //free(trial[FIXME][FIXME]);
    //free(data[0]);
    
    return 0;
}

void ptmcmc(struct Model ***model, struct Chain *chain, struct Flags *flags)
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
        
        logL1 = 0.0;
        logL2 = 0.0;
        for(int i=0; i<flags->NDATA; i++)
        {
            logL1 += model[olda][i]->logL + model[olda][i]->logLnorm;
            logL2 += model[oldb][i]->logL + model[oldb][i]->logLnorm;
        }
        
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
    //double t0=100;
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
        //S[ic] += (A[ic][0] - A[ic][1])/nu;
        
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
        model_y->logL     = gaussian_log_likelihood(orbit, data, model_y);
        model_y->logLnorm = gaussian_log_likelihood_constant_norm(data, model_y);
        
        /*
         H = [p(d|y)/p(d|x)]/T x p(y)/p(x) x q(x|y)/q(y|x)
         */
        logH += ( (model_y->logL+model_y->logLnorm) - (model_x->logL+model_x->logLnorm) )/chain->temperature[ic]; //delta logL
        if(flags->burnin) logH /= chain->annealing;
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
            model_y->logL = gaussian_log_likelihood(orbit, data, model_y);
            
            /*
             H = [p(d|y)/p(d|x)]/T x p(y)/p(x) x q(x|y)/q(y|x)
             */
            logH += (model_y->logL - model_x->logL)/chain->temperature[ic]; //delta logL
            if(flags->burnin) logH /= chain->annealing;
        }
        logH += logPy  - logPx;  //priors
        logH += logQxy - logQyx; //proposals
        
        loga = log(gsl_rng_uniform(chain->r[ic]));
        
        
        if(isfinite(logH) && logH > loga)
        {
            
            //KAL? Print something for cov draw here??
            
            if(ic==0  && model_y->logL - model_x->logL < -20.)
            {
                
                printf("%s logH=%g, logLx=%g, logLy=%g\n",proposal[nprop]->name,logH,model_x->logL , model_y->logL);
                printf("   dlogQ=%g, logQxy=%g, logQyx=%g\n",logQxy - logQyx,logQxy,logQyx);
                
                /*
                 FILE *fptr = fopen("crazyhop.dat","a");
                 for(int k=0; k<data->NP; k++)fprintf(fptr,"%.12g ",source_x->params[k]);
                 for(int k=0; k<data->NP; k++)fprintf(fptr,"%.12g ",source_y->params[k]);
                 fprintf(fptr,"\n");
                 fclose(fptr);
                 */
                
                //exit(1);
            }
            
            
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
            //TODO: insert draw fro galaxy prior into draw_from_prior()
            logQyx = (*proposal[nprop]->function)(data, model_y, model_y->source[create], proposal[nprop], model_y->source[create]->params, chain->r[ic]);
            logQxy = 0;
            
            map_array_to_params(model_y->source[create], model_y->source[create]->params, data->T);
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
        model_y->logL = gaussian_log_likelihood(orbit, data, model_y);
        
        /*
         H = [p(d|y)/p(d|x)]/T x p(y)/p(x) x q(x|y)/q(y|x)
         */
        logH += (model_y->logL - model_x->logL)/chain->temperature[ic]; //delta logL
        if(flags->burnin) logH /= chain->annealing;
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
        
        if(ic==0  && model_y->logL - model_x->logL < -20.)
        {
            
            printf("RJ%s logH=%g, logLx=%g, logLy=%g\n",proposal[nprop]->name,logH,model_x->logL , model_y->logL);
            printf("   dlogQ=%g, logQxy=%g, logQyx=%g\n",logQxy - logQyx,logQxy,logQyx);
            //exit(1);
        }
        
        proposal[nprop]->accept[ic]++;
        copy_model(model_y,model_x);
    }
    
}

void data_mcmc(struct Orbit *orbit, struct Data **data, struct Model **model, struct Chain *chain, struct Flags *flags, struct Proposal **proposal, int ic)
{
    double logH  = 0.0; //(log) Hastings ratio
    double loga  = 1.0; //(log) transition probability
    double logQ  = 0.0;
    
    struct Model **trial = malloc(sizeof(struct Model *) * flags->NDATA);
    
    for(int i=0; i<flags->NDATA; i++)
    {
        trial[i] = malloc(sizeof(struct Model));
        
        
        alloc_model(trial[i],model[i]->Nmax,data[i]->N,data[i]->Nchannel, data[i]->NP, flags->NT);
        
        set_uniform_prior(flags, trial[i], data[i], 0);
        
        copy_model(model[i],trial[i]);
    }
    
    logQ += t0_shift(data[0], trial[0], trial[0]->source[0], proposal[0], trial[0]->source[0]->params, chain->r[ic]);
    
    for(int j=1; j<flags->NDATA; j++)
    {
        for(int i=0; i<flags->NT; i++)
        {
            trial[j]->t0[i] = trial[0]->t0[i];
            //      for(int n=0; n<trial[0]->Nlive; n++)
            //      {
            //        double dt = trial[j]->t0[i] - model[j]->t0[i];
            //        trial[0]->source[n]->params[0] += -(1.e-7/5.)*dt;
            //      }
        }
    }
    
    for(int j=0; j<flags->NDATA; j++)
    {
        // Form master template
        /*
         passing generate_signal_model -1 results in full recalculation of waveform model
         */
        generate_signal_model(orbit, data[j], trial[j], -1);
        //generate_signal_model(orbit, data[j], model[j], -1);
        
        /*
         H = [p(d|y)/p(d|x)]/T x p(y)/p(x) x q(x|y)/q(y|x)
         */
        if(!flags->prior)
        {
            // get likelihood for y
            trial[j]->logL = gaussian_log_likelihood(orbit, data[j], trial[j]);
            //model[j]->logL = gaussian_log_likelihood(orbit, data[j], model[j]);
            
            logH += (trial[j]->logL - model[j]->logL)/chain->temperature[ic];
            if(flags->burnin) logH /= chain->annealing;
        }
        logH += logQ; //delta logL
        
    }
    
    loga = log(gsl_rng_uniform(chain->r[ic]));
    
    if(logH > loga)
    {
        for(int j=0; j<flags->NDATA; j++)
        {
            copy_model(trial[j],model[j]);
        }
    }
    
    for(int i=0; i<flags->NDATA; i++)
    {
        free_model(trial[i]);
    }
    free(trial);
    
}




