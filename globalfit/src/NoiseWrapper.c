//
//  NoiseWrapper.c
//  
//
//  Created by Tyson Littenberg on 2/5/21.
//

#include <sys/stat.h>
#include <string.h>

#include <mpi.h>
#include <omp.h>

#include <LISA.h>
#include <GalacticBinary.h>
#include <GalacticBinaryIO.h>
#include <GalacticBinaryData.h>
#include <GalacticBinaryPrior.h>
#include <GalacticBinaryModel.h>
#include <GalacticBinaryProposal.h>
#include <GalacticBinaryMCMC.h>

#include "GalacticBinaryWrapper.h"
#include "NoiseWrapper.h"

void alloc_noise_data(struct NoiseData *noise_data, int procID)
{
    noise_data->status = 0;
    noise_data->procID = procID;
    noise_data->flags = NULL;//malloc(sizeof(struct Flags));
    noise_data->orbit = NULL;//malloc(sizeof(struct Orbit));
    noise_data->chain = malloc(sizeof(struct Chain));
    noise_data->data  = malloc(sizeof(struct Data));
}

void setup_noise_data(struct NoiseData *noise_data, struct GBMCMCData *gbmcmc_data, struct TDI *tdi_full)
{
    int procID = noise_data->procID;
    char dirname[MAXSTRINGSIZE];

    noise_data->data->downsample = gbmcmc_data->data->downsample;
    noise_data->data->Nwave      = gbmcmc_data->data->Nwave;
    
    
    noise_data->chain->NC      = gbmcmc_data->chain->NC;
    noise_data->data->Nchannel = gbmcmc_data->data->Nchannel;
    noise_data->data->NP       = gbmcmc_data->data->NP;
    strcpy(noise_data->data->format,gbmcmc_data->data->format);

    int qpad = gbmcmc_data->data->qpad;
    int N = gbmcmc_data->data->N - 2*qpad;
    int Nseg = gbmcmc_data->procID_max - gbmcmc_data->procID_min + 1;
    double T = gbmcmc_data->data->T;
    
    noise_data->data->T = T;
    noise_data->data->N = N*Nseg;
    noise_data->data->qpad = 0.0;

    //TODO: Fix this hacky noise model frequency alignment
    noise_data->data->fmin = gbmcmc_data->data->fmin + (double)qpad/T; //start noise model after padding
    noise_data->data->fmin += (double)N/T; //shift start by one segment (root doesn't have gbmcmc)


    noise_data->flags = malloc(sizeof(struct Flags));
    memcpy(noise_data->flags, gbmcmc_data->flags, sizeof(struct Flags));

    sprintf(noise_data->flags->runDir,"noise");
    mkdir(noise_data->flags->runDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    sprintf(dirname,"%s/chains",noise_data->flags->runDir);
    mkdir(dirname,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    sprintf(dirname,"%s/data",noise_data->flags->runDir);
    mkdir(dirname,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    noise_data->orbit = gbmcmc_data->orbit;
    
    alloc_data(noise_data->data, noise_data->flags);
    
    noise_data->model = malloc(sizeof(struct Model*)*gbmcmc_data->chain->NC);
    noise_data->trial = malloc(sizeof(struct Model*)*gbmcmc_data->chain->NC);
    
    select_frequency_segment(noise_data->data, tdi_full, procID);
}


void initialize_noise_sampler(struct NoiseData *noise_data)
{
    /* Aliases to gbmcmc structures */
    struct Flags *flags = noise_data->flags;
    struct Orbit *orbit = noise_data->orbit;
    struct Chain *chain = noise_data->chain;
    struct Data *data   = noise_data->data;
    struct Model **model = noise_data->model;
    struct Model **trial = noise_data->trial;
    
    /* Get noise spectrum for data segment */
    GalacticBinaryGetNoiseModel(data,orbit,flags);
    
    
    /* Initialize parallel chain */
    if(flags->resume)
        initialize_chain(chain, flags, &data->cseed, "a");
    else
        initialize_chain(chain, flags, &data->cseed, "w");
    
    /* Initialize GBMCMC sampler state */
    initialize_noise_state(data, orbit, flags, chain, model, trial);
    
    /* Set sampler counter */
    noise_data->mcmc_step = -flags->NBURN;
}

void initialize_noise_state(struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Chain *chain, struct Model **model, struct Model **trial)
{
    int NC = chain->NC;
    int DMAX = flags->DMAX;
    for(int ic=0; ic<NC; ic++)
    {
        
        trial[ic] = malloc(sizeof(struct Model));
        model[ic] = malloc(sizeof(struct Model));
        
        alloc_model(trial[ic],DMAX,data->N,data->Nchannel, data->NP, data->NT);
        alloc_model(model[ic],DMAX,data->N,data->Nchannel, data->NP, flags->NT);
                        
        //set noise model
        for(int j=0; j<flags->NT; j++) copy_noise(data->noise[j], model[ic]->noise[j]);
        
        // Form master model & compute likelihood of starting position
        generate_noise_model(data, model[ic]);
        generate_signal_model(orbit, data, model[ic], -1);
        
        if(!flags->prior)
        {
<<<<<<< HEAD
            model[ic]->logL     = gaussian_log_likelihood(data, model[ic]);
=======
            model[ic]->logL     = gaussian_log_likelihood(orbit, data, model[ic]);
>>>>>>> 3848e2f8e40eaadda41c5966f615a685041c9b4e
            model[ic]->logLnorm = gaussian_log_likelihood_constant_norm(data, model[ic]);
        }
        else model[ic]->logL = model[ic]->logLnorm = 0.0;
                
    }//end loop over chains

}

int update_noise_sampler(struct NoiseData *noise_data)
{
    /* Aliases to gbmcmc structures */
    struct Flags *flags = noise_data->flags;
    struct Orbit *orbit = noise_data->orbit;
    struct Chain *chain = noise_data->chain;
    struct Data *data   = noise_data->data;
    struct Model **model = noise_data->model;
    struct Model **trial = noise_data->trial;

    int NC = chain->NC;
    
    //For saving the number of threads actually given
    int numThreads;
    
#pragma omp parallel num_threads(flags->threads)
    {
        int threadID;
        //Save individual thread number
        threadID = omp_get_thread_num();
        
        //Only one thread runs this section
        if(threadID==0)  numThreads = omp_get_num_threads();
        
#pragma omp barrier
        /* The MCMC loop */
                
#pragma omp barrier
        // (parallel) loop over chains
        for(int ic=threadID; ic<NC; ic+=numThreads)
        {
            
            //loop over frequency segments
            struct Model *model_ptr = model[chain->index[ic]];
            struct Model *trial_ptr = trial[chain->index[ic]];
            
            //explicitly kill the signal model
            //loop over time segments
            for(int n=0; n<model_ptr->NT; n++)
            {
                for(int i=0; i<data->N*2; i++)
                {
                    model_ptr->tdi[n]->A[i] = 0.0;
                    model_ptr->tdi[n]->E[i] = 0.0;
                }
            }
            
            //update likelihood
<<<<<<< HEAD
            model_ptr->logL     = gaussian_log_likelihood(data, model_ptr);
=======
            model_ptr->logL     = gaussian_log_likelihood(orbit, data, model_ptr);
>>>>>>> 3848e2f8e40eaadda41c5966f615a685041c9b4e
            model_ptr->logLnorm = gaussian_log_likelihood_constant_norm(data, model_ptr);

            //evolve sampler
            for(int steps=0; steps < 100; steps++)
            {
                noise_model_mcmc(orbit, data, model_ptr, trial_ptr, chain, flags, ic);
            }
        }// end (parallel) loop over chains
         
        
        //Next section is single threaded. Every thread must get here before continuing
#pragma omp barrier
        if(threadID==0){
            ptmcmc(model,chain,flags);
            adapt_temperature_ladder(chain, noise_data->mcmc_step+flags->NBURN);
            noise_data->mcmc_step++;
        }
        //Can't continue MCMC until single thread is finished
#pragma omp barrier
        
    }// End of parallelization
    
    print_noise_state(data, model[chain->index[0]], chain->noiseFile[0], noise_data->mcmc_step);

    if(noise_data->mcmc_step>=0 && noise_data->mcmc_step%data->downsample==0 && noise_data->mcmc_step/data->downsample < data->Nwave)
    {
        struct Model *model_ptr = model[chain->index[0]];

        for(int n=0; n<data->N; n++)
        {
            data->S_pow[n][0][0][noise_data->mcmc_step/data->downsample] = model_ptr->noise[0]->SnA[n];
            data->S_pow[n][1][0][noise_data->mcmc_step/data->downsample] = model_ptr->noise[0]->SnE[n];
        }
    }
    return 1;
}
