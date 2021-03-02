//
//  GalacticBinaryWrapper.c
//  
//
//  Created by Tyson Littenberg on 1/28/21.
//

#include <mpi.h>
#include <omp.h>
#include <string.h>

#include <stdio.h>

#include <LISA.h>

#include <GalacticBinary.h>
#include <GalacticBinaryIO.h>
#include <GalacticBinaryData.h>
#include <GalacticBinaryPrior.h>
#include <GalacticBinaryModel.h>
#include <GalacticBinaryProposal.h>
#include <GalacticBinaryWaveform.h>
#include <GalacticBinaryMCMC.h>

#include "GalacticBinaryWrapper.h"

#define N_TDI_CHANNELS 2

void alloc_gbmcmc_data(struct GBMCMCData *gbmcmc_data, int procID, int procID_min, int procID_max)
{
    gbmcmc_data->status = 0;
    gbmcmc_data->procID = procID;
    gbmcmc_data->procID_min = procID_min;
    gbmcmc_data->procID_max = procID_max;
    gbmcmc_data->flags = malloc(sizeof(struct Flags));
    gbmcmc_data->orbit = malloc(sizeof(struct Orbit));
    gbmcmc_data->chain = malloc(sizeof(struct Chain));
    gbmcmc_data->data  = malloc(sizeof(struct Data));
    gbmcmc_data->prior = malloc(sizeof(struct Prior));
}

void select_frequency_segment(struct Data *data, struct TDI *tdi_full, int procID)
{
    //get max and min samples
    data->fmin = data->fmin + (double)(procID*(data->N - 2*data->qpad))/data->T;
    data->fmax = data->fmin + data->N/data->T;
    data->qmin = (int)(data->fmin*data->T);
    data->qmax = data->qmin+data->N;
    
    //store frequency segment in TDI structure
    struct TDI *tdi = data->tdi[0];
    struct TDI *raw = data->raw[0];
    for(int n=0; n<2*data->N; n++)
    {
        int m = data->qmin*2+n;
        
        /* data to be used by sampler */
        tdi->X[n] = tdi_full->X[m];
        tdi->Y[n] = tdi_full->Y[m];
        tdi->Z[n] = tdi_full->Z[m];
        tdi->A[n] = tdi_full->A[m];
        tdi->E[n] = tdi_full->E[m];
        tdi->T[n] = tdi_full->T[m];
        
        /* raw data to be used as reference */
        raw->X[n] = tdi_full->X[m];
        raw->Y[n] = tdi_full->Y[m];
        raw->Z[n] = tdi_full->Z[m];
        raw->A[n] = tdi_full->A[m];
        raw->E[n] = tdi_full->E[m];
        raw->T[n] = tdi_full->T[m];
    }
}

void get_frequency_segment(struct Data *data, struct TDI *tdi_full, int Nsamples, int root, int procID)
{
    //first tell all processes how large the dataset is
    MPI_Bcast(&Nsamples, 1, MPI_INT, root, MPI_COMM_WORLD);
    //MPI_Bcast(&data->T, 1, MPI_DOUBLE, root, MPI_COMM_WORLD); //only needed if read data maps to 2^N
    //MPI_Bcast(&data->sqT, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);

    //all but root process need to allocate memory for TDI structure
    if(procID!=root) alloc_tdi(tdi_full, Nsamples, N_TDI_CHANNELS);
    
    //now broadcast contents of TDI structure
    MPI_Bcast(&tdi_full->delta, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Bcast(tdi_full->X, 2*Nsamples, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Bcast(tdi_full->Y, 2*Nsamples, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Bcast(tdi_full->Z, 2*Nsamples, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Bcast(tdi_full->A, 2*Nsamples, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Bcast(tdi_full->E, 2*Nsamples, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Bcast(tdi_full->T, 2*Nsamples, MPI_DOUBLE, root, MPI_COMM_WORLD);
    
    /* select frequency segment for each process */
    select_frequency_segment(data, tdi_full, procID-1);
}

void broadcast_cache(struct Data *data, int root, int procID)
{

    /* broadcast number of each lines in the cache */
    MPI_Bcast(&data->Ncache, 1, MPI_INT, root, MPI_COMM_WORLD);

    /* all but root process need to allocate memory for cache sructure */
    if(procID!=root)
    {
        data->cache = malloc(data->Ncache*sizeof(char *));
        for(int n=0; n<data->Ncache; n++)
            data->cache[n] = (char *) malloc(MAXSTRINGSIZE);
    }
    
    /* broadcast each line */
    for(int n=0; n<data->Ncache; n++)
        MPI_Bcast(data->cache[n], MAXSTRINGSIZE, MPI_CHAR, root, MPI_COMM_WORLD);

}

void initialize_gbmcmc_sampler(struct GBMCMCData *gbmcmc_data)
{
    /* Aliases to gbmcmc structures */
    struct Flags *flags = gbmcmc_data->flags;
    struct Orbit *orbit = gbmcmc_data->orbit;
    struct Chain *chain = gbmcmc_data->chain;
    struct Data *data   = gbmcmc_data->data;
    struct Prior *prior = gbmcmc_data->prior;
    struct Proposal **proposal = gbmcmc_data->proposal;
    struct Model **model = gbmcmc_data->model;
    struct Model **trial = gbmcmc_data->trial;
    
    /* Lowest rank process has extra IO */
    if(gbmcmc_data->procID==gbmcmc_data->procID_min) flags->quiet=0;
    
    /* Get noise spectrum for data segment */
    GalacticBinaryGetNoiseModel(data,orbit,flags);
    
    /* Initialize parallel chain */
    if(flags->resume)
        initialize_chain(chain, flags, &data->cseed, "a");
    else
        initialize_chain(chain, flags, &data->cseed, "w");
    
    /* Initialize priors */
    if(flags->galaxyPrior) set_galaxy_prior(flags, prior);
    if(flags->update) set_gmm_prior(flags, data, prior);
    
    /* Initialize MCMC proposals */
    initialize_proposal(orbit, data, prior, chain, flags, proposal, flags->DMAX);
    
    /* Initialize GBMCMC sampler state */
    initialize_gbmcmc_state(data, orbit, flags, chain, proposal, model, trial);
    
    /* Set sampler counter */
    gbmcmc_data->mcmc_step = -flags->NBURN;
    
    /* Store data segment in working directory */
    print_data(data, data->tdi[0], flags, 0);
}

static void print_sampler_state(struct GBMCMCData *gbmcmc_data)
{
    struct Model *model = gbmcmc_data->model[gbmcmc_data->chain->index[0]];
    fprintf(stdout,"GBMCMC Process %i on step %i: sources = %i, logL = %g\n",
            gbmcmc_data->procID,
            gbmcmc_data->mcmc_step,
            model->Nlive,
            model->logL+model->logLnorm);
}

void run_gbmcmc_sampler(struct GBMCMCData *gbmcmc_data)
{
    /* Aliases to gbmcmc structures */
    struct Flags *flags = gbmcmc_data->flags;
    struct Orbit *orbit = gbmcmc_data->orbit;
    struct Chain *chain = gbmcmc_data->chain;
    struct Data *data   = gbmcmc_data->data;
    struct Prior *prior = gbmcmc_data->prior;
    struct Proposal **proposal = gbmcmc_data->proposal;
    struct Model **model = gbmcmc_data->model;
    struct Model **trial = gbmcmc_data->trial;
    
    int NC = chain->NC;
    int mcmc_start = -flags->NBURN;
    
    //For saving the number of threads actually given
    int numThreads;
    
    print_sampler_state(gbmcmc_data);
#pragma omp parallel num_threads(flags->threads)
    {
        int threadID;
        //Save individual thread number
        threadID = omp_get_thread_num();
        
        //Only one thread runs this section
        if(threadID==0)  numThreads = omp_get_num_threads();
        
#pragma omp barrier
        /* The MCMC loop */
        do
        {
            if(threadID==0)
            {
                if(gbmcmc_data->mcmc_step<0) flags->burnin=1;
                else       flags->burnin=0;
            }
            
#pragma omp barrier
            // (parallel) loop over chains
            for(int ic=threadID; ic<NC; ic+=numThreads)
            {
                
                //loop over frequency segments
                struct Model *model_ptr = model[chain->index[ic]];
                struct Model *trial_ptr = trial[chain->index[ic]];
                
                //update likelihood
<<<<<<< HEAD
                model_ptr->logL     = gaussian_log_likelihood(data, model_ptr);
=======
                model_ptr->logL     = gaussian_log_likelihood(orbit, data, model_ptr);
>>>>>>> 3848e2f8e40eaadda41c5966f615a685041c9b4e
                model_ptr->logLnorm = gaussian_log_likelihood_constant_norm(data, model_ptr);

                for(int steps=0; steps < 100; steps++)
                {
                    galactic_binary_mcmc(orbit, data, model_ptr, trial_ptr, chain, flags, prior, proposal, ic);
                }
                
                //reverse jump birth/death move
                if(flags->rj) galactic_binary_rjmcmc(orbit, data, model_ptr, trial_ptr, chain, flags, prior, proposal, ic);
                
                //update fisher matrix for each chain
                if(gbmcmc_data->mcmc_step%100==0)
                {
                    for(int n=0; n<model_ptr->Nlive; n++)
                    {
                        galactic_binary_fisher(orbit, data, model_ptr->source[n], data->noise[FIXME]);
                    }
                }
                
            }// end (parallel) loop over chains
             //Next section is single threaded. Every thread must get here before continuing
#pragma omp barrier
            
            if(threadID==0){
                ptmcmc(model,chain,flags);
                adapt_temperature_ladder(chain, gbmcmc_data->mcmc_step+flags->NBURN);
                
                print_chain_files(data, model, chain, flags, gbmcmc_data->mcmc_step);
                
                //track maximum log Likelihood
                if(gbmcmc_data->mcmc_step%100==0)
                {
                    if(update_max_log_likelihood(model, chain, flags))
                    {
                        gbmcmc_data->mcmc_step = -flags->NBURN;
                        //MPI_Bcast(&mcmc, 1, MPI_INT, gbmcmc_data->procID, MPI_COMM_WORLD);
                    }
                }
                
                //update run status
                if(gbmcmc_data->mcmc_step%data->downsample==0 && gbmcmc_data->mcmc_step>mcmc_start)
                {
                    
                    //minimal screen output
                    print_sampler_state(gbmcmc_data);
                    
                    //save chain state to resume sampler
                    save_chain_state(data, model, chain, flags, gbmcmc_data->mcmc_step);
                    
                }
                
                //dump waveforms to file, update avgLogL for thermodynamic integration
                if(gbmcmc_data->mcmc_step>0 && gbmcmc_data->mcmc_step%data->downsample==0)
                {
                    save_waveforms(data, model[chain->index[0]], gbmcmc_data->mcmc_step/data->downsample);
                    
                    for(int ic=0; ic<NC; ic++)
                    {
                        chain->dimension[ic][model[chain->index[ic]]->Nlive]++;
                        for(int i=0; i<flags->NDATA; i++)
                        chain->avgLogL[ic] += model[chain->index[ic]]->logL + model[chain->index[ic]]->logLnorm;
                    }
                }
                
                if(gbmcmc_data->mcmc_step%100==0)
                {
                    exchange_gbmcmc_source_params(gbmcmc_data);
                }
                gbmcmc_data->mcmc_step++;
            }
            //Can't continue MCMC until single thread is finished
#pragma omp barrier
            
        }while(gbmcmc_data->mcmc_step < flags->NMCMC); // end MCMC loop
        
    }// End of parallelization
}

int update_gbmcmc_sampler(struct GBMCMCData *gbmcmc_data)
{
    /* Aliases to gbmcmc structures */
    struct Flags *flags = gbmcmc_data->flags;
    struct Orbit *orbit = gbmcmc_data->orbit;
    struct Chain *chain = gbmcmc_data->chain;
    struct Data *data   = gbmcmc_data->data;
    struct Prior *prior = gbmcmc_data->prior;
    struct Proposal **proposal = gbmcmc_data->proposal;
    struct Model **model = gbmcmc_data->model;
    struct Model **trial = gbmcmc_data->trial;
    
    int NC = chain->NC;
    int mcmc_start = -flags->NBURN;
    
    /* Exchange parameters with neighbors */
    exchange_gbmcmc_source_params(gbmcmc_data);

    /* exit if this segment is finished */
    if(gbmcmc_data->mcmc_step >= flags->NMCMC) return 0;
    
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
        
        if(threadID==0)
        {
            flags->burnin   = (gbmcmc_data->mcmc_step<0) ? 1 : 0;
            flags->maximize = (gbmcmc_data->mcmc_step<-flags->NBURN/2) ? 1 : 0;
        }
        
#pragma omp barrier
        // (parallel) loop over chains
        for(int ic=threadID; ic<NC; ic+=numThreads)
        {
            
            //loop over frequency segments
            struct Model *model_ptr = model[chain->index[ic]];
            struct Model *trial_ptr = trial[chain->index[ic]];
            
            for(int steps=0; steps < 100; steps++)
            {
                galactic_binary_mcmc(orbit, data, model_ptr, trial_ptr, chain, flags, prior, proposal, ic);
            }
            
            //reverse jump birth/death move
            if(flags->rj) galactic_binary_rjmcmc(orbit, data, model_ptr, trial_ptr, chain, flags, prior, proposal, ic);
            
            //update fisher matrix for each chain
            if(gbmcmc_data->mcmc_step%100==0)
            {
                for(int n=0; n<model_ptr->Nlive; n++)
                {
                    galactic_binary_fisher(orbit, data, model_ptr->source[n], data->noise[FIXME]);
                }
            }
            
        }// end (parallel) loop over chains
         //Next section is single threaded. Every thread must get here before continuing
#pragma omp barrier
        
        if(threadID==0){
            ptmcmc(model,chain,flags);
            adapt_temperature_ladder(chain, gbmcmc_data->mcmc_step+flags->NBURN);
            
            print_chain_files(data, model, chain, flags, gbmcmc_data->mcmc_step);
            
            //track maximum log Likelihood
            if(gbmcmc_data->mcmc_step%100==0)
            {
                if(update_max_log_likelihood(model, chain, flags))
                {
                    gbmcmc_data->mcmc_step = -flags->NBURN;
                }
            }
            
            //update run status
            if(gbmcmc_data->mcmc_step%data->downsample==0 && gbmcmc_data->mcmc_step>mcmc_start)
            {
                
                //minimal screen output
                print_sampler_state(gbmcmc_data);
                
                //save chain state to resume sampler
                save_chain_state(data, model, chain, flags, gbmcmc_data->mcmc_step);
                
            }
            
            //dump waveforms to file, update avgLogL for thermodynamic integration
            if(gbmcmc_data->mcmc_step>0 && gbmcmc_data->mcmc_step%data->downsample==0)
            {
                save_waveforms(data, model[chain->index[0]], gbmcmc_data->mcmc_step/data->downsample);
                
                for(int ic=0; ic<NC; ic++)
                {
                    chain->dimension[ic][model[chain->index[ic]]->Nlive]++;
                    for(int i=0; i<flags->NDATA; i++)
                    chain->avgLogL[ic] += model[chain->index[ic]]->logL + model[chain->index[ic]]->logLnorm;
                }
            }
            
            gbmcmc_data->mcmc_step++;
        }
        //Can't continue MCMC until single thread is finished
#pragma omp barrier
        
    }// End of parallelization
    
    return 1;
}

void exchange_gbmcmc_source_params(struct GBMCMCData *gbmcmc_data)
{
    /* aliases for contents of GBMCMCdata structure */
    int procID     = gbmcmc_data->procID;
    int procID_min = gbmcmc_data->procID_min;
    int procID_max = gbmcmc_data->procID_max;
    
    struct Flags *flags = gbmcmc_data->flags;
    struct Orbit *orbit = gbmcmc_data->orbit;
    struct Chain *chain = gbmcmc_data->chain;
    struct Data  *data  = gbmcmc_data->data;
    struct Model *model = gbmcmc_data->model[chain->index[0]];
    
    int tag = 1; /* MPI tag for parameter exchange */
    
    //aliases for needed contents of model structure
    int Nshare = 0;
    int Nlive  = model->Nlive;
    int NP = model->NP;
    int share_flag[Nlive];
    struct Source **source = model->source;
    
    /* count how many sources are shareable */
    for(int n=0; n<Nlive; n++)
    {
        share_flag[n]=0;
        //find central bin of source
        double q_sample = model->source[n]->f0*data->T;
        
        //only share if the source is not in the padding region
        if(q_sample > data->qmin+data->qpad && q_sample < data->qmax-data->qpad)
        {
            share_flag[n]=1;
            Nshare++;
        }
    }
    
    //build array of all source parameters to ship
    int Nparams = Nshare*NP;
    double params[Nparams]; //TODO: Allow for different number of parameters for each source
    Nshare = 0;
    for(int i=0; i<Nlive; i++)
    {
        if(share_flag[i])
        {
            for(int n=0; n<NP; n++)
            {
                params[Nshare*NP+n] = source[i]->params[n];
            }
            Nshare++;
            
        }
    }
    
    /* send to either side */
    int left_neighbor = procID-1;
    int right_neighbor = procID+1;
    if(left_neighbor>=procID_min)
        MPI_Send(&params, Nparams, MPI_DOUBLE, left_neighbor, tag, MPI_COMM_WORLD);
    if(right_neighbor<=procID_max)
        MPI_Send(&params, Nparams, MPI_DOUBLE, right_neighbor, tag, MPI_COMM_WORLD);
    
    /* probe for messages from either side*/
    MPI_Status status;
    MPI_Status status_left;
    MPI_Status status_right;
    
    int Nparams_left =0;
    int Nparams_right=0;
    
    if(left_neighbor>=procID_min)
    {   MPI_Probe(left_neighbor, tag, MPI_COMM_WORLD, &status_left);
        MPI_Get_count(&status_left, MPI_DOUBLE, &Nparams_left);
    }
    if(right_neighbor<=procID_max)
    {
        MPI_Probe(right_neighbor, tag, MPI_COMM_WORLD, &status_right);
        MPI_Get_count(&status_right, MPI_DOUBLE, &Nparams_right);
    }
    
    /* create arrays to receive messages */
    double params_right[Nparams_right];
    double params_left[Nparams_left];
    
    /* recieve */
    if(left_neighbor>=procID_min)
        MPI_Recv(&params_left, Nparams_left, MPI_DOUBLE, left_neighbor, tag, MPI_COMM_WORLD,&status);
    if(right_neighbor<=procID_max)
        MPI_Recv(&params_right, Nparams_right, MPI_DOUBLE, right_neighbor, tag, MPI_COMM_WORLD,&status);
    
    /* populate new model structure with neighboring waveforms */
    int Nparams_new = Nparams_left + Nparams_right;
    int Nlive_new = Nparams_new/NP;
    
    if(Nlive_new > 0)
    {
        struct Model *new_model = malloc(sizeof(struct Model));
        alloc_model(new_model,Nlive_new,data->N,data->Nchannel,NP,flags->NT);
        new_model->Nlive = Nlive_new;
        //set noise model
        for(int j=0; j<flags->NT; j++) copy_noise(data->noise[j], new_model->noise[j]);
        
        /* unpack recieved parameter vectors into source structures and generate waveforms */
        int m=0;
        for(int n=0; n<Nparams_left/NP; n++)
        {
            for(int i=0; i<NP; i++)
            {
                new_model->source[n]->params[i] = params_left[n*NP+i];
            }
        }
        for(int n=0; n<Nparams_right/NP; n++)
        {
            for(int i=0; i<NP; i++)
            {
                m = Nparams_left/NP + n;
                new_model->source[m]->params[i] = params_right[n*NP+i];
            }
        }
        
        for(int n=0; n<Nlive_new; n++)
        {
            map_array_to_params(new_model->source[n], new_model->source[n]->params, data->T);
        }
        
        /* generate signal model for received parameters */
        generate_signal_model(orbit, data, new_model, -1);
        
        /* form residual */
        for(int i=0; i<data->NT; i++)
        {
            for(int n=0; n<2*data->N; n++)
            {
                //TODO: need support for X,Y,Z,T channels
                data->tdi[i]->A[n] = data->raw[i]->A[n] - new_model->tdi[i]->A[n];
                data->tdi[i]->E[n] = data->raw[i]->E[n] - new_model->tdi[i]->E[n];
            }
        }
        
        /* update likelihoods */
        for(int ic=0; ic<chain->NC; ic++)
        {
            model = gbmcmc_data->model[chain->index[ic]];
            if(!flags->prior)
            {
<<<<<<< HEAD
                model->logL     = gaussian_log_likelihood(data, model);
=======
                model->logL     = gaussian_log_likelihood(orbit, data, model);
>>>>>>> 3848e2f8e40eaadda41c5966f615a685041c9b4e
                model->logLnorm = gaussian_log_likelihood_constant_norm(data, model);
            }
            else model->logL = model->logLnorm = 0.0;
        }
        
        /* clean up */
        free_model(new_model);
    }
}

int get_gbmcmc_status(struct GBMCMCData *gbmcmc_data, int Nproc, int root, int procID)
{
    int GBMCMC_Status = 0;
    int PIDmin = gbmcmc_data->procID_min;
    int PIDmax = gbmcmc_data->procID_max;
    
    /* worker nodes report status of their sampler */
    if(procID >= PIDmin && procID <= PIDmax)
        MPI_Send(&gbmcmc_data->status, 1, MPI_INT, root, 0, MPI_COMM_WORLD);
    
    /* root node receives worker status, sums to find global status */
    if(procID==root)
    {
        GBMCMC_Status = 0;
        for(int i=PIDmin; i<=PIDmax; i++)
        {
            MPI_Status status;
            MPI_Recv(&gbmcmc_data->status, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
            GBMCMC_Status += gbmcmc_data->status;
        }
    }
    
    /* root node shares global status with worker nodes */
    MPI_Bcast(&GBMCMC_Status, 1, MPI_INT, root, MPI_COMM_WORLD);
    
    return GBMCMC_Status;
}

