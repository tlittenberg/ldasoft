//
//  UCBWrapper.c
//
//
//  Created by Tyson Littenberg on 1/28/21.
//

#include <mpi.h>
#include <omp.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdio.h>

#include <lisa.h>
#include <data.h>
#include <ucb.h>

#include "ucb_wrapper.h"



void alloc_ucb_data(struct UCBData *ucb_data, int procID)
{
    ucb_data->status  = 0;
    ucb_data->procID  = procID;
    ucb_data->flags   = malloc(sizeof(struct Flags));
    ucb_data->orbit   = malloc(sizeof(struct Orbit));
    ucb_data->chain   = malloc(sizeof(struct Chain));
    ucb_data->data    = malloc(sizeof(struct Data));
    ucb_data->prior   = malloc(sizeof(struct Prior));
    ucb_data->catalog = malloc(sizeof(struct Catalog));
}

void setup_ucb_data(struct UCBData *ucb_data, struct TDI *tdi_full)
{
    int procID     = ucb_data->procID;
    int procID_min = ucb_data->procID_min;
    int procID_max = ucb_data->procID_max;
    
    /* don't let procID go below procID_min (for frequency spacing) */
    if(procID<procID_min) procID = procID_min;
    
    /* Aliases to ucb structures */
    struct Flags *flags     = ucb_data->flags;
    struct Orbit *orbit     = ucb_data->orbit;
    struct Chain *chain     = ucb_data->chain;
    struct Data *data       = ucb_data->data;
    struct Catalog *catalog = ucb_data->catalog;
    
    /* silence output except for in highest-ranking process */
    flags->quiet = 1;
    if(procID==procID_max) flags->quiet = 0;
    
    /* Finish allocating UCB structures now that we know the number of PT chains */
    ucb_data->proposal = malloc(chain->NProp*sizeof(struct Proposal*));
    ucb_data->model = malloc(sizeof(struct Model*)*chain->NC);
    ucb_data->trial = malloc(sizeof(struct Model*)*chain->NC);

    /* Initialize LISA orbit model */
    initialize_orbit(data,orbit,flags);

    select_frequency_segment(data, tdi_full);
    
    /* Load gb catalog cache file for proposals/priors */
    if(flags->catalog)
        GalacticBinaryLoadCatalogCache(data, flags, catalog);
    
    /*
     Initialize measured time of model update.
     Used to determin number of steps relative to mbh model
     */
    ucb_data->cpu_time = 1.0;
}

#define GRIDFILE "ucb_frequency_spacing.dat"

void setup_frequency_segment(struct UCBData *ucb_data)
{
    int procID = ucb_data->procID;
    int procID_min = ucb_data->procID_min;
    int procID_max = ucb_data->procID_max;
    struct Data *data = ucb_data->data;

    if(access(GRIDFILE,F_OK)==0)
    {
        FILE *fgrid=fopen(GRIDFILE,"r");
        int n,id;
        double fstart,fstop;
        int fgridsize=0;
        while(!feof(fgrid))
        {
            fscanf(fgrid,"%i%i%lg%lg",&id,&n,&fstart,&fstop);
            fgridsize++;
        }
        rewind(fgrid);
        fgridsize--;
        
        double fmin[fgridsize];
        double fmax[fgridsize];
        
        for(int i=0; i<fgridsize; i++)
        {
            fscanf(fgrid,"%i%i%lg%lg",&id,&n,&fstart,&fstop);
            fmin[i] = fstart;
            fmax[i] = fstop;;
        }
        fclose(fgrid);
        
        data->fmin = fmin[procID-procID_min];
        data->fmax = fmax[procID-procID_min];
        data->N = (int)round((data->fmax-data->fmin)*data->T);
        data->qmin = (int)(data->fmin*data->T);
        data->qmax = data->qmin+data->N;
    }
    else
    {
        if(procID==procID_min) fprintf(stdout,"Did not find %s, manually setting up frequency grid\n",GRIDFILE);
        
        //how many ucb nodes
        int N_node = procID_max - procID_min + 1;
        
        //how many section sizes?
        int Smin =  (int)round(log(data->N-2*data->qpad)/log(2.));
        int Smax =  (int)round(log(data->Nmax-2*data->qpad)/log(2.));
        int N_seg = Smax - Smin + 1;//(int)round(log((double)Smax - (double)Smin + 1.)/log(2.));
        
        //integer part of nodes per section
        int n = (int)floor((double)N_node/(double)N_seg);
        
        //which section am I in?
        int k = (procID - procID_min)/n;
        if(k>N_seg-1) k = N_seg-1;
        
        //size of nodes in my section
        data->N = (int)round(pow(2,Smin+k));
        
        //start bin of my node?
        int Nsum=0;
        for(int node=procID_min; node<procID; node++)
        {
            k = (node - procID_min)/n;
            if(k>N_seg-1) k = N_seg-1;
            Nsum += (int)round(pow(2,Smin+k));
        }
        data->fmin = data->fmin + Nsum/data->T;
        data->fmax = data->fmin + data->N/data->T;
        data->qmin = (int)(data->fmin*data->T);
        data->qmax = data->qmin+data->N;
    }

    //add padding
    data->N += 2*data->qpad;
    data->qmin -= data->qpad;
    data->qmax += data->qpad;
    data->fmin = (double)data->qmin/data->T;
    data->fmax = (double)data->qmax/data->T;

}

void select_frequency_segment(struct Data *data, struct TDI *tdi_full)
{
    //store frequency segment in TDI structure
    struct TDI *tdi = data->tdi;
    struct TDI *raw = data->raw;
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
        
        /* raw data to be used for inter-segment ucb swaps */
        raw->X[n] = tdi_full->X[m];
        raw->Y[n] = tdi_full->Y[m];
        raw->Z[n] = tdi_full->Z[m];
        raw->A[n] = tdi_full->A[m];
        raw->E[n] = tdi_full->E[m];
        raw->T[n] = tdi_full->T[m];
    }
}

void initialize_ucb_sampler(struct UCBData *ucb_data)
{
    /* Aliases to ucb structures */
    struct Flags *flags = ucb_data->flags;
    struct Orbit *orbit = ucb_data->orbit;
    struct Chain *chain = ucb_data->chain;
    struct Data *data   = ucb_data->data;
    struct Prior *prior = ucb_data->prior;
    struct Proposal **proposal = ucb_data->proposal;
    struct Model **model = ucb_data->model;
    struct Model **trial = ucb_data->trial;
    struct Catalog *catalog = ucb_data->catalog;
    
    /* Lowest rank process has extra IO */
    if(ucb_data->procID==ucb_data->procID_min) flags->quiet=1;
    
    /* Get noise spectrum for data segment */
    GetNoiseModel(data,orbit,flags);
    
    /* Initialize parallel chain */
    if(flags->resume)
        initialize_chain(chain, flags, &data->cseed, "a");
    else
        initialize_chain(chain, flags, &data->cseed, "w");
    
    /* Initialize priors */
    if(flags->galaxyPrior) set_galaxy_prior(flags, prior);
    if(flags->update) set_gmm_prior(flags, data, prior, catalog);
    
    /* Initialize MCMC proposals */
    initialize_proposal(orbit, data, prior, chain, flags, catalog, proposal, flags->DMAX);
    
    /* Initialize UCB sampler state */
    struct Source *inj = NULL;
    initialize_ucb_state(data, orbit, flags, chain, proposal, model, trial, inj);
        
    /* Set sampler counter */
    ucb_data->mcmc_step = -flags->NBURN;
    
    /* Start analysis from saved chain state */
    if(flags->resume)
    {
        char filename[MAXSTRINGSIZE];
        
        //check for files needed to resume
        FILE *fptr = NULL;
        int file_error = 0;
        
        for(int ic=0; ic<chain->NC; ic++)
        {
            sprintf(filename,"%s/checkpoint/chain_state_%i.dat",flags->runDir,ic);
            
            if( (fptr = fopen(filename,"r")) == NULL )
            {
                fprintf(stderr,"Warning: Could not checkpoint run state for segment %i\n",ucb_data->procID);
                fprintf(stderr,"         Parameter file %s does not exist\n",filename);
                file_error++;
                break;
            }
        }
        
        //if all of the files exist resume run from checkpointed state
        if(!file_error) restore_chain_state(orbit, data, model, chain, flags, &ucb_data->mcmc_step);
    }
    
    /* Store data segment in working directory */
    print_data(data, data->tdi, flags);

    /* Store post-processing script */
    print_gb_catalog_script(flags, data, orbit);

}

static void print_sampler_state(struct UCBData *ucb_data)
{
    struct Model *model = ucb_data->model[ucb_data->chain->index[0]];
    fprintf(stdout,"UCB Process %i on step %i: sources = %i, logL = %g\n",
            ucb_data->procID,
            ucb_data->mcmc_step,
            model->Nlive,
            model->logL+model->logLnorm);
}

int update_ucb_sampler(struct UCBData *ucb_data)
{
    clock_t start = clock();

    /* Aliases to UCB structures */
    struct Flags *flags = ucb_data->flags;
    struct Orbit *orbit = ucb_data->orbit;
    struct Chain *chain = ucb_data->chain;
    struct Data *data   = ucb_data->data;
    struct Prior *prior = ucb_data->prior;
    struct Proposal **proposal = ucb_data->proposal;
    struct Model **model = ucb_data->model;
    struct Model **trial = ucb_data->trial;
    
    int NC = chain->NC;
    int mcmc_start = -flags->NBURN;
        
    /* exit if this segment is finished */
    if(ucb_data->mcmc_step >= flags->NMCMC) return 0;

    /* set flags based on current state of sampler */
    flags->burnin   = (ucb_data->mcmc_step<0) ? 1 : 0;
    flags->maximize = 0;//(ucb_data->mcmc_step<-flags->NBURN/2) ? 1 : 0;
    

    /* The MCMC loop */
    int numThreads;
    int numSteps = 10;
#pragma omp parallel num_threads(flags->threads)
    {
        //Save individual thread number
        int threadID = omp_get_thread_num();
        
        //Only one thread runs this section
        if(threadID==0) numThreads = omp_get_num_threads();
        
#pragma omp barrier
        // (parallel) loop over chains
        for(int ic=threadID; ic<NC; ic+=numThreads)
        {
            
            //loop over frequency segments
            struct Model *model_ptr = model[chain->index[ic]];
            struct Model *trial_ptr = trial[chain->index[ic]];
            
            //update model likelihood using new residual
            model_ptr->logL = gaussian_log_likelihood(data, model_ptr);
            model_ptr->logLnorm = gaussian_log_likelihood_model_norm(data, model_ptr);

            //sync up model and trial pointers
            copy_model(model_ptr,trial_ptr);

            for(int steps=0; steps<numSteps; steps++)
            {
                for(int m=0; m<100; m++)
                {
                    //reverse jump birth/death or split/merge moves
                    if(gsl_rng_uniform(chain->r[ic])<0.25 && flags->rj)
                        galactic_binary_rjmcmc(orbit, data, model_ptr, trial_ptr, chain, flags, prior, proposal, ic);
                    
                    //fixed dimension parameter updates
                    else
                        galactic_binary_mcmc(orbit, data, model_ptr, trial_ptr, chain, flags, prior, proposal, ic);
                }
            }
            
            //update fisher matrix for each chain
            for(int n=0; n<model_ptr->Nlive; n++)
            {
                galactic_binary_fisher(orbit, data, model_ptr->source[n], model_ptr->noise);
            }

            
        }// end (parallel) loop over chains
    }//end parallel section
#pragma omp barrier
    
    ptmcmc(model,chain,flags);
    adapt_temperature_ladder(chain, ucb_data->mcmc_step+flags->NBURN);
                    
    print_chain_files(data, model, chain, flags, ucb_data->mcmc_step);
        
    //track maximum log Likelihood
    if(update_max_log_likelihood(model, chain, flags))
        ucb_data->mcmc_step = -flags->NBURN;
    
    //update run status
    if(ucb_data->mcmc_step%data->downsample==0 && ucb_data->mcmc_step>mcmc_start)
    {
        
        //minimal screen output
        print_sampler_state(ucb_data);
        
        //save chain state to resume sampler
        save_chain_state(data, model, chain, flags, ucb_data->mcmc_step);
    }
    
    //dump waveforms to file, update avgLogL for thermodynamic integration
    if(ucb_data->mcmc_step>0 && ucb_data->mcmc_step%data->downsample==0)
    {
        save_waveforms(data, model[chain->index[0]], ucb_data->mcmc_step/data->downsample);
        
        for(int ic=0; ic<NC; ic++)
        {
            chain->dimension[ic][model[chain->index[ic]]->Nlive]++;
            for(int i=0; i<flags->NDATA; i++)
                chain->avgLogL[ic] += model[chain->index[ic]]->logL + model[chain->index[ic]]->logLnorm;
        }
    }
    
    ucb_data->mcmc_step+=numSteps;
    
    clock_t stop = clock();
    ucb_data->cpu_time = (double)(stop-start);
    
    return 1;
}

void exchange_ucb_source_params(struct UCBData *ucb_data)
{
    /* aliases for contents of UCBData structure */
    int procID     = ucb_data->procID;
    int procID_min = ucb_data->procID_min;
    int procID_max = ucb_data->procID_max;
    
    struct Flags *flags = ucb_data->flags;
    struct Orbit *orbit = ucb_data->orbit;
    struct Chain *chain = ucb_data->chain;
    struct Data  *data  = ucb_data->data;
    struct Model *model = ucb_data->model[chain->index[0]];
    
    int tag = 1; /* MPI tag for parameter exchange */
    
    //aliases for needed contents of model structure
    int Nshare = 0;
    int Nlive  = model->Nlive;
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
        alloc_model(new_model,Nlive_new,data->N,data->Nchannel);
        new_model->Nlive = Nlive_new;
        //set noise model
        copy_noise(data->noise, new_model->noise);
        
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
        for(int n=0; n<2*data->N; n++)
        {
            data->tdi->X[n] = data->raw->X[n] - new_model->tdi->X[n];
            data->tdi->Y[n] = data->raw->Y[n] - new_model->tdi->Y[n];
            data->tdi->Z[n] = data->raw->Z[n] - new_model->tdi->Z[n];
        }
                
        /* clean up */
        free_model(new_model);
    }
    
    /* update likelihoods */
    for(int ic=0; ic<chain->NC; ic++)
    {
        model = ucb_data->model[chain->index[ic]];
        if(!flags->prior)
        {
            model->logL     = gaussian_log_likelihood(data, model);
            model->logLnorm = gaussian_log_likelihood_constant_norm(data, model);
        }
        else model->logL = model->logLnorm = 0.0;
    }

}

int get_ucb_status(struct UCBData *ucb_data, int Nproc, int root, int procID)
{
    int UCB_Status = 0;
    int PIDmin = ucb_data->procID_min;
    int PIDmax = ucb_data->procID_max;
    
    /* worker nodes report status of their sampler */
    if(procID >= PIDmin && procID <= PIDmax)
        MPI_Send(&ucb_data->status, 1, MPI_INT, root, 0, MPI_COMM_WORLD);
    
    /* root node receives worker status, sums to find global status */
    if(procID==root)
    {
        UCB_Status = 0;
        for(int i=PIDmin; i<=PIDmax; i++)
        {
            MPI_Status status;
            MPI_Recv(&ucb_data->status, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
            UCB_Status += ucb_data->status;
        }
    }
    
    /* root node shares global status with worker nodes */
    MPI_Bcast(&UCB_Status, 1, MPI_INT, root, MPI_COMM_WORLD);
    
    return UCB_Status;
}

void print_ucb_state(struct UCBData *ucb_data, FILE *fptr, int counter)
{
    struct Data *data   = ucb_data->data;
    struct Model **model = ucb_data->model;
    struct Chain *chain = ucb_data->chain;
    
    int n = chain->index[0];

    fprintf(fptr,"%i %i ", counter, model[n]->Nlive);
    for(int i=0; i<model[n]->Nlive; i++)
    {
        print_source_params(data, model[n]->source[i], fptr);
    }
    fprintf(fptr,"\n");
}

