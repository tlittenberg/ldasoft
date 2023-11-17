//
//  VerificationBinaryWrapper.c
//  ldasoft
//
//  Created by Tyson Littenberg on 8/11/21.
//

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>

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
#include "VerificationBinaryWrapper.h"

void alloc_vbmcmc_data(struct VBMCMCData *vbmcmc_data, struct GBMCMCData *gbmcmc_data, int procID)
{
    vbmcmc_data->status = 0;
    vbmcmc_data->procID = procID;
    vbmcmc_data->flags = malloc(sizeof(struct Flags));
    vbmcmc_data->orbit = malloc(sizeof(struct Orbit));

    memcpy(vbmcmc_data->flags, gbmcmc_data->flags, sizeof(struct Flags));

    vbmcmc_data->chain_vec = malloc(vbmcmc_data->flags->NVB*sizeof(struct Chain*));
    vbmcmc_data->data_vec  = malloc(vbmcmc_data->flags->NVB*sizeof(struct Data*));
    vbmcmc_data->prior_vec = malloc(vbmcmc_data->flags->NVB*sizeof(struct Prior*));
    vbmcmc_data->proposal_vec = malloc(vbmcmc_data->flags->NVB*sizeof(struct Proposal**));
    vbmcmc_data->trial_vec = malloc(vbmcmc_data->flags->NVB*sizeof(struct Model**));
    vbmcmc_data->model_vec = malloc(vbmcmc_data->flags->NVB*sizeof(struct Model**));
    
    for(int n=0; n<vbmcmc_data->flags->NVB; n++)
    {
        vbmcmc_data->data_vec[n] = malloc(sizeof(struct Data));
        copy_data(gbmcmc_data->data,vbmcmc_data->data_vec[n]);
        
        vbmcmc_data->chain_vec[n] = malloc(sizeof(struct Chain));
        vbmcmc_data->chain_vec[n]->NP = 7;
        vbmcmc_data->chain_vec[n]->NC = gbmcmc_data->chain->NC;
    }
}


void setup_vbmcmc_data(struct VBMCMCData *vbmcmc_data, struct GBMCMCData *gbmcmc_data, struct TDI *tdi_full)
{
    /* Aliases to gbmcmc structures */
    struct Flags *flags = vbmcmc_data->flags;
    struct Orbit *orbit = vbmcmc_data->orbit;
    struct Chain **chain_vec = vbmcmc_data->chain_vec;
    struct Data **data_vec   = vbmcmc_data->data_vec;
    
    /*
     * Set custom flags for verification binary analysis
     *
     * We are using the injection infrastructure to keep track
     * of the known binary parameters
     */
    flags->knownSource = 1;
    flags->snrPrior    = 0;
    flags->galaxyPrior = 0;
    flags->fixSky      = 1;
    flags->fixFreq     = 1;
    flags->cheat       = 1; //initializes chain at injection values
    flags->NINJ        = flags->NVB;
    flags->NBURN       = 0; //no burn in for vbmcmc
    
    /* parse verification binary files */
    FILE *vbFile = fopen(flags->vbFile,"r");
    
    //strip off header
    char header[MAXSTRINGSIZE];
    if(fgets(header, MAXSTRINGSIZE, vbFile)==NULL)
    {
        fprintf(stderr,"Error reading %s\n",flags->vbFile);
        exit(1);
    }
    
    /* Initialize LISA orbit model */
    initialize_orbit(gbmcmc_data->data, orbit, flags);
    
    /* initialize all data structures */
    for(int n=0; n<flags->NVB; n++)
    {
        
        struct Chain *chain=chain_vec[n];
        struct Data *data=data_vec[n];
        
        copy_data(gbmcmc_data->data,data);
        chain->NC = chain_vec[0]->NC; //number of chains
        
        
        /* Initialize data structures */
        alloc_data(data, flags);
        
        /* Get source from verification binary file */
        GetVerificationBinary(data, flags, vbFile);
        
        /* Pull out the right strain data */
        select_frequency_segment(data, tdi_full);
        
        /* set approximate f/fstar for segment */
        data->sine_f_on_fstar = sin((data->fmin + (data->fmax-data->fmin)/2.)/orbit->fstar);
    }
    
    /* Setup the rest of the model */
    for(int n=0; n<flags->NVB; n++)
    {
        vbmcmc_data->prior_vec[n] = malloc(sizeof(struct Prior));
        vbmcmc_data->proposal_vec[n] = malloc(vbmcmc_data->chain_vec[n]->NP*sizeof(struct Proposal*));
        vbmcmc_data->trial_vec[n] = malloc(sizeof(struct Model*)*vbmcmc_data->chain_vec[n]->NC);
        vbmcmc_data->model_vec[n] = malloc(sizeof(struct Model*)*vbmcmc_data->chain_vec[n]->NC);
    }
    
    /*
     Initialize measured time of model update.
     Used to determine number of steps relative to mbh model
     */
    vbmcmc_data->cpu_time = 1.0;
}


void initialize_vbmcmc_sampler(struct VBMCMCData *vbmcmc_data)
{
    struct Flags *flags = vbmcmc_data->flags;
    struct Orbit *orbit = vbmcmc_data->orbit;
    flags->quiet=1;
    
    for(int n=0; n<vbmcmc_data->flags->NVB; n++)
    {
        struct Chain *chain = vbmcmc_data->chain_vec[n];
        struct Data *data   = vbmcmc_data->data_vec[n];
        struct Prior *prior = vbmcmc_data->prior_vec[n];
        struct Proposal **proposal = vbmcmc_data->proposal_vec[n];
        struct Model **model = vbmcmc_data->model_vec[n];
        struct Model **trial = vbmcmc_data->trial_vec[n];
        
        /* Initialize parallel chain */
        if(flags->resume)
            initialize_chain(chain, flags, &data->cseed, "a");
        else
            initialize_chain(chain, flags, &data->cseed, "w");
        
        /* Initialize MCMC proposals */
        initialize_vb_proposal(orbit, data, prior, chain, flags, proposal, flags->DMAX);
        
        /* Initialize data models */
        initialize_gbmcmc_state(data, orbit, flags, chain, proposal, model, trial);
        
        /* Store data segment in working directory */
        if(vbmcmc_data->procID==1) print_data(data, data->tdi[0], flags, 0);
        
        /* Store post-processing script */
        print_gb_catalog_script(flags, data, orbit);
        
    }
    
    /* Set sampler counter */
    vbmcmc_data->mcmc_step = -flags->NBURN;
    
    
}

int update_vbmcmc_sampler(struct VBMCMCData *vbmcmc_data)
{
    clock_t start = clock();
    
    /* Aliases to gbmcmc structures */
    struct Flags *flags = vbmcmc_data->flags;
    struct Orbit *orbit = vbmcmc_data->orbit;
    struct Chain **chain_vec = vbmcmc_data->chain_vec;
    struct Data **data_vec   = vbmcmc_data->data_vec;
    struct Prior **prior_vec = vbmcmc_data->prior_vec;
    struct Proposal ***proposal_vec = vbmcmc_data->proposal_vec;
    struct Model ***model_vec = vbmcmc_data->model_vec;
    struct Model ***trial_vec = vbmcmc_data->trial_vec;
    
    struct Data *data = NULL;
    struct Chain *chain = NULL;
    struct Proposal **proposal = NULL;
    struct Model **model = NULL;
    
    int NC = chain_vec[0]->NC;
    int mcmc_start = -flags->NBURN;
    
    /* exit if this segment is finished */
    if(vbmcmc_data->mcmc_step >= flags->NMCMC) return 0;
    
    /* set flags based on current state of sampler */
    flags->burnin   = (vbmcmc_data->mcmc_step<0) ? 1 : 0;
    flags->maximize = 0;//(vbmcmc_data->mcmc_step<-flags->NBURN/2) ? 1 : 0;
    
    //For saving the number of threads actually given
    int numThreads;
    int numSteps=100;
#pragma omp parallel num_threads(flags->threads)
    {
        //Save individual thread number
        int threadID = omp_get_thread_num();;
        
        //Only one thread runs this section
        if(threadID==0)  numThreads = omp_get_num_threads();
        
#pragma omp barrier
        // (parallel) loop over chains
        for(int ic=threadID; ic<NC; ic+=numThreads)
        {
            
            for(int n=0; n<flags->NVB; n++)
            {
                
                struct Model *model_ptr = model_vec[n][chain_vec[n]->index[ic]];
                struct Model *trial_ptr = trial_vec[n][chain_vec[n]->index[ic]];
                
                model_ptr->logL = gaussian_log_likelihood(data_vec[n], model_ptr);
                model_ptr->logLnorm = gaussian_log_likelihood_model_norm(data_vec[n], model_ptr);

                //sync up model and trial pointers
                copy_model(model_ptr,trial_ptr);

                for(int steps=0; steps<numSteps; steps++)
                        galactic_binary_mcmc(orbit, data_vec[n], model_ptr, trial_ptr, chain_vec[n], flags, prior_vec[n], proposal_vec[n], ic);
            }
        }// end (parallel) loop over chains
    }//end parallel section
#pragma omp barrier
    
    for(int n=0; n<flags->NVB; n++)
    {
        model = model_vec[n];
        data = data_vec[n];
        proposal = proposal_vec[n];
        chain = chain_vec[n];
        
        //update fisher matrix for each chain
        for(int ic=0; ic<NC; ic++)
        {
            galactic_binary_fisher(orbit, data, model[chain->index[ic]]->source[0], model[chain->index[ic]]->noise[FIXME]);
        }

        ptmcmc(model_vec[n],chain_vec[n],flags);
        adapt_temperature_ladder(chain_vec[n], vbmcmc_data->mcmc_step+flags->NBURN);
        
        print_chain_files(data_vec[n], model_vec[n], chain_vec[n], flags, vbmcmc_data->mcmc_step);
        
        //track maximum log Likelihood
        if(update_max_log_likelihood(model, chain, flags))
            vbmcmc_data->mcmc_step = -flags->NBURN;
        
        //store reconstructed waveform
        if(!flags->quiet) print_waveform_draw(data, model[chain->index[0]], flags);
        
        //update run status
        if(vbmcmc_data->mcmc_step%data->downsample==0 && vbmcmc_data->mcmc_step>mcmc_start)
        {
            
            if(!flags->quiet)
            {
                print_chain_state(data, chain, model[chain->index[0]], flags, stdout, vbmcmc_data->mcmc_step); //writing to file
                fprintf(stdout,"Sources: %i\n",model[chain->index[0]]->Nlive);
                print_acceptance_rates(proposal, chain->NP, 0, stdout);
            }
            
            //save chain state to resume sampler
            save_chain_state(data, model, chain, flags, vbmcmc_data->mcmc_step);
            
        }
        
        //dump waveforms to file, update avgLogL for thermodynamic integration
        if(vbmcmc_data->mcmc_step%data->downsample==0)
        {
            save_waveforms(data, model[chain->index[0]], vbmcmc_data->mcmc_step/data->downsample);
            
            for(int ic=0; ic<NC; ic++)
            {
                chain->dimension[ic][model[chain->index[ic]]->Nlive]++;
                for(int i=0; i<flags->NDATA; i++)
                    chain->avgLogL[ic] += model[chain->index[ic]]->logL + model[chain->index[ic]]->logLnorm;
            }
        }
    }
    vbmcmc_data->mcmc_step+=numSteps;
    
    clock_t stop = clock();
    vbmcmc_data->cpu_time = (double)(stop-start);
    
    return 1;
}

void select_vbmcmc_segments(struct VBMCMCData *vbmcmc_data, struct TDI *tdi)
{
    for(int n=0; n<vbmcmc_data->flags->NVB; n++)
    {
        struct Data *data = vbmcmc_data->data_vec[n];
        select_frequency_segment(data,tdi);
        //struct Flags *flags = vbmcmc_data->flags;
        //print_data(data, data->tdi[0], flags, 0);
    }
}

void print_vbmcmc_state(struct VBMCMCData *vbmcmc_data, FILE *fptr, int counter)
{
    struct Data **data_vec   = vbmcmc_data->data_vec;
    struct Model ***model_vec = vbmcmc_data->model_vec;
    struct Chain **chain_vec = vbmcmc_data->chain_vec;
    struct Flags *flags = vbmcmc_data->flags;
    
    fprintf(fptr,"%i %i " ,counter, flags->NVB);
    
    for(int n=0; n<flags->NVB; n++)
    {
        struct Model **model = model_vec[n];
        struct Data *data = data_vec[n];
        struct Chain *chain = chain_vec[n];
        
        int m = chain->index[0];

        for(int i=0; i<model[m]->Nlive; i++)
        {
            print_source_params(data, model[m]->source[i], fptr);
        }
    }
    fprintf(fptr,"\n");

}
