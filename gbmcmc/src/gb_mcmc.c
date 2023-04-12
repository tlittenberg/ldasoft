
/*
 *  Copyright (C) 2021 Tyson B. Littenberg (MSFC-ST12), Neil J. Cornish
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
 @file gb_mcmc.c
 \brief Main function for stand-alone GBMCMC sampler
 */

/*  REQUIRED LIBRARIES  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <sys/stat.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <omp.h>

#include <LISA.h>

#include "gbmcmc.h"

/**
 * This is the main function
 *
 */
int main(int argc, char *argv[])
{
    
    time_t start, stop;
    start = time(NULL);
    
    int NMAX = 10;   //max number of frequency & time segments

    char filename[PATH_BUFSIZE];

    /* check arguments */
    print_LISA_ASCII_art(stdout);
    print_version(stdout);
    if(argc==1) print_usage();
    
    
    /* Allocate data structures */
    struct Flags *flags = malloc(sizeof(struct Flags));
    struct Orbit *orbit = malloc(sizeof(struct Orbit));
    struct Chain *chain = malloc(sizeof(struct Chain));
    struct Data  *data = malloc(sizeof(struct Data));
        
    /* Parse command line and set defaults/flags */
    data->t0   = calloc( NMAX , sizeof(double) );
    data->tgap = calloc( NMAX , sizeof(double) );
    
    parse(argc,argv,data,orbit,flags,chain,NMAX);
    int NC = chain->NC;
    int DMAX = flags->DMAX;
    int mcmc_start = -flags->NBURN;
    
    /* Setup output directories for chain and data structures */
    pathprintf(data->dataDir,"%s/data",flags->runDir);
    pathprintf(chain->chainDir,"%s/chains",flags->runDir);
    pathprintf(chain->chkptDir,"%s/checkpoint",flags->runDir);
    
    mkdir(flags->runDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(data->dataDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain->chainDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain->chkptDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    /* Initialize data structures */
    alloc_data(data, flags);
    
    /* Initialize LISA orbit model */
    initialize_orbit(data, orbit, flags);

    /* Inject strain data */
    if(flags->strainData)
    {
        GalacticBinaryReadData(data,orbit,flags);
    }
    
    if(flags->NINJ>0)
    {
        /* Inject gravitational wave signal */
        if(flags->knownSource)
            GalacticBinaryInjectVerificationSource(data,orbit,flags);
        else
            GalacticBinaryInjectSimulatedSource(data,orbit,flags);
        
        /* set approximate f/fstar for segment */
        data->sine_f_on_fstar = sin((data->fmin + (data->fmax-data->fmin)/2.)/orbit->fstar);

    }
    
    
    /* Load catalog cache file for proposals/priors */
    if(flags->catalog)
    {
        GalacticBinaryLoadCatalogCache(data, flags);
        GalacticBinaryParseCatalogCache(data);
        GalacticBinaryLoadCatalog(data);
    }

    
    /* Initialize data-dependent proposal */
    setup_frequency_proposal(data, flags);
    
    /* Initialize parallel chain */
    if(flags->resume)
        initialize_chain(chain, flags, &data->cseed, "a");
    else
        initialize_chain(chain, flags, &data->cseed, "w");
    
    /* Initialize priors */
    struct Prior *prior = malloc(sizeof(struct Prior));
    if(flags->galaxyPrior) set_galaxy_prior(flags, prior);
    if(flags->update) set_gmm_prior(flags, data, prior);
    
    /* Initialize MCMC proposals */
    struct Proposal **proposal = malloc(chain->NP*sizeof(struct Proposal*));
    initialize_proposal(orbit, data, prior, chain, flags, proposal, DMAX);
    
    /* Test noise model */
    //test_noise_model(orbit);
    
    /* Initialize data models */
    struct Model **trial = malloc(sizeof(struct Model*)*NC);
    struct Model **model = malloc(sizeof(struct Model*)*NC);
    initialize_gbmcmc_state(data, orbit, flags, chain, proposal, model, trial);

    
    /* Start analysis from saved chain state */
    if(flags->resume)
    {
        fprintf(stdout,"\n=============== Checkpointing ===============\n");
        
        //check for files needed to resume
        FILE *fptr = NULL;
        int file_error = 0;
        
        for(int ic=0; ic<chain->NC; ic++)
        {
            pathprintf(filename,"%s/checkpoint/chain_state_%i.dat",flags->runDir,ic);
            
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
    
    /*test proposals
     FILE *test=fopen("proposal_test.dat","w");
     for(int i=0; i<100000; i++)
     {
     double logP = draw_from_gmm_prior(data, model[0], model[0]->source[0], proposal[7], model[0]->source[0]->params, chain->r[0]);
     print_source_params(data, model[0]->source[0], test);
     fprintf(test,"%lg\n",logP);
     }
     fclose(test);*/
    //exit(1);
    
    //test covariance proposal
    if(flags->updateCov) test_covariance_proposal(data, flags, model[0], prior, proposal[8], chain->r[0]);
    
    
    /* Write example gb_catalog bash script in run directory */
    print_gb_catalog_script(flags, data, orbit);
    
    //For saving the number of threads actually given
    int numThreads;
    int mcmc = mcmc_start;
    #pragma omp parallel num_threads(flags->threads)
    {
        int threadID;
        //Save individual thread number
        threadID = omp_get_thread_num();
        
        //Only one thread runs this section
        if(threadID==0)  numThreads = omp_get_num_threads();
        
        #pragma omp barrier
        
        /* The MCMC loop */
        for(; mcmc < flags->NMCMC;)
        {
            if(threadID==0)
            {
                flags->burnin   = (mcmc<0) ? 1 : 0;
                flags->maximize = 0;//(mcmc<-flags->NBURN/2) ? 1 : 0;
            }
            
            #pragma omp barrier
            // (parallel) loop over chains
            for(int ic=threadID; ic<NC; ic+=numThreads)
            {
                
                //loop over frequency segments
                struct Model *model_ptr = model[chain->index[ic]];
                struct Model *trial_ptr = trial[chain->index[ic]];
                copy_model(model_ptr,trial_ptr);
                
                for(int steps=0; steps < 100; steps++)
                {
                    //for(int j=0; j<model_ptr->Nlive; j++)
                    galactic_binary_mcmc(orbit, data, model_ptr, trial_ptr, chain, flags, prior, proposal, ic);
                                        
                }//loop over MCMC steps
                                
                //reverse jump birth/death move
                if(flags->rj)galactic_binary_rjmcmc(orbit, data, model_ptr, trial_ptr, chain, flags, prior, proposal, ic);
                
                if( (flags->strainData || flags->simNoise) && !flags->psd)
                    noise_model_mcmc(orbit, data, model_ptr, trial_ptr, chain, flags, ic);

                //update fisher matrix for each chain
                if(mcmc%100==0)
                {
                    for(int n=0; n<model_ptr->Nlive; n++)
                    {
                        galactic_binary_fisher(orbit, data, model_ptr->source[n], data->noise[FIXME]);
                    }
                }
                
                //update start time for data segments
                if(flags->gap) data_mcmc(orbit, data, model[chain->index[ic]], chain, flags, proposal, ic);
                
            }// end (parallel) loop over chains
            
            //Next section is single threaded. Every thread must get here before continuing
            #pragma omp barrier
            if(threadID==0){
                ptmcmc(model,chain,flags);
                adapt_temperature_ladder(chain, mcmc+flags->NBURN);
                
                print_chain_files(data, model, chain, flags, mcmc);
                
                //track maximum log Likelihood
                if(mcmc%100)
                {
                    if(update_max_log_likelihood(model, chain, flags)) mcmc = -flags->NBURN;
                }
                
                //store reconstructed waveform
                if(!flags->quiet) print_waveform_draw(data, model[chain->index[0]], flags);
                
                //update run status
                if(mcmc%data->downsample==0)
                {
                    
                    if(!flags->quiet)
                    {
                        print_chain_state(data, chain, model[chain->index[0]], flags, stdout, mcmc); //writing to file
                        fprintf(stdout,"Sources: %i\n",model[chain->index[0]]->Nlive);
                        print_acceptance_rates(proposal, chain->NP, 0, stdout);
                    }
                    
                    //save chain state to resume sampler
                    save_chain_state(data, model, chain, flags, mcmc);
                    
                }
                
                //dump waveforms to file, update avgLogL for thermodynamic integration
                if(mcmc>0 && mcmc%data->downsample==0)
                {
                    save_waveforms(data, model[chain->index[0]], mcmc/data->downsample);
                    
                    for(int ic=0; ic<NC; ic++)
                    {
                        chain->dimension[ic][model[chain->index[ic]]->Nlive]++;
                        for(int i=0; i<flags->NDATA; i++)
                        chain->avgLogL[ic] += model[chain->index[ic]]->logL + model[chain->index[ic]]->logLnorm;
                    }
                }
                mcmc++;
            }
            //Can't continue MCMC until single thread is finished
            #pragma omp barrier
            
        }// end MCMC loop
        
    }// End of parallelization
    
    //store final state of sampler
    save_chain_state(data, model, chain, flags, mcmc);

    //print aggregate run files/results
    print_waveforms_reconstruction(data,flags);
    print_noise_reconstruction(data,flags);
    print_evidence(chain,flags);
    
    pathprintf(filename,"%s/data/waveform_strain.dat",flags->runDir);
    FILE *waveFile = fopen(filename,"w");
    print_waveform_strain(data,model[chain->index[0]],waveFile);
    fclose(waveFile);


    pathprintf(filename,"%s/avg_log_likelihood.dat",flags->runDir);
    FILE *chainFile = fopen(filename,"w");
    for(int ic=0; ic<NC; ic++) fprintf(chainFile,"%lg %lg\n",1./chain->temperature[ic],chain->avgLogL[ic]/(double)(flags->NMCMC/data->downsample));
    fclose(chainFile);
    
    //print total run time
    stop = time(NULL);
    
    printf(" ELAPSED TIME = %g seconds on %i thread(s)\n",(double)(stop-start),numThreads);
    pathprintf(filename,"%s/gb_mcmc.log",flags->runDir);
    FILE *runlog = fopen(filename,"a");
    fprintf(runlog," ELAPSED TIME = %g seconds on %i thread(s)\n",(double)(stop-start),numThreads);
    fclose(runlog);
    
    //free memory and exit cleanly
    for(int ic=0; ic<NC; ic++)
    {
        free_model(model[ic]);
        free_model(trial[ic]);
    }
    if(flags->orbit)free_orbit(orbit);
    //free_noise(data->noise[FIXME]);
    //free_tdi(data->tdi[FIXME]);
    free_chain(chain,flags);
    //free(model[FIXME][FIXME]);
    //free(trial[FIXME][FIXME]);
    //free(data);
    
    return 0;
}

