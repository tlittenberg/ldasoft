
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
 @file vb_mcmc.c
 \brief Main function for dedicated verification binary sampler
 */

/*  REQUIRED LIBRARIES  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <getopt.h>

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
    char filename[MAXSTRINGSIZE];

    /* check arguments */
    print_LISA_ASCII_art(stdout);
    print_version(stdout);
    if(argc==1) print_usage();

    
    /* Allocate data structures */
    struct Flags *flags = malloc(sizeof(struct Flags));
    struct Orbit *orbit = malloc(sizeof(struct Orbit));

    /* get vbmcmc-specific arguments before allocating Data structure*/
    parse_vb_list(argc,argv,flags);
    if(flags->NVB==0)
    {
        fprintf(stdout,"Verification binary list required\n");
        fprintf(stdout, "Example ./vb_mcmc --known-sources /path/to/full_list.txt\n");
        return 0;
    }

    struct Data *data;
    struct Chain *chain;
    struct Data  **data_vec = malloc(flags->NVB*sizeof(struct Data*));
    struct Chain **chain_vec = malloc(flags->NVB*sizeof(struct Chain*));
    for(int n=0; n<flags->NVB; n++)
    {
        data_vec[n] = malloc(sizeof(struct Data));
        chain_vec[n] = malloc(sizeof(struct Chain));
    }
        
    data=data_vec[0];
    chain=chain_vec[0];
    parse(argc,argv,data,orbit,flags,chain,NMAX);
    int NC = chain->NC;
    int DMAX = flags->DMAX;
    int mcmc_start = -flags->NBURN;
    
    /*
     * Set custom flags for verification binary analysis
     *
     * We are using the injection infrastructure to keep track
     * of the known binary parameters
     */
    flags->knownSource = 1;
    flags->snrPrior    = 0;
    flags->fixSky      = 1;
    flags->fixFreq     = 1;
    flags->cheat       = 1; //initializes chain at injection values
    flags->NINJ        = flags->NVB;
    
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
    initialize_orbit(data, orbit, flags);

    
    /* initialize all data structures */

    /* Setup output directories for data and chain files */
    mkdir(flags->runDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    for(int n=0; n<flags->NVB; n++)
    {
        chain=chain_vec[n];
        data=data_vec[n];
        data->nseed+=n;

        char subDir[PATH_BUFSIZE];
        pathprintf(subDir,"%s/seg%02d",flags->runDir,n);
        mkdir(subDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        pathprintf(data->dataDir,"%s/data",subDir);
        pathprintf(chain->chainDir,"%s/chains",subDir);
        pathprintf(chain->chkptDir,"%s/checkpoint",subDir);
        
        mkdir(data->dataDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        mkdir(chain->chainDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        mkdir(chain->chkptDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        if(n>0)
        {
            copy_data(data_vec[0],data);
            chain->NP = chain_vec[0]->NP; //number of proposals
            chain->NC = chain_vec[0]->NC; //number of chains

        }
        
        /* Initialize data structures */
        alloc_data(data, flags);
        
        /* Get source from verification binary file */
        GetVerificationBinary(data, flags, vbFile);
        
        /* Read strain data */
        if(flags->hdf5Data)
            GalacticBinaryReadData(data,orbit,flags);

        /* Simulate strain data */
        else
            GalacticBinaryInjectVerificationSet(data, orbit, flags);
        
        /* set approximate f/fstar for segment */
        data->sine_f_on_fstar = sin((data->fmin + (data->fmax-data->fmin)/2.)/orbit->fstar);

        //print various data products for plotting
        print_data(data, data->tdi[0], flags, 0);
        
        
        /* Initialize parallel chain */
        if(flags->resume)
            initialize_chain(chain, flags, &data->cseed, "a");
        else
            initialize_chain(chain, flags, &data->cseed, "w");
        
    }

    /* Setup the rest of the model */
    struct Proposal **proposal = NULL;
    struct Model **model = NULL;
    struct Prior **prior_vec = malloc(flags->NVB*sizeof(struct Prior *));
    struct Proposal ***proposal_vec = malloc(flags->NVB*sizeof(struct Proposal **));
    struct Model ***trial_vec = malloc(flags->NVB*sizeof(struct Model**));
    struct Model ***model_vec = malloc(flags->NVB*sizeof(struct Model**));

    for(int n=0; n<flags->NVB; n++)
    {
        /* Initialize priors */
        prior_vec[n] = malloc(sizeof(struct Prior));
        
        /* Initialize MCMC proposals */
        proposal_vec[n] = malloc(chain->NP*sizeof(struct Proposal*));
        initialize_vb_proposal(orbit, data_vec[n], prior_vec[n], chain_vec[n], flags, proposal_vec[n], DMAX);
        
        /* Initialize data models */
        trial_vec[n] = malloc(sizeof(struct Model*)*NC);
        model_vec[n] = malloc(sizeof(struct Model*)*NC);
        initialize_gbmcmc_state(data_vec[n], orbit, flags, chain_vec[n], proposal_vec[n], model_vec[n], trial_vec[n]);
    }
    
    /* Start analysis from saved chain state */
    if(flags->resume)
    {
        fprintf(stdout,"\n=============== Checkpointing ===============\n");
        for(int n=0; n<flags->NVB; n++)
        {
            //check for files needed to resume
            FILE *fptr = NULL;
            int file_error = 0;
            
            for(int ic=0; ic<chain_vec[n]->NC; ic++)
            {
                pathprintf(filename,"%s/chain_state_%i.dat",chain_vec[n]->chkptDir,ic);
                
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
                restore_chain_state(orbit, data_vec[n], model_vec[n], chain_vec[n], flags, &mcmc_start);
            }
        }
        fprintf(stdout,"============================================\n\n");
    }
            
    /* Write example gb_catalog bash script in run directory */
    print_gb_catalog_script(flags, data_vec[0], orbit);
    
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
                flags->maximize = (mcmc<-flags->NBURN/2) ? 1 : 0;
            }
            
            #pragma omp barrier
            // (parallel) loop over chains
            for(int ic=threadID; ic<NC; ic+=numThreads)
            {
                
                //loop over verification binary segments
                for(int n=0; n<flags->NVB; n++)
                {
                    
                    struct Model *model_ptr = model_vec[n][chain_vec[n]->index[ic]];
                    struct Model *trial_ptr = trial_vec[n][chain_vec[n]->index[ic]];
                    
                    for(int steps=0; steps < 100; steps++)
                    {
                        galactic_binary_mcmc(orbit, data_vec[n], model_ptr, trial_ptr, chain_vec[n], flags, prior_vec[n], proposal_vec[n], ic);
                    }//loop over MCMC steps
                    
                    //update fisher matrix for each chain
                    if(mcmc%100==0)
                    {
                        for(int i=0; i<model_ptr->Nlive; i++)
                        {
                            galactic_binary_fisher(orbit, data_vec[n], model_ptr->source[i], data_vec[n]->noise[FIXME]);
                        }
                    }
                }
                
            }// end (parallel) loop over chains
            
            //Next section is single threaded. Every thread must get here before continuing
            #pragma omp barrier
            if(threadID==0){
                
                for(int n=0; n<flags->NVB; n++)
                {
                    model = model_vec[n];
                    data = data_vec[n];
                    proposal = proposal_vec[n];
                    chain = chain_vec[n];
                    
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
                }
            mcmc++;
            }
            //Can't continue MCMC until single thread is finished
#pragma omp barrier
            
        }// end MCMC loop
        
    }// End of parallelization
    
    //print aggregate run files/results
    for(int n=0; n<flags->NVB; n++)
    {
        print_waveforms_reconstruction(data_vec[n],flags);
        print_noise_reconstruction(data_vec[n],flags);
    }
    
    
    //print total run time
    stop = time(NULL);
    
    printf(" ELAPSED TIME = %g seconds on %i thread(s)\n",(double)(stop-start),numThreads);
    pathprintf(filename,"%s/vb_mcmc.log",flags->runDir);
    FILE *runlog = fopen(filename,"a");
    fprintf(runlog," ELAPSED TIME = %g seconds on %i thread(s)\n",(double)(stop-start),numThreads);
    fclose(runlog);
    
    
    return 0;
}



