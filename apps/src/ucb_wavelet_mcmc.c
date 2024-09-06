/*
 *  Copyright (C) 2024 Tyson B. Littenberg (MSFC-ST12)
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
 @file ucb_wavelet_mcmc.c
 \brief Main function for wavelet-domain UCB sampler app `ucb_wavelet_mcmc` 
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

#include <glass_utils.h>
#include <glass_noise.h>
#include <glass_ucb.h>

static void print_usage()
{
    print_glass_usage();
    print_ucb_usage();
    
    fprintf(stdout,"EXAMPLE:\n");
    fprintf(stdout,"ucb_wavelet_mcmc --inj [path to]/ldasoft/ucb/etc/sources/precision/PrecisionSource_0.txt --cheat\n");
    fprintf(stdout,"\n");

    exit(0);
}

/**
 * This is the main function
 *
 */
int main(int argc, char *argv[])
{
    fprintf(stdout, "\n============= Wavelet-domain UCB MCMC ============\n");

    time_t start, stop;
    start = time(NULL);
    char filename[MAXSTRINGSIZE];

    /* check arguments */
    print_LISA_ASCII_art(stdout);
    print_version(stdout);
    if(argc==1) print_usage();
    
    /* Allocate data structures */
    struct Flags  *flags = malloc(sizeof(struct Flags));
    struct Orbit  *orbit = malloc(sizeof(struct Orbit));
    struct Chain  *chain = malloc(sizeof(struct Chain));
    struct Data   *data  = malloc(sizeof(struct Data));
    
    /* Parse command line and set defaults/flags */
    parse_data_args(argc,argv,data,orbit,flags,chain,"wavelet");
    parse_ucb_args(argc,argv,flags);
    if(flags->help) print_usage();

    int NC = chain->NC;
    int DMAX = flags->DMAX;
    int mcmc_start = -flags->NBURN;

    /* Setup output directories for chain and data structures */
    sprintf(data->dataDir,"%s/data",flags->runDir);
    sprintf(chain->chainDir,"%s/chains",flags->runDir);
    sprintf(chain->chkptDir,"%s/checkpoint",flags->runDir);
    
    mkdir(flags->runDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(data->dataDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain->chainDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain->chkptDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    /* Initialize data structures */
    alloc_data(data, flags);

    /*
    Define and set up Orbit structure which contains spacecraft ephemerides
    */
    initialize_orbit(data, orbit, flags);
    initialize_interpolated_analytic_orbits(orbit, data->T, data->t0);

    /*
    Only supporting software injections r/n
    */
    struct Source **inj=malloc(DMAX*sizeof(struct Source*));
    for(int n=0; n<DMAX; n++) inj[n] = malloc(sizeof(struct Source));
    
    UCBInjectSimulatedSource(data,orbit,flags,inj);

    /* Load catalog cache file for proposals/priors */
    struct Catalog *catalog=malloc(sizeof(struct Catalog));
    if(flags->catalog)
        UCBLoadCatalogCache(data, flags, catalog);

    

    /* Initialize parallel chain */
    if(flags->resume)
        initialize_chain(chain, flags, &data->cseed, "a");
    else
        initialize_chain(chain, flags, &data->cseed, "w");
    
    /* Initialize priors */
    struct Prior *prior = malloc(sizeof(struct Prior));
    if(flags->galaxyPrior) set_galaxy_prior(flags, prior);
    if(flags->update) set_gmm_prior(flags, data, prior, catalog);

    /* Initialize MCMC proposals */
    struct Proposal **proposal = malloc(UCB_PROPOSAL_NPROP*sizeof(struct Proposal*));
    initialize_proposal(orbit, data, prior, chain, flags, catalog, proposal, DMAX);
    
    
    /* Initialize data models */
    struct Model **trial = malloc(sizeof(struct Model*)*NC);
    struct Model **model = malloc(sizeof(struct Model*)*NC);
    initialize_ucb_state(data, orbit, flags, chain, proposal, model, trial, inj);

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
                    ucb_mcmc(orbit, data, model_ptr, trial_ptr, chain, flags, prior, proposal, ic);
                                                
                //update fisher matrix for each chain
                for(int n=0; n<model_ptr->Nlive; n++)
                    ucb_fisher_wavelet(orbit, data, model_ptr->source[n], data->noise);
                
            }// end (parallel) loop over chains
            
            //Next section is single threaded. Every thread must get here before continuing
            #pragma omp barrier
            if(threadID==0)
            {
                ptmcmc(model,chain,flags);
                adapt_temperature_ladder(chain, mcmc+flags->NBURN);
                
                print_chain_files(data, model, chain, flags, mcmc);
                
                //track maximum log Likelihood
                if(mcmc%100)
                {
                    if(update_max_log_likelihood(model, chain, flags)) mcmc = -flags->NBURN;
                }
                                
                //update run status
                if(mcmc%data->downsample==0)
                {
                    
                    if(!flags->quiet)
                    {
                        print_chain_state(data, chain, model[chain->index[0]], flags, stdout, mcmc); //writing to file
                        fprintf(stdout,"Sources: %i/%i\n",model[chain->index[0]]->Nlive,model[chain->index[0]]->Neff-1);
                        print_acceptance_rates(proposal, UCB_PROPOSAL_NPROP, 0, stdout);
                    }
                    
                    //save chain state to resume sampler
                    save_chain_state(data, model, chain, flags, mcmc);
                    
                }
                
                if(mcmc>-flags->NBURN+flags->NBURN/10. && model[0]->Neff < model[0]->Nmax && flags->rj)
                {
                    for(int ic=0; ic<NC; ic++) model[ic]->Neff++;
                    mcmc = -flags->NBURN;
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
    print_evidence(chain,flags);
    
    sprintf(filename,"%s/data/waveform_strain.dat",flags->runDir);
    FILE *waveFile = fopen(filename,"w");
    print_waveform_strain(data,model[chain->index[0]],waveFile);
    fclose(waveFile);


    sprintf(filename,"%s/avg_log_likelihood.dat",flags->runDir);
    FILE *chainFile = fopen(filename,"w");
    for(int ic=0; ic<NC; ic++) fprintf(chainFile,"%lg %lg\n",1./chain->temperature[ic],chain->avgLogL[ic]/(double)(flags->NMCMC/data->downsample));
    fclose(chainFile);
    
    //print total run time
    stop = time(NULL);
    
    printf(" ELAPSED TIME = %g seconds on %i thread(s)\n",(double)(stop-start),numThreads);
    sprintf(filename,"%s/ucb_mcmc.log",flags->runDir);
    FILE *runlog = fopen(filename,"a");
    fprintf(runlog," ELAPSED TIME = %g seconds on %i thread(s)\n",(double)(stop-start),numThreads);
    fclose(runlog);
        
    return 0;
}

