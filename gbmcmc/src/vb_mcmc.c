
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

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <omp.h>

#include <LISA.h>

#include "GalacticBinary.h"
#include "GalacticBinaryIO.h"
#include "GalacticBinaryData.h"
#include "GalacticBinaryMath.h"
#include "GalacticBinaryPrior.h"
#include "GalacticBinaryModel.h"
#include "GalacticBinaryProposal.h"
#include "GalacticBinaryWaveform.h"
#include "GalacticBinaryCatalog.h"
#include "GalacticBinaryMCMC.h"

/**
 * vbmcmc prototypes and static functions
 *
 */

void initialize_vb_proposal(struct Orbit *orbit, struct Data *data, struct Prior *prior, struct Chain *chain, struct Flags *flags, struct Proposal **proposal, int NMAX);
void GetVerificationBinary(struct Data *data, struct Flags *flags, FILE *vbFile);


static void **copy_argv(int argc, char **argv, char **new_argv)
{
    for(int i = 0; i < argc; ++i)
    {
        size_t length = strlen(argv[i])+1;
        new_argv[i] = malloc(length);
        memcpy(new_argv[i], argv[i], length);
    }
    new_argv[argc] = NULL;
}

static void parse_vb_list(int argc, char **argv, struct Flags *flags)
{
    
    int vb_list_flag = 0;
    
    static struct option long_options[] =
    {
        /* These options set a flag. */
        {"known-sources", required_argument, 0, 0},
        {0, 0, 0, 0}
    };
    
    opterr = 0;
    int opt=0;
    int long_index=0;
    
    //copy argv since getopt permutes order
    char **argv_copy=malloc((argc+1) * sizeof *argv_copy);
    copy_argv(argc,argv,argv_copy);

    
    //Loop through argv string and find argument for verification binaries
    while ((opt = getopt_long_only(argc, argv_copy,"apl:b:", long_options, &long_index )) != -1)
    {
        
        switch (opt)
        {
            case 0:
                if(strcmp("known-sources", long_options[long_index].name) == 0)
                {
                    strcpy(flags->vbFile,optarg);
                    vb_list_flag=1;
                }

                break;
            default:
                break;
                //print_usage();
                //exit(EXIT_FAILURE);
        }
    }
    
    if(!vb_list_flag)
    {
        fprintf(stdout,"Verification binary list required\n");
        fprintf(stdout, "Example ./vb_mcmc --known-sources /path/to/full_list.txt\n");
        exit(1);
    }

    //count lines in the file (one source per line)
    flags->NVB=0;    //repurpose NINJ flag for number of VBs
    char *line;
    char buffer[MAXSTRINGSIZE];
    
    FILE *sourceFile = fopen(flags->vbFile,"r");
    while( (line=fgets(buffer, MAXSTRINGSIZE, sourceFile)) != NULL) flags->NVB++;
    fclose(sourceFile);
    
    //reset opt counter
    optind = 0;
}

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
    parse(argc,argv,data,orbit,flags,chain,NMAX,0,0);
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

    /* Initialize LISA orbit model */
    initialize_orbit(data, orbit, flags);

    
    /* initialize all data structures */
    for(int n=0; n<flags->NVB; n++)
    {
        /* Setup output directories for data and chain files */
        sprintf(data_vec[n]->dataDir,"%s/data_%i",flags->runDir,n);
        sprintf(chain_vec[n]->chainDir,"%s/chains_%i",flags->runDir,n);

        chain=chain_vec[n];
        data=data_vec[n];
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
        GalacticBinaryReadData(data,orbit,flags);
        
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
    struct Prior *prior = NULL;
    struct Proposal **proposal = NULL;
    struct Model **trial = NULL;
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
//    if(flags->resume)
//    {
//        fprintf(stdout,"\n=============== Checkpointing ===============\n");
//
//        //check for files needed to resume
//        FILE *fptr = NULL;
//        int file_error = 0;
//
//        for(int ic=0; ic<chain->NC; ic++)
//        {
//            sprintf(filename,"%s/checkpoint/chain_state_%i.dat",flags->runDir,ic);
//
//            if( (fptr = fopen(filename,"r")) == NULL )
//            {
//                fprintf(stderr,"Warning: Could not checkpoint run state\n");
//                fprintf(stderr,"         Parameter file %s does not exist\n",filename);
//                file_error++;
//                break;
//            }
//        }
//
//        //if all of the files exist resume run from checkpointed state
//        if(!file_error)
//        {
//            fprintf(stdout,"   Checkpoint files found. Resuming chain\n");
//            restore_chain_state(orbit, data[FIXME], model, chain, flags, &mcmc_start);
//        }
//        fprintf(stdout,"============================================\n\n");
//    }
            
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
            //loop over frequency segments
            
            if(threadID==0)
            {
                flags->burnin   = (mcmc<0) ? 1 : 0;
                flags->maximize = (mcmc<-flags->NBURN/2) ? 1 : 0;
            }
            
            #pragma omp barrier
            // (parallel) loop over chains
            for(int ic=threadID; ic<NC; ic+=numThreads)
            {
                
                for(int n=0; n<flags->NVB; n++)
                {
                    /*
                    model = model_vec[n];
                    trial = trial_vec[n];
                    data = data_vec[n];
                    prior = prior_vec[n];
                    proposal = proposal_vec[n];
                    chain = chain_vec[n];
                     */
                    
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
                    trial = trial_vec[n];
                    data = data_vec[n];
                    prior = prior_vec[n];
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
    sprintf(filename,"%s/vb_mcmc.log",flags->runDir);
    FILE *runlog = fopen(filename,"a");
    fprintf(runlog," ELAPSED TIME = %g seconds on %i thread(s)\n",(double)(stop-start),numThreads);
    fclose(runlog);
    
    
    return 0;
}

void initialize_vb_proposal(struct Orbit *orbit, struct Data *data, struct Prior *prior, struct Chain *chain, struct Flags *flags, struct Proposal **proposal, int NMAX)
{
    int NC = chain->NC;
    double check  =0.0;
    double rjcheck=0.0;
    
    for(int i=0; i<chain->NP; i++)
    {
        proposal[i] = malloc(sizeof(struct Proposal));

        proposal[i]->trial  = malloc(NC*sizeof(int));
        proposal[i]->accept = malloc(NC*sizeof(int));
        
        for(int ic=0; ic<NC; ic++)
        {
            proposal[i]->trial[ic]  = 1;
            proposal[i]->accept[ic] = 0;
        }
        
        switch(i)
        {
            case 0:
                sprintf(proposal[i]->name,"prior");
                proposal[i]->function = &draw_from_uniform_prior;
                proposal[i]->density  = &prior_density;
                proposal[i]->weight   = 0.0;
                proposal[i]->rjweight = 0.0;
                setup_prior_proposal(flags, prior, proposal[i]);
                check   += proposal[i]->weight;
                rjcheck += proposal[i]->rjweight;
                break;
            case 1:
                sprintf(proposal[i]->name,"fstat draw");
                //setup_fstatistic_proposal(orbit, data, flags, proposal[i]);
                proposal[i]->function = &draw_from_fstatistic;
                proposal[i]->density  = &evaluate_fstatistic_proposal;
                proposal[i]->weight   = 0.0;
                proposal[i]->rjweight = 0.0; //that's a 1 all right.  don't panic
                check += proposal[i]->weight;
                break;
            case 2:
                sprintf(proposal[i]->name,"fstat jump");
                
                //re-use setup of fstat proposal from case 0
                proposal[i]->size   = proposal[1]->size;
                proposal[i]->norm   = proposal[1]->norm;
                proposal[i]->maxp   = proposal[1]->maxp;
                proposal[i]->vector = proposal[1]->vector;
                proposal[i]->matrix = proposal[1]->matrix;
                proposal[i]->tensor = proposal[1]->tensor;
                
                proposal[i]->function = &jump_from_fstatistic;
                proposal[i]->density  = &evaluate_fstatistic_proposal;
                proposal[i]->weight   = 0.0;
                proposal[i]->rjweight = 0.0;
                check   += proposal[i]->weight;
                rjcheck += proposal[i]->rjweight;
                break;
            case 3:
                sprintf(proposal[i]->name,"extrinsic prior");
                proposal[i]->function = &draw_from_extrinsic_prior;
                proposal[i]->density  = &prior_density;
                proposal[i]->weight   = 0.0;
                proposal[i]->rjweight = 0.0;
                check   += proposal[i]->weight;
                rjcheck += proposal[i]->rjweight;
                break;
            case 4:
                sprintf(proposal[i]->name,"fisher");
                proposal[i]->function = &draw_from_fisher;
                proposal[i]->density  = &symmetric_density;
                proposal[i]->weight = 1.0; //that's a 1 all right.  don't panic
                proposal[i]->rjweight = 0.0;
                //check   += proposal[i]->weight;
                rjcheck += proposal[i]->rjweight;
                break;
            case 5:
                sprintf(proposal[i]->name,"fm shift");
                proposal[i]->function = &fm_shift;
                proposal[i]->density  = &symmetric_density;
                proposal[i]->weight   = 0.0;
                proposal[i]->rjweight = 0.0;
                check   += proposal[i]->weight;
                rjcheck += proposal[i]->rjweight;
                break;
            case 6:
                sprintf(proposal[i]->name,"psi-phi jump");
                proposal[i]->function = &psi_phi_jump;
                proposal[i]->density  = &symmetric_density;
                proposal[i]->weight   = 0.2;
                proposal[i]->rjweight = 0.0;
                check   += proposal[i]->weight;
                rjcheck += proposal[i]->rjweight;
                break;
            case 7:
                sprintf(proposal[i]->name,"gmm draw");
                proposal[i]->function = &draw_from_gmm_prior;
                proposal[i]->density = &gmm_prior_density;
                proposal[i]->weight   = 0.0;
                proposal[i]->rjweight = 0.0;
                if(flags->catalog)
                {
                    setup_gmm_proposal(data, proposal[i]);
                    proposal[i]->weight   = 0.0;
                    proposal[i]->rjweight = 0.0;
                }
                check   += proposal[i]->weight;
                rjcheck += proposal[i]->rjweight;
                break;
            case 8:
                sprintf(proposal[i]->name,"cov draw");
                proposal[i]->function = &draw_from_cov;
                proposal[i]->density  = &cov_density;
                proposal[i]->weight  = 0.0;
                proposal[i]->rjweight = 0.0;
                if(flags->updateCov)
                {
                    setup_covariance_proposal(data, flags, proposal[i]);
                    proposal[i]->weight   = 0.0;
                    proposal[i]->rjweight = 0.0;
                }
                check   += proposal[i]->weight;
                rjcheck += proposal[i]->rjweight;
                break;
            default:
                break;
        }
    }
    //Fisher proposal fills in the cracks for fixed D moves
    proposal[4]->weight -= check;
    
    //Fstat proposal fills in the cracks for trans D moves
    proposal[1]->rjweight -= rjcheck;
    
    if(proposal[4]->weight<0.0 || proposal[1]->rjweight < 0.0)
    {
        fprintf(stderr,"Proposal weights not normalized (line %d of file %s)\n",__LINE__,__FILE__);
        exit(1);
    }
    
    if(!flags->quiet)
    {
        fprintf(stdout,"\n============== Proposal Cocktail ==============\n");
        fprintf(stdout,"   MCMC proposals:\n");
        for(int i=0; i<chain->NP; i++)
        {
            if(proposal[i]->weight>0.0)fprintf(stdout,"     %i) %s %lg\n",i,proposal[i]->name,proposal[i]->weight);
        }
        fprintf(stdout,"   RJMCMC proposals:\n");
        for(int i=0; i<chain->NP; i++)
        {
            if(proposal[i]->rjweight)fprintf(stdout,"     %i) %s %lg\n",i,proposal[i]->name,proposal[i]->rjweight);
        }
        fprintf(stdout,"===============================================\n");
    }
}

void GetVerificationBinary(struct Data *data, struct Flags *flags, FILE *vbFile)
{
    /* Get injection parameters */
    double f0,dfdt,costheta,phi,m1,m2,D; //read from injection file
    double cosi,phi0,psi;                //drawn from prior
    double Mc,amp;                       //calculated
    
    int check = fscanf(vbFile,"%lg %lg %lg %lg %lg %lg %lg %lg",&f0,&dfdt,&costheta,&phi,&m1,&m2,&cosi,&D);
    if(!check)
    {
        fprintf(stderr,"Error reading %s\n",flags->vbFile);
        exit(1);
    }
    
    //incoming distance in kpc, function expects pc
    D *= 1000.0;
    
    //compute derived parameters
    Mc  = chirpmass(m1,m2);
    amp = galactic_binary_Amp(Mc, f0, D);
    
    //initialize extrinsic parameters
    phi0 = 0.0;
    psi  = 0.0;
    
    struct Source *inj = data->inj;
    
    //set bandwidth of data segment centered on injection
    data->fmin = f0 - (data->N/2)/data->T;
    data->fmax = f0 + (data->N/2)/data->T;
    data->qmin = (int)(data->fmin*data->T);
    data->qmax = data->qmin+data->N;
    
    //recompute fmin and fmax so they align with a bin
    data->fmin = data->qmin/data->T;
    data->fmax = data->qmax/data->T;
    
    //map parameters to vector
    inj->f0       = f0;
    inj->dfdt     = dfdt;
    inj->costheta = costheta;
    inj->phi      = phi;
    inj->amp      = amp;
    inj->cosi     = cosi;
    inj->phi0     = phi0;
    inj->psi      = psi;
    map_params_to_array(inj, inj->params, data->T);
}


