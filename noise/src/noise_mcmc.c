//
//  noise_mcmc.c
//
//
//  Created by Tyson Littenberg on 4/06/21.
//

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
#include <gbmcmc.h>
#include <Noise.h>

int main(int argc, char *argv[])
{
    time_t start, stop;
    start = time(NULL);
    char filename[MAXSTRINGSIZE];
    
    /* check arguments */
    print_LISA_ASCII_art(stdout);
    print_version(stdout);
    if(argc==1) print_usage();

    /* Allocate data structures */
    struct Data *data   = malloc(sizeof(struct Data));
    struct Flags *flags = malloc(sizeof(struct Flags));
    struct Orbit *orbit = malloc(sizeof(struct Orbit));
    struct Chain *chain = malloc(sizeof(struct Chain));
    
    parse(argc,argv,data,orbit,flags,chain,1);
        
    /* Setup output directories for data and chain files */
    sprintf(data->dataDir,"%s/data",flags->runDir);
    sprintf(chain->chainDir,"%s/chains",flags->runDir);
    sprintf(chain->chkptDir,"%s/checkpoint",flags->runDir);
    
    mkdir(flags->runDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(data->dataDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain->chainDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain->chkptDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    /* Initialize data structures */
    alloc_data(data, flags);
    
    /* Initialize LISA orbit model */
    initialize_orbit(data, orbit, flags);
    
    /* Initialize chain structure and files */
    initialize_chain(chain, flags, &data->cseed, "a");
    
    /* read data */
    if(flags->strainData)
        GalacticBinaryReadData(data,orbit,flags);
    else if (flags->simNoise)
        GalacticBinarySimulateData(data, orbit, flags);
    
    /* Initialize Instrument Noise Model */
    printf("   ...initialize instrument noise model\n");
    struct InstrumentModel **model = malloc(chain->NC*sizeof(struct InstrumentModel *));
    for(int ic=0; ic<chain->NC; ic++)
    {
        model[ic] = malloc(sizeof(struct InstrumentModel));
        initialize_instrument_model(orbit, data, model[ic]);
    }
    sprintf(filename,"%s/instrument_noise_model.dat",data->dataDir);
    print_noise_model(model[0]->psd, filename);


    /* Initialize Galactic Foreground Model */
    if(flags->confNoise)printf("   ...initialize foreground noise model\n");
    struct ForegroundModel **galaxy = malloc(chain->NC*sizeof(struct ForegroundModel *));
    for(int ic=0; ic<chain->NC; ic++)
    {
        galaxy[ic] = malloc(sizeof(struct ForegroundModel));
        if(flags->confNoise) initialize_foreground_model(orbit, data, galaxy[ic]);
    }
    if(flags->confNoise)
    {
        sprintf(filename,"%s/foreground_noise_model.dat",data->dataDir);
        print_noise_model(galaxy[0]->psd, filename);
    }

    /* Combine noise components to form covariance matrix */
    for(int ic=0; ic<chain->NC; ic++)
        if(flags->confNoise) generate_full_covariance_matrix(model[ic]->psd, galaxy[ic]->psd, data->Nchannel);

    /* get initial likelihood */
    for(int ic=0; ic<chain->NC; ic++)
    {
        invert_noise_covariance_matrix(model[ic]->psd);
        model[ic]->logL = noise_log_likelihood(data, model[ic]->psd);
    }

    sprintf(filename,"%s/full_noise_model.dat",data->dataDir);
    print_noise_model(model[0]->psd, filename);

    
    
    //MCMC
    printf("\n==== Noise MCMC Sampler ====\n");
    
    sprintf(filename,"%s/noise_chain.dat",chain->chainDir);
    FILE *noiseChainFile = fopen(filename,"w");

    FILE *foregroundChainFile = NULL;
    if(flags->confNoise)
    {
        sprintf(filename,"%s/foreground_chain.dat",chain->chainDir);
        foregroundChainFile = fopen(filename,"w");
    }
    
    int numThreads;
    int step = 0;
    int NC = chain->NC;
    
    #pragma omp parallel num_threads(flags->threads)
    {
        int threadID;
        
        //Save individual thread number
        threadID = omp_get_thread_num();
        
        //Only one thread runs this section
        if(threadID==0)  numThreads = omp_get_num_threads();
        
        #pragma omp barrier
        
        /* The MCMC loop */
        for(; step<flags->NMCMC;)
        {
            #pragma omp barrier
            
            // (parallel) loop over chains
            for(int ic=threadID; ic<NC; ic+=numThreads)
            {
                struct InstrumentModel *model_ptr = model[chain->index[ic]];
                struct ForegroundModel *galaxy_ptr = galaxy[chain->index[ic]];
                for(int mc=0; mc<10; mc++)
                {
                    noise_instrument_model_mcmc(orbit, data, model_ptr, galaxy_ptr, chain, flags, ic);
                    if(flags->confNoise) noise_foreground_model_mcmc(orbit, data, model_ptr, galaxy_ptr, chain, flags, ic);
                }
            }// end (parallel) loop over chains
            
            //Next section is single threaded. Every thread must get here before continuing
            
            #pragma omp barrier
            
            if(threadID==0)
            {
                noise_ptmcmc(model, chain, flags);
                
                if(step%(flags->NMCMC/10)==0)printf("noise_mcmc at step %i\n",step);

                print_instrument_state(model[chain->index[0]], noiseChainFile, step);
                if(flags->confNoise) print_foreground_state(galaxy[chain->index[0]], foregroundChainFile, step);

                if(step%(flags->NMCMC/10)==0)
                {
                    generate_instrument_noise_model(data,orbit,model[chain->index[0]]);
                    sprintf(filename,"%s/current_instrument_noise_model.dat",data->dataDir);
                    print_noise_model(model[chain->index[0]]->psd, filename);

                    if(flags->confNoise)
                    {
                        generate_galactic_foreground_model(data,orbit,galaxy[chain->index[0]]);
                        sprintf(filename,"%s/current_foreground_noise_model.dat",data->dataDir);
                        print_noise_model(galaxy[chain->index[0]]->psd, filename);
                    }
                }
                
                if(step%data->downsample==0 && step/data->downsample < data->Nwave)
                {
                    generate_instrument_noise_model(data,orbit,model[chain->index[0]]);
                    if(flags->confNoise)
                    {
                        generate_galactic_foreground_model(data,orbit,galaxy[chain->index[0]]);
                        generate_full_covariance_matrix(model[chain->index[0]]->psd,galaxy[chain->index[0]]->psd, data->Nchannel);
                    }
                    
                    for(int n=0; n<data->N; n++)
                        for(int i=0; i<data->Nchannel; i++)
                            data->S_pow[n][i][0][step/data->downsample] = model[chain->index[0]]->psd->C[i][i][n];
                }

                step++;
                
                
            }
            //Can't continue MCMC until single thread is finished
            #pragma omp barrier
            
        }// end of MCMC loop
        
    }// End of parallelization
    
    fclose(noiseChainFile);
    if(flags->confNoise)fclose(foregroundChainFile);

    generate_instrument_noise_model(data,orbit,model[chain->index[0]]);
    sprintf(filename,"%s/final_instrument_noise_model.dat",data->dataDir);
    print_noise_model(model[chain->index[0]]->psd, filename);

    if(flags->confNoise)
    {
        generate_galactic_foreground_model(data,orbit,galaxy[chain->index[0]]);
        sprintf(filename,"%s/final_foreground_noise_model.dat",data->dataDir);
        print_noise_model(galaxy[chain->index[0]]->psd, filename);
        
        generate_full_covariance_matrix(model[chain->index[0]]->psd, galaxy[chain->index[0]]->psd, data->Nchannel);
        sprintf(filename,"%s/final_full_noise_model.dat",data->dataDir);
        print_noise_model(model[chain->index[0]]->psd, filename);
    }
    
    print_noise_reconstruction(data, flags);

    sprintf(filename,"%s/whitened_data.dat",data->dataDir);
    if(flags->confNoise)generate_full_covariance_matrix(model[chain->index[0]]->psd,galaxy[chain->index[0]]->psd, data->Nchannel);
    print_whitened_data(data, model[chain->index[0]]->psd, filename);

    for(int ic=0; ic<chain->NC; ic++) free_instrument_model(model[ic]);
    free(model);
    
    //print total run time
    stop = time(NULL);
    
    printf(" ELAPSED TIME = %g seconds\n",(double)(stop-start));
    
    
    return 0;
}

