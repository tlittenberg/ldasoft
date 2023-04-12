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
#include <GalacticBinary.h>
#include <GalacticBinaryIO.h>
#include <GalacticBinaryMath.h>
#include <GalacticBinaryData.h>
#include <GalacticBinaryModel.h>

#include "Noise.h"

int main(int argc, char *argv[])
{
    time_t start, stop;
    start = time(NULL);
    char filename[PATH_BUFSIZE];
    
    print_LISA_ASCII_art(stdout);
    
    struct Data *data   = malloc(sizeof(struct Data));
    struct Flags *flags = malloc(sizeof(struct Flags));
    struct Orbit *orbit = malloc(sizeof(struct Orbit));
    struct Chain *chain = malloc(sizeof(struct Chain));
    
    parse(argc,argv,data,orbit,flags,chain,1);
    
    /*
     * Get Data
     */
    
    /* Setup output directories for data and chain files */
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
    
    /* Initialize chain structure and files */
    initialize_chain(chain, flags, &data->cseed, "a");
    
    /* read data */
    if(flags->strainData)
        GalacticBinaryReadData(data,orbit,flags);
    else if (flags->simNoise)
        GalacticBinarySimulateData(data, orbit, flags);
    
    /*
     * Initialize Spline Model
     */
    int Nspline = 32+1;
    struct SplineModel **model = malloc(chain->NC*sizeof(struct SplineModel *));
    for(int ic=0; ic<chain->NC; ic++)
    {
        model[ic] = malloc(sizeof(struct SplineModel));
        initialize_spline_model(orbit, data, model[ic], Nspline);
    }
    
    pathprintf(filename,"%s/initial_spline_points.dat",data->dataDir);
    print_noise_model(model[0]->spline, filename);
    
    pathprintf(filename,"%s/interpolated_spline_points.dat",data->dataDir);
    print_noise_model(model[0]->psd, filename);
    
    
    //MCMC
    pathprintf(filename,"%s/chain_file.dat",chain->chainDir);
    FILE *chainFile = fopen(filename,"w");

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
        for(; step<100000;)
        {
            
            #pragma omp barrier
            
            // (parallel) loop over chains
            for(int ic=threadID; ic<NC; ic+=numThreads)
            {
                struct SplineModel *model_ptr = model[chain->index[ic]];
                for(int mc=0; mc<10; mc++)
                {
                    
                    if(gsl_rng_uniform(chain->r[ic])<0.9)
                        noise_spline_model_mcmc(orbit, data, model_ptr, chain, flags, ic);
                    else
                        noise_spline_model_rjmcmc(orbit, data, model_ptr, chain, flags, ic);
                }
            }// end (parallel) loop over chains
            
            //Next section is single threaded. Every thread must get here before continuing
            
            #pragma omp barrier
            
            if(threadID==0)
            {
                spline_ptmcmc(model, chain, flags);
                
                if(step%10000==0)printf("noise_mcmc at step %i\n",step);
                
                if(step%100==0)
                {
                    print_spline_state(model[chain->index[0]], chainFile, step);
                    
                    sprintf(filename,"%s/current_interpolated_spline_points.dat",data->dataDir);
                    print_noise_model(model[chain->index[0]]->psd, filename);
                    
                    sprintf(filename,"%s/current_spline_points.dat",data->dataDir);
                    print_noise_model(model[chain->index[0]]->spline, filename);


                }
                
                if(step%data->downsample==0 && step/data->downsample < data->Nwave)
                {
                    for(int n=0; n<data->N; n++)
                    {
                        data->S_pow[n][0][0][step/data->downsample] = model[chain->index[0]]->psd->SnA[n];
                        data->S_pow[n][1][0][step/data->downsample] = model[chain->index[0]]->psd->SnE[n];
                    }
                    

                }

                step++;
                
                
            }
            //Can't continue MCMC until single thread is finished
            #pragma omp barrier
            
        }// end of MCMC loop
        
    }// End of parallelization
    
    fclose(chainFile);
    
    pathprintf(filename,"%s/final_spline_points.dat",data->dataDir);
    print_noise_model(model[chain->index[0]]->spline, filename);
    
    pathprintf(filename,"%s/final_interpolated_spline_points.dat",data->dataDir);
    print_noise_model(model[chain->index[0]]->psd, filename);
    
    print_noise_reconstruction(data, flags);

    for(int ic=0; ic<chain->NC; ic++) free_spline_model(model[ic]);
    free(model);
    
    //print total run time
    stop = time(NULL);
    
    printf(" ELAPSED TIME = %g seconds\n",(double)(stop-start));
    
    
    return 0;
}

