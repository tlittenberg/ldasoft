//
//  NoiseWrapper.c
//  
//
//  Created by Tyson Littenberg on 2/5/21.
//

#include <mpi.h>

#include <glass_utils.h>
#include <glass_noise.h>
#include <glass_ucb.h>
#include <mbh.h>

#include "glass_ucb_wrapper.h"
#include "glass_vgb_wrapper.h"
#include "glass_mbh_wrapper.h"
#include "glass_noise_wrapper.h"

void alloc_noise_data(struct NoiseData *noise_data, struct UCBData *ucb_data, int procID, int nProc)
{
    noise_data->status = 0;
    noise_data->procID = procID;
    noise_data->nProc = nProc;
    noise_data->flags = malloc(sizeof(struct Flags));
    noise_data->chain = malloc(sizeof(struct Chain));
    noise_data->data  = malloc(sizeof(struct Data));
    
    noise_data->orbit = ucb_data->orbit;
    memcpy(noise_data->flags, ucb_data->flags, sizeof(struct Flags));

}



void select_noise_segment(struct Noise *psd_full, struct Data *data, struct Chain *chain, struct Model **model)
{
    double Tobs = data->T;
    int qstart = (int)(psd_full->f[0]*Tobs);
    int q0 = (int)(model[0]->noise->f[0]*Tobs);
    int dq = q0 - qstart;
    
    for(int n=0; n<chain->NC; n++)
    {
        for(int i=0; i<model[n]->noise->Nchannel; i++)
        {
            for(int j=0; j<model[n]->noise->Nchannel; j++)
            {
                
                memcpy(model[n]->noise->C[i][j], psd_full->C[i][j]+dq, data->N*sizeof(double));
                memcpy(model[n]->noise->invC[i][j], psd_full->invC[i][j]+dq, data->N*sizeof(double));
            }
        }
        memcpy(model[n]->noise->detC, psd_full->detC+dq, data->N*sizeof(double));

    }

}

void setup_noise_data(struct NoiseData *noise_data, struct UCBData *ucb_data, struct VGBData *vgb_data, struct MBHData *mbh_data, struct TDI *tdi_full, int procID)
{
    noise_data->data->downsample = ucb_data->data->downsample;
    noise_data->data->Nwave      = 100;
    
    noise_data->chain->NC      = ucb_data->chain->NC;
    noise_data->data->Nchannel = ucb_data->data->Nchannel;
    strcpy(noise_data->data->format,ucb_data->data->format);

    double T = ucb_data->data->T;
    
    noise_data->data->T = T;
    
    //set noise model to cover ucb segment
    noise_data->data->fmin = ucb_data->data->fmin;
    noise_data->data->fmax = ucb_data->data->fmax;
    
    //set all processes noise models based on max/min ucb segment
    MPI_Bcast(&noise_data->data->fmin, 1, MPI_DOUBLE, ucb_data->procID_min, MPI_COMM_WORLD);
    MPI_Bcast(&noise_data->data->fmax, 1, MPI_DOUBLE, ucb_data->procID_max, MPI_COMM_WORLD);

    //adjust noise model bandwidth to account for VBs
    for(int n=0; n<vgb_data->flags->NVB; n++)
    {
        noise_data->data->fmin = (vgb_data->data_vec[n]->fmin < noise_data->data->fmin ) ? vgb_data->data_vec[n]->fmin : noise_data->data->fmin;
        noise_data->data->fmax = (vgb_data->data_vec[n]->fmax > noise_data->data->fmax ) ? vgb_data->data_vec[n]->fmax : noise_data->data->fmax;
    }

    //pad noise model
    noise_data->data->fmin -= 1./T;
    noise_data->data->fmax += 1./T;

    //adjust noise model bandwidth to account for MBHs
    if(mbh_data->NMBH>0)
    {
        //set limits of noise model to cover both models
        noise_data->data->fmin = 2./T;//(mbh_data->data->fmin < noise_data->data->fmin ) ? mbh_data->data->fmin : noise_data->data->fmin;
        noise_data->data->fmax = (mbh_data->data->fmax > noise_data->data->fmax ) ? mbh_data->data->fmax : noise_data->data->fmax;

        //pad noise model even more (MBH bandwidth fluctuates)
    }
        
    noise_data->data->N = (int)((noise_data->data->fmax - noise_data->data->fmin)*T);
    
    alloc_data(noise_data->data, noise_data->flags);
    
    noise_data->psd = malloc(sizeof(struct Noise *)*ucb_data->chain->NC);
    noise_data->inst_model = malloc(sizeof(struct InstrumentModel*)*ucb_data->chain->NC);
    noise_data->inst_trial = malloc(sizeof(struct InstrumentModel*)*ucb_data->chain->NC);
    noise_data->conf_model = malloc(sizeof(struct ForegroundModel*)*ucb_data->chain->NC);
    noise_data->conf_trial = malloc(sizeof(struct ForegroundModel*)*ucb_data->chain->NC);

    //get max and min samples
    noise_data->data->qmin = (int)(noise_data->data->fmin*noise_data->data->T);
    noise_data->data->qmax = noise_data->data->qmin+noise_data->data->N;
    noise_data->data->fmax = (double)noise_data->data->qmax/T;
    
    //store max and min frequency in MBH structure
    mbh_data->data->fmin = noise_data->data->fmin;
    mbh_data->data->fmax = noise_data->data->fmax;
    
    select_frequency_segment(noise_data->data, tdi_full);
    
    /*
     Initialize measured time of model update.
     Used to determine number of steps relative to mbh model
     */
    noise_data->cpu_time = 1.0;

}


void initialize_noise_sampler(struct NoiseData *noise_data)
{
    /* Aliases to ucb structures */
    struct Flags *flags = noise_data->flags;
    struct Chain *chain = noise_data->chain;
    struct Data *data   = noise_data->data;
    
    //first check if file exists
    char filename[MAXSTRINGSIZE];
    sprintf(filename,"%s/current_spline_points.dat",data->dataDir);
    int check=0;
    FILE *test=NULL;
    if( (test=fopen(filename,"r")) )
    {
        fclose(test);
        check=1;
    }
    
    /* Initialize parallel chain & sampler state */
    if(flags->resume && check)
    {
        initialize_chain(chain, flags, &data->cseed, "a");
        resume_noise_state(noise_data);
    }
    else
    {
        initialize_chain(chain, flags, &data->cseed, "w");
        initialize_noise_state(noise_data);
    }
        
    /* Set sampler counter */
    noise_data->mcmc_step = -flags->NBURN;
    
    /* Store data segment in working directory */
    print_data(data, data->tdi, flags);

}

void initialize_noise_state(struct NoiseData *noise_data)
{
    /* Aliases to ucb structures */
    struct Orbit *orbit = noise_data->orbit;
    struct Chain *chain = noise_data->chain;
    struct Data *data   = noise_data->data;
    struct Flags *flags = noise_data->flags;
    struct Noise **psd = noise_data->psd;
    struct InstrumentModel **inst_model = noise_data->inst_model;
    struct InstrumentModel **inst_trial = noise_data->inst_trial;
    struct ForegroundModel **conf_model = noise_data->conf_model;
    struct ForegroundModel **conf_trial = noise_data->conf_trial;

    int NC = chain->NC;
    
    //populate spline model
    for(int ic=0; ic<NC; ic++)
    {
        /* Initialize work space for assembling full Noise model */
        psd[ic] = malloc(sizeof(struct Noise));
        alloc_noise(psd[ic], data->NFFT, data->Nchannel);

        /* Initialize Instrument Noise Model */
        inst_model[ic] = malloc(sizeof(struct InstrumentModel));
        inst_trial[ic] = malloc(sizeof(struct InstrumentModel));
        initialize_instrument_model(orbit, data, inst_model[ic]);
        initialize_instrument_model(orbit, data, inst_trial[ic]);

        /* Initialize Galactic Foreground Model */
        if(flags->confNoise) 
        {
            conf_model[ic] = malloc(sizeof(struct ForegroundModel));
            conf_trial[ic] = malloc(sizeof(struct ForegroundModel));
            initialize_foreground_model(orbit, data, conf_model[ic]);
            initialize_foreground_model(orbit, data, conf_trial[ic]);
        }

    }
    
    char filename[128];
    sprintf(filename,"%s/instrument_noise_model.dat",data->dataDir);
    print_noise_model(inst_model[0]->psd, filename);
    
    if(flags->confNoise)
    {
        sprintf(filename,"%s/foreground_noise_model.dat",data->dataDir);
        print_noise_model(conf_model[0]->psd, filename);
    }


}

void resume_noise_state(struct NoiseData *noise_data)
{
    /* Aliases to ucb structures */
    struct Orbit *orbit = noise_data->orbit;
    struct Chain *chain = noise_data->chain;
    struct Data *data   = noise_data->data;
    struct InstrumentModel **inst_model = noise_data->inst_model;
    struct InstrumentModel **inst_trial = noise_data->inst_trial;

    int NC = chain->NC;    
    
    //count lines in file
    char filename[MAXSTRINGSIZE];
    sprintf(filename,"%s/noise_chain.dat",data->dataDir);
    FILE *noiseFile = fopen(filename,"r");
    int i;
    double junk;
    int Nstep=0;
    while(!feof(noiseFile))
    {
        fscanf(noiseFile,"%i %lg",&i,&junk);//iteration and logL
        for(int n=0; n<inst_model[0]->Nlink; n++) fscanf(noiseFile,"%lg",&junk);//acceleration noise parameters
        for(int n=0; n<inst_model[0]->Nlink; n++) fscanf(noiseFile,"%lg",&junk);//OMS noise parameters
        Nstep++;
    }
    rewind(noiseFile);
    Nstep--;
    
    //initialize instrument noise model
    for(int ic=0; ic<NC; ic++)
    {
        inst_model[ic] = malloc(sizeof(struct InstrumentModel));
        inst_trial[ic] = malloc(sizeof(struct InstrumentModel));
        initialize_instrument_model(orbit, data, inst_model[ic]);
        initialize_instrument_model(orbit, data, inst_trial[ic]);
    }
    
    //set instrument model to stored values
    for(int n=0; n<Nstep; n++)
    {
        fscanf(noiseFile,"%i %lg",&i,&junk);//iteration and logL
        for(int n=0; n<inst_model[0]->Nlink; n++) fscanf(noiseFile,"%lg",&inst_model[0]->sacc[n]);//acceleration noise parameters
        for(int n=0; n<inst_model[0]->Nlink; n++) fscanf(noiseFile,"%lg",&inst_model[0]->soms[n]);//OMS noise parameters
    }
    fclose(noiseFile);

    generate_instrument_noise_model(data,orbit,inst_model[0]);
    invert_noise_covariance_matrix(inst_model[0]->psd);
    inst_model[0]->logL = noise_log_likelihood(data, inst_model[0]->psd);
    
    for(int ic=1; ic<NC; ic++) copy_instrument_model(inst_model[0], inst_model[ic]);
}

int update_noise_sampler(struct NoiseData *noise_data)
{
    clock_t start = clock();
    
    /* Aliases to ucb structures */
    struct Flags *flags = noise_data->flags;
    struct Orbit *orbit = noise_data->orbit;
    struct Chain *chain = noise_data->chain;
    struct Data *data   = noise_data->data;
    struct Noise **psd  = noise_data->psd;
    struct InstrumentModel **inst_model = noise_data->inst_model;
    struct InstrumentModel **inst_trial = noise_data->inst_trial;
    struct ForegroundModel **conf_model = noise_data->conf_model;
    struct ForegroundModel **conf_trial = noise_data->conf_trial;

    int NC = chain->NC;
    
    //For saving the number of threads actually given
    int numThreads;
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
            
            //loop over frequency segments
            struct Noise *psd_ptr = psd[chain->index[ic]];
            struct InstrumentModel *inst_model_ptr = inst_model[chain->index[ic]];
            struct InstrumentModel *inst_trial_ptr = inst_trial[chain->index[ic]];
            struct ForegroundModel *conf_model_ptr = conf_model[chain->index[ic]];
            struct ForegroundModel *conf_trial_ptr = conf_trial[chain->index[ic]];

            //update log likelihood (data may have changed)
            inst_model_ptr->logL = noise_log_likelihood(data, inst_model_ptr->psd);
            
            //evolve fixed dimension sampler
            for(int steps=0; steps<10; steps++)
            {
                noise_instrument_model_mcmc(orbit, data, inst_model_ptr, inst_trial_ptr, conf_model_ptr, psd_ptr, chain, flags, ic);
                if(flags->confNoise) noise_foreground_model_mcmc(orbit, data, inst_model_ptr, conf_model_ptr, conf_trial_ptr, psd_ptr, chain, flags, ic);
            }
            
            
        }// end (parallel) loop over chains
        
    }// End of parallelization
#pragma omp barrier
    
    noise_ptmcmc(inst_model,chain,flags);
    adapt_temperature_ladder(chain, noise_data->mcmc_step+flags->NBURN);
    
    //print chain files
    fprintf(chain->noiseFile[0],"%i %.12g ",noise_data->mcmc_step,inst_model[chain->index[0]]->logL);
    print_instrument_state(inst_model[chain->index[0]], chain->noiseFile[0]);
    fprintf(chain->noiseFile[0],"\n");

    if(flags->confNoise)
    {
        fprintf(chain->foregroundFile[0],"%i %.12g ",noise_data->mcmc_step, conf_model[chain->index[0]]->logL);
        print_foreground_state(conf_model[chain->index[0]], chain->foregroundFile[0]);
        fprintf(chain->foregroundFile[0],"\n");
    }

    //save point estimate of noise model
    int i = (noise_data->mcmc_step+flags->NBURN)%data->Nwave;
    generate_instrument_noise_model(data,orbit,inst_model[chain->index[0]]);
    if(flags->confNoise)
    {
        generate_galactic_foreground_model(data,orbit,conf_model[chain->index[0]]);
        generate_full_covariance_matrix(inst_model[chain->index[0]]->psd,conf_model[chain->index[0]]->psd, data->Nchannel);
    }
    invert_noise_covariance_matrix(inst_model[chain->index[0]]->psd);

    for(int n=0; n<data->N; n++)
    {
        for(int m=0; m<data->Nchannel; m++)
            data->S_pow[n][m][i] = inst_model[chain->index[0]]->psd->C[m][m][n];
    }
    
    noise_data->mcmc_step++;

    clock_t stop = clock();
    noise_data->cpu_time = (double)(stop-start);

    return 1;
}

void print_noise_state(struct NoiseData *noise_data, FILE *fptr, int counter)
{
    struct Chain *chain = noise_data->chain;
    int ic = chain->index[0];

    struct InstrumentModel *inst_model = noise_data->inst_model[ic];
    struct ForegroundModel *conf_model = noise_data->conf_model[ic];

    fprintf(fptr,"%i ",counter);
    print_instrument_state(inst_model, fptr);
    if(noise_data->flags->confNoise) print_foreground_state(conf_model, fptr);
    fprintf(fptr,"\n");
}
