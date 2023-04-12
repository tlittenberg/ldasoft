//
//  MBHWrapper.c
//  global_fit
//
//  Created by Tyson Littenberg on 9/7/21.
//

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <time.h>
#include <math.h>
#include <omp.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

#include <LISA.h>
#include <gbmcmc.h>
#include <mbh.h>

#include "GalacticBinaryWrapper.h"
#include "MBHWrapper.h"

void parse_mbh_args(int argc, char **argv, struct MBHData *data)
{
    static struct option long_options[] =
    {
        /* These options set a flag. */
        {"mbh-search-path", required_argument, 0, 0},
        {"duration", required_argument, 0, 0},
        {"start-time", required_argument, 0, 0},
        {0,0,0,0}
    };
    
    opterr = 0;
    int opt = 0;
    int long_index=0;

    //copy argv since getopt permutes order
    char **argv_copy=malloc((argc+1) * sizeof *argv_copy);
    copy_argv(argc,argv,argv_copy); //defined in GalacticBinaryIO.c

    //flags to check that all flags are set
    double Tobs=0.0;
    double Tstart=0.0;
    int pathFlag = 0;
    data->NMBH = 0;

    //Loop through argv string and find argument for mbh binaries
    while ((opt = getopt_long_only(argc, argv_copy,"apl:b:", long_options, &long_index )) != -1)
    {
        
        switch (opt)
        {
            case 0:
                if(strcmp("duration", long_options[long_index].name) == 0) Tobs = atof(optarg);
                if(strcmp("start-time", long_options[long_index].name) == 0) Tstart = atof(optarg);
                if(strcmp("mbh-search-path", long_options[long_index].name) == 0)
                {
                    strcpy(data->searchDir,optarg);
                    pathFlag = 1;
                }
                
                break;
            default:
                break;
        }
    }
        
    //reset opt counter
    optind = 0;
    
    /* Parse the search parameter file */
    if(pathFlag)
    {
        int NMBH_search = 0;
        
        char filename[PATH_BUFSIZE];
        pathprintf(filename,"%s/search_sources.dat",data->searchDir);
        FILE *searchFile = fopen(filename,"r");
        
        
        //determine the number of sources found in the search
        char *line;
        char buffer[MAXSTRINGSIZE];
        while( (line=fgets(buffer, MAXSTRINGSIZE, searchFile)) != NULL) NMBH_search++;

        rewind(searchFile);
        
        //store search parameters
        double **searchParams = malloc(NMBH_search*sizeof(double *));
        for(int i=0; i<NMBH_search; i++)
        {
            searchParams[i] = malloc(NParams*sizeof(double));
            
            //burn off first two columns (what are they?)
            double x;
            ufscanf(searchFile,"%lg %lg",&x,&x);
            for(int j=0; j<NParams; j++) ufscanf(searchFile,"%lg",&searchParams[i][j]);
            
            //merger time is parameter 5
            double t_merge = searchParams[i][5];
            if(t_merge > Tstart && t_merge < Tstart + Tobs) data->NMBH++;
        }
        rewind(searchFile);
        
        //store search parameters in segment
        data->segParams = malloc(data->NMBH*sizeof(double *));
        for(int i=0; i<data->NMBH; i++) data->segParams[i] = malloc(NParams*sizeof(double));
        
        int counter = 0;
        for(int i=0; i<NMBH_search; i++)
        {
            //merger time is parameter 5
            double t_merge = searchParams[i][5];
            if(t_merge > Tstart && t_merge < Tstart + Tobs)
            {
                memcpy(data->segParams[counter], searchParams[i], NParams*sizeof(double));
                counter++;
            }
        }
        
        for(int i=0; i<NMBH_search; i++) free(searchParams[i]);
        free(searchParams);

        fclose(searchFile);

    }

}

void alloc_mbh_data(struct MBHData *mbh_data, struct GBMCMCData *gbmcmc_data, int procID)
{
    mbh_data->status = 0;
    mbh_data->procID = procID;
    mbh_data->data = malloc(sizeof(struct MBH_Data));
    mbh_data->het = malloc(sizeof(struct Het));
    mbh_data->flags = malloc(sizeof(struct Flags));
    mbh_data->tdi = malloc(sizeof(struct TDI));
    
    memcpy(mbh_data->flags, gbmcmc_data->flags, sizeof(struct Flags));
}



void select_mbh_segment(struct MBHData *mbh_data, struct TDI *tdi_full)
{
    struct MBH_Data *data = mbh_data->data;
    /*
     fill data matrix with FD TDI data
      MBH sampler uses GSL FFT conventions for complex arrays
     */
    data->data[0][0] = tdi_full->A[0];
    data->data[1][0] = tdi_full->E[0];
    data->data[0][1] = 0.0;
    data->data[1][1] = 0.0;
    for(int i=1; i<data->N/2; i++)
    {
        int j = i;
        int k = data->N-i;
        
        /* data to be used by mbh sampler */
        data->data[0][j] = tdi_full->A[2*i];   //re
        data->data[0][k] = tdi_full->A[2*i+1]; //im
        data->data[1][j] = tdi_full->E[2*i];   //re
        data->data[1][k] = tdi_full->E[2*i+1]; //im
    }
}

void select_mbh_noise(struct MBHData *mbh_data, struct Noise *psd)
{
    double Tobs = mbh_data->data->Tobs;
    int nstart = (int)(psd->f[0]*Tobs);
    
    for(int i=0; i<nstart; i++)
    {
        mbh_data->data->SN[0][i] = 1.0;
        mbh_data->data->SN[1][i] = 1.0;
        mbh_data->data->SM[0][i] = 1.0;
        mbh_data->data->SM[1][i] = 1.0;
    }
    for(int i=0; i<psd->N; i++)
    {
        //TODO: Check this factor of 2!
        mbh_data->data->SN[0][i+nstart] = 2*psd->SnA[i];
        mbh_data->data->SN[1][i+nstart] = 2*psd->SnE[i];
        mbh_data->data->SM[0][i+nstart] = 2*psd->SnA[i];
        mbh_data->data->SM[1][i+nstart] = 2*psd->SnE[i];
    }
}

void setup_mbh_data(struct MBHData *mbh_data, struct GBMCMCData *gbmcmc_data, struct TDI *tdi_full, int procID)
{

    mbh_data->data->N = 2*tdi_full->N; //tdi->N is # of bins, mbh->N is # of time samples
    mbh_data->data->Nch = 2;
    mbh_data->data->Tobs = gbmcmc_data->data->T;
    mbh_data->data->dt = mbh_data->data->Tobs/(double)mbh_data->data->N;
    mbh_data->data->sqrtTobs = sqrt(mbh_data->data->Tobs);
    mbh_data->data->Tstart = gbmcmc_data->data->t0[0];
    mbh_data->data->Tend = mbh_data->data->Tstart + mbh_data->data->Tobs;
    
    //allocate data matrix
    mbh_data->data->data = malloc(mbh_data->data->Nch*sizeof(double *));
    mbh_data->data->SN = malloc(mbh_data->data->Nch*sizeof(double *));//double_matrix(dat->Nch,dat->N/2);
    mbh_data->data->SM = malloc(mbh_data->data->Nch*sizeof(double *));//double_matrix(dat->Nch,dat->N/2);

    for(int n=0; n<mbh_data->data->Nch; n++)
    {
        mbh_data->data->data[n] = malloc(mbh_data->data->N*sizeof(double));
        mbh_data->data->SN[n] = malloc(mbh_data->data->N/2 * sizeof(double));
        mbh_data->data->SM[n] = malloc(mbh_data->data->N/2 * sizeof(double));
    }
    
    /*
     fill data matrix with FD TDI data
      MBH sampler uses GSL FFT conventions for complex arrays
     */
    select_mbh_segment(mbh_data, tdi_full);

    
    /* set number of chains to fill available number of threads */
    mbh_data->NC = 24; //minimum chain count
    
    //Chains should be a multiple of threads for best usage of cores
    if(mbh_data->NC % mbh_data->flags->threads !=0){
        mbh_data->NC += mbh_data->flags->threads - (mbh_data->NC % mbh_data->flags->threads);
    }
    
    /* allocate the memory needed for the MBH update() function */
    mbh_data->NH = 1000;
    mbh_data->heat = double_vector(mbh_data->NC);
    mbh_data->logLx = double_vector(mbh_data->NC);
    mbh_data->paramx = double_matrix(mbh_data->NC,NParams);
    mbh_data->paramy = double_matrix(mbh_data->NC,NParams);
    mbh_data->history = double_tensor(mbh_data->NC,mbh_data->NH,NParams);
    mbh_data->ejump = double_matrix(mbh_data->NC,NParams);
    mbh_data->Fisher = double_tensor(mbh_data->NC,NParams,NParams);
    mbh_data->evec = double_tensor(mbh_data->NC,NParams,NParams);
    mbh_data->sx = double_matrix(mbh_data->NC,mbh_data->data->Nch);
    mbh_data->sy = double_matrix(mbh_data->NC,mbh_data->data->Nch);
    mbh_data->max = (double*)malloc(sizeof(double)* (NParams));
    mbh_data->min = (double*)malloc(sizeof(double)* (NParams));
    mbh_data->who = int_vector(mbh_data->NC);
    mbh_data->av = int_matrix(5,mbh_data->NC);
    mbh_data->cv = int_matrix(5,mbh_data->NC);
    
    for (int i=0; i< mbh_data->NC; i++) mbh_data->who[i] = i;

    
    /* Set up RNG for MBH sampler */
    const gsl_rng_type * P;
    gsl_rng_env_setup();
    P = gsl_rng_default;

    mbh_data->rvec = (gsl_rng **)malloc(sizeof(gsl_rng *) * (mbh_data->NC+1));

    for(int i = 0 ; i<= mbh_data->NC; i++){
        mbh_data->rvec[i] = gsl_rng_alloc(P);
        gsl_rng_set(mbh_data->rvec[i] , i);
    }

    set_mbh_priors(mbh_data->data,2,mbh_data->min,mbh_data->max);

    //initialize PSD model
    for(int i=0; i<mbh_data->data->N/2; i++)
    {
        double f = (double)i/mbh_data->data->Tobs;
        mbh_data->data->SM[0][i] = AEnoise_FF(Larm, fstar, f)/sqrt(2.);
        mbh_data->data->SM[1][i] = AEnoise_FF(Larm, fstar, f)/sqrt(2.);
        mbh_data->data->SN[0][i] = AEnoise_FF(Larm, fstar, f)/sqrt(2.);
        mbh_data->data->SN[1][i] = AEnoise_FF(Larm, fstar, f)/sqrt(2.);
    }
    
    //get max and min frequency extent of MBH models
    double fmin=1.;
    double fmax=0.;
    double fstart;
    double fstop;
    double *params=malloc(NParams*sizeof(double));
    for(int n=0; n<mbh_data->NMBH; n++)
    {
        for(int i=0; i<NParams; i++) params[i] = mbh_data->segParams[n][i];
        map_params(2, params);
        
        het_space(mbh_data->data, mbh_data->het, 2, params, mbh_data->min, mbh_data->max);
        
        fstart = mbh_data->het->MN/mbh_data->data->Tobs;
        fstop  = mbh_data->het->MM/mbh_data->data->Tobs;
        
        if(fstart < fmin) fmin = fstart;
        if(fstop  > fmax) fmax = fstop;
        
        //het_space allocates memory!
        free(mbh_data->het->fgrid);
        free(mbh_data->het->freq);

    }
    mbh_data->data->fmin = fmin;
    mbh_data->data->fmax = fmax;
        
    
    /*
     create storage for full TDI of waveform model
     */
    alloc_tdi(mbh_data->tdi, tdi_full->N, tdi_full->Nchannel);
    mbh_data->tdi->delta = tdi_full->delta;

    /*
     Initialize measured time of model update.
     Used to scale number of steps for other models
     */
    mbh_data->cpu_time = 1.0;

    free(params);
}

void initialize_mbh_sampler(struct MBHData *mbh_data)
{
    
    /* set up chain directory and file */
    char filename[PATH_BUFSIZE];
    pathprintf(mbh_data->chainDir,"%s",mbh_data->flags->runDir);
    pathprintf(filename,"%s/chain.dat",mbh_data->chainDir);
    
    if(mbh_data->flags->resume)
    {
        //first count samples in the chain file
        char* line;
        char lineBuffer[MAXSTRINGSIZE];
        mbh_data->chainFile = fopen(filename,"r");
        int MBHSTEPS = 0;
        while((line = fgets(lineBuffer, MAXSTRINGSIZE, mbh_data->chainFile)) != NULL) MBHSTEPS++;
        MBHSTEPS--;
        rewind(mbh_data->chainFile);
        
        //copy all but last (corrupted) line of chain.dat to new file
        char newfilename[PATH_BUFSIZE];
        pathprintf(newfilename,"%s/chain.dat.new",mbh_data->chainDir);
        FILE *newChainFile = fopen(newfilename,"w");
        for(int i=0; i<MBHSTEPS-1; i++)
        {
            line=fgets(lineBuffer,MAXSTRINGSIZE,mbh_data->chainFile);
            fprintf(newChainFile,"%s",line);
        }
        
        //parse last (complete) line of chain.dat to initialize sampler state
        char *column = strtok(line," ");
        double contents[18];
        for(int n=0; n<18; n++)
        {
            sscanf(column, "%lg", &contents[n]);
            column=strtok(NULL," ");
        }
        for(int i=0; i<NParams; i++) mbh_data->paramx[0][i] = contents[i+2];
        map_params(2, mbh_data->paramx[0]);
        
        
        fclose(mbh_data->chainFile);
        fclose(newChainFile);

        //replace chain.dat with chain.dat.new
        if(!rename(newfilename, filename)) fprintf(stdout,"%s is renamed %s successfully\n",newfilename,filename);
        else fprintf(stdout,"%s could not be renamed.\n",newfilename);
        fflush(stdout);

        //reopen chain.dat to append new samples
        pathprintf(filename,"%s/chain.dat",mbh_data->chainDir);
        mbh_data->chainFile = fopen(filename,"a");

    }
    else
    {
        mbh_data->chainFile = fopen(filename,"w");
        for(int i=0; i<NParams; i++) mbh_data->paramx[0][i] = mbh_data->segParams[mbh_data->procID - mbh_data->procID_min][i];
        map_params(2, mbh_data->paramx[0]);
    }
    for(int i=1; i<mbh_data->NC; i++) for(int j=0; j<NParams; j++) mbh_data->paramx[i][j] = mbh_data->paramx[0][j];

    het_space(mbh_data->data, mbh_data->het, 2, mbh_data->paramx[0], mbh_data->min, mbh_data->max);
    heterodyne(mbh_data->data, mbh_data->het, 2, mbh_data->paramx[0]);

    //overwrite stored global max and min f w/ local range for this MBH
    mbh_data->data->fmin = mbh_data->het->MN/mbh_data->data->Tobs;
    mbh_data->data->fmax = mbh_data->het->MM/mbh_data->data->Tobs;

    //initialize noise model
    for (int i=0; i< mbh_data->NC; i++)
    {
        for (int j=0; j< mbh_data->data->Nch; j++) mbh_data->sx[i][j] = 1.0;
    }
    
    //initialize parameter chain buffer
    for (int i=0; i< mbh_data->NC; i++)
    {
        for (int k=0; k< mbh_data->NH; k++)
        {
            for (int j=0; j< NParams; j++) mbh_data->history[i][k][j] = mbh_data->paramx[i][j];
        }
    }
    
    //initialize likelihood for each chain
    for (int i=0; i< mbh_data->NC; i++)
    {
        mbh_data->logLx[i] = log_likelihood_het(mbh_data->data, mbh_data->het, 2, mbh_data->paramx[i], mbh_data->sx[i]);
    }
    

    /* set up temperature ladder */
    
    // run NCC cold chains
    for (int i=0; i< NCC; i++) mbh_data->heat[i] = 1.0;
    
    
    double SNR = mbh_data->het->SNR;
    double x = pow((SNR/5.0),1.0/(double)(mbh_data->NC-NCC));
    if(x > 1.3) x = 1.3;
    for (int i=NCC; i< mbh_data->NC; i++) mbh_data->heat[i] = mbh_data->heat[i-1]*x;

    
    /* store data segment in ASCII format */
    /*
    FILE *tempFile;
        
    sprintf(filename,"%s/power_data_0.dat",mbh_data->flags->runDir);
    tempFile = fopen(filename,"w");
    for(int i=mbh_data->het->MN; i<mbh_data->het->MM; i++)
    {
        int re = i;
        int im = mbh_data->data->N-i;
        double f = (double)i/mbh_data->data->Tobs;
        fprintf(tempFile,"%lg %lg %lg\n",f,
                mbh_data->data->data[0][re]*mbh_data->data->data[0][re]+mbh_data->data->data[0][im]*mbh_data->data->data[0][im],
                mbh_data->data->data[1][re]*mbh_data->data->data[1][re]+mbh_data->data->data[1][im]*mbh_data->data->data[1][im]);
    }
    fclose(tempFile);
    */
    
    
    /* set up Fisher matrix proposal for each chain */
    FisherHet(mbh_data->data, mbh_data->het, 2, mbh_data->paramx[0], mbh_data->Fisher[0]);
    FisherEvec(mbh_data->Fisher[0], mbh_data->ejump[0], mbh_data->evec[0], NParams);
    efix(mbh_data->data, mbh_data->het, 1, 2, mbh_data->paramx[0], mbh_data->min,mbh_data-> max, mbh_data->ejump[0], mbh_data->evec[0], 1.0);

    for(int i=1; i<mbh_data->NC; i++)
    {
        memcpy(mbh_data->ejump[i],mbh_data->ejump[0],NParams*sizeof(double));
        for(int j=0; j<NParams; j++)
        {
            memcpy(mbh_data->Fisher[i][j],mbh_data->Fisher[0][j],NParams*sizeof(double));
            memcpy(mbh_data->evec[i][j],mbh_data->evec[0][j],NParams*sizeof(double));
        }
    }
    
    freehet(mbh_data->het);//TODO: Temporarily disable updates to het_space

    mbh_data->mcmc_step=0;
}

int update_mbh_sampler(struct MBHData *mbh_data)
{
    clock_t start = clock();

    //aliases to MBH_Data structure members
    struct MBH_Data *dat = mbh_data->data;
    struct Het *het = mbh_data->het;
    struct Flags *flags = mbh_data->flags;
    double *logLx = mbh_data->logLx;
    double **paramx = mbh_data->paramx;
    double **paramy = mbh_data->paramy;
    double **sx = mbh_data->sx;
    double **sy = mbh_data->sy;
    double *min = mbh_data->min;
    double *max = mbh_data->max;
    int *who = mbh_data->who;
    double *heat = mbh_data->heat;
    double ***history = mbh_data->history;
    double ***Fisher = mbh_data->Fisher;
    double **ejump = mbh_data->ejump;
    double ***evec = mbh_data->evec;
    int **cv = mbh_data->cv;
    int **av = mbh_data->av;
    int NH = mbh_data->NH;
    int NC = mbh_data->NC;
    gsl_rng **rvec = mbh_data->rvec;
    FILE *chain = mbh_data->chainFile;
    
    //some extra stuff
    int *m = calloc(NC,sizeof(int));
    
    //proposal schedule
    int typ;
    double a = 0.5;
    double b = 0.2;
    double c = 0.0;

    
    // use one of the cold chains to produce the reference waveform
    het_space(dat, het, 2, paramx[who[0]], min, max);
    heterodyne(dat, het, 2, paramx[who[0]]);

    
    //initialize likelihood for each chain
    for(int i=0; i< NC; i++)  logLx[i] = log_likelihood_het(dat, het, 2, paramx[i], sx[i]);

    //compute Fisher matrices for current parameters
    //#pragma omp parallel for
    for(int i=0; i<NC; i++)
    {
        FisherHet(dat, het, 2, paramx[i], Fisher[i]);
        FisherEvec(Fisher[i], ejump[i], evec[i], NParams);
        efix(dat, het, 1, 2, paramx[i], min, max, ejump[i], evec[i], 1.0);
    }

    double alpha = gsl_rng_uniform(rvec[0]);
    double beta;
    
    for(int cycle = 0; cycle<1000; cycle++)
    {
        // decide if we are doing a MCMC update of all the chains or a PT swap
        if((NC > 1) && (alpha < 0.2)) // chain swap
        {
            int hold; //hold on to current chain index
            alpha = (double)(NC-1)*gsl_rng_uniform(rvec[0]);
            int j = (int)(alpha);
            beta = exp((logLx[who[j]]-logLx[who[j+1]])/heat[j+1] - (logLx[who[j]]-logLx[who[j+1]])/heat[j]);
            alpha = gsl_rng_uniform(rvec[0]);
            if(beta > alpha)
            {
                hold = who[j];
                who[j] = who[j+1];
                who[j+1] = hold;
            }
        }
        else // MCMC update
        {
            for(int j = 0; j < NC; j++)  for(int i = 0; i < NParams; i++) paramy[j][i] = paramx[j][i];
            
            // all chains do the same type of update since some (especially type 2) are much slower than the others. Saves them waiting on others to finish
            alpha = gsl_rng_uniform(rvec[0]);
            
            if(alpha > a)      typ = 0;
            else if(alpha > b) typ = 1;
            else if(alpha > c) typ = 2;
            else               typ = 3;
            
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
                for(int ic=threadID; ic<NC; ic+=numThreads) update(dat, het, typ, ic, 2, logLx, paramx, paramy, sx, sy, min, max, who, heat, history, NH, ejump, evec, cv, av, rvec[ic]);
                #pragma omp barrier
            }// End of parallelization
        }
        
        // add to the history file
        for(int k=0; k < NC; k++)
        {
            int q = who[k];
            int i = m[k]%NH;
            // the history file is kept for each temperature
            for(int j=0; j<NParams; j++) history[k][i][j] = paramx[q][j];
            m[k]++;
        }

        if(cycle%100==0) print_mbh_chain_file(dat, het, who, paramx, logLx, sx, 2, 1, chain);

    }//end cycle
    
    free(m);
    freehet(mbh_data->het);
    
    clock_t stop = clock();
    mbh_data->cpu_time = (double)(stop-start);
    mbh_data->mcmc_step++;
    return 1;
}

void get_mbh_waveform(struct MBHData *mbh_data)
{
    int NF = mbh_data->het->MM - mbh_data->het->MN;
    int index = (int)(mbh_data->data->fmin*mbh_data->data->Tobs);

    /* temporary storage for phase and amplitude */
    double *f = malloc(NF*sizeof(double));
    double *A_amp = malloc(NF*sizeof(double));
    double *E_amp = malloc(NF*sizeof(double));
    double *A_phi = malloc(NF*sizeof(double));
    double *E_phi = malloc(NF*sizeof(double));

    for(int n=0; n<NF; n++) f[n] = mbh_data->data->fmin + (double)n/mbh_data->data->Tobs;

    /* returns phase and amplitude for A and E channels on input frequency grid */
    int *who = mbh_data->who;
    fullphaseamp(mbh_data->data, 2, NF, mbh_data->paramx[who[0]], f, A_amp, E_amp, A_phi, E_phi);

    /* insert TDI(f) into correct section of TDI structure */
    for(int n=0; n<2*mbh_data->tdi->N; n++)
    {
        mbh_data->tdi->A[n] = 0.0;
        mbh_data->tdi->E[n] = 0.0;
    }
    for(int n=0; n<NF; n++)
    {
        int re = 2*(n+index);
        int im = re+1;
        
        mbh_data->tdi->A[re] = A_amp[n]*cos(A_phi[n]);
        mbh_data->tdi->A[im] = A_amp[n]*sin(A_phi[n]);
        mbh_data->tdi->E[re] = E_amp[n]*cos(E_phi[n]);
        mbh_data->tdi->E[im] = E_amp[n]*sin(E_phi[n]);

    }
    

    free(f);
    free(A_amp);
    free(E_amp);
    free(A_phi);
    free(E_phi);

}

