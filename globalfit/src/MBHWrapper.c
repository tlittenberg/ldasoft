//
//  MBHWrapper.c
//  global_fit
//
//  Created by Tyson Littenberg on 9/7/21.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

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
        
        char filename[MAXSTRINGSIZE];
        sprintf(filename,"%s/search_sources.dat",data->searchDir);
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
            fscanf(searchFile,"%lg %lg",&x,&x);
            for(int j=0; j<NParams; j++) fscanf(searchFile,"%lg",&searchParams[i][j]);
            
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
    
    memcpy(mbh_data->flags, gbmcmc_data->flags, sizeof(struct Flags));
}

static void set_mbh_priors(struct MBH_Data *dat, double *min, double *max)
{
    max[0] = log(0.44*5.0e8);
    max[1] = log(5.0e8);
    min[0] = log(1.0e2);
    min[1] = log(1.0e3);
    
    
    max[2] = 0.999;
    max[3] = 0.999;
    max[4] = PI;
    max[5] = 2.0*dat->Tend;
    max[6] = log(1.0e3);
    max[7] = 1.0;
    max[8] = 2.0*PI;
    max[9] = PI;
    max[10] = 1.0;
    
    
    min[2] = -0.999;
    min[3] = -0.999;
    min[4] = 0.0;
    min[5] = 1.01*dat->Tstart;
    min[6] = log(0.1);
    min[7] = -1.0;
    min[8] = 0.0;
    min[9] = 0.0;
    min[10] = -1.0;
}

void setup_mbh_data(struct MBHData *mbh_data, struct GBMCMCData *gbmcmc_data, struct TDI *tdi_full, int procID)
{
    int NC = 24; //TODO: Find somewhere to put NC

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
    mbh_data->data->data[0][0] = tdi_full->A[0];
    mbh_data->data->data[1][0] = tdi_full->E[0];
    mbh_data->data->data[0][1] = 0.0;
    mbh_data->data->data[1][1] = 0.0;
    for(int i=1; i<mbh_data->data->N/2; i++)
    {
        int j = i;
        int k = mbh_data->data->N-i;
        
        /* data to be used by mbh sampler */
        mbh_data->data->data[0][j] = tdi_full->A[2*i];   //re
        mbh_data->data->data[0][k] = tdi_full->A[2*i+1]; //im
        mbh_data->data->data[1][j] = tdi_full->E[2*i];   //re
        mbh_data->data->data[1][k] = tdi_full->E[2*i+1]; //im
    
    }
    
    /* allocate the memory needed for the MBH update() function */
    mbh_data->NH = 1000;
    mbh_data->heat = double_vector(NC);
    mbh_data->logLx = double_vector(NC);
    mbh_data->paramx = double_matrix(NC,NParams);
    mbh_data->paramy = double_matrix(NC,NParams);
    mbh_data->history = double_tensor(NC,mbh_data->NH,NParams);
    mbh_data->ejump = double_matrix(NC,NParams);
    mbh_data->evec = double_tensor(NC,NParams,NParams);
    mbh_data->sx = double_matrix(NC,mbh_data->data->Nch);
    mbh_data->sy = double_matrix(NC,mbh_data->data->Nch);
    mbh_data->max = (double*)malloc(sizeof(double)* (NParams));
    mbh_data->min = (double*)malloc(sizeof(double)* (NParams));
    mbh_data->who = int_vector(NC);
    mbh_data->av = int_matrix(5,NC);
    mbh_data->cv = int_matrix(5,NC);
    
    for (int i=0; i< NC; i++) mbh_data->who[i] = i;

    
    /* Set up RNG for MBH sampler */
    const gsl_rng_type * P;
    gsl_rng_env_setup();
    P = gsl_rng_default;

    mbh_data->rvec = (gsl_rng **)malloc(sizeof(gsl_rng *) * (NC+1));

    for(int i = 0 ; i<= NC; i++){
        mbh_data->rvec[i] = gsl_rng_alloc(P);
        gsl_rng_set(mbh_data->rvec[i] , i);
    }

    set_mbh_priors(mbh_data->data,mbh_data->min,mbh_data->max);

    //initialize PSD model
    for(int i=0; i<mbh_data->data->N/2; i++)
    {
        double f = (double)i/mbh_data->data->Tobs;
        mbh_data->data->SM[0][i] = AEnoise_FF(Larm, fstar, f)/sqrt(2.);
        mbh_data->data->SM[1][i] = AEnoise_FF(Larm, fstar, f)/sqrt(2.);
        mbh_data->data->SN[0][i] = AEnoise_FF(Larm, fstar, f)/sqrt(2.);
        mbh_data->data->SN[1][i] = AEnoise_FF(Larm, fstar, f)/sqrt(2.);
    }
    
    //load reference parameters
    //TODO: ACTUALLY READ FILES!
    //1.288603445394820e+03
    //7.820288034997643e+05
    mbh_data->paramx[0][0]=1.012977502365184e+06;
    mbh_data->paramx[0][1]=7.988754990395494e+05;
    mbh_data->paramx[0][2]=5.864208707957034e-01;
    mbh_data->paramx[0][3]=3.858559251489717e-01;
    mbh_data->paramx[0][4]=2.299319320412110e+00;
    mbh_data->paramx[0][5]=4.799245290472143e+06;
    mbh_data->paramx[0][6]=1.736720156296295e+01;
    mbh_data->paramx[0][7]=4.718350322394367e-01;
    mbh_data->paramx[0][8]=4.514867444831509e+00;
    mbh_data->paramx[0][9]=2.046608869505739e-01;
    mbh_data->paramx[0][10]=2.401959040445284e-02;
    map_params(2, mbh_data->paramx[0]);

    for(int i=1; i<NC; i++) for(int j=0; j<NParams; j++) mbh_data->paramx[i][j] = mbh_data->paramx[0][j];

    //set up heterodyne likelihood
    het_space(mbh_data->data, mbh_data->het, 2, mbh_data->paramx[0], mbh_data->min, mbh_data->max);
    heterodyne(mbh_data->data, mbh_data->het, 2, mbh_data->paramx[0]);

    //initialize noise model
    for (int i=0; i< NC; i++)
    {
        for (int j=0; j< mbh_data->data->Nch; j++) mbh_data->sx[i][j] = 1.0;
    }
    
    //initialize parameter chain buffer
    for (int i=0; i< NC; i++)
    {
        for (int k=0; k< mbh_data->NH; k++)
        {
            for (int j=0; j< NParams; j++) mbh_data->history[i][k][j] = mbh_data->paramx[i][j];
        }
    }
    
    //initialize likelihood for each chain
    for (int i=0; i< NC; i++)
    {
        mbh_data->logLx[i] = log_likelihood_het(mbh_data->data, mbh_data->het, 2, mbh_data->paramx[i], mbh_data->sx[i]);
    }
    

    /* set up temperature ladder */
    
    // run NCC cold chains
    for (int i=0; i< NCC; i++) mbh_data->heat[i] = 1.0;
    
    
    double SNR = mbh_data->het->SNR;
    double x = pow((SNR/5.0),1.0/(double)(NC-NCC));
    if(x > 1.3) x = 1.3;
    for (int i=NCC; i< NC; i++) mbh_data->heat[i] = mbh_data->heat[i-1]*x;
    
    
    /* check that data got filled correctly */
    if(procID>=mbh_data->procID_min && procID<=mbh_data->procID_max)
    {
        char filename[1024];
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
    }
}
