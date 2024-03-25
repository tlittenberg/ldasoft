/*
 *  Copyright (C) 2023 Tyson B. Littenberg (MSFC-ST12)
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


#include "glass_utils.h"
#include "gitversion.h"

#define FIXME 0

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60


void printProgress (double percentage)
{
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    fprintf(stdout, "\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush (stdout);
}


void print_version(FILE *fptr)
{
    fprintf(fptr, "\n");
    fprintf(fptr, "=============== GLASS Version: ==============\n\n");
    //fprintf(fptr, "  Git remote origin: %s\n", GIT_URL);
    //fprintf(fptr, "  Git version: %s\n", GIT_VER);
    fprintf(fptr, "  Git commit: %s\n", GITVERSION);
    //fprintf(fptr, "  Git commit author: %s\n",GIT_AUTHOR);
    //fprintf(fptr, "  Git commit date: %s\n", GIT_DATE);
    fprintf(fptr, "\n=============================================\n\n");
}

void setup_run_directories(struct Flags *flags, struct Data *data, struct Chain *chain)
{
    
    sprintf(data->dataDir,"%s/data",flags->runDir);
    sprintf(chain->chainDir,"%s/chains",flags->runDir);
    sprintf(chain->chkptDir,"%s/checkpoint",flags->runDir);

    mkdir(flags->runDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(data->dataDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain->chainDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain->chkptDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

}

void initialize_orbit(struct Data *data, struct Orbit *orbit, struct Flags *flags)
{
    /* Load spacecraft ephemerides */
    switch(flags->orbit)
    {
        case 0:
            initialize_analytic_orbit(orbit);
            break;
        case 1:
            initialize_numeric_orbit(orbit);
            break;
        default:
            fprintf(stderr,"unsupported orbit type\n");
            exit(1);
            break;
    }
    
    /* set approximate f/fstar for segment */
    data->sine_f_on_fstar = sin((data->fmin + (data->fmax-data->fmin)/2.)/orbit->fstar);
}

//TODO: Move file pointers to Model instead of Chain structs
void initialize_chain(struct Chain *chain, struct Flags *flags, long *seed, const char *mode)
{
    int ic;
    int NC = chain->NC;
    char filename[MAXSTRINGSIZE];

    chain->index = calloc(NC,sizeof(int));
    chain->acceptance = calloc(NC,sizeof(double));
    chain->temperature = calloc(NC,sizeof(double));
    chain->avgLogL     = calloc(NC,sizeof(double));
    chain->dimension   = calloc(NC,sizeof(int *));
    for(ic=0; ic<NC; ic++)
    {
        chain->index[ic]=ic;
        chain->acceptance[ic] = 1.0;
        chain->temperature[ic] = pow(1.2,(double)ic);
        chain->avgLogL[ic] = 0.0;
        chain->dimension[ic] = calloc(flags->DMAX,sizeof(int));
        for(int id=0; id<flags->DMAX; id++) chain->dimension[ic][id] = 0;
    }
    //set hottest chain to ~infinite temperature
    if(NC>1) chain->temperature[NC-1] = 1e12;
    chain->logLmax = 0.0;
    
    chain->r = malloc(NC*sizeof(gsl_rng *));
    chain->T = malloc(NC*sizeof(const gsl_rng_type *));
    
    for(ic=0; ic<NC; ic++)
    {
        chain->T[ic] = gsl_rng_default;
        chain->r[ic] = gsl_rng_alloc(chain->T[ic]);
        gsl_rng_env_setup();
        gsl_rng_set (chain->r[ic], *seed);
        *seed = (long)gsl_rng_get(chain->r[ic]);
    }
    
    if(!flags->quiet)
    {
        sprintf(filename,"%s/log_likelihood_chain.dat",chain->chainDir);
        chain->likelihoodFile = fopen(filename,mode);
        
        sprintf(filename,"%s/temperature_chain.dat",chain->chainDir);
        chain->temperatureFile = fopen(filename,mode);
    }
    
    chain->chainFile = malloc(NC*sizeof(FILE *));
    sprintf(filename,"%s/model_chain.dat.0",chain->chainDir);
    chain->chainFile[0] = fopen(filename,mode);
    
    chain->parameterFile = malloc(NC*sizeof(FILE *));
    sprintf(filename,"%s/parameter_chain.dat.0",chain->chainDir);
    chain->parameterFile[0] = fopen(filename,mode);
    
    chain->dimensionFile = malloc(flags->DMAX*sizeof(FILE *));
    for(int i=0; i<flags->DMAX; i++)
    {
        /* only create these files when needed */
        chain->dimensionFile[i]=NULL;
    }
    
    chain->noiseFile = malloc(NC*sizeof(FILE *));
    sprintf(filename,"%s/noise_chain.dat.0",chain->chainDir);
    chain->noiseFile[0] = fopen(filename,mode);
    
    if(flags->confNoise)
    {
        chain->foregroundFile = malloc(NC*sizeof(FILE *));
        sprintf(filename,"%s/foreground_chain.dat.0",chain->chainDir);
        chain->foregroundFile[0] = fopen(filename,mode);
    }
    
    if(flags->calibration)
    {
        chain->calibrationFile = malloc(NC*sizeof(FILE *));
        sprintf(filename,"%s/calibration_chain.dat.0",chain->chainDir);
        chain->calibrationFile[0] = fopen(filename,mode);
    }
    
    if(flags->verbose)
    {
        for(ic=1; ic<NC; ic++)
        {
            sprintf(filename,"%s/parameter_chain.dat.%i",chain->chainDir,ic);
            chain->parameterFile[ic] = fopen(filename,mode);
            
            sprintf(filename,"%s/model_chain.dat.%i",chain->chainDir,ic);
            chain->chainFile[ic] = fopen(filename,mode);
            
            sprintf(filename,"%s/noise_chain.dat.%i",chain->chainDir,ic);
            chain->noiseFile[ic] = fopen(filename,mode);
        }
    }
}

void alloc_data(struct Data *data, struct Flags *flags)
{
    int NMCMC = flags->NMCMC;
        
    data->logN = log((double)(2*data->N*data->Nchannel));
    
    //data->inj = malloc(sizeof(struct Source));
    //alloc_source(data->inj,data->N,data->Nchannel,data->NP);
    
    data->tdi   = malloc(sizeof(struct TDI));
    data->raw   = malloc(sizeof(struct TDI));
    data->noise = malloc(sizeof(struct Noise));
            
    alloc_tdi(data->tdi, data->N, data->Nchannel);
    alloc_tdi(data->raw, data->N, data->Nchannel);
    alloc_noise(data->noise, data->N, data->Nchannel);
    
    //reconstructed signal model
    int i_re,i_im;
    data->h_rec = malloc(data->N*2*sizeof(double **));
    data->h_res = malloc(data->N*2*sizeof(double **));
    data->r_pow = malloc(data->N*sizeof(double **));
    data->h_pow = malloc(data->N*sizeof(double **));
    data->S_pow = malloc(data->N*sizeof(double **));
    
    //number of waveform samples to save
    data->Nwave=100;
    
    //downsampling rate of post-burn-in samples
    data->downsample = NMCMC/data->Nwave;
    
    for(int i=0; i<data->N; i++)
    {
        i_re = i*2;
        i_im = i_re+1;
        
        data->S_pow[i]    = malloc(data->Nchannel*sizeof(double *));
        data->h_pow[i]    = malloc(data->Nchannel*sizeof(double *));
        data->r_pow[i]    = malloc(data->Nchannel*sizeof(double *));
        data->h_rec[i_re] = malloc(data->Nchannel*sizeof(double *));
        data->h_rec[i_im] = malloc(data->Nchannel*sizeof(double *));
        data->h_res[i_re] = malloc(data->Nchannel*sizeof(double *));
        data->h_res[i_im] = malloc(data->Nchannel*sizeof(double *));
        for(int n=0; n<data->Nchannel; n++)
        {
            data->S_pow[i][n]    = calloc(data->Nwave,sizeof(double));
            data->h_pow[i][n]    = calloc(data->Nwave,sizeof(double));
            data->r_pow[i][n]    = calloc(data->Nwave,sizeof(double));
            data->h_rec[i_re][n] = calloc(data->Nwave,sizeof(double));
            data->h_rec[i_im][n] = calloc(data->Nwave,sizeof(double));
            data->h_res[i_re][n] = calloc(data->Nwave,sizeof(double));
            data->h_res[i_im][n] = calloc(data->Nwave,sizeof(double));
        }
    }
    
    //Spectrum proposal
    data->p = calloc(data->N,sizeof(double));

    //catalog of previously detected sources
    //data->catalog = malloc(sizeof(struct Catalog));

}

void alloc_noise(struct Noise *noise, int NFFT, int Nchannel)
{
    noise->N = NFFT;
    noise->Nchannel = Nchannel;
    
    noise->eta = calloc(Nchannel,sizeof(double));

    noise->f   = calloc(NFFT,sizeof(double));

    noise->C    = malloc(Nchannel*sizeof(double **));
    noise->invC = malloc(Nchannel*sizeof(double **));
    
    for(int i=0; i<Nchannel; i++)
    {
        noise->eta[i] = 1.0;
        noise->C[i]    = malloc(Nchannel*sizeof(double *));
        noise->invC[i] = malloc(Nchannel*sizeof(double *));
        
        for(int j=0; j<Nchannel; j++)
        {
            noise->C[i][j]    = calloc(NFFT,sizeof(double));
            noise->invC[i][j] = calloc(NFFT,sizeof(double));
        }
    }

    noise->detC = calloc(NFFT,sizeof(double));
    noise->transfer = calloc(NFFT,sizeof(double));
    
    int n;
    for(n=0; n<NFFT; n++)
    {
        for(int i=0; i<Nchannel; i++) noise->C[i][i][n] = 1.0;
        for(int i=0; i<Nchannel; i++)
        {
            for(int j=i+1; i<Nchannel; i++)
            {
                noise->C[i][j][n] = 0.0;
                noise->C[j][i][n] = 0.0;
            }
        }
        noise->transfer[n] = 1.0;
    }
}

void alloc_calibration(struct Calibration *calibration)
{
    calibration->dampA = 0.0;
    calibration->dampE = 0.0;
    calibration->dampX = 0.0;
    calibration->dphiA = 0.0;
    calibration->dphiE = 0.0;
    calibration->dphiX = 0.0;
    calibration->real_dphiA = 1.0;
    calibration->real_dphiE = 1.0;
    calibration->real_dphiX = 1.0;
    calibration->imag_dphiA = 0.0;
    calibration->imag_dphiE = 0.0;
    calibration->imag_dphiX = 0.0;
}

//TODO: Expand copy_data() to include everything, and then replace where needed (NoiseWrapper.c, ...)
void copy_data(struct Data *origin, struct Data *copy)
{
    memcpy(copy->format, origin->format, sizeof(origin->format));
    memcpy(copy->fileName, origin->fileName, sizeof(origin->fileName));
    copy->T=origin->T;
    copy->sqT=origin->sqT;
    copy->N=origin->N;
    copy->Nchannel=origin->Nchannel;
    copy->qpad=origin->qpad;
    copy->cseed=origin->cseed;
    copy->nseed=origin->nseed;
    copy->iseed=origin->iseed;
    copy->t0   = origin->t0;
}

void copy_noise(struct Noise *origin, struct Noise *copy)
{
    memcpy(copy->eta,origin->eta,origin->Nchannel*sizeof(double));

    memcpy(copy->f, origin->f, origin->N*sizeof(double));

    copy_Cij(origin->C, copy->C, origin->Nchannel, origin->N);
    copy_Cij(origin->invC, copy->invC, origin->Nchannel, origin->N);

    memcpy(copy->detC, origin->detC, origin->N*sizeof(double));
    memcpy(copy->transfer, origin->transfer, origin->N*sizeof(double));
}

void copy_Cij(double ***origin, double ***copy, int M, int N)
{
    for(int i=0; i<M; i++)
        for(int j=0; j<M; j++)
            memcpy(copy[i][j], origin[i][j], N*sizeof(double));
}

void copy_calibration(struct Calibration *origin, struct Calibration *copy)
{
    copy=origin;
    /*
    copy->dampA   = origin->dampA;
    copy->dampE   = origin->dampE;
    copy->dampX   = origin->dampX;
    copy->dphiA = origin->dphiA;
    copy->dphiE = origin->dphiE;
    copy->dphiX = origin->dphiX;
    copy->real_dphiA = origin->real_dphiA;
    copy->real_dphiE = origin->real_dphiE;
    copy->real_dphiX = origin->real_dphiX;
    copy->imag_dphiA = origin->imag_dphiA;
    copy->imag_dphiE = origin->imag_dphiE;
    copy->imag_dphiX = origin->imag_dphiX;
     */
}

void free_noise(struct Noise *noise)
{
    free(noise->eta);
    free(noise->f);
    for(int i=0; i<noise->Nchannel; i++)
    {
        for(int j=0; j<noise->Nchannel; j++)
        {
            free(noise->C[i][j]);
            free(noise->invC[i][j]);
        }
        free(noise->C[i]);
        free(noise->invC[i]);
    }
    free(noise->C);
    free(noise->invC);
    free(noise->detC);
    free(noise->transfer);
    free(noise);
}

void free_chain(struct Chain *chain, struct Flags *flags)
{
    free(chain->index);
    free(chain->acceptance);
    free(chain->temperature);
    free(chain->avgLogL);
    for(int ic=0; ic<chain->NC; ic++)
    {
        gsl_rng_free(chain->r[ic]);
        free(chain->dimension[ic]);
    }
    free(chain->dimension);
    free(chain->r);
    free(chain->T);
    
    if(!flags->quiet)
    {
        fclose(chain->likelihoodFile);
        fclose(chain->temperatureFile);
    }

    fclose(chain->chainFile[0]);

    fclose(chain->parameterFile[0]);
    
    for(int i=0; i<flags->DMAX; i++)
    {
        /* only create these files when needed */
        if(chain->dimensionFile[i]!=NULL) fclose(chain->dimensionFile[i]);
    }
    free(chain->dimensionFile);

    fclose(chain->noiseFile[0]);
    
    if(flags->calibration)
    {
        fclose(chain->calibrationFile[0]);
        free(chain->calibrationFile);
    }

    
    if(flags->verbose)
    {
        for(int ic=1; ic<chain->NC; ic++)
        {
            fclose(chain->chainFile[ic]);
            fclose(chain->parameterFile[ic]);
            fclose(chain->noiseFile[ic]);
        }
    }
    free(chain->chainFile);
    free(chain->parameterFile);
    free(chain->noiseFile);
    
    free(chain);
}

void free_calibration(struct Calibration *calibration)
{
    free(calibration);
}

void ReadHDF5(struct Data *data, struct TDI *tdi, struct Flags *flags)
{
    /* LDASOFT-formatted structure for TDI data */
    struct TDI *tdi_td = malloc(sizeof(struct TDI));
        
    if(!strcmp(data->format,"frequency"))  LISA_Read_HDF5_LDC_RADLER_TDI(tdi_td, data->fileName);
    if(!strcmp(data->format,"sangria")) LISA_Read_HDF5_LDC_TDI(tdi_td, data->fileName, "/obs/tdi");
    
    
    /* Select time segment of full data set */
    double start_time = data->t0;
    double stop_time = start_time + data->T;
    double dt = tdi_td->delta;
    double Tobs = stop_time - start_time;
    int N = (int)floor(Tobs/dt);

    /* work space for selecting and FT'ing time series */
    double *X = malloc(N*sizeof(double));
    double *Y = malloc(N*sizeof(double));
    double *Z = malloc(N*sizeof(double));
    double *A = malloc(N*sizeof(double));
    double *E = malloc(N*sizeof(double));
    double *T = malloc(N*sizeof(double));

    /* Allocate data->tdi structure for Fourier transform output */
    alloc_tdi(tdi, N/2, N_TDI_CHANNELS);
    tdi->delta = 1./Tobs;

    /* Select requested time segment */
    int n_start = (int)floor(start_time/dt); // first sample of time segment
    
    for(int n=0; n<N; n++)
    {
        int m = n_start+n;
        X[n] = tdi_td->X[m];
        Y[n] = tdi_td->Y[m];
        Z[n] = tdi_td->Z[m];
        A[n] = tdi_td->A[m];
        E[n] = tdi_td->E[m];
        T[n] = tdi_td->T[m];
    }
    
    /* lets get rid of those black holes */
    if(flags->no_mbh)
    {
        struct TDI *tdi_td_mbhb = malloc(sizeof(struct TDI));
        LISA_Read_HDF5_LDC_TDI(tdi_td_mbhb, data->fileName, "/sky/mbhb/tdi");
        for(int n=0; n<N; n++)
        {
            int m = n_start+n;
            X[n] -= tdi_td_mbhb->X[m];
            Y[n] -= tdi_td_mbhb->Y[m];
            Z[n] -= tdi_td_mbhb->Z[m];
            A[n] -= tdi_td_mbhb->A[m];
            E[n] -= tdi_td_mbhb->E[m];
            T[n] -= tdi_td_mbhb->T[m];
        }
        free_tdi(tdi_td_mbhb);
    }
    
    /* lets get rid of the galaxy */
    if(flags->no_ucb)
    {
        struct TDI *tdi_td_dgb = malloc(sizeof(struct TDI));
        LISA_Read_HDF5_LDC_TDI(tdi_td_dgb, data->fileName, "/sky/dgb/tdi");
        for(int n=0; n<N; n++)
        {
            int m = n_start+n;
            X[n] -= tdi_td_dgb->X[m];
            Y[n] -= tdi_td_dgb->Y[m];
            Z[n] -= tdi_td_dgb->Z[m];
            A[n] -= tdi_td_dgb->A[m];
            E[n] -= tdi_td_dgb->E[m];
            T[n] -= tdi_td_dgb->T[m];
        }
        free_tdi(tdi_td_dgb);
    
        struct TDI *tdi_td_igb = malloc(sizeof(struct TDI));
        LISA_Read_HDF5_LDC_TDI(tdi_td_igb, data->fileName, "/sky/igb/tdi");
        for(int n=0; n<N; n++)
        {
            int m = n_start+n;
            X[n] -= tdi_td_igb->X[m];
            Y[n] -= tdi_td_igb->Y[m];
            Z[n] -= tdi_td_igb->Z[m];
            A[n] -= tdi_td_igb->A[m];
            E[n] -= tdi_td_igb->E[m];
            T[n] -= tdi_td_igb->T[m];
        }
        free_tdi(tdi_td_igb);
    }

    /* lets get rid of the verification binaries */
    if(flags->no_vgb)
    {
        struct TDI *tdi_td_vgb = malloc(sizeof(struct TDI));
        LISA_Read_HDF5_LDC_TDI(tdi_td_vgb, data->fileName, "/sky/vgb/tdi");
        for(int n=0; n<N; n++)
        {
            int m = n_start+n;
            X[n] -= tdi_td_vgb->X[m];
            Y[n] -= tdi_td_vgb->Y[m];
            Z[n] -= tdi_td_vgb->Z[m];
            A[n] -= tdi_td_vgb->A[m];
            E[n] -= tdi_td_vgb->E[m];
            T[n] -= tdi_td_vgb->T[m];
        }
        free_tdi(tdi_td_vgb);
    }
    
    /* Tukey window time-domain TDI channels tdi_td */
    double alpha = (2.0*FILTER_LENGTH/Tobs);
    
    tukey(X, alpha, N);
    tukey(Y, alpha, N);
    tukey(Z, alpha, N);
    tukey(A, alpha, N);
    tukey(E, alpha, N);
    tukey(T, alpha, N);
    
    
    /* Fourier transform time-domain TDI channels */
    gsl_fft_real_wavetable * real = gsl_fft_real_wavetable_alloc (N);
    gsl_fft_real_workspace * work = gsl_fft_real_workspace_alloc (N);

    gsl_fft_real_transform (X, 1, N, real, work);
    gsl_fft_real_transform (Y, 1, N, real, work);
    gsl_fft_real_transform (Z, 1, N, real, work);
    gsl_fft_real_transform (A, 1, N, real, work);
    gsl_fft_real_transform (E, 1, N, real, work);
    gsl_fft_real_transform (T, 1, N, real, work);

    /* Normalize FD data */
    double rft_norm = sqrt(Tobs)/(double)N;
    
    /* Account for losses from windowing
    double tukey_norm = tukey_scale(alpha, N);
    rft_norm /= tukey_norm;
    */
    
    for(int n=0; n<N; n++)
    {
        X[n] *= rft_norm;
        Y[n] *= rft_norm;
        Z[n] *= rft_norm;
        A[n] *= rft_norm;
        E[n] *= rft_norm;
        T[n] *= rft_norm;
    }
        
    /* unpack GSL-formatted arrays to the way GLASS expects them */
    unpack_gsl_rft_output(tdi->X, X, N);
    unpack_gsl_rft_output(tdi->Y, Y, N);
    unpack_gsl_rft_output(tdi->Z, Z, N);
    unpack_gsl_rft_output(tdi->A, A, N);
    unpack_gsl_rft_output(tdi->E, E, N);
    unpack_gsl_rft_output(tdi->T, T, N);
    
    /* Free memory */
    gsl_fft_real_wavetable_free (real);
    gsl_fft_real_workspace_free (work);
    free_tdi(tdi_td);
    free(X);
    free(Y);
    free(Z);
    free(A);
    free(E);
    free(T);
    
}

void ReadASCII(struct Data *data, struct TDI *tdi)
{
    double f;
    double junk;
    
    FILE *fptr = fopen(data->fileName,"r");
    
    //count number of samples
    int Nsamples = 0;
    while(!feof(fptr))
    {
        int check = fscanf(fptr,"%lg %lg %lg %lg %lg",&f,&junk,&junk,&junk,&junk);
        if(!check)
        {
            fprintf(stderr,"Error reading %s\n",data->fileName);
            exit(1);
        }
        Nsamples++;
    }
    rewind(fptr);
    Nsamples--;
    
    //load full dataset into TDI structure
    alloc_tdi(tdi, Nsamples, 3);
    
    for(int n=0; n<Nsamples; n++)
    {
        int check = fscanf(fptr,"%lg %lg %lg %lg %lg",&f,&tdi->A[2*n],&tdi->A[2*n+1],&tdi->E[2*n],&tdi->E[2*n+1]);
        if(!check)
        {
            fprintf(stderr,"Error reading %s\n",data->fileName);
            exit(1);
        }
        
    }
    fclose(fptr);
}

void ReadData(struct Data *data, struct Orbit *orbit, struct Flags *flags)
{
    if(!flags->quiet) fprintf(stdout,"\n==== ReadData ====\n");
    
    struct TDI *tdi = data->tdi;
    
    
    /* load full dataset */
    struct TDI *tdi_full = malloc(sizeof(struct TDI));
    if(flags->hdf5Data)
        ReadHDF5(data,tdi_full,flags);
    else
        ReadASCII(data,tdi_full);
    
    
    /* select frequency segment */
    
    //get max and min samples
    data->fmax = data->fmin + data->N/data->T;
    data->qmin = (int)(data->fmin*data->T);
    data->qmax = data->qmin+data->N;
    
    //store frequency segment in TDI structure
    for(int n=0; n<2*data->N; n++)
    {
        int m = data->qmin*2+n;
        tdi->X[n] = tdi_full->X[m];
        tdi->Y[n] = tdi_full->Y[m];
        tdi->Z[n] = tdi_full->Z[m];
        tdi->A[n] = tdi_full->A[m];
        tdi->E[n] = tdi_full->E[m];
        tdi->T[n] = tdi_full->T[m];
    }
    
    //Get noise spectrum for data segment
    GetNoiseModel(data,orbit,flags);
    
    //Add Gaussian noise to injection
    if(flags->simNoise) AddNoise(data,tdi);
    
    //print various data products for plotting
    print_data(data, tdi, flags);
    
    //free memory
    free_tdi(tdi_full);
}

void GetNoiseModel(struct Data *data, struct Orbit *orbit, struct Flags *flags)
{
    double Spm, Sop;
    
    //if you are simulating/fitting the noise
    if(!flags->psd)
    {
        for(int n=0; n<data->N; n++)
        {
            double f = data->fmin + (double)(n)/data->T;
            data->noise->f[n] = f;
            data->noise->transfer[n] = noise_transfer_function(f/orbit->fstar);

            if(strcmp(data->format,"phase")==0)
            {
                data->noise->C[0][0][n] = AEnoise(orbit->L, orbit->fstar, f);
                data->noise->C[1][1][n] = AEnoise(orbit->L, orbit->fstar, f);
                data->noise->C[0][1][n] = 0.0;
                if(flags->confNoise)
                {
                    data->noise->C[0][0][n] += GBnoise(data->T,f);
                    data->noise->C[1][1][n] += GBnoise(data->T,f);
                }
            }
            else if(strcmp(data->format,"frequency")==0)
            {
                get_noise_levels("radler",f,&Spm,&Sop);
                data->noise->C[0][0][n] = AEnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);
                data->noise->C[1][1][n] = AEnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);
                data->noise->C[0][1][n] = 0.0;
                if(flags->confNoise)
                {
                    data->noise->C[0][0][n] += GBnoise_FF(data->T, orbit->fstar, f);
                    data->noise->C[1][1][n] += GBnoise_FF(data->T, orbit->fstar, f);
                }
            }
            else if(strcmp(data->format,"sangria")==0)
            {
                //TODO: Need a sqrt(2) to match Sangria data/noise
                get_noise_levels("sangria",f,&Spm,&Sop);
                data->noise->C[0][0][n] = AEnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);
                data->noise->C[1][1][n] = AEnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);
                data->noise->C[0][1][n] = 0.0;
                if(flags->confNoise)
                {
                    data->noise->C[0][0][n] += GBnoise_FF(data->T, orbit->fstar, f);
                    data->noise->C[1][1][n] += GBnoise_FF(data->T, orbit->fstar, f);
                }
            }
            else
            {
                fprintf(stderr,"Unsupported data format %s\n",data->format);
                exit(1);
            }

            //TODO: 3-channel model only has support for Sangria data conventions
            if(data->Nchannel==3)
            {
                //TODO: Need a sqrt(2) to match Sangria data/noise
                get_noise_levels("sangria",f,&Spm,&Sop);
                data->noise->C[0][0][n] = XYZnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);
                data->noise->C[1][1][n] = XYZnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);
                data->noise->C[2][2][n] = XYZnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);

                data->noise->C[0][1][n] = data->noise->C[1][0][n] = XYZcross_FF(orbit->L, orbit->fstar, f, Spm, Sop);
                data->noise->C[0][2][n] = data->noise->C[2][0][n] = XYZcross_FF(orbit->L, orbit->fstar, f, Spm, Sop);
                data->noise->C[1][2][n] = data->noise->C[2][1][n] = XYZcross_FF(orbit->L, orbit->fstar, f, Spm, Sop);

                if(flags->confNoise)
                {
                    double GBnoise=GBnoise_FF(data->T, orbit->fstar, f)/1.5; //GBnoise_FF() is hard-coded for AE channels
                    data->noise->C[0][0][n] += GBnoise;
                    data->noise->C[1][1][n] += GBnoise;
                    data->noise->C[2][2][n] += GBnoise;
                    data->noise->C[0][1][n] += -0.5*GBnoise;
                    data->noise->C[0][2][n] += -0.5*GBnoise;
                    data->noise->C[1][2][n] += -0.5*GBnoise;
                    data->noise->C[1][0][n] += -0.5*GBnoise;
                    data->noise->C[2][0][n] += -0.5*GBnoise;
                    data->noise->C[2][1][n] += -0.5*GBnoise;
                }

            }
        }
        
        invert_noise_covariance_matrix(data->noise);

    }
    //use PSD from file
    else
    {
        
        //parse input PSD file
        FILE *psdFile = fopen(flags->psdFile,"r");
        int lines=0;
        double f_temp, SnA_temp, SnE_temp;
        while(!feof(psdFile))
        {
            fscanf(psdFile,"%lg %lg %lg",&f_temp,&SnA_temp,&SnE_temp);
            lines++;
        }
        rewind(psdFile);
        lines--;
        
        double *f   = malloc(lines*sizeof(double));
        double *SnA = malloc(lines*sizeof(double));
        double *SnE = malloc(lines*sizeof(double));
        
        for(int l=0; l<lines; l++) fscanf(psdFile,"%lg %lg %lg",&f[l],&SnA[l],&SnE[l]);
        
        //interpolate input psd onto segment grid
        double *fint = malloc(data->N*sizeof(double));
        for(int n=0; n<data->N; n++) fint[n] = data->fmin + (double)(n)/data->T;

        CubicSplineGSL(lines, f, SnA, data->N, fint, data->noise->C[0][0]);
        CubicSplineGSL(lines, f, SnE, data->N, fint, data->noise->C[1][1]);
        
        free(f);
        free(SnA);
        free(SnE);
        fclose(psdFile);
    }
}

void AddNoise(struct Data *data, struct TDI *tdi)
{
    
    printf("   ...adding Gaussian noise realization\n");
    
    //set RNG for noise
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc(T);
    gsl_rng_env_setup();
    gsl_rng_set (r, data->nseed);
    
    double n_re[data->Nchannel];
    double n_im[data->Nchannel];
    double u_re[data->Nchannel];
    double u_im[data->Nchannel];
    
    //get LU decomposition of covariance matrix
    double **L = malloc(data->Nchannel*sizeof(double*));
    double **C = malloc(data->Nchannel*sizeof(double*));
    for(int i=0; i<data->Nchannel; i++)
    {
        L[i] = malloc(data->Nchannel*sizeof(double));
        C[i] = malloc(data->Nchannel*sizeof(double));
    }
    
    
    
    for(int n=0; n<data->N; n++)
    {
        for(int i=0; i<data->Nchannel; i++)
        {
            u_re[i] = gsl_ran_gaussian (r,1);
            u_im[i] = gsl_ran_gaussian (r,1);
            n_re[i] = n_im[i] = 0.0;
        }
 
        // make sure both diagonals of the covariance matrix are filled
        for(int i=0; i<data->Nchannel; i++)
            for(int j=i; j<data->Nchannel; j++)
                C[i][j] = C[j][i] = data->noise->C[i][j][n];

        cholesky_decomp(C, L, data->Nchannel);

        // n = Lu
        for(int i=0; i<data->Nchannel; i++)
        {
            for(int j=0; j<data->Nchannel; j++)
            {
                    n_re[i] += L[i][j]*u_re[j]/sqrt(2.);
                    n_im[i] += L[i][j]*u_im[j]/sqrt(2.);
            }
        }
        
        switch(data->Nchannel)
        {
            case 1:
                tdi->X[2*n]   += n_re[0];
                tdi->X[2*n+1] += n_im[0];
                break;
            case 2:
                tdi->A[2*n]   += n_re[0];
                tdi->A[2*n+1] += n_im[0];
                tdi->E[2*n]   += n_re[1];
                tdi->E[2*n+1] += n_im[1];
                break;
            case 3:
                tdi->X[2*n]   += n_re[0];
                tdi->X[2*n+1] += n_im[0];
                tdi->Y[2*n]   += n_re[1];
                tdi->Y[2*n+1] += n_im[1];
                tdi->Z[2*n]   += n_re[2];
                tdi->Z[2*n+1] += n_im[2];
                break;
        }
    }

    gsl_rng_free(r);
    for(int i=0; i<data->Nchannel; i++)
    {
        free(L[i]);
        free(C[i]);
    }
    free(L);
    free(C);
}

void SimulateData(struct Data *data, struct Orbit *orbit, struct Flags *flags)
{
    if(!flags->quiet) fprintf(stdout,"\n==== SimulateData ====\n");
    struct TDI *tdi = data->tdi;

    //get max and min samples
    data->fmax = data->fmin + data->N/data->T;
    data->qmin = (int)(data->fmin*data->T);
    data->qmax = data->qmin+data->N;

    //Get noise spectrum for data segment
    GetNoiseModel(data,orbit,flags);
    
    //Add Gaussian noise to injection
    if(flags->simNoise) AddNoise(data,tdi);
    
    //print various data products for plotting
    print_data(data, tdi, flags);

}

void print_data(struct Data *data, struct TDI *tdi, struct Flags *flags)
{
    char filename[128];
    FILE *fptr;
    
    sprintf(filename,"%s/power_data.dat",data->dataDir);
    fptr=fopen(filename,"w");
    
    for(int i=0; i<data->N; i++)
    {
        double f = (double)(i+data->qmin)/data->T;
        switch(data->Nchannel)
        {
            case 1:
                fprintf(fptr,"%.12g %lg\n", f, tdi->X[2*i]*tdi->X[2*i]+tdi->X[2*i+1]*tdi->X[2*i+1]);
                break;
            case 2:
                fprintf(fptr,"%.12g %lg %lg\n", f, tdi->A[2*i]*tdi->A[2*i]+tdi->A[2*i+1]*tdi->A[2*i+1], tdi->E[2*i]*tdi->E[2*i]+tdi->E[2*i+1]*tdi->E[2*i+1]);
                break;
            case 3:
                fprintf(fptr,"%.12g %lg %lg %lg\n", f, tdi->X[2*i]*tdi->X[2*i]+tdi->X[2*i+1]*tdi->X[2*i+1], tdi->Y[2*i]*tdi->Y[2*i]+tdi->Y[2*i+1]*tdi->Y[2*i+1], tdi->Z[2*i]*tdi->Z[2*i]+tdi->Z[2*i+1]*tdi->Z[2*i+1]);
                break;
        }
    }
    fclose(fptr);
    
    sprintf(filename,"%s/data.dat",data->dataDir);
    fptr=fopen(filename,"w");
    
    for(int i=0; i<data->N; i++)
    {
        double f = (double)(i+data->qmin)/data->T;
        switch(data->Nchannel)
        {
            case 1:
                fprintf(fptr,"%.12g %lg %lg\n", f, tdi->X[2*i],tdi->X[2*i+1]);
                break;
            case 2:
                fprintf(fptr,"%.12g %lg %lg %lg %lg\n", f, tdi->A[2*i],tdi->A[2*i+1], tdi->E[2*i],tdi->E[2*i+1]);
                break;
            case 3:
                fprintf(fptr,"%.12g %lg %lg %lg %lg %lg %lg\n", f, tdi->X[2*i],tdi->X[2*i+1], tdi->Y[2*i],tdi->Y[2*i+1], tdi->Z[2*i],tdi->Z[2*i+1]);
                break;
        }
    }
    fclose(fptr);
    
    sprintf(filename,"%s/power_noise.dat",data->dataDir);
    fptr=fopen(filename,"w");
    
    for(int i=0; i<data->N; i++)
    {
        double f = (double)(i+data->qmin)/data->T;
        switch(data->Nchannel)
        {
            case 1:
                fprintf(fptr,"%.12g %lg\n", f, data->noise->C[0][0][i]);
                break;
            case 2:
                fprintf(fptr,"%.12g %lg %lg\n", f, data->noise->C[0][0][i], data->noise->C[1][1][i]);
                break;
            case 3:
                fprintf(fptr,"%.12g %lg %lg %lg\n", f, data->noise->C[0][0][i], data->noise->C[1][1][i], data->noise->C[2][2][i]);
                break;
        }
    }
    fclose(fptr);
}

void print_glass_usage()
{
    fprintf(stdout,"\n");
    fprintf(stdout,"=============== GLASS Usage: ============== \n");
    fprintf(stdout,"REQUIRED:\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"OPTIONAL:\n");
    fprintf(stdout,"  -h | --help        : print help message and exit         \n");
    fprintf(stdout,"  -v | --verbose     : enable verbose output               \n");
    fprintf(stdout,"  -q | --quiet       : restrict output                     \n");
    fprintf(stdout,"  -d | --debug       : leaner settings for quick running   \n");
    fprintf(stdout,"\n");
    
    //LISA
    fprintf(stdout,"       =========== LISA =========== \n");
    fprintf(stdout,"       --orbit       : orbit ephemerides file (2.5 GM MLDC)\n");
    fprintf(stdout,"       --channels    : # of channels [1->X,2->AE,3->XYZ](3)\n");
    fprintf(stdout,"       --phase       : phase data (fractional frequency)   \n");
    fprintf(stdout,"       --sangria     : use LDC Sangria TDI conventions     \n");
    fprintf(stdout,"\n");
    
    //Data
    fprintf(stdout,"       =========== Data =========== \n");
    fprintf(stdout,"       --data        : strain data file (ASCII)            \n");
    fprintf(stdout,"       --h5-data     : strain data file (HDF5)             \n");
    fprintf(stdout,"       --h5-no-mbh   : remove mbhs from HDF5 data          \n");
    fprintf(stdout,"       --h5-no-ucb   : remove ucbs from HDF5 data          \n");
    fprintf(stdout,"       --h5-no-vgb   : remove vgbs from HDF5 data          \n");
    fprintf(stdout,"       --h5-no-noise : remove noise from HDF5 data (TODO)  \n");
    fprintf(stdout,"       --psd         : psd data file (ASCII)               \n");
    fprintf(stdout,"       --samples     : number of frequency bins (512)      \n");
    fprintf(stdout,"       --padding     : number of bins padded on segment (0)\n");
    fprintf(stdout,"       --start-time  : initial time of epoch  (0)          \n");
    fprintf(stdout,"       --fmin        : minimum frequency                   \n");
    fprintf(stdout,"       --fmax        : maximum frequency (overrides --samples)\n");
    fprintf(stdout,"       --duration    : duration of epoch (31457280)        \n");
    fprintf(stdout,"       --sim-noise   : data w/out noise realization        \n");
    fprintf(stdout,"       --conf-noise  : include model for confusion noise   \n");
    fprintf(stdout,"       --noiseseed   : seed for noise RNG                  \n");
    fprintf(stdout,"\n");
    
    //Chain
    fprintf(stdout,"       ========== Chains ========== \n");
    fprintf(stdout,"       --steps       : number of mcmc steps (10000)        \n");
    fprintf(stdout,"       --chainseed   : seed for MCMC RNG                   \n");
    fprintf(stdout,"       --chains      : number of parallel chains (20)      \n");
    fprintf(stdout,"       --no-burnin   : skip burn in steps                  \n");
    fprintf(stdout,"       --resume      : restart from checkpoint             \n");
    fprintf(stdout,"       --threads     : number of parallel threads (max)    \n");
    fprintf(stdout,"       --prior       : sample from prior                   \n");
    fprintf(stdout,"       --no-rj       : turn off RJMCMC                     \n");
    fprintf(stdout,"\n");
    
    //Misc.
    fprintf(stdout,"       =========== Misc =========== \n");
    fprintf(stdout,"       --rundir      : top level run directory ['./']\n");
    fprintf(stdout,"       --match-in1   : input paramaters for overlap [filename] \n");
    fprintf(stdout,"       --match-in2   : output match values [filename] \n");
    fprintf(stdout,"\n");

    /*
    fprintf(stdout,"USEFUL TOBS:\n");
    fprintf(stdout,"   1 wk: %.0f\n",62914560./2./52.);
    fprintf(stdout,"   1 mo: %.0f \n",62914560./2./12.);
    fprintf(stdout,"   2 mo: %.0f \n",62914560./2./6.);
    fprintf(stdout,"   1 yr: %.0f (default)\n",62914560./2.);
    fprintf(stdout,"   2 yr: %.0f \n",62914560.);
    fprintf(stdout,"   5 yr: %.0f \n",5.*62914560./2.);
    fprintf(stdout,"  10 yr: %.0f \n",10.*62914560./2.);
     */
    fprintf(stdout,"\n");
}

void parse_data_args(int argc, char **argv, struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Chain *chain)
{
    //copy argv since getopt permutes order
    char **argv_copy=malloc((argc+1) * sizeof *argv_copy);
    copy_argv(argc,argv,argv_copy);
    opterr=0; //suppress warnings about unknown arguments
    
    //Set defaults
    flags->help        = 0;
    flags->calibration = 0;
    flags->verbose     = 0;
    flags->quiet       = 0;
    flags->simNoise    = 0;
    flags->confNoise   = 0;
    flags->burnin      = 1;
    flags->debug       = 0;
    flags->strainData  = 0;
    flags->hdf5Data    = 0;
    flags->psd         = 0;
    flags->orbit       = 0;
    flags->prior       = 0;
    flags->resume      = 0;
    flags->NMCMC       = 10000;
    flags->NBURN       = 10000;
    flags->threads     = omp_get_max_threads();
    sprintf(flags->runDir,"./");
    chain->NC          = 12;//number of chains
    int set_fmax_flag  = 0; //flag watching for if fmax is set by CLI
    
    /* Simulated data building blocks */
    flags->no_mbh = 0;
    flags->no_ucb = 0;
    flags->no_vgb = 0;
    flags->no_noise = 0;
    
    /*
     default data format is 'phase'
     optional support for 'frequency' a la LDCs
     */
    sprintf(data->format,"sangria");
        
    data->T        = 31457280; /* one "mldc years" at 15s sampling */
    data->t0       = 0.0; /* start time of data segment in seconds */
    data->sqT      = sqrt(data->T);
    data->N        = 512;
    data->Nchannel = 3; //1=X, 2=AE, 3=XYZ
    data->qpad     = 0;
    data->fmin     = 1e-4; //Hz
    
    data->cseed = 150914;
    data->nseed = 151226;
    data->iseed = 151012;
        
    //Specifying the expected options
    static struct option long_options[] =
    {
        /* These options set a flag. */
        {"samples",    required_argument, 0, 0},
        {"padding",    required_argument, 0, 0},
        {"duration",   required_argument, 0, 0},
        {"start-time", required_argument, 0, 0},
        {"orbit",      required_argument, 0, 0},
        {"chains",     required_argument, 0, 0},
        {"chainseed",  required_argument, 0, 0},
        {"noiseseed",  required_argument, 0, 0},
        {"data",       required_argument, 0, 0},
        {"h5-data",    required_argument, 0, 0},
        {"psd",        required_argument, 0, 0},
        {"fmin",       required_argument, 0, 0},
        {"fmax",       required_argument, 0, 0},
        {"channels",   required_argument, 0, 0},
        {"steps",      required_argument, 0, 0},
        {"threads",    required_argument, 0, 0},
        {"rundir",     required_argument, 0, 0},
        
        /* These options don’t set a flag.
         We distinguish them by their indices. */
        {"help",        no_argument, 0,'h'},
        {"verbose",     no_argument, 0,'v'},
        {"quiet",       no_argument, 0,'q'},
        {"debug",       no_argument, 0,'d'},
        {"resume",      no_argument, 0, 0 },
        {"sim-noise",   no_argument, 0, 0 },
        {"conf-noise",  no_argument, 0, 0 },
        {"phase",       no_argument, 0, 0 },
        {"sangria",     no_argument, 0, 0 },
        {"prior",       no_argument, 0, 0 },
        {"no-burnin",   no_argument, 0, 0 },
        {"no-rj",       no_argument, 0, 0 },
        {"calibration", no_argument, 0, 0 },
        {"h5-no-mbh",   no_argument, 0, 0 },
        {"h5-no-ucb",   no_argument, 0, 0 },
        {"h5-no-vgb",   no_argument, 0, 0 },
        {"h5-no-noise", no_argument, 0, 0 },
        {0, 0, 0, 0}
    };
    
    int opt=0;
    int long_index=0;
    
    //Loop through argv string and pluck out arguments
    while ((opt = getopt_long_only(argc, argv_copy,"apl:b:", long_options, &long_index )) != -1)
    {
        switch (opt)
        {
                
            case 0:
                if(strcmp("samples",     long_options[long_index].name) == 0) data->N           = atoi(optarg);
                if(strcmp("padding",     long_options[long_index].name) == 0) data->qpad        = atoi(optarg);
                if(strcmp("start-time",  long_options[long_index].name) == 0) data->t0       = (double)atof(optarg);
                if(strcmp("chains",      long_options[long_index].name) == 0) chain->NC         = atoi(optarg);
                if(strcmp("chainseed",   long_options[long_index].name) == 0) data->cseed       = (long)atoi(optarg);
                if(strcmp("noiseseed",   long_options[long_index].name) == 0) data->nseed       = (long)atoi(optarg);
                if(strcmp("injseed",     long_options[long_index].name) == 0) data->iseed       = (long)atoi(optarg);
                if(strcmp("sim-noise",   long_options[long_index].name) == 0) flags->simNoise   = 1;
                if(strcmp("conf-noise",  long_options[long_index].name) == 0) flags->confNoise  = 1;
                if(strcmp("prior",       long_options[long_index].name) == 0) flags->prior      = 1;
                if(strcmp("no-burnin",   long_options[long_index].name) == 0) flags->burnin     = 0;
                if(strcmp("no-rj",       long_options[long_index].name) == 0) flags->rj         = 0;
                if(strcmp("calibration", long_options[long_index].name) == 0) flags->calibration= 1;
                if(strcmp("resume",      long_options[long_index].name) == 0) flags->resume     = 1;
                if(strcmp("h5-no-mbh",   long_options[long_index].name) == 0) flags->no_mbh     = 1;
                if(strcmp("h5-no-ucb",   long_options[long_index].name) == 0) flags->no_ucb     = 1;
                if(strcmp("h5-no-vgb",   long_options[long_index].name) == 0) flags->no_vgb     = 1;
                if(strcmp("h5-no-noise", long_options[long_index].name) == 0) flags->no_noise   = 1;
                if(strcmp("threads",     long_options[long_index].name) == 0) flags->threads    = atoi(optarg);
                if(strcmp("rundir",      long_options[long_index].name) == 0) strcpy(flags->runDir,optarg);
                if(strcmp("phase",       long_options[long_index].name) == 0) sprintf(data->format,"phase");
                if(strcmp("sangria",     long_options[long_index].name) == 0) sprintf(data->format,"sangria");
                if(strcmp("fmin",        long_options[long_index].name) == 0) sscanf(optarg, "%lg", &data->fmin);
                if(strcmp("fmax",        long_options[long_index].name) == 0)
                {
                    set_fmax_flag = 1;
                    sscanf(optarg, "%lg", &data->fmax);
                }
                if(strcmp("duration",    long_options[long_index].name) == 0)
                {   
                    data->T   = (double)atof(optarg);
                    data->sqT = sqrt(data->T);
                }
                if(strcmp("steps",       long_options[long_index].name) == 0)
                {
                    flags->NMCMC = atoi(optarg);
                    flags->NBURN = flags->NMCMC;
                }
                if(strcmp("data", long_options[long_index].name) == 0)
                {
                    checkfile(optarg);
                    flags->strainData = 1;
                    sprintf(data->fileName,"%s",optarg);
                }
                if(strcmp("h5-data", long_options[long_index].name) == 0)
                {
                    checkfile(optarg);
                    flags->hdf5Data = 1;
                    flags->strainData = 1;
                    sprintf(data->fileName,"%s",optarg);
                }
                if(strcmp("psd", long_options[long_index].name) == 0)
                {
                    checkfile(optarg);
                    flags->psd = 1;
                    sprintf(flags->psdFile,"%s",optarg);
                }
                if(strcmp("orbit", long_options[long_index].name) == 0)
                {
                    checkfile(optarg);
                    flags->orbit = 1;
                    sprintf(orbit->OrbitFileName,"%s",optarg);
                }
                if(strcmp("channels",long_options[long_index].name) == 0)
                {
                    data->Nchannel = (int)atoi(optarg);
                    if(data->Nchannel<1 || data->Nchannel>3)
                    {
                        fprintf(stderr,"Requested umber of channels (%i) not supported\n",data->Nchannel);
                        fprintf(stderr,"Use --channels 1 for X (Michelson) data\n");
                        fprintf(stderr,"    --channels 2 for AE data\n");
                        fprintf(stderr,"    --channels 3 for XYZ data\n");
                        exit(1);
                    }
                }
                break;
            case 'd' : flags->debug = 1;
                break;
            case 'h' : flags->help = 1;
                break;
            case 'v' : flags->verbose = 1;
                break;
            case 'q' : flags->quiet = 1;
                break;
            default:
                break;
        }
    }
    if(flags->cheat || !flags->burnin) flags->NBURN = 0;
    
    if(flags->verbose && flags->quiet)
    {
        fprintf(stderr,"--verbose and --quiet flags are in conflict\n");
        exit(1);
    }
    
    //Chains should be a multiple of threads for best usage of cores
    if(chain->NC % flags->threads !=0){
        chain->NC += flags->threads - (chain->NC % flags->threads);
    }
    
    //override size of data if fmax was requested
    if(set_fmax_flag) data->N = (int)floor((data->fmax - data->fmin)*data->T);
    
    //pad data
    data->N += 2*data->qpad;
    data->fmin -= data->qpad/data->T;
        
    //map fmin to nearest bin
    data->fmin = floor(data->fmin*data->T)/data->T;
    data->fmax = data->fmin + (double)data->N/data->T;

    //calculate helper quantities for likelihood normalizations
    data->logfmin   = log(data->fmin);
    data->sum_log_f = 0.0;
    for(int n=0; n<data->N; n++)
    {
        data->sum_log_f += log(data->fmin + (double)n/data->T);
    }
    
    //Print version control
//    sprintf(filename,"glass.log");
//    FILE *runlog = fopen(filename,"w");
//    print_version(runlog);
    
    //Report on set parameters
//    if(!flags->quiet) print_run_settings(argc, argv, data, orbit, flags, stdout);
//    print_run_settings(argc, argv, data, orbit, flags, runlog);
    
//    fclose(runlog);
    
    //reset opt counter
    optind = 0;

    //free placeholder for argvs
    for(int i=0; i<=argc; i++)free(argv_copy[i]);
    free(argv_copy);

}

void copy_argv(int argc, char **argv, char **new_argv)
{
    for(int i = 0; i < argc; ++i)
    {
        size_t length = strlen(argv[i])+1;
        new_argv[i] = malloc(length);
        memcpy(new_argv[i], argv[i], length);
    }
    new_argv[argc] = NULL;
}

int checkfile(char filename[])
{
    FILE *fptr = fopen(filename, "r");
    if(fptr)
    {
        fclose(fptr);
        return 1;
    }
    else
    {
        fprintf(stderr,"File %s does not exist\n",filename);
        fprintf(stderr,"\n");
        exit(1);
    }
}
