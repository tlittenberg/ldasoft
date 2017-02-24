//
//  GalacticBinaryIO.c
//
//
//  Created by Littenberg, Tyson B. (MSFC-ZP12) on 1/15/17.
//
//

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

#include "LISA.h"
#include "GalacticBinary.h"
#include "GalacticBinaryIO.h"
#include "GalacticBinaryModel.h"

void print_usage()
{
  fprintf(stdout,"\n");
  fprintf(stdout,"Usage: \n");
  fprintf(stdout,"REQUIRED:\n");
  fprintf(stdout,"       --orbit       : orbit ephemerides file             \n");
  fprintf(stdout,"OPTIONAL:\n");
  fprintf(stdout,"  -h | --help        : print help message and exit        \n");
  fprintf(stdout,"  -v | --verbose     : enable verbose output              \n");
  fprintf(stdout,"       --samples     : number of frequency bins (2048)    \n");
  fprintf(stdout,"       --start-time  : initial time of segment  (0)       \n");
  fprintf(stdout,"       --duration    : duration of time segment (62914560)\n");
  fprintf(stdout,"       --noiseseed   : seed for noise RNG                 \n");
  fprintf(stdout,"       --chainseed   : seed for MCMC RNG                  \n");
  fprintf(stdout,"       --injseed     : seed for injection parameters      \n");
  fprintf(stdout,"       --inj         : inject signal                      \n");
  fprintf(stdout,"       --zero-noise  : data w/out noise realization       \n");
  fprintf(stdout,"       --links       : number of links [4->X,6->AE] (6)   \n");
  fprintf(stdout,"--\n");
  fprintf(stdout,"EXAMPLE:\n");
  fprintf(stdout,"./gb_mcmc --orbit ../config/OrbitConfig1.txt --verbose --inj ../data/sources/RXJ0806.dat\n");
  fprintf(stdout,"\n");
  exit(EXIT_FAILURE);
}

void parse(int argc, char **argv, struct Data **data, struct Orbit *orbit, struct Flags *flags, int Nmax)
{
  //Set defaults
  flags->verbose   = 0;
  flags->injection = 0;
  flags->zeroNoise = 0;
  
  for(int i=0; i<Nmax; i++)
  {
    data[i]->t0       = 0.0;
    data[i]->T        = 62914560.0; /* two "mldc years" at 15s sampling */
    data[i]->N        = 256;
    data[i]->Nchannel = 2; //1=X, 2=AE
    
    data[i]->cseed = 150914+i;
    data[i]->nseed = 151226+i;
    data[i]->iseed = 151012+i;
  }
  
  flags->injFile = malloc(10*sizeof(char *));
  for(int n=0; n<10; n++) flags->injFile[n] = malloc(1024*sizeof(char));
  
  if(argc==1) print_usage();
  
  //Specifying the expected options
  static struct option long_options[] =
  {
    /* These options set a flag. */
    {"samples",   required_argument, 0, 0},
    {"duration",  required_argument, 0, 0},
    {"start-time",required_argument, 0, 0},
    {"orbit",     required_argument, 0, 0},
    {"chainseed", required_argument, 0, 0},
    {"noiseseed", required_argument, 0, 0},
    {"injseed",   required_argument, 0, 0},
    {"inj",       required_argument, 0, 0},
    {"links",     required_argument, 0, 0},
    
    /* These options don’t set a flag.
     We distinguish them by their indices. */
    {"help",    no_argument,       0,  'h' },
    {"verbose", no_argument,       0,  'v' },
    {"zero-noise", no_argument,    0,   0  },
    {0,         0,                 0,   0  }
  };
  
  int opt=0;
  int long_index=0;
  
  //Print command line
  FILE *out = fopen("gb_mcmc.sh","w");
  fprintf(out,"#!/bin/sh\n\n");
  for(opt=0; opt<argc; opt++) fprintf(out,"%s ",argv[opt]);
  fprintf(out,"\n\n");
  fclose(out);
  
  //Loop through argv string and pluck out arguments
  while ((opt = getopt_long_only(argc, argv,"apl:b:", long_options, &long_index )) != -1)
  {
    switch (opt)
    {
        
      case 0:
        if(strcmp("samples",   long_options[long_index].name) == 0) data[0]->N     = (long)atof(optarg);
        if(strcmp("duration",  long_options[long_index].name) == 0) data[0]->T     = (double)atof(optarg);
        if(strcmp("start-time",long_options[long_index].name) == 0) data[0]->t0    = (double)atof(optarg);
        if(strcmp("chainseed", long_options[long_index].name) == 0) data[0]->cseed = (long)atoi(optarg);
        if(strcmp("noiseseed", long_options[long_index].name) == 0) data[0]->nseed = (long)atoi(optarg);
        if(strcmp("injseed",   long_options[long_index].name) == 0) data[0]->iseed = (long)atoi(optarg);
        if(strcmp("zero-noise",long_options[long_index].name) == 0) flags->zeroNoise = 1;
        if(strcmp("orbit",     long_options[long_index].name) == 0) sprintf(orbit->OrbitFileName,"%s",optarg);
        if(strcmp("inj",       long_options[long_index].name) == 0)
        {
          sprintf(flags->injFile[flags->injection],"%s",optarg);
          flags->injection++;
          if(flags->injection>Nmax)
          {
            fprintf(stderr,"Requested number of injections is too large (%i/%i)\n",flags->injection,Nmax);
            fprintf(stderr,"Remove at least %i --inj arguments\n",flags->injection-Nmax);
            fprintf(stderr,"Now exiting to system\n");
            exit(1);
          }
        }
        if(strcmp("links",      long_options[long_index].name) == 0)
        {
          int Nlinks = (int)atoi(optarg);
          switch(Nlinks)
          {
            case 4:
              data[0]->Nchannel=1;
              break;
            case 6:
              data[0]->Nchannel=2;
              break;
            default:
              fprintf(stderr,"Requested umber of links (%i) not supported\n",Nlinks);
              fprintf(stderr,"Use --links 4 for X (Michelson) data\n");
              fprintf(stderr,"    --links 6 for AE data\n");
              exit(1);
          }
        }
        break;
      case 'h' :
        print_usage();
        exit(EXIT_FAILURE);
        break;
      case 'v' : flags->verbose = 1;
        break;
      default: print_usage();
        exit(EXIT_FAILURE);
    }
  }
  
  // copy command line args to other data structures
  for(int i=1; i<Nmax; i++)
  {
    data[i]->t0       = data[0]->t0;
    data[i]->T        = data[0]->T;
    data[i]->N        = data[0]->N;
    data[i]->Nchannel = data[0]->Nchannel;
    
    data[i]->cseed = data[0]->cseed+i;
    data[i]->nseed = data[0]->nseed+i;
    data[i]->iseed = data[0]->iseed+i;
  }

  
  // check for required arguments
  int abort=0;
  
//  if(data->duration[0]=='\0')
//  {
//    printf("Missing required argument: --duration\n");
//    abort++;
//  }
//  else data->T = (double)atof(data->duration);
  
  if(abort>0)exit(EXIT_FAILURE);
  
  
  //Report on set parameters
  fprintf(stdout,"\n");
  fprintf(stdout,"=============== RUN SETTINGS ===============\n");
  fprintf(stdout,"\n");
  fprintf(stdout,"  GalacticBinary version: %s\n", VERSION);
  fprintf(stdout,"  Command Line: ");
  for(opt=0; opt<argc; opt++) fprintf(stdout,"%s ",argv[opt]);
  fprintf(stdout,"\n");
  fprintf(stdout,"\n");
  fprintf(stdout,"  Orbit file .......... %s   \n",orbit->OrbitFileName);
  fprintf(stdout,"  Data channels ........");
  switch(data[0]->Nchannel)
  {
    case 1:
      fprintf(stdout,"X\n");
      break;
    case 2:
      fprintf(stdout,"AE\n");
      break;
  }
  fprintf(stdout,"  Data sample size .... %i   \n",data[0]->N);
  fprintf(stdout,"  Data start time ..... %.0f \n",data[0]->t0);
  fprintf(stdout,"  Data duration ....... %.0f \n",data[0]->T);
  fprintf(stdout,"  MCMC chain seed ..... %li  \n",data[0]->cseed);
  fprintf(stdout,"\n");
  fprintf(stdout,"================= RUN FLAGS ================\n");
  if(flags->verbose)  fprintf(stdout,"  Verbose flag ........ ENABLED \n");
  else                fprintf(stdout,"  Verbose flag ........ DISABLED\n");
  if(flags->injection>0)
  {
    fprintf(stdout,"  Injected sources..... %i\n",flags->injection);
    fprintf(stdout,"     seed ............. %li\n",data[0]->iseed);
    for(int i=0; i<flags->injection; i++)
    {
      fprintf(stdout,"     source ........... %s\n",flags->injFile[i]);
    }
  }
  else                fprintf(stdout,"  Injection is ........ DISABLED\n");
  if(flags->zeroNoise)fprintf(stdout,"  Noise realization is. DISABLED\n");
  else
  {
    fprintf(stdout,"  Noise realization is. ENABLED\n");
    fprintf(stdout,"  Noise seed .......... %li  \n",data[0]->nseed);
  }
  fprintf(stdout,"\n");
  fprintf(stdout,"\n");
  
}

void print_chain_files(struct Data *data, struct Model **model, struct Chain *chain, struct Flags *flags, int step)
{
  int n,ic;
  
  //Always print logL & temperature chains
  fprintf(chain->likelihoodFile,  "%i ",step);
  fprintf(chain->temperatureFile, "%i ",step);
  for(ic=0; ic<chain->NC; ic++)
  {
    n = chain->index[ic];
    fprintf(chain->likelihoodFile,  "%lg ",model[n]->logL+model[n]->logLnorm);
    fprintf(chain->temperatureFile, "%lg ",1./chain->temperature[ic]);
  }
  fprintf(chain->likelihoodFile, "\n");
  fprintf(chain->temperatureFile,"\n");
  
  //Always print cold chain
  n = chain->index[0];
  print_chain_state(data, chain, model[n], chain->parameterFile[0], step);
  print_noise_state(data, model[n], chain->noiseFile[0], step);

  //Print hot chains if verbose flag
  if(flags->verbose)
  {
    for(ic=1; ic<chain->NC; ic++)
    {
      n = chain->index[ic];
      print_chain_state(data, chain, model[n], chain->parameterFile[ic], step);
      print_noise_state(data, model[n], chain->noiseFile[ic], step);
    }//loop over chains
  }//verbose flag
}

void print_chain_state(struct Data *data, struct Chain *chain, struct Model *model, FILE *fptr, int step)
{
  for(int i=0; i<model->Nlive; i++)
  {
    fprintf(fptr, "%i ",step);
    fprintf(fptr, "%lg ",model->logL+model->logLnorm);
    print_source_params(data,model->source[i],fptr);
    fprintf(fptr, "\n");
  }
}

void print_noise_state(struct Data *data, struct Model *model, FILE *fptr, int step)
{
  fprintf(fptr, "%i ",step);
  fprintf(fptr, "%lg ",model->logL+model->logLnorm);

  switch(data->Nchannel)
  {
    case 1:
      fprintf(fptr, "%lg ", model->noise->etaX);
      break;
    case 2:
      fprintf(fptr, "%lg ", model->noise->etaA);
      fprintf(fptr, "%lg ", model->noise->etaE);
      break;
  }
  fprintf(fptr, "\n");
}

void print_source_params(struct Data *data, struct Source *source, FILE *fptr)
{
  //map to parameter names (just to make code readable)
  map_array_to_params(source, source->params, data->T);
  
  fprintf(fptr,"%.12g ",source->f0);
  fprintf(fptr,"%.12g ",source->dfdt);
  fprintf(fptr,"%lg ",source->amp);
  fprintf(fptr,"%lg ",source->phi);
  fprintf(fptr,"%lg ",source->costheta);
  fprintf(fptr,"%lg ",source->cosi);
  fprintf(fptr,"%lg ",source->psi);
  fprintf(fptr,"%lg ",source->phi0);
}

void save_waveforms(struct Data *data, struct Model *model, int mcmc)
{
  int n_re,n_im;
  double A_re,A_im,E_re,E_im,X_re,X_im,R_re,R_im;
  
  switch(data->Nchannel)
  {
    case 1:
      for(int n=0; n<data->N; n++)
      {
        n_re = 2*n;
        n_im = n_re++;
        
        X_re = model->tdi->X[n_re];
        X_im = model->tdi->X[n_im];
        
        data->h_rec[n_re][0][mcmc] = X_re;
        data->h_rec[n_im][0][mcmc] = X_im;
        
        R_re = data->tdi->X[n_re] - X_re;
        R_im = data->tdi->X[n_im] - X_im;
        
        data->h_res[n][0][mcmc] = R_re*R_re + R_im*R_im;
        data->h_pow[n][0][mcmc] = X_re*X_re + X_im*X_im;
        
        data->S_pow[n][0][mcmc] = model->noise->SnX[n];
      }
      break;
    case 2:
      for(int n=0; n<data->N; n++)
      {
        n_re = 2*n;
        n_im = n_re++;

        A_re = model->tdi->A[n_re];
        A_im = model->tdi->A[n_im];
        E_re = model->tdi->E[n_re];
        E_im = model->tdi->E[n_im];
        
        data->h_rec[n_re][0][mcmc] = A_re;
        data->h_rec[n_im][0][mcmc] = A_im;
        data->h_rec[n_re][1][mcmc] = E_re;
        data->h_rec[n_im][1][mcmc] = E_im;
        
        R_re = data->tdi->A[n_re] - A_re;
        R_im = data->tdi->A[n_im] - A_im;
        data->h_res[n][0][mcmc] = R_re*R_re + R_im*R_im;
        
        R_re = data->tdi->E[n_re] - E_re;
        R_im = data->tdi->E[n_im] - E_im;
        data->h_res[n][1][mcmc] = R_re*R_re + R_im*R_im;
        
        data->h_pow[n][0][mcmc] = A_re*A_re + A_im*A_im;
        data->h_pow[n][1][mcmc] = E_re*E_re + E_im*E_im;
        
        data->S_pow[n][0][mcmc] = model->noise->SnA[n];
        data->S_pow[n][1][mcmc] = model->noise->SnE[n];
      }
      break;
  }
}

void print_reconstructed_waveforms(struct Data *data)
{
  //sort h reconstructions
  for(int n=0; n<data->N*2; n++)
  {
    for(int m=0; m<data->Nchannel; m++)
    {
      gsl_sort(data->h_rec[n][m],1,data->Nwave);
    }
  }
  for(int n=0; n<data->N; n++)
  {
    for(int m=0; m<data->Nchannel; m++)
    {
      gsl_sort(data->h_res[n][m],1,data->Nwave);
      gsl_sort(data->h_pow[n][m],1,data->Nwave);
      gsl_sort(data->S_pow[n][m],1,data->Nwave);
    }
  }
  
  FILE *fptr_rec=fopen("power_reconstruction.dat","w");
  FILE *fptr_res=fopen("power_residual.dat","w");
  FILE *fptr_Snf=fopen("power_noise.dat","w");
  //double X_med,X_lo_50,X_hi_50,X_lo_90,X_hi_90;
  double A_med,A_lo_50,A_hi_50,A_lo_90,A_hi_90;
  double E_med,E_lo_50,E_hi_50,E_lo_90,E_hi_90;
  for(int i=0; i<data->N; i++)
  {
    double f = (double)(i+data->qmin)/data->T;
    
    A_med   = gsl_stats_median_from_sorted_data   (data->h_res[i][0], 1, data->Nwave);
    A_lo_50 = gsl_stats_quantile_from_sorted_data (data->h_res[i][0], 1, data->Nwave, 0.25);
    A_hi_50 = gsl_stats_quantile_from_sorted_data (data->h_res[i][0], 1, data->Nwave, 0.75);
    A_lo_90 = gsl_stats_quantile_from_sorted_data (data->h_res[i][0], 1, data->Nwave, 0.05);
    A_hi_90 = gsl_stats_quantile_from_sorted_data (data->h_res[i][0], 1, data->Nwave, 0.95);
    
    E_med   = gsl_stats_median_from_sorted_data   (data->h_res[i][1], 1, data->Nwave);
    E_lo_50 = gsl_stats_quantile_from_sorted_data (data->h_res[i][1], 1, data->Nwave, 0.25);
    E_hi_50 = gsl_stats_quantile_from_sorted_data (data->h_res[i][1], 1, data->Nwave, 0.75);
    E_lo_90 = gsl_stats_quantile_from_sorted_data (data->h_res[i][1], 1, data->Nwave, 0.05);
    E_hi_90 = gsl_stats_quantile_from_sorted_data (data->h_res[i][1], 1, data->Nwave, 0.95);
    
    fprintf(fptr_res,"%lg ",f);
    fprintf(fptr_res,"%lg ",A_med);
    fprintf(fptr_res,"%lg ",A_lo_50);
    fprintf(fptr_res,"%lg ",A_hi_50);
    fprintf(fptr_res,"%lg ",A_lo_90);
    fprintf(fptr_res,"%lg ",A_hi_90);
    fprintf(fptr_res,"%lg ",E_med);
    fprintf(fptr_res,"%lg ",E_lo_50);
    fprintf(fptr_res,"%lg ",E_hi_50);
    fprintf(fptr_res,"%lg ",E_lo_90);
    fprintf(fptr_res,"%lg ",E_hi_90);
    fprintf(fptr_res,"\n");
    
    A_med   = gsl_stats_median_from_sorted_data   (data->h_pow[i][0], 1, data->Nwave);
    A_lo_50 = gsl_stats_quantile_from_sorted_data (data->h_pow[i][0], 1, data->Nwave, 0.25);
    A_hi_50 = gsl_stats_quantile_from_sorted_data (data->h_pow[i][0], 1, data->Nwave, 0.75);
    A_lo_90 = gsl_stats_quantile_from_sorted_data (data->h_pow[i][0], 1, data->Nwave, 0.05);
    A_hi_90 = gsl_stats_quantile_from_sorted_data (data->h_pow[i][0], 1, data->Nwave, 0.95);
    
    E_med   = gsl_stats_median_from_sorted_data   (data->h_pow[i][1], 1, data->Nwave);
    E_lo_50 = gsl_stats_quantile_from_sorted_data (data->h_pow[i][1], 1, data->Nwave, 0.25);
    E_hi_50 = gsl_stats_quantile_from_sorted_data (data->h_pow[i][1], 1, data->Nwave, 0.75);
    E_lo_90 = gsl_stats_quantile_from_sorted_data (data->h_pow[i][1], 1, data->Nwave, 0.05);
    E_hi_90 = gsl_stats_quantile_from_sorted_data (data->h_pow[i][1], 1, data->Nwave, 0.95);
    
    
    fprintf(fptr_rec,"%lg ",f);
    fprintf(fptr_rec,"%lg ",A_med);
    fprintf(fptr_rec,"%lg ",A_lo_50);
    fprintf(fptr_rec,"%lg ",A_hi_50);
    fprintf(fptr_rec,"%lg ",A_lo_90);
    fprintf(fptr_rec,"%lg ",A_hi_90);
    fprintf(fptr_rec,"%lg ",E_med);
    fprintf(fptr_rec,"%lg ",E_lo_50);
    fprintf(fptr_rec,"%lg ",E_hi_50);
    fprintf(fptr_rec,"%lg ",E_lo_90);
    fprintf(fptr_rec,"%lg ",E_hi_90);
    fprintf(fptr_rec,"\n");

    
    A_med   = gsl_stats_median_from_sorted_data   (data->S_pow[i][0], 1, data->Nwave);
    A_lo_50 = gsl_stats_quantile_from_sorted_data (data->S_pow[i][0], 1, data->Nwave, 0.25);
    A_hi_50 = gsl_stats_quantile_from_sorted_data (data->S_pow[i][0], 1, data->Nwave, 0.75);
    A_lo_90 = gsl_stats_quantile_from_sorted_data (data->S_pow[i][0], 1, data->Nwave, 0.05);
    A_hi_90 = gsl_stats_quantile_from_sorted_data (data->S_pow[i][0], 1, data->Nwave, 0.95);
    
    E_med   = gsl_stats_median_from_sorted_data   (data->S_pow[i][1], 1, data->Nwave);
    E_lo_50 = gsl_stats_quantile_from_sorted_data (data->S_pow[i][1], 1, data->Nwave, 0.25);
    E_hi_50 = gsl_stats_quantile_from_sorted_data (data->S_pow[i][1], 1, data->Nwave, 0.75);
    E_lo_90 = gsl_stats_quantile_from_sorted_data (data->S_pow[i][1], 1, data->Nwave, 0.05);
    E_hi_90 = gsl_stats_quantile_from_sorted_data (data->S_pow[i][1], 1, data->Nwave, 0.95);
    
    
    fprintf(fptr_Snf,"%lg ",f);
    fprintf(fptr_Snf,"%lg ",A_med);
    fprintf(fptr_Snf,"%lg ",A_lo_50);
    fprintf(fptr_Snf,"%lg ",A_hi_50);
    fprintf(fptr_Snf,"%lg ",A_lo_90);
    fprintf(fptr_Snf,"%lg ",A_hi_90);
    fprintf(fptr_Snf,"%lg ",E_med);
    fprintf(fptr_Snf,"%lg ",E_lo_50);
    fprintf(fptr_Snf,"%lg ",E_hi_50);
    fprintf(fptr_Snf,"%lg ",E_lo_90);
    fprintf(fptr_Snf,"%lg ",E_hi_90);
    fprintf(fptr_Snf,"\n");

  }
  fclose(fptr_res);
  fclose(fptr_rec);
  fclose(fptr_Snf);
}


