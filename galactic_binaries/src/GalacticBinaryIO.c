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

#define FIXME 0

void print_usage()
{
  fprintf(stdout,"\n");
  fprintf(stdout,"Usage: \n");
  fprintf(stdout,"REQUIRED:\n");
  fprintf(stdout,"\n");
  fprintf(stdout,"OPTIONAL:\n");
  fprintf(stdout,"  -h | --help        : print help message and exit         \n");
  fprintf(stdout,"  -v | --verbose     : enable verbose output               \n");
  fprintf(stdout,"       --orbit       : orbit ephemerides file (2.5 GM MLDC)\n");
  fprintf(stdout,"       --samples     : number of frequency bins (2048)     \n");
  fprintf(stdout,"       --segments    : number of data segments (1)         \n");
  fprintf(stdout,"       --start-time  : initial time of segment  (0)        \n");
  fprintf(stdout,"       --gap-time    : duration of data gaps (0)           \n");
  fprintf(stdout,"       --duration    : duration of time segment (62914560) \n");
  fprintf(stdout,"       --noiseseed   : seed for noise RNG                  \n");
  fprintf(stdout,"       --chainseed   : seed for MCMC RNG                   \n");
  fprintf(stdout,"       --chains      : number of parallel chains (20)      \n");
  fprintf(stdout,"       --injseed     : seed for injection parameters       \n");
  fprintf(stdout,"       --inj         : inject signal                       \n");
  fprintf(stdout,"       --fix-sky     : pin sky params to injection         \n");
  fprintf(stdout,"       --known-source: injection is VB (draw orientation)  \n");
  fprintf(stdout,"       --cheat       : start chain at injection parameters \n");
  fprintf(stdout,"       --zero-noise  : data w/out noise realization        \n");
  fprintf(stdout,"       --f-double-dot: include f double dot in model       \n");
  fprintf(stdout,"       --links       : number of links [4->X,6->AE] (6)    \n");
  fprintf(stdout,"       --prior       : sample from prior                   \n");
  fprintf(stdout,"--\n");
  fprintf(stdout,"EXAMPLE:\n");
  fprintf(stdout,"./gb_mcmc --orbit ../config/OrbitConfig1.txt --verbose --inj ../data/sources/RXJ0806.dat\n");
  fprintf(stdout,"\n");
  exit(EXIT_FAILURE);
}

void parse(int argc, char **argv, struct Data ***data, struct Orbit *orbit, struct Flags *flags, struct Chain *chain, int Nmax)
{
  //Set defaults
  flags->verbose     = 0;
  flags->injection   = 0;
  flags->zeroNoise   = 0;
  flags->fixSky      = 0;
  flags->cheat       = 0;
  flags->knownSource = 0;
  flags->segment     = 1;
  flags->orbit       = 0;
  flags->prior       = 0;
  chain->NC          = 20;

  for(int i=0; i<Nmax; i++)
  {
    for(int j=0; j<Nmax; j++)
    {
      data[i][j]->t0       = 0.0;
      data[i][j]->tgap     = 0.0;
      data[i][j]->T        = 62914560.0; /* two "mldc years" at 15s sampling */
      data[i][j]->N        = 1024;
      data[i][j]->NP       = 8; //default includes fdot
      data[i][j]->Nchannel = 2; //1=X, 2=AE
      
      data[i][j]->cseed = 150914+i*Nmax+j;
      data[i][j]->nseed = 151226+i*Nmax+j;
      data[i][j]->iseed = 151012+i*Nmax+j;
    }
  }
  
  flags->injFile = malloc(10*sizeof(char *));
  for(int n=0; n<10; n++) flags->injFile[n] = malloc(1024*sizeof(char));
  
  //if(argc==1) print_usage();
  
  //Specifying the expected options
  static struct option long_options[] =
  {
    /* These options set a flag. */
    {"samples",   required_argument, 0, 0},
    {"duration",  required_argument, 0, 0},
    {"segments",  required_argument, 0, 0},
    {"start-time",required_argument, 0, 0},
    {"gap-time",  required_argument, 0, 0},
    {"orbit",     required_argument, 0, 0},
    {"chains",    required_argument, 0, 0},
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
    {"fix-sky", no_argument,       0,   0  },
    {"known-source",no_argument,   0,   0  },
    {"f-double-dot",no_argument,   0,   0  },
    {"prior",   no_argument,       0,   0  },
    {"cheat",   no_argument,       0,   0  },
    {0,         0,                 0,   0  }
  };
  
  int opt=0;
  int long_index=0;
  
  //Print command line
  FILE *out = fopen("run.sh","w");
  fprintf(out,"#!/bin/sh\n\n");
  for(opt=0; opt<argc; opt++) fprintf(out,"%s ",argv[opt]);
  fprintf(out,"\n\n");
  fclose(out);
  
  //Loop through argv string and pluck out arguments
  struct Data *data_ptr = data[0][0];
  while ((opt = getopt_long_only(argc, argv,"apl:b:", long_options, &long_index )) != -1)
  {
    switch (opt)
    {
        
      case 0:
        if(strcmp("samples",     long_options[long_index].name) == 0) data_ptr->N       = atoi(optarg);
        if(strcmp("segments",    long_options[long_index].name) == 0) flags->segment    = atoi(optarg);
        if(strcmp("duration",    long_options[long_index].name) == 0) data_ptr->T       = (double)atof(optarg);
        if(strcmp("start-time",  long_options[long_index].name) == 0) data_ptr->t0      = (double)atof(optarg);
        if(strcmp("gap-time",    long_options[long_index].name) == 0) data_ptr->tgap    = (double)atof(optarg);
        if(strcmp("chains",      long_options[long_index].name) == 0) chain->NC         = atoi(optarg);
        if(strcmp("chainseed",   long_options[long_index].name) == 0) data_ptr->cseed   = (long)atoi(optarg);
        if(strcmp("noiseseed",   long_options[long_index].name) == 0) data_ptr->nseed   = (long)atoi(optarg);
        if(strcmp("injseed",     long_options[long_index].name) == 0) data_ptr->iseed   = (long)atoi(optarg);
        if(strcmp("zero-noise",  long_options[long_index].name) == 0) flags->zeroNoise  = 1;
        if(strcmp("fix-sky",     long_options[long_index].name) == 0) flags->fixSky     = 1;
        if(strcmp("prior",       long_options[long_index].name) == 0) flags->prior      = 1;
        if(strcmp("cheat",       long_options[long_index].name) == 0) flags->cheat      = 1;
        if(strcmp("f-double-dot",long_options[long_index].name) == 0) data_ptr->NP      = 9;
        if(strcmp("known-source",long_options[long_index].name) == 0)
        {
          flags->knownSource = 1;
          flags->fixSky      = 1;
        }
        if(strcmp("orbit",       long_options[long_index].name) == 0)
        {
          flags->orbit = 1;
          sprintf(orbit->OrbitFileName,"%s",optarg);
        }
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
        if(strcmp("links",long_options[long_index].name) == 0)
        {
          int Nlinks = (int)atoi(optarg);
          switch(Nlinks)
          {
            case 4:
              data_ptr->Nchannel=1;
              break;
            case 6:
              data_ptr->Nchannel=2;
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
  for(int i=0; i<flags->injection; i++)
  {
    for(int j=0; j<flags->segment; j++)
    {
      data[i][j]->t0       = data[0][0]->t0 + j*(data[0][0]->T + data[0][0]->tgap);
      data[i][j]->tgap     = data[0][0]->tgap;
      data[i][j]->T        = data[0][0]->T;
      data[i][j]->N        = data[0][0]->N;
      data[i][j]->NP       = data[0][0]->NP;
      data[i][j]->Nchannel = data[0][0]->Nchannel;
      
      data[i][j]->cseed = data[0][0]->cseed+i*flags->injection + j;
      data[i][j]->nseed = data[0][0]->nseed+i*flags->injection + j;
      data[i][j]->iseed = data[0][0]->iseed+i*flags->injection + j;      
    }
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
  switch(flags->orbit)
  {
    case 0:
      fprintf(stdout,"  Orbit model is ...... EccentricInclined \n");
      break;
    case 1:
      fprintf(stdout,"  Orbit file .......... %s   \n",orbit->OrbitFileName);
      break;
  }
  fprintf(stdout,"  Data channels ........");
  switch(data_ptr->Nchannel)
  {
    case 1:
      fprintf(stdout,"X\n");
      break;
    case 2:
      fprintf(stdout,"AE\n");
      break;
  }
  fprintf(stdout,"  Data sample size .... %i   \n",data_ptr->N);
  fprintf(stdout,"  Data start time ..... %.0f \n",data_ptr->t0);
  fprintf(stdout,"  Data duration ....... %.0f \n",data_ptr->T);
  fprintf(stdout,"  Data segments ....... %i   \n",flags->segment);
  fprintf(stdout,"  Data gap duration.....%.0f \n",data_ptr->tgap);
  fprintf(stdout,"  MCMC chain seed ..... %li  \n",data_ptr->cseed);
  fprintf(stdout,"\n");
  fprintf(stdout,"================= RUN FLAGS ================\n");
  if(flags->verbose)  fprintf(stdout,"  Verbose flag ........ ENABLED \n");
  else                fprintf(stdout,"  Verbose flag ........ DISABLED\n");
  if(flags->injection>0)
  {
    fprintf(stdout,"  Injected sources..... %i\n",flags->injection);
    fprintf(stdout,"     seed ............. %li\n",data_ptr->iseed);
    for(int i=0; i<flags->injection; i++)
    {
      fprintf(stdout,"     source ........... %s\n",flags->injFile[i]);
    }
  }
  else                fprintf(stdout,"  Injection is ........ DISABLED\n");
  if(flags->fixSky)   fprintf(stdout,"  Sky parameters are... DISABLED\n");
  else                fprintf(stdout,"  Sky parameters are... ENABLED\n");
  if(flags->zeroNoise)fprintf(stdout,"  Noise realization is. DISABLED\n");
  else
  {
    fprintf(stdout,"  Noise realization is. ENABLED\n");
    fprintf(stdout,"  Noise seed .......... %li  \n",data_ptr->nseed);
  }
  fprintf(stdout,"\n");
  fprintf(stdout,"\n");
  
}

void print_chain_files(struct Data *data, struct Model ****model, struct Chain *chain, struct Flags *flags, int step)
{
  int i,j,n,ic;
  
  //Always print logL & temperature chains
  fprintf(chain->likelihoodFile,  "%i ",step);
  fprintf(chain->temperatureFile, "%i ",step);
  double logL;
  for(ic=0; ic<chain->NC; ic++)
  {
    n = chain->index[ic];
    logL=0.0;
    for(i=0; i<flags->injection; i++) for(j=0; j<flags->segment; j++)logL += model[n][i][j]->logL+model[n][i][j]->logLnorm;
    fprintf(chain->likelihoodFile,  "%lg ",logL);
    fprintf(chain->temperatureFile, "%lg ",1./chain->temperature[ic]);
  }
  fprintf(chain->likelihoodFile, "\n");
  fprintf(chain->temperatureFile,"\n");
  
  //Always print cold chain
  n = chain->index[0];
  for(i=0; i<flags->injection; i++)
  {
    print_chain_state(data, chain, model[n][i], flags, chain->parameterFile[0], step);
    print_noise_state(data, model[n][i][FIXME], chain->noiseFile[0], step);
  }
  
  //Print hot chains if verbose flag
  if(flags->verbose)
  {
    for(ic=1; ic<chain->NC; ic++)
    {
      n = chain->index[ic];
      print_chain_state(data, chain, model[n][0], flags, chain->parameterFile[ic], step);
      print_noise_state(data, model[n][0][FIXME], chain->noiseFile[ic], step);
    }//loop over chains
  }//verbose flag
}

void print_chain_state(struct Data *data, struct Chain *chain, struct Model **model, struct Flags *flags, FILE *fptr, int step)
{
  double logL=0.0;
  for(int i=0; i<flags->segment; i++) logL += model[i]->logL+model[i]->logLnorm;
  for(int i=0; i<model[0]->Nlive; i++)
  {
    fprintf(fptr, "%i ",step);
    fprintf(fptr, "%lg ",logL);
    for(int j=0; j<flags->segment; j++)fprintf(fptr, "%.12g ",model[j]->t0);
    //fprintf(fptr, "%lg ",data->tgap);
    print_source_params(data,model[0]->source[i],fptr);
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
  fprintf(fptr,"%.12g ",source->amp);
  fprintf(fptr,"%.12g ",source->phi);
  fprintf(fptr,"%.12g ",source->costheta);
  fprintf(fptr,"%.12g ",source->cosi);
  fprintf(fptr,"%.12g ",source->psi);
  fprintf(fptr,"%.12g ",source->phi0);
  if(source->NP>8)
    fprintf(fptr,"%.12g ",source->d2fdt2);
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

void print_waveform(struct Data *data, struct Model *model, FILE *fptr)
{
  for(int n=0; n<data->N; n++)
  {
    int re = 2*n;
    int im = re+1;
    double f = data->fmin + (double)n/data->T;
    fprintf(fptr,"%.12g ",f);
    fprintf(fptr,"%.12g ",data->tdi->A[re]);
    fprintf(fptr,"%.12g ",data->tdi->A[im]);
    fprintf(fptr,"%.12g ",data->tdi->E[re]);
    fprintf(fptr,"%.12g ",data->tdi->E[im]);
    fprintf(fptr,"%.12g ",model->tdi->A[re]);
    fprintf(fptr,"%.12g ",model->tdi->A[im]);
    fprintf(fptr,"%.12g ",model->tdi->E[re]);
    fprintf(fptr,"%.12g ",model->tdi->E[im]);
    fprintf(fptr,"\n");
  }

}

void print_waveform_draw(struct Data ***data, struct Model ***model, struct Flags *flags)
{
  FILE *fptr;
  char filename[128];
  
  for(int i=0; i<flags->injection; i++)
  {
    for(int j=0; j<flags->segment; j++)
    {
      sprintf(filename,"waveform_draw_%i_%i.dat",i,j);
      fptr=fopen(filename,"w");
      print_waveform(data[i][j], model[i][j], fptr);
      fclose(fptr);
    }
  }
}

void print_waveforms_reconstruction(struct Data *data, int seg)
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
  
  char filename[1024];
  sprintf(filename,"power_reconstruction_%i.dat",seg);
  FILE *fptr_rec=fopen(filename,"w");
  sprintf(filename,"power_residual_%i.dat",seg);
  FILE *fptr_res=fopen(filename,"w");
  sprintf(filename,"power_noise_%i.dat",seg);
  FILE *fptr_Snf=fopen(filename,"w");
  
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
    
    fprintf(fptr_res,"%.12g ",f);
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
    
    
    fprintf(fptr_rec,"%.12g ",f);
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
    
    
    fprintf(fptr_Snf,"%.12g ",f);
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


