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
#include <sys/stat.h>

#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

#include "LISA.h"
#include "GalacticBinary.h"
#include "GalacticBinaryIO.h"
#include "GalacticBinaryModel.h"

#define FIXME 0

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60
//void printProgress (double percentage)
//{
//  int val = (int) (percentage * 100);
//  int lpad = (int) (percentage * PBWIDTH);
//  int rpad = PBWIDTH - lpad;
//  printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
//  fflush (stdout);
//}

void printProgress (double percentage)
{
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush (stdout);
}

static int checkfile(char filename[])
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
  fprintf(stdout,"       --steps       : number of mcmc steps (10000)        \n");
  fprintf(stdout,"       --samples     : number of frequency bins (2048)     \n");
  fprintf(stdout,"       --segments    : number of data segments (1)         \n");
  fprintf(stdout,"       --start-time  : initial time of segment  (0)        \n");
  fprintf(stdout,"       --gap-time    : duration of data gaps (0)           \n");
  fprintf(stdout,"       --fmin        : minimum frequency                   \n");
  fprintf(stdout,"       --duration    : duration of time segment (62914560) \n");
  fprintf(stdout,"       --noiseseed   : seed for noise RNG                  \n");
  fprintf(stdout,"       --chainseed   : seed for MCMC RNG                   \n");
  fprintf(stdout,"       --chains      : number of parallel chains (20)      \n");
  fprintf(stdout,"       --injseed     : seed for injection parameters       \n");
  fprintf(stdout,"       --inj         : inject signal                       \n");
  fprintf(stdout,"       --data        : strain data file                    \n");
  fprintf(stdout,"       --fix-sky     : pin sky params to injection         \n");
  fprintf(stdout,"       --sky-prior   : use galaxy model for sky prior      \n");
  fprintf(stdout,"       --snr-prior   : use SNR-based amplitude prior       \n");
  fprintf(stdout,"       --known-source: injection is VB (draw orientation)  \n");
  fprintf(stdout,"       --detached    : detached binary(i.e., use Mc prior) \n");
  fprintf(stdout,"       --cheat       : start chain at injection parameters \n");
  fprintf(stdout,"       --update      : use chain as proposal [filename]    \n");
  fprintf(stdout,"       --zero-noise  : data w/out noise realization        \n");
  fprintf(stdout,"       --f-double-dot: include f double dot in model       \n");
  fprintf(stdout,"       --links       : number of links [4->X,6->AE] (6)    \n");
  fprintf(stdout,"       --no-rj       : used fixed dimension                \n");
  fprintf(stdout,"       --fit-gap     : fit for time gaps between segments  \n");
  fprintf(stdout,"       --prior       : sample from prior                   \n");
  fprintf(stdout,"       --debug       : leaner settings for quick running   \n");
  fprintf(stdout,"--\n");
  fprintf(stdout,"EXAMPLE:\n");
  fprintf(stdout,"./gb_mcmc --orbit ../config/OrbitConfig1.txt --verbose --inj ../data/sources/RXJ0806.dat\n");
  fprintf(stdout,"\n");
  fprintf(stdout,"USEFUL TOBS:\n");
  fprintf(stdout,"   1 wk: %.0f\n",62914560./2./52.);
  fprintf(stdout,"   1 mo: %.0f \n",62914560./2./12.);
  fprintf(stdout,"   2 mo: %.0f \n",62914560./2./6.);
  fprintf(stdout,"   1 yr: %.0f \n",62914560./2.);
  fprintf(stdout,"   2 yr: %.0f (default)\n",62914560.);
  fprintf(stdout,"   5 yr: %.0f \n",5.*62914560./2.);
  fprintf(stdout,"  10 yr: %.0f \n",10.*62914560./2.);
  fprintf(stdout,"\n");
  exit(EXIT_FAILURE);
}

void parse(int argc, char **argv, struct Data **data, struct Orbit *orbit, struct Flags *flags, struct Chain *chain, int Nmax)
{
  if(argc==1) print_usage();
  
  //Set defaults
  flags->rj          = 0;
  flags->rj          = 1;
  flags->verbose     = 0;
  flags->NF          = 0;
  flags->zeroNoise   = 0;
  flags->fixSky      = 0;
  flags->skyPrior    = 0;
  flags->snrPrior    = 0;
  flags->cheat       = 0;
  flags->debug       = 0;
  flags->detached    = 0;
  flags->strainData  = 0;
  flags->knownSource = 0;
  flags->NT          = 1;
  flags->orbit       = 0;
  flags->prior       = 0;
  flags->update      = 0;
  flags->NMAX        = Nmax;
  flags->NMCMC       = 10000;
  flags->NBURN       = 10000;
  chain->NP          = 5; //number of proposals
  chain->NC          = 12;//number of chains

  for(int i=0; i<Nmax; i++)
  {
    data[i]->t0   = malloc(sizeof(double)*Nmax);
    data[i]->tgap = malloc(sizeof(double)*Nmax);
    
    for(int j=0; j<Nmax; j++)
    {
      data[i]->t0[j]   = 0.0;
      data[i]->tgap[j] = 0.0;
    }
    
    data[i]->T        = 62914560.0; /* two "mldc years" at 15s sampling */
    data[i]->N        = 1024;
    data[i]->NP       = 8; //default includes fdot
    data[i]->Nchannel = 2; //1=X, 2=AE
    
    data[i]->cseed = 150914+i*Nmax;
    data[i]->nseed = 151226+i*Nmax;
    data[i]->iseed = 151012+i*Nmax;
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
    {"data",      required_argument, 0, 0},
    {"fmin",      required_argument, 0, 0},
    {"links",     required_argument, 0, 0},
    {"update",    required_argument, 0, 0},
    {"steps",     required_argument, 0, 0},
    
    /* These options don’t set a flag.
     We distinguish them by their indices. */
    {"help",        no_argument, 0,'h'},
    {"verbose",     no_argument, 0,'v'},
    {"zero-noise",  no_argument, 0, 0 },
    {"fix-sky",     no_argument, 0, 0 },
    {"sky-prior",   no_argument, 0, 0 },
    {"known-source",no_argument, 0, 0 },
    {"f-double-dot",no_argument, 0, 0 },
    {"detached",    no_argument, 0, 0 },
    {"prior",       no_argument, 0, 0 },
    {"cheat",       no_argument, 0, 0 },
    {"debug",       no_argument, 0, 0 },
    {"no-rj",       no_argument, 0, 0 },
    {"fit-gap",     no_argument, 0, 0 },
    {0, 0, 0, 0}
  };
  
  int opt=0;
  int long_index=0;
  
  //Loop through argv string and pluck out arguments
  struct Data *data_ptr = data[0];
  while ((opt = getopt_long_only(argc, argv,"apl:b:", long_options, &long_index )) != -1)
  {
    switch (opt)
    {
        
      case 0:
        if(strcmp("samples",     long_options[long_index].name) == 0) data_ptr->N       = atoi(optarg);
        if(strcmp("segments",    long_options[long_index].name) == 0) flags->NT         = atoi(optarg);
        if(strcmp("duration",    long_options[long_index].name) == 0) data_ptr->T       = (double)atof(optarg);
        if(strcmp("start-time",  long_options[long_index].name) == 0) data_ptr->t0[0]   = (double)atof(optarg);
        if(strcmp("fmin",        long_options[long_index].name) == 0) data_ptr->fmin    = (double)atof(optarg);
        if(strcmp("gap-time",    long_options[long_index].name) == 0) data_ptr->tgap[0] = (double)atof(optarg);
        if(strcmp("chains",      long_options[long_index].name) == 0) chain->NC         = atoi(optarg);
        if(strcmp("chainseed",   long_options[long_index].name) == 0) data_ptr->cseed   = (long)atoi(optarg);
        if(strcmp("noiseseed",   long_options[long_index].name) == 0) data_ptr->nseed   = (long)atoi(optarg);
        if(strcmp("injseed",     long_options[long_index].name) == 0) data_ptr->iseed   = (long)atoi(optarg);
        if(strcmp("zero-noise",  long_options[long_index].name) == 0) flags->zeroNoise  = 1;
        if(strcmp("fix-sky",     long_options[long_index].name) == 0) flags->fixSky     = 1;
        if(strcmp("sky-prior",   long_options[long_index].name) == 0) flags->skyPrior   = 1;
        if(strcmp("snr-prior",   long_options[long_index].name) == 0) flags->snrPrior   = 1;
        if(strcmp("prior",       long_options[long_index].name) == 0) flags->prior      = 1;
        if(strcmp("f-double-dot",long_options[long_index].name) == 0) data_ptr->NP      = 9;
        if(strcmp("detached",    long_options[long_index].name) == 0) flags->detached   = 1;
        if(strcmp("cheat",       long_options[long_index].name) == 0) flags->cheat      = 1;
        if(strcmp("debug",       long_options[long_index].name) == 0) flags->debug      = 1;
        if(strcmp("no-rj",       long_options[long_index].name) == 0) flags->rj         = 0;
        if(strcmp("fit-gap",     long_options[long_index].name) == 0) flags->gap        = 1;
        if(strcmp("steps",       long_options[long_index].name) == 0)
        {
          flags->NMCMC = atoi(optarg);
          flags->NBURN = flags->NMCMC;
        }
        if(strcmp("known-source",long_options[long_index].name) == 0)
        {
          flags->knownSource = 1;
          flags->fixSky      = 1;
        }
        if(strcmp("data", long_options[long_index].name) == 0)
        {
          checkfile(optarg);
          flags->NF++;
          flags->strainData = 1;
          sprintf(data_ptr->fileName,"%s",optarg);
        }
        if(strcmp("orbit", long_options[long_index].name) == 0)
        {
          checkfile(optarg);
          flags->orbit = 1;
          sprintf(orbit->OrbitFileName,"%s",optarg);
        }
        if(strcmp("inj", long_options[long_index].name) == 0)
        {
          checkfile(optarg);
          sprintf(flags->injFile[flags->NF],"%s",optarg);
          flags->NF++;
          if(flags->NF>Nmax)
          {
            fprintf(stderr,"Requested number of injections is too large (%i/%i)\n",flags->NF,Nmax);
            fprintf(stderr,"Remove at least %i --inj arguments\n",flags->NF-Nmax);
            fprintf(stderr,"Now exiting to system\n");
            exit(1);
          }
        }
        if(strcmp("update", long_options[long_index].name) == 0)
        {
          checkfile(optarg);
          flags->update=1;
          sprintf(flags->cdfFile,"%s",optarg);
          chain->NP++;
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
  if(flags->cheat) flags->NBURN = 0;

  // copy command line args to other data structures
  for(int i=0; i<flags->NF; i++)
  {
    data[i]->NT = flags->NT;
    for(int j=0; j<flags->NT; j++)
    {
      data[i]->t0[j]   = data[0]->t0[0] + j*(data[0]->T + data[0]->tgap[0]);
      data[i]->tgap[j] = data[0]->tgap[0];
    }
    data[i]->T        = data[0]->T;
    data[i]->N        = data[0]->N;
    data[i]->NT       = data[0]->N;
    data[i]->NP       = data[0]->NP;
    data[i]->Nchannel = data[0]->Nchannel;
    
    data[i]->cseed = data[0]->cseed+i*flags->NF;
    data[i]->nseed = data[0]->nseed+i*flags->NF;
    data[i]->iseed = data[0]->iseed+i*flags->NF;
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

  // run looks good to go, make directories and save command line
  mode_t process_mask = umask(0);
  mkdir("chains",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  mkdir("data",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  umask(process_mask);
  
  //Print command line
  FILE *out = fopen("run.sh","w");
  fprintf(out,"#!/bin/sh\n\n");
  for(opt=0; opt<argc; opt++) fprintf(out,"%s ",argv[opt]);
  fprintf(out,"\n\n");
  fclose(out);
  
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
  fprintf(stdout,"  Data start time ..... %.0f \n",data_ptr->t0[0]);
  fprintf(stdout,"  Data duration ....... %.0f \n",data_ptr->T);
  fprintf(stdout,"  Data segments ....... %i   \n",flags->NT);
  fprintf(stdout,"  Data gap duration.....%.0f \n",data_ptr->tgap[0]);
  fprintf(stdout,"  MCMC steps............%i   \n",flags->NMCMC);
  fprintf(stdout,"  MCMC burnin steps.....%i   \n",flags->NBURN);
  fprintf(stdout,"  MCMC chain seed ..... %li  \n",data_ptr->cseed);
  fprintf(stdout,"\n");
  fprintf(stdout,"================= RUN FLAGS ================\n");
  if(flags->verbose)  fprintf(stdout,"  Verbose flag ........ ENABLED \n");
  else                fprintf(stdout,"  Verbose flag ........ DISABLED\n");
  if(flags->NF>0)
  {
    fprintf(stdout,"  Injected sources..... %i\n",flags->NF);
    fprintf(stdout,"     seed ............. %li\n",data_ptr->iseed);
    for(int i=0; i<flags->NF; i++)
    {
      fprintf(stdout,"     source ........... %s\n",flags->injFile[i]);
    }
  }
  else                fprintf(stdout,"  Injection is ........ DISABLED\n");
  if(flags->fixSky)   fprintf(stdout,"  Sky parameters are... DISABLED\n");
  else                fprintf(stdout,"  Sky parameters are... ENABLED\n");
  if(flags->skyPrior) fprintf(stdout,"  Galaxy prior is ..... ENABLED\n");
  else                fprintf(stdout,"  Galaxy prior is ..... DISABLED\n");
  if(flags->snrPrior) fprintf(stdout,"  SNR prior is ........ ENABLED\n");
  else                fprintf(stdout,"  SNR prior is ........ DISABLED\n");
  if(flags->zeroNoise)fprintf(stdout,"  Noise realization is. DISABLED\n");
  else
  {
    fprintf(stdout,"  Noise realization is. ENABLED\n");
    fprintf(stdout,"  Noise seed .......... %li  \n",data_ptr->nseed);
  }
  if(flags->rj)       fprintf(stdout,"  RJMCMC is ........... ENABLED\n");
  else                fprintf(stdout,"  RJMCMC is ........... DISABLED\n");
  if(flags->detached) fprintf(stdout,"  Mchirp prior is...... ENABLED\n");
  else                fprintf(stdout,"  Mchirp prior is...... DISABLED\n");
  fprintf(stdout,"\n");
  fprintf(stdout,"\n");
  
}

void print_chain_files(struct Data *data, struct Model ***model, struct Chain *chain, struct Flags *flags, int step)
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
    for(i=0; i<flags->NF; i++) logL += model[n][i]->logL+model[n][i]->logLnorm;
    fprintf(chain->likelihoodFile,  "%lg ",logL);
    fprintf(chain->temperatureFile, "%lg ",1./chain->temperature[ic]);
  }
  fprintf(chain->likelihoodFile, "\n");
  fprintf(chain->temperatureFile,"\n");
  
  //Always print cold chain
  n = chain->index[0];
  for(i=0; i<flags->NF; i++)
  {
    print_chain_state(data, chain, model[n][i], flags, chain->chainFile[0], step);
    print_noise_state(data, model[n][i], chain->noiseFile[0], step);
  }
  
  //Print sampling parameters
  for(j=0; j<flags->NF; j++)
  {
    int D = model[n][j]->Nlive;
    for(i=0; i<D; i++)
    {
      print_source_params(data,model[n][j]->source[i],chain->parameterFile[0]);
      fprintf(chain->parameterFile[0],"\n");

      if(step>0)
      {
        print_source_params(data,model[n][j]->source[i],chain->dimensionFile[D]);
        fprintf(chain->dimensionFile[D],"\n");
      }
    }
  }
  
  //Print hot chains if verbose flag
  if(flags->verbose)
  {
    for(j=0; j<flags->NF; j++)
    {
      for(ic=1; ic<chain->NC; ic++)
      {
        n = chain->index[ic];
        print_chain_state(data, chain, model[n][j], flags, chain->chainFile[ic], step);
        print_noise_state(data, model[n][j], chain->noiseFile[ic], step);
      }//loop over chains
    }//end loop over
  }//verbose flag
}

void print_chain_state(struct Data *data, struct Chain *chain, struct Model *model, struct Flags *flags, FILE *fptr, int step)
{
  double logL=0.0;
  logL = model->logL+model->logLnorm;
  fprintf(fptr, "%i ",step);
  fprintf(fptr, "%i ",model->Nlive);
  fprintf(fptr, "%lg ",logL);
  for(int j=0; j<flags->NT; j++)fprintf(fptr, "%.12g ",model->t0[j]);
  if(flags->verbose)
  {
    for(int i=0; i<model->Nlive; i++)
    {
      print_source_params(data,model->source[i],fptr);
    }
  }
  fprintf(fptr, "\n");
}

void print_noise_state(struct Data *data, struct Model *model, FILE *fptr, int step)
{
  fprintf(fptr, "%i ",step);
  fprintf(fptr, "%lg ",model->logL+model->logLnorm);

  for(int i=0; i<model->NT; i++)
  {
    switch(data->Nchannel)
    {
      case 1:
        fprintf(fptr, "%lg ", model->noise[i]->etaX);
        break;
      case 2:
        fprintf(fptr, "%lg ", model->noise[i]->etaA);
        fprintf(fptr, "%lg ", model->noise[i]->etaE);
        break;
    }
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

void scan_source_params(struct Data *data, struct Source *source, FILE *fptr)
{
  
  fscanf(fptr,"%lg",&source->f0);
  fscanf(fptr,"%lg",&source->dfdt);
  fscanf(fptr,"%lg",&source->amp);
  fscanf(fptr,"%lg",&source->phi);
  fscanf(fptr,"%lg",&source->costheta);
  fscanf(fptr,"%lg",&source->cosi);
  fscanf(fptr,"%lg",&source->psi);
  fscanf(fptr,"%lg",&source->phi0);
  if(source->NP>8)
    fscanf(fptr,"%lg",&source->d2fdt2);

  //map to parameter names (just to make code readable)
  map_params_to_array(source, source->params, data->T);

}

void save_waveforms(struct Data *data, struct Model *model, int mcmc)
{
  int n_re,n_im;
  double A_re,A_im,E_re,E_im,X_re,X_im,R_re,R_im;
  
  for(int i=0; i<model->NT; i++)
  {
    switch(data->Nchannel)
    {
      case 1:
        for(int n=0; n<data->N; n++)
        {
          n_re = 2*n;
          n_im = n_re++;
          
          X_re = model->tdi[i]->X[n_re];
          X_im = model->tdi[i]->X[n_im];
          
          data->h_rec[n_re][0][i][mcmc] = X_re;
          data->h_rec[n_im][0][i][mcmc] = X_im;
          
          R_re = data->tdi[i]->X[n_re] - X_re;
          R_im = data->tdi[i]->X[n_im] - X_im;
          
          data->h_res[n][0][i][mcmc] = R_re*R_re + R_im*R_im;
          data->h_pow[n][0][i][mcmc] = X_re*X_re + X_im*X_im;
          
          data->S_pow[n][0][i][mcmc] = model->noise[i]->SnX[n];
        }
        break;
      case 2:
        for(int n=0; n<data->N; n++)
        {
          n_re = 2*n;
          n_im = n_re++;
          
          A_re = model->tdi[i]->A[n_re];
          A_im = model->tdi[i]->A[n_im];
          E_re = model->tdi[i]->E[n_re];
          E_im = model->tdi[i]->E[n_im];
          
          data->h_rec[n_re][0][i][mcmc] = A_re;
          data->h_rec[n_im][0][i][mcmc] = A_im;
          data->h_rec[n_re][1][i][mcmc] = E_re;
          data->h_rec[n_im][1][i][mcmc] = E_im;
          
          R_re = data->tdi[i]->A[n_re] - A_re;
          R_im = data->tdi[i]->A[n_im] - A_im;
          data->h_res[n][0][i][mcmc] = R_re*R_re + R_im*R_im;
          
          R_re = data->tdi[i]->E[n_re] - E_re;
          R_im = data->tdi[i]->E[n_im] - E_im;
          data->h_res[n][1][i][mcmc] = R_re*R_re + R_im*R_im;
          
          data->h_pow[n][0][i][mcmc] = A_re*A_re + A_im*A_im;
          data->h_pow[n][1][i][mcmc] = E_re*E_re + E_im*E_im;
          
          data->S_pow[n][0][i][mcmc] = model->noise[i]->SnA[n];
          data->S_pow[n][1][i][mcmc] = model->noise[i]->SnE[n];
        }
        break;
    }
  }
}

void print_waveform(struct Data *data, struct Model *model, FILE *fptr)
{
  for(int n=0; n<data->N; n++)
  {
    int re = 2*n;
    int im = re+1;
    double f = data->fmin + (double)n/data->T;
//    for(int i=0; i<model->NT; i++)
//    {
    int i = 0;
      fprintf(fptr,"%.12g ",f);
      fprintf(fptr,"%.12g ",data->tdi[i]->A[re]*data->tdi[i]->A[re] + data->tdi[i]->A[im]*data->tdi[i]->A[im]);
      fprintf(fptr,"%.12g ",data->tdi[i]->E[re]*data->tdi[i]->E[re] + data->tdi[i]->E[im]*data->tdi[i]->E[im]);
      fprintf(fptr,"%.12g ",model->tdi[i]->A[re]*model->tdi[i]->A[re] + model->tdi[i]->A[im]*model->tdi[i]->A[im]);
      fprintf(fptr,"%.12g ",model->tdi[i]->E[re]*model->tdi[i]->E[re] + model->tdi[i]->E[im]*model->tdi[i]->E[im]);
      fprintf(fptr,"%.12g ",(data->tdi[i]->A[re]-model->tdi[i]->A[re])*(data->tdi[i]->A[re]-model->tdi[i]->A[re]) + (data->tdi[i]->A[im]-model->tdi[i]->A[im])*(data->tdi[i]->A[im]-model->tdi[i]->A[im]) );
      fprintf(fptr,"%.12g ",(data->tdi[i]->E[re]-model->tdi[i]->E[re])*(data->tdi[i]->E[re]-model->tdi[i]->E[re]) + (data->tdi[i]->E[im]-model->tdi[i]->E[im])*(data->tdi[i]->E[im]-model->tdi[i]->E[im]) );
      fprintf(fptr,"\n");
//    }
  }
}

void print_waveform_draw(struct Data **data, struct Model **model, struct Flags *flags)
{
  FILE *fptr;
  char filename[128];
  
  for(int i=0; i<flags->NF; i++)
  {
      sprintf(filename,"data/waveform_draw_%i.dat",i);
      fptr=fopen(filename,"w");
      print_waveform(data[i], model[i], fptr);
      fclose(fptr);
  }
}

void print_waveforms_reconstruction(struct Data *data, int seg)
{
  char filename[1024];
  FILE *fptr_rec;
  FILE *fptr_res;
  FILE *fptr_Snf;

  //sort h reconstructions
  for(int k=0; k<data->NT; k++)
  {
    for(int n=0; n<data->N*2; n++)
    {
      for(int m=0; m<data->Nchannel; m++)
      {
        gsl_sort(data->h_rec[n][m][k],1,data->Nwave);
      }
    }
    
    printf("sort h_res etc.\n");
    for(int n=0; n<data->N; n++)
    {
      for(int m=0; m<data->Nchannel; m++)
      {
        gsl_sort(data->h_res[n][m][k],1,data->Nwave);
        gsl_sort(data->h_pow[n][m][k],1,data->Nwave);
        gsl_sort(data->S_pow[n][m][k],1,data->Nwave);
      }
    }
    
    sprintf(filename,"data/power_reconstruction_t%i_f%i.dat",k,seg);
    fptr_rec=fopen(filename,"w");
    sprintf(filename,"data/power_residual_t%i_f%i.dat",k,seg);
    fptr_res=fopen(filename,"w");
    sprintf(filename,"data/power_noise_t%i_f%i.dat",k,seg);
    fptr_Snf=fopen(filename,"w");
    
    //double X_med,X_lo_50,X_hi_50,X_lo_90,X_hi_90;
    double A_med,A_lo_50,A_hi_50,A_lo_90,A_hi_90;
    double E_med,E_lo_50,E_hi_50,E_lo_90,E_hi_90;
    
    for(int i=0; i<data->N; i++)
    {
      double f = (double)(i+data->qmin)/data->T;
      
      A_med   = gsl_stats_median_from_sorted_data   (data->h_res[i][0][k], 1, data->Nwave);
      A_lo_50 = gsl_stats_quantile_from_sorted_data (data->h_res[i][0][k], 1, data->Nwave, 0.25);
      A_hi_50 = gsl_stats_quantile_from_sorted_data (data->h_res[i][0][k], 1, data->Nwave, 0.75);
      A_lo_90 = gsl_stats_quantile_from_sorted_data (data->h_res[i][0][k], 1, data->Nwave, 0.05);
      A_hi_90 = gsl_stats_quantile_from_sorted_data (data->h_res[i][0][k], 1, data->Nwave, 0.95);
      
      E_med   = gsl_stats_median_from_sorted_data   (data->h_res[i][1][k], 1, data->Nwave);
      E_lo_50 = gsl_stats_quantile_from_sorted_data (data->h_res[i][1][k], 1, data->Nwave, 0.25);
      E_hi_50 = gsl_stats_quantile_from_sorted_data (data->h_res[i][1][k], 1, data->Nwave, 0.75);
      E_lo_90 = gsl_stats_quantile_from_sorted_data (data->h_res[i][1][k], 1, data->Nwave, 0.05);
      E_hi_90 = gsl_stats_quantile_from_sorted_data (data->h_res[i][1][k], 1, data->Nwave, 0.95);
      
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
      
      A_med   = gsl_stats_median_from_sorted_data   (data->h_pow[i][0][k], 1, data->Nwave);
      A_lo_50 = gsl_stats_quantile_from_sorted_data (data->h_pow[i][0][k], 1, data->Nwave, 0.25);
      A_hi_50 = gsl_stats_quantile_from_sorted_data (data->h_pow[i][0][k], 1, data->Nwave, 0.75);
      A_lo_90 = gsl_stats_quantile_from_sorted_data (data->h_pow[i][0][k], 1, data->Nwave, 0.05);
      A_hi_90 = gsl_stats_quantile_from_sorted_data (data->h_pow[i][0][k], 1, data->Nwave, 0.95);
      
      E_med   = gsl_stats_median_from_sorted_data   (data->h_pow[i][1][k], 1, data->Nwave);
      E_lo_50 = gsl_stats_quantile_from_sorted_data (data->h_pow[i][1][k], 1, data->Nwave, 0.25);
      E_hi_50 = gsl_stats_quantile_from_sorted_data (data->h_pow[i][1][k], 1, data->Nwave, 0.75);
      E_lo_90 = gsl_stats_quantile_from_sorted_data (data->h_pow[i][1][k], 1, data->Nwave, 0.05);
      E_hi_90 = gsl_stats_quantile_from_sorted_data (data->h_pow[i][1][k], 1, data->Nwave, 0.95);
      
      
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
      
      
      A_med   = gsl_stats_median_from_sorted_data   (data->S_pow[i][0][k], 1, data->Nwave);
      A_lo_50 = gsl_stats_quantile_from_sorted_data (data->S_pow[i][0][k], 1, data->Nwave, 0.25);
      A_hi_50 = gsl_stats_quantile_from_sorted_data (data->S_pow[i][0][k], 1, data->Nwave, 0.75);
      A_lo_90 = gsl_stats_quantile_from_sorted_data (data->S_pow[i][0][k], 1, data->Nwave, 0.05);
      A_hi_90 = gsl_stats_quantile_from_sorted_data (data->S_pow[i][0][k], 1, data->Nwave, 0.95);
      
      E_med   = gsl_stats_median_from_sorted_data   (data->S_pow[i][1][k], 1, data->Nwave);
      E_lo_50 = gsl_stats_quantile_from_sorted_data (data->S_pow[i][1][k], 1, data->Nwave, 0.25);
      E_hi_50 = gsl_stats_quantile_from_sorted_data (data->S_pow[i][1][k], 1, data->Nwave, 0.75);
      E_lo_90 = gsl_stats_quantile_from_sorted_data (data->S_pow[i][1][k], 1, data->Nwave, 0.05);
      E_hi_90 = gsl_stats_quantile_from_sorted_data (data->S_pow[i][1][k], 1, data->Nwave, 0.95);
      
      
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
}


