/*
*  Copyright (C) 2019 Tyson B. Littenberg (MSFC-ST12), Neil J. Cornish
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



#include <math.h>
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
#include "GalacticBinaryWaveform.h"
#include "gitversion.h"

#define FIXME 0

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60


void printProgress (double percentage)
{
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush (stdout);
}


void print_version(FILE *fptr)
{
  fprintf(fptr, "\n");
  fprintf(fptr, "============== GBMCMC Version: =============\n\n");
  //fprintf(fptr, "  Git remote origin: %s\n", GIT_URL);
  //fprintf(fptr, "  Git version: %s\n", GIT_VER);
  fprintf(fptr, "  Git commit: %s\n", gitversion);
  //fprintf(fptr, "  Git commit author: %s\n",GIT_AUTHOR);
  //fprintf(fptr, "  Git commit date: %s\n", GIT_DATE);
}

void print_run_settings(int argc, char **argv, struct Data *data_ptr, struct Orbit *orbit, struct Flags *flags, FILE *fptr)
{
  fprintf(fptr,"\n");
  fprintf(fptr,"=============== RUN SETTINGS ===============\n");
  fprintf(fptr,"\n");
  fprintf(fptr,"  Command Line: ");
  for(int opt=0; opt<argc; opt++) fprintf(fptr,"%s ",argv[opt]);
  fprintf(fptr,"\n");
  fprintf(fptr,"\n");
  switch(flags->orbit)
  {
    case 0:
      fprintf(fptr,"  Orbit model is ...... EccentricInclined \n");
      break;
    case 1:
      fprintf(fptr,"  Orbit file .......... %s   \n",orbit->OrbitFileName);
      break;
  }
  fprintf(fptr,"  Data channels ........");
  switch(data_ptr->Nchannel)
  {
    case 1:
      fprintf(fptr,"X\n");
      break;
    case 2:
      fprintf(fptr,"AE\n");
      break;
  }
  fprintf(fptr,"  Data sample size .... %i   \n",data_ptr->N);
  fprintf(fptr,"  Data padding size ... %i   \n",data_ptr->qpad);
  fprintf(fptr,"  Data start time ..... %.0f \n",data_ptr->t0[0]);
  fprintf(fptr,"  Data start frequency. %.16g\n",data_ptr->fmin);
  fprintf(fptr,"  Data duration ....... %.0f \n",data_ptr->T);
  fprintf(fptr,"  Data segments ....... %i   \n",flags->NT);
  fprintf(fptr,"  Data gap duration.....%.0f \n",data_ptr->tgap[0]);
  fprintf(fptr,"  Data format is........%s   \n",data_ptr->format);
  fprintf(fptr,"  Max # of sources......%i   \n",flags->DMAX);
  fprintf(fptr,"  MCMC steps............%i   \n",flags->NMCMC);
  fprintf(fptr,"  MCMC burnin steps.....%i   \n",flags->NBURN);
  fprintf(fptr,"  MCMC chain seed ..... %li  \n",data_ptr->cseed);
  fprintf(fptr,"\n");
  fprintf(fptr,"================= RUN FLAGS ================\n");
  if(flags->verbose)  fprintf(fptr,"  Verbose flag ........ ENABLED \n");
  else                fprintf(fptr,"  Verbose flag ........ DISABLED\n");
  if(flags->quiet)    fprintf(fptr,"  Quiet flag .......... ENABLED \n");
  else                fprintf(fptr,"  Quiet flag .......... DISABLED\n");
  if(flags->NINJ>0)
  {
    fprintf(fptr,"  Injected sources..... %i\n",flags->NINJ);
    fprintf(fptr,"     seed ............. %li\n",data_ptr->iseed);
    for(int i=0; i<flags->NINJ; i++)
    {
      fprintf(fptr,"     source ........... %s\n",flags->injFile[i]);
    }
  }
  else                   fprintf(fptr,"  Injection is ........ DISABLED\n");
  if(flags->fixSky)      fprintf(fptr,"  Sky parameters are... DISABLED\n");
  if(flags->fixFreq)     fprintf(fptr,"  Freq paramseters are. DISABLED\n");
  else                   fprintf(fptr,"  Sky parameters are... ENABLED\n");
  if(flags->calibration) fprintf(fptr,"  Calibration is....... ENABLED\n");
  else                   fprintf(fptr,"  Calibration is....... DISABLED\n");
  if(flags->galaxyPrior) fprintf(fptr,"  Galaxy prior is ..... ENABLED\n");
  else                   fprintf(fptr,"  Galaxy prior is ..... DISABLED\n");
  if(flags->snrPrior)    fprintf(fptr,"  SNR prior is ........ ENABLED\n");
  else                   fprintf(fptr,"  SNR prior is ........ DISABLED\n");
  if(flags->simNoise)
  {
    fprintf(fptr,"  Noise simulation is.. ENABLED\n");
    fprintf(fptr,"  Noise seed .......... %li  \n",data_ptr->nseed);
  }
  else                fprintf(fptr,"  Noise simulation is.. DISABLED\n");
  if(flags->rj)       fprintf(fptr,"  RJMCMC is ........... ENABLED\n");
  else                fprintf(fptr,"  RJMCMC is ........... DISABLED\n");
  if(flags->detached) fprintf(fptr,"  Mchirp prior is...... ENABLED\n");
  else                fprintf(fptr,"  Mchirp prior is...... DISABLED\n");
  fprintf(fptr,"\n");
  fprintf(fptr,"\n");
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

void print_usage()
{
  fprintf(stdout,"\n");
  fprintf(stdout,"=============== GBMCMC Usage: ============== \n");
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
  fprintf(stdout,"       --links       : number of links [4->X,6->AE] (6)    \n");
  fprintf(stdout,"       --frac-freq   : fractional frequency data (phase)   \n");
  fprintf(stdout,"\n");

  //Data
  fprintf(stdout,"       =========== Data =========== \n");
  fprintf(stdout,"       --data        : strain data file                    \n");
  fprintf(stdout,"       --samples     : number of frequency bins (2048)     \n");
  fprintf(stdout,"       --padding     : number of bins padded on segment (0)\n");
  fprintf(stdout,"       --segments    : number of data segments (1)         \n");
  fprintf(stdout,"       --start-time  : initial time of segment  (0)        \n");
  fprintf(stdout,"       --gap-time    : duration of data gaps (0)           \n");
  fprintf(stdout,"       --fmin        : minimum frequency                   \n");
  fprintf(stdout,"       --duration    : duration of time segment (62914560) \n");
  fprintf(stdout,"       --sim-noise   : data w/out noise realization        \n");
  fprintf(stdout,"       --conf-noise  : include model for confusion noise   \n");
  fprintf(stdout,"       --noiseseed   : seed for noise RNG                  \n");
  fprintf(stdout,"       --inj         : inject signal                       \n");
  fprintf(stdout,"       --injseed     : seed for injection parameters       \n");
  fprintf(stdout,"       --catalog     : list of known sources               \n");
  fprintf(stdout,"\n");

  //Chain
  fprintf(stdout,"       ========== Chains ========== \n");
  fprintf(stdout,"       --steps       : number of mcmc steps (100000)       \n");
  fprintf(stdout,"       --chainseed   : seed for MCMC RNG                   \n");
  fprintf(stdout,"       --chains      : number of parallel chains (20)      \n");
  fprintf(stdout,"       --no-burnin   : skip burn in steps                  \n");
  fprintf(stdout,"\n");

  //Model
  fprintf(stdout,"       ========== Model =========== \n");
  fprintf(stdout,"       --sources     : maximum number of sources (10)      \n");
  fprintf(stdout,"       --cheat       : start chain at injection parameters \n");
  fprintf(stdout,"       --f-double-dot: include f double dot in model       \n");
  fprintf(stdout,"       --prior       : sample from prior                   \n");
  fprintf(stdout,"       --no-rj       : used fixed dimension                \n");
  fprintf(stdout,"       --calibration : marginalize over calibration errors \n");
  fprintf(stdout,"       --fit-gap     : fit for time gaps between segments  \n");
  fprintf(stdout,"\n");

  //Priors & Proposals
  fprintf(stdout,"       ==== Priors & Proposals ==== \n");
  fprintf(stdout,"       --fix-sky     : pin sky params to injection         \n");
  fprintf(stdout,"       --galaxy-prior: use galaxy model for sky prior      \n");
  fprintf(stdout,"       --snr-prior   : use SNR-based amplitude prior       \n");
  fprintf(stdout,"       --em-prior    : update prior ranges from other obs  \n");
  fprintf(stdout,"       --known-source: injection is VB (draw orientation)  \n");
  fprintf(stdout,"       --detached    : detached binary(i.e., use Mc prior) \n");
  fprintf(stdout,"       --update      : use chain as proposal [filename]    \n");
  fprintf(stdout,"       --update-cov  : use cov mtrx proposal [filename]    \n");
  fprintf(stdout,"\n");
  
  //Misc.
  fprintf(stdout,"       =========== Misc =========== \n");
  fprintf(stdout,"       --match-in1      : input paramaters for overlap [filename] \n");
  fprintf(stdout,"       --match-in2      : output match values [filename] \n");
  
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
  print_LISA_ASCII_art(stdout);
  print_version(stdout);

  if(argc==1) print_usage();
  
  int DMAX_default = 10;
  
  //Set defaults
  flags->calibration = 0;
  flags->rj          = 1;
  flags->verbose     = 0;
  flags->quiet       = 0;
  flags->NDATA       = 1;
  flags->NINJ        = 0;
  flags->simNoise    = 0;
  flags->confNoise   = 0;
  flags->fixSky      = 0;
  flags->fixFreq     = 0;
  flags->galaxyPrior = 0;
  flags->snrPrior    = 0;
  flags->emPrior     = 0;
  flags->cheat       = 0;
  flags->burnin      = 1;
  flags->debug       = 0;
  flags->detached    = 0;
  flags->strainData  = 0;
  flags->knownSource = 0;
  flags->catalog     = 0;
  flags->NT          = 1;
  flags->orbit       = 0;
  flags->prior       = 0;
  flags->update      = 0;
  flags->updateCov   = 0;
  flags->match       = 0;
  flags->resume      = 0;
  flags->DMAX        = DMAX_default;
  flags->NMCMC       = 100000;
  flags->NBURN       = 100000;
  chain->NP          = 9; //number of proposals
  chain->NC          = 12;//number of chains
  
  for(int i=0; i<Nmax; i++)
  {
    /*
     default data format is 'phase' 
     optional support for 'frequency' a la LDCs
    */
    sprintf(data[i]->format,"phase");
    
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
    data[i]->DMAX     = DMAX_default;//maximum number of sources
    data[i]->qpad     = 0;
    
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
    {"padding",   required_argument, 0, 0},
    {"duration",  required_argument, 0, 0},
    {"segments",  required_argument, 0, 0},
    {"sources",   required_argument, 0, 0},
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
    {"update-cov",required_argument, 0, 0},
    {"match-in1", required_argument, 0, 0},
    {"match-in2", required_argument, 0, 0},
    {"steps",     required_argument, 0, 0},
    {"em-prior",  required_argument, 0, 0},
    {"catalog",   required_argument, 0, 0},
    
    /* These options don’t set a flag.
     We distinguish them by their indices. */
    {"help",        no_argument, 0,'h'},
    {"verbose",     no_argument, 0,'v'},
    {"quiet",       no_argument, 0,'q'},
    {"debug",       no_argument, 0,'d'},
    {"resume",      no_argument, 0, 0 },
    {"sim-noise",   no_argument, 0, 0 },
    {"conf-noise",  no_argument, 0, 0 },
    {"frac-freq",   no_argument, 0, 0 },
    {"fix-sky",     no_argument, 0, 0 },
    {"fix-freq",    no_argument, 0, 0 },
    {"galaxy-prior",no_argument, 0, 0 },
    {"snr-prior",   no_argument, 0, 0 },
    {"known-source",no_argument, 0, 0 },
    {"f-double-dot",no_argument, 0, 0 },
    {"detached",    no_argument, 0, 0 },
    {"prior",       no_argument, 0, 0 },
    {"cheat",       no_argument, 0, 0 },
    {"no-burnin",   no_argument, 0, 0 },
    {"no-rj",       no_argument, 0, 0 },
    {"fit-gap",     no_argument, 0, 0 },
    {"calibration", no_argument, 0, 0 },
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
        if(strcmp("padding",     long_options[long_index].name) == 0) data_ptr->qpad    = atoi(optarg);
        if(strcmp("segments",    long_options[long_index].name) == 0) flags->NT         = atoi(optarg);
        if(strcmp("duration",    long_options[long_index].name) == 0) data_ptr->T       = (double)atof(optarg);
        if(strcmp("start-time",  long_options[long_index].name) == 0) data_ptr->t0[0]   = (double)atof(optarg);
        if(strcmp("fmin",        long_options[long_index].name) == 0) sscanf(optarg, "%lg", &data_ptr->fmin);
        if(strcmp("gap-time",    long_options[long_index].name) == 0) data_ptr->tgap[0] = (double)atof(optarg);
        if(strcmp("chains",      long_options[long_index].name) == 0) chain->NC         = atoi(optarg);
        if(strcmp("chainseed",   long_options[long_index].name) == 0) data_ptr->cseed   = (long)atoi(optarg);
        if(strcmp("noiseseed",   long_options[long_index].name) == 0) data_ptr->nseed   = (long)atoi(optarg);
        if(strcmp("injseed",     long_options[long_index].name) == 0) data_ptr->iseed   = (long)atoi(optarg);
        if(strcmp("sim-noise",   long_options[long_index].name) == 0) flags->simNoise   = 1;
        if(strcmp("conf-noise",  long_options[long_index].name) == 0) flags->confNoise  = 1;
        if(strcmp("fix-sky",     long_options[long_index].name) == 0) flags->fixSky     = 1;
        if(strcmp("fix-freq",    long_options[long_index].name) == 0) flags->fixFreq    = 1;
        if(strcmp("galaxy-prior",long_options[long_index].name) == 0) flags->galaxyPrior= 1;
        if(strcmp("snr-prior",   long_options[long_index].name) == 0) flags->snrPrior   = 1;
        if(strcmp("prior",       long_options[long_index].name) == 0) flags->prior      = 1;
        if(strcmp("f-double-dot",long_options[long_index].name) == 0) data_ptr->NP      = 9;
        if(strcmp("detached",    long_options[long_index].name) == 0) flags->detached   = 1;
        if(strcmp("cheat",       long_options[long_index].name) == 0) flags->cheat      = 1;
        if(strcmp("no-burnin",   long_options[long_index].name) == 0) flags->burnin     = 0;
        if(strcmp("no-rj",       long_options[long_index].name) == 0) flags->rj         = 0;
        if(strcmp("fit-gap",     long_options[long_index].name) == 0) flags->gap        = 1;
        if(strcmp("calibration", long_options[long_index].name) == 0) flags->calibration= 1;
        if(strcmp("resume",      long_options[long_index].name) == 0)
          flags->resume=1;
        if(strcmp("sources",     long_options[long_index].name) == 0)
        {
          data_ptr->DMAX    = atoi(optarg);
          flags->DMAX       = atoi(optarg);
        }
        if(strcmp("em-prior",    long_options[long_index].name) == 0)
        {
          flags->emPrior = 1;
          sprintf(flags->pdfFile,"%s",optarg);
        }
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
        if(strcmp("catalog",long_options[long_index].name) == 0)
        {
            flags->catalog = 1;
            sprintf(flags->catalogFile,"%s",optarg);
        }
        if(strcmp("data", long_options[long_index].name) == 0)
        {
          checkfile(optarg);
          //flags->NDATA++;
          flags->strainData = 1;
          sprintf(data_ptr->fileName,"%s",optarg);
        }
        if(strcmp("frac-freq",   long_options[long_index].name) == 0)
        {
          for(int i=0; i<Nmax; i++) sprintf(data[i]->format,"frequency");
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
          sprintf(flags->injFile[flags->NINJ],"%s",optarg);
          flags->NINJ++;
          if(flags->NINJ>Nmax)
          {
            fprintf(stderr,"WARNING: Requested number of injections is too large (%i/%i)\n",flags->NINJ,Nmax);
            fprintf(stderr,"Should you remove at least %i --inj arguments?\n",flags->NINJ-Nmax);
            fprintf(stderr,"Now exiting to system\n");
            exit(1);
          }
        }
        if(strcmp("update", long_options[long_index].name) == 0)
        {
          checkfile(optarg);
          flags->update=1;
          sprintf(flags->cdfFile,"%s",optarg);
        }
        if(strcmp("update-cov", long_options[long_index].name) == 0)
        {
          checkfile(optarg);
          flags->updateCov=1;
          sprintf(flags->covFile,"%s",optarg);
        }
        if(strcmp("match-in1", long_options[long_index].name) == 0)
        {
            checkfile(optarg);
            flags->match=1;
            sprintf(flags->matchInfile1,"%s",optarg);
        }
            if(strcmp("match-in2", long_options[long_index].name) == 0)
        {
            sprintf(flags->matchInfile2,"%s",optarg);
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
      case 'd' : flags->debug = 1;
        break;
      case 'h' :
        print_usage();
        exit(EXIT_FAILURE);
        break;
      case 'v' : flags->verbose = 1;
        break;
      case 'q' : flags->quiet = 1;
        break;
      default: print_usage();
        exit(EXIT_FAILURE);
    }
  }
  if(flags->cheat || !flags->burnin) flags->NBURN = 0;

  if(flags->verbose && flags->quiet)
  {
    fprintf(stderr,"--verbose and --quiet flags are in conflict\n");
    exit(1);
  }
  
  //pad data
  data[0]->N += 2*data[0]->qpad;
  data[0]->fmin -= data[0]->qpad/data[0]->T;

    
  // copy command line args to other data structures
  for(int i=0; i<flags->NDATA; i++)
  {
    data[i]->NT = flags->NT;
    for(int j=0; j<flags->NT; j++)
    {
      data[i]->t0[j]   = data[0]->t0[0] + j*(data[0]->T + data[0]->tgap[0]);
      data[i]->tgap[j] = data[0]->tgap[0];
    }
    data[i]->T        = data[0]->T;
    data[i]->qpad     = data[0]->qpad;
    data[i]->N        = data[0]->N;
    data[i]->NT       = data[0]->N;
    data[i]->NP       = data[0]->NP;
    data[i]->Nchannel = data[0]->Nchannel;
    data[i]->DMAX     = data[0]->DMAX;
      
    
    data[i]->cseed = data[0]->cseed+i*flags->NDATA;
    data[i]->nseed = data[0]->nseed+i*flags->NDATA;
    data[i]->iseed = data[0]->iseed+i*flags->NDATA;

    //map fmin to nearest bin
    data[i]->fmin = floor(data[i]->fmin*data[i]->T)/data[i]->T;
    
    //calculate helper quantities for likelihood normalizations
    data[i]->logfmin   = log(data[i]->fmin);
    data[i]->sum_log_f = 0.0;
    for(int n=0; n<data[i]->N; n++)
    {
      data[i]->sum_log_f += log(data[i]->fmin + (double)n/data[i]->T);
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

  // run looks good to go, make directories and save command line
  mode_t process_mask = umask(0);
  mkdir("checkpoint",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  mkdir("chains",    S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  mkdir("data",      S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  umask(process_mask);
  
  //Print command line
  FILE *out = fopen("run.sh","w");
  fprintf(out,"#!/bin/sh\n\n");
  for(opt=0; opt<argc; opt++) fprintf(out,"%s ",argv[opt]);
  fprintf(out,"\n\n");
  fclose(out);
  
  //Print version control
  FILE *runlog = fopen("gb_mcmc.log","w");
  print_version(runlog);
  
  //Report on set parameters
  print_run_settings(argc, argv, data_ptr, orbit, flags, stdout);
  print_run_settings(argc, argv, data_ptr, orbit, flags, runlog);
  
  fclose(runlog);

}

void save_chain_state(struct Data **data, struct Model ***model, struct Chain *chain, struct Flags *flags, int step)
{
  char filename[128];
  FILE *stateFile;
  for(int ic=0; ic<chain->NC; ic++)
  {
    sprintf(filename,"checkpoint/chain_state_%i.dat",ic);
    stateFile = fopen(filename,"w");

    int n = chain->index[ic];
    
    fprintf(stateFile,"%.12g\n",chain->logLmax);
    
    for(int j=0; j<flags->NDATA; j++)
    {
      print_chain_state(data[j], chain, model[n][j], flags, stateFile, step);
      print_noise_state(data[j], model[n][j], stateFile, step);
      if(flags->calibration)
        print_calibration_state(data[j], model[n][j], stateFile, step);
      
      int D = model[n][j]->Nlive;
      for(int i=0; i<D; i++)
      {
        print_source_params(data[j],model[n][j]->source[i],stateFile);
        fprintf(stateFile,"\n");
      }
    }
    
    fclose(stateFile);
  }
}

void restore_chain_state(struct Orbit *orbit, struct Data **data, struct Model ***model, struct Chain *chain, struct Flags *flags, int *step)
{
  char filename[128];
  FILE *stateFile;
  chain->logLmax=0.0;
  for(int ic=0; ic<chain->NC; ic++)
  {
    sprintf(filename,"checkpoint/chain_state_%i.dat",ic);
    stateFile = fopen(filename,"r");

    int n = chain->index[ic];
    
    fscanf(stateFile,"%lg",&chain->logLmax);
    
    for(int j=0; j<flags->NDATA; j++)
    {
      scan_chain_state(data[j], chain, model[n][j], flags, stateFile, step);
      scan_noise_state(data[j], model[n][j], stateFile, step);
      if(flags->calibration)
        scan_calibration_state(data[j], model[n][j], stateFile, step);
      
      int D = model[n][j]->Nlive;
      for(int i=0; i<D; i++)
      {
        scan_source_params(data[j],model[n][j]->source[i], stateFile);
        galactic_binary_fisher(orbit, data[j], model[n][j]->source[i], data[j]->noise[0]);
      }
      
      generate_noise_model(data[j], model[n][j]);
      generate_signal_model(orbit, data[j], model[n][j], -1);
      
      if(!flags->prior)
      {
        model[n][j]->logL = gaussian_log_likelihood(orbit, data[j], model[n][j]);
        model[n][j]->logLnorm = gaussian_log_likelihood_constant_norm(data[j], model[n][j]);
      }
      else model[n][j]->logL = model[n][j]->logLnorm = 0.0;

    }
    
    fclose(stateFile);
  }
}

void print_chain_files(struct Data *data, struct Model ***model, struct Chain *chain, struct Flags *flags, int step)
{
  int i,j,n,ic;
  
  //Print logL & temperature chains
  if(!flags->quiet)
  {
    fprintf(chain->likelihoodFile,  "%i ",step);
    fprintf(chain->temperatureFile, "%i ",step);
    double logL;
    for(ic=0; ic<chain->NC; ic++)
    {
      n = chain->index[ic];
      logL=0.0;
      for(i=0; i<flags->NDATA; i++) logL += model[n][i]->logL+model[n][i]->logLnorm;
      fprintf(chain->likelihoodFile,  "%lg ",logL);
      fprintf(chain->temperatureFile, "%lg ",1./chain->temperature[ic]);
    }
    fprintf(chain->likelihoodFile, "\n");
    fprintf(chain->temperatureFile,"\n");
  }
  
  //Print cold chains
  n = chain->index[0];
  for(i=0; i<flags->NDATA; i++)
  {
    print_chain_state(data, chain, model[n][i], flags, chain->chainFile[0], step);
    if(!flags->quiet || step>0)
      print_noise_state(data, model[n][i], chain->noiseFile[0], step);
    if(flags->calibration)
      print_calibration_state(data, model[n][i], chain->calibrationFile[0], step);
  }
    if(flags->verbose)
    {
        fflush(chain->chainFile[0]);
        fflush(chain->noiseFile[0]);
      if(flags->calibration) fflush(chain->calibrationFile[0]);
    }
  
  //Print sampling parameters
  for(j=0; j<flags->NDATA; j++)
  {
    int D = model[n][j]->Nlive;
    for(i=0; i<D; i++)
    {
      if(!flags->quiet || step>0)
      {
        print_source_params(data,model[n][j]->source[i],chain->parameterFile[0]);
        fprintf(chain->parameterFile[0],"\n");
        if(flags->verbose)fflush(chain->parameterFile[0]);
      }
      if(step>0)
      {
        print_source_params(data,model[n][j]->source[i],chain->dimensionFile[D]);
        fprintf(chain->dimensionFile[D],"\n");
      }
    }
  }
  
  //Print calibration parameters
  for(j=0; j<flags->NDATA; j++)
  {
    
  }
  
  //Print hot chains if verbose flag
  if(flags->verbose)
  {
    for(j=0; j<flags->NDATA; j++)
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

void scan_chain_state(struct Data *data, struct Chain *chain, struct Model *model, struct Flags *flags, FILE *fptr, int *step)
{
  fscanf(fptr, "%i",step);
  fscanf(fptr, "%i",&model->Nlive);
  fscanf(fptr, "%lg",&model->logL);
  fscanf(fptr, "%lg",&model->logLnorm);
  for(int j=0; j<flags->NT; j++)fscanf(fptr, "%lg",&model->t0[j]);
  if(flags->verbose)
  {
    for(int i=0; i<model->Nlive; i++)
    {
      scan_source_params(data,model->source[i],fptr);
    }
  }
}

void print_chain_state(struct Data *data, struct Chain *chain, struct Model *model, struct Flags *flags, FILE *fptr, int step)
{
  fprintf(fptr, "%i ",step);
  fprintf(fptr, "%i ",model->Nlive);
  fprintf(fptr, "%lg ",model->logL);
  fprintf(fptr, "%lg ",model->logLnorm);
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

void scan_calibration_state(struct Data *data, struct Model *model, FILE *fptr, int *step)
{
  fscanf(fptr, "%i %lg %lg",step, &model->logL,&model->logLnorm);
  
  for(int i=0; i<model->NT; i++)
  {
    switch(data->Nchannel)
    {
      case 1:
        fscanf(fptr, "%lg", &model->calibration[i]->dampX);
        fscanf(fptr, "%lg", &model->calibration[i]->dphiX);
        break;
      case 2:
        fscanf(fptr, "%lg", &model->calibration[i]->dampA);
        fscanf(fptr, "%lg", &model->calibration[i]->dphiA);
        fscanf(fptr, "%lg", &model->calibration[i]->dampE);
        fscanf(fptr, "%lg", &model->calibration[i]->dphiE);
        break;
    }
  }
}
void print_calibration_state(struct Data *data, struct Model *model, FILE *fptr, int step)
{
  fprintf(fptr, "%i ",step);
  fprintf(fptr, "%lg %lg ",model->logL, model->logLnorm);
  
  for(int i=0; i<model->NT; i++)
  {
    switch(data->Nchannel)
    {
      case 1:
        fprintf(fptr, "%lg ", model->calibration[i]->dampX);
        fprintf(fptr, "%lg ", model->calibration[i]->dphiX);
        break;
      case 2:
        fprintf(fptr, "%lg ", model->calibration[i]->dampA);
        fprintf(fptr, "%lg ", model->calibration[i]->dphiA);
        fprintf(fptr, "%lg ", model->calibration[i]->dampE);
        fprintf(fptr, "%lg ", model->calibration[i]->dphiE);
        break;
    }
  }
  fprintf(fptr, "\n");
}

void scan_noise_state(struct Data *data, struct Model *model, FILE *fptr, int *step)
{
  fscanf(fptr, "%i ",step);
  fscanf(fptr, "%lg %lg ", &model->logL, &model->logLnorm);

  for(int i=0; i<model->NT; i++)
  {
    switch(data->Nchannel)
    {
      case 1:
        fscanf(fptr, "%lg", &model->noise[i]->etaX);
        break;
      case 2:
        fscanf(fptr, "%lg", &model->noise[i]->etaA);
        fscanf(fptr, "%lg", &model->noise[i]->etaE);
        break;
    }
  }
}

void print_noise_state(struct Data *data, struct Model *model, FILE *fptr, int step)
{
  fprintf(fptr, "%i ",step);
  fprintf(fptr, "%lg %lg ",model->logL, model->logLnorm);

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
  
  fprintf(fptr,"%.16g ",source->f0);
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

          data->h_res[n_re][0][i][mcmc] = R_re;
          data->h_res[n_im][0][i][mcmc] = R_im;

          data->r_pow[n][0][i][mcmc] = R_re*R_re + R_im*R_im;
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

          data->h_res[n_re][0][i][mcmc] = R_re;
          data->h_res[n_im][0][i][mcmc] = R_im;

          data->r_pow[n][0][i][mcmc] = R_re*R_re + R_im*R_im;
          
          R_re = data->tdi[i]->E[n_re] - E_re;
          R_im = data->tdi[i]->E[n_im] - E_im;

          data->h_res[n_re][1][i][mcmc] = R_re;
          data->h_res[n_im][1][i][mcmc] = R_im;

          data->r_pow[n][1][i][mcmc] = R_re*R_re + R_im*R_im;
          
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
  
  int N = 1;
  if(flags->NINJ>1) N = flags->NINJ;
  for(int i=0; i<N; i++)
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
  FILE *fptr_var;
  FILE *fptr_Snf;

  //get variance of residual
  double ***res_var = malloc(data->N*sizeof(double **));
  for(int n=0; n<data->N; n++)
  {
    res_var[n] = malloc(data->Nchannel*sizeof(double *));
    for(int m=0; m<data->Nchannel; m++)
    {
      res_var[n][m] =  malloc(data->NT*sizeof(double));
      for(int k=0; k<data->NT; k++) res_var[n][m][k] = 0.0;
    }
  }

  
  for(int k=0; k<data->NT; k++)
  {
    for(int n=0; n<data->N*2; n++)
    {
      for(int m=0; m<data->Nchannel; m++)
      {
        gsl_sort(data->h_rec[n][m][k],1,data->Nwave);
      }
    }

    for(int n=0; n<data->N; n++)
    {
      for(int m=0; m<data->Nchannel; m++)
      {
        gsl_sort(data->r_pow[n][m][k],1,data->Nwave);
        gsl_sort(data->h_pow[n][m][k],1,data->Nwave);
        gsl_sort(data->S_pow[n][m][k],1,data->Nwave);
        res_var[n][m][k] = gsl_stats_variance(data->h_rec[2*n][m][k], 1, data->Nwave)+gsl_stats_variance(data->h_rec[2*n+1][m][k], 1, data->Nwave);
      }
    }
    
    sprintf(filename,"data/power_reconstruction_t%i_f%i.dat",k,seg);
    fptr_rec=fopen(filename,"w");
    sprintf(filename,"data/power_residual_t%i_f%i.dat",k,seg);
    fptr_res=fopen(filename,"w");
    sprintf(filename,"data/power_noise_t%i_f%i.dat",k,seg);
    fptr_Snf=fopen(filename,"w");
    sprintf(filename,"data/variance_residual_t%i_f%i.dat",k,seg);
    fptr_var=fopen(filename,"w");

    //double X_med,X_lo_50,X_hi_50,X_lo_90,X_hi_90;
    double A_med,A_lo_50,A_hi_50,A_lo_90,A_hi_90;
    double E_med,E_lo_50,E_hi_50,E_lo_90,E_hi_90;
    
    for(int i=0; i<data->N; i++)
    {
      double f = (double)(i+data->qmin)/data->T;
      fprintf(fptr_var,"%.12g %.12g %.12g\n",f,res_var[i][0][k],res_var[i][1][k]);
      
      A_med   = gsl_stats_median_from_sorted_data   (data->r_pow[i][0][k], 1, data->Nwave);
      A_lo_50 = gsl_stats_quantile_from_sorted_data (data->r_pow[i][0][k], 1, data->Nwave, 0.25);
      A_hi_50 = gsl_stats_quantile_from_sorted_data (data->r_pow[i][0][k], 1, data->Nwave, 0.75);
      A_lo_90 = gsl_stats_quantile_from_sorted_data (data->r_pow[i][0][k], 1, data->Nwave, 0.05);
      A_hi_90 = gsl_stats_quantile_from_sorted_data (data->r_pow[i][0][k], 1, data->Nwave, 0.95);
      
      E_med   = gsl_stats_median_from_sorted_data   (data->r_pow[i][1][k], 1, data->Nwave);
      E_lo_50 = gsl_stats_quantile_from_sorted_data (data->r_pow[i][1][k], 1, data->Nwave, 0.25);
      E_hi_50 = gsl_stats_quantile_from_sorted_data (data->r_pow[i][1][k], 1, data->Nwave, 0.75);
      E_lo_90 = gsl_stats_quantile_from_sorted_data (data->r_pow[i][1][k], 1, data->Nwave, 0.05);
      E_hi_90 = gsl_stats_quantile_from_sorted_data (data->r_pow[i][1][k], 1, data->Nwave, 0.95);
      
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
    
    fclose(fptr_var);
    fclose(fptr_res);
    fclose(fptr_rec);
    fclose(fptr_Snf);
  }
}


