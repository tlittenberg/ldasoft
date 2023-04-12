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
#include <omp.h>

#include <sys/stat.h>

#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

#include <util.h>
#include <LISA.h>

#include "GalacticBinary.h"
#include "GalacticBinaryIO.h"
#include "GalacticBinaryMath.h"
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
    fprintf(stdout, "\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush (stdout);
}


void print_version(FILE *fptr)
{
    fprintf(fptr, "\n");
    fprintf(fptr, "============== LDASOFT Version: =============\n\n");
    //fprintf(fptr, "  Git remote origin: %s\n", GIT_URL);
    //fprintf(fptr, "  Git version: %s\n", GIT_VER);
    fprintf(fptr, "  Git commit: %s\n", GITVERSION);
    //fprintf(fptr, "  Git commit author: %s\n",GIT_AUTHOR);
    //fprintf(fptr, "  Git commit date: %s\n", GIT_DATE);
}

void setup_run_directories(struct Flags *flags, struct Data *data, struct Chain *chain)
{
    
    pathprintf(data->dataDir,"%s/data",flags->runDir);
    pathprintf(chain->chainDir,"%s/chains",flags->runDir);
    pathprintf(chain->chkptDir,"%s/checkpoint",flags->runDir);

    mkdir(flags->runDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(data->dataDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain->chainDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain->chkptDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

}

void print_gb_catalog_script(struct Flags *flags, struct Data *data, struct Orbit *orbit)
{
    
    //back out original input f & N
    int samples = data->N - 2*data->qpad;
    double fmin = data->fmin + data->qpad/data->T;
    
    char filename[PATH_BUFSIZE];
    pathprintf(filename,"%s/example_gb_catalog.sh",flags->runDir);
    FILE *fptr = fopen(filename,"w");
    
    fprintf(fptr,"#!/bin/sh\n\n");
    fprintf(fptr,"if [ \"$#\" -ne 1 ]; then\n");
    fprintf(fptr,"\t echo \"You must enter model dimension\"\n");
    fprintf(fptr,"fi\n\n");
    
    fprintf(fptr,"gb_catalog ");
    
    //Required
    if(flags->NVB>0) fprintf(fptr,"--fmin $2 ");
    else fprintf(fptr,"--fmin %.12g ", fmin);
    fprintf(fptr,"--samples %i ", samples);
    fprintf(fptr,"--padding %i ",data->qpad);
    fprintf(fptr,"--duration %f ",data->T);
    fprintf(fptr,"--start-time %f ",data->t0[0]);
    fprintf(fptr,"--sources $1 --chain-file chains/dimension_chain.dat.$1 ");
    
    //Optional
    if(strcmp(data->format,"frequency")==0)
        fprintf(fptr,"--frac-freq ");
    if(strcmp(data->format,"sangria")==0)
        fprintf(fptr,"--sangria ");
    if(flags->orbit)
        fprintf(fptr,"--orbit %s ",orbit->OrbitFileName);
    if(data->NP==9)
        fprintf(fptr,"--f-double-dot ");
    if(data->Nchannel==1)
        fprintf(fptr,"--links 4 ");
    
    fprintf(fptr,"\n\n");
    
    //Recommendations
    fprintf(fptr,"# Consider including the following options:\n");
    fprintf(fptr,"#\t--match       : match threshold for waveforms (0.8)\n");
    fprintf(fptr,"#\t--noise-file  : reconstructed noise model\n");
    fprintf(fptr,"#\t\t e.g., data/power_noise_0.dat\n");
    fprintf(fptr,"#\t--catalog     : list of known sources\n");
    fprintf(fptr,"#\t--Tcatalog    : observing time of previous catalog\n");
    
    fclose(fptr);
}

void print_run_settings(int argc, char **argv, struct Data *data, struct Orbit *orbit, struct Flags *flags, FILE *fptr)
{
    fprintf(fptr,"\n");
    fprintf(fptr,"=============== RUN SETTINGS ===============\n");
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
    switch(data->Nchannel)
    {
        case 1:
            fprintf(fptr,"X\n");
            break;
        case 2:
            fprintf(fptr,"AE\n");
            break;
    }
    fprintf(fptr,"  Data sample size .... %i   \n",data->N);
    fprintf(fptr,"  Data padding size ... %i   \n",data->qpad);
    fprintf(fptr,"  Data start time ..... %.0f \n",data->t0[0]);
    fprintf(fptr,"  Data start frequency. %.16g\n",data->fmin);
    fprintf(fptr,"  Data duration ....... %.0f \n",data->T);
    fprintf(fptr,"  Data epochs ......... %i   \n",flags->NT);
    fprintf(fptr,"  Data gap duration.....%.0f \n",data->tgap[0]);
    fprintf(fptr,"  Data format is........%s   \n",data->format);
    fprintf(fptr,"  Max # of sources......%i   \n",flags->DMAX-1);
    fprintf(fptr,"  MCMC steps............%i   \n",flags->NMCMC);
    fprintf(fptr,"  MCMC burnin steps.....%i   \n",flags->NBURN);
    fprintf(fptr,"  MCMC chain seed ..... %li  \n",data->cseed);
    fprintf(fptr,"  Number of threads ... %i   \n",flags->threads);
    fprintf(fptr,"  Run Directory is .... %s\n",flags->runDir);
    fprintf(fptr,"\n");
    fprintf(fptr,"================= RUN FLAGS ================\n");
    if(flags->verbose)  fprintf(fptr,"  Verbose flag ........ ENABLED \n");
    else                fprintf(fptr,"  Verbose flag ........ DISABLED\n");
    if(flags->quiet)    fprintf(fptr,"  Quiet flag .......... ENABLED \n");
    else                fprintf(fptr,"  Quiet flag .......... DISABLED\n");
    if(flags->NINJ>0)
    {
        fprintf(fptr,"  Injected sources..... %i\n",flags->NINJ);
        fprintf(fptr,"     seed ............. %li\n",data->iseed);
        for(int i=0; i<flags->NINJ; i++)
        {
            fprintf(fptr,"     source ........... %s\n",flags->injFile[i]);
        }
    }
    else                   fprintf(fptr,"  Injection is ........ DISABLED\n");
    if(flags->fixSky)      fprintf(fptr,"  Sky parameters are... DISABLED\n");
    if(flags->fixFreq)     fprintf(fptr,"  Freq parameters are.. DISABLED\n");
    else                   fprintf(fptr,"  Freq parameters are.. ENABLED\n");
    if(flags->fixFdot)     fprintf(fptr,"  Fdot parameters are.. DISABLED\n");
    else                   fprintf(fptr,"  Fdot parameters are.. ENABLED\n");
    if(flags->calibration) fprintf(fptr,"  Calibration is....... ENABLED\n");
    else                   fprintf(fptr,"  Calibration is....... DISABLED\n");
    if(flags->galaxyPrior) fprintf(fptr,"  Galaxy prior is ..... ENABLED\n");
    else                   fprintf(fptr,"  Galaxy prior is ..... DISABLED\n");
    if(flags->snrPrior)    fprintf(fptr,"  SNR prior is ........ ENABLED\n");
    else                   fprintf(fptr,"  SNR prior is ........ DISABLED\n");
    if(flags->simNoise)
    {
        fprintf(fptr,"  Noise simulation is.. ENABLED\n");
        fprintf(fptr,"  Noise seed .......... %li  \n",data->nseed);
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
    fprintf(stdout,"       --sangria     : use LDC Sangria TDI conventions     \n");
    fprintf(stdout,"\n");
    
    //Data
    fprintf(stdout,"       =========== Data =========== \n");
    fprintf(stdout,"       --data        : strain data file (ASCII)            \n");
    fprintf(stdout,"       --h5-data     : strain data file (HDF5)             \n");
    fprintf(stdout,"       --psd         : psd data file (ASCII)               \n");
    fprintf(stdout,"       --samples     : number of frequency bins (2048)     \n");
    fprintf(stdout,"       --samples_max : max size of segment (2048)          \n");
    fprintf(stdout,"       --padding     : number of bins padded on segment (0)\n");
    fprintf(stdout,"       --epochs      : number of time segments (1)         \n");
    fprintf(stdout,"       --start-time  : initial time of epoch  (0)          \n");
    fprintf(stdout,"       --gap-time    : duration of data gaps (0)           \n");
    fprintf(stdout,"       --fmin        : minimum frequency                   \n");
    fprintf(stdout,"       --duration    : duration of epoch (62914560)        \n");
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
    fprintf(stdout,"       --resume      : restart from checkpoint             \n");
    fprintf(stdout,"       --threads     : number of parallel threads (max)    \n");
    fprintf(stdout,"\n");
    
    //Model
    fprintf(stdout,"       ========== Model =========== \n");
    fprintf(stdout,"       --sources     : maximum number of sources (10)      \n");
    fprintf(stdout,"       --cheat       : start chain at injection parameters \n");
    fprintf(stdout,"       --f-double-dot: include f double dot in model       \n");
    fprintf(stdout,"       --prior       : sample from prior                   \n");
    fprintf(stdout,"       --no-rj       : used fixed dimension                \n");
    fprintf(stdout,"       --calibration : marginalize over calibration errors \n");
    fprintf(stdout,"       --fit-gap     : fit for time gaps between epochs    \n");
    fprintf(stdout,"\n");
    
    //Priors & Proposals
    fprintf(stdout,"       ==== Priors & Proposals ==== \n");
    fprintf(stdout,"       --fix-sky     : pin sky params to injection         \n");
    fprintf(stdout,"       --fix-freq    : pin frequency to injection          \n");
    fprintf(stdout,"       --fix-fdot    : pin f && fdot to injection          \n");
    fprintf(stdout,"       --galaxy-prior: use galaxy model for sky prior      \n");
    fprintf(stdout,"       --no-snr-prior: don't use SNR-based amplitude prior \n");
    fprintf(stdout,"       --em-prior    : update prior ranges from other obs  \n");
    fprintf(stdout,"       --known-source: injection is VB (draw orientation)  \n");
    fprintf(stdout,"       --detached    : detached binary(i.e., use Mc prior) \n");
    fprintf(stdout,"       --update      : gmm of posterior as prior [filename]\n");
    fprintf(stdout,"       --update-cov  : use cov matrix proposal [filename]  \n");
    fprintf(stdout,"\n");
    
    //Misc.
    fprintf(stdout,"       =========== Misc =========== \n");
    fprintf(stdout,"       --rundir      : top level run directory ['./']\n");
    fprintf(stdout,"       --match-in1   : input paramaters for overlap [filename] \n");
    fprintf(stdout,"       --match-in2   : output match values [filename] \n");
    
    fprintf(stdout,"--\n");
    fprintf(stdout,"EXAMPLE:\n");
    fprintf(stdout,"gb_mcmc --inj [path to]/ldasoft/gbmcmc/etc/sources/precision/PrecisionSource_0.txt\n");
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
    exit(0);
}

void parse(int argc, char **argv, struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Chain *chain, int Nmax)
{
    
    int DMAX_default = 10;
    int NmaxFlag = 0; //flag if Nmax is set at command line
    
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
    flags->fixFdot     = 0;
    flags->galaxyPrior = 0;
    flags->snrPrior    = 1;
    flags->emPrior     = 0;
    flags->cheat       = 0;
    flags->burnin      = 1;
    flags->debug       = 0;
    flags->detached    = 0;
    flags->strainData  = 0;
    flags->hdf5Data    = 0;
    flags->psd         = 0;
    flags->knownSource = 0;
    flags->catalog     = 0;
    flags->NT          = 1;
    flags->gap         = 0;
    flags->orbit       = 0;
    flags->prior       = 0;
    flags->update      = 0;
    flags->updateCov   = 0;
    flags->match       = 0;
    flags->resume      = 0;
    flags->DMAX        = DMAX_default;
    flags->NMCMC       = 100000;
    flags->NBURN       = 100000;
    flags->threads     = omp_get_max_threads();
    sprintf(flags->runDir,"./");
    chain->NP          = 9; //number of proposals
    chain->NC          = 12;//number of chains
        
    /*
     default data format is 'phase'
     optional support for 'frequency' a la LDCs
     */
    sprintf(data->format,"phase");
    
    data->t0   = calloc(Nmax,sizeof(double));
    data->tgap = calloc(Nmax,sizeof(double));
    
    data->T        = 62914560.0; /* two "mldc years" at 15s sampling */
    data->sqT      = sqrt(data->T);
    data->N        = 1024;
    data->NP       = 8; //default includes fdot
    data->Nchannel = 2; //1=X, 2=AE
    data->DMAX     = DMAX_default;//maximum number of sources
    data->qpad     = 0;
    
    data->cseed = 150914;
    data->nseed = 151226;
    data->iseed = 151012;
    
    
    flags->injFile = malloc(10*sizeof(char *));
    for(int n=0; n<10; n++) flags->injFile[n] = malloc(1024*sizeof(char));
        
    //Specifying the expected options
    static struct option long_options[] =
    {
        /* These options set a flag. */
        {"samples",    required_argument, 0, 0},
        {"samples_max",required_argument, 0, 0},
        {"padding",    required_argument, 0, 0},
        {"duration",   required_argument, 0, 0},
        {"epochs",     required_argument, 0, 0},
        {"sources",    required_argument, 0, 0},
        {"start-time", required_argument, 0, 0},
        {"gap-time",   required_argument, 0, 0},
        {"orbit",      required_argument, 0, 0},
        {"chains",     required_argument, 0, 0},
        {"chainseed",  required_argument, 0, 0},
        {"noiseseed",  required_argument, 0, 0},
        {"injseed",    required_argument, 0, 0},
        {"inj",        required_argument, 0, 0},
        {"data",       required_argument, 0, 0},
        {"h5-data",    required_argument, 0, 0},
        {"psd",        required_argument, 0, 0},
        {"fmin",       required_argument, 0, 0},
        {"links",      required_argument, 0, 0},
        {"update-cov", required_argument, 0, 0},
        {"match-in1",  required_argument, 0, 0},
        {"match-in2",  required_argument, 0, 0},
        {"steps",      required_argument, 0, 0},
        {"em-prior",   required_argument, 0, 0},
        {"catalog",    required_argument, 0, 0},
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
        {"frac-freq",   no_argument, 0, 0 },
        {"sangria",     no_argument, 0, 0 },
        {"update",      no_argument, 0, 0 },
        {"fix-sky",     no_argument, 0, 0 },
        {"fix-freq",    no_argument, 0, 0 },
        {"fix-fdot",    no_argument, 0, 0 },
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
    
    //Print command line
//    char filename[MAXSTRINGSIZE];
//    sprintf(filename,"run.sh");
//    FILE *out = fopen(filename,"w");
//    fprintf(out,"#!/bin/sh\n\n");
//    for(opt=0; opt<argc; opt++) fprintf(out,"%s ",argv[opt]);
//    fprintf(out,"\n\n");
//    fclose(out);

    //Loop through argv string and pluck out arguments
    while ((opt = getopt_long_only(argc, argv,"apl:b:", long_options, &long_index )) != -1)
    {
        switch (opt)
        {
                
            case 0:
                if(strcmp("samples",     long_options[long_index].name) == 0) data->N           = atoi(optarg);
                if(strcmp("padding",     long_options[long_index].name) == 0) data->qpad        = atoi(optarg);
                if(strcmp("epochs",      long_options[long_index].name) == 0) flags->NT         = atoi(optarg);
                if(strcmp("start-time",  long_options[long_index].name) == 0) data->t0[0]       = (double)atof(optarg);
                if(strcmp("fmin",        long_options[long_index].name) == 0) sscanf(optarg, "%lg", &data->fmin);
                if(strcmp("gap-time",    long_options[long_index].name) == 0) data->tgap[0]     = (double)atof(optarg);
                if(strcmp("chains",      long_options[long_index].name) == 0) chain->NC         = atoi(optarg);
                if(strcmp("chainseed",   long_options[long_index].name) == 0) data->cseed       = (long)atoi(optarg);
                if(strcmp("noiseseed",   long_options[long_index].name) == 0) data->nseed       = (long)atoi(optarg);
                if(strcmp("injseed",     long_options[long_index].name) == 0) data->iseed       = (long)atoi(optarg);
                if(strcmp("sim-noise",   long_options[long_index].name) == 0) flags->simNoise   = 1;
                if(strcmp("conf-noise",  long_options[long_index].name) == 0) flags->confNoise  = 1;
                if(strcmp("fix-sky",     long_options[long_index].name) == 0) flags->fixSky     = 1;
                if(strcmp("fix-freq",    long_options[long_index].name) == 0) flags->fixFreq    = 1;
                if(strcmp("update",      long_options[long_index].name) == 0) flags->update     = 1;
                if(strcmp("galaxy-prior",long_options[long_index].name) == 0) flags->galaxyPrior= 1;
                if(strcmp("no-snr-prior",long_options[long_index].name) == 0) flags->snrPrior   = 0;
                if(strcmp("prior",       long_options[long_index].name) == 0) flags->prior      = 1;
                if(strcmp("f-double-dot",long_options[long_index].name) == 0) data->NP          = 9;
                if(strcmp("detached",    long_options[long_index].name) == 0) flags->detached   = 1;
                if(strcmp("cheat",       long_options[long_index].name) == 0) flags->cheat      = 1;
                if(strcmp("no-burnin",   long_options[long_index].name) == 0) flags->burnin     = 0;
                if(strcmp("no-rj",       long_options[long_index].name) == 0) flags->rj         = 0;
                if(strcmp("fit-gap",     long_options[long_index].name) == 0) flags->gap        = 1;
                if(strcmp("calibration", long_options[long_index].name) == 0) flags->calibration= 1;
                if(strcmp("resume",      long_options[long_index].name) == 0) flags->resume     = 1;
                if(strcmp("threads",     long_options[long_index].name) == 0) flags->threads    = atoi(optarg);
                if(strcmp("rundir",      long_options[long_index].name) == 0) strcpy(flags->runDir,optarg);
                if(strcmp("duration",    long_options[long_index].name) == 0)
                {   data->T   = (double)atof(optarg);
                    data->sqT = sqrt(data->T);
                }
                if(strcmp("sources",     long_options[long_index].name) == 0)
                {
                    data->DMAX  = atoi(optarg)+1;
                    flags->DMAX = atoi(optarg)+1;
                }
                if(strcmp("em-prior",    long_options[long_index].name) == 0)
                {
                    flags->emPrior = 1;
                    sprintf(flags->pdfFile,"%s",optarg);
                }
                if(strcmp("samples_max", long_options[long_index].name) == 0)
                {
                    NmaxFlag = 1;
                    data->Nmax = atoi(optarg);
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
                    flags->fixFreq     = 1;
                }
                if(strcmp("fix-fdot",long_options[long_index].name) == 0)
                {
                    flags->fixFdot = 1;
                    flags->fixFreq = 1; //if you are fixing fdot you ought to be fixing f as well.
                }
                if(strcmp("catalog",long_options[long_index].name) == 0)
                {
                    checkfile(optarg);
                    flags->catalog = 1;
                    sprintf(flags->catalogFile,"%s",optarg);
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
                if(strcmp("frac-freq",   long_options[long_index].name) == 0)
                {
                    for(int i=0; i<Nmax; i++) sprintf(data->format,"frequency");
                }
                if(strcmp("sangria",   long_options[long_index].name) == 0)
                {
                    for(int i=0; i<Nmax; i++) sprintf(data->format,"sangria");
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
                            data->Nchannel=1;
                            break;
                        case 6:
                            data->Nchannel=2;
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
    
    //pad data
    if(!NmaxFlag) data->Nmax = data->N;
    data->N += 2*data->qpad;
    data->Nmax += 2*data->qpad;
    data->fmin -= data->qpad/data->T;
    
    
    // copy command line args to other data structures
    data->NT = flags->NT;
    for(int j=0; j<flags->NT; j++)
    {
        data->t0[j]   = data->t0[0] + j*(data->T + data->tgap[0]);
        data->tgap[j] = data->tgap[0];
    }
    
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
//    sprintf(filename,"gb_mcmc.log");
//    FILE *runlog = fopen(filename,"w");
//    print_version(runlog);
    
    //Report on set parameters
//    if(!flags->quiet) print_run_settings(argc, argv, data, orbit, flags, stdout);
//    print_run_settings(argc, argv, data, orbit, flags, runlog);
    
//    fclose(runlog);
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

void parse_vb_list(int argc, char **argv, struct Flags *flags)
{
    flags->NVB=0;
    int vb_list_flag = 0;
    
    static struct option long_options[] =
    {
        /* These options set a flag. */
        {"known-sources", required_argument, 0, 0},
        {0, 0, 0, 0}
    };
    
    opterr = 0;
    int opt=0;
    int long_index=0;
    
    //copy argv since getopt permutes order
    char **argv_copy=malloc((argc+1) * sizeof *argv_copy);
    copy_argv(argc,argv,argv_copy);

    
    //Loop through argv string and find argument for verification binaries
    while ((opt = getopt_long_only(argc, argv_copy,"apl:b:", long_options, &long_index )) != -1)
    {
        
        switch (opt)
        {
            case 0:
                if(strcmp("known-sources", long_options[long_index].name) == 0)
                {
                    strcpy(flags->vbFile,optarg);
                    vb_list_flag=1;
                }

                break;
            default:
                break;
                //print_usage();
                //exit(EXIT_FAILURE);
        }
    }
    
    if(vb_list_flag)
    {
        //count lines in the file (one source per line)
        char *line;
        char buffer[MAXSTRINGSIZE];
        
        FILE *sourceFile = fopen(flags->vbFile,"r");
        
        //strip off header
        line = fgets(buffer, MAXSTRINGSIZE, sourceFile);
        if(line==NULL)
        {
            fprintf(stderr,"Error reading %s\n",flags->vbFile);
            exit(1);
        }
        
        while( (line=fgets(buffer, MAXSTRINGSIZE, sourceFile)) != NULL) flags->NVB++;
        fclose(sourceFile);
    }

    //reset opt counter
    optind = 0;
}

void save_chain_state(struct Data *data, struct Model **model, struct Chain *chain, struct Flags *flags, int step)
{
    char filename[PATH_BUFSIZE];
    FILE *stateFile;
    for(int ic=0; ic<chain->NC; ic++)
    {
        pathprintf(filename,"%s/chain_state_%i.dat",chain->chkptDir,ic);
        stateFile = fopen(filename,"w");
        
        int n = chain->index[ic];
        
        fprintf(stateFile,"%.12g\n",chain->logLmax);
        
        print_chain_state(data, chain, model[n], flags, stateFile, step);
        print_noise_state(data, model[n], stateFile, step);
        if(flags->calibration)
            print_calibration_state(data, model[n], stateFile, step);
        
        int D = model[n]->Nlive;
        for(int i=0; i<D; i++)
        {
            print_source_params(data,model[n]->source[i],stateFile);
            fprintf(stateFile,"\n");
        }
        
        fclose(stateFile);
    }
}

void restore_chain_state(struct Orbit *orbit, struct Data *data, struct Model **model, struct Chain *chain, struct Flags *flags, int *step)
{
    char filename[PATH_BUFSIZE];
    FILE *stateFile;
    chain->logLmax=0.0;
    for(int ic=0; ic<chain->NC; ic++)
    {
        pathprintf(filename,"%s/checkpoint/chain_state_%i.dat",flags->runDir,ic);
        stateFile = fopen(filename,"r");
        
        int n = chain->index[ic];
        
        int check = fscanf(stateFile,"%lg",&chain->logLmax);
        if(!check)
        {
            fprintf(stderr,"Error reading checkpoint file\n");
            exit(1);
        }
        
        scan_chain_state(data, chain, model[n], flags, stateFile, step);
        scan_noise_state(data, model[n], stateFile, step);
        if(flags->calibration)
            scan_calibration_state(data, model[n], stateFile, step);
        
        int D = model[n]->Nlive;
        for(int i=0; i<D; i++)
        {
            scan_source_params(data,model[n]->source[i], stateFile);
            galactic_binary_fisher(orbit, data, model[n]->source[i], data->noise[0]);
        }
        
        generate_noise_model(data, model[n]);
        generate_signal_model(orbit, data, model[n], -1);
        
        if(!flags->prior)
        {
            model[n]->logL = gaussian_log_likelihood(data, model[n]);
            model[n]->logLnorm = gaussian_log_likelihood_constant_norm(data, model[n]);
        }
        else model[n]->logL = model[n]->logLnorm = 0.0;
        
        
        fclose(stateFile);
    }
}

void print_chain_files(struct Data *data, struct Model **model, struct Chain *chain, struct Flags *flags, int step)
{
    int i,n,ic;
    
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
            logL += model[n]->logL+model[n]->logLnorm;
            fprintf(chain->likelihoodFile,  "%lg ",logL);
            fprintf(chain->temperatureFile, "%lg ",1./chain->temperature[ic]);
        }
        fprintf(chain->likelihoodFile, "\n");
        fprintf(chain->temperatureFile,"\n");
    }
    
    //Print cold chains
    n = chain->index[0];
    
    print_chain_state(data, chain, model[n], flags, chain->chainFile[0], step);
    if(!flags->quiet || step>0)
        print_noise_state(data, model[n], chain->noiseFile[0], step);
    if(flags->calibration)
        print_calibration_state(data, model[n], chain->calibrationFile[0], step);
    
    if(flags->verbose)
    {
        fflush(chain->chainFile[0]);
        fflush(chain->noiseFile[0]);
        if(flags->calibration) fflush(chain->calibrationFile[0]);
    }
    
    //Print sampling parameters
    int D = model[n]->Nlive;
    for(i=0; i<D; i++)
    {
        print_source_params(data,model[n]->source[i],chain->parameterFile[0]);
        if(flags->verbose)
        {
            //numerical SNR
            double snr_n = snr(model[n]->source[i], data->noise[0]);
            //analytic SNR
            double snr_a = analytic_snr(exp(model[n]->source[i]->params[3]), data->noise[0]->SnA[0], data->sine_f_on_fstar, data->sqT);
            
            fprintf(chain->parameterFile[0],"%lg %lg ",snr_a,snr_n);
        }
        fprintf(chain->parameterFile[0],"\n");
        if(flags->verbose)fflush(chain->parameterFile[0]);
        
        if(step>0)
        {
            if(chain->dimensionFile[D]==NULL)
            {
                char filename[PATH_BUFSIZE];
                pathprintf(filename,"%s/dimension_chain.dat.%i",chain->chainDir,D);
                if(flags->resume)chain->dimensionFile[D] = fopen(filename,"a");
                else             chain->dimensionFile[D] = fopen(filename,"w");
            }
            print_source_params(data,model[n]->source[i],chain->dimensionFile[D]);
            fprintf(chain->dimensionFile[D],"\n");
        }
    }
    
    //Print calibration parameters
    
    //Print hot chains if verbose flag
    if(flags->verbose)
    {
        for(ic=1; ic<chain->NC; ic++)
        {
            n = chain->index[ic];
            print_chain_state(data, chain, model[n], flags, chain->chainFile[ic], step);
            print_noise_state(data, model[n], chain->noiseFile[ic], step);
        }//loop over chains
    }//verbose flag
}

void scan_chain_state(struct Data *data, struct Chain *chain, struct Model *model, struct Flags *flags, FILE *fptr, int *step)
{
    int check = 0;
    check += fscanf(fptr, "%i",step);
    check += fscanf(fptr, "%i",&model->Nlive);
    check += fscanf(fptr, "%lg",&model->logL);
    check += fscanf(fptr, "%lg",&model->logLnorm);
    for(int j=0; j<flags->NT; j++) check += fscanf(fptr, "%lg",&model->t0[j]);
    if(!check)
    {
        fprintf(stderr,"Error reading checkpoint files\n");
        exit(1);
    }
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
    int check = fscanf(fptr, "%i %lg %lg",step, &model->logL,&model->logLnorm);
    
    for(int i=0; i<model->NT; i++)
    {
        switch(data->Nchannel)
        {
            case 1:
                check += fscanf(fptr, "%lg", &model->calibration[i]->dampX);
                check += fscanf(fptr, "%lg", &model->calibration[i]->dphiX);
                break;
            case 2:
                check += fscanf(fptr, "%lg", &model->calibration[i]->dampA);
                check += fscanf(fptr, "%lg", &model->calibration[i]->dphiA);
                check += fscanf(fptr, "%lg", &model->calibration[i]->dampE);
                check += fscanf(fptr, "%lg", &model->calibration[i]->dphiE);
                break;
        }
    }
    
    if(!check)
    {
        fprintf(stderr,"Error reading calibration files\n");
        exit(1);
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
    int check=0;
    check+=fscanf(fptr, "%i ",step);
    check+=fscanf(fptr, "%lg %lg ", &model->logL, &model->logLnorm);
    
    for(int i=0; i<model->NT; i++)
    {
        switch(data->Nchannel)
        {
            case 1:
                check+=fscanf(fptr, "%lg", &model->noise[i]->etaX);
                break;
            case 2:
                check+=fscanf(fptr, "%lg", &model->noise[i]->etaA);
                check+=fscanf(fptr, "%lg", &model->noise[i]->etaE);
                break;
        }
    }
    if(!check)
    {
        fprintf(stderr,"Error reading noise file\n");
        exit(1);
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
    int check = 0;
    check+=fscanf(fptr,"%lg",&source->f0);
    check+=fscanf(fptr,"%lg",&source->dfdt);
    check+=fscanf(fptr,"%lg",&source->amp);
    check+=fscanf(fptr,"%lg",&source->phi);
    check+=fscanf(fptr,"%lg",&source->costheta);
    check+=fscanf(fptr,"%lg",&source->cosi);
    check+=fscanf(fptr,"%lg",&source->psi);
    check+=fscanf(fptr,"%lg",&source->phi0);
    if(source->NP>8)
        check+=fscanf(fptr,"%lg",&source->d2fdt2);
    
    if(!check)
    {
        fprintf(stdout,"Error reading source file\n");
        exit(1);
    }
    
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

void print_waveform_strain(struct Data *data, struct Model *model, FILE *fptr)
{
    for(int n=0; n<data->N; n++)
    {
        int re = 2*n;
        int im = re+1;
        double f = data->fmin + (double)n/data->T;

        int i = 0;
        fprintf(fptr,"%.12g ",f);
        fprintf(fptr,"%.12g ",model->tdi[i]->A[re]);
        fprintf(fptr,"%.12g ",model->tdi[i]->A[im]);
        fprintf(fptr,"%.12g ",model->tdi[i]->E[re]);
        fprintf(fptr,"%.12g\n",model->tdi[i]->E[im]);
    }
}


void print_waveform_draw(struct Data *data, struct Model *model, struct Flags *flags)
{
    FILE *fptr;
    char filename[PATH_BUFSIZE];
    
    pathprintf(filename,"%s/waveform_draw.dat",data->dataDir);
    fptr=fopen(filename,"w");
    print_waveform(data, model, fptr);
    fclose(fptr);
}

void print_noise_reconstruction(struct Data *data, struct Flags *flags)
{
    FILE *fptr_Snf;
    char filename[PATH_BUFSIZE];
    
    for(int k=0; k<data->NT; k++)
    {
        pathprintf(filename,"%s/power_noise_t%i.dat",data->dataDir,k);
        fptr_Snf=fopen(filename,"w");

        for(int i=0; i<data->N; i++)
        {
            gsl_sort(data->S_pow[i][0][k],1,data->Nwave);
            gsl_sort(data->S_pow[i][1][k],1,data->Nwave);

            
            double f = (double)(i+data->qmin)/data->T;
            
            double A_med   = gsl_stats_median_from_sorted_data   (data->S_pow[i][0][k], 1, data->Nwave);
            double A_lo_50 = gsl_stats_quantile_from_sorted_data (data->S_pow[i][0][k], 1, data->Nwave, 0.25);
            double A_hi_50 = gsl_stats_quantile_from_sorted_data (data->S_pow[i][0][k], 1, data->Nwave, 0.75);
            double A_lo_90 = gsl_stats_quantile_from_sorted_data (data->S_pow[i][0][k], 1, data->Nwave, 0.05);
            double A_hi_90 = gsl_stats_quantile_from_sorted_data (data->S_pow[i][0][k], 1, data->Nwave, 0.95);
            
            double E_med   = gsl_stats_median_from_sorted_data   (data->S_pow[i][1][k], 1, data->Nwave);
            double E_lo_50 = gsl_stats_quantile_from_sorted_data (data->S_pow[i][1][k], 1, data->Nwave, 0.25);
            double E_hi_50 = gsl_stats_quantile_from_sorted_data (data->S_pow[i][1][k], 1, data->Nwave, 0.75);
            double E_lo_90 = gsl_stats_quantile_from_sorted_data (data->S_pow[i][1][k], 1, data->Nwave, 0.05);
            double E_hi_90 = gsl_stats_quantile_from_sorted_data (data->S_pow[i][1][k], 1, data->Nwave, 0.95);
                        
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
        fclose(fptr_Snf);
    }
    
}

void print_waveforms_reconstruction(struct Data *data, struct Flags *flags)
{
    char filename[PATH_BUFSIZE];
    FILE *fptr_rec;
    FILE *fptr_res;
    FILE *fptr_var;
    
    //get variance of residual
    double ***res_var = malloc(data->N*sizeof(double **));
    for(int n=0; n<data->N; n++)
    {
        res_var[n] = malloc(data->Nchannel*sizeof(double *));
        for(int m=0; m<data->Nchannel; m++)
        {
            res_var[n][m] =  calloc(data->NT,sizeof(double));
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
        
        pathprintf(filename,"%s/power_reconstruction_t%i.dat",data->dataDir,k);
        fptr_rec=fopen(filename,"w");
        pathprintf(filename,"%s/power_residual_t%i.dat",data->dataDir,k);
        fptr_res=fopen(filename,"w");
        pathprintf(filename,"%s/variance_residual_t%i.dat",data->dataDir,k);
        fptr_var=fopen(filename,"w");
        
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
            
        }
        
        fclose(fptr_var);
        fclose(fptr_res);
        fclose(fptr_rec);
    }
    
    for(int n=0; n<data->N; n++)
    {
        for(int m=0; m<data->Nchannel; m++)
        {
            free(res_var[n][m]);
        }
        free(res_var[n]);
    }
    free(res_var);
}

void print_data(struct Data *data, struct TDI *tdi, struct Flags *flags, int t_index)
{
    char filename[PATH_BUFSIZE];
    FILE *fptr;

    pathprintf(filename,"%s/waveform_injection_%i.dat",data->dataDir,t_index);
    fptr=fopen(filename,"w");
    for(int i=0; i<data->N; i++)
    {
        double f = (double)(i+data->qmin)/data->T;
        fprintf(fptr,"%lg %lg %lg %lg %lg",
                f,
                tdi->A[2*i],tdi->A[2*i+1],
                tdi->E[2*i],tdi->E[2*i+1]);
        fprintf(fptr,"\n");
    }
    fclose(fptr);
    
    pathprintf(filename,"%s/power_injection_%i.dat",data->dataDir,t_index);
    fptr=fopen(filename,"w");
    for(int i=0; i<data->N; i++)
    {
        double f = (double)(i+data->qmin)/data->T;
        fprintf(fptr,"%.12g %lg %lg ",
                f,
                tdi->A[2*i]*tdi->A[2*i]+tdi->A[2*i+1]*tdi->A[2*i+1],
                tdi->E[2*i]*tdi->E[2*i]+tdi->E[2*i+1]*tdi->E[2*i+1]);
        fprintf(fptr,"\n");
    }
    fclose(fptr);
    
    pathprintf(filename,"%s/power_data_%i.dat",data->dataDir,t_index);
    fptr=fopen(filename,"w");
    
    for(int i=0; i<data->N; i++)
    {
        double f = (double)(i+data->qmin)/data->T;
        fprintf(fptr,"%.12g %lg %lg ",
                f,
                tdi->A[2*i]*tdi->A[2*i]+tdi->A[2*i+1]*tdi->A[2*i+1],
                tdi->E[2*i]*tdi->E[2*i]+tdi->E[2*i+1]*tdi->E[2*i+1]);
        fprintf(fptr,"\n");
    }
    fclose(fptr);
    
    pathprintf(filename,"%s/data_%i.dat",data->dataDir,t_index);
    fptr=fopen(filename,"w");
    
    for(int i=0; i<data->N; i++)
    {
        double f = (double)(i+data->qmin)/data->T;
        fprintf(fptr,"%.12g %lg %lg %lg %lg",
                f,
                tdi->A[2*i],tdi->A[2*i+1],
                tdi->E[2*i],tdi->E[2*i+1]);
        fprintf(fptr,"\n");
    }
    fclose(fptr);
    
    pathprintf(filename,"%s/power_noise_%i.dat",data->dataDir,t_index);
    fptr=fopen(filename,"w");
    
    for(int i=0; i<data->N; i++)
    {
        double f = (double)(i+data->qmin)/data->T;
        fprintf(fptr,"%.12g %lg %lg ",
                f,
                data->noise[t_index]->SnA[i],
                data->noise[t_index]->SnE[i]);
        fprintf(fptr,"\n");
    }
    fclose(fptr);
}

void print_evidence(struct Chain *chain,struct Flags *flags)
{
    char filename[PATH_BUFSIZE];
    pathprintf(filename,"%s/evidence.dat",flags->runDir);
    FILE *zFile = fopen(filename,"w");
    for(int i=0; i<flags->DMAX; i++) fprintf(zFile,"%i %i\n",i,chain->dimension[0][i]);
    fclose(zFile);
}


