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

#include <glass_utils.h>

#include "glass_ucb_model.h"
#include "glass_ucb_waveform.h"
#include "glass_ucb_io.h"

void print_ucb_usage()
{
    
    fprintf(stdout,"\n");
    fprintf(stdout,"================ UCB Usage: =============== \n");
    fprintf(stdout,"REQUIRED:\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"OPTIONAL:\n");

    //Priors & Proposals
    fprintf(stdout,"       ==== Priors & Proposals ==== \n");
    fprintf(stdout,"       --sources     : max # of sources for UCB model (10) \n");
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
    fprintf(stdout,"       --catalog     : list of known sources               \n");
    fprintf(stdout,"       --ucb-grid    : ucb frequency grid [filename]       \n");
    fprintf(stdout,"\n");

    //Injections
    fprintf(stdout,"       ======== Injections ======== \n");
    fprintf(stdout,"       --inj         : inject signal                       \n");
    fprintf(stdout,"       --injseed     : seed for injection parameters       \n");
    fprintf(stdout,"\n");
}

void parse_ucb_args(int argc, char **argv, struct Flags *flags)
{
    //copy argv since getopt permutes order
    char **argv_copy=malloc((argc+1) * sizeof *argv_copy);
    copy_argv(argc,argv,argv_copy);
    opterr=0; //suppress warnings about unknown arguments

    //Set defaults
    flags->NINJ        = 0;
    flags->fixSky      = 0;
    flags->fixFreq     = 0;
    flags->fixFdot     = 0;
    flags->galaxyPrior = 0;
    flags->snrPrior    = 1;
    flags->emPrior     = 0;
    flags->cheat       = 0;
    flags->detached    = 0;
    flags->knownSource = 0;
    flags->catalog     = 0;
    flags->grid        = 0;
    flags->update      = 0;
    flags->updateCov   = 0;
    flags->match       = 0;
    flags->DMAX        = 10;
    
    flags->injFile = malloc(flags->DMAX *sizeof(char *));
    for(int n=0; n<flags->DMAX ; n++) flags->injFile[n] = malloc(1024*sizeof(char));
    
    //Specifying the expected options
    static struct option long_options[] =
    {
        /* These options set a flag. */
        {"sources",    required_argument, 0, 0},
        {"injseed",    required_argument, 0, 0},
        {"inj",        required_argument, 0, 0},
        {"update-cov", required_argument, 0, 0},
        {"match-in1",  required_argument, 0, 0},
        {"match-in2",  required_argument, 0, 0},
        {"em-prior",   required_argument, 0, 0},
        {"catalog",    required_argument, 0, 0},
        {"ucb-grid",   required_argument, 0, 0},
        
        /* These options don’t set a flag.
         We distinguish them by their indices. */
        {"update",      no_argument, 0, 0 },
        {"fix-sky",     no_argument, 0, 0 },
        {"fix-freq",    no_argument, 0, 0 },
        {"fix-fdot",    no_argument, 0, 0 },
        {"galaxy-prior",no_argument, 0, 0 },
        {"snr-prior",   no_argument, 0, 0 },
        {"known-source",no_argument, 0, 0 },
        {"f-double-dot",no_argument, 0, 0 },
        {"detached",    no_argument, 0, 0 },
        {"cheat",       no_argument, 0, 0 },
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
                if(strcmp("fix-sky",     long_options[long_index].name) == 0) flags->fixSky     = 1;
                if(strcmp("fix-freq",    long_options[long_index].name) == 0) flags->fixFreq    = 1;
                if(strcmp("update",      long_options[long_index].name) == 0) flags->update     = 1;
                if(strcmp("galaxy-prior",long_options[long_index].name) == 0) flags->galaxyPrior= 1;
                if(strcmp("no-snr-prior",long_options[long_index].name) == 0) flags->snrPrior   = 0;
                //TODO: get NP back on the CLI if(strcmp("f-double-dot",long_options[long_index].name) == 0) data->NP          = 9;
                if(strcmp("detached",    long_options[long_index].name) == 0) flags->detached   = 1;
                if(strcmp("cheat",       long_options[long_index].name) == 0) flags->cheat      = 1;
                if(strcmp("sources",     long_options[long_index].name) == 0)
                {
                    flags->DMAX = atoi(optarg);
                }
                if(strcmp("em-prior",    long_options[long_index].name) == 0)
                {
                    flags->emPrior = 1;
                    sprintf(flags->pdfFile,"%s",optarg);
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
                if(strcmp("ucb-grid",long_options[long_index].name) == 0)
                {
                    checkfile(optarg);
                    flags->grid = 1;
                    sprintf(flags->ucbGridFile,"%s",optarg);
                }
                if(strcmp("inj", long_options[long_index].name) == 0)
                {
                    checkfile(optarg);
                    sprintf(flags->injFile[flags->NINJ],"%s",optarg);
                    flags->NINJ++;
                    if(flags->NINJ>flags->DMAX )
                    {
                        fprintf(stderr,"WARNING: Requested number of injections is too large (%i/%i)\n",flags->NINJ,flags->DMAX );
                        fprintf(stderr,"Should you remove at least %i --inj arguments?\n",flags->NINJ-flags->DMAX );
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
                break;
            case 'h' :
                print_ucb_usage();
                break;
            default:
                break;
        }
    }
    
    //Print version control
    //    sprintf(filename,"glass.log");
    //    FILE *runlog = fopen(filename,"w");
    //    print_version(runlog);
    
    //Report on set parameters
    //    if(!flags->quiet) print_run_settings(argc, argv_copy, data, orbit, flags, stdout);
    //    print_run_settings(argc, argv_copy, data, orbit, flags, runlog);
    
    //    fclose(runlog);

    //reset opt counter
    optind = 0;

    //free placeholder for argvs
    for(int i=0; i<=argc; i++)free(argv_copy[i]);
    free(argv_copy);
}

void parse_vgb_args(int argc, char **argv, struct Flags *flags)
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

    //free placeholder for argvs
    for(int i=0; i<=argc; i++)free(argv_copy[i]);
    free(argv_copy);
}


void print_ucb_catalog_script(struct Flags *flags, struct Data *data, struct Orbit *orbit)
{
    
    //back out original input f & N
    int samples = data->NFFT - 2*data->qpad;
    double fmin = data->fmin + data->qpad/data->T;
    
    char filename[MAXSTRINGSIZE];
    sprintf(filename,"%s/example_ucb_catalog.sh",flags->runDir);
    FILE *fptr = fopen(filename,"w");
    
    fprintf(fptr,"#!/bin/sh\n\n");
    fprintf(fptr,"if [ \"$#\" -ne 1 ]; then\n");
    fprintf(fptr,"\t echo \"You must enter model dimension\"\n");
    fprintf(fptr,"fi\n\n");
    
    fprintf(fptr,"ucb_catalog ");
    
    //Required
    if(flags->NVB>0) fprintf(fptr,"--fmin $2 ");
    else fprintf(fptr,"--fmin %.12g ", fmin);
    fprintf(fptr,"--samples %i ", samples);
    fprintf(fptr,"--padding %i ",data->qpad);
    fprintf(fptr,"--duration %f ",data->T);
    fprintf(fptr,"--start-time %f ",data->t0);
    fprintf(fptr,"--sources $1 --chain-file chains/dimension_chain.dat.$1 ");
    
    //Optional
    if(strcmp(data->format,"phase")==0)
        fprintf(fptr,"--phase ");
    if(strcmp(data->format,"sangria")==0)
        fprintf(fptr,"--sangria ");
    if(flags->orbit)
        fprintf(fptr,"--orbit %s ",orbit->OrbitFileName);
    if(UCB_MODEL_NP==9)
        fprintf(fptr,"--f-double-dot ");
    if(data->Nchannel==1)
        fprintf(fptr,"--channels 1 ");
    
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
    fprintf(fptr,"  Data start time ..... %.0f \n",data->t0);
    fprintf(fptr,"  Data start frequency. %.16g\n",data->fmin);
    fprintf(fptr,"  Data duration ....... %.0f \n",data->T);
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

void save_chain_state(struct Data *data, struct Model **model, struct Chain *chain, struct Flags *flags, int step)
{
    char filename[128];
    FILE *stateFile;
    for(int ic=0; ic<chain->NC; ic++)
    {
        sprintf(filename,"%s/chain_state_%i.dat",chain->chkptDir,ic);
        stateFile = fopen(filename,"w");
        
        int n = chain->index[ic];
        
        fprintf(stateFile,"%.12g\n",chain->logLmax);
        
        print_chain_state(data, chain, model[n], flags, stateFile, step);
        print_psd_state(data, model[n], stateFile, step);
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
    char filename[128];
    FILE *stateFile;
    chain->logLmax=0.0;
    for(int ic=0; ic<chain->NC; ic++)
    {
        sprintf(filename,"%s/checkpoint/chain_state_%i.dat",flags->runDir,ic);
        stateFile = fopen(filename,"r");
        
        int n = chain->index[ic];
        
        int check = fscanf(stateFile,"%lg",&chain->logLmax);
        if(!check)
        {
            fprintf(stderr,"Error reading checkpoint file\n");
            exit(1);
        }
        
        scan_chain_state(data, chain, model[n], flags, stateFile, step);
        scan_psd_state(data, model[n], stateFile, step);
        if(flags->calibration)
            scan_calibration_state(data, model[n], stateFile, step);
        
        int D = model[n]->Nlive;
        for(int i=0; i<D; i++)
        {
            scan_source_params(data,model[n]->source[i], stateFile);
            ucb_fisher(orbit, data, model[n]->source[i], data->noise);
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
        print_psd_state(data, model[n], chain->noiseFile[0], step);
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
            double snr_n;
            if(!strcmp("fourier",data->basis)) snr_n = snr(model[n]->source[i], data->noise);
            if(!strcmp("wavelet",data->basis)) snr_n = snr_wavelet(model[n]->source[i], data->noise);
            //analytic SNR
            double snr_a = analytic_snr(exp(model[n]->source[i]->params[3]), data->noise->C[0][0][0], data->sine_f_on_fstar, data->sqT);
            
            fprintf(chain->parameterFile[0],"%lg %lg ",snr_a,snr_n);
        }
        fprintf(chain->parameterFile[0],"\n");
        if(flags->verbose)fflush(chain->parameterFile[0]);
        
        if(step>0)
        {
            if(chain->dimensionFile[D]==NULL)
            {
                char filename[MAXSTRINGSIZE];
                sprintf(filename,"%s/dimension_chain.dat.%i",chain->chainDir,D);
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
            print_psd_state(data, model[n], chain->noiseFile[ic], step);
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
    check += fscanf(fptr, "%lg",&model->t0);
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
    fprintf(fptr, "%.12g ",model->t0);
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
    
    switch(data->Nchannel)
    {
        case 1:
            check += fscanf(fptr, "%lg", &model->calibration->dampX);
            check += fscanf(fptr, "%lg", &model->calibration->dphiX);
            break;
        case 2:
            check += fscanf(fptr, "%lg", &model->calibration->dampA);
            check += fscanf(fptr, "%lg", &model->calibration->dphiA);
            check += fscanf(fptr, "%lg", &model->calibration->dampE);
            check += fscanf(fptr, "%lg", &model->calibration->dphiE);
            break;
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
    
    switch(data->Nchannel)
    {
        case 1:
            fprintf(fptr, "%lg ", model->calibration->dampX);
            fprintf(fptr, "%lg ", model->calibration->dphiX);
            break;
        case 2:
            fprintf(fptr, "%lg ", model->calibration->dampA);
            fprintf(fptr, "%lg ", model->calibration->dphiA);
            fprintf(fptr, "%lg ", model->calibration->dampE);
            fprintf(fptr, "%lg ", model->calibration->dphiE);
            break;
    }

    fprintf(fptr, "\n");
}

void scan_psd_state(struct Data *data, struct Model *model, FILE *fptr, int *step)
{
    int check=0;
    check+=fscanf(fptr, "%i ",step);
    check+=fscanf(fptr, "%lg %lg ", &model->logL, &model->logLnorm);
    
    for(int n=0; n<data->Nchannel; n++)
        check+=fscanf(fptr, "%lg", &model->noise->eta[n]);

    if(!check)
    {
        fprintf(stderr,"Error reading noise file\n");
        exit(1);
    }
}

void print_psd_state(struct Data *data, struct Model *model, FILE *fptr, int step)
{
    fprintf(fptr, "%i ",step);
    fprintf(fptr, "%lg %lg ",model->logL, model->logLnorm);
    
    for(int n=0; n<data->Nchannel; n++)
        fprintf(fptr, "%lg ", model->noise->eta[n]);

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
    if(UCB_MODEL_NP>8)
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
    if(UCB_MODEL_NP>8)
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
    double A_re,A_im,E_re,E_im,X_re,X_im,Y_re,Y_im,Z_re,Z_im,R_re,R_im;
    
    if(!strcmp(data->format,"fourier"))
    {
        switch(data->Nchannel)
        {
            case 1:
                for(int n=0; n<data->NFFT; n++)
                {
                    n_re = 2*n;
                    n_im = n_re++;
                    
                    X_re = model->tdi->X[n_re];
                    X_im = model->tdi->X[n_im];
                    
                    data->h_rec[n_re][0][mcmc] = X_re;
                    data->h_rec[n_im][0][mcmc] = X_im;
                    
                    R_re = data->tdi->X[n_re] - X_re;
                    R_im = data->tdi->X[n_im] - X_im;
                    
                    data->h_res[n_re][0][mcmc] = R_re;
                    data->h_res[n_im][0][mcmc] = R_im;
                    
                    data->r_pow[n][0][mcmc] = R_re*R_re + R_im*R_im;
                    data->h_pow[n][0][mcmc] = X_re*X_re + X_im*X_im;
                    
                    data->S_pow[n][0][mcmc] = 1./model->noise->invC[0][0][n];
                }
                break;
            case 2:
                for(int n=0; n<data->NFFT; n++)
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
                    
                    data->h_res[n_re][0][mcmc] = R_re;
                    data->h_res[n_im][0][mcmc] = R_im;
                    
                    data->r_pow[n][0][mcmc] = R_re*R_re + R_im*R_im;
                    
                    R_re = data->tdi->E[n_re] - E_re;
                    R_im = data->tdi->E[n_im] - E_im;
                    
                    data->h_res[n_re][1][mcmc] = R_re;
                    data->h_res[n_im][1][mcmc] = R_im;
                    
                    data->r_pow[n][1][mcmc] = R_re*R_re + R_im*R_im;
                    
                    data->h_pow[n][0][mcmc] = A_re*A_re + A_im*A_im;
                    data->h_pow[n][1][mcmc] = E_re*E_re + E_im*E_im;
                    
                    data->S_pow[n][0][mcmc] = 1./model->noise->invC[0][0][n];
                    data->S_pow[n][1][mcmc] = 1./model->noise->invC[1][1][n];
                }
                break;
            case 3:
                for(int n=0; n<data->NFFT; n++)
                {
                    n_re = 2*n;
                    n_im = n_re++;
                    
                    X_re = model->tdi->X[n_re];
                    X_im = model->tdi->X[n_im];
                    Y_re = model->tdi->Y[n_re];
                    Y_im = model->tdi->Y[n_im];
                    Z_re = model->tdi->Z[n_re];
                    Z_im = model->tdi->Z[n_im];
                    
                    data->h_rec[n_re][0][mcmc] = X_re;
                    data->h_rec[n_im][0][mcmc] = X_im;
                    data->h_rec[n_re][1][mcmc] = Y_re;
                    data->h_rec[n_im][1][mcmc] = Y_im;
                    data->h_rec[n_re][2][mcmc] = Z_re;
                    data->h_rec[n_im][2][mcmc] = Z_im;
                    
                    R_re = data->tdi->X[n_re] - X_re;
                    R_im = data->tdi->X[n_im] - X_im;
                    
                    data->h_res[n_re][0][mcmc] = R_re;
                    data->h_res[n_im][0][mcmc] = R_im;
                    
                    data->r_pow[n][0][mcmc] = R_re*R_re + R_im*R_im;
                    
                    R_re = data->tdi->Y[n_re] - Y_re;
                    R_im = data->tdi->Y[n_im] - Y_im;
                    
                    data->h_res[n_re][1][mcmc] = R_re;
                    data->h_res[n_im][1][mcmc] = R_im;
                    
                    data->r_pow[n][1][mcmc] = R_re*R_re + R_im*R_im;
                    
                    R_re = data->tdi->Z[n_re] - Z_re;
                    R_im = data->tdi->Z[n_im] - Z_im;
                    
                    data->h_res[n_re][2][mcmc] = R_re;
                    data->h_res[n_im][2][mcmc] = R_im;
                    
                    data->r_pow[n][2][mcmc] = R_re*R_re + R_im*R_im;
                    
                    
                    data->h_pow[n][0][mcmc] = X_re*X_re + X_im*X_im;
                    data->h_pow[n][1][mcmc] = Y_re*Y_re + Y_im*Y_im;
                    data->h_pow[n][2][mcmc] = Z_re*Z_re + Z_im*Z_im;
                    
                    data->S_pow[n][0][mcmc] = 1./model->noise->invC[0][0][n];
                    data->S_pow[n][1][mcmc] = 1./model->noise->invC[1][1][n];
                    data->S_pow[n][2][mcmc] = 1./model->noise->invC[2][2][n];
                    
                }
                break;
        }
    }//end if Fourier
    if(!strcmp(data->basis,"wavelet"))
    {
        for(int n=0; n<data->N; n++)
        {
            data->h_rec[n][0][mcmc] = model->tdi->X[n];
            data->h_rec[n][1][mcmc] = model->tdi->Y[n];
            data->h_rec[n][2][mcmc] = model->tdi->Z[n];
            
            data->h_res[n][0][mcmc] = data->tdi->X[n] - model->tdi->X[n];
            data->h_res[n][1][mcmc] = data->tdi->Y[n] - model->tdi->Y[n];
            data->h_res[n][2][mcmc] = data->tdi->Z[n] - model->tdi->Z[n];
            
            data->r_pow[n][0][mcmc] = data->h_res[n][0][mcmc]*data->h_res[n][0][mcmc];
            data->r_pow[n][1][mcmc] = data->h_res[n][1][mcmc]*data->h_res[n][1][mcmc];
            data->r_pow[n][2][mcmc] = data->h_res[n][2][mcmc]*data->h_res[n][2][mcmc];
            
            
            data->h_rec[n][0][mcmc] = model->tdi->X[n]*model->tdi->X[n];
            data->h_rec[n][1][mcmc] = model->tdi->Y[n]*model->tdi->Y[n];
            data->h_rec[n][2][mcmc] = model->tdi->Z[n]*model->tdi->Z[n];
            
            data->S_pow[n][0][mcmc] = 1./model->noise->invC[0][0][n];
            data->S_pow[n][1][mcmc] = 1./model->noise->invC[1][1][n];
            data->S_pow[n][2][mcmc] = 1./model->noise->invC[2][2][n];
        }
    }
}

void print_waveform(struct Data *data, struct Model *model, FILE *fptr)
{
    for(int n=0; n<data->NFFT; n++)
    {
        int re = 2*n;
        int im = re+1;
        double f = data->fmin + (double)n/data->T;

        fprintf(fptr,"%.12g ",f);
        switch(data->Nchannel)
        {
            case 2:
                fprintf(fptr,"%.12g ",data->tdi->A[re]*data->tdi->A[re] + data->tdi->A[im]*data->tdi->A[im]);
                fprintf(fptr,"%.12g ",data->tdi->E[re]*data->tdi->E[re] + data->tdi->E[im]*data->tdi->E[im]);
                
                fprintf(fptr,"%.12g ",model->tdi->A[re]*model->tdi->A[re] + model->tdi->A[im]*model->tdi->A[im]);
                fprintf(fptr,"%.12g ",model->tdi->E[re]*model->tdi->E[re] + model->tdi->E[im]*model->tdi->E[im]);
                
                fprintf(fptr,"%.12g ",(data->tdi->A[re]-model->tdi->A[re])*(data->tdi->A[re]-model->tdi->A[re]) + (data->tdi->A[im]-model->tdi->A[im])*(data->tdi->A[im]-model->tdi->A[im]) );
                fprintf(fptr,"%.12g ",(data->tdi->E[re]-model->tdi->E[re])*(data->tdi->E[re]-model->tdi->E[re]) + (data->tdi->E[im]-model->tdi->E[im])*(data->tdi->E[im]-model->tdi->E[im]) );
                
                break;
            case 3:
                fprintf(fptr,"%.12g ",data->tdi->X[re]*data->tdi->X[re] + data->tdi->X[im]*data->tdi->X[im]);
                fprintf(fptr,"%.12g ",data->tdi->Y[re]*data->tdi->Y[re] + data->tdi->Y[im]*data->tdi->Y[im]);
                fprintf(fptr,"%.12g ",data->tdi->Z[re]*data->tdi->Z[re] + data->tdi->Z[im]*data->tdi->Z[im]);

                fprintf(fptr,"%.12g ",model->tdi->X[re]*model->tdi->X[re] + model->tdi->X[im]*model->tdi->X[im]);
                fprintf(fptr,"%.12g ",model->tdi->Y[re]*model->tdi->Y[re] + model->tdi->Y[im]*model->tdi->Y[im]);
                fprintf(fptr,"%.12g ",model->tdi->Z[re]*model->tdi->Z[re] + model->tdi->Z[im]*model->tdi->Z[im]);
                
                fprintf(fptr,"%.12g ",(data->tdi->X[re]-model->tdi->X[re])*(data->tdi->X[re]-model->tdi->X[re]) + (data->tdi->X[im]-model->tdi->X[im])*(data->tdi->X[im]-model->tdi->X[im]) );
                fprintf(fptr,"%.12g ",(data->tdi->Y[re]-model->tdi->Y[re])*(data->tdi->Y[re]-model->tdi->Y[re]) + (data->tdi->Y[im]-model->tdi->Y[im])*(data->tdi->Y[im]-model->tdi->Y[im]) );
                fprintf(fptr,"%.12g ",(data->tdi->Z[re]-model->tdi->Z[re])*(data->tdi->Z[re]-model->tdi->Z[re]) + (data->tdi->Z[im]-model->tdi->Z[im])*(data->tdi->Z[im]-model->tdi->Z[im]) );

                break;
        }
        fprintf(fptr,"\n");
        //    }
    }
}

void print_waveform_strain(struct Data *data, struct Model *model, FILE *fptr)
{
    if(!strcmp(data->basis,"fourier"))
    {
        for(int n=0; n<data->NFFT; n++)
        {
            int re = 2*n;
            int im = re+1;
            double f = data->fmin + (double)n/data->T;
            
            fprintf(fptr,"%.12g ",f);
            switch(data->Nchannel)
            {
                case 2:
                    fprintf(fptr,"%.12g ",model->tdi->A[re]);
                    fprintf(fptr,"%.12g ",model->tdi->A[im]);
                    fprintf(fptr,"%.12g ",model->tdi->E[re]);
                    fprintf(fptr,"%.12g\n",model->tdi->E[im]);
                    break;
                case 3:
                    fprintf(fptr,"%.12g ",model->tdi->X[re]);
                    fprintf(fptr,"%.12g ",model->tdi->X[im]);
                    fprintf(fptr,"%.12g ",model->tdi->Y[re]);
                    fprintf(fptr,"%.12g ",model->tdi->Y[im]);
                    fprintf(fptr,"%.12g ",model->tdi->Z[re]);
                    fprintf(fptr,"%.12g\n",model->tdi->Z[im]);
                    break;
                    
            }
        }
    }
    if(!strcmp(data->basis,"wavelet"))
    {
        for(int j=0; j<data->wdm->NF; j++)
        {
            for(int i=0; i<data->wdm->NT; i++)
            {
                int k;
                wavelet_pixel_to_index(data->wdm, i, j, &k);
                fprintf(fptr,"%.12g %.12g ",i*WAVELET_DURATION,j*WAVELET_BANDWIDTH + WAVELET_BANDWIDTH/2);
                fprintf(fptr,"%.12g ",model->tdi->X[k]);
                fprintf(fptr,"%.12g ",model->tdi->Y[k]);
                fprintf(fptr,"%.12g ",model->tdi->Z[k]);
                fprintf(fptr,"\n");
            }
            fprintf(fptr,"\n");
        }
    }
}


void print_waveform_draw(struct Data *data, struct Model *model, struct Flags *flags)
{
    FILE *fptr;
    char filename[128];
    
    sprintf(filename,"%s/waveform_draw.dat",data->dataDir);
    fptr=fopen(filename,"w");
    print_waveform(data, model, fptr);
    fclose(fptr);
}
void print_waveforms_reconstruction(struct Data *data, struct Flags *flags)
{
    char filename[1024];
    FILE *fptr_rec;
    FILE *fptr_res;
    FILE *fptr_var;
    
    int Nsamples;
    if(!strcmp(data->basis,"fourier")) Nsamples = data->NFFT;
    if(!strcmp(data->basis,"wavelet")) Nsamples = data->N;
    
    //get variance of residual
    double **res_var = malloc(Nsamples*sizeof(double *));
    for(int n=0; n<Nsamples; n++)
        res_var[n] = calloc(data->Nchannel,sizeof(double));
    
    for(int n=0; n<Nsamples; n++)
    {
        for(int m=0; m<data->Nchannel; m++)
        {
            gsl_sort(data->h_rec[n][m],1,data->Nwave);
        }
    }
    
    for(int n=0; n<Nsamples; n++)
    {
        for(int m=0; m<data->Nchannel; m++)
        {
            gsl_sort(data->r_pow[n][m],1,data->Nwave);
            gsl_sort(data->h_pow[n][m],1,data->Nwave);
            gsl_sort(data->S_pow[n][m],1,data->Nwave);
            res_var[n][m] = gsl_stats_variance(data->h_rec[2*n][m], 1, data->Nwave)+gsl_stats_variance(data->h_rec[2*n+1][m], 1, data->Nwave);
        }
    }
    
    sprintf(filename,"%s/power_reconstruction.dat",data->dataDir);
    fptr_rec=fopen(filename,"w");
    sprintf(filename,"%s/power_residual.dat",data->dataDir);
    fptr_res=fopen(filename,"w");
    sprintf(filename,"%s/variance_residual.dat",data->dataDir);
    fptr_var=fopen(filename,"w");
    
    double med,lo_50,hi_50,lo_90,hi_90;
    
    for(int i=0; i<Nsamples; i++)
    {
        double f = (double)(i+data->qmin)/data->T;
        fprintf(fptr_var,"%.12g ",f);
        fprintf(fptr_res,"%.12g ",f);
        fprintf(fptr_rec,"%.12g ",f);

        for(int n=0; n<data->Nchannel; n++)
        {
            fprintf(fptr_var,"%.12g ",res_var[i][n]);
            
            med   = gsl_stats_median_from_sorted_data   (data->r_pow[i][n], 1, data->Nwave);
            lo_50 = gsl_stats_quantile_from_sorted_data (data->r_pow[i][n], 1, data->Nwave, 0.25);
            hi_50 = gsl_stats_quantile_from_sorted_data (data->r_pow[i][n], 1, data->Nwave, 0.75);
            lo_90 = gsl_stats_quantile_from_sorted_data (data->r_pow[i][n], 1, data->Nwave, 0.05);
            hi_90 = gsl_stats_quantile_from_sorted_data (data->r_pow[i][n], 1, data->Nwave, 0.95);
            
            fprintf(fptr_res,"%lg ",med);
            fprintf(fptr_res,"%lg ",lo_50);
            fprintf(fptr_res,"%lg ",hi_50);
            fprintf(fptr_res,"%lg ",lo_90);
            fprintf(fptr_res,"%lg ",hi_90);
            
            med   = gsl_stats_median_from_sorted_data   (data->h_pow[i][n], 1, data->Nwave);
            lo_50 = gsl_stats_quantile_from_sorted_data (data->h_pow[i][n], 1, data->Nwave, 0.25);
            hi_50 = gsl_stats_quantile_from_sorted_data (data->h_pow[i][n], 1, data->Nwave, 0.75);
            lo_90 = gsl_stats_quantile_from_sorted_data (data->h_pow[i][n], 1, data->Nwave, 0.05);
            hi_90 = gsl_stats_quantile_from_sorted_data (data->h_pow[i][n], 1, data->Nwave, 0.95);
            
            fprintf(fptr_rec,"%lg ",med);
            fprintf(fptr_rec,"%lg ",lo_50);
            fprintf(fptr_rec,"%lg ",hi_50);
            fprintf(fptr_rec,"%lg ",lo_90);
            fprintf(fptr_rec,"%lg ",hi_90);
        }
        fprintf(fptr_var,"\n");
        fprintf(fptr_res,"\n");
        fprintf(fptr_rec,"\n");

    }
    
    fclose(fptr_var);
    fclose(fptr_res);
    fclose(fptr_rec);
    
    for(int n=0; n<Nsamples; n++)
    {
        free(res_var[n]);
    }
    free(res_var);
}

void print_evidence(struct Chain *chain,struct Flags *flags)
{
    char filename[MAXSTRINGSIZE];
    sprintf(filename,"%s/evidence.dat",flags->runDir);
    FILE *zFile = fopen(filename,"w");
    for(int i=0; i<flags->DMAX; i++) fprintf(zFile,"%i %i\n",i,chain->dimension[0][i]);
    fclose(zFile);
}


