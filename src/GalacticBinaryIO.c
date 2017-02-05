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

#include "LISA.h"
#include "GalacticBinary.h"

void print_usage()
{
  fprintf(stdout,"\n");
  fprintf(stdout,"Usage: \n");
  fprintf(stdout,"REQUIRED:\n");
  fprintf(stdout,"       --orbit       : orbit ephemerides file             \n");
  fprintf(stdout,"OPTIONAL:\n");
  fprintf(stdout,"  -h | --help        : print help message and exit        \n");
  fprintf(stdout,"  -v | --verbose     : enable verbose output              \n");
  fprintf(stdout,"       --start-time  : initial time of segment  (0)       \n");
  fprintf(stdout,"       --duration    : duration of time segment (62914560)\n");
  fprintf(stdout,"--\n");
  fprintf(stdout,"EXAMPLE:\n");
  fprintf(stdout,"./gb_mcmc --orbit ../config/OrbitConfig1.txt --verbose --inj ../data/sources/RXJ0806.dat\n");
  fprintf(stdout,"\n");
  exit(EXIT_FAILURE);
}

void parse(int argc, char **argv, struct Data *data, struct Orbit *orbit, struct Flags *flags)
{
  //Set defaults
  flags->verbose   = 0;
  flags->injection = 0;
  flags->zeroNoise = 0;
  data->t0 = 0.0;
  data->T  = 62914560.0; /* two "mldc years" at 15s sampling */
  
  data->cseed = 150914;
  data->nseed = 151226;
  data->iseed = 151012;
  
  if(argc==1) print_usage();
  
  //Specifying the expected options
  static struct option long_options[] =
  {
    /* These options set a flag. */
    {"duration",  required_argument, 0, 0},
    {"start-time",required_argument, 0, 0},
    {"orbit",     required_argument, 0, 0},
    {"chainseed", required_argument, 0, 0},
    {"noiseseed", required_argument, 0, 0},
    {"injseed",   required_argument, 0, 0},
    {"inj",       required_argument, 0, 0},
    
    /* These options don’t set a flag.
     We distinguish them by their indices. */
    {"help",    no_argument,       0,  'h' },
    {"verbose", no_argument,       0,  'v' },
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
        if(strcmp("duration",  long_options[long_index].name) == 0) data->T     = (double)atof(optarg);
        if(strcmp("start-time",long_options[long_index].name) == 0) data->t0    = (double)atof(optarg);
        if(strcmp("chainseed", long_options[long_index].name) == 0) data->cseed = (long)atoi(optarg);
        if(strcmp("noiseseed", long_options[long_index].name) == 0) data->nseed = (long)atoi(optarg);
        if(strcmp("injseed",   long_options[long_index].name) == 0) data->iseed = (long)atoi(optarg);
        if(strcmp("zero-noise",long_options[long_index].name) == 0) flags->zeroNoise = 1;
        if(strcmp("orbit",     long_options[long_index].name) == 0) sprintf(orbit->OrbitFileName,"%s",optarg);
        if(strcmp("inj",       long_options[long_index].name) == 0)
        {
          sprintf(data->injFile,"%s",optarg);
          flags->injection=1;
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
  fprintf(stdout,"  Command Line: ");
  for(opt=0; opt<argc; opt++) fprintf(stdout,"%s ",argv[opt]);
  fprintf(stdout,"\n");
  fprintf(stdout,"\n");
  fprintf(stdout,"  Orbit file is ....... %s \n",orbit->OrbitFileName);
  fprintf(stdout,"  Data start time is .. %.0f \n",data->t0);
  fprintf(stdout,"  Data duration is .... %.0f \n",data->T);
  fprintf(stdout,"  MCMC chain seed is .. %li  \n",data->cseed);
  fprintf(stdout,"  Noise seed is ....... %li  \n",data->nseed);
  fprintf(stdout,"  Injection seed is ... %li  \n",data->iseed);
  fprintf(stdout,"\n");
  fprintf(stdout,"================= RUN FLAGS ================\n");
  if(flags->verbose)  fprintf(stdout,"  Verbose flag ........ ENABLED \n");
  else                fprintf(stdout,"  Verbose flag ........ DISABLED\n");
  if(flags->injection)fprintf(stdout,"  Injection is ........ %s\n",data->injFile);
  else                fprintf(stdout,"  Injection is ........ DISABLED\n");
  if(flags->zeroNoise)fprintf(stdout,"  Noise realization is. DISABLED\n");
  else                fprintf(stdout,"  Noise realization is. ENABLED\n");
  
  fprintf(stdout,"\n");
  fprintf(stdout,"\n");
  
}



