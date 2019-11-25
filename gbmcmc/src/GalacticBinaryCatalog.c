/*
*  Copyright (C) 2019 Tyson B. Littenberg (MSFC-ST12), Kristen Lackeos, Neil J. Cornish
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



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <getopt.h>


#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* ===============  PROTOTYPE DECLARATIONS FOR INTERNAL FUNCTIONS ========= */

#include "LISA.h"
#include "Constants.h"
#include "GalacticBinary.h"
#include "GalacticBinaryIO.h"
#include "GalacticBinaryMath.h"
#include "GalacticBinaryModel.h"
#include "GalacticBinaryWaveform.h"
#include "GalacticBinaryCatalog.h"

/* ============================  MAIN PROGRAM  ============================ */
static void print_usage_catalog();
static void parse_catalog(int argc, char **argv, struct Data **data, struct Orbit *orbit, struct Flags *flags, int Nmax);

int main(int argc, char *argv[])
{

  /* ************************************************************** */
  /*             Allocate & Initialize Data Structures              */
  /* ************************************************************** */

    int NTEMP = 1;   //needed size of data structure

    
    /* Allocate data structures */
    struct Flags *flags = malloc(sizeof(struct Flags));
    struct Orbit *orbit = malloc(sizeof(struct Orbit));
    struct Data  **data_tmp = malloc(sizeof(struct Data*)*NTEMP); //data[NF]
    
    
    /* Parse command line and set defaults/flags */
    for(int i=0; i<NTEMP; i++)
    {
        data_tmp[i] = malloc(sizeof(struct Data));
        data_tmp[i]->t0   = malloc( NTEMP * sizeof(double) );
    }
    parse_catalog(argc,argv,data_tmp,orbit,flags,NTEMP);
    alloc_data(data_tmp, flags);
    struct Data *data  = data_tmp[0];
    data->qmin = (int)(data->fmin*data->T);
  
    //File containing chain samples
    FILE *chain_file = fopen(data->fileName,"r");

    //Orbits
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
            return(1);
            break;
    }
        

  //alias for catalog->entry pointers used later on
  struct Entry *entry = NULL;

  struct Source *sample = NULL;
  sample = malloc(sizeof *sample);
  alloc_source(sample, data->N, data->Nchannel, data->NP);

  
    //count lines in chain file
    int N=0;
    while(!feof(chain_file))
    {
        scan_source_params(data, sample, chain_file);
        N++;
    }
    N--;
    rewind(chain_file);
    
  
  //selection criteria for catalog entries
    int matchFlag;          //track if sample matches any entries
    double Match;     //match between pairs of waveforms
    double tolerance = 0.8; //tolerance on match to be considered associated
    double dqmax = 10;      //maximum frequency separation to try match calculation (in frequency bins)
  

  //Book-keeping
  int DMAX = data->DMAX; //maximum number of sources per chain sample (needs to be read in from above)
  int IMAX = N/DMAX;     //maximum number of chain samples (needs to be read in from above)
  int NMAX = N; //(absurd) upper limit on total number of catalog entries
      
    struct Catalog *catalog = NULL;
    catalog = malloc(sizeof(struct Catalog));
    catalog->N = 0; //start with 0 sources in catalog
    catalog->entry = malloc(NMAX*sizeof(struct Entry*));
  
    /* ************************************************************** */
    /*             Allocate & Initialize Instrument Model             */
    /* ************************************************************** */
    
    struct Noise *noise = NULL;
    noise = malloc(flags->NT*sizeof(struct Noise));
    alloc_noise(noise, data->N);
        
    //Noise model
    //Get noise spectrum for data segment
    for(int n=0; n<data->N; n++)
    {
        double f = data->fmin + (double)(n)/data->T;
        if(strcmp(data->format,"phase")==0)
        {
            noise->SnA[n] = AEnoise(orbit->L, orbit->fstar, f);
            noise->SnE[n] = AEnoise(orbit->L, orbit->fstar, f);
        }
        else if(strcmp(data->format,"frequency")==0)
        {
            noise->SnA[n] = AEnoise_FF(orbit->L, orbit->fstar, f);
            noise->SnE[n] = AEnoise_FF(orbit->L, orbit->fstar, f);
        }
        else
        {
            fprintf(stderr,"Unsupported data format %s",data->format);
            exit(1);
        }
    }

    
  /* **************************************************************** */
  /*        First sample of the chain initializes entry list          */
  /* **************************************************************** */
  for(int d=0; d<DMAX; d++)
  {
    
    //parse source in first sample of chain file
    scan_source_params(data, sample, chain_file);
    
    //Book-keeping of waveform in time-frequency volume
    galactic_binary_alignment(orbit, data, sample);
    
    //calculate waveform model of sample
    galactic_binary(orbit, data->format, data->T, data->t0[0], sample->params, data->NP, sample->tdi->X, sample->tdi->A, sample->tdi->E, sample->BW, data->Nchannel);
        
    //add new source to catalog
    create_new_source(catalog, sample, noise, IMAX, data->N, sample->tdi->Nchannel, data->NP);
  }


  /* ****************************************************************/
  /*            Now loop over the rest of the chain file            */
  /* ****************************************************************/
  
  fprintf(stdout,"\nLooping over chain file\n");
  for(int i=1; i<IMAX; i++)
  {
    if(i%(IMAX/100)==0)printProgress((double)i/(double)IMAX);

    //check each source in chain sample
    for(int d=0; d<DMAX; d++)
    {
        //parse source parameters
        scan_source_params(data, sample, chain_file);

        //calculate waveform model of sample
        galactic_binary_alignment(orbit, data, sample);

        //calculate waveform model of sample
        galactic_binary(orbit, data->format, data->T, data->t0[0], sample->params, data->NP, sample->tdi->X, sample->tdi->A, sample->tdi->E, sample->BW, data->Nchannel);
        
      //calculate match of sample and all entries
      matchFlag = 0;
      for(int n=0; n<catalog->N; n++)
      {
        entry = catalog->entry[n];
        
        //check frequency separation
        double q_entry  = entry->source[0]->f0 * data->T;
        double q_sample = sample->f0 * data->T;
        
        if( fabs(q_entry-q_sample) > dqmax ) Match = -1.0;
        
        //calculate match
        else Match = waveform_match(sample, entry->source[0], noise);

        if(Match > tolerance)
        {
          matchFlag = 1;
        
          //append sample to entry
          entry->match[entry->I] = Match;
          append_sample_to_entry(entry, sample, IMAX, data->N, data->Nchannel, data->NP);

          //stop looping over entries in catalog
          break;
        }

      }//end loop over catalog entries

      //if the match tolerence is never met, add as new source
      if(!matchFlag) create_new_source(catalog, sample, noise, IMAX, data->N, data->Nchannel, data->NP);

    }//end loop over sources in chain sample

  }//end loop over chain
  fprintf(stdout,"\n");
  fclose(chain_file);
  

  /* ****************************************************************/
  /*             Select entries that have enough weight             */
  /* ****************************************************************/
  double weight;
  double weight_threshold = 0.2; //what fraction of samples must an entry have to count?

  int detections = 0;            //number of entries that meet the weight threshold
  int detection_index[NMAX];     //list of entry indicies for detected sources

  for(int n=0; n<catalog->N; n++)
  {

    weight = (double)catalog->entry[n]->I/(double)IMAX;
    if(weight > weight_threshold)
    {
      //record index of catalog entry for detected source
      detection_index[detections] = n;

      //increment the number of detections
      detections++;
    }
  }
  fprintf(stdout,"\nNumber of discrete sources is %d.\n",detections);

  /* *************************************************************** */
  /*           Format selected entries as L3 data products           */
  /* *************************************************************** */
  
  char outdir[MAXSTRINGSIZE];
  sprintf(outdir,"catalog_%i",DMAX);
  mkdir(outdir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  char filename[128];
  sprintf(filename,"%s/entries.dat",outdir);
  FILE *catalogFile = fopen(filename,"w");

  double f_med;
  double *f_vec;
  fprintf(stdout,"\nPost processing events\n");
  for(int d=0; d<detections; d++)
  {
    printProgress((double)(d+1)/(double)detections);

    int n = detection_index[d];
    entry = catalog->entry[n];
    
    //get median frequency as identifier of source
    f_vec = malloc(entry->I*sizeof(double));
    for(int i=0; i<entry->I; i++) f_vec[i] = entry->source[i]->f0;
    gsl_sort(f_vec, 1, entry->I);
    f_med = gsl_stats_median_from_sorted_data(f_vec, 1, entry->I);
    free(f_vec);

    //find sample containing median frequency
    for(int i=0; i<entry->I; i++)
    {
      if(f_med == entry->source[i]->f0)
      {
        //store index of median sample
        entry->i = i;
        //replace stored SNR with median sample
        entry->SNR = snr(entry->source[i],noise);
        break;
      }
    }

    //name source based on median frequency
    sprintf(entry->name,"GW%010li",(long)(f_med*1e10));
    
    //evidence for source related to number of visits in the chain
    entry->evidence = (double)entry->I/(double)IMAX;
    
    fprintf(catalogFile,"%s %lg %lg\n",entry->name, entry->SNR, entry->evidence);
  }
  fprintf(stdout,"\n");
  
    
  /* *************************************************************** */
  /*           Save source detection parameters to file              */
  /* *************************************************************** */
  
  for(int d=0; d<detections; d++)
  {
    //open file for detection
    
    FILE *out;
    
    int n = detection_index[d];
    entry = catalog->entry[n];
    
    sprintf(filename, "%s/%s_chain.dat", outdir,entry->name);
    out = fopen( filename, "w");
    
    //add parameters to file
    for(int k=0; k<entry->I; k++)
    {
      print_source_params(data,entry->source[k],out);
      fprintf(out,"%lg\n",entry->match[k]);
    }
    
    // close detection file
    fclose(out);
  }

  if(flags->orbit)free_orbit(orbit);
  
  return 0;
}


/* ============================  SUBROUTINES  ============================ */

void alloc_entry(struct Entry *entry, int IMAX)
{
  entry->I = 0;
  entry->source = malloc(IMAX*sizeof(struct Source*));
  entry->match  = malloc(IMAX*sizeof(double));

}

void create_new_source(struct Catalog *catalog, struct Source *sample, struct Noise *noise, int IMAX, int NFFT, int Nchannel, int NP)
{
  int N = catalog->N;
  
  //allocate memory for new entry in catalog
  catalog->entry[N] = malloc(sizeof(struct Entry));
  struct Entry *entry = catalog->entry[N];
  
  alloc_entry(entry,IMAX);
  entry->source[entry->I] = malloc(sizeof(struct Source));
  alloc_source(entry->source[entry->I], NFFT, Nchannel, NP);

  //add sample to the catalog as the new entry
  copy_source(sample, entry->source[entry->I]);

  //store SNR of reference sample to set match criteria
  entry->SNR = snr(sample,noise);
  
  entry->match[entry->I] = 1.0;

  entry->I++; //increment number of samples for entry
  catalog->N++;//increment number of entries for catalog
}

void append_sample_to_entry(struct Entry *entry, struct Source *sample, int IMAX, int NFFT, int Nchannel, int NP)
{
  //malloc source structure for this sample's entry
  entry->source[entry->I] = malloc(sizeof(struct Source));
  alloc_source(entry->source[entry->I], NFFT, Nchannel, NP);
  
  //copy chain sample into catalog entry
  copy_source(sample, entry->source[entry->I]);
  
  //increment number of stored samples for this entry
  entry->I++;
}

static void print_usage_catalog()
{
  fprintf(stdout,"\n");
  fprintf(stdout,"============== GBCATALOG Usage: ============ \n");
  fprintf(stdout,"REQUIRED:\n");
  fprintf(stdout,"       --chain-file  : chain file to be sorted into catalog\n");
  fprintf(stdout,"       --sources     : maximum number of sources (10)      \n");
  fprintf(stdout,"       --fmin        : minimum frequency                   \n");
  fprintf(stdout,"\n");
  fprintf(stdout,"OPTIONAL:\n");
  fprintf(stdout,"  -h | --help        : print help message and exit         \n");
  fprintf(stdout,"       --orbit       : orbit ephemerides file (2.5 GM MLDC)\n");
  fprintf(stdout,"       --samples     : number of frequency bins (2048)     \n");
  fprintf(stdout,"       --duration    : duration of time segment (62914560) \n");
  fprintf(stdout,"       --frac-freq   : fractional frequency data (phase)   \n");
  fprintf(stdout,"       --f-double-dot: include f double dot in model       \n");
  fprintf(stdout,"       --links       : number of links [4->X,6->AE] (6)    \n");
  fprintf(stdout,"--\n");
  fprintf(stdout,"EXAMPLE:\n");
  fprintf(stdout,"./gb_catalog --fmin 0.004 --samples 256 --duration 31457280 --sources 5 --chain-file chains/dimension_chain.dat.5");
  fprintf(stdout,"\n");
  fprintf(stdout,"\n");
  exit(EXIT_FAILURE);
}

static void parse_catalog(int argc, char **argv, struct Data **data, struct Orbit *orbit, struct Flags *flags, int Nmax)
{
  print_LISA_ASCII_art(stdout);
  print_version(stdout);

  if(argc==1) print_usage_catalog();

  int DMAX_default = 10;

  //Set defaults
  flags->orbit       = 0;
  flags->match       = 0;
  flags->NT          = 1;
  flags->NDATA       = 1;
  flags->DMAX        = DMAX_default;


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

    data[i]->cseed = 150914+i*Nmax;
    data[i]->nseed = 151226+i*Nmax;
    data[i]->iseed = 151012+i*Nmax;
  }



  //Specifying the expected options
  static struct option long_options[] =
  {
    /* These options set a flag. */
    {"samples",   required_argument, 0, 0},
    {"duration",  required_argument, 0, 0},
    {"segments",  required_argument, 0, 0},
    {"sources",   required_argument, 0, 0},
    {"orbit",     required_argument, 0, 0},
    {"fmin",      required_argument, 0, 0},
    {"links",     required_argument, 0, 0},
    {"chain-file",required_argument, 0, 0},

    /* These options don’t set a flag.
     We distinguish them by their indices. */
    {"help",        no_argument, 0,'h'},
    {"verbose",     no_argument, 0,'v'},
    {"frac-freq",   no_argument, 0, 0 },
    {"f-double-dot",no_argument, 0, 0 },
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
        if(strcmp("duration",    long_options[long_index].name) == 0) data_ptr->T       = (double)atof(optarg);
        if(strcmp("fmin",        long_options[long_index].name) == 0) data_ptr->fmin    = (double)atof(optarg);
        if(strcmp("f-double-dot",long_options[long_index].name) == 0) data_ptr->NP      = 9;
        if(strcmp("sources",     long_options[long_index].name) == 0)
        {
          data_ptr->DMAX    = atoi(optarg);
          flags->DMAX       = atoi(optarg);
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
        if(strcmp("chain-file", long_options[long_index].name) == 0)
        {
            checkfile(optarg);
            sprintf(data_ptr->fileName,"%s",optarg);
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
        print_usage_catalog();
        exit(EXIT_FAILURE);
        break;
      case 'v' : flags->verbose = 1;
        break;
      default: print_usage_catalog();
        exit(EXIT_FAILURE);
    }
  }

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
  }

  //Print command line
  char filename[128];
  sprintf(filename,"catalog_%i.sh",data[0]->DMAX);
  FILE *out = fopen(filename,"w");
  fprintf(out,"#!/bin/sh\n\n");
  for(opt=0; opt<argc; opt++) fprintf(out,"%s ",argv[opt]);
  fprintf(out,"\n\n");
  fclose(out);

  //Print version control
  FILE *runlog = fopen("gb_catalog.log","w");
  print_version(runlog);

  fclose(runlog);
}
