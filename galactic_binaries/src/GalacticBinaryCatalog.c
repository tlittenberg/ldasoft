//
//  GalacticBinaryCatalog.c
//  
//
//  Created by Littenberg, Tyson B. (MSFC-ST12) on 9/17/19.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>

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

int main(int argc, char *argv[])
{

    
    FILE *chain_file1;
    int NTEMP = 1;   //needed size of data structure

    
    /* Allocate data structures */
    struct Flags *flags = malloc(sizeof(struct Flags));
    struct Orbit *orbit = malloc(sizeof(struct Orbit));
    struct Chain *chain = malloc(sizeof(struct Chain));
    struct Data  **data_tmp = malloc(sizeof(struct Data*)*NTEMP); //data[NF]
    
    
    /* Parse command line and set defaults/flags */
    for(int i=0; i<NTEMP; i++)
    {
        data_tmp[i] = malloc(sizeof(struct Data));
        data_tmp[i]->t0   = malloc( NTEMP * sizeof(double) );
    }
    parse(argc,argv,data_tmp,orbit,flags,chain,NTEMP);
    alloc_data(data_tmp, flags);
    struct Data *data  = data_tmp[0];
    data->qmin = (int)(data->fmin*data->T);
  
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
    
    chain_file1 = fopen(flags->matchInfile1,"r");

    if ( chain_file1 == NULL )
    {
        printf("match-in1 is null\n");
    }
    

    //count lines in chain file
    double junk;
    int N=0;
    while(!feof(chain_file1))
    {
        fscanf(chain_file1,"%lg %lg %lg %lg %lg %lg %lg %lg",&junk,&junk,&junk,&junk,&junk,&junk,&junk,&junk);
        N++;
    }
    N--;
    rewind(chain_file1);
    
  
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
  
    //alias for catalog->entry pointers used later on
    struct Entry *entry = NULL;

    struct Source *sample = NULL;
    sample = malloc(sizeof *sample);
    alloc_source(sample, data->N, data->Nchannel, data->NP);

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
    scan_source_params(data, sample, chain_file1);
    
    //Book-keeping of waveform in time-frequency volume
    galactic_binary_alignment(orbit, data, sample);
    
    //calculate waveform model of sample
    galactic_binary(orbit, data->format, data->T, data->t0[0], sample->params, data->NP, sample->tdi->X, sample->tdi->A, sample->tdi->E, sample->BW, data->Nchannel);
        
    //add new source to catalog
    create_new_source(catalog, sample, IMAX, data->N, sample->tdi->Nchannel, data->NP);
  }


  /* ****************************************************************/
  /*            Now loop over the rest of the chain file            */
  /* ****************************************************************/
  
  printf("Looping over chain file\n");
  for(int i=1; i<IMAX; i++)
  {
    if(i%(IMAX/100)==0)printProgress((double)i/(double)IMAX);

    //check each source in chain sample
    for(int d=0; d<DMAX; d++)
    {
        //parse source parameters
        scan_source_params(data, sample, chain_file1);

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
          append_sample_to_entry(entry, sample, IMAX, data->N, data->Nchannel, data->NP);

          //stop looping over entries in catalog
          break;
        }

      }//end loop over catalog entries

      //if the match tolerence is never met, add as new source
      if(!matchFlag)create_new_source(catalog, sample, IMAX, data->N, data->Nchannel, data->NP);

    }//end loop over sources in chain sample

  }//end loop over chain
  printf("\n");
  fclose(chain_file1);
  

  /* ****************************************************************/
  /*             Select entries that have enough weight             */
  /* ****************************************************************/
  double weight;
  double weight_threshold = 0.2; //how many samples must an entry have to any hope of counting?

  int detections = 0;           //number of entries that meet the weight threshold
  int detection_index[NMAX];    //list of entry indicies for detected sources

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
  /*           Format selected entries as L4 data products           */
  /* *************************************************************** */
  
  char outdir[MAXSTRINGSIZE];
  sprintf(outdir,"catalog");
  mkdir(outdir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  FILE *catalogFile = fopen("post/entries.dat","w");

  double *f_vec;
  double f_med;
  for(int d=0; d<detections; d++)
  {
    int n = detection_index[d];
    entry = catalog->entry[n];
    
    //get median frequency as identifier of source
    f_vec = malloc(entry->I*sizeof(double));
    for(int i=0; i<entry->I; i++) f_vec[i] = entry->source[i]->f0;
    gsl_sort(f_vec, 1, entry->I);
    f_med = gsl_stats_median_from_sorted_data(f_vec, 1, entry->I);

    //name source based on median frequency
    sprintf(entry->name,"GW%08d",(int)(f_med*1e8));
    
    //evidence for source related to number of visits in the chain
    entry->evidence = (double)entry->I/(double)IMAX;

    fprintf(catalogFile,"%s %lg\n",entry->name, entry->evidence);
  }
  
    
  /* *************************************************************** */
  /*           Save source detection parameters to file              */
  /* *************************************************************** */
  
  for(int d=0; d<detections; d++)
  {
    //open file for detection
    
    char filename[100];
    FILE *out;
    
    int n = detection_index[d];
    entry = catalog->entry[n];
    
    sprintf(filename, "catalog/%s_chain.dat", entry->name);
    out = fopen( filename, "w");
    
    
    //add parameters to file
    
    for(int k=0; k<entry->I; k++)
    {
      fprintf(out,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",entry->source[k]->f0,entry->source[k]->dfdt,entry->source[k]->amp,entry->source[k]->phi,entry->source[k]->costheta,entry->source[k]->cosi,entry->source[k]->psi,entry->source[k]->phi0);
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
}

void create_new_source(struct Catalog *catalog, struct Source *sample, int IMAX, int NFFT, int Nchannel, int NP)
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
