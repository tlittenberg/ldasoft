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

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* ===============  PROTOTYPE DECLARATIONS FOR INTERNAL FUNCTIONS ========= */

#include "LISA.h"
#include "Constants.h"
#include "GalacticBinary.h"
#include "GalacticBinaryMath.h"
#include "GalacticBinaryModel.h"
#include "GalacticBinaryWaveform.h"
#include "GalacticBinaryCatalog.h"

/* ============================  MAIN PROGRAM  ============================ */

int main(int argc, char *argv[])
{

  //Load posterior samples
  /*
   -get chain file names
   -open chain files
   -figure out size of things
   */
  
  //selection criteria for catalog entries
  int matchFlag;          //track if sample matches any entries
  double Match = 0.0;     //match between pairs of waveforms
  double tolerance = 0.9; //tolerance on match to be considered associated
  double dqmax = 10;      //maximum frequency separation to try match calculation (in frequency bins)
  
  //data volume (needs to be read in somehow)
  int NFFT     = 2048;  //maximum number of samples for waveform model
  int Nchannel = 2;     //number of data channels (2 == A,e)
  int NP       = 8;     //number of parameters

  //Book-keeping
  int DMAX = 9;         //maximum number of sources per chain sample (needs to be read in from above)
  int IMAX = 10000;     //maximum number of chain samples (needs to be read in from above)
  int NMAX = DMAX*IMAX; //(absurd) upper limit on total number of catalog entries
  
  //Allocate memory and initialize structures
  struct Catalog *catalog = NULL;
  catalog = malloc(sizeof(struct Catalog));
  
  catalog->N = 0; //start with 0 sources in catalog
  catalog->entry = malloc(NMAX*sizeof(struct Entry*));
  
  //for holding current state of chain
  struct Source *sample = malloc(sizeof(struct Source));
  alloc_source(sample, NFFT, Nchannel, NP);

  //shortcut to entry pointer in catalog structure
  struct Entry *entry = NULL;

  
  /******************************************************************/
  /*             Allocate & Initialize Instrument Model             */
  /******************************************************************/

  //Orbits
  
  //Noise model
  
  //What else?
  
  
  /******************************************************************/
  /*        First sample of the chain initializes entry list        */
  /******************************************************************/
  for(int d=0; d<DMAX; d++)
  {
    //parse source in first sample of chain file
    
    //populate source structure
    
    //calculate waveform model of sample
    
    //add new source to catalog
    create_new_source(catalog, sample, IMAX, NFFT, Nchannel, NP);
  }
  

  /******************************************************************/
  /*            Now loop over the rest of the chain file            */
  /******************************************************************/
  for(int i=1; i<IMAX; i++)
  {
    //check each source in chain sample
    for(int d=0; d<DMAX; d++)
    {
      //parse source parameters
      
      //populate sample structure with current parameters
      
      //calculate waveform model of sample
      
      //calculate match of sample and all entries
      matchFlag = 0;
      for(int n=0; n<catalog->N; n++)
      {
        entry = catalog->entry[n];
       
        //check frequency separation
        double q_entry  = entry->source[0]->params[0];
        double q_sample = sample->params[0];
        if( fabs(q_entry-q_sample) < dqmax ) continue;
        
        //calculate match
        //Match = galactic_binary_match(sample, entry, noise);
      
        if(Match > tolerance)
        {
          matchFlag = 1;
          
          //append sample to entry
          append_sample_to_entry(entry, sample, IMAX, NFFT, Nchannel, NP);
          
          //stop looping over entries in catalog
          break;
        }
        
      }//end loop over catalog entries
      
      //if the match tolerence is never met, add as new source
      if(!matchFlag)create_new_source(catalog, sample, IMAX, NFFT, Nchannel, NP);
    
    }//end loop over sources in chain sample
  
  }//end loop over chain
  
  
  /******************************************************************/
  /*             Select entries that have enough weight             */
  /******************************************************************/
  double weight;
  double weight_threshold = .2; //how many samples must an entry have to any hope of counting?
  
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
  
  
  /******************************************************************/
  /*           Format selected entries as L4 data products          */
  /******************************************************************/
  for(int d=0; d<detections; d++)
  {
    int n = detection_index[d];
    entry = catalog->entry[n];
  }
  
  
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
  struct Entry *entry = catalog->entry[N];
  entry = malloc(sizeof(struct Entry));
  alloc_entry(entry,IMAX);
  entry->source[entry->I] = malloc(sizeof(struct Source));
  alloc_source(entry->source[entry->I], NFFT, Nchannel, NP);
  
  //add sample to the catalog as the new entry
  copy_source(sample, catalog->entry[N]->source[0]);
  catalog->entry[N]->I++; //increment number of samples for entry
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
