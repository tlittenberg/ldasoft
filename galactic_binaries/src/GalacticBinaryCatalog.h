//
//  GalacticBinaryCatalog.h
//  
//
//  Created by Littenberg, Tyson B. (MSFC-ST12) on 9/17/19.
//

#ifndef GalacticBinaryCatalog_h
#define GalacticBinaryCatalog_h

#include <stdio.h>

struct Catalog
{
  int N; //number of discrete sources in catalog
  struct Entry **entry; //discrete catalog entry
};

struct Entry
{
  int I;                  //number of chain samples
  char name[128];         //source name
  struct Source **source; //source structure contains parameters, defined in GalacticBinary.h
  double *match;          //match between sample and ref. source
  double evidence;        //source evidence
  double SNR;             //reference SNR of source
  int i;                  //sample containing med. freq.
};

void alloc_entry(struct Entry *entry, int IMAX);
void create_new_source(struct Catalog *catalog, struct Source *sample, struct Noise *noise, int IMAX, int NFFT, int Nchannel, int NP);
void append_sample_to_entry(struct Entry *entry, struct Source *sample, int IMAX, int NFFT, int Nchannel, int NP);


#endif /* GalacticBinaryCatalog_h */
