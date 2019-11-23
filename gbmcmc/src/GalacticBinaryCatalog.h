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
