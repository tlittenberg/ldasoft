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
  struct Source **source; //source structure contains parameters, defined in GalacticBinary.h
};



#endif /* GalacticBinaryCatalog_h */
