//
//  GalacticBinaryData.h
//  
//
//  Created by Littenberg, Tyson B. (MSFC-ZP12) on 2/3/17.
//
//

#ifndef GalacticBinaryData_h
#define GalacticBinaryData_h

#include <stdio.h>

void GalacticBinaryReadData(struct Data **data_vec, struct Orbit *orbit, struct Flags *flags);
void GalacticBinarySimulateData(struct Data *data);
void GalacticBinaryInjectVerificationSource(struct Data **data_vec, struct Orbit *orbit, struct Flags *flags);
void GalacticBinaryInjectSimulatedSource(struct Data **data_vec, struct Orbit *orbit, struct Flags *flags);
void GalacticBinaryCatalogSNR(struct Data *data, struct Orbit *orbit, struct Flags *flags);

#endif /* GalacticBinaryData_h */
