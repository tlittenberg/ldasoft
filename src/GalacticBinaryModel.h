//
//  GalacticBinaryModel.h
//  
//
//  Created by Littenberg, Tyson B. (MSFC-ZP12) on 1/15/17.
//
//

#ifndef GalacticBinaryModel_h
#define GalacticBinaryModel_h

#include <stdio.h>

void map_params_to_array(struct Source *source, double *params);
void map_array_to_params(struct Source *source, double *params);


void initialize_chain(struct Chain *chain, long seed, int i);
void initialize_model(struct Model *model, int Nmax, int NFFT, int Nchannel);
void initialize_noise(struct Noise *noise, int NFFT);
void initialize_tdi(struct TDI *tdi, int NFFT, int Nchannel);
void initialize_source(struct Source *source, int NFFT, int Nchannel);

void free_tdi(struct TDI *tdi);
void free_noise(struct Noise *noise);
void free_model(struct Model *model);
void free_source(struct Source *source);

#endif /* GalacticBinaryModel_h */
