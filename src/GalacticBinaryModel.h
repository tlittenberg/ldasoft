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

void simualte_data(struct Data *data, struct Flags *flags, struct Source **injections, int Ninj);
void generate_signal_model(struct Orbit *orbit, struct Data *data, struct Model *model);
void generate_noise_model(struct Data *data, struct Model *model);

double gaussian_log_likelihood(struct Orbit *orbit, struct Data *data, struct Model *model);
double gaussian_log_likelihood_constant_norm(struct Data *data, struct Model *model);
double gaussian_log_likelihood_model_norm(struct Data *data, struct Model *model);

void update_max_log_likelihood(struct Model ***model, struct Chain *chain, struct Flags *flags);

void map_params_to_array(struct Source *source, double *params, double T);
void map_array_to_params(struct Source *source, double *params, double T);

void alloc_data(struct Data **data_vec, int NMCMC, int NMAX);

void initialize_chain(struct Chain *chain, struct Flags *flags, long *seed, int NC);
void alloc_model(struct Model *model, int Nmax, int NFFT, int Nchannel);
void alloc_noise(struct Noise *noise, int NFFT);
void alloc_tdi(struct TDI *tdi, int NFFT, int Nchannel);
void alloc_source(struct Source *source, int NFFT, int Nchannel);

void copy_source(struct Source *origin, struct Source *copy);
void copy_model(struct Model *origin, struct Model *copy);
void copy_tdi(struct TDI *origin, struct TDI *copy);
void copy_noise(struct Noise *origin, struct Noise *copy);

void free_tdi(struct TDI *tdi);
void free_noise(struct Noise *noise);
void free_model(struct Model *model);
void free_source(struct Source *source);
void free_chain(struct Chain *chain, struct Flags *flags);


#endif /* GalacticBinaryModel_h */
