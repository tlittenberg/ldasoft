/*
*  Copyright (C) 2019 Tyson B. Littenberg (MSFC-ST12), Neil J. Cornish
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


#ifndef GalacticBinaryModel_h
#define GalacticBinaryModel_h

#include <stdio.h>

void simualte_data(struct Data *data, struct Flags *flags, struct Source **injections, int Ninj);
void generate_signal_model(struct Orbit *orbit, struct Data *data, struct Model *model, int index);
void generate_noise_model(struct Data *data, struct Model *model);
void generate_calibration_model(struct Data *data, struct Model *model);
void apply_calibration_model(struct Data *data, struct Model *model);

double gaussian_log_likelihood(struct Orbit *orbit, struct Data *data, struct Model *model);
double gaussian_log_likelihood_constant_norm(struct Data *data, struct Model *model);
double gaussian_log_likelihood_model_norm(struct Data *data, struct Model *model);

int update_max_log_likelihood(struct Model ***model, struct Chain *chain, struct Flags *flags);

void map_params_to_array(struct Source *source, double *params, double T);
void map_array_to_params(struct Source *source, double *params, double T);

void alloc_data(struct Data **data_vec, struct Flags *flags);

void initialize_chain(struct Chain *chain, struct Flags *flags, long *seed, const char *mode);
void alloc_model(struct Model *model, int Nmax, int NFFT, int Nchannel, int NP, int NT);
void alloc_noise(struct Noise *noise, int NFFT);
void alloc_tdi(struct TDI *tdi, int NFFT, int Nchannel);
void alloc_source(struct Source *source, int NFFT, int Nchannel, int NP);
void alloc_calibration(struct Calibration *calibration);

int compare_model(struct Model *a, struct Model *b);

void copy_source(struct Source *origin, struct Source *copy);
void copy_model(struct Model *origin, struct Model *copy);
void copy_tdi(struct TDI *origin, struct TDI *copy);
void copy_noise(struct Noise *origin, struct Noise *copy);
void copy_calibration(struct Calibration *origin, struct Calibration *copy);

void free_tdi(struct TDI *tdi);
void free_noise(struct Noise *noise);
void free_model(struct Model *model);
void free_source(struct Source *source);
void free_chain(struct Chain *chain, struct Flags *flags);
void free_calibration(struct Calibration *calibration);

#endif /* GalacticBinaryModel_h */
