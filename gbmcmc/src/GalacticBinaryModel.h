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

/**
 @file GalacticBinaryModel.h
 \brief Functions for defining, manipulating, and evaluating Galactic Binary model
 */


#ifndef GalacticBinaryModel_h
#define GalacticBinaryModel_h

#include <stdio.h>

/**
\brief Create galactic binary model waveform

 Computes LISA response to signal with parameters indexed by `source_id`
 in the Model::Source, and creates meta template of all sources in the model.
 */
void generate_signal_model(struct Orbit *orbit, struct Data *data, struct Model *model, int source_id);

/**
\brief Modify galactic binary model waveform

 Computes LISA response to signal with parameters indexed by `source_id`
 in the Model::Source, and updates meta template of all sources in the model.
 */
void update_signal_model(struct Orbit *orbit, struct Data *data, struct Model *model_x, struct Model *model_y, int source_id);

/**
\brief F-statistic maximization of galactic binary parameters

 Wrapper of GalacticBinaryFstatistic.c functions for maximizing waveform
 over \f$(\mathcal{A},\cos\iota,\psi,\varphi_0)\f$ by filtering on
 original data.
 
 \todo test filtering on residuals
 */
void maximize_signal_model(struct Orbit *orbit, struct Data *data, struct Model *model, int source_id);

/**
 \brief Create LISA instrument noise model
 
 Computes \f$S_n(f)\f$ from Model::noise.
 */
void generate_noise_model(struct Data *data, struct Model *model);

/**
 \brief Analytic in-place inversion of noise covariance matrix
 
 Substitutes contents of noise->Cij[n] with inverse, and sets deteriminent noise->detC[n]
 */
void invert_noise_covariance_matrix(struct Noise *noise);

/**
 \brief Create LISA instrument calibration model
 
 Computes phase and amplitude corrections from Model::calibration parameters.
 */
void generate_calibration_model(struct Data *data, struct Model *model);

/**
\brief Apply amplitude and phase corrections.

Computes new LISA instrument response Model::tdi after applying
 amplitude and phase corrections from calibration parameters.
 */
void apply_calibration_model(struct Data *data, struct Model *model);

/**
 \brief Compute argument of Gaussian likelihood
 
 Computes residual of data and meta-template from Model and noise weighted inner product.
 @return \f$ -\frac{1}{2}(d-h|d-h) \f$
 */
double gaussian_log_likelihood(struct Data *data, struct Model *model);

/**
 \brief Compute normalization of Gaussian likelihood for constant noise level
 
 For noise models that are (approximately) a constant over the band
 of \f$N\f$ Fourier bins, parameterized by a multiplyer \f$\eta\f$
 
 @return \f$ \sum_{\rm TDI} N\log\eta_{\rm TDI} \f$
 */
double gaussian_log_likelihood_constant_norm(struct Data *data, struct Model *model);

/**
 \brief Compute normalization of Gaussian likelihood for arbitrary noise level
 
 For noise models that are free to vary over the band
 of \f$N\f$ Fourier bins
 
 @return \f$ \sum_{\rm TDI} \sum_f \log S_{n,{\rm TDI}}(f) \f$
 */
double gaussian_log_likelihood_model_norm(struct Data *data, struct Model *model);

/**
 \brief Compute difference in log Likelihood from changing parameters of one source
 
 Updates residuals from changing single source `source_id` and computes change in likelihood sum only over frequency range where waveforms were non-zero.
 @return \f$ -\frac{1}{2}\left[(d-h_{\rm new}|d-h__{\rm new}) - (d-h_{\rm old}|d-h_{\rm old})\right] \f$
 */
double delta_log_likelihood(struct Data *data, struct Model *model_x, struct Model *model_y, int source_id);


/**
 \brief Check for increase in maximum log likelihood
 */
int update_max_log_likelihood(struct Model **model, struct Chain *chain, struct Flags *flags);

/**
 \brief Converts physical UCB parameters to array expected by GalacticBinaryWaveform.c
 */
void map_params_to_array(struct Source *source, double *params, double T);

/**
 \brief Converts array expected by GalacticBinaryWaveform.c to
 physical UCB parameters
 */
void map_array_to_params(struct Source *source, double *params, double T);

/**
 \brief Loads spacecraft orbits, either analytically or from tabulated ephemerides.
 */
void initialize_orbit(struct Data *data, struct Orbit *orbit, struct Flags *flags);

/**
 \brief Allocates and initializes Chain structure and prepares output files.
 */
void initialize_chain(struct Chain *chain, struct Flags *flags, long *seed, const char *mode);

/** @name Allocate memory for structures */
///@{
void alloc_data(struct Data *data, struct Flags *flags);
void alloc_model(struct Model *model, int Nmax, int NFFT, int Nchannel, int NP, int NT);
void alloc_noise(struct Noise *noise, int NFFT, int Nchannel);
void alloc_source(struct Source *source, int NFFT, int Nchannel, int NP);
void alloc_calibration(struct Calibration *calibration);
///@}

/**
 \brief Shallow copy of Data structure
 */
void copy_data(struct Data *origin, struct Data *copy);

/** @name Deep copy structure contents */
///@{
void copy_source(struct Source *origin, struct Source *copy);
void copy_model(struct Model *origin, struct Model *copy);
void copy_noise(struct Noise *origin, struct Noise *copy);
void copy_model_lite(struct Model *origin, struct Model *copy);
void copy_calibration(struct Calibration *origin, struct Calibration *copy);
///@}

/** @name Free memory for structures */
///@{
void free_noise(struct Noise *noise);
void free_model(struct Model *model);
void free_source(struct Source *source);
void free_chain(struct Chain *chain, struct Flags *flags);
void free_calibration(struct Calibration *calibration);
///@}

/**
 \brief Deep comparison of Model contents
 
 @return 0 if models are identical
 @return 1 if models are different
 */
int compare_model(struct Model *a, struct Model *b);


#endif /* GalacticBinaryModel_h */
