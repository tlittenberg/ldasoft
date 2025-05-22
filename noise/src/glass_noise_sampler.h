/*
 * Copyright 2023 Tyson B. Littenberg
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 @file glass_noise_sampler.h
 \brief Sampling routines for NOISE module
 
 Including
 - Trans-dimension spline noise model sampler
 - Fixed-dimension link-level instrument noise sampler
 - Fixed-dimension phenomenological ucb foreground sampler
 */
#ifndef glass_noise_sampler_h
#define glass_noise_sampler_h

/**
 \brief In-place parallel tempering exchange of instrument noise `model` states
 */
void noise_ptmcmc(struct InstrumentModel **model, struct Chain *chain, struct Flags *flags);

/**
 \brief In-place parallel tempering exchange of spline noise `model` states
 */
void spline_ptmcmc(struct SplineModel **model, struct Chain *chain, struct Flags *flags);

/**
 \brief Fixed-dimension update of each parallel tempered spline `model` state
 */
void noise_spline_model_mcmc(struct Orbit *orbit, struct Data *data, struct SplineModel *model, struct Chain *chain, struct Flags *flags, int ic);

/**
 \brief Trans-dimension update of each parallel tempered spline `model` state
 */
void noise_spline_model_rjmcmc(struct Orbit *orbit, struct Data *data, struct SplineModel *model, struct Chain *chain, struct Flags *flags, int ic);

/**
 \brief Fixed-dimension update of each parallel tempered instrument noise `model` state
 */
void noise_instrument_model_mcmc(struct Orbit *orbit, struct Data *data, struct InstrumentModel *model, struct InstrumentModel *trial, struct ForegroundModel *galaxy, struct Noise *psd, struct Chain *chain, struct Flags *flags, int ic);

/**
 \brief Fixed-dimension update of each parallel tempered galactic foreground noise `model` state
 */
void noise_foreground_model_mcmc(struct Data *data, struct InstrumentModel *noise, struct ForegroundModel *model, struct ForegroundModel *trial, struct Noise *psd, struct Chain *chain, struct Flags *flags, int ic);


#endif /* glass_noise_sampler_h */
