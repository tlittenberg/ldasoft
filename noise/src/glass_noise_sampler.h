/*
 *  Copyright (C) 2023 Tyson B. Littenberg (MSFC-ST12)
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
