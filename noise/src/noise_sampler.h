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

#ifndef noise_sampler_h
#define noise_sampler_h

#include <stdio.h>
#include "noise.h"

/**
 \brief In-place parallel tempering exchange of `model` states
 */
void spline_ptmcmc(struct SplineModel **model, struct Chain *chain, struct Flags *flags);
void noise_ptmcmc(struct InstrumentModel **model, struct Chain *chain, struct Flags *flags);

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
void noise_instrument_model_mcmc(struct Orbit *orbit, struct Data *data, struct InstrumentModel *model, struct ForegroundModel *galaxy, struct Chain *chain, struct Flags *flags, int ic);

/**
 \brief Fixed-dimension update of each parallel tempered galactic foreground noise `model` state
 */
void noise_foreground_model_mcmc(struct Orbit *orbit, struct Data *data, struct InstrumentModel *noise, struct ForegroundModel *model, struct Chain *chain, struct Flags *flags, int ic);

/**
 \brief Set initial state of spline `model`
 */
void initialize_spline_model(struct Orbit *orbit, struct Data *data, struct SplineModel *model, int Nspline);

/**
 \brief Set initial state of instrument noise `model`
 */
void initialize_instrument_model(struct Orbit *orbit, struct Data *data, struct InstrumentModel *model);

/**
 \brief Set initial state of instrument noise `model`
 */
void initialize_foreground_model(struct Orbit *orbit, struct Data *data, struct ForegroundModel *model);


#endif /* noise_sampler_h */
