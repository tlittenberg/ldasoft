/*
 *  Copyright (C) 2021 Tyson B. Littenberg (MSFC-ST12), Neil J. Cornish
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
 @file glass_ucb_sampler.h
 \brief Sampling routines for UCB module
 
 Including
 - Fixed dimension galactic binary MCMC
 - Trans-dimension galactic binary RJMCMC
 - Parallel tempering chain exchanges ptmcmc()
 - Fixed dimension noise model MCMC
 - Fixed dimension data model MCMC
 */

#ifndef ucb_sampler_h
#define ucb_sampler_h

/** \brief Parallel tempering exchange
 
 Cycles through all chains and proposes swaps between adjacent pairs.
 */
void ptmcmc(struct Model **model, struct Chain *chain, struct Flags *flags);

/**
 \brief Adaptive temperature spacing
 
 Adjusts temperature spacing between tempered chains according with the goal of having an even acceptance rate between all pairs. The sensitivity of the adjustment asymptotically goes to zero as the sampler approaches the end of the burn-in phase.
 */
void adapt_temperature_ladder(struct Chain *chain, int mcmc);

/**
 \brief Fixed dimension galactic binary MCMC
 
 One fixed-dimension MCMC step for galactic binary model.
 The sampler chooses a source at random from the full model to update.
 @param[in] ic index for which chain is being updated
 */
void ucb_mcmc(struct Orbit *orbit, struct Data *data, struct Model *model, struct Model *trial, struct Chain *chain, struct Flags *flags, struct Prior *prior, struct Proposal **proposal, int ic);

/**
 \brief Trans-dimension galactic binary RJMCMC
 
 One trans-dimension RJMCMC step for galactic binary model.
 The sampler chooses to either add or remove a galactic binary signal to the model.
 @param[in] ic index for which chain is being updated
 */
void ucb_rjmcmc(struct Orbit *orbit, struct Data *data, struct Model *model, struct Model *trial, struct Chain *chain, struct Flags *flags, struct Prior *prior, struct Proposal **proposal, int ic);

/**
 \brief Noise model MCMC
 
 One fixed-dimension MCMC step to update the noise model.
 Currently this only the variance of the orthogonal A and E TDI channels. This will be generalized to include more complete covariance matrices, including non-stationary and non-orthogonal noise.
 @param[in] ic index for which chain is being updated
 */
void noise_model_mcmc(struct Orbit *orbit, struct Data *data, struct Model *model, struct Model *trial, struct Chain *chain, struct Flags *flags, int ic);

/**
 \brief Setup UCB sampler
 
 Get all UCB structures into their initial state.
 */
void initialize_ucb_state(struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Chain *chain, struct Proposal **proposal, struct Model **model, struct Model **trial, struct Source **inj_vec);

#endif /* ucb_sampler_h */
