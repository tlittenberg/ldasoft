/*
 * Copyright 2019 Tyson B. Littenberg & Neil J. Cornish
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
