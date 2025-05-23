/*
 * Copyright 2019 Tyson B. Littenberg, Neil J. Cornish
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
 @file glass_ucb_data.h
 \brief Functions for data handling, including reading external data, simulating internal data, and signal injections.
 
 \todo There is **a lot** of redundant code in these functions.  That needs to be fixed.
 */

#ifndef ucb_data_h
#define ucb_data_h

#include <stdio.h>


/**
 \brief Injection routine for ENSEMBLE of EM-known binaries
 
  **THIS FUNCTION AND UCBInjectVerificationSource() NEED TO BE CONSOLIDATED**
 EM observations do not provide information about polarization angle \f$\psi\f$ or (currently) initial phase \f$\varphi_0\f$. Those two parameters are missing from the injection files and drawn from their priors.
 
 */
void UCBInjectVerificationSet(struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Source *inj);

/**
 \brief Injection routine for EM-known binaries
 
 EM observations do not provide information about polarization angle \f$\psi\f$ or (currently) initial phase \f$\varphi_0\f$. Those two parameters are missing from the injection files and drawn from their priors.
 
 */
void UCBInjectVerificationSource(struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Source *inj);

/**
 \brief Injection routine for generic binaries
 
 Unlike verification sources, this code expects the polarization angle \f$\psi\f$ or (currently) initial phase \f$\varphi_0\f$ to be included in the parameter files.
 */
void UCBInjectSimulatedSource(struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Source **inj_vec);

/**
 \brief parse VB parameter files and convert got GBMCMC parameters
 */
void GetVerificationBinary(struct Data *data, struct Flags *flags, struct Source *inj, FILE *vbFile);

/**
 \brief Store full contents of input cache file from `--catalog` argument
 */
void UCBLoadCatalogCache(struct Data *data, struct Flags *flags, struct Catalog *catalog);



#endif /* ucb_data_h */
