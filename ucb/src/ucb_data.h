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
 @file ucb_data.h
 \brief Codes for data handling, including reading external data, simulating internal data, and signal injections.
 
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
void UCBInjectSimulatedSource(struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Source *inj);

/**
 \brief parse VB parameter files and convert got GBMCMC parameters
 */
void GetVerificationBinary(struct Data *data, struct Flags *flags, struct Source *inj, FILE *vbFile);

/**
 \brief Store full contents of input cache file from `--catalog` argument
 */
void UCBLoadCatalogCache(struct Data *data, struct Flags *flags, struct Catalog *catalog);



#endif /* ucb_data_h */
