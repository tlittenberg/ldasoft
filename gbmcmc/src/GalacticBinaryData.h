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
 @file GalacticBinaryData.h
 \brief Codes for data handling, including reading external data, simulating internal data, and signal injections.
 
 \todo There is **a lot** of redundant code in these functions.  That needs to be fixed.
 */

#ifndef GalacticBinaryData_h
#define GalacticBinaryData_h

#include <stdio.h>

/**
 \brief Reads data from external source input using `--data` flag
 */
void GalacticBinaryReadData(struct Data *data, struct Orbit *orbit, struct Flags *flags);

/**
 \brief Reads LDC-formatted HDF5 data using `--h5-data` flag
 */
void GalacticBinaryReadHDF5(struct Data *data, struct TDI *tdi);

/**
 \brief Reads ASCII data using `--data` flag
 */
void GalacticBinaryReadASCII(struct Data *data, struct TDI *tdi);

/**
 \brief Get theoretical noise PSDs for TDI channels
 */
void GalacticBinaryGetNoiseModel(struct Data *data, struct Orbit *orbit, struct Flags *flags);

/**
 \brief Add simulated Gaussian noise realization to data
 */
void GalacticBinaryAddNoise(struct Data *data, struct TDI *tdi);

/**
 \brief Generate noise-only simulated data
 */
void GalacticBinarySimulateData(struct Data *data, struct Orbit *orbit, struct Flags *flags);

/**
 \brief Injection routine for EM-known binaries
 
 EM observations do not provide information about polarization angle \f$\psi\f$ or (currently) initial phase \f$\varphi_0\f$. Those two parameters are missing from the injection files and drawn from their priors.
 
 */
void GalacticBinaryInjectVerificationSource(struct Data *data, struct Orbit *orbit, struct Flags *flags);

/**
 \brief Injection routine for generic binaries
 
 Unlike verification sources, this code expects the polarization angle \f$\psi\f$ or (currently) initial phase \f$\varphi_0\f$ to be included in the parameter files.
 */
void GalacticBinaryInjectSimulatedSource(struct Data *data, struct Orbit *orbit, struct Flags *flags);

/**
 \brief parse VB parameter files and convert got GBMCMC parameters
 */
void GetVerificationBinary(struct Data *data, struct Flags *flags, FILE *vbFile);

/**
 \brief UNUSED: Computes SNR of all injected sources and saves to file.
 */
void GalacticBinaryCatalogSNR(struct Data *data, struct Orbit *orbit, struct Flags *flags);

/**
 \brief Store full contents of input cache file from `--catalog` argument
 */
void GalacticBinaryLoadCatalogCache(struct Data *data, struct Flags *flags);

/**
 \brief Select sources in `--catalog` that are in frequency segment and store metadata
 */
void GalacticBinaryParseCatalogCache(struct Data *data);

/**
 \brief Get data for sources in `--catalog`.
 */
void GalacticBinaryLoadCatalog(struct Data *data);


/**
 \brief Removes known sources which are outside of analysis window but whose signal power may extend into the region of interest.
 
 Signal parameters are read from file picked input using the `--catalog` command line option. Sources with \f$f_0<f_{\rm min}\f$ but \f$f_0 + (B/2)/T > f_{\rm min}\f$ are removed from the data.
 
 Any sources in the catalog with \f$ f_{\rm min}\leq f_0 \leq f_{\rm max} \f$ (including window padding) are not used in the subtraction.
 */
void GalacticBinaryCleanEdges(struct Data *data, struct Orbit *orbit, struct Flags *flags);



#endif /* GalacticBinaryData_h */
