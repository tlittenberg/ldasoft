/*
 * Copyright 2019 Tyson B. Littenberg, Kristen Lackeos, Neil J. Cornish
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
 @file glass_ucb_catalog.h
 \brief Codes for postprocessing `ucb_mcmc` results and formatting for catalog production
 
 */


#ifndef ucb_catalog_h
#define ucb_catalog_h

#include <stdio.h>

/*!
 \brief Prototype structure for catalog of detected sources.
 
 Contains size of the recovered catalog and structures for each sources.
*/
struct Catalog
{
    int N; //!<number of discrete sources in catalog
    struct Entry **entry; //!<discrete catalog entries
};

/*!
 \brief Prototype structure for individual source entries in the catalog.
 
 Contains metadata describing/labeling the source,
 and the full posterior reconstruction (parameters, waveforms, etc.).
*/
struct Entry
{
    int Nchain;             //!<number of chain samples
    char name[128];         //!<source name
    char parent[128];       //!<source parent name
    char path[1024];        //!<path to catalog entry
    struct Source **source; //!<source structure contains parameters, defined in ucb.h
    double *match;          //!<match between sample and ref. source
    double *distance;       //!<metric distance between sample and ref. source
    double evidence;        //!<source evidence
    double SNR;             //!<reference SNR of source
    int i;                  //!<sample containing med. freq.
    int *stepFlag;          //!<array containing on/off for entry at each MCMC step.
    struct GMM *gmm;        //!<Gaussian Mixture Model representation of posterior.
};

/**
 \brief Allocates memory for catalog entry (i.e. indivudual source) and initializes with input source.
 */
void alloc_entry(struct Entry *entry, int IMAX);

/**
 \brief Frees memory for catalog entry (i.e. indivudual source)
 */
void free_entry(struct Entry *entry, int IMAX);

/**
 \brief Allocates memory for new catalog entry (i.e. individual source) without initializing contents.
 */
void create_empty_source(struct Catalog *catalog, int NFFT, int Nchannel);

/**
 \brief Allocates memory for catalog entry (i.e. indivudual source) and initializes with input `sample`.
 */
void create_new_source(struct Catalog *catalog, struct Source *sample, struct Noise *noise, int i, int IMAX, int NFFT, int Nchannel);

/**
\brief Adds input `sample` to existing catalog entry and increments counters.
*/
void append_sample_to_entry(struct Entry *entry, struct Source *sample, int IMAX, int NFFT, int Nchannel);

/**
 \brief computes correlation matrix of all parameters \f$ \rho_{ij} = \frac{\langle (\theta_i - \bar\theta_i)(\theta_j - \bar\theta_j) \rangle}{\sigma_i \sigma_j} \f$
 */
void get_correlation_matrix(struct Data *data, struct Catalog *catalog, int *detection_index, int detections, int IMAX, double **corr);

/**
 \brief Wrapper for using functions in GMM_with_EM.c to represent posterior samples of `entry` as a Gaussian Mixture Model.
 */
int gaussian_mixture_model_wrapper(double **ranges, struct Flags *flags, struct Entry *entry, char *outdir, size_t NMODE, size_t NTHIN, unsigned int *seed, double *BIC);



#endif /* ucb_catalog_h */
