/*
 * Copyright 2021 Tyson B. Littenberg (MSFC-ST12)
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
 @file glass_mbh_wrapper.h
 \brief Wrapper functions to MBH sampler.

 Functions for handling the memory, data sharing, and sampling of the MBH block for the global fit analysis.
 */

#ifndef MBHWrapper_h
#define MBHWrapper_h

struct MBHData
{
    int mcmc_step;
    int status;
    int flag;

    int procID; //!<MPI process identifier
    int procID_min; //!<lowest rank MPI process for MBH block
    int procID_max; //!<highest rank MPI process for MBH block
    
    struct Flags *flags; //!<Command line flags
    struct TDI *tdi; //!<TDI representation of waveform models
    
    /*!
     LMBH-specific types
     */
    struct MBH_Data *data; //!< Structure containing TDI data and metadata
    struct Het *het; //!< Structure containing everything needed for the heterodyne likelihood
    double *logLx; //!< Current state of likelihood for each chain
    double **paramx; //!< Current state of parameters for each chain
    double **paramy; //!< Proposed state of parameters for each chain
    double **sx; //!< Current state of noise scale factor UNUSED
    double **sy; //!< Proposed state of noise scale factor UNUSED
    double *min; //!< Lower end of uniform priors
    double *max; //!< Upper end of uniform priors
    int *who; //!< Keep track of which chain is at which temperature
    double *heat; //!< Temperature ladder for PTMCMC
    double ***history; //!< Buffer of chain samples for Differential Evolution
    double ***Fisher; //!< Fisher Information Matrix
    double **ejump; //!< Scale for Fisher matrix proposals
    double ***evec; //!< Directions for Fisher matrix proposals
    int **cv; //!< Unknown
    int **av; //!< Unknown
    gsl_rng **rvec; //!< Vector of RNG seeds for multithreading PTMCMC

    int NC; //!< Number of chains
    int NH; //!< Some global variable?
    int NMBH; //!< Number of MBHs found from search in data segment
    double logLmax; //!< Max log likelihood
    
    char searchDir[MAXSTRINGSIZE]; //!<store `DIRNAME` containing output from mbh search
    double **segParams; //!<store parameters from search phase
    
    char chainDir[MAXSTRINGSIZE]; //!<location to print chain files
    FILE *chainFile;
    
    double cpu_time; //!<CPU time for single update
};

void parse_mbh_args(int argc, char **argv, struct MBHData *data);

void alloc_mbh_data(struct MBHData *mbh_data, struct UCBData *gbmcmc_data, int procID);

void select_mbh_segment(struct MBHData *mbh_data, struct TDI *tdi_full);

void select_mbh_noise(struct MBHData *mbh_data, struct Noise *psd);

void setup_mbh_data(struct MBHData *mbh_data, struct UCBData *gbmcmc_data, struct TDI *tdi_full, int procID);

void initialize_mbh_sampler(struct MBHData *mbh_data);

int update_mbh_sampler(struct MBHData *mbh_data);

void get_mbh_waveform(struct MBHData *mbh_data);

void print_mbh_state(struct MBHData *mbh_data, FILE *fptr, int counter);

#endif /* MBHWrapper_h */
