//
//  GalacticBinaryWrapper.h
//  
//
//  Created by Tyson Littenberg on 1/28/21.
//

/**
 @file GalacticBinaryWrapper.h
 \brief Wrapper functions to GBMCMC sampler.

 Functions for handling the memory, data sharing, and sampling of the GBMCMC block for the global fit analysis.
 */

#ifndef GalacticBinaryWrapper_h
#define GalacticBinaryWrapper_h

/*!
 * \brief Collection of all structures needed for GBMCMC sampler plus metadata.
 *
 * The GBMCMCData structure is a wrapper for all data structures needed for the GBMCMC sampler,
 * as well as metdata tracking the current state of the GBMCMC block in the global fit
 *
 */
struct GBMCMCData
{
    
    int mcmc_step; //!<current step number of GBMCMC sampler
    int status; //!<flag indicating if sampler is still working (0) or finished (1)
    
    int procID; //!<MPI process identifier
    int procID_min; //!<lowest rank MPI process for GBMCMC block
    int procID_max; //!<highest rank MPI process for GBMCMC block
    
    /** @name GBMCMC structures defined in GalacticBinary.h */
     ///@{
    struct Flags *flags;
    struct Orbit *orbit;
    struct Chain *chain;
    struct Data  *data;
    struct Prior *prior;
    struct Proposal **proposal;
    struct Model **trial;
    struct Model **model;
    ///@}

    double cpu_time; //!<CPU time for single block update
};

/**
 \brief Allocate memory for GBMCMC structures
 */
void alloc_gbmcmc_data(struct GBMCMCData *gbmcmc_data, int procID);

/**
 \brief Initialize memory for GB structures
 
 Allocates memory for internals of the data structures,
 initializes values, assigns frequency segment, and
 loads/shares saved ucb catalog file if one is specified
 at the command line with the `--catalog` flag.
 */
void setup_gbmcmc_data(struct GBMCMCData *gbmcmc_data, struct TDI *tdi_full);

/**
 \brief Assigns frequency segment to GBMCMC process
 
 Sets contents in `data` structure determining the bandwidth of the
 analysis window for the GBMCMC process. The settings are determined by the
 `--fmin` and `--samples` command line args and the proces ID number.
 
 The `--fmin` and `--samples` arguments are overridden if the file
 `ucb_frequency_spacing.dat` is found in the run directory.
 */
void setup_frequency_segment(struct GBMCMCData *gbmcmc_data);

/**
 \brief Loads data from analysis bandwidth
 
 Copies subset of TDI data in `tdi_full` to `data` corresponding to analysis window
 */
void select_frequency_segment(struct Data *data, struct TDI *tdi_full);

/**
 \brief Share contents of catalog file to all nodes
 
 Copies contents of the cache file from the `--catalog` argument from the `root` node to all other worker nodes using MPI.
 */
void broadcast_cache(struct Data *data, int root, int procID);

/**
 \brief Initialize state of GBMCMC model
 
 Prepares GBMCMC data structures for sampling, including setting the parallel tempering scheme, proposal distributions, and initial state of the model.
 */
void initialize_gbmcmc_sampler(struct GBMCMCData *gbmcmc_data);

/**
 \brief Evolve state of GBMCMC model
 
 Executes a fixed number of MCMC, RJMCMC, and PTMCMC steps for each of the parallel tempered chains.
 Each chain executes `numSteps` cycles of 100 model updates (90% MCMC, 10% RJCMMC0.
 The chains are updated in parallel on different threads using `OpenMP`. The state of the
 After all of those cycles are completed the PTMCMC swap is tried, the current state is written to file, and the model iteration counter is advanced by `numSteps`.
 
 @returns 0 if sampler is completed
 @returns 1 else
 */
int update_gbmcmc_sampler(struct GBMCMCData *gbmcmc_data);

/**
 \brief Send/Receive GBMCMC params from neighboring segments

 */
void exchange_gbmcmc_source_params(struct GBMCMCData *gbmcmc_data);

/**
 \brief Send/Receive GBMCMC sampler status

 @returns 0 when all GBMCMC samplers are complete
 */
int get_gbmcmc_status(struct GBMCMCData *gbmcmc_data, int Nproc, int root, int procID);

/**
 \brief Print current state of sampler in single line indexed by overall iteration number
 
 */
void print_gbmcmc_state(struct GBMCMCData *gbmcmc_data, FILE *fptr, int counter);

#endif /* GalacticBinaryWrapper_h */

