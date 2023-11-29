//
//  ucb_wrapper.h
//
//
//  Created by Tyson Littenberg on 1/28/21.
//

/**
 @file GalacticBinaryWrapper.h
 \brief Wrapper functions to UCB sampler.

 Functions for handling the memory, data sharing, and sampling of the UCB block for the global fit analysis.
 */

#ifndef UltraCompactBinaryWrapper_h
#define UltraCompactBinaryWrapper_h

/*!
 * \brief Collection of all structures needed for UCB sampler plus metadata.
 *
 * The UCBData structure is a wrapper for all data structures needed for the UCB sampler,
 * as well as metdata tracking the current state of the UCB block in the global fit
 *
 */
struct UCBData
{
    
    int mcmc_step; //!<current step number of UCB sampler
    int status; //!<flag indicating if sampler is still working (0) or finished (1)
    
    int procID; //!<MPI process identifier
    int procID_min; //!<lowest rank MPI process for UCB block
    int procID_max; //!<highest rank MPI process for UCB block
    
    /** @name UCB structures defined in GalacticBinary.h */
     ///@{
    struct Flags *flags;
    struct Orbit *orbit;
    struct Chain *chain;
    struct Data  *data;
    struct Prior *prior;
    struct Proposal **proposal;
    struct Model **trial;
    struct Model **model;
    struct Catalog *catalog;
    ///@}

    double cpu_time; //!<CPU time for single block update
};

/**
 \brief Allocate memory for UCB structures
 */
void alloc_ucb_data(struct UCBData *ucb_data, int procID);

/**
 \brief Initialize memory for GB structures
 
 Allocates memory for internals of the data structures,
 initializes values, assigns frequency segment, and
 loads/shares saved ucb catalog file if one is specified
 at the command line with the `--catalog` flag.
 */
void setup_ucb_data(struct UCBData *ucb_data, struct TDI *tdi_full);

/**
 \brief Assigns frequency segment to UCB process
 
 Sets contents in `data` structure determining the bandwidth of the
 analysis window for the UCB process. The settings are determined by the
 `--fmin` and `--samples` command line args and the proces ID number.
 
 The `--fmin` and `--samples` arguments are overridden if the file
 `ucb_frequency_spacing.dat` is found in the run directory.
 */
void setup_frequency_segment(struct UCBData *ucb_data);

/**
 \brief Loads data from analysis bandwidth
 
 Copies subset of TDI data in `tdi_full` to `data` corresponding to analysis window
 */
void select_frequency_segment(struct Data *data, struct TDI *tdi_full);

/**
 \brief Initialize state of UCB model
 
 Prepares UCB data structures for sampling, including setting the parallel tempering scheme, proposal distributions, and initial state of the model.
 */
void initialize_ucb_sampler(struct UCBData *ucb_data);

/**
 \brief Evolve state of UCB model
 
 Executes a fixed number of MCMC, RJMCMC, and PTMCMC steps for each of the parallel tempered chains.
 Each chain executes `numSteps` cycles of 100 model updates (90% MCMC, 10% RJCMMC0.
 The chains are updated in parallel on different threads using `OpenMP`. The state of the
 After all of those cycles are completed the PTMCMC swap is tried, the current state is written to file, and the model iteration counter is advanced by `numSteps`.
 
 @returns 0 if sampler is completed
 @returns 1 else
 */
int update_ucb_sampler(struct UCBData *ucb_data);

/**
 \brief Send/Receive UCB params from neighboring segments

 */
void exchange_ucb_source_params(struct UCBData *ucb_data);

/**
 \brief Send/Receive UCB sampler status

 @returns 0 when all UCB samplers are complete
 */
int get_ucb_status(struct UCBData *ucb_data, int Nproc, int root, int procID);

/**
 \brief Print current state of sampler in single line indexed by overall iteration number
 
 */
void print_ucb_state(struct UCBData *ucb_data, FILE *fptr, int counter);

#endif /* UltraCompactBinaryWrapper_h */

