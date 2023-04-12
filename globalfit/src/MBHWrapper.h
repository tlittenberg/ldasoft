//
//  MBHWrapper.h
//  global_fit
//
//  Created by Tyson Littenberg on 9/7/21.
//

#ifndef MBHWrapper_h
#define MBHWrapper_h

#include <stdio.h>

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
    
    char searchDir[PATH_BUFSIZE]; //!<store `DIRNAME` containing output from mbh search
    double **segParams; //!<store parameters from search phase
    
    char chainDir[PATH_BUFSIZE]; //!<location to print chain files
    FILE *chainFile;
    
    double cpu_time; //!<CPU time for single update
};

void parse_mbh_args(int argc, char **argv, struct MBHData *data);

void alloc_mbh_data(struct MBHData *mbh_data, struct GBMCMCData *gbmcmc_data, int procID);

void select_mbh_segment(struct MBHData *mbh_data, struct TDI *tdi_full);

void select_mbh_noise(struct MBHData *mbh_data, struct Noise *psd);

void setup_mbh_data(struct MBHData *mbh_data, struct GBMCMCData *gbmcmc_data, struct TDI *tdi_full, int procID);

void initialize_mbh_sampler(struct MBHData *mbh_data);

int update_mbh_sampler(struct MBHData *mbh_data);

void get_mbh_waveform(struct MBHData *mbh_data);

#endif /* MBHWrapper_h */
