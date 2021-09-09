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
    double **ejump; //!< Scale for Fisher matrix proposals
    double ***evec; //!< Directions for Fisher matrix proposals
    int **cv; //!< Unknown
    int **av; //!< Unknown
    gsl_rng **rvec; //!< Vector of RNG seeds for multithreading PTMCMC

    int NH; //!< Some global variable?
    double logLmax; //!< Max log likelihood
    
    char searchDir[MAXSTRINGSIZE]; //!<store `DIRNAME` containing output from mbh search
    int searchSegment;
    int searchSource;
};

void parse_mbh_args(int argc, char **argv, struct MBHData *data);

void alloc_mbh_data(struct MBHData *mbh_data, struct GBMCMCData *gbmcmc_data, int procID, int procID_min, int procID_max);

void setup_mbh_data(struct MBHData *mbh_data, struct GBMCMCData *gbmcmc_data, struct TDI *tdi_full, int procID);

#endif /* MBHWrapper_h */
