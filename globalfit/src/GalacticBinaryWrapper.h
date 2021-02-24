//
//  GalacticBinaryWrapper.h
//  
//
//  Created by Tyson Littenberg on 1/28/21.
//

#ifndef GalacticBinaryWrapper_h
#define GalacticBinaryWrapper_h

struct GBMCMCData
{
    int mcmc_step;
    int status;
    
    int procID; //!<MPI process identifier
    int procID_min; //!<lowest rank MPI process for GBMCMC block
    int procID_max; //!<highest rank MPI process for GBMCMC block
    
    struct Flags *flags;
    struct Orbit *orbit;
    struct Chain *chain;
    struct Data  *data;
    struct Prior *prior;
    struct Proposal **proposal;
    struct Model **trial;
    struct Model **model;
};

void alloc_gbmcmc_data(struct GBMCMCData *gbmcmc_data, int procID, int procID_min, int procID_max);

void select_frequency_segment(struct Data *data, struct TDI *tdi_full, int procID);

void get_frequency_segment(struct Data *data, struct TDI *tdi_full, int Nsamples, int root, int procID);

void broadcast_cache(struct Data *data, int root, int procID);

void initialize_gbmcmc_sampler(struct GBMCMCData *gbmcmc_data);

void run_gbmcmc_sampler(struct GBMCMCData *gbmcmc_data);
int update_gbmcmc_sampler(struct GBMCMCData *gbmcmc_data);

void exchange_gbmcmc_source_params(struct GBMCMCData *gbmcmc_data);

int get_gbmcmc_status(struct GBMCMCData *gbmcmc_data, int Nproc, int root, int procID);

#endif /* GalacticBinaryWrapper_h */

