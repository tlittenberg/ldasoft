//
//  VerificationBinaryWrapper.h
//  ldasoft
//
//  Created by Tyson Littenberg on 8/11/21.
//

#ifndef VerificationBinaryWrapper_h
#define VerificationBinaryWrapper_h

struct VGBData
{
    int mcmc_step;
    int status;
    
    int procID; //!<MPI process identifier
    int nProc;  //!<Number of processing segments
    int procID_min; //!<lowest rank MPI process for VGB block
    int procID_max; //!<highest rank MPI process for VGB block

    struct Flags *flags;
    struct Orbit *orbit;
    struct Chain **chain_vec;
    struct Data  **data_vec;
    struct Prior **prior_vec;
    struct Proposal ***proposal_vec;
    struct Model ***trial_vec;
    struct Model ***model_vec;
    struct Source **vgb_vec;
    
    double cpu_time; //!<CPU time for single update
};

void alloc_vgb_data(struct VGBData *vgb_data, struct UCBData *ucb_data, int procID);

void setup_vgb_data(struct VGBData *vgb_data, struct UCBData *ucb_data, struct TDI *tdi_full);

void initialize_vgb_sampler(struct VGBData *vgb_data);

int update_vgb_sampler(struct VGBData *vgb_data);

void select_vgb_segments(struct VGBData *vgb_data, struct TDI *tdi);

void print_vgb_state(struct VGBData *vgb_data, FILE *fptr, int counter);

#endif /* VerificationBinaryWrapper_h */
