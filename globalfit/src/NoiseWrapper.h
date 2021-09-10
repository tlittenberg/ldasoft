//
//  NoiseWrapper.h
//  
//
//  Created by Tyson Littenberg on 2/5/21.
//

#ifndef NoiseWrapper_h
#define NoiseWrapper_h

struct NoiseData
{
    int mcmc_step;
    int status;
    
    int procID; //!<MPI process identifier
    int nProc;  //!<Number of processing segments
    int procID_min; //!<lowest rank MPI process for Noise block
    int procID_max; //!<highest rank MPI process for Noise block

    struct Flags *flags;
    struct Orbit *orbit;
    struct Chain *chain;
    struct Data  *data;
    struct SplineModel **model;
};

void alloc_noise_data(struct NoiseData *noise_data, struct GBMCMCData *gbmcmc_data, int procID, int nProc);
void setup_noise_data(struct NoiseData *noise_data, struct GBMCMCData *gbmcmc_data, struct MBHData *mbh_data, struct TDI *tdi_full, int procID);
void select_noise_segment(struct Noise *psd_full, struct Data *data, struct Chain *chain, struct Model **model);

void initialize_noise_sampler(struct NoiseData *noise_data);
void initialize_noise_state(struct NoiseData *noise_data);

int update_noise_sampler(struct NoiseData *noise_data);

#endif /* NoiseWrapper_h */
