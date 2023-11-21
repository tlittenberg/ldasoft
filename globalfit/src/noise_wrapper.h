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
    struct InstrumentModel **inst_model;
    struct ForegroundModel **conf_model;

    double cpu_time; //!<CPU time for single update
};

void alloc_noise_data(struct NoiseData *noise_data, struct UCBData *ucb_data, int procID, int nProc);
void setup_noise_data(struct NoiseData *noise_data, struct UCBData *ucb_data, struct VGBData *vbmcmc_data, struct MBHData *mbh_data, struct TDI *tdi_full, int procID);
void select_noise_segment(struct Noise *psd_full, struct Data *data, struct Chain *chain, struct Model **model);

void initialize_noise_sampler(struct NoiseData *noise_data);
void initialize_noise_state(struct NoiseData *noise_data);
void resume_noise_state(struct NoiseData *noise_data);

int update_noise_sampler(struct NoiseData *noise_data);

/* bad name because of naming conflict in GalacticBinaryIO.c */
void print_noise_state(struct NoiseData *noise_data, FILE *fptr, int counter);

#endif /* NoiseWrapper_h */
