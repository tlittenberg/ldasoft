//
//  GlobalFit.c
//  
//
//  Created by Tyson Littenberg on 1/26/21.
//

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>

#include <mpi.h>
#include <omp.h>

#include <LISA.h>

#include <GalacticBinary.h>
#include <GalacticBinaryIO.h>
#include <GalacticBinaryData.h>
#include <GalacticBinaryPrior.h>
#include <GalacticBinaryModel.h>
#include <GalacticBinaryProposal.h>
#include <GalacticBinaryWaveform.h>
#include <GalacticBinaryMCMC.h>

#include <Noise.h>

#include "GalacticBinaryWrapper.h"
#include "NoiseWrapper.h"

#define NMAX 10

static void share_gbmcmc_residual(struct GBMCMCData *gbmcmc_data, struct NoiseData *noise_data, int GBMCMC_Flag, int Noise_Flag, int root)
{
    if(GBMCMC_Flag)
    {
        struct Data *data = gbmcmc_data->data;
        struct Chain *chain = gbmcmc_data->chain;
        struct Model *model = gbmcmc_data->model[chain->index[0]];
        MPI_Send(model->residual[0]->A, 2*data->N, MPI_DOUBLE, root, 0, MPI_COMM_WORLD);
        MPI_Send(model->residual[0]->E, 2*data->N, MPI_DOUBLE, root, 1, MPI_COMM_WORLD);
    }
    if(Noise_Flag)
    {
        MPI_Status status;
        int index = 0;
        int N = gbmcmc_data->data->N*2;
        int procID_min = gbmcmc_data->procID_min;
        int procID_max = gbmcmc_data->procID_max;
        struct Data *data = noise_data->data;
        for(int n=procID_min; n<=procID_max; n++)
        {
            index = 2*(n - procID_min)*(gbmcmc_data->data->N - 2*gbmcmc_data->data->qpad);
            MPI_Recv(data->tdi[0]->A+index, N, MPI_DOUBLE, n, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(data->tdi[0]->E+index, N, MPI_DOUBLE, n, 1, MPI_COMM_WORLD, &status);
        }
    }
}

static void share_noise_model(struct GBMCMCData *gbmcmc_data, struct NoiseData *noise_data, int GBMCMC_Flag, int Noise_Flag, int root)
{
    if(Noise_Flag)
    {
        int index = 0;
        int procID_min = gbmcmc_data->procID_min;
        int procID_max = gbmcmc_data->procID_max;
        
        struct Chain *chain = noise_data->chain;
        struct SplineModel *model = noise_data->model[chain->index[0]];
        for(int n=procID_min; n<=procID_max; n++)
        {
            index = (n - procID_min)*(gbmcmc_data->data->N - 2*gbmcmc_data->data->qpad);
            MPI_Send(model->psd->SnA+index, gbmcmc_data->data->N, MPI_DOUBLE, n, 0, MPI_COMM_WORLD);
            MPI_Send(model->psd->SnE+index, gbmcmc_data->data->N, MPI_DOUBLE, n, 1, MPI_COMM_WORLD);
        }
    }
    if(GBMCMC_Flag)
    {
        MPI_Status status;
        struct Data *data = gbmcmc_data->data;
        struct Chain *chain = gbmcmc_data->chain;
        struct Model *model = gbmcmc_data->model[chain->index[0]];
        
        MPI_Recv(model->noise[0]->SnA, data->N, MPI_DOUBLE, root, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(model->noise[0]->SnE, data->N, MPI_DOUBLE, root, 1, MPI_COMM_WORLD, &status);
        
        //copy new noise parameters to each chain & update PSD
        for(int i=1; i<chain->NC; i++)
        {
            memcpy(gbmcmc_data->model[chain->index[i]]->noise[0]->SnA,model->noise[0]->SnA, data->N*sizeof(double));
            memcpy(gbmcmc_data->model[chain->index[i]]->noise[0]->SnE,model->noise[0]->SnE, data->N*sizeof(double));
        }
    }
}


int main(int argc, char *argv[])
{
    time_t start, stop;
    start = time(NULL);

    int Nproc, procID;
    int root = 0; //root process

    MPI_Init(&argc, &argv); //start parallelization
    
    /* get process ID, and total number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &Nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &procID);

    /* root process check command line */
    if(procID==root)
    {
        //check command line format
        print_LISA_ASCII_art(stdout);
        if(argc==1) print_usage();
        
    } MPI_Barrier(MPI_COMM_WORLD);

    /* Allocate data structures */
    struct GBMCMCData *gbmcmc_data = malloc(sizeof(struct GBMCMCData));
    alloc_gbmcmc_data(gbmcmc_data, procID, 1, Nproc-1);
    
    struct NoiseData *noise_data = malloc(sizeof(struct NoiseData));
    alloc_noise_data(noise_data, procID, Nproc-1);
    
    /* Aliases to gbmcmc structures */
    struct Flags *flags = gbmcmc_data->flags;
    struct Orbit *orbit = gbmcmc_data->orbit;
    struct Chain *chain = gbmcmc_data->chain;
    struct Data *data   = gbmcmc_data->data;

    /* all processes parse command line and set defaults/flags */
    parse(argc,argv,data,orbit,flags,chain,NMAX,Nproc,procID);
    
    /* Finish allocating GBMCMC structures now that we know the number of PT chains */
    gbmcmc_data->proposal = malloc(chain->NP*sizeof(struct Proposal*));
    gbmcmc_data->model = malloc(sizeof(struct Model*)*chain->NC);
    gbmcmc_data->trial = malloc(sizeof(struct Model*)*chain->NC);

    /* Initialize data structures */
    alloc_data(data, flags);
        
    /* TDI structure to hold full dataset */
    struct TDI *tdi_full = malloc(sizeof(struct TDI));

    /* root process reads data */
    if(procID==root) GalacticBinaryReadHDF5(data,tdi_full);

    /* broadcast data to all gbmcmc processes and select frequency segment */
    get_frequency_segment(gbmcmc_data->data, tdi_full, tdi_full->N, root, procID);

    /* set up data for noise model processes */
    if(procID==root) setup_noise_data(noise_data, gbmcmc_data, tdi_full);

    /* Initialize LISA orbit model */
    initialize_orbit(data,orbit,flags);
    
    /* Load gb catalog cache file for proposals/priors */
    if(flags->catalog)
    {
        if(procID==root) GalacticBinaryLoadCatalogCache(data, flags);
        
        broadcast_cache(data, root, procID);
        
        GalacticBinaryParseCatalogCache(data);
        GalacticBinaryLoadCatalog(data);
    }


    /*
     * Initialize all of the samplers
     *
     */

    //choose which sampler to run based on procID
    int GBMCMC_Flag = 0;
    int Noise_Flag = 0;

    /* Assign processes to GBMCMC model */
    if(procID >= gbmcmc_data->procID_min && procID <= gbmcmc_data->procID_max)
    {
        GBMCMC_Flag = 1;
        initialize_gbmcmc_sampler(gbmcmc_data);
        print_gb_catalog_script(flags, data, orbit);
    }
    
    /* ssign processes to Noise model */
    if(procID == root) //TODO: generalize so that noise model doesn't have to be root
    {
        Noise_Flag = 1;
        initialize_noise_sampler(noise_data);
    }
    
    
    /*
     * Master Blocked Gibbs sampler
     *
     */
    do
    {
        /* ============================= */
        /*     ULTRACOMPACT BINARIES     */
        /* ============================= */

        /* gbmcmc sampler gibbs update */
        if(GBMCMC_Flag) gbmcmc_data->status = update_gbmcmc_sampler(gbmcmc_data);

        /* get global status of gbmcmc samplers */
        gbmcmc_data->status = get_gbmcmc_status(gbmcmc_data,Nproc,root,procID);

        /* share gbmcmc residual with other worker nodes */
        share_gbmcmc_residual(gbmcmc_data, noise_data, GBMCMC_Flag, Noise_Flag, root);
        
        /* ============================= */
        /*    INSTRUMENT NOISE MODEL     */
        /* ============================= */

        /* noise model update */
        if(Noise_Flag)
            noise_data->status = update_noise_sampler(noise_data);
        
        /* share noise model with other worker nodes */
        share_noise_model(gbmcmc_data, noise_data, GBMCMC_Flag, Noise_Flag, root);
        
        /* ============================= */
        /*  MASSIVE BLACK HOLE BINARIES  */
        /* ============================= */

        /* mbh model update */

        /* share mbh residual with other worker nodes */

    }while(gbmcmc_data->status!=0);
    
    /*
     * Post processing model components
     *
     */
    if(GBMCMC_Flag)
    {
        /* waveform reconstructions */
        print_waveforms_reconstruction(gbmcmc_data->data, gbmcmc_data->flags);
                
        /* evidence results */
        print_evidence(chain,flags);

        /* flush & close chain files */
        free_chain(chain,flags);
    }
    if(Noise_Flag)  print_noise_reconstruction(noise_data->data, noise_data->flags);
    
    /* clean up */
    free_tdi(tdi_full);

    //print total run time
    stop = time(NULL);

    if(procID==root) printf(" ELAPSED TIME = %g seconds on %i processes\n",(double)(stop-start),Nproc);

    MPI_Finalize();//ends the parallelization

    return 0;
}


