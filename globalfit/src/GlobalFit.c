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
#include "VerificationBinaryWrapper.h"
#include "NoiseWrapper.h"

#define NMAX 10

struct GlobalFitData
{
    struct TDI *tdi_full;
    struct TDI *tdi_store;
    struct TDI *tdi_ucb;
    struct TDI *tdi_vgb;
    struct Noise *psd;
};

static void setup_gf_data(struct GlobalFitData *global_fit)
{
    alloc_tdi(global_fit->tdi_store, global_fit->tdi_full->N, global_fit->tdi_full->Nchannel);
    alloc_tdi(global_fit->tdi_ucb, global_fit->tdi_full->N, global_fit->tdi_full->Nchannel);
    alloc_tdi(global_fit->tdi_vgb, global_fit->tdi_full->N, global_fit->tdi_full->Nchannel);
    global_fit->tdi_store->delta = global_fit->tdi_full->delta;
    global_fit->tdi_vgb->delta = global_fit->tdi_full->delta;
    global_fit->tdi_ucb->delta = global_fit->tdi_full->delta;
    
    copy_tdi(global_fit->tdi_full, global_fit->tdi_store);
}

static void share_vbmcmc_model(struct VBMCMCData *vbmcmc_data, struct GlobalFitData *gf, int root, int procID)
{

    /* get waveforms from vbmcmc sampler and send to root */
    if(procID==1)
    {
        struct Flags *flags = vbmcmc_data->flags;

        for(int i=0; i<gf->tdi_vgb->N; i++)
        {
            gf->tdi_vgb->A[i] = 0.0;
            gf->tdi_vgb->E[i] = 0.0;
        }
        for(int n=0; n<flags->NVB; n++)
        {
            struct Data *data = vbmcmc_data->data_vec[n];
            struct Chain *chain = vbmcmc_data->chain_vec[n];
            struct Model *model = vbmcmc_data->model_vec[n][chain->index[0]];
            
            int index = 2*data->qmin;
            
            for(int i=0; i<data->N; i++)
            {
                gf->tdi_vgb->A[i+index] += model->tdi[0]->A[i];
                gf->tdi_vgb->E[i+index] += model->tdi[0]->E[i];
            }
        }
    }
    
    /* broadcast vbmcmc model to all worker nodes */
    MPI_Bcast(gf->tdi_vgb->A, gf->tdi_vgb->N, MPI_DOUBLE, 1, MPI_COMM_WORLD);
    MPI_Bcast(gf->tdi_vgb->E, gf->tdi_vgb->N, MPI_DOUBLE, 1, MPI_COMM_WORLD);

}

static void share_gbmcmc_model(struct GBMCMCData *gbmcmc_data, int GBMCMC_Flag, struct GlobalFitData *gf, int root, int procID)
{
    if(GBMCMC_Flag)
    {
        struct Data *data = gbmcmc_data->data;
        struct Chain *chain = gbmcmc_data->chain;
        struct Model *model = gbmcmc_data->model[chain->index[0]];
        MPI_Send(model->tdi[0]->A, 2*data->N, MPI_DOUBLE, root, 0, MPI_COMM_WORLD);
        MPI_Send(model->tdi[0]->E, 2*data->N, MPI_DOUBLE, root, 1, MPI_COMM_WORLD);
    }
    
    if(procID==root)
    {
        MPI_Status status;
        int index = 0;
        
        int N = gbmcmc_data->data->N*2;
        
        int procID_min = gbmcmc_data->procID_min;
        int procID_max = gbmcmc_data->procID_max;
        
        int qmin = gbmcmc_data->data->qmin + gbmcmc_data->data->N;
        
        double *A = malloc(N*sizeof(double));
        double *E = malloc(N*sizeof(double));
        
        //zero out ucb model
        for(int i=0; i<gf->tdi_ucb->N; i++)
        {
            gf->tdi_ucb->A[i] = 0.0;
            gf->tdi_ucb->E[i] = 0.0;
        }
        
        for(int n=procID_min; n<=procID_max; n++)
        {
            index = 2*(qmin + (n - procID_min)*(gbmcmc_data->data->N - 2*gbmcmc_data->data->qpad));
            MPI_Recv(A, N, MPI_DOUBLE, n, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(E, N, MPI_DOUBLE, n, 1, MPI_COMM_WORLD, &status);
            
            for(int i=0; i<N; i++)
            {
                gf->tdi_ucb->A[index+i] += A[i];
                gf->tdi_ucb->E[index+i] += E[i];
            }
        }
        
        free(A);
        free(E);
    }
    
    /* broadcast gbmcmc model to all worker nodes */
    MPI_Bcast(gf->tdi_ucb->A, gf->tdi_vgb->N, MPI_DOUBLE, 1, MPI_COMM_WORLD);
    MPI_Bcast(gf->tdi_ucb->E, gf->tdi_vgb->N, MPI_DOUBLE, 1, MPI_COMM_WORLD);

}

static void share_noise_model(struct NoiseData *noise_data, struct GlobalFitData *global_fit, int root, int procID)
{

    int ic = 0;
    if(procID==root)
    {
        ic = noise_data->chain->index[0];
        struct Noise *model_psd = noise_data->model[ic]->psd;
        copy_noise(model_psd,global_fit->psd);

    }
    MPI_Bcast(global_fit->psd->f, global_fit->psd->N, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Bcast(global_fit->psd->SnX, global_fit->psd->N, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Bcast(global_fit->psd->SnA, global_fit->psd->N, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Bcast(global_fit->psd->SnE, global_fit->psd->N, MPI_DOUBLE, root, MPI_COMM_WORLD);
}

static void create_residual(struct GlobalFitData *global_fit, int GBMCMC_Flag, int VBMCMC_Flag)
{
    for(int i=0; i<global_fit->tdi_vgb->N; i++)
    {
        global_fit->tdi_full->A[i] = global_fit->tdi_store->A[i];
        global_fit->tdi_full->E[i] = global_fit->tdi_store->E[i];
        
        if(!GBMCMC_Flag)
        {
            global_fit->tdi_full->A[i] -= global_fit->tdi_ucb->A[i];
            global_fit->tdi_full->E[i] -= global_fit->tdi_ucb->E[i];
        }

        if(!VBMCMC_Flag)
        {
            global_fit->tdi_full->A[i] -= global_fit->tdi_vgb->A[i];
            global_fit->tdi_full->E[i] -= global_fit->tdi_vgb->E[i];
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

    /* Allocate GBMCMC data structures */
    struct GBMCMCData *gbmcmc_data = malloc(sizeof(struct GBMCMCData));
    alloc_gbmcmc_data(gbmcmc_data, procID, 2, Nproc-1);
        
    /* Aliases to gbmcmc structures */
    struct Flags *flags = gbmcmc_data->flags;
    struct Orbit *orbit = gbmcmc_data->orbit;
    struct Chain *chain = gbmcmc_data->chain;
    struct Data *data   = gbmcmc_data->data;

    /* all processes parse command line and set defaults/flags */
    parse_vb_list(argc, argv, flags);
    parse(argc,argv,data,orbit,flags,chain,NMAX,procID);
    
    /* Allocate remaining data structures */
    struct VBMCMCData *vbmcmc_data = malloc(sizeof(struct VBMCMCData));
    alloc_vbmcmc_data(vbmcmc_data, gbmcmc_data, procID);

    struct NoiseData *noise_data = malloc(sizeof(struct NoiseData));
    alloc_noise_data(noise_data, gbmcmc_data, procID, Nproc-1);

    /* Setup output directories for chain and data structures */
    if(procID==0)
    {
        sprintf(noise_data->flags->runDir,"%s_noise",noise_data->flags->runDir);
        sprintf(noise_data->data->dataDir,"%s/data",noise_data->flags->runDir);
        sprintf(noise_data->chain->chainDir,"%s/chains",noise_data->flags->runDir);
        sprintf(noise_data->chain->chkptDir,"%s/checkpoint",noise_data->flags->runDir);

        setup_run_directories(noise_data->flags, noise_data->data, noise_data->chain);
    }
    else if(procID==1)
    {
        sprintf(vbmcmc_data->flags->runDir,"%s_vgb",vbmcmc_data->flags->runDir);
        for(int n=0; n<vbmcmc_data->flags->NVB; n++)
        {
            sprintf(vbmcmc_data->data_vec[n]->dataDir,"%s/data_%i",vbmcmc_data->flags->runDir,n);
            sprintf(vbmcmc_data->chain_vec[n]->chainDir,"%s/chains_%i",vbmcmc_data->flags->runDir,n);
            sprintf(vbmcmc_data->chain_vec[n]->chkptDir,"%s/checkpoint_%i",vbmcmc_data->flags->runDir,n);
            setup_run_directories(vbmcmc_data->flags, vbmcmc_data->data_vec[n], vbmcmc_data->chain_vec[n]);
        }
    }
    else
    {
        sprintf(gbmcmc_data->flags->runDir,"%s_ucb_%i",gbmcmc_data->flags->runDir, procID-1);
        sprintf(gbmcmc_data->data->dataDir,"%s/data",gbmcmc_data->flags->runDir);
        sprintf(gbmcmc_data->chain->chainDir,"%s/chains",gbmcmc_data->flags->runDir);
        sprintf(gbmcmc_data->chain->chkptDir,"%s/checkpoint",gbmcmc_data->flags->runDir);
        setup_run_directories(gbmcmc_data->flags, gbmcmc_data->data, gbmcmc_data->chain);
    }
    /* Finish allocating GBMCMC structures now that we know the number of PT chains */
    gbmcmc_data->proposal = malloc(chain->NP*sizeof(struct Proposal*));
    gbmcmc_data->model = malloc(sizeof(struct Model*)*chain->NC);
    gbmcmc_data->trial = malloc(sizeof(struct Model*)*chain->NC);

    /* Initialize data structures */
    alloc_data(data, flags);
        
    /* TDI structure to hold full dataset & composite models */
    struct GlobalFitData *global_fit = malloc(sizeof(struct GlobalFitData));
    global_fit->tdi_full = malloc(sizeof(struct TDI));
    global_fit->tdi_store = malloc(sizeof(struct TDI));
    global_fit->tdi_vgb = malloc(sizeof(struct TDI));
    global_fit->tdi_ucb = malloc(sizeof(struct TDI));
    global_fit->psd = malloc(sizeof(struct Noise));
    
    /* root process reads data */
    if(procID==root) GalacticBinaryReadHDF5(data,global_fit->tdi_full);
    MPI_Barrier(MPI_COMM_WORLD);

    /* alias of full TDI data */
    struct TDI *tdi_full = global_fit->tdi_full;

    /* broadcast data to all gbmcmc processes and select frequency segment */
    get_frequency_segment(gbmcmc_data->data, tdi_full, tdi_full->N, root, procID, gbmcmc_data->procID_min);

    /* composite models are allocated to be the same size as tdi_full */
    setup_gf_data(global_fit);

    /* set up data for noise model processes */
    setup_noise_data(noise_data, gbmcmc_data, tdi_full, procID);
    
    /* allocate global fit noise model for all processes */
    alloc_noise(global_fit->psd, noise_data->data->N);

    /* set up data for verification binary model processes */
    setup_vbmcmc_data(vbmcmc_data, gbmcmc_data, tdi_full);

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
    int VBMCMC_Flag = 0;
    int Noise_Flag = 0;
    
    /* Assign processes to Noise model */
    if(procID == 0) //TODO: generalize so that noise model doesn't have to be root
    {
        Noise_Flag = 1;
        initialize_noise_sampler(noise_data);
    }
    share_noise_model(noise_data,global_fit,root,procID);

    /* Assign processes to VB model */
    if(procID == 1 && vbmcmc_data->flags->NVB>0)
    {
        VBMCMC_Flag = 1;
        initialize_vbmcmc_sampler(vbmcmc_data);
    }
    
    /* Assign processes to UCB model */
    if(procID >= gbmcmc_data->procID_min && procID <= gbmcmc_data->procID_max)
    {
        GBMCMC_Flag = 1;
        initialize_gbmcmc_sampler(gbmcmc_data);

        print_gb_catalog_script(flags, data, orbit);
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
        if(GBMCMC_Flag)
        {
            
            create_residual(global_fit, GBMCMC_Flag, VBMCMC_Flag);
            
            select_frequency_segment(gbmcmc_data->data, tdi_full);
            
            select_noise_segment(global_fit->psd, gbmcmc_data->model[0]->noise[0]);
            
            
            for(int i=1; i<gbmcmc_data->chain->NC; i++)
            {
                memcpy(gbmcmc_data->model[gbmcmc_data->chain->index[i]]->noise[0]->SnA,gbmcmc_data->model[0]->noise[0]->SnA, gbmcmc_data->data->N*sizeof(double));
                memcpy(gbmcmc_data->model[gbmcmc_data->chain->index[i]]->noise[0]->SnE,gbmcmc_data->model[0]->noise[0]->SnE, gbmcmc_data->data->N*sizeof(double));
            }
            
            
            gbmcmc_data->status = update_gbmcmc_sampler(gbmcmc_data);
        }

        /* get global status of gbmcmc samplers */
        gbmcmc_data->status = get_gbmcmc_status(gbmcmc_data,Nproc,root,procID);

        /* share gbmcmc residual with root node */
        share_gbmcmc_model(gbmcmc_data, GBMCMC_Flag, global_fit, root, procID);
           
        /* ============================= */
        /*   VERIFICATION BINARY MODEL   */
        /* ============================= */

        /* vbmcmc sampler gibbs update */
        if(VBMCMC_Flag)
        {
            create_residual(global_fit, GBMCMC_Flag, VBMCMC_Flag);

            for(int n=0; n<vbmcmc_data->flags->NVB; n++)
            {
                select_noise_segment(global_fit->psd, vbmcmc_data->model_vec[n][0]->noise[0]);
                for(int i=1; i<vbmcmc_data->chain_vec[n]->NC; i++)
                {
                    memcpy(vbmcmc_data->model_vec[n][vbmcmc_data->chain_vec[n]->index[i]]->noise[0]->SnA,vbmcmc_data->model_vec[n][0]->noise[0]->SnA, vbmcmc_data->data_vec[n]->N*sizeof(double));
                    memcpy(vbmcmc_data->model_vec[n][vbmcmc_data->chain_vec[n]->index[i]]->noise[0]->SnE,vbmcmc_data->model_vec[n][0]->noise[0]->SnE, vbmcmc_data->data_vec[n]->N*sizeof(double));
                }
            }

            select_vbmcmc_segments(vbmcmc_data, tdi_full);
            vbmcmc_data->status = update_vbmcmc_sampler(vbmcmc_data);
        }
        
        /* share vbmcmc residual with other worker nodes */
        share_vbmcmc_model(vbmcmc_data, global_fit, root, procID);

        /* ============================= */
        /*    INSTRUMENT NOISE MODEL     */
        /* ============================= */

        /* noise model update */
        if(Noise_Flag)
        {
            create_residual(global_fit, GBMCMC_Flag, VBMCMC_Flag);

            select_frequency_segment(noise_data->data, tdi_full);
            noise_data->status = update_noise_sampler(noise_data);
        }
        /* share noise model with other worker nodes */
        share_noise_model(noise_data, global_fit, root, procID);

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
    }
    if(VBMCMC_Flag)
    {
        for(int n=0; n<vbmcmc_data->flags->NVB; n++)
        {
            /* waveform reconstructions */
            print_waveforms_reconstruction(vbmcmc_data->data_vec[n], vbmcmc_data->flags);
        }
    }
    if(Noise_Flag)
    {
        char filename[128];
        sprintf(filename,"%s/data/final_spline_points.dat",noise_data->flags->runDir);
        print_noise_model(noise_data->model[noise_data->chain->index[0]]->spline, filename);

        sprintf(filename,"%s/data/final_interpolated_spline_points.dat",noise_data->flags->runDir);
        print_noise_model(noise_data->model[noise_data->chain->index[0]]->psd, filename);

        print_noise_reconstruction(noise_data->data, noise_data->flags);
    }

    
    //print total run time
    stop = time(NULL);

    if(procID==root) printf(" ELAPSED TIME = %g seconds on %i processes\n",(double)(stop-start),Nproc);

    MPI_Finalize();//ends the parallelization

    return 0;
}


