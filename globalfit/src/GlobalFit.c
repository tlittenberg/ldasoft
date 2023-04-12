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
#include <gbmcmc.h>
#include <Noise.h>
#include <mbh.h>

#include "GalacticBinaryWrapper.h"
#include "VerificationBinaryWrapper.h"
#include "MBHWrapper.h"
#include "NoiseWrapper.h"

#define N_TDI_CHANNELS 2
#define NMAX 10

struct GlobalFitData
{
    int nUCB;
    int nMBH;
    int nVGB;
    
    struct TDI *tdi_full;
    struct TDI *tdi_store;
    struct TDI *tdi_ucb;
    struct TDI *tdi_vgb;
    struct TDI *tdi_mbh;
    struct Noise *psd;
    
    double block_time;
    double max_block_time;
};

static void alloc_gf_data(struct GlobalFitData *global_fit)
{
    global_fit->tdi_full = malloc(sizeof(struct TDI));
    global_fit->tdi_store = malloc(sizeof(struct TDI));
    global_fit->tdi_vgb = malloc(sizeof(struct TDI));
    global_fit->tdi_ucb = malloc(sizeof(struct TDI));
    global_fit->tdi_mbh = malloc(sizeof(struct TDI));
    global_fit->psd = malloc(sizeof(struct Noise));
}

static void setup_gf_data(struct GlobalFitData *global_fit)
{
    alloc_tdi(global_fit->tdi_store, global_fit->tdi_full->N, global_fit->tdi_full->Nchannel);
    alloc_tdi(global_fit->tdi_ucb, global_fit->tdi_full->N, global_fit->tdi_full->Nchannel);
    alloc_tdi(global_fit->tdi_vgb, global_fit->tdi_full->N, global_fit->tdi_full->Nchannel);
    alloc_tdi(global_fit->tdi_mbh, global_fit->tdi_full->N, global_fit->tdi_full->Nchannel);
    global_fit->tdi_store->delta = global_fit->tdi_full->delta;
    global_fit->tdi_vgb->delta = global_fit->tdi_full->delta;
    global_fit->tdi_ucb->delta = global_fit->tdi_full->delta;
    global_fit->tdi_mbh->delta = global_fit->tdi_full->delta;

    copy_tdi(global_fit->tdi_full, global_fit->tdi_store);
}

static void share_data(struct TDI *tdi_full, int root, int procID)
{
    //first tell all processes how large the dataset is
    MPI_Bcast(&tdi_full->N, 1, MPI_INT, root, MPI_COMM_WORLD);

    //all but root process need to allocate memory for TDI structure
    if(procID!=root) alloc_tdi(tdi_full, tdi_full->N, N_TDI_CHANNELS);
    
    //now broadcast contents of TDI structure
    MPI_Bcast(&tdi_full->delta, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Bcast(tdi_full->X, 2*tdi_full->N, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Bcast(tdi_full->Y, 2*tdi_full->N, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Bcast(tdi_full->Z, 2*tdi_full->N, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Bcast(tdi_full->A, 2*tdi_full->N, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Bcast(tdi_full->E, 2*tdi_full->N, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Bcast(tdi_full->T, 2*tdi_full->N, MPI_DOUBLE, root, MPI_COMM_WORLD);
    
}


static void share_vbmcmc_model(struct GBMCMCData *gbmcmc_data,
                               struct VBMCMCData *vbmcmc_data,
                               struct MBHData *mbh_data,
                               struct GlobalFitData *gf,
                               int root, int procID)
{

    /* get waveforms from vbmcmc sampler and send to root */
    if(procID==1)
    {
        struct Flags *flags = vbmcmc_data->flags;

        for(int i=0; i<gf->tdi_vgb->N*2; i++)
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
            
            for(int i=0; i<data->N*2; i++)
            {
                gf->tdi_vgb->A[i+index] += model->tdi[0]->A[i];
                gf->tdi_vgb->E[i+index] += model->tdi[0]->E[i];
            }
        }
        
        // Send full model back to root
        MPI_Send(gf->tdi_vgb->A, gf->tdi_vgb->N*2, MPI_DOUBLE, root, 0, MPI_COMM_WORLD);
        MPI_Send(gf->tdi_vgb->E, gf->tdi_vgb->N*2, MPI_DOUBLE, root, 1, MPI_COMM_WORLD);
    }
    
    
    /* Root sends vb model to worker nodes */
    if(procID==root)
    {
        MPI_Status status;

        //receive full model from VBMCMC node
        MPI_Recv(gf->tdi_vgb->A, gf->tdi_vgb->N*2, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(gf->tdi_vgb->E, gf->tdi_vgb->N*2, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);
        
        //send full model to GBMCMC nodes
        for(int id=gbmcmc_data->procID_min; id<=gbmcmc_data->procID_max; id++)
        {
            MPI_Send(gf->tdi_vgb->A, gf->tdi_vgb->N*2, MPI_DOUBLE, id, 0, MPI_COMM_WORLD);
            MPI_Send(gf->tdi_vgb->E, gf->tdi_vgb->N*2, MPI_DOUBLE, id, 1, MPI_COMM_WORLD);
        }
        
        //send full model to MBH nodes
        for(int id=mbh_data->procID_min; id<=mbh_data->procID_max; id++)
        {
            MPI_Send(gf->tdi_vgb->A, gf->tdi_vgb->N*2, MPI_DOUBLE, id, 0, MPI_COMM_WORLD);
            MPI_Send(gf->tdi_vgb->E, gf->tdi_vgb->N*2, MPI_DOUBLE, id, 1, MPI_COMM_WORLD);
        }

    }
    
    /* Receive vgb segment at gbmcmc models */
    if(procID>=gbmcmc_data->procID_min && procID<=gbmcmc_data->procID_max)
    {
        MPI_Status status;

        MPI_Recv(gf->tdi_vgb->A, gf->tdi_vgb->N*2, MPI_DOUBLE, root, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(gf->tdi_vgb->E, gf->tdi_vgb->N*2, MPI_DOUBLE, root, 1, MPI_COMM_WORLD, &status);
    }
    
    /* Recieve vgb segment at MBH models */
    if(procID>=mbh_data->procID_min && procID<=mbh_data->procID_max)
    {
        MPI_Status status;

        MPI_Recv(gf->tdi_vgb->A, gf->tdi_vgb->N*2, MPI_DOUBLE, root, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(gf->tdi_vgb->E, gf->tdi_vgb->N*2, MPI_DOUBLE, root, 1, MPI_COMM_WORLD, &status);
    }

}

static void share_gbmcmc_model(struct GBMCMCData *gbmcmc_data,
                               struct VBMCMCData *vbmcmc_data,
                               struct MBHData *mbh_data,
                               struct GlobalFitData *gf,
                               int root, int procID)
{
    if(procID>=gbmcmc_data->procID_min && procID<=gbmcmc_data->procID_max)
    {
        struct Data *data = gbmcmc_data->data;
        struct Chain *chain = gbmcmc_data->chain;
        struct Model *model = gbmcmc_data->model[chain->index[0]];
        
        int N = 2*data->N;
        MPI_Send(&data->qmin, 1, MPI_INT, root, 0, MPI_COMM_WORLD);
        MPI_Send(&N, 1, MPI_INT, root, 1, MPI_COMM_WORLD);
        MPI_Send(model->tdi[0]->A, N, MPI_DOUBLE, root, 2, MPI_COMM_WORLD);
        MPI_Send(model->tdi[0]->E, N, MPI_DOUBLE, root, 3, MPI_COMM_WORLD);
    }
    
    if(procID==root)
    {
        MPI_Status status;
                
        //zero out ucb model
        for(int i=0; i<gf->tdi_ucb->N*2; i++)
        {
            gf->tdi_ucb->A[i] = 0.0;
            gf->tdi_ucb->E[i] = 0.0;
        }

        int procID_min = gbmcmc_data->procID_min;
        int procID_max = gbmcmc_data->procID_max;

        for(int n=procID_min; n<=procID_max; n++)
        {
            int N,qmin;
            MPI_Recv(&qmin, 1, MPI_INT, n, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&N, 1, MPI_INT, n, 1, MPI_COMM_WORLD, &status);
            
            double *A = malloc(N*sizeof(double));
            double *E = malloc(N*sizeof(double));
            
            MPI_Recv(A, N, MPI_DOUBLE, n, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(E, N, MPI_DOUBLE, n, 3, MPI_COMM_WORLD, &status);
            
            for(int i=0; i<N; i++)
            {
                gf->tdi_ucb->A[2*qmin+i] += A[i];
                gf->tdi_ucb->E[2*qmin+i] += E[i];
            }

            free(A);
            free(E);
        }
    }
    
    /* broadcast gbmcmc model to all non-UCB worker nodes */
    if(procID==root)
    {
        //send to VBMCMC node
        for(int id=vbmcmc_data->procID_min; id<=vbmcmc_data->procID_max; id++)
        {
            MPI_Send(gf->tdi_ucb->A, gf->tdi_ucb->N*2, MPI_DOUBLE, id, 0, MPI_COMM_WORLD);
            MPI_Send(gf->tdi_ucb->E, gf->tdi_ucb->N*2, MPI_DOUBLE, id, 1, MPI_COMM_WORLD);
        }
        
        //send to MBH node
        for(int id=mbh_data->procID_min; id<=mbh_data->procID_max; id++)
        {
            MPI_Send(gf->tdi_ucb->A, gf->tdi_ucb->N*2, MPI_DOUBLE, id, 0, MPI_COMM_WORLD);
            MPI_Send(gf->tdi_ucb->E, gf->tdi_ucb->N*2, MPI_DOUBLE, id, 1, MPI_COMM_WORLD);
        }
    }
    
    //recieve
    if(procID>=vbmcmc_data->procID_min && procID<=vbmcmc_data->procID_max)
    {
        MPI_Status status;

        //receive UCB model from root node
        MPI_Recv(gf->tdi_ucb->A, gf->tdi_ucb->N*2, MPI_DOUBLE, root, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(gf->tdi_ucb->E, gf->tdi_ucb->N*2, MPI_DOUBLE, root, 1, MPI_COMM_WORLD, &status);
    }
    
    if(procID>=mbh_data->procID_min && procID<=mbh_data->procID_max)
    {
        MPI_Status status;

        //receive UCB model from root node
        MPI_Recv(gf->tdi_ucb->A, gf->tdi_ucb->N*2, MPI_DOUBLE, root, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(gf->tdi_ucb->E, gf->tdi_ucb->N*2, MPI_DOUBLE, root, 1, MPI_COMM_WORLD, &status);
    }

}

static void share_mbh_model(struct GBMCMCData *gbmcmc_data,
                            struct VBMCMCData *vbmcmc_data,
                            struct MBHData *mbh_data,
                            struct GlobalFitData *global_fit,
                            int root, int procID)
{
    int NF;
    int index;
    
    /* send MBH models to root */
    if(procID >= mbh_data->procID_min && procID <= mbh_data->procID_max)
    {
        NF = mbh_data->het->MM - mbh_data->het->MN;
        index = 2*(int)(mbh_data->data->fmin*mbh_data->data->Tobs);
        
        get_mbh_waveform(mbh_data);

        MPI_Send(&NF, 1, MPI_INT, root, 0, MPI_COMM_WORLD); //number of frequency bins for MBH model
        MPI_Send(&index, 1, MPI_INT, root, 1, MPI_COMM_WORLD); //first bin for MBH model
        MPI_Send(mbh_data->tdi->A+index, 2*NF, MPI_DOUBLE, root, 2, MPI_COMM_WORLD);
        MPI_Send(mbh_data->tdi->E+index, 2*NF, MPI_DOUBLE, root, 3, MPI_COMM_WORLD);
    }
    
    /* combine all incoming MBH waveforms */
    if(procID==root)
    {
        MPI_Status status;
        int N = global_fit->tdi_full->N;
                
        double *A = malloc(N*sizeof(double));
        double *E = malloc(N*sizeof(double));
        
        //zero out mbh model
        for(int i=0; i<global_fit->tdi_mbh->N*2; i++)
        {
            global_fit->tdi_mbh->A[i] = 0.0;
            global_fit->tdi_mbh->E[i] = 0.0;
        }

        for(int n=mbh_data->procID_min; n<=mbh_data->procID_max; n++)
        {
            //zero out the workspace for collecting the worker nodes' waveforms
            for(int i=0; i<N; i++) A[i] = E[i] = 0.0;
            
            //get model from worker node
            MPI_Recv(&NF, 1, MPI_INT, n, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&index, 1, MPI_INT, n, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(A+index, 2*NF, MPI_DOUBLE, n, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(E+index, 2*NF, MPI_DOUBLE, n, 3, MPI_COMM_WORLD, &status);
            
            //remove from joint residual
            for(int i=0; i<2*NF; i++)
            {
                global_fit->tdi_mbh->A[index+i] += A[i+index];
                global_fit->tdi_mbh->E[index+i] += E[i+index];
            }            
        }
        
        free(A);
        free(E);
    }


    
    /* Broadcast joint MBH model to all worker nodes */
    MPI_Bcast(global_fit->tdi_mbh->A, global_fit->tdi_mbh->N*2, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Bcast(global_fit->tdi_mbh->E, global_fit->tdi_mbh->N*2, MPI_DOUBLE, root, MPI_COMM_WORLD);
    
    /* Remove current state of MBH model from joint fit */
    if(procID >= mbh_data->procID_min && procID <= mbh_data->procID_max)
    {
        for(int i=0; i<2*NF; i++)
        {
            global_fit->tdi_mbh->A[index+i] -= mbh_data->tdi->A[index+i];
            global_fit->tdi_mbh->E[index+i] -= mbh_data->tdi->E[index+i];
        }
    }

    /* Broadcast run time of first MBH source (used to scale other updates) */
    MPI_Bcast(&mbh_data->cpu_time, 1, MPI_DOUBLE, mbh_data->procID_min, MPI_COMM_WORLD);


}

static void share_noise_model(struct NoiseData *noise_data, struct GBMCMCData *gbmcmc_data, struct VBMCMCData *vbmcmc_data, struct MBHData *mbh_data, struct GlobalFitData *global_fit, int root, int procID)
{
    
    int ic = 0;
    if(procID==root)
    {
        ic = noise_data->chain->index[0];
        struct Noise *model_psd = noise_data->model[ic]->psd;
        copy_noise(model_psd,global_fit->psd);
    }
    
    MPI_Bcast(global_fit->psd->SnA, global_fit->psd->N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(global_fit->psd->SnE, global_fit->psd->N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

static void create_residual(struct GlobalFitData *global_fit, int GBMCMC_Flag, int VBMCMC_Flag, int MBH_Flag)
{
    
    memcpy(global_fit->tdi_full->A, global_fit->tdi_store->A, 2*global_fit->tdi_full->N*sizeof(double));
    memcpy(global_fit->tdi_full->E, global_fit->tdi_store->E, 2*global_fit->tdi_full->N*sizeof(double));
 
    for(int i=0; i<2*global_fit->tdi_full->N; i++)
    {
        if(global_fit->nUCB>0 && !GBMCMC_Flag)
        {
            global_fit->tdi_full->A[i] -= global_fit->tdi_ucb->A[i];
            global_fit->tdi_full->E[i] -= global_fit->tdi_ucb->E[i];
        }

        if(global_fit->nVGB>0 && !VBMCMC_Flag)
        {
            global_fit->tdi_full->A[i] -= global_fit->tdi_vgb->A[i];
            global_fit->tdi_full->E[i] -= global_fit->tdi_vgb->E[i];
        }

        /*
         MBH nodes hold all other MBH states in their global_fit structure
         */
        if(global_fit->nMBH>0)
        {
            global_fit->tdi_full->A[i] -= global_fit->tdi_mbh->A[i];
            global_fit->tdi_full->E[i] -= global_fit->tdi_mbh->E[i];
        }
    }

}

static void print_data_state(struct NoiseData *noise_data, struct GBMCMCData *gbmcmc_data, struct VBMCMCData *vbmcmc_data, struct MBHData *mbh_data, int GBMCMC_Flag, int VBMCMC_Flag, int Noise_Flag, int MBH_Flag)
{
    if(Noise_Flag)
    {
        print_data(noise_data->data, noise_data->data->tdi[0], noise_data->flags, 0);
        char filename[PATH_BUFSIZE];
        pathprintf(filename,"%s/data/current_spline_points.dat",noise_data->flags->runDir);
        print_noise_model(noise_data->model[noise_data->chain->index[0]]->spline, filename);
        
        pathprintf(filename,"%s/data/current_interpolated_spline_points.dat",noise_data->flags->runDir);
        print_noise_model(noise_data->model[noise_data->chain->index[0]]->psd, filename);
    }
    if(VBMCMC_Flag)
    {
        for(int n=0; n<vbmcmc_data->flags->NVB; n++)print_data(vbmcmc_data->data_vec[n], vbmcmc_data->data_vec[n]->tdi[0], vbmcmc_data->flags, 0);
    }
    if(GBMCMC_Flag)
    {
        copy_noise(gbmcmc_data->model[gbmcmc_data->chain->index[0]]->noise[0], gbmcmc_data->data->noise[0]);
        print_data(gbmcmc_data->data, gbmcmc_data->data->tdi[0], gbmcmc_data->flags, 0);
    }
    if(MBH_Flag)
    {
        char tempFileName[PATH_BUFSIZE];
        pathprintf(tempFileName,"%s/power_data_0.dat.temp",mbh_data->flags->runDir);
        FILE *tempFile = fopen(tempFileName,"w");
        for(int i=mbh_data->het->MN; i<mbh_data->het->MM; i++)
        {
            int re = i;
            int im = mbh_data->data->N-i;
            double f = (double)i/mbh_data->data->Tobs;
            fprintf(tempFile,"%lg %lg %lg\n",f,
                    mbh_data->data->data[0][re]*mbh_data->data->data[0][re]+mbh_data->data->data[0][im]*mbh_data->data->data[0][im],
                    mbh_data->data->data[1][re]*mbh_data->data->data[1][re]+mbh_data->data->data[1][im]*mbh_data->data->data[1][im]);
        }
        fclose(tempFile);
        
        
        pathprintf(tempFileName,"%s/data_0.dat.temp",mbh_data->flags->runDir);
        tempFile = fopen(tempFileName,"w");
        for(int i=0; i<mbh_data->het->MN; i++)
        {
            double f = (double)i/mbh_data->data->Tobs;
            fprintf(tempFile,"%lg %lg %lg %lg %lg\n",f, 0.0,0.0,0.0,0.0);
            
        }
        for(int i=mbh_data->het->MN; i<mbh_data->het->MM; i++)
        {
            int re = 2*(i-mbh_data->het->MN);
            int im = re+1;
            double f = (double)i/mbh_data->data->Tobs;
            fprintf(tempFile,"%lg %lg %lg %lg %lg\n",f,
                    mbh_data->tdi->A[re],mbh_data->tdi->A[im],
                    mbh_data->tdi->E[re],mbh_data->tdi->E[im]);
        }
        for(int i=mbh_data->het->MM; i<mbh_data->data->N; i++)
        {
            double f = (double)i/mbh_data->data->Tobs;
            fprintf(tempFile,"%lg %lg %lg %lg %lg\n",f, 0.0,0.0,0.0,0.0);
            
        }
        
        fclose(tempFile);
        
    }
}

static void blocked_gibbs_load_balancing(struct GlobalFitData *global_fit, int root, int procID, int Nproc)
{
    
    double *block_time_vec=malloc(sizeof(double)*Nproc);
    
    if(procID == root) MPI_Gather(&global_fit->block_time, 1, MPI_DOUBLE, block_time_vec, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
    else               MPI_Gather(&global_fit->block_time, 1, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, root, MPI_COMM_WORLD);
    
    
    if(procID==root)
    {
        global_fit->max_block_time=0;
        for(int n=0; n<Nproc; n++) if(block_time_vec[n]>global_fit->max_block_time) global_fit->max_block_time=block_time_vec[n];
    }
    MPI_Bcast(&global_fit->max_block_time, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
    free(block_time_vec);
}


int main(int argc, char *argv[])
{
    time_t start, stop;
    start = time(NULL);

    int Nproc, procID;
    int root = 0; //root process

    MPI_Init(&argc, &argv); //start parallelization
    MPI_Barrier(MPI_COMM_WORLD);
    
    /* get process ID, and total number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &Nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &procID);

    /* root process check command line */
    if(procID==root)
    {
        //check command line format
        print_LISA_ASCII_art(stdout);
        if(argc==1) print_usage();
    }

    /* Allocate structures to hold global model */
    struct GlobalFitData *global_fit  = malloc(sizeof(struct GlobalFitData));
    struct GBMCMCData    *gbmcmc_data = malloc(sizeof(struct GBMCMCData));
    struct VBMCMCData    *vbmcmc_data = malloc(sizeof(struct VBMCMCData));
    struct NoiseData     *noise_data  = malloc(sizeof(struct NoiseData));
    struct MBHData       *mbh_data    = malloc(sizeof(struct MBHData));

    alloc_gf_data(global_fit);

    /* Allocate GBMCMC data structures */
    alloc_gbmcmc_data(gbmcmc_data, procID);
        
    /* Aliases to gbmcmc structures */
    struct Flags *flags = gbmcmc_data->flags;
    struct Orbit *orbit = gbmcmc_data->orbit;
    struct Chain *chain = gbmcmc_data->chain;
    struct Data  *data  = gbmcmc_data->data;

    /* all processes parse command line and set defaults/flags */
    parse_vb_list(argc, argv, flags);
    parse_mbh_args(argc, argv, mbh_data);
    parse(argc,argv,data,orbit,flags,chain,NMAX);
    
    /* Allocate remaining data structures */
    alloc_noise_data(noise_data, gbmcmc_data, procID, Nproc-1); //noise runs on root process
    alloc_vbmcmc_data(vbmcmc_data, gbmcmc_data, procID); //vbs run on process 1
    alloc_mbh_data(mbh_data, gbmcmc_data, procID); //next are the MBH processes
    
    /* Store size of each model component */
    global_fit->nMBH = mbh_data->NMBH;
    global_fit->nVGB = vbmcmc_data->flags->NVB;
    global_fit->nUCB = Nproc - 1; //room for noise model
    if(global_fit->nMBH>0) global_fit->nUCB -= global_fit->nMBH; //room for MBH mdodels
    if(global_fit->nVGB>0) global_fit->nUCB -= 1; //room for VGB mdodels

    /* Assign processes to models */
    int pid_counter = 0;
    
    //noise model takes one node
    noise_data->procID_min = pid_counter;
    noise_data->procID_max = pid_counter;
    pid_counter++;
    
    //vb model takes one node
    vbmcmc_data->procID_min =  0;
    vbmcmc_data->procID_max = -1;
    if(global_fit->nVGB>0)
    {
        vbmcmc_data->procID_min = pid_counter;
        vbmcmc_data->procID_max = pid_counter;
        pid_counter++;
    }
    
    //mbh model takes one node/source
    mbh_data->procID_min =  0;
    mbh_data->procID_max = -1;
    if(global_fit->nMBH>0)
    {
        mbh_data->procID_min = pid_counter;
        mbh_data->procID_max = pid_counter+mbh_data->NMBH-1;
        pid_counter+=mbh_data->NMBH;
    }
    
    //ucb model takes remaining nodes
    gbmcmc_data->procID_min = pid_counter;
    gbmcmc_data->procID_max = Nproc-1;

    /* Tell the user how the resources are allocated */
    if(procID==root)
    {
        
        fprintf(stdout,"\n =============== Global Fit Analysis ============== \n");
        fprintf(stdout,"  %i noise  processes (pid %i)\n",1+noise_data->procID_max-noise_data->procID_min,noise_data->procID_min);
        if(global_fit->nVGB>0) fprintf(stdout,"  %i vbmcmc processes (pid %i)\n",1+vbmcmc_data->procID_max-vbmcmc_data->procID_min,vbmcmc_data->procID_min);
        if(global_fit->nMBH>0) fprintf(stdout,"  %i mbh    processes (pid %i-%i)\n",1+mbh_data->procID_max-mbh_data->procID_min,mbh_data->procID_min,mbh_data->procID_max);
        fprintf(stdout,"  %i gbmcmc processes (pid %i-%i)\n",1+gbmcmc_data->procID_max-gbmcmc_data->procID_min,gbmcmc_data->procID_min,gbmcmc_data->procID_max);
        fprintf(stdout," ================================================== \n");
   }

    /* Setup output directories for chain and data structures */
    if(procID==0) mkdir(gbmcmc_data->flags->runDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(procID==0)
    {
        //pathprintf(noise_data->flags->runDir,"%s/noise",noise_data->flags->runDir);
        pathappendprintf(noise_data->flags->runDir, "/noise");
        setup_run_directories(noise_data->flags, noise_data->data, noise_data->chain);
        
    }
    else if(procID>=vbmcmc_data->procID_min && procID<=vbmcmc_data->procID_max)
    {
        pathappendprintf(vbmcmc_data->flags->runDir,"/vgb");
        mkdir(vbmcmc_data->flags->runDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        //save original runDir
        char runDir[PATH_BUFSIZE];
        pathprintf(runDir,"%s",vbmcmc_data->flags->runDir);

        for(int n=0; n<vbmcmc_data->flags->NVB; n++)
        {
            //temporarily assign runDir to seg subdir
            pathprintf(vbmcmc_data->flags->runDir,"%s/seg_%04d",runDir,n);
            setup_run_directories(vbmcmc_data->flags, vbmcmc_data->data_vec[n], vbmcmc_data->chain_vec[n]);
        }
        
        //restore original runDir
        pathprintf(vbmcmc_data->flags->runDir,"%s",runDir);

    }
    else if(procID>=gbmcmc_data->procID_min && procID<=gbmcmc_data->procID_max)
    {
        pathappendprintf(gbmcmc_data->flags->runDir,"/ucb");
        mkdir(gbmcmc_data->flags->runDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        pathappendprintf(gbmcmc_data->flags->runDir,"/seg_%04d",procID-gbmcmc_data->procID_min);
        setup_run_directories(gbmcmc_data->flags, gbmcmc_data->data, gbmcmc_data->chain);
    }
    else if(procID>=mbh_data->procID_min && procID <=mbh_data->procID_max)
    {
        pathappendprintf(mbh_data->flags->runDir,"/mbh");
        mkdir(mbh_data->flags->runDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        pathappendprintf(mbh_data->flags->runDir,"/src%04d",procID-mbh_data->procID_min);
        mkdir(mbh_data->flags->runDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    }
       
    if(procID>=gbmcmc_data->procID_min && procID<=gbmcmc_data->procID_max)
    {
        setup_frequency_segment(gbmcmc_data);
    }
    
    /* Initialize gbmcmc data structures */
    alloc_data(data, flags);
            
    /* root process reads data */
    if(procID==root) GalacticBinaryReadHDF5(data,global_fit->tdi_full);
    MPI_Barrier(MPI_COMM_WORLD);

    /* alias of full TDI data */
    struct TDI *tdi_full = global_fit->tdi_full;

    /* broadcast data to all gbmcmc processes and select frequency segment */
    share_data(tdi_full, root, procID);
    
    /* composite models are allocated to be the same size as tdi_full */
    setup_gf_data(global_fit);

    /* set up data for ucb model processes */
    setup_gbmcmc_data(gbmcmc_data, tdi_full);
    
    /* set up data for verification binary model processes */
    if(vbmcmc_data->flags->NVB>0)setup_vbmcmc_data(vbmcmc_data, gbmcmc_data, tdi_full);

    /* set up data for mbh model processes */
    if(mbh_data->NMBH>0)setup_mbh_data(mbh_data, gbmcmc_data, tdi_full, procID);

    /* set up data for noise model processes */
    setup_noise_data(noise_data, gbmcmc_data, vbmcmc_data, mbh_data, tdi_full, procID);

    /* allocate global fit noise model for all processes */
    alloc_noise(global_fit->psd, noise_data->data->N);
    
    /*
     * Initialize all of the samplers
     *
     */
    //choose which sampler to run based on procID
    int GBMCMC_Flag = 0;
    int VBMCMC_Flag = 0;
    int Noise_Flag = 0;
    int MBH_Flag = 0;

    /* Assign processes to Noise model */
    if(procID == 0)
    {
        Noise_Flag = 1;
        initialize_noise_sampler(noise_data);

        int ic = noise_data->chain->index[0];
        struct Noise *model_psd = noise_data->model[ic]->psd;
        copy_noise(model_psd,global_fit->psd);

    }
    MPI_Bcast(global_fit->psd->f, global_fit->psd->N, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Bcast(global_fit->psd->SnX, global_fit->psd->N, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Bcast(global_fit->psd->SnA, global_fit->psd->N, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Bcast(global_fit->psd->SnE, global_fit->psd->N, MPI_DOUBLE, root, MPI_COMM_WORLD);
    share_noise_model(noise_data,gbmcmc_data,vbmcmc_data,mbh_data,global_fit,root,procID);

    /* Assign processes to VB model */
    if(procID>=vbmcmc_data->procID_min && procID<=vbmcmc_data->procID_max)
    {
        VBMCMC_Flag = 1;
        initialize_vbmcmc_sampler(vbmcmc_data);
    }
    
    /* Assign processes to UCB model */
    if(procID >= gbmcmc_data->procID_min && procID <= gbmcmc_data->procID_max)
    {
        GBMCMC_Flag = 1;
        initialize_gbmcmc_sampler(gbmcmc_data);
    }

    /* Assign processes to MBH model */
    if(procID >= mbh_data->procID_min && procID <= mbh_data->procID_max)
    {
        MBH_Flag = 1;
        initialize_mbh_sampler(mbh_data);
    }

    /* distribute current state of models to worker nodes */
    share_noise_model (noise_data, gbmcmc_data, vbmcmc_data, mbh_data, global_fit, root, procID);
    if(global_fit->nUCB>0)share_gbmcmc_model(gbmcmc_data, vbmcmc_data, mbh_data, global_fit, root, procID);
    if(global_fit->nVGB>0)share_vbmcmc_model(gbmcmc_data, vbmcmc_data, mbh_data, global_fit, root, procID);
    if(global_fit->nMBH>0)share_mbh_model   (gbmcmc_data, vbmcmc_data, mbh_data, global_fit, root, procID);

    /* number of update steps for each module (scaled to MBH model update) */
    int cycle;
    global_fit->max_block_time=1;
    /*
     * Master Blocked Gibbs sampler
     *
     */
    do
    {
        cycle=1;
        /* ============================= */
        /*     ULTRACOMPACT BINARIES     */
        /* ============================= */

        /* gbmcmc sampler gibbs update */
        if(GBMCMC_Flag)
        {
            create_residual(global_fit, GBMCMC_Flag, VBMCMC_Flag, MBH_Flag);
            
            select_frequency_segment(gbmcmc_data->data, tdi_full);

            select_noise_segment(global_fit->psd, gbmcmc_data->data, gbmcmc_data->chain, gbmcmc_data->model);
            

            cycle = (int)round(global_fit->max_block_time/gbmcmc_data->cpu_time);

            exchange_gbmcmc_source_params(gbmcmc_data);
            if(procID%2==0)
            {
                for(int i=0; i<((cycle > 1 ) ? cycle : 1); i++)
                    gbmcmc_data->status = update_gbmcmc_sampler(gbmcmc_data);
            }
            
            exchange_gbmcmc_source_params(gbmcmc_data);
            if(procID%2!=0)
            {
                for(int i=0; i<((cycle > 1 ) ? cycle : 1); i++)
                    gbmcmc_data->status = update_gbmcmc_sampler(gbmcmc_data);
            }
            
            global_fit->block_time = gbmcmc_data->cpu_time;
        }

        /* ============================= */
        /*  MASSIVE BLACK HOLE BINARIES  */
        /* ============================= */

        /* mbh model update */
        if(MBH_Flag)
        {
            create_residual(global_fit, GBMCMC_Flag, VBMCMC_Flag, MBH_Flag);
            
            select_mbh_segment(mbh_data, tdi_full);
            
            select_mbh_noise(mbh_data, global_fit->psd);
            
            cycle = (int)round(global_fit->max_block_time/mbh_data->cpu_time);
            for(int i=0; i<((cycle > 1 ) ? cycle : 1); i++)
                mbh_data->status = update_mbh_sampler(mbh_data);
            
            global_fit->block_time = mbh_data->cpu_time;
        }
        
        
        /* ========================================= */
        /* MPI EXCHANGES OF UCB & MBH MODEL STATES */
        /* ========================================= */

        if(global_fit->nUCB>0)share_gbmcmc_model(gbmcmc_data, vbmcmc_data, mbh_data, global_fit, root, procID);
        if(global_fit->nMBH>0)share_mbh_model   (gbmcmc_data, vbmcmc_data, mbh_data, global_fit, root, procID);

        
        /* ============================= */
        /*   VERIFICATION BINARY MODEL   */
        /* ============================= */

        /* vbmcmc sampler gibbs update */
        if(VBMCMC_Flag)
        {
            create_residual(global_fit, GBMCMC_Flag, VBMCMC_Flag, MBH_Flag);

            select_vbmcmc_segments(vbmcmc_data, tdi_full);

            for(int n=0; n<vbmcmc_data->flags->NVB; n++)
                select_noise_segment(global_fit->psd, vbmcmc_data->data_vec[n], vbmcmc_data->chain_vec[n], vbmcmc_data->model_vec[n]);

            //cycle = (int)round(global_fit->max_block_time/vbmcmc_data->cpu_time);
            //for(int i=0; i<((cycle > 1 ) ? cycle : 1); i++)
            vbmcmc_data->status = update_vbmcmc_sampler(vbmcmc_data);
            
            global_fit->block_time = vbmcmc_data->cpu_time;
        }

        /* ============================= */
        /*    INSTRUMENT NOISE MODEL     */
        /* ============================= */

        /* noise model update */
        if(Noise_Flag)
        {
            create_residual(global_fit, GBMCMC_Flag, VBMCMC_Flag, MBH_Flag);

            select_frequency_segment(noise_data->data, tdi_full);

            //cycle = (int)round(global_fit->max_block_time/noise_data->cpu_time);
            //for(int i=0; i<((cycle > 1 ) ? cycle : 1); i++)
            noise_data->status = update_noise_sampler(noise_data);
            
            global_fit->block_time = noise_data->cpu_time;
        }

        /* get global status of gbmcmc samplers */
        gbmcmc_data->status = get_gbmcmc_status(gbmcmc_data,Nproc,root,procID);

        /* ========================================= */
        /* MPI EXCHANGES OF NOISE & VGB MODEL STATES */
        /* ========================================= */
        
        share_noise_model (noise_data, gbmcmc_data, vbmcmc_data, mbh_data, global_fit, root, procID);
        if(global_fit->nVGB>0)share_vbmcmc_model(gbmcmc_data, vbmcmc_data, mbh_data, global_fit, root, procID);

        /* send time spent in each block to root for load balancing */
        blocked_gibbs_load_balancing(global_fit, root, procID, Nproc);

        /* DEBUG */
        print_data_state(noise_data,gbmcmc_data,vbmcmc_data,mbh_data,GBMCMC_Flag,VBMCMC_Flag,Noise_Flag,MBH_Flag);

        
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
        char filename[PATH_BUFSIZE];
        pathprintf(filename,"%s/data/final_spline_points.dat",noise_data->flags->runDir);
        print_noise_model(noise_data->model[noise_data->chain->index[0]]->spline, filename);

        pathprintf(filename,"%s/data/final_interpolated_spline_points.dat",noise_data->flags->runDir);
        print_noise_model(noise_data->model[noise_data->chain->index[0]]->psd, filename);

        print_noise_reconstruction(noise_data->data, noise_data->flags);
    }
    if(MBH_Flag)
    {
        /* waveform reconstructions */
        //print_mbh_waveform_reconstruction(mbh_data);
    }

    
    //print total run time
    stop = time(NULL);

    if(procID==root) printf(" ELAPSED TIME = %g seconds on %i processes\n",(double)(stop-start),Nproc);

    MPI_Finalize();//ends the parallelization

    return 0;
}


