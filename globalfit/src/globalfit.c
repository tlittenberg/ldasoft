//
//  GlobalFit.c
//  
//
//  Created by Tyson Littenberg on 1/26/21.
//

#include <mpi.h>

#include <glass_utils.h>
#include <glass_ucb.h>
#include <glass_noise.h>
#include <mbh.h>

#include "glass_ucb_wrapper.h"
#include "glass_vgb_wrapper.h"
#include "glass_mbh_wrapper.h"
#include "glass_noise_wrapper.h"

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
    
    FILE *chainFile;
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

static void dump_data(struct Data *data, struct Flags *flags)
{
    FILE *fptr;
    char filename[MAXSTRINGSIZE];
    struct TDI *tdi = data->tdi;

    sprintf(filename,"%s/data/power_data.dat",flags->runDir);
    fptr = fopen(filename,"w");
    for(int i=0; i<data->N; i++)
    {
        double f = (double)(i+data->qmin)/data->T;
        switch(tdi->Nchannel)
        {
            case 1:
                fprintf(fptr,"%.12g %lg\n", f, tdi->X[2*i]*tdi->X[2*i]+tdi->X[2*i+1]*tdi->X[2*i+1]);
                break;
            case 2:
                fprintf(fptr,"%.12g %lg %lg\n", f, tdi->A[2*i]*tdi->A[2*i]+tdi->A[2*i+1]*tdi->A[2*i+1], tdi->E[2*i]*tdi->E[2*i]+tdi->E[2*i+1]*tdi->E[2*i+1]);
                break;
            case 3:
                fprintf(fptr,"%.12g %lg %lg %lg\n", f, tdi->X[2*i]*tdi->X[2*i]+tdi->X[2*i+1]*tdi->X[2*i+1], tdi->Y[2*i]*tdi->Y[2*i]+tdi->Y[2*i+1]*tdi->Y[2*i+1], tdi->Z[2*i]*tdi->Z[2*i]+tdi->Z[2*i+1]*tdi->Z[2*i+1]);
                break;
        }
    }
    fclose(fptr);

    sprintf(filename,"%s/data/data.dat",flags->runDir);
    fptr = fopen(filename,"w");
    for(int i=0; i<tdi->N; i++)
    {
        double f = (double)(i+data->qmin)/data->T;
        switch(data->Nchannel)
        {
            case 1:
                fprintf(fptr,"%.12g %lg %lg\n", f, tdi->X[2*i],tdi->X[2*i+1]);
                break;
            case 2:
                fprintf(fptr,"%.12g %lg %lg %lg %lg\n", f, tdi->A[2*i],tdi->A[2*i+1], tdi->E[2*i],tdi->E[2*i+1]);
                break;
            case 3:
                fprintf(fptr,"%.12g %lg %lg %lg %lg %lg %lg\n", f, tdi->X[2*i],tdi->X[2*i+1], tdi->Y[2*i],tdi->Y[2*i+1], tdi->Z[2*i],tdi->Z[2*i+1]);
                break;
        }
    }
    fclose(fptr);
}

static void share_vgb_model(struct UCBData *ucb_data,
                               struct VGBData *vgb_data,
                               struct MBHData *mbh_data,
                               struct GlobalFitData *gf,
                               int root, int procID)
{

    /* get waveforms from vgb sampler and send to root */
    if(procID==1)
    {
        struct Flags *flags = vgb_data->flags;

        for(int i=0; i<gf->tdi_vgb->N*2; i++)
        {
            gf->tdi_vgb->X[i] = 0.0;
            gf->tdi_vgb->Y[i] = 0.0;
            gf->tdi_vgb->Z[i] = 0.0;
        }
        for(int n=0; n<flags->NVB; n++)
        {
            struct Data *data = vgb_data->data_vec[n];
            struct Chain *chain = vgb_data->chain_vec[n];
            struct Model *model = vgb_data->model_vec[n][chain->index[0]];
            
            int index = 2*data->qmin;
            
            for(int i=0; i<data->N*2; i++)
            {
                gf->tdi_vgb->X[i+index] += model->tdi->X[i];
                gf->tdi_vgb->Y[i+index] += model->tdi->Y[i];
                gf->tdi_vgb->Z[i+index] += model->tdi->Z[i];
            }
        }
        
        // Send full model back to root
        MPI_Send(gf->tdi_vgb->X, gf->tdi_vgb->N*2, MPI_DOUBLE, root, 0, MPI_COMM_WORLD);
        MPI_Send(gf->tdi_vgb->Y, gf->tdi_vgb->N*2, MPI_DOUBLE, root, 1, MPI_COMM_WORLD);
        MPI_Send(gf->tdi_vgb->Z, gf->tdi_vgb->N*2, MPI_DOUBLE, root, 2, MPI_COMM_WORLD);
    }
    
    
    /* Root sends vb model to worker nodes */
    if(procID==root)
    {
        MPI_Status status;

        //receive full model from VGB node
        MPI_Recv(gf->tdi_vgb->X, gf->tdi_vgb->N*2, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(gf->tdi_vgb->Y, gf->tdi_vgb->N*2, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(gf->tdi_vgb->Z, gf->tdi_vgb->N*2, MPI_DOUBLE, 1, 2, MPI_COMM_WORLD, &status);

        //send full model to UCB nodes
        for(int id=ucb_data->procID_min; id<=ucb_data->procID_max; id++)
        {
            MPI_Send(gf->tdi_vgb->X, gf->tdi_vgb->N*2, MPI_DOUBLE, id, 0, MPI_COMM_WORLD);
            MPI_Send(gf->tdi_vgb->Y, gf->tdi_vgb->N*2, MPI_DOUBLE, id, 1, MPI_COMM_WORLD);
            MPI_Send(gf->tdi_vgb->Z, gf->tdi_vgb->N*2, MPI_DOUBLE, id, 2, MPI_COMM_WORLD);
        }
        
        //send full model to MBH nodes
        for(int id=mbh_data->procID_min; id<=mbh_data->procID_max; id++)
        {
            MPI_Send(gf->tdi_vgb->X, gf->tdi_vgb->N*2, MPI_DOUBLE, id, 0, MPI_COMM_WORLD);
            MPI_Send(gf->tdi_vgb->Y, gf->tdi_vgb->N*2, MPI_DOUBLE, id, 1, MPI_COMM_WORLD);
            MPI_Send(gf->tdi_vgb->Z, gf->tdi_vgb->N*2, MPI_DOUBLE, id, 2, MPI_COMM_WORLD);
        }

    }
    
    /* Receive vgb segment at ucb models */
    if(procID>=ucb_data->procID_min && procID<=ucb_data->procID_max)
    {
        MPI_Status status;

        MPI_Recv(gf->tdi_vgb->X, gf->tdi_vgb->N*2, MPI_DOUBLE, root, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(gf->tdi_vgb->Y, gf->tdi_vgb->N*2, MPI_DOUBLE, root, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(gf->tdi_vgb->Z, gf->tdi_vgb->N*2, MPI_DOUBLE, root, 2, MPI_COMM_WORLD, &status);
    }
    
    /* Recieve vgb segment at MBH models */
    if(procID>=mbh_data->procID_min && procID<=mbh_data->procID_max)
    {
        MPI_Status status;

        MPI_Recv(gf->tdi_vgb->X, gf->tdi_vgb->N*2, MPI_DOUBLE, root, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(gf->tdi_vgb->Y, gf->tdi_vgb->N*2, MPI_DOUBLE, root, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(gf->tdi_vgb->Z, gf->tdi_vgb->N*2, MPI_DOUBLE, root, 2, MPI_COMM_WORLD, &status);
    }

}

static void share_ucb_model(struct UCBData *ucb_data,
                               struct VGBData *vgb_data,
                               struct MBHData *mbh_data,
                               struct GlobalFitData *gf,
                               int root, int procID)
{
    if(procID>=ucb_data->procID_min && procID<=ucb_data->procID_max)
    {
        struct Data *data = ucb_data->data;
        struct Chain *chain = ucb_data->chain;
        struct Model *model = ucb_data->model[chain->index[0]];
        
        int N = 2*data->N;
        MPI_Send(&data->qmin, 1, MPI_INT, root, 0, MPI_COMM_WORLD);
        MPI_Send(&N, 1, MPI_INT, root, 1, MPI_COMM_WORLD);
        MPI_Send(model->tdi->X, N, MPI_DOUBLE, root, 2, MPI_COMM_WORLD);
        MPI_Send(model->tdi->Y, N, MPI_DOUBLE, root, 3, MPI_COMM_WORLD);
        MPI_Send(model->tdi->Z, N, MPI_DOUBLE, root, 4, MPI_COMM_WORLD);
    }
    
    if(procID==root)
    {
        MPI_Status status;
                
        //zero out ucb model
        for(int i=0; i<gf->tdi_ucb->N*2; i++)
        {
            gf->tdi_ucb->X[i] = 0.0;
            gf->tdi_ucb->Y[i] = 0.0;
            gf->tdi_ucb->Z[i] = 0.0;
        }

        int procID_min = ucb_data->procID_min;
        int procID_max = ucb_data->procID_max;

        for(int n=procID_min; n<=procID_max; n++)
        {
            int N,qmin;
            MPI_Recv(&qmin, 1, MPI_INT, n, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&N, 1, MPI_INT, n, 1, MPI_COMM_WORLD, &status);
            
            double *X = malloc(N*sizeof(double));
            double *Y = malloc(N*sizeof(double));
            double *Z = malloc(N*sizeof(double));

            MPI_Recv(X, N, MPI_DOUBLE, n, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(Y, N, MPI_DOUBLE, n, 3, MPI_COMM_WORLD, &status);
            MPI_Recv(Z, N, MPI_DOUBLE, n, 4, MPI_COMM_WORLD, &status);

            for(int i=0; i<N; i++)
            {
                gf->tdi_ucb->X[2*qmin+i] += X[i];
                gf->tdi_ucb->Y[2*qmin+i] += Y[i];
                gf->tdi_ucb->Z[2*qmin+i] += Z[i];
            }

            free(X);
            free(Y);
            free(Z);
        }
    }
    
    /* broadcast ucb model to all non-UCB worker nodes */
    if(procID==root)
    {
        //send to VGB node
        for(int id=vgb_data->procID_min; id<=vgb_data->procID_max; id++)
        {
            MPI_Send(gf->tdi_ucb->X, gf->tdi_ucb->N*2, MPI_DOUBLE, id, 0, MPI_COMM_WORLD);
            MPI_Send(gf->tdi_ucb->Y, gf->tdi_ucb->N*2, MPI_DOUBLE, id, 1, MPI_COMM_WORLD);
            MPI_Send(gf->tdi_ucb->Z, gf->tdi_ucb->N*2, MPI_DOUBLE, id, 2, MPI_COMM_WORLD);
        }
        
        //send to MBH node
        for(int id=mbh_data->procID_min; id<=mbh_data->procID_max; id++)
        {
            MPI_Send(gf->tdi_ucb->X, gf->tdi_ucb->N*2, MPI_DOUBLE, id, 0, MPI_COMM_WORLD);
            MPI_Send(gf->tdi_ucb->Y, gf->tdi_ucb->N*2, MPI_DOUBLE, id, 1, MPI_COMM_WORLD);
            MPI_Send(gf->tdi_ucb->Z, gf->tdi_ucb->N*2, MPI_DOUBLE, id, 2, MPI_COMM_WORLD);
        }
    }
    
    //recieve
    if(procID>=vgb_data->procID_min && procID<=vgb_data->procID_max)
    {
        MPI_Status status;

        //receive UCB model from root node
        MPI_Recv(gf->tdi_ucb->X, gf->tdi_ucb->N*2, MPI_DOUBLE, root, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(gf->tdi_ucb->Y, gf->tdi_ucb->N*2, MPI_DOUBLE, root, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(gf->tdi_ucb->Z, gf->tdi_ucb->N*2, MPI_DOUBLE, root, 2, MPI_COMM_WORLD, &status);
    }
    
    if(procID>=mbh_data->procID_min && procID<=mbh_data->procID_max)
    {
        MPI_Status status;

        //receive UCB model from root node
        MPI_Recv(gf->tdi_ucb->X, gf->tdi_ucb->N*2, MPI_DOUBLE, root, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(gf->tdi_ucb->Y, gf->tdi_ucb->N*2, MPI_DOUBLE, root, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(gf->tdi_ucb->Z, gf->tdi_ucb->N*2, MPI_DOUBLE, root, 2, MPI_COMM_WORLD, &status);
    }

}

static void share_mbh_model(struct UCBData *ucb_data,
                            struct VGBData *vgb_data,
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
        MPI_Send(mbh_data->tdi->X+index, 2*NF, MPI_DOUBLE, root, 2, MPI_COMM_WORLD);
        MPI_Send(mbh_data->tdi->Y+index, 2*NF, MPI_DOUBLE, root, 3, MPI_COMM_WORLD);
        MPI_Send(mbh_data->tdi->Z+index, 2*NF, MPI_DOUBLE, root, 4, MPI_COMM_WORLD);
    }
    
    /* combine all incoming MBH waveforms */
    if(procID==root)
    {
        MPI_Status status;
        int N = global_fit->tdi_full->N;
                
        double *X = malloc(N*sizeof(double));
        double *Y = malloc(N*sizeof(double));
        double *Z = malloc(N*sizeof(double));

        //zero out mbh model
        for(int i=0; i<global_fit->tdi_mbh->N*2; i++)
        {
            global_fit->tdi_mbh->X[i] = 0.0;
            global_fit->tdi_mbh->Y[i] = 0.0;
            global_fit->tdi_mbh->Z[i] = 0.0;
        }

        for(int n=mbh_data->procID_min; n<=mbh_data->procID_max; n++)
        {
            //zero out the workspace for collecting the worker nodes' waveforms
            for(int i=0; i<N; i++) X[i] = Y[i] = Z[i] = 0.0;
            
            //get model from worker node
            MPI_Recv(&NF, 1, MPI_INT, n, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&index, 1, MPI_INT, n, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(X+index, 2*NF, MPI_DOUBLE, n, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(Y+index, 2*NF, MPI_DOUBLE, n, 3, MPI_COMM_WORLD, &status);
            MPI_Recv(Z+index, 2*NF, MPI_DOUBLE, n, 4, MPI_COMM_WORLD, &status);

            //remove from joint residual
            for(int i=0; i<2*NF; i++)
            {
                global_fit->tdi_mbh->X[index+i] += X[i+index];
                global_fit->tdi_mbh->Y[index+i] += Y[i+index];
                global_fit->tdi_mbh->Z[index+i] += Z[i+index];
            }
        }
        
        free(X);
        free(Y);
        free(Z);
    }


    
    /* Broadcast joint MBH model to all worker nodes */
    MPI_Bcast(global_fit->tdi_mbh->X, global_fit->tdi_mbh->N*2, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Bcast(global_fit->tdi_mbh->Y, global_fit->tdi_mbh->N*2, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Bcast(global_fit->tdi_mbh->Z, global_fit->tdi_mbh->N*2, MPI_DOUBLE, root, MPI_COMM_WORLD);

    /* Remove current state of MBH model from joint fit */
    if(procID >= mbh_data->procID_min && procID <= mbh_data->procID_max)
    {
        for(int i=0; i<2*NF; i++)
        {
            global_fit->tdi_mbh->X[index+i] -= mbh_data->tdi->X[index+i];
            global_fit->tdi_mbh->Y[index+i] -= mbh_data->tdi->Y[index+i];
            global_fit->tdi_mbh->Z[index+i] -= mbh_data->tdi->Z[index+i];

        }
    }

    /* Broadcast run time of first MBH source (used to scale other updates) */
    MPI_Bcast(&mbh_data->cpu_time, 1, MPI_DOUBLE, mbh_data->procID_min, MPI_COMM_WORLD);


}

static void share_noise_model(struct NoiseData *noise_data, struct UCBData *ucb_data, struct VGBData *vgb_data, struct MBHData *mbh_data, struct GlobalFitData *global_fit, int root, int procID)
{
    
    int ic = 0;
    if(procID==root)
    {
        ic = noise_data->chain->index[0];
        struct Noise *model_psd = noise_data->inst_model[ic]->psd;
        copy_noise(model_psd,global_fit->psd);
    }
    
    for(int i=0; i<global_fit->psd->Nchannel; i++)
    {
        for(int j=0; j<global_fit->psd->Nchannel; j++)
        {
            MPI_Bcast(global_fit->psd->C[i][j], global_fit->psd->N, MPI_DOUBLE, root, MPI_COMM_WORLD);
            MPI_Bcast(global_fit->psd->invC[i][j], global_fit->psd->N, MPI_DOUBLE, root, MPI_COMM_WORLD);
        }
    }
    MPI_Bcast(global_fit->psd->detC, global_fit->psd->N, MPI_DOUBLE, root, MPI_COMM_WORLD);
}

static void create_residual(struct GlobalFitData *global_fit, int UCB_Flag, int VGB_Flag, int MBH_Flag)
{
    
    memcpy(global_fit->tdi_full->X, global_fit->tdi_store->X, 2*global_fit->tdi_full->N*sizeof(double));
    memcpy(global_fit->tdi_full->Y, global_fit->tdi_store->Y, 2*global_fit->tdi_full->N*sizeof(double));
    memcpy(global_fit->tdi_full->Z, global_fit->tdi_store->Z, 2*global_fit->tdi_full->N*sizeof(double));

    for(int i=0; i<2*global_fit->tdi_full->N; i++)
    {
        if(global_fit->nUCB>0 && !UCB_Flag)
        {
            global_fit->tdi_full->X[i] -= global_fit->tdi_ucb->X[i];
            global_fit->tdi_full->Y[i] -= global_fit->tdi_ucb->Y[i];
            global_fit->tdi_full->Z[i] -= global_fit->tdi_ucb->Z[i];
        }

        if(global_fit->nVGB>0 && !VGB_Flag)
        {
            global_fit->tdi_full->X[i] -= global_fit->tdi_vgb->X[i];
            global_fit->tdi_full->Y[i] -= global_fit->tdi_vgb->Y[i];
            global_fit->tdi_full->Z[i] -= global_fit->tdi_vgb->Z[i];
        }

        /*
         MBH nodes hold all other MBH states in their global_fit structure
         */
        if(global_fit->nMBH>0)
        {
            global_fit->tdi_full->X[i] -= global_fit->tdi_mbh->X[i];
            global_fit->tdi_full->Y[i] -= global_fit->tdi_mbh->Y[i];
            global_fit->tdi_full->Z[i] -= global_fit->tdi_mbh->Z[i];
        }
    }

}

static void print_data_state(struct NoiseData *noise_data, struct UCBData *ucb_data, struct VGBData *vgb_data, struct MBHData *mbh_data, int UCB_Flag, int VGB_Flag, int Noise_Flag, int MBH_Flag)
{
    if(Noise_Flag)
    {
        print_data(noise_data->data, noise_data->data->tdi, noise_data->flags);
        char filename[256];
        sprintf(filename,"%s/data/current_instrument_noise_model.dat",noise_data->flags->runDir);
        generate_instrument_noise_model(noise_data->data,noise_data->orbit,noise_data->inst_model[noise_data->chain->index[0]]);
        print_noise_model(noise_data->inst_model[noise_data->chain->index[0]]->psd, filename);
        
        if(noise_data->flags->confNoise)
        {
            sprintf(filename,"%s/data/current_foreground_noise_model.dat",noise_data->flags->runDir);
            print_noise_model(noise_data->conf_model[noise_data->chain->index[0]]->psd, filename);
        }

    }
    if(VGB_Flag)
    {
        for(int n=0; n<vgb_data->flags->NVB; n++)print_data(vgb_data->data_vec[n], vgb_data->data_vec[n]->tdi, vgb_data->flags);
    }
    if(UCB_Flag)
    {
        print_data(ucb_data->data, ucb_data->data->tdi, ucb_data->flags);
    }
    if(MBH_Flag)
    {
        char tempFileName[MAXSTRINGSIZE];
        sprintf(tempFileName,"%s/data/power_data.dat.temp",mbh_data->flags->runDir);
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
        
        
        sprintf(tempFileName,"%s/data_0.dat.temp",mbh_data->flags->runDir);
        tempFile = fopen(tempFileName,"w");
        for(int i=0; i<mbh_data->het->MN; i++)
        {
            double f = (double)i/mbh_data->data->Tobs;
            fprintf(tempFile,"%lg %lg %lg %lg %lg %lg %lg\n",f, 0.0,0.0,0.0,0.0,0.0,0.0);
            
        }
        for(int i=mbh_data->het->MN; i<mbh_data->het->MM; i++)
        {
            int re = 2*(i-mbh_data->het->MN);
            int im = re+1;
            double f = (double)i/mbh_data->data->Tobs;
            fprintf(tempFile,"%lg %lg %lg %lg %lg %lg %lg\n",f,
                    mbh_data->tdi->X[re],mbh_data->tdi->X[im],
                    mbh_data->tdi->Y[re],mbh_data->tdi->Y[im],
                    mbh_data->tdi->Z[re],mbh_data->tdi->Z[im]);
        }
        for(int i=mbh_data->het->MM; i<mbh_data->data->N; i++)
        {
            double f = (double)i/mbh_data->data->Tobs;
            fprintf(tempFile,"%lg %lg %lg %lg %lg %lg %lg\n",f, 0.0,0.0,0.0,0.0,0.0,0.0);
            
        }
        
        fclose(tempFile);
        
    }
}

static void print_globalfit_state(struct NoiseData *noise_data, struct UCBData *ucb_data, struct VGBData *vgb_data, struct MBHData *mbh_data, int UCB_Flag, int VGB_Flag, int Noise_Flag, int MBH_Flag, FILE *fptr, int counter)
{
    if(Noise_Flag)
        print_noise_state(noise_data, fptr, counter);
    
    if(VGB_Flag)
        print_vgb_state(vgb_data, fptr, counter);
    
    if(UCB_Flag)
        print_ucb_state(ucb_data, fptr, counter);
    
    if(MBH_Flag)
        print_mbh_state(mbh_data, fptr, counter);
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

static void print_usage()
{
    print_glass_usage();
    print_ucb_usage();
    exit(0);
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

    //check command line format
    if(procID==root) 
    {
        fprintf(stdout, "\n================= GLOBAL FIT ================\n");
        print_LISA_ASCII_art(stdout);
    }
    if(argc==1) print_usage();

    /* Allocate structures to hold global model */
    struct GlobalFitData *global_fit = malloc(sizeof(struct GlobalFitData));
    struct NoiseData     *noise_data = malloc(sizeof(struct NoiseData));
    struct UCBData       *ucb_data   = malloc(sizeof(struct UCBData));
    struct VGBData       *vgb_data   = malloc(sizeof(struct VGBData));
    struct MBHData       *mbh_data   = malloc(sizeof(struct MBHData));

    alloc_gf_data(global_fit);

    /* Allocate UCB data structures */
    alloc_ucb_data(ucb_data, procID);
        
    /* Aliases to ucb structures */
    struct Flags *flags = ucb_data->flags;
    struct Orbit *orbit = ucb_data->orbit;
    struct Chain *chain = ucb_data->chain;
    struct Data  *data  = ucb_data->data;

    /* all processes parse command line and set defaults/flags */
    parse_data_args(argc,argv,data,orbit,flags,chain);
    parse_vgb_args(argc, argv, flags);
    parse_mbh_args(argc, argv, mbh_data);
    parse_ucb_args(argc, argv, flags);
    if(procID==0 && flags->help) print_usage();
    
    /* Allocate remaining data structures */
    alloc_noise_data(noise_data, ucb_data, procID, Nproc-1); //noise runs on root process
    alloc_vgb_data(vgb_data, ucb_data, procID); //vbs run on process 1
    alloc_mbh_data(mbh_data, ucb_data, procID); //next are the MBH processes
    
    /* Store size of each model component */
    global_fit->nMBH = mbh_data->NMBH;
    global_fit->nVGB = vgb_data->flags->NVB;
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
    vgb_data->procID_min =  0;
    vgb_data->procID_max = -1;
    if(global_fit->nVGB>0)
    {
        vgb_data->procID_min = pid_counter;
        vgb_data->procID_max = pid_counter;
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
    ucb_data->procID_min = pid_counter;
    ucb_data->procID_max = Nproc-1;

    /* Tell the user how the resources are allocated */
    if(procID==root)
    {
        
        fprintf(stdout,"\n =============== Global Fit Analysis ============== \n");
        fprintf(stdout,"  %i noise  processes (pid %i)\n",1+noise_data->procID_max-noise_data->procID_min,noise_data->procID_min);
        if(global_fit->nVGB>0) fprintf(stdout,"  %i vgb processes (pid %i)\n",1+vgb_data->procID_max-vgb_data->procID_min,vgb_data->procID_min);
        if(global_fit->nMBH>0) fprintf(stdout,"  %i mbh    processes (pid %i-%i)\n",1+mbh_data->procID_max-mbh_data->procID_min,mbh_data->procID_min,mbh_data->procID_max);
        fprintf(stdout,"  %i ucb processes (pid %i-%i)\n",1+ucb_data->procID_max-ucb_data->procID_min,ucb_data->procID_min,ucb_data->procID_max);
        fprintf(stdout," ================================================== \n");
   }

    /* Setup output directories for chain and data structures */
    if(procID==0)
    {
        /* top level run directory */
        mkdir(ucb_data->flags->runDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        
        /* chains directory for global fit joint samples */
        char dirname[MAXSTRINGSIZE];
        sprintf(dirname,"%s/samples",ucb_data->flags->runDir);
        mkdir(dirname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        /* data directory for raw data dump */
        sprintf(dirname,"%s/data",ucb_data->flags->runDir);
        mkdir(dirname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(procID==0)
    {
        /* joint chain file for noise parameters */
        char filename[MAXSTRINGSIZE];
        sprintf(filename,"%s/samples/noise_samples_%04d.dat",noise_data->flags->runDir,procID);
        global_fit->chainFile=fopen(filename,"w");

        /* output directory for noise module */
        sprintf(noise_data->flags->runDir,"%s/noise",noise_data->flags->runDir);
        setup_run_directories(noise_data->flags, noise_data->data, noise_data->chain);
    }
    else if(procID>=vgb_data->procID_min && procID<=vgb_data->procID_max)
    {
        
        /* joint chain file for noise parameters */
        char filename[MAXSTRINGSIZE];
        sprintf(filename,"%s/samples/vgb_samples_%04d.dat",vgb_data->flags->runDir,procID);
        global_fit->chainFile=fopen(filename,"w");

        sprintf(vgb_data->flags->runDir,"%s/vgb",vgb_data->flags->runDir);
        mkdir(vgb_data->flags->runDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        //save original runDir
        char runDir[MAXSTRINGSIZE];
        sprintf(runDir,"%s",vgb_data->flags->runDir);

        for(int n=0; n<vgb_data->flags->NVB; n++)
        {
            //temporarily assign runDir to seg subdir
            sprintf(vgb_data->flags->runDir,"%s/seg_%04d",runDir,n);
            setup_run_directories(vgb_data->flags, vgb_data->data_vec[n], vgb_data->chain_vec[n]);
        }
        
        //restore original runDir
        sprintf(vgb_data->flags->runDir,"%s",runDir);
    }
    else if(procID>=ucb_data->procID_min && procID<=ucb_data->procID_max)
    {
        /* joint chain file for noise parameters */
        char filename[MAXSTRINGSIZE];
        sprintf(filename,"%s/samples/ucb_samples_%04d.dat",ucb_data->flags->runDir,procID);
        global_fit->chainFile=fopen(filename,"w");

        sprintf(ucb_data->flags->runDir,"%s/ucb",ucb_data->flags->runDir);
        mkdir(ucb_data->flags->runDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        sprintf(ucb_data->flags->runDir,"%s/seg_%04d",ucb_data->flags->runDir, procID-ucb_data->procID_min);
        setup_run_directories(ucb_data->flags, ucb_data->data, ucb_data->chain);

    }
    else if(procID>=mbh_data->procID_min && procID <=mbh_data->procID_max)
    {
        
        /* joint chain file for noise parameters */
        char filename[MAXSTRINGSIZE];
        sprintf(filename,"%s/samples/mbh_samples_%04d.dat",mbh_data->flags->runDir,procID);
        global_fit->chainFile=fopen(filename,"w");

        sprintf(mbh_data->flags->runDir,"%s/mbh",mbh_data->flags->runDir);
        mkdir(mbh_data->flags->runDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        sprintf(mbh_data->flags->runDir,"%s/src%04d",mbh_data->flags->runDir, procID-mbh_data->procID_min);
        mkdir(mbh_data->flags->runDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        
    }
       
    if(procID>=ucb_data->procID_min && procID<=ucb_data->procID_max)
    {
        setup_frequency_segment(ucb_data);
    }
    
    /* Initialize ucb data structures */
    alloc_data(data, flags);
            
    /* root process reads data */
    if(procID==root) ReadHDF5(data,global_fit->tdi_full,flags);
    MPI_Barrier(MPI_COMM_WORLD);

    /* alias of full TDI data */
    struct TDI *tdi_full = global_fit->tdi_full;

    /* broadcast data to all ucb processes and select frequency segment */
    share_data(tdi_full, root, procID);
    
    /* composite models are allocated to be the same size as tdi_full */
    setup_gf_data(global_fit);

    /* set up data for ucb model processes */
    setup_ucb_data(ucb_data, tdi_full);
    
    /* set up data for verification binary model processes */
    if(vgb_data->flags->NVB>0)setup_vgb_data(vgb_data, ucb_data, tdi_full);

    /* set up data for mbh model processes */
    if(mbh_data->NMBH>0)setup_mbh_data(mbh_data, ucb_data, tdi_full, procID);

    /* set up data for noise model processes */
    setup_noise_data(noise_data, ucb_data, vgb_data, mbh_data, tdi_full, procID);

    /* allocate global fit noise model for all processes */
    alloc_noise(global_fit->psd, noise_data->data->N, noise_data->data->Nchannel);
    
    /* print ASCII copy of raw data for visualizations */
    if(procID == 0) dump_data(noise_data->data, flags);
    
    /*
     * Initialize all of the samplers
     *
     */
    //choose which sampler to run based on procID
    int UCB_Flag = 0;
    int VGB_Flag = 0;
    int Noise_Flag = 0;
    int MBH_Flag = 0;

    /* Assign processes to Noise model */
    if(procID == 0)
    {
        Noise_Flag = 1;
        initialize_noise_sampler(noise_data);

        int ic = noise_data->chain->index[0];
        struct Noise *model_psd = noise_data->inst_model[ic]->psd;
        invert_noise_covariance_matrix(model_psd);
        copy_noise(model_psd,global_fit->psd);

    }
    MPI_Bcast(global_fit->psd->f, global_fit->psd->N, MPI_DOUBLE, root, MPI_COMM_WORLD);
    for(int i=0; i<global_fit->psd->Nchannel; i++)
    {
        for(int j=0; j<global_fit->psd->Nchannel; j++)
        {
            MPI_Bcast(global_fit->psd->C[i][j], global_fit->psd->N, MPI_DOUBLE, root, MPI_COMM_WORLD);
            MPI_Bcast(global_fit->psd->invC[i][j], global_fit->psd->N, MPI_DOUBLE, root, MPI_COMM_WORLD);
        }
    }
    MPI_Bcast(global_fit->psd->detC, global_fit->psd->N, MPI_DOUBLE, root, MPI_COMM_WORLD);

    share_noise_model(noise_data,ucb_data,vgb_data,mbh_data,global_fit,root,procID);

    /* Assign processes to VB model */
    if(procID>=vgb_data->procID_min && procID<=vgb_data->procID_max)
    {
        VGB_Flag = 1;
        initialize_vgb_sampler(vgb_data);
    }
    
    /* Assign processes to UCB model */
    if(procID >= ucb_data->procID_min && procID <= ucb_data->procID_max)
    {
        UCB_Flag = 1;
        initialize_ucb_sampler(ucb_data);
    }

    /* Assign processes to MBH model */
    if(procID >= mbh_data->procID_min && procID <= mbh_data->procID_max)
    {
        MBH_Flag = 1;
        initialize_mbh_sampler(mbh_data);
    }

    /* distribute current state of models to worker nodes */
    share_noise_model (noise_data, ucb_data, vgb_data, mbh_data, global_fit, root, procID);
    if(global_fit->nUCB>0)share_ucb_model(ucb_data, vgb_data, mbh_data, global_fit, root, procID);
    if(global_fit->nVGB>0)share_vgb_model(ucb_data, vgb_data, mbh_data, global_fit, root, procID);
    if(global_fit->nMBH>0)share_mbh_model   (ucb_data, vgb_data, mbh_data, global_fit, root, procID);

    /* number of update steps for each module (scaled to MBH model update) */
    int cycle;
    global_fit->max_block_time=1;
    /*
     * Master Blocked Gibbs sampler
     *
     */
    int global_fit_counter = 0;
    do
    {
        cycle=1;
        /* ============================= */
        /*     ULTRACOMPACT BINARIES     */
        /* ============================= */

        /* ucb sampler gibbs update */
        if(UCB_Flag)
        {
            create_residual(global_fit, UCB_Flag, VGB_Flag, MBH_Flag);
            
            select_frequency_segment(ucb_data->data, tdi_full);

            select_noise_segment(global_fit->psd, ucb_data->data, ucb_data->chain, ucb_data->model);
            

            cycle = (int)round(global_fit->max_block_time/ucb_data->cpu_time);

            exchange_ucb_source_params(ucb_data);
            if(procID%2==0)
            {
                for(int i=0; i<((cycle > 1 ) ? cycle : 1); i++)
                    ucb_data->status = update_ucb_sampler(ucb_data);
            }
            
            exchange_ucb_source_params(ucb_data);
            if(procID%2!=0)
            {
                for(int i=0; i<((cycle > 1 ) ? cycle : 1); i++)
                    ucb_data->status = update_ucb_sampler(ucb_data);
            }
            
            global_fit->block_time = ucb_data->cpu_time;
        }

        /* ============================= */
        /*  MASSIVE BLACK HOLE BINARIES  */
        /* ============================= */

        /* mbh model update */
        if(MBH_Flag)
        {
            create_residual(global_fit, UCB_Flag, VGB_Flag, MBH_Flag);
            
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

        if(global_fit->nUCB>0)share_ucb_model(ucb_data, vgb_data, mbh_data, global_fit, root, procID);
        if(global_fit->nMBH>0)share_mbh_model   (ucb_data, vgb_data, mbh_data, global_fit, root, procID);

        
        /* ============================= */
        /*   VERIFICATION BINARY MODEL   */
        /* ============================= */

        /* vgb sampler gibbs update */
        if(VGB_Flag)
        {
            create_residual(global_fit, UCB_Flag, VGB_Flag, MBH_Flag);

            select_vgb_segments(vgb_data, tdi_full);

            for(int n=0; n<vgb_data->flags->NVB; n++)
                select_noise_segment(global_fit->psd, vgb_data->data_vec[n], vgb_data->chain_vec[n], vgb_data->model_vec[n]);

            //cycle = (int)round(global_fit->max_block_time/vgb_data->cpu_time);
            //for(int i=0; i<((cycle > 1 ) ? cycle : 1); i++)
            vgb_data->status = update_vgb_sampler(vgb_data);
            
            global_fit->block_time = vgb_data->cpu_time;
        }

        /* ============================= */
        /*    INSTRUMENT NOISE MODEL     */
        /* ============================= */

        /* noise model update */
        if(Noise_Flag)
        {
            create_residual(global_fit, UCB_Flag, VGB_Flag, MBH_Flag);

            select_frequency_segment(noise_data->data, tdi_full);

            //cycle = (int)round(global_fit->max_block_time/noise_data->cpu_time);
            //for(int i=0; i<((cycle > 1 ) ? cycle : 1); i++)
            noise_data->status = update_noise_sampler(noise_data);
            
            global_fit->block_time = noise_data->cpu_time;
        }

        /* get global status of ucb samplers */
        ucb_data->status = get_ucb_status(ucb_data,Nproc,root,procID);

        /* ========================================= */
        /* MPI EXCHANGES OF NOISE & VGB MODEL STATES */
        /* ========================================= */
        
        share_noise_model (noise_data, ucb_data, vgb_data, mbh_data, global_fit, root, procID);
        if(global_fit->nVGB>0)share_vgb_model(ucb_data, vgb_data, mbh_data, global_fit, root, procID);

        /* send time spent in each block to root for load balancing */
        blocked_gibbs_load_balancing(global_fit, root, procID, Nproc);

        /* DEBUG */
        //print_data_state(noise_data,ucb_data,vgb_data,mbh_data,UCB_Flag,VGB_Flag,Noise_Flag,MBH_Flag);

        /* save state of global model */
        print_globalfit_state(noise_data,ucb_data,vgb_data,mbh_data,UCB_Flag,VGB_Flag,Noise_Flag,MBH_Flag, global_fit->chainFile, global_fit_counter);
        
        global_fit_counter++;

    }while(ucb_data->status!=0);
    
    /*
     * Post processing model components
     *
     */
    if(UCB_Flag)
    {
        /* waveform reconstructions */
        print_waveforms_reconstruction(ucb_data->data, ucb_data->flags);
                
        /* evidence results */
        print_evidence(chain,flags);
    }
    if(VGB_Flag)
    {
        for(int n=0; n<vgb_data->flags->NVB; n++)
        {
            /* waveform reconstructions */
            print_waveforms_reconstruction(vgb_data->data_vec[n], vgb_data->flags);
        }
    }
    if(Noise_Flag)
    {
        char filename[128];
        sprintf(filename,"%s/data/final_instrument_noise_model.dat",noise_data->flags->runDir);
        generate_instrument_noise_model(noise_data->data,noise_data->orbit,noise_data->inst_model[noise_data->chain->index[0]]);
        print_noise_model(noise_data->inst_model[noise_data->chain->index[0]]->psd, filename);

        print_noise_reconstruction(noise_data->data, noise_data->flags);
    }
    if(MBH_Flag)
    {
        /* waveform reconstructions */
        //print_mbh_waveform_reconstruction(mbh_data);
    }
    fclose(global_fit->chainFile);
    
    print_data_state(noise_data,ucb_data,vgb_data,mbh_data,UCB_Flag,VGB_Flag,Noise_Flag,MBH_Flag);

    //print total run time
    stop = time(NULL);

    if(procID==root) printf(" ELAPSED TIME = %g seconds on %i processes\n",(double)(stop-start),Nproc);

    MPI_Finalize();//ends the parallelization

    return 0;
}


