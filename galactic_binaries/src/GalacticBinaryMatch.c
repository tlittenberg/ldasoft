
/***************************  REQUIRED LIBRARIES  ***************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//#include "omp.h"

/*************  PROTOTYPE DECLARATIONS FOR INTERNAL FUNCTIONS  **************/

#include "LISA.h"
#include "Constants.h"
#include "BayesLine.h"
#include "GalacticBinary.h"
#include "GalacticBinaryIO.h"
#include "GalacticBinaryData.h"
#include "GalacticBinaryMath.h"
#include "GalacticBinaryPrior.h"
#include "GalacticBinaryModel.h"
#include "GalacticBinaryProposal.h"
#include "GalacticBinaryWaveform.h"



void hallowelt(struct Flags *flags);
/* ============================  MAIN PROGRAM  ============================ */

int main(int argc, char *argv[])
{
    

//    int ic;
//    FILE *match_file;
    FILE *chain_file;
//    FILE *chain_file2;
    int NMAX = 10;   //max number of frequency & time segments
    int DMAX = 30;   //100; //max number of GB waveforms
    
    /* Allocate data structures */
    struct Flags *flags = malloc(sizeof(struct Flags));
    struct Orbit *orbit = malloc(sizeof(struct Orbit));
    struct Chain *chain = malloc(sizeof(struct Chain));
    struct Data  **data = malloc(sizeof(struct Data*)*NMAX); //data[NF]
    
    
    /* Parse command line and set defaults/flags */
    for(int i=0; i<NMAX; i++)
    {
        data[i] = malloc(sizeof(struct Data));
        data[i]->t0   = malloc( NMAX * sizeof(double) );
        data[i]->tgap = malloc( NMAX * sizeof(double) );
    }
    parse(argc,argv,data,orbit,flags,chain,NMAX,DMAX);
    
    
    /* Load spacecraft ephemerides */
    switch(flags->orbit)
    {
        case 0:
            initialize_analytic_orbit(orbit);
            break;
        case 1:
            initialize_numeric_orbit(orbit);
            break;
        default:
            fprintf(stderr,"unsupported orbit type\n");
            return(1);
            break;
    }
    
    /* Initialize data structures */
    alloc_data(data, flags);
    
    
//    FILE *chain_file = fopen(argv[1],"r");
    chain_file = fopen(flags->matchInfile,"r");
//    chain_file2 = fopen(flags->matchInfile,"r");
//    FILE *match_file = fopen(argv[2],"w");
//    match_file = fopen("/Users/klackeos/Desktop/match_out.dat","w");

    if ( chain_file == NULL )
    {
        printf("infile is null\n");
    }

//    if ( match_file == NULL )
//    {
//        printf("outfile is null\n");
//    }
    int i;
    double *f0,*dfdt,*costheta,*phi,*amp,*cosi,*phi0,*psi,junk;

    int N=0;
    while(!feof(chain_file))
    {
        fscanf(chain_file,"%lg %lg %lg %lg %lg %lg %lg %lg",&junk,&junk,&junk,&junk,&junk,&junk,&junk,&junk);
        N++;
    }
    
    rewind(chain_file);
    
    N--;
    
    f0 = malloc(N * sizeof(double));
    dfdt = malloc(N * sizeof(double));
    costheta = malloc(N * sizeof(double));
    phi = malloc(N * sizeof(double));
    amp = malloc(N * sizeof(double));
    cosi = malloc(N * sizeof(double));
    phi0 = malloc(N * sizeof(double));
    psi = malloc(N * sizeof(double));
    
    for(i=0; i < N; i++)
        fscanf(chain_file,"%lg %lg %lg %lg %lg %lg %lg %lg",&f0[i],&dfdt[i],&amp[i],&phi[i],&costheta[i],&cosi[i],&psi[i],&phi0[i]);
    
    
    
//    double *f0,*dfdt,*costheta,*phi,*amp,*cosi,*phi0,*psi,junk;//read parameters from chain file
//    int n;
    //count sources in file
//    for(n=0; n<N; n++)
//    {
//        fprintf(match_file,"%lg %lg %lg %lg %lg %lg %lg %lg\n",f0,dfdt,amp,phi,costheta,cosi,psi,phi0);
//    }
    
    
    
    
    for(int kk=0; kk<N; kk++)
    {
//    while(!feof(chain_file))
//    {
//      scanf(chain_file,"%lg %lg %lg %lg %lg %lg %lg %lg",&f0[kk],&dfdt[kk],&amp[kk],&phi[kk],&costheta[kk],&cosi[kk],&psi[kk],&phi0[kk]);
        fprintf(stdout,"\n%d %lg %lg %lg %lg %lg %lg %lg %lg\n",kk, f0[kk],dfdt[kk],amp[kk],phi[kk],costheta[kk],cosi[kk],psi[kk],phi0[kk]);
        struct Data **data_vec = data;
        struct Data *data2  = data_vec[0];
        
        const gsl_rng_type *T = gsl_rng_default;
        gsl_rng *r = gsl_rng_alloc(T);
        gsl_rng_env_setup();
        gsl_rng_set (r, data_vec[0]->iseed);

    
    
        
//        for(int jj=0; jj<flags->NT; jj++)
//        {
            int jj=0;
            struct TDI *tdi = data2->tdi[jj];
            
            
            //set bandwidth of data segment centered on injection
                data2->fmin = f0[kk] - (data2->N/2)/data2->T;
                data2->fmax = f0[kk] + (data2->N/2)/data2->T;
                data2->qmin = (int)(data2->fmin*data2->T);
                data2->qmax = data2->qmin+data2->N;
                
                //recompute fmin and fmax so they align with a bin
                data2->fmin = data2->qmin/data2->T;
                data2->fmax = data2->qmax/data2->T;
        
            
//            struct Source *src1 = data2->inj;
            struct Source *src1 = NULL;
            alloc_source(...);
            struct Source *src2 = NULL;
            alloc_source(...);

            for(int n=0; n<2*data2->N; n++)
            {
                src1->tdi->A[n] = 0.0;
                src1->tdi->E[n] = 0.0;
                src1->tdi->X[n] = 0.0;
            }
            
            //map polarization angle into [0:pi], preserving relation to phi0
            if(psi[kk]>M_PI) psi[kk]  -= M_PI;
            if(phi0[kk]>PI2) phi0[kk] -= PI2;
            
            //map parameters to vector
            src1->f0       = f0[kk];
            src1->dfdt     = dfdt[kk];
            src1->costheta = costheta[kk];
            src1->phi      = phi[kk];
            src1->amp      = amp[kk];
            src1->cosi     = cosi[kk];
            src1->phi0     = phi0[kk];
            src1->psi      = psi[kk];
            if(data2->NP>8)
            src1->d2fdt2 = 11.0/3.0*dfdt[kk]*dfdt[kk]/f0[kk];
            //src1->d2fdt2 = fddot;
        
            map_params_to_array(src1, src1->params, data2->T);
        
            //Book-keeping of injection time-frequency volume
            galactic_binary_alignment(orbit, data2, src1);
        
            galactic_binary(orbit, data2->format, data2->T, data2->t0[jj], src1->params, data2->NP, src1->tdi->X, src1->tdi->A, src1->tdi->E, src1->BW, 2);
            
            
            //Add waveform to data TDI channels
            for(int n=0; n<src1->BW; n++)
            {
                int i = n+src1->imin;

                tdi->X[2*i]=0.0;
                tdi->X[2*i+1]=0.0;
                tdi->A[2*i]=0.0;
                tdi->A[2*i+1]=0.0;
                tdi->E[2*i]=0.0;
                tdi->E[2*i+1]=0.0;
                
                tdi->X[2*i]   += src1->tdi->X[2*n];
                tdi->X[2*i+1] += src1->tdi->X[2*n+1];
                
                tdi->A[2*i]   += src1->tdi->A[2*n];
                tdi->A[2*i+1] += src1->tdi->A[2*n+1];
                
                tdi->E[2*i]   += src1->tdi->E[2*n];
                tdi->E[2*i+1] += src1->tdi->E[2*n+1];
            }
            
            
            //Get noise spectrum for data segment
            for(int n=0; n<data2->N; n++)
            {
                double f = data2->fmin + (double)(n)/data2->T;
                if(strcmp(data2->format,"phase")==0)
                {
                    data2->noise[jj]->SnA[n] = AEnoise(orbit->L, orbit->fstar, f);
                    data2->noise[jj]->SnE[n] = AEnoise(orbit->L, orbit->fstar, f);
                }
                else if(strcmp(data2->format,"frequency")==0)
                {
                    data2->noise[jj]->SnA[n] = AEnoise_FF(orbit->L, orbit->fstar, f);
                    data2->noise[jj]->SnE[n] = AEnoise_FF(orbit->L, orbit->fstar, f);
                }
                else
                {
                    fprintf(stderr,"Unsupported data format %s",data2->format);
                    exit(1);
                }
            }
            
//            Get injected SNR
            fprintf(stdout,"   ...injected SNR=%g\n",snr(src1, data2->noise[jj]));
        
//            struct Data  **data_b = malloc(sizeof(struct Data*)*NMAX); //data[NF]
            for(int k=0; k<N; k++)
            {
//                scanf(chain_file2,"%lg %lg %lg %lg %lg %lg %lg %lg",&f0[k],&dfdt[k],&amp[k],&phi[k],&costheta[k],&cosi[k],&psi[k],&phi0[k]);
//                fprintf(stdout,"\n%lg %lg %lg %lg %lg %lg %lg %lg %d\n",f0,dfdt,amp,phi,costheta,cosi,psi,phi0,k);
                fprintf(stdout,"\n%lg %lg %lg %lg %lg %lg %lg %lg\n",f0[k],dfdt[k],amp[k],phi[k],costheta[k],cosi[k],psi[k],phi0[k]);
                
                fprintf(stdout,"\nf01=%lg\n",f0[kk]);
                fprintf(stdout,"\nf02=%lg\n",f0[k]);

//                /* Parse command line and set defaults/flags */
//                for(int i=0; i<NMAX; i++)
//                {
//                    data_b[i] = malloc(sizeof(struct Data));
//                    data_b[i]->t0   = malloc( NMAX * sizeof(double) );
//                    data_b[i]->tgap = malloc( NMAX * sizeof(double) );
//                }
//
                
                struct Data **data_vec2 = data;
                struct Data *data3  = data_vec2[0];
                struct Source *src2 = data3->inj;
                int jj=0;
                struct TDI *tdi2 = data3->tdi[jj];

                const gsl_rng_type *T = gsl_rng_default;
                gsl_rng *r = gsl_rng_alloc(T);
                gsl_rng_env_setup();
                gsl_rng_set (r, data_vec2[0]->iseed);
                

                
                //        for(int jj=0; jj<flags->NT; jj++)
                //        {
                
                
                //set bandwidth of data segment centered on injection
                data3->fmin = f0[k] - (data3->N/2)/data3->T;
                data3->fmax = f0[k] + (data3->N/2)/data3->T;
                data3->qmin = (int)(data3->fmin*data3->T);
                data3->qmax = data3->qmin+data3->N;
                
                
                //recompute fmin and fmax so they align with a bin
                data3->fmin = data3->qmin/data3->T;
                data3->fmax = data3->qmax/data3->T;
                
                
                
                
                for(int n=0; n<2*data3->N; n++)
                {
                    src2->tdi->A[n] = 0.0;
                    src2->tdi->E[n] = 0.0;
                    src2->tdi->X[n] = 0.0;
                }
                
                //map polarization angle into [0:pi], preserving relation to phi0
                if(psi[k]>M_PI) psi[k]  -= M_PI;
                if(phi0[k]>PI2) phi0[k] -= PI2;
                
                //map parameters to vector
                src2->f0       = f0[k];
                src2->dfdt     = dfdt[k];
                src2->costheta = costheta[k];
                src2->phi      = phi[k];
                src2->amp      = amp[k];
                src2->cosi     = cosi[k];
                src2->phi0     = phi0[k];
                src2->psi      = psi[k];
                
                if(data3->NP>8)
                    src2->d2fdt2 = 11.0/3.0*dfdt[k]*dfdt[k]/f0[k];
                //src2->d2fdt2 = fddot;
                
                map_params_to_array(src2, src2->params, data3->T);
                
                //Book-keeping of injection time-frequency volume
                galactic_binary_alignment(orbit, data3, src2);
                
                galactic_binary(orbit, data3->format, data3->T, data3->t0[jj], src2->params, data3->NP, src2->tdi->X, src2->tdi->A, src2->tdi->E, src2->BW, 2);
                
                
                //Add waveform to data TDI channels
                for(int n=0; n<src2->BW; n++)
                {
                    int i = n+src2->imin;
                    
                    tdi2->X[2*i]=0.0;
                    tdi2->X[2*i+1]=0.0;
                    tdi2->A[2*i]=0.0;
                    tdi2->A[2*i+1]=0.0;
                    tdi2->E[2*i]=0.0;
                    tdi2->E[2*i+1]=0.0;
                    
                    tdi2->X[2*i]   += src2->tdi->X[2*n];
                    tdi2->X[2*i+1] += src2->tdi->X[2*n+1];
                    
                    tdi2->A[2*i]   += src2->tdi->A[2*n];
                    tdi2->A[2*i+1] += src2->tdi->A[2*n+1];
                    
                    tdi2->E[2*i]   += src2->tdi->E[2*n];
                    tdi2->E[2*i+1] += src2->tdi->E[2*n+1];
                }
                
                
                //Get noise spectrum for data segment
                for(int n=0; n<data3->N; n++)
                {
                    double f = data3->fmin + (double)(n)/data3->T;
                    if(strcmp(data3->format,"phase")==0)
                    {
                        data3->noise[jj]->SnA[n] = AEnoise(orbit->L, orbit->fstar, f);
                        data3->noise[jj]->SnE[n] = AEnoise(orbit->L, orbit->fstar, f);
                    }
                    else if(strcmp(data3->format,"frequency")==0)
                    {
                        data3->noise[jj]->SnA[n] = AEnoise_FF(orbit->L, orbit->fstar, f);
                        data3->noise[jj]->SnE[n] = AEnoise_FF(orbit->L, orbit->fstar, f);
                    }
                    else
                    {
                        fprintf(stderr,"Unsupported data format %s",data3->format);
                        exit(1);
                    }
                }
                
                
                double match;
                //COMPUTE MATCH:
                
                double snr2=0;
                fprintf(stdout,"\nf01=%lg\n",src1->f0);
                fprintf(stdout,"\nf02=%lg\n",src2->f0);
                snr2 = fourier_nwip(tdi->A,tdi2->A,data2->noise[jj]->SnA,tdi->N)/(snr(src1, data2->noise[jj])*snr(src2, data3->noise[jj]));
                //            fprintf(stdout,"\nsnr2=%lg\n",fourier_nwip(tdi->A,tdi->A,data2->noise[jj]->SnA,tdi->N));
                snr2 += fourier_nwip(tdi->E,tdi2->E,data2->noise[jj]->SnA,tdi->N)/(snr(src1, data2->noise[jj])*snr(src2, data3->noise[jj]));
                //            fprintf(stdout,"\nsnr2=%lg\n",fourier_nwip(tdi->E,tdi->E,data2->noise[jj]->SnA,tdi->N));
                
                match=snr2;
                fprintf(stdout,"\noverlap=%lg\n",match);
//                free_source(src2);
            }
//            free_source(src1);
        
        
        
//        }//end jj loop over segments
    }//end nn loop over sources in file
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    printf("\n");
    fclose(chain_file);
//    fclose(chain_file2);
//    fclose(match_file);
    
    
    
//    hallowelt(flags);
    if(flags->orbit)free_orbit(orbit);
    free(f0);
    free(dfdt);
    free(costheta);
    free(phi);
    free(amp);
    free(cosi);
    free(phi0);
    free(psi);
    
    return 0;
}

void hallowelt(struct Flags *flags)
{
    fprintf(stdout,"Hallo welt!\n");
}
