
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
    FILE *chain_file1;
    FILE *chain_file2;
    int NMAX = 10;   //max number of frequency & time segments
    int DMAX = 30;   //100; //max number of GB waveforms
    
    /* Allocate data structures */
    struct Flags *flags = malloc(sizeof(struct Flags));
    struct Orbit *orbit = malloc(sizeof(struct Orbit));
    struct Chain *chain = malloc(sizeof(struct Chain));
    struct Data  **data = malloc(sizeof(struct Data*)*NMAX); //data[NF]
    
    
//     Parse command line and set defaults/flags
    for(int i=0; i<NMAX; i++)
    {
        data[i] = malloc(sizeof(struct Data));
        data[i]->t0   = malloc( NMAX * sizeof(double) );
    }
    parse(argc,argv,data,orbit,flags,chain,NMAX,DMAX);
    alloc_data(data, flags);
    struct Data **data_vec = data;
    struct Data *data2  = data_vec[0];

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
    
    
    chain_file1 = fopen(flags->matchInfile1,"r");
    chain_file2 = fopen(flags->matchInfile2,"r");

    if ( chain_file1 == NULL )
    {
        printf("match-in1 is null\n");
    }
    
    if ( chain_file2 == NULL )
    {
        printf("match-in2 is null\n");
    }

//    int i;
    double f0,dfdt,costheta,phi,amp,cosi,phi0,psi;
    double f02,dfdt2,costheta2,phi2,amp2,cosi2,phi02,psi2;
    double junk;

    int N=0;
    while(!feof(chain_file1))
    {
        fscanf(chain_file1,"%lg %lg %lg %lg %lg %lg %lg %lg",&junk,&junk,&junk,&junk,&junk,&junk,&junk,&junk);
        N++;
    }
    
    rewind(chain_file1);
    
    N--;
    
    for(int i=0; i < N; i++)
        fscanf(chain_file1,"%lg %lg %lg %lg %lg %lg %lg %lg",&f0,&dfdt,&amp,&phi,&costheta,&cosi,&psi,&phi0);
    
    for(int i=0; i < N; i++)
        fscanf(chain_file2,"%lg %lg %lg %lg %lg %lg %lg %lg",&f02,&dfdt2,&amp2,&phi2,&costheta2,&cosi2,&psi2,&phi02);
    
    
            struct Source *src1 = NULL;
            src1 = malloc(sizeof(struct Source));
            alloc_source(src1, data2->N,2,data2->NP);
            double *SnA1 = NULL;
            SnA1 = malloc(sizeof(data2->N));
            double *SnE1 = NULL;
            SnE1 = malloc(sizeof(data2->N));
            struct Source *src2 = NULL;
            src2 = malloc(sizeof(struct Source));
            alloc_source(src2, data2->N,2,data2->NP);
            double *SnA2 = NULL;
            SnA2 = malloc(sizeof(data2->N));
            double *SnE2 = NULL;
            SnE2 = malloc(sizeof(data2->N));
////
            for(int n=0; n<2*data2->N; n++)
            {
                src1->tdi->A[n] = 0.0;
                src1->tdi->E[n] = 0.0;
                src1->tdi->X[n] = 0.0;
                src2->tdi->A[n] = 0.0;
                src2->tdi->E[n] = 0.0;
                src2->tdi->X[n] = 0.0;
            }
////
//            //map polarization angle into [0:pi], preserving relation to phi0
            if(psi>M_PI) psi  -= M_PI;
            if(phi0>PI2) phi0 -= PI2;
            if(psi2>M_PI) psi2  -= M_PI;
            if(phi02>PI2) phi02 -= PI2;
            //map parameters to vector
            src1->f0       = f0;
            src1->dfdt     = dfdt;
            src1->costheta = costheta;
            src1->phi      = phi;
            src1->amp      = amp;
            src1->cosi     = cosi;
            src1->phi0     = phi0;
            src1->psi      = psi;
            src2->f0       = f02;
            src2->dfdt     = dfdt2;
            src2->costheta = costheta2;
            src2->phi      = phi2;
            src2->amp      = amp2;
            src2->cosi     = cosi2;
            src2->phi0     = phi02;
            src2->psi      = psi2;
//            if(src1->NP>8)src1->d2fdt2 = 11.0/3.0*dfdt*dfdt/f0;
            //src1->d2fdt2 = fddot;

    
    

            map_params_to_array(src1, src1->params, data2->T);
            map_params_to_array(src2, src2->params, data2->T);
//
//            //Book-keeping of injection time-frequency volume
            galactic_binary_alignment(orbit, data2, src1);
//
            galactic_binary(orbit, data2->format, data2->T, data2->t0[0], src1->params, data2->NP, src1->tdi->X, src1->tdi->A, src1->tdi->E, src1->BW, 2);
    
            galactic_binary_alignment(orbit, data2, src2);
            //
            galactic_binary(orbit, data2->format, data2->T, data2->t0[0], src2->params, data2->NP, src2->tdi->X, src2->tdi->A, src2->tdi->E, src2->BW, 2);

//
//            //Add waveform to data TDI channels
            for(int n=0; n<data2->N; n++)
            {
                int i = n+src1->imin;

                src1->tdi->X[2*i]=0.0;
                src1->tdi->X[2*i+1]=0.0;
                src1->tdi->A[2*i]=0.0;
                src1->tdi->A[2*i+1]=0.0;
                src1->tdi->E[2*i]=0.0;
                src1->tdi->E[2*i+1]=0.0;

                src1->tdi->X[2*i]   += src1->tdi->X[2*n];
                src1->tdi->X[2*i+1] += src1->tdi->X[2*n+1];

                src1->tdi->A[2*i]   += src1->tdi->A[2*n];
                src1->tdi->A[2*i+1] += src1->tdi->A[2*n+1];

                src1->tdi->E[2*i]   += src1->tdi->E[2*n];
                src1->tdi->E[2*i+1] += src1->tdi->E[2*n+1];
                
                src2->tdi->X[2*i]=0.0;
                src2->tdi->X[2*i+1]=0.0;
                src2->tdi->A[2*i]=0.0;
                src2->tdi->A[2*i+1]=0.0;
                src2->tdi->E[2*i]=0.0;
                src2->tdi->E[2*i+1]=0.0;

                src2->tdi->X[2*i]   += src2->tdi->X[2*n];
                src2->tdi->X[2*i+1] += src2->tdi->X[2*n+1];

                src2->tdi->A[2*i]   += src2->tdi->A[2*n];
                src2->tdi->A[2*i+1] += src2->tdi->A[2*n+1];

                src2->tdi->E[2*i]   += src2->tdi->E[2*n];
                src2->tdi->E[2*i+1] += src2->tdi->E[2*n+1];
            }
//
//
//            //Get noise spectrum for data segment
            for(int n=0; n<data2->N; n++)
            {
                double f = data2->fmin + (double)(n)/data2->T;
                if(strcmp(data2->format,"phase")==0)
                {
                    SnA1[n] = AEnoise(orbit->L, orbit->fstar, f);
                    SnE1[n] = AEnoise(orbit->L, orbit->fstar, f);
                    SnA2[n] = AEnoise(orbit->L, orbit->fstar, f);
                    SnE2[n] = AEnoise(orbit->L, orbit->fstar, f);
                }
                else if(strcmp(data2->format,"frequency")==0)
                {
                    SnA1[n] = AEnoise_FF(orbit->L, orbit->fstar, f);
                    SnE1[n] = AEnoise_FF(orbit->L, orbit->fstar, f);
                    SnA2[n] = AEnoise_FF(orbit->L, orbit->fstar, f);
                    SnE2[n] = AEnoise_FF(orbit->L, orbit->fstar, f);
                }
                else
                {
                    fprintf(stderr,"Unsupported data format %s",data2->format);
                    exit(1);
                }
            }
    
//
                double match;
//                //COMPUTE MATCH:
//
                double snr1=0;
                double snrx=0;
                double snr2=0;
    
                snr1 += fourier_nwip(src2->tdi->A,src2->tdi->A,SnA2,src2->tdi->N);
                snr1 += fourier_nwip(src2->tdi->E,src2->tdi->E,SnE2,src2->tdi->N);
                snr2 += fourier_nwip(src1->tdi->A,src1->tdi->A,SnA2,src2->tdi->N);
                snr2 += fourier_nwip(src1->tdi->E,src1->tdi->E,SnE2,src2->tdi->N);

    
                snrx = fourier_nwip(src1->tdi->A,src2->tdi->A,SnA2,src2->tdi->N)/sqrt(snr1*snr2);
                snrx += fourier_nwip(src1->tdi->E,src2->tdi->E,SnE2,src2->tdi->N)/sqrt(snr1*snr2);
                match=snrx;
                fprintf(stdout,"\noverlap=%lg\n",match);
    
    
    
    hallowelt(flags);
    if(flags->orbit)free_orbit(orbit);
    printf("\n");
    fclose(chain_file1);
    fclose(chain_file2);

    
    return 0;
}

void hallowelt(struct Flags *flags)
{
    fprintf(stdout,"Done!\n");
}
