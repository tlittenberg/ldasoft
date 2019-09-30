
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
#include "GalacticBinary.h"
#include "GalacticBinaryIO.h"
#include "GalacticBinaryMath.h"
#include "GalacticBinaryModel.h"
#include "GalacticBinaryWaveform.h"

// CODE USAGE:
// ~/gb_match --match-in1 ~/input1.dat --match-in2 ~/input2.dat --frac-freq --fmin 0.001249 --samples 512 --duration 62914560

void print_match(struct Flags *flags, double match);
/* ============================  MAIN PROGRAM  ============================ */

int main(int argc, char *argv[])
{

    FILE *chain_file1;
    FILE *chain_file2;
    int NMAX = 10;   //max number of frequency & time segments

    
    /* Allocate data structures */
    struct Flags *flags = malloc(sizeof(struct Flags));
    struct Orbit *orbit = malloc(sizeof(struct Orbit));
    struct Chain *chain = malloc(sizeof(struct Chain));
    struct Data  **data = malloc(sizeof(struct Data*)*NMAX); //data[NF]
    
    
//   Parse command line and set defaults/flags
    for(int i=0; i<NMAX; i++)
    {
        data[i] = malloc(sizeof(struct Data));
        data[i]->t0   = malloc( NMAX * sizeof(double) );
    }
    parse(argc,argv,data,orbit,flags,chain,NMAX);
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
    
    // Read data from input files
    for(int i=0; i < N; i++)
        fscanf(chain_file1,"%lg %lg %lg %lg %lg %lg %lg %lg",&f0,&dfdt,&amp,&phi,&costheta,&cosi,&psi,&phi0);
    
    for(int i=0; i < N; i++)
        fscanf(chain_file2,"%lg %lg %lg %lg %lg %lg %lg %lg",&f02,&dfdt2,&amp2,&phi2,&costheta2,&cosi2,&psi2,&phi02);
    
            //allocate memory for two sources and noise
            struct Source *src1 = NULL;
            src1 = malloc(sizeof(struct Source));
            alloc_source(src1, data2->N,2,data2->NP);
            double *SnA1 = NULL;
            SnA1 = malloc(data2->N*sizeof(double));
            double *SnE1 = NULL;
            SnE1 = malloc(data2->N*sizeof(double));
            struct Source *src2 = NULL;
            src2 = malloc(sizeof(struct Source));
            alloc_source(src2, data2->N,2,data2->NP);

            for(int n=0; n<2*data2->N; n++)
            {
                src1->tdi->A[n] = 0.0;
                src1->tdi->E[n] = 0.0;
                src1->tdi->X[n] = 0.0;
                src2->tdi->A[n] = 0.0;
                src2->tdi->E[n] = 0.0;
                src2->tdi->X[n] = 0.0;
            }

            //map polarization angle into [0:pi], preserving relation to phi0
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

            //Book-keeping of injection time-frequency volume
            galactic_binary_alignment(orbit, data2, src1);

            galactic_binary(orbit, data2->format, data2->T, data2->t0[0], src1->params, data2->NP, src1->tdi->X, src1->tdi->A, src1->tdi->E, src1->BW, 2);
    
            galactic_binary_alignment(orbit, data2, src2);
            
            galactic_binary(orbit, data2->format, data2->T, data2->t0[0], src2->params, data2->NP, src2->tdi->X, src2->tdi->A, src2->tdi->E, src2->BW, 2);




            //Get noise spectrum for data segment
            for(int n=0; n<data2->N; n++)
            {
                double f = data2->fmin + (double)(n)/data2->T;
                if(strcmp(data2->format,"phase")==0)
                {
                    SnA1[n] = AEnoise(orbit->L, orbit->fstar, f);
                    SnE1[n] = AEnoise(orbit->L, orbit->fstar, f);
                }
                else if(strcmp(data2->format,"frequency")==0)
                {
                    SnA1[n] = AEnoise_FF(orbit->L, orbit->fstar, f);
                    SnE1[n] = AEnoise_FF(orbit->L, orbit->fstar, f);
                }
                else
                {
                    fprintf(stderr,"Unsupported data format %s",data2->format);
                    exit(1);
                }
            }
    //SHOULD BE A FUNCTION: INPUT SNA and SRC1 AND SRC2, return match double. lives in gb_waveform or *math*... `waveform_match'
    //COMPUTE OVERLAP:
                double match;
                double snr1=0;
                double snr2=0;
                double snrx=0;
    
                snr1 += fourier_nwip(src1->tdi->A,src1->tdi->A,SnA1,src1->tdi->N);
                snr1 += fourier_nwip(src1->tdi->E,src1->tdi->E,SnE1,src1->tdi->N);
                snr2 += fourier_nwip(src2->tdi->A,src2->tdi->A,SnA1,src2->tdi->N);
                snr2 += fourier_nwip(src2->tdi->E,src2->tdi->E,SnE1,src2->tdi->N);

                snrx = fourier_nwip(src1->tdi->A,src2->tdi->A,SnA1,src2->tdi->N)/sqrt(snr1*snr2);
                snrx += fourier_nwip(src1->tdi->E,src2->tdi->E,SnE1,src2->tdi->N)/sqrt(snr1*snr2);
                match=snrx;
    
    //Print results of overlap calc.
    print_match(flags,match);
    if(flags->orbit)free_orbit(orbit);
    printf("\n");
    fclose(chain_file1);
    fclose(chain_file2);

    return 0;
}

void print_match(struct Flags *flags, double match)
{
    if(match > 0.999)
    {
        fprintf(stdout,"These are the same source.\n");
        fprintf(stdout,"overlap = %lg > 0.999\n",match);
    }
    else
    {
        fprintf(stdout,"Different sources.\n");
        fprintf(stdout,"overlap = %lg < 0.999\n",match);
    }
}
