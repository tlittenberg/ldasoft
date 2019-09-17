
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
    
//    galactic_binary(orbit, data->format, data->T, data->t0[jj], inj->params, 8, inj->tdi->X, inj->tdi->A, inj->tdi->E, inj->BW, 2);

//    int ic;
    FILE *match_file;
    FILE *chain_file;
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
    
//    if(argc!=3)
//    {
//        fprintf(stdout,"Usage: gb_match chain_file match_file\n");
//        return 0;
//    }
    
//    FILE *chain_file = fopen(argv[1],"r");
    chain_file = fopen(flags->matchInfile,"r");
//    FILE *match_file = fopen(argv[2],"w");
    match_file = fopen("/Users/klackeos/Desktop/match_out.dat","w");

    if ( chain_file == NULL )
    {
        printf("infile is null\n");
    }

    if ( match_file == NULL )
    {
        printf("outfile is null\n");
    }
    
    
    double f0,dfdt,costheta,phi,amp,cosi,phi0,psi;//read parameters from chain file
    int n;
    //count sources in file
    int N=0;
    while(!feof(chain_file))
    {
        fscanf(chain_file,"%lg %lg %lg %lg %lg %lg %lg %lg",&f0,&dfdt,&amp,&phi,&costheta,&cosi,&psi,&phi0);
        
        N++;
    }
    
    rewind(chain_file);
    
    N--;
    
    for(n=0; n<N; n++)
    {
        fprintf(match_file,"%lg %lg %lg %lg %lg %lg %lg %lg\n",f0,dfdt,amp,phi,costheta,cosi,psi,phi0);
    }
    
    
    
    
    
    
    
    struct Data **data_vec = data;
    
    struct Data *data2  = data_vec[0];
    
//    sprintf(data2->format,"frequency");
    
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc(T);
    gsl_rng_env_setup();
    gsl_rng_set (r, data_vec[0]->iseed);

    
    for(int nn=0; nn<N; nn++)
    {
//        fscanf(injectionFile,"%lg %lg %lg %lg %lg %lg %lg %lg",&f0,&dfdt,&costheta,&phi,&amp,&cosi,&psi,&phi0);
//        //fscanf(injectionFile,"%lg %lg %lg %lg %lg %lg %lg %lg %lg",&f0,&dfdt,&costheta,&phi,&amp,&cosi,&psi,&phi0,&fddot);
        
        for(int jj=0; jj<flags->NT; jj++)
        {
            
            struct TDI *tdi = data2->tdi[jj];
            
            
            //set bandwidth of data segment centered on injection
            if(nn==0)
            {
                data2->fmin = f0 - (data2->N/2)/data2->T;
                data2->fmax = f0 + (data2->N/2)/data2->T;
                data2->qmin = (int)(data2->fmin*data2->T);
                data2->qmax = data2->qmin+data2->N;
                
                //recompute fmin and fmax so they align with a bin
                data2->fmin = data2->qmin/data2->T;
                data2->fmax = data2->qmax/data2->T;
                
                if(jj==0)fprintf(stdout,"Frequency bins for segment [%i,%i]\n",data2->qmin,data2->qmax);
                fprintf(stdout,"   ...start time: %g\n",data2->t0[jj]);
            }
            
            
            struct Source *inj = data2->inj;
            
            for(int n=0; n<2*data2->N; n++)
            {
                inj->tdi->A[n] = 0.0;
                inj->tdi->E[n] = 0.0;
                inj->tdi->X[n] = 0.0;
            }
            
            //map polarization angle into [0:pi], preserving relation to phi0
            if(psi>M_PI) psi  -= M_PI;
            if(phi0>PI2) phi0 -= PI2;
            
            //map parameters to vector
            inj->f0       = f0;
            inj->dfdt     = dfdt;
            inj->costheta = costheta;
            inj->phi      = phi;
            inj->amp      = amp;
            inj->cosi     = cosi;
            inj->phi0     = phi0;
            inj->psi      = psi;
            if(data2->NP>8)
            inj->d2fdt2 = 11.0/3.0*dfdt*dfdt/f0;
            //inj->d2fdt2 = fddot;
            
            map_params_to_array(inj, inj->params, data2->T);
            
//            //save parameters to file
//            sprintf(filename,"injection_parameters_%i_%i.dat",ii,jj);
//            if(nn==0)paramFile=fopen(filename,"w");
//            else     paramFile=fopen(filename,"a");
//            fprintf(paramFile,"%lg ",data->t0[jj]);
//            print_source_params(data, inj, paramFile);
//            fprintf(paramFile,"\n");
//            fclose(paramFile);
            
            //Book-keeping of injection time-frequency volume
            galactic_binary_alignment(orbit, data2, inj);
            
            printf("   ...bandwidth : %i\n",inj->BW);
            printf("   ...fdot      : %g\n",inj->dfdt*data2->T*data2->T);
            printf("   ...fddot     : %g\n",inj->d2fdt2*data2->T*data2->T*data2->T);
            
            //Simulate gravitational wave signal
            //double t0 = data->t0 + jj*(data->T + data->tgap);
            printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! t0 = %g\n",data2->t0[jj]);
            galactic_binary(orbit, data2->format, data2->T, data2->t0[jj], inj->params, data2->NP, inj->tdi->X, inj->tdi->A, inj->tdi->E, inj->BW, 2);
//            fprintf(stdout,"%s",data2->format);

            //Add waveform to data TDI channels
            for(int n=0; n<inj->BW; n++)
            {
                int i = n+inj->imin;
                
                tdi->X[2*i]   += inj->tdi->X[2*n];
                tdi->X[2*i+1] += inj->tdi->X[2*n+1];
                
                tdi->A[2*i]   += inj->tdi->A[2*n];
                tdi->A[2*i+1] += inj->tdi->A[2*n+1];
                
                tdi->E[2*i]   += inj->tdi->E[2*n];
                tdi->E[2*i+1] += inj->tdi->E[2*n+1];
            }
            
//            sprintf(filename,"data/waveform_injection_%i_%i.dat",ii,jj);
//            fptr=fopen(filename,"w");
//            for(int i=0; i<data->N; i++)
//            {
//                double f = (double)(i+data->qmin)/data->T;
//                fprintf(fptr,"%lg %lg %lg %lg %lg",
//                        f,
//                        tdi->A[2*i],tdi->A[2*i+1],
//                        tdi->E[2*i],tdi->E[2*i+1]);
//                fprintf(fptr,"\n");
//            }
//            fclose(fptr);
            
//            sprintf(filename,"data/power_injection_%i_%i.dat",ii,jj);
//            fptr=fopen(filename,"w");
//            for(int i=0; i<data->N; i++)
//            {
//                double f = (double)(i+data->qmin)/data->T;
//                fprintf(fptr,"%.12g %lg %lg ",
//                        f,
//                        tdi->A[2*i]*tdi->A[2*i]+tdi->A[2*i+1]*tdi->A[2*i+1],
//                        tdi->E[2*i]*tdi->E[2*i]+tdi->E[2*i+1]*tdi->E[2*i+1]);
//                fprintf(fptr,"\n");
//            }
//            fclose(fptr);
            
            //Get noise spectrum for data segment
            for(int n=0; n<data2->N; n++)
            {
                double f = data2->fmin + (double)(n)/data2->T;
                if(strcmp(data2->format,"phase")==0)
                {
                    data2->noise[0]->SnA[n] = AEnoise(orbit->L, orbit->fstar, f);
                    data2->noise[0]->SnE[n] = AEnoise(orbit->L, orbit->fstar, f);
                    if(flags->confNoise)
                    {
                        data2->noise[0]->SnA[n] += GBnoise(data2->T,f);
                        data2->noise[0]->SnE[n] += GBnoise(data2->T,f);
                    }
                }
                else if(strcmp(data2->format,"frequency")==0)
                {
                    data2->noise[0]->SnA[n] = AEnoise_FF(orbit->L, orbit->fstar, f);
                    data2->noise[0]->SnE[n] = AEnoise_FF(orbit->L, orbit->fstar, f);
                    if(flags->confNoise)
                    {
                        data2->noise[0]->SnA[n] += GBnoise_FF(data2->T, orbit->fstar, f);
                        data2->noise[0]->SnE[n] += GBnoise_FF(data2->T, orbit->fstar, f);
                    }
                }
                else
                {
                    fprintf(stderr,"Unsupported data format %s",data2->format);
                    exit(1);
                }
            }
            
            //Get injected SNR
            fprintf(stdout,"   ...injected SNR=%g\n",snr(inj, data2->noise[jj]));
            
            //Add Gaussian noise to injection
            gsl_rng_set (r, data2->nseed+jj);
            
            if(!flags->zeroNoise && nn==0)
            {
                printf("   ...adding Gaussian noise realization\n");
                
                for(int n=0; n<data2->N; n++)
                {
                    tdi->A[2*n]   += gsl_ran_gaussian (r, sqrt(data2->noise[0]->SnA[n])/2.);
                    tdi->A[2*n+1] += gsl_ran_gaussian (r, sqrt(data2->noise[0]->SnA[n])/2.);
                    
                    tdi->E[2*n]   += gsl_ran_gaussian (r, sqrt(data2->noise[0]->SnE[n])/2.);
                    tdi->E[2*n+1] += gsl_ran_gaussian (r, sqrt(data2->noise[0]->SnE[n])/2.);
                }
            }
            
            //Compute fisher information matrix of injection
//            printf("   ...computing Fisher Information Matrix of injection\n");
//
//            galactic_binary_fisher(orbit, data2, inj, data2->noise[jj]);
            
            
//            printf("\n Fisher Matrix:\n");
//            for(int i=0; i<data2->NP; i++)
//            {
//                fprintf(stdout," ");
//                for(int j=0; j<data2->NP; j++)
//                {
//                    if(inj->fisher_matrix[i][j]<0)fprintf(stdout,"%.2e ", inj->fisher_matrix[i][j]);
//                    else                          fprintf(stdout,"+%.2e ",inj->fisher_matrix[i][j]);
//                }
//                fprintf(stdout,"\n");
//            }
            
//            printf("\n Fisher std. errors:\n");
//            for(int j=0; j<data2->NP; j++)  fprintf(stdout," %.4e\n", sqrt(inj->fisher_matrix[j][j]));
            
            
            
//            sprintf(filename,"data/power_data_%i_%i.dat",ii,jj);
//            fptr=fopen(filename,"w");
//
//            for(int i=0; i<data->N; i++)
//            {
//                double f = (double)(i+data->qmin)/data->T;
//                fprintf(fptr,"%.12g %lg %lg ",
//                        f,
//                        tdi->A[2*i]*tdi->A[2*i]+tdi->A[2*i+1]*tdi->A[2*i+1],
//                        tdi->E[2*i]*tdi->E[2*i]+tdi->E[2*i+1]*tdi->E[2*i+1]);
//                fprintf(fptr,"\n");
//            }
//            fclose(fptr);
            
//            sprintf(filename,"data/data_%i_%i.dat",ii,jj);
//            fptr=fopen(filename,"w");
//
//            for(int i=0; i<data->N; i++)
//            {
//                double f = (double)(i+data->qmin)/data->T;
//                fprintf(fptr,"%.12g %lg %lg %lg %lg",
//                        f,
//                        tdi->A[2*i],tdi->A[2*i+1],
//                        tdi->E[2*i],tdi->E[2*i+1]);
//                fprintf(fptr,"\n");
//            }
//            fclose(fptr);
            
            //TODO: fill X vectors with A channel for now
            for(int n=0; n<data2->N; n++)
            {
                tdi->X[2*n]   = tdi->A[2*n];
                tdi->X[2*n+1] = tdi->A[2*n+1];
            }
            
        }//end jj loop over segments
    }//end nn loop over sources in file

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    printf("\n");
    fclose(chain_file);
    fclose(match_file);
    
    
    
    hallowelt(flags);
    
    if(flags->orbit)free_orbit(orbit);
    
    return 0;
}

void hallowelt(struct Flags *flags)
{
    fprintf(stdout,"Hallo welt!\n");
}
