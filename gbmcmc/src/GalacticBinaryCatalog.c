/*
 *  Copyright (C) 2019 Tyson B. Littenberg (MSFC-ST12), Kristen Lackeos, Neil J. Cornish
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <getopt.h>


#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort_vector.h>

#include <LISA.h>
#include <GMM_with_EM.h>


#include "GalacticBinary.h"
#include "GalacticBinaryIO.h"
#include "GalacticBinaryMath.h"
#include "GalacticBinaryModel.h"
#include "GalacticBinaryPrior.h"
#include "GalacticBinaryWaveform.h"
#include "GalacticBinaryCatalog.h"


void alloc_entry(struct Entry *entry, int IMAX)
{
    entry->I = 0;
    entry->source = malloc(IMAX*sizeof(struct Source*));
    entry->match  = malloc(IMAX*sizeof(double));
    entry->distance  = malloc(IMAX*sizeof(double));
    entry->stepFlag = calloc(IMAX,sizeof(int));
    entry->gmm = malloc(sizeof(struct GMM));
}

void create_empty_source(struct Catalog *catalog, int NFFT, int Nchannel, int NP)
{
    int N = catalog->N;
    
    //allocate memory for new entry in catalog
    catalog->entry[N] = malloc(sizeof(struct Entry));
    struct Entry *entry = catalog->entry[N];
    
    alloc_entry(entry,1);
    entry->source[entry->I] = malloc(sizeof(struct Source));
    alloc_source(entry->source[entry->I], NFFT, Nchannel, NP);

    entry->match[entry->I] = 1.0;
    entry->distance[entry->I] = 0.0;
    
    entry->I++; //increment number of samples for entry
    catalog->N++;//increment number of entries for catalog
}

void create_new_source(struct Catalog *catalog, struct Source *sample, struct Noise *noise, int i, int IMAX, int NFFT, int Nchannel, int NP)
{
    int N = catalog->N;
    
    //allocate memory for new entry in catalog
    catalog->entry[N] = malloc(sizeof(struct Entry));
    struct Entry *entry = catalog->entry[N];
    
    alloc_entry(entry,IMAX);
    entry->source[entry->I] = malloc(sizeof(struct Source));
    alloc_source(entry->source[entry->I], NFFT, Nchannel, NP);
    
    //add sample to the catalog as the new entry
    copy_source(sample, entry->source[entry->I]);
    
    //store SNR of reference sample to set match criteria
    entry->SNR = snr(sample,noise);
    
    entry->match[entry->I] = 1.0;
    entry->distance[entry->I] = 0.0;
    entry->stepFlag[i] = 1;
    
    entry->I++; //increment number of samples for entry
    catalog->N++;//increment number of entries for catalog
}

void append_sample_to_entry(struct Entry *entry, struct Source *sample, int IMAX, int NFFT, int Nchannel, int NP)
{
    //malloc source structure for this sample's entry
    entry->source[entry->I] = malloc(sizeof(struct Source));
    /* This is taking up too much memory
    alloc_source(entry->source[entry->I], NFFT, Nchannel, NP);
    
    //copy chain sample into catalog entry
    copy_source(sample, entry->source[entry->I]);
    */

    /* leaner way of storing source info */
    
    //only copy source parameters
    entry->source[entry->I]->NP = NP;
    entry->source[entry->I]->params=calloc(NP,sizeof(double));
    memcpy(entry->source[entry->I]->params, sample->params, sample->NP*sizeof(double));
    
    //need Tobs which isn't stored in entries
    double T = sample->params[0]/sample->f0;
    
    //get physical parameters for waveform calculations later
    map_array_to_params(entry->source[entry->I], entry->source[entry->I]->params, T);

    //increment number of stored samples for this entry
    entry->I++;
}

void get_correlation_matrix(struct Data *data, struct Catalog *catalog, int *detection_index, int detections, int IMAX, double **corr)
{
    int NP = data->NP;
    struct Entry *entry = NULL;
    
    /*
     compute mean and variance for each parameter
     */
    double **mean = malloc(detections*sizeof(double *));
    double **var = malloc(detections*sizeof(double *));
    
    for(int d=0; d<detections; d++)
    {
        mean[d] = calloc(NP,sizeof(double));
        var[d] = calloc(NP,sizeof(double));
        
        entry = catalog->entry[detection_index[d]];

        for(int n=0; n<NP; n++)
        {
            double x;
            
            for(int i=0; i<entry->I; i++)
            {
                x = entry->source[i]->params[n];
                mean[d][n] += x;
            }
            mean[d][n] /= (double)entry->I;

            /*
             computing variance after mean (instead of using 1-pass method)
             because rounding error was non-negligible for frequency parameter
             std << mu
             */
            for(int i=0; i<entry->I; i++)
            {
                x = entry->source[i]->params[n];
                var[d][n] += (x - mean[d][n])*(x - mean[d][n]);
            }
            
            var[d][n] /=(double)(entry->I);
        }
    }
    
    /*
     compute correlation matrix
     */
    int N = detections*NP;
    struct Entry *n_entry=NULL;
    struct Entry *m_entry=NULL;
    
    for(int n=0; n<N; n++)
    {
        for(int m=0; m<N; m++)
        {
            //which source row?
            int nd = n/NP;

            //which source column?
            int md = m/NP;

            //which parameter row?
            int nx = n - nd*NP;
            
            //which parameter column?
            int mx = m - md*NP;
                        
            //which entries?
            n_entry = catalog->entry[detection_index[nd]];
            m_entry = catalog->entry[detection_index[md]];
            

            /*
             for each element of the correlation matrix,
             sum over non-null chain samples
             */
            int corr_count=0;
            int ncount=0;
            int mcount=0;
                                     
            for(int i=0; i<IMAX; i++)
            {
                //is this a valid pairing?
                if(n_entry->stepFlag[i] && m_entry->stepFlag[i])
                {
                    double X = n_entry->source[ncount]->params[nx];
                    double Y = m_entry->source[mcount]->params[mx];
                    corr[n][m] += (X - mean[nd][nx]) * (Y - mean[md][mx]);
                    corr_count++;
                    
                }
                
                //advance n-counter
                if(n_entry->stepFlag[i]) ncount++;
                
                //advance m-counter
                if(m_entry->stepFlag[i]) mcount++;
                
            }
            corr[n][m] /= (double)corr_count;
            corr[n][m] /= sqrt(var[nd][nx]*var[md][mx]);;
        }
    }
}

int gaussian_mixture_model_wrapper(double **ranges, struct Flags *flags, struct Entry *entry, char *outdir, size_t NP, size_t NMODE, size_t NTHIN, gsl_rng *seed, double *BIC)
{
    if(flags->verbose)fprintf(stdout,"Event %s, NMODE=%i\n",entry->name,(int)NMODE);
    
    // number of samples
    size_t NMCMC = entry->I;
    
    // number of EM iterations
    size_t NSTEP = 100;
    
    // thin chain
    NMCMC /= NTHIN;
    
    struct Sample **samples = malloc(NMCMC*sizeof(struct Sample*));
    for(size_t n=0; n<NMCMC; n++)
    {
        samples[n] = malloc(sizeof(struct Sample));
        samples[n]->x = gsl_vector_alloc(NP);
        samples[n]->p = gsl_vector_alloc(NMODE);
        samples[n]->w = gsl_vector_alloc(NMODE);
    }
    
    // covariance matrices for different modes
    struct MVG **modes = malloc(NMODE*sizeof(struct MVG*));
    for(size_t n=0; n<NMODE; n++)
    {
        modes[n] = malloc(sizeof(struct MVG));
        alloc_MVG(modes[n],NP);
    }
    
    // Logistic mapping of samples onto R
    double y;
    double pmin,pmax;
    gsl_vector *y_vec = gsl_vector_alloc(NMCMC);
    
    /* parse chain file */
    gsl_vector **params = malloc(NP*sizeof(gsl_vector *));
    for(size_t n=0; n<NP; n++) params[n] = gsl_vector_alloc(NMCMC);
    double value[NP];
    for(size_t i=0; i<NMCMC; i++)
    {
        value[0] = entry->source[i*NTHIN]->f0;
        value[1] = entry->source[i*NTHIN]->costheta;
        value[2] = entry->source[i*NTHIN]->phi;
        value[3] = log(entry->source[i*NTHIN]->amp);
        value[4] = entry->source[i*NTHIN]->cosi;
        value[5] = entry->source[i*NTHIN]->psi;
        value[6] = entry->source[i*NTHIN]->phi0;
        if(NP>7)
            value[7] = entry->source[i*NTHIN]->dfdt;
        if(NP>8)
            value[8] = entry->source[i*NTHIN]->d2fdt2;
        
        for(size_t n=0; n<NP; n++)
        {
            //gsl_vector_set(samples[i]->x,n,value[n]);
            gsl_vector_set(params[n],i,value[n]);
        }
    }
    
    /* Use priors to set min and max of each parameter*/
    for(size_t n=0; n<NP; n++)
    {
        // copy max and min into each MVG structure
        for(size_t k=0; k<NMODE; k++)
        {
            gsl_matrix_set(modes[k]->minmax,n,0,ranges[n][0]);
            gsl_matrix_set(modes[k]->minmax,n,1,ranges[n][1]);
        }
    }
    
    
    /* map params to R with logit function */
    for(size_t n=0; n<NP; n++)
    {
        pmin = gsl_matrix_get(modes[0]->minmax,n,0);
        pmax = gsl_matrix_get(modes[0]->minmax,n,1);
        logit_mapping(params[n], y_vec, pmin, pmax);
        
        for(size_t i=0; i<NMCMC; i++)
        {
            y = gsl_vector_get(y_vec,i);
            gsl_vector_set(samples[i]->x,n,y);
        }
    }
    
    /* The main Gaussian Mixture Model with Expectation Maximization function */
    double logL;
    if(GMM_with_EM(modes,samples,NMCMC,NSTEP,seed,&logL,BIC)) return 1;
    
    /* Write GMM results to binary for pick up by other processes */
    char filename[BUFFER_SIZE];
    sprintf(filename,"%s/%s_gmm.bin",outdir,entry->name);
    FILE *fptr = fopen(filename,"wb");
    fwrite(&NMODE, sizeof NMODE, 1, fptr);
    for(size_t n=0; n<NMODE; n++) write_MVG(modes[n],fptr);
    fclose(fptr);
    
    /* print 1D PDFs and 2D contours of GMM model */
    if(flags->verbose) print_model(modes, samples, NMCMC, logL, *BIC, NMODE);
    
    /* clean up */
    for(size_t n=0; n<NMCMC; n++)
    {
        gsl_vector_free(samples[n]->x);
        gsl_vector_free(samples[n]->p);
        gsl_vector_free(samples[n]->w);
        free(samples[n]);
    }
    free(samples);
    
    // covariance matrices for different modes
    for(size_t n=0; n<NMODE; n++) free_MVG(modes[n]);
    free(modes);
    
    return 0;
}

