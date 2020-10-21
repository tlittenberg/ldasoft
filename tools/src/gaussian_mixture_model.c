/*
 *  Author: Tyson B. Littenberg (MSFC-ST12)
 *  Created: 07.27.2020
 *
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
 *
 * Gaussian Mixture Model using the Expectation-Maximization Algorithm
 *
 * Compile: gcc gaussian_mixture_model.c -lm -lgsl -lgslcblas -o gaussian_mixture_model
 * Usage: gaussian_mixture_model -h
 * Example: Fit 4 gaussians to the first two columns of `chain.dat`
 *          gaussian_mixture_model --file chain.dat --modes 4 --nparams 2
 *
 */

/**
@file gaussian_mixture_model.c
\brief Gaussian Mixture Model fit to input data samples using Expectation-Maximization Algorithm.
 
 \todo More generalization is possible.
*/


/*
 REQUIRED LIBRARIES
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <getopt.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort_vector.h>

#include "GMM_with_EM.h"




/* map parameters from R to [xmin,xmax] with sigmoid function
for(size_t n=0; n<NP; n++)
{
    logit(params[n], y_vec, pxminmax);

    for(size_t i=0; i<NMCMC; i++)
    {
        y = gsl_vector_get(y_vec,i);
        gsl_vector_set(samples[i]->x,n,y);
    }
}*/

/**
 * \brief Main function for data handling and calling GM_with_EM() algorithm
 *
 */
int main(int argc, char* argv[])
{
    // chain file
    FILE *chainFile=NULL;
    
    // dimension of model
    size_t NP = 0;
    
    // maximum number of modes
    size_t NMODE = 2;
    
    // number of samples
    size_t NMCMC = 0;
        
    // number of EM iterations
    size_t NSTEP = 500;

    // thinning rate of input chain
    size_t NTHIN = 1;
    
    // array for flagging parameters to log
    size_t *LFLAG = calloc(100,sizeof(size_t));
    
    // rng
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc(T);
    gsl_rng_env_setup();
    gsl_rng_set (r, 190521);
    
    
    /* parse command line */
    struct option cmd_opt[] =
    {
        { "help",     no_arg,  0, 'h' },
        { "file",     req_arg, 0, 'f' },
        { "log",      req_arg, 0, 'l' },
        { "modes",    req_arg, 0, 'm' },
        { "nparams",  req_arg, 0, 'n' },
        { "seed",     req_arg, 0, 's' },
        { "thin",     req_arg, 0, 't' },
        { 0, 0, 0, 0 }
    };
    
    char args[] = "hl:f:m:n:s:t:";
    char *program = argv[0];
    
    if(argc==1)
    {
        printUsage(program);
        exit(0);
    }
    while(1)
    {
        int opt_indx = 0;
        int c;
        
        c = getopt_long(argc, argv, args, cmd_opt, &opt_indx );
        if(c==-1) // end of options
            break;
        
        switch(c)
        {
            case 'f':
                chainFile = fopen(optarg,"r");
                break;
            case 'h': // help
                printUsage(program);
                exit(0);
            case 'l':
                LFLAG[(size_t)atoi(optarg)] = 1;
                break;
            case 'm':
                NMODE = (size_t)atoi(optarg);
                break;
            case 'n':
                NP = (size_t)atoi(optarg);
                break;
            case 's':
                gsl_rng_set(r,(long)atoi(optarg));
                break;
            case 't':
                NTHIN = (size_t)atoi(optarg);
                break;
            default:
                fprintf(stderr,"unknown error while parsing options\n" );
                exit(1);
        }
    }
    
    /* count lines in file */
    char* line;
    char lineBuffer[BUFFER_SIZE];
    while((line = fgets(lineBuffer, BUFFER_SIZE, chainFile)) != NULL) NMCMC++;
    rewind(chainFile);
    
    //thin chain
    NMCMC /= NTHIN;
    
    
    /* allocate memory in data structures*/
    
    // raw chain samples
    gsl_vector **params = malloc(NP*sizeof(gsl_vector *));
    for(size_t n=0; n<NP; n++) params[n] = gsl_vector_alloc(NMCMC);
    
    // mapped chain samples
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
    double x,y,p;
    double pmin,pmax;
    gsl_vector **pminmax = malloc(NP*sizeof(gsl_vector *));
    for(size_t i=0; i<NP; i++) pminmax[i] = gsl_vector_alloc(2);
    gsl_vector *x_vec = gsl_vector_alloc(NMCMC);
    gsl_vector *y_vec = gsl_vector_alloc(NMCMC);

    /* parse chain file */
    double value;
    char *column;
    for(size_t i=0; i<NMCMC; i++)
    {
        for(size_t j=0; j<NTHIN; j++) line = fgets(lineBuffer, BUFFER_SIZE, chainFile);
        
        column=strtok(line," ");
        
        for(size_t n=0; n<NP; n++)
        {
            sscanf(column, "%lg", &value);
            if(LFLAG[n]) value = log(value);
            //gsl_vector_set(samples[i]->x,n,value);
            gsl_vector_set(params[n],i,value);
            column=strtok(NULL," ");
        }
    }
    
    
    
    
    /* Get max and min for each parameter */
    for(size_t n=0; n<NP; n++)
    {
        int err = 0;
        double *temp = malloc(10*sizeof(double));
        err = gsl_sort_vector_smallest(temp, 10, params[n]);
        pmin = temp[0] - (temp[9]-temp[0]); //pad min to avoid infs in mapping
        
        err = gsl_sort_vector_largest(temp, 10, params[n]);
        pmax = temp[0] + (temp[0]-temp[9]); //pad max to avoid infs in mapping
        
        /* cpopy max and min into each MVG structure */
        for(size_t k=0; k<NMODE; k++)
        {
            gsl_matrix_set(modes[k]->minmax,n,0,pmin);
            gsl_matrix_set(modes[k]->minmax,n,1,pmax);
        }
    }
    
    
    /* map params to R with logit function */
    for(size_t n=0; n<NP; n++)
    {
        double pmin = gsl_matrix_get(modes[0]->minmax,n,0);
        double pmax = gsl_matrix_get(modes[0]->minmax,n,1);
        logit_mapping(params[n], y_vec, pmin, pmax);

        for(size_t i=0; i<NMCMC; i++)
        {
            y = gsl_vector_get(y_vec,i);
            gsl_vector_set(samples[i]->x,n,y);
        }
    }


    /* The main Gaussian Mixture Model with Expectation Maximization function */
    double logL, BIC;
    if(GMM_with_EM(modes,samples,NMCMC,NSTEP,r,&logL,&BIC)) return 1;
    

    /* Write GMM results to binary for pick up by other processes */
    char filename[BUFFER_SIZE];
    sprintf(filename,"gmm_%i.bin",(int)NMODE);
    FILE *fptr = fopen(filename,"wb");
    for(size_t n=0; n<NMODE; n++) write_MVG(modes[n],fptr);
    fclose(fptr);

    /* Read GMM results to binary for pick up by other processes */
    sprintf(filename,"gmm_%i.bin",(int)NMODE);
    fptr = fopen(filename,"rb");
    for(size_t n=0; n<NMODE; n++) read_MVG(modes[n],fptr);
    fclose(fptr);

    
    print_model(modes, samples, NMCMC, logL, BIC, NMODE);

    
    for(size_t n=0; n<NMODE; n++) fprintf(stdout,"%i %g\n",(int)n, modes[n]->p);

    
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

    free(LFLAG);

    
    return 0;
}



