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

/** @name Flags for command line parser */
///@{
#define no_arg 0   //!< no argument
#define req_arg 1  //!< required argument
#define opt_arg 2  //!< optional argument
///@}

/** @name Boundaries for variables */
///@{
#define BUFFER_SIZE 1024 //!< max size of `char` buffers
#define PMIN 1.e-16 //!< floor on probabilitiy densities (avoid singularities at \f$p=0\f$
///@}

/** @name Progress bar settings */
///@{
#define PBSTR "||||||||||||||||||||||||||||||||||||||||" //!< displayed progress
#define PBWIDTH 40 //!<size of progress bar
///@}



/**
 *\brief print progress bar
 */
static void printProgress (double percentage)
{
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush (stdout);
}

/**
 * \brief Show usage
 */
void printUsage(const char *program)
{
    fprintf( stdout,"\n");
    fprintf( stdout, "Gaussian Mixture Model (GMM):\n");
    fprintf( stdout, "  Expectation-Maximization algorithm to fit GMM to \n");
    fprintf( stdout, "  input list of data samples e.g. from MCMC\n");
    fprintf( stdout, "\n");
    fprintf( stdout, "Usage: %s required [optional] \n", program );
    fprintf( stdout, "\n");
    fprintf( stdout, "  Required:\n");
    fprintf( stdout, "     -f, --file=FILE        filename for input chain\n" );
    fprintf( stdout, "     -n, --nparams=INT      number of model parameters\n" );
    fprintf( stdout, "                            extra columns are ignored\n");
    fprintf( stdout, "\n");
    fprintf( stdout, "  Optional:\n");
    fprintf( stdout, "    [-h, --help]            print this message and exit\n" );
    fprintf( stdout, "    [-l, --log=INT]         use log(params) in specified column number [0,m-1]\n" );
    fprintf( stdout, "                            multiple arguments can be given.\n");
    fprintf( stdout, "    [-m, --modes=INT]       number of modes in fit (2)\n");
    fprintf( stdout, "    [-s, --seed=LONG]       RNG seed\n");
    fprintf( stdout, "    [-t, --thin=INT]        downsample rate for chains (1)\n");
    fprintf( stdout,"\n");
}

/**
 * \brief Data structure for individual samples,
 * along with their relative weights and probabilites in each mode.
 */
struct Sample
{
    gsl_vector *x; //!< location in parameter space
    gsl_vector *p; //!< p(x) for each mode
    gsl_vector *w; //!< weight sample for mode
};

/**
 * \brief Data structure for each mode,
 * including parameters and covariance matrix products used in calculations
 */
struct Mode
{
    gsl_vector *mu; //!< means
    gsl_matrix *C; //!< covariance matrix
    gsl_matrix *Cinv; //!< inverse covariance matrix
    gsl_matrix *evectors; //!< eigenvectors
    gsl_vector *evalues; //!< eigenvalues
    double detC; //!< determinant of covariance matrix
    double p; //!< prior for Mode (i.e. weighting)
    double Neff; //!< effective number of samples in mode
};


/**
 * \brief Evaluates the probability density of a multviariate Gaussian with input mean \f$\mu\f$
 * and covariance matrix \f$ C \f$.
 * \param[in] x vector of parameters
 * \param[in] mu vector of means
 * \param[in] Cinv inverse covariance matrix
 * \param[out] probability \f$ \frac{\exp -\frac{1}{2}\left(x - \mu\right)^T C^{-1} \left(x - \mu\right)}{\sqrt{(2\pi)^N \det C}} \f$
 */
double multivariate_gaussian(gsl_vector *x, gsl_vector *mu, gsl_matrix *Cinv, double detC)
{
    size_t N = x->size;
    gsl_vector *dx = gsl_vector_alloc(N);
    
    // x-mu
    for(size_t n=0; n<N; n++)
    {
        double xi = gsl_vector_get(x,n);
        double mi = gsl_vector_get(mu,n);
        double d  = xi-mi;
        
        gsl_vector_set(dx,n,d);
    }
    
    
    /*
     (x-mu)^T C^-1 (x-mu)
     */
    double CdotX;
    double chi2 = 0.0;
    for(size_t m=0; m<N; m++)
    {
        CdotX = 0.0;
        for(size_t n=0; n<N; n++)
        {
            CdotX += gsl_matrix_get(Cinv,m,n)*gsl_vector_get(dx,n);
        }
        chi2 += CdotX*gsl_vector_get(dx,m);
    }
    gsl_vector_free(dx);
    
    double p = exp(-0.5*chi2)/sqrt(pow(2.0*M_PI,N)*detC);
    
    return p > PMIN ? p : PMIN;
}

/**
 * \brief Wrapper for computing matrix inversion using
 * `GSL` functions.
 * \param[in] A input symmetric matrix
 * \param[out] Ainv inverse of A
 * \param[out] detA determinant of A
 * \param[out] R recipricol condition number [0,1]
 */
void invert_matrix(gsl_matrix *A, gsl_matrix *Ainv, double *detA, double *R)
{
    gsl_set_error_handler_off();
    
    //error catchers
    int err = 0;
    
    //get size of matrix (assumed to be NxN)
    size_t N = A->size1;
    
    //some workspace
    gsl_permutation * permutation = gsl_permutation_alloc(N);
    
    //copy A into Ainv because LU decomposition destroys the matrix
    gsl_matrix_memcpy(Ainv,A);
    
    //cholesky decomposition
    int i;
    err += gsl_linalg_cholesky_decomp(Ainv);
    
    //get condition number
    gsl_vector *work = gsl_vector_alloc(3*N);
    err += gsl_linalg_cholesky_rcond(Ainv, R, work);

    //inverse of A
    err += gsl_linalg_cholesky_invert(Ainv);
    
    //get deteriminant, need LU decomposition
    gsl_matrix *L = gsl_matrix_calloc(N,N);
    gsl_matrix_memcpy(L,A);
    gsl_linalg_LU_decomp(L,permutation,&i);
    *detA = gsl_linalg_LU_det(L,i);
    
    //clean up
    gsl_matrix_free(L);
    gsl_vector_free (work);
    gsl_permutation_free (permutation);
}

/**
 * \brief Wrapper for computing matrix eigenvalue decomposition using
 * `GSL` functions.
 * \param[in] A input symmetric matrix
 * \param[out] evec matrix where each row has components of eigenvector
 * \param[out] eval vector where each element is the eigenvalue associated with the same row of `evec`.
 */
void decompose_matrix(gsl_matrix *A, gsl_matrix *evec, gsl_vector *eval)
{
    //error catchers
    int err = 0;
    
    //get size of matrix (assumed to be NxN)
    size_t N = A->size1;
    
    //get deteriminant, need LU decomposition
    gsl_matrix *Atemp = gsl_matrix_calloc(N,N);
    gsl_eigen_symmv_workspace * workspace = gsl_eigen_symmv_alloc (N);
    gsl_permutation * permutation = gsl_permutation_alloc(N);
    
    //copy A into Atemp because eigen_symmv destroys the matrix
    gsl_matrix_memcpy(Atemp,A);
    
    //the reason we're all here...
    err += gsl_eigen_symmv (Atemp, eval, evec, workspace);
    
    gsl_matrix_free (Atemp);
    gsl_eigen_symmv_free (workspace);
    gsl_permutation_free (permutation);
    
}

/**
 * \brief Log-likelihood of Gaussian Mixture Model
 * \param[in] modes parameters of each Gaussian including relative weights \f$\alpha_{k}\f$
 * \param[in] samples data points \f$x_i\f$ and probabilities \f$p(x_i|k)\f$ for each Gaussian \f$k\f$
 * \param[in] NMCMC number of samples
 * \param[in] NMODE number of modes
 * \param[out] log-likelihood \f$\log L = \sum_k^{\rm NMCMC} \log \sum_i^{\rm NMODE} \alpha_k p(x_i | k)\f$
 *
 */
double log_likelihood(struct Mode **modes, struct Sample **samples, int NMCMC, int NMODE)
{
    
    double logL = 0.0;
    for(size_t i=0; i<NMCMC; i++)
    {
        double P = 0.0;
        for(size_t k=0; k<NMODE; k++)
        {
            P += modes[k]->p*gsl_vector_get(samples[i]->p,k);
        }
        if(P==0) exit(1);
        logL += log(P);
    }
    return logL;
}

/**
 * \brief Print joint 1D distributions for each parameter to file
 */
void print_1D_pdfs(struct Mode **modes, struct Sample **samples, size_t NMCMC, char root[], size_t ix)
{
    char filename[128];
    sprintf(filename,"%s_%i.dat",root,(int)ix);
    FILE *fptr = fopen(filename,"w");
    
    size_t NMODE = samples[0]->p->size;
    
    double *xvec = malloc(NMCMC*sizeof(double));
    double xmin,xmax;
    double x0 =  1e60;
    double xf = -1e60;
    for(size_t k=0; k<NMODE; k++)
    {
        for(size_t i=0; i<NMCMC; i++) xvec[i] = gsl_vector_get(samples[i]->x,ix);
        gsl_stats_minmax(&xmin,&xmax,xvec,1,NMCMC);
        if(xmin<x0) x0 = xmin;
        if(xmax>xf) xf = xmax;
    }
    
    double p;
    double x;
    double dx = (xf-x0)/100.;
    for(int n=0; n<=100; n++)
    {
        p = 0.0;
        x = xmin + (double)n*dx;
        for(size_t k=0; k<NMODE; k++)
        {
            double mean = gsl_vector_get(modes[k]->mu,ix);
            double var  = gsl_matrix_get(modes[k]->C,ix,ix);
            p += modes[k]->p*exp( -0.5*(x-mean)*(x-mean)/var )/sqrt(2*M_PI*var);
        }
        fprintf(fptr,"%.16g %.16g\n",x,p);
    }
    
    fclose(fptr);
    
    
    free(xvec);
}

/**
 * \brief Print 1,2, and 3\f$\sigma\f$ contours of each individual
 *  Gaussian for in the model for each parameter pair
 */
void print_2D_contours(struct Mode **modes, size_t NMODE, char root[], size_t x1, size_t x2)
{
    char filename[128];
    FILE *fptr = NULL;
    
    struct Mode **submodes = malloc(NMODE*sizeof(struct Mode*));
    for(size_t k=0; k<NMODE; k++)
    {
        submodes[k] = malloc(sizeof(struct Mode));
        submodes[k]->mu = gsl_vector_alloc(2);
        submodes[k]->C = gsl_matrix_calloc(2,2);
        submodes[k]->Cinv = gsl_matrix_calloc(2,2);
        submodes[k]->evectors = gsl_matrix_calloc(2,2);
        submodes[k]->evalues = gsl_vector_calloc(2);
    }
    
    //Pick parameters
    size_t X[2] = {x1,x2};
    
    for(size_t k=0; k<NMODE; k++)
    {
        for(size_t n=0; n<2; n++)
        {
            gsl_vector_set(submodes[k]->mu,n,gsl_vector_get(modes[k]->mu,X[n]));
            for(size_t m=0; m<2; m++) gsl_matrix_set(submodes[k]->C,m,n,gsl_matrix_get(modes[k]->C,X[m],X[n]));
        }
        
        decompose_matrix(submodes[k]->C, submodes[k]->evectors, submodes[k]->evalues);
    }
    
    
    double x,y;
    for(size_t k=0; k<NMODE; k++)
    {
        sprintf(filename,"%s_%i_%i_%i.dat",root,(int)x1,(int)x2,(int)k);
        fptr=fopen(filename,"w");
        for(int n=0; n<=100; n++)
        {
            double theta = atan2(gsl_matrix_get(submodes[k]->evectors,0,1),gsl_matrix_get(submodes[k]->evectors,0,0));
            double Rx = sqrt(gsl_vector_get(submodes[k]->evalues,0));
            double Ry = sqrt(gsl_vector_get(submodes[k]->evalues,1));
            double Cx = gsl_vector_get(submodes[k]->mu,0);
            double Cy = gsl_vector_get(submodes[k]->mu,1);
            
            double angle = n*(2.*M_PI/100.);
            x = 1.*( Rx*cos(angle)*cos(theta) + Ry*sin(angle)*sin(theta) ) + Cx;
            y = 1.*(-Rx*cos(angle)*sin(theta) + Ry*sin(angle)*cos(theta) ) + Cy;
            fprintf(fptr,"%.16lg %.16lg ",x,y);
            
            x = 2.*( Rx*cos(angle)*cos(theta) + Ry*sin(angle)*sin(theta) ) + Cx;
            y = 2.*(-Rx*cos(angle)*sin(theta) + Ry*sin(angle)*cos(theta) ) + Cy;
            fprintf(fptr,"%.16lg %.16lg ",x,y);
            
            x = 3.*( Rx*cos(angle)*cos(theta) + Ry*sin(angle)*sin(theta) ) + Cx;
            y = 3.*(-Rx*cos(angle)*sin(theta) + Ry*sin(angle)*cos(theta) ) + Cy;
            fprintf(fptr,"%.16lg %.16lg ",x,y);
            
            fprintf(fptr,"\n");
        }
        fclose(fptr);
    }
    
    for(size_t k=0; k<NMODE; k++)
    {
        gsl_vector_free(submodes[k]->mu);
        gsl_matrix_free(submodes[k]->C);
        gsl_matrix_free(submodes[k]->Cinv);
        gsl_matrix_free(submodes[k]->evectors);
        gsl_vector_free(submodes[k]->evalues);
    }
    free(submodes);
}

/**
 * \brief Wrapper for print_1D_pdfs() and print_2D_contours()
 */
void print_model(struct Mode **modes, struct Sample **samples, size_t NMCMC, double logL, double BIC, size_t step)
{
    size_t NP = samples[0]->x->size;
    size_t NMODE = samples[0]->w->size;
    char filename[128];
    
    sprintf(filename,"gmm.dat");
    FILE *fptr = fopen(filename,"a");
    fprintf(fptr,"%i %g %g\n",(int)step, logL, BIC);
    fclose(fptr);
    
    for(size_t m=0; m<NP; m++)
    {
        for(size_t n=m; n<NP; n++)
        {
            /* Get 1D PDFs for plotting */
            if(m==n)
            {
                sprintf(filename,"pdf_%i",(int)step);
                print_1D_pdfs(modes, samples, NMCMC, filename, m);
            }
            
            /* Get Eigenvectors & Eigenvalues of Covariance matrix for plotting */
            else
            {
                sprintf(filename,"contours_%i",(int)step);
                print_2D_contours(modes, NMODE, filename, m, n);
            }
        }
    }
}


/**
 * \brief The Expectation-Maximization (EM) Algorithm
 *
 * Execulte one iteration of the EM algorithm to update fit to multivariate gaussian model (could be generalized)
 *
 * *E-step*: Compute probability for each sample to belong to each mode
 *
 * *M-step*: Recompute mean, covariance, and relative weight of newly weighted samples
 *
 * \param[in] samples data points \f$x_i\f$ and probabilities \f$p(x_i|k)\f$ for each Gaussian \f$k\f$
 * \param[in] NMCMC number of samples
 * \param[out] logL log-likelihood of input model
 * \param[out] BIC Bayesian Information Criteria (BIC) for input model
 * \param[in,out] modes parameters of each Gaussian including relative weights \f$\alpha_{k}\f$, updated by M-step.
 */
void expectation_maximization(struct Sample **samples, struct Mode **modes, size_t NMCMC, double *logL, double *BIC)
{
    size_t NP    = samples[0]->x->size;
    size_t NMODE = samples[0]->p->size;
    
    // aliases to structures
    struct Mode *M = NULL;
    struct Sample *s = NULL;
    
    // helper quantities for building sums etc.
    double norm;
    double mu=0.0;
    double C=0.0;
    double R;
    double p;
    
    /*
     E-step:
     compute probability for each sample to belong to each mode
     */
    
    /* compute p(x|mode) for each sample for each mode */
    for(size_t i=0; i<NMCMC; i++)
    {
        
        norm=0.0;
        s = samples[i];
        
        //loop over modes
        for(size_t k=0; k<NMODE; k++)
        {
            M = modes[k];
            
            //compute p(x|mode)
            p = multivariate_gaussian(s->x, M->mu, M->Cinv, M->detC);
            
            gsl_vector_set(s->p,k,p);
            gsl_vector_set(s->w,k,M->p*p);
            norm += M->p*p;
        }
        gsl_vector_scale(s->w,1./norm);
    }
    
    /* weigh the number of samples in each mode */
    for(size_t k=0; k<NMODE; k++)
    {
        modes[k]->Neff = 0;
        for(size_t i=0; i<NMCMC; i++) modes[k]->Neff += gsl_vector_get(samples[i]->w,k);
        if(modes[k]->Neff < 1.0)
        {
            exit(1);
        }
        modes[k]->p = modes[k]->Neff/(double)NMCMC;
    }
    
    /* check convergence with log likelihood & BIC */
    *logL = log_likelihood(modes, samples, NMCMC, NMODE);
    *BIC = -2.*(*logL) + (double)NMODE*((double)NP*((double)NP+3.)/2. + 1)*log((double)NMCMC);
    printf(" logL = %g,  BIC = %g     ",*logL, *BIC);
    
    
    /*
     M-Step:
     recompute mean, covariance, and relative weight of newly weighted samples
     */
    
    /* compute weighted mean and covariance for each mode */
    for(size_t k=0; k<NMODE; k++)
    {
        M = modes[k];
        
        //get new means
        for(size_t n=0; n<NP; n++)
        {
            mu = 0.0;
            
            //mu is a weighted average for each mode
            for(size_t i=0; i<NMCMC; i++)
                mu += (gsl_vector_get(samples[i]->w,k)/modes[k]->Neff) * gsl_vector_get(samples[i]->x,n);
            
            gsl_vector_set(M->mu,n,mu);
        }
        
        //get new covariance
        for(size_t m=0; m<NP; m++)
        {
            //loop over parameters again (wasteful, not using symmetry
            for(size_t n=0; n<NP; n++)
            {
                C = 0.0;
                
                //loop over samples
                for(size_t i=0; i<NMCMC; i++)
                {
                    double dx_m = gsl_vector_get(samples[i]->x,m) - gsl_vector_get(M->mu,m);
                    double dx_n = gsl_vector_get(samples[i]->x,n) - gsl_vector_get(M->mu,n);
                    
                    C += (gsl_vector_get(samples[i]->w,k)/modes[k]->Neff)*(dx_m)*(dx_n);
                }
                gsl_matrix_set(M->C,m,n,C);
            }
        }
        
        //invert new matrix to evaluate the probabilities
        invert_matrix(modes[k]->C, modes[k]->Cinv, &modes[k]->detC, &R);
    }
}


/**
 * \brief Main function for data handling and iterating through EM algorithm
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
    
    // chain samples
    struct Sample **samples = malloc(NMCMC*sizeof(struct Sample*));
    for(size_t n=0; n<NMCMC; n++)
    {
        samples[n] = malloc(sizeof(struct Sample));
        samples[n]->x = gsl_vector_alloc(NP);
        samples[n]->p = gsl_vector_alloc(NMODE);
        samples[n]->w = gsl_vector_alloc(NMODE);
    }
    
    // covariance matrices for different modes
    struct Mode **modes = malloc(NMODE*sizeof(struct Mode*));
    for(size_t n=0; n<NMODE; n++)
    {
        modes[n] = malloc(sizeof(struct Mode));
        modes[n]->mu = gsl_vector_alloc(NP);
        modes[n]->C = gsl_matrix_calloc(NP,NP);
        modes[n]->Cinv = gsl_matrix_calloc(NP,NP);
        modes[n]->evectors = gsl_matrix_calloc(NP,NP);
        modes[n]->evalues = gsl_vector_calloc(NP);
    }
    
    
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
            gsl_vector_set(samples[i]->x,n,value);
            column=strtok(NULL," ");
        }
    }
    
    
    /*
     Initialize Model
     */
    
    /* construct diagonal covariance matrix of full sample variances */
    double x_temp[NMCMC];
    double mean_temp, var_temp;
    for(size_t i=0; i<NP; i++)
    {
        for(size_t n=0; n<NMCMC; n++) x_temp[n] = gsl_vector_get(samples[n]->x,i);
        mean_temp = gsl_stats_mean(x_temp,1,NMCMC);
        var_temp  = gsl_stats_variance_m(x_temp,1,NMCMC, mean_temp);
        
        //set diagonals of C
        for(size_t n=0; n<NMODE; n++) gsl_matrix_set(modes[n]->C,i,i,var_temp);
        
    }
    
    /* place covariance matrices at random draws from the chain file */
    double R; //condition number of matrix
    for(size_t k=0; k<NMODE; k++)
    {
        //pick a sample from the chain to anchor each covariance matrix
        int fair_draw = (int)gsl_ran_flat(r,0,NMCMC);
        for(size_t n=0; n<NP; n++) gsl_vector_set(modes[k]->mu,n,gsl_vector_get(samples[fair_draw]->x,n));
        
        //set priors for each model
        modes[k]->p = (double)1./(double)NMODE;
        
        //get inverset, determinant, etc.
        invert_matrix(modes[k]->C, modes[k]->Cinv, &modes[k]->detC, &R);
    }
    
    /* EM Algorithm for Gaussian Mixture Models */
    size_t step=0;
    double logL,BIC;
    double BICmin = 1e60;
    while(step<NSTEP)
    {
        printProgress((double)(step+1)/NSTEP);
        expectation_maximization(samples, modes, NMCMC, &logL, &BIC);
        if(floor(BIC) < floor(BICmin))
        {
            BICmin = BIC;
            step=0;
        }
        step++;
    }
    printf("\n");

    print_model(modes, samples, NMCMC, logL, BIC, NMODE);

    
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
    for(size_t n=0; n<NMODE; n++)
    {
        gsl_vector_free(modes[n]->mu);
        gsl_matrix_free(modes[n]->C);
        gsl_matrix_free(modes[n]->Cinv);
        gsl_matrix_free(modes[n]->evectors);
        gsl_vector_free(modes[n]->evalues);
        free(modes[n]);
    }
    free(modes);

    free(LFLAG);

    
    return 0;
}



