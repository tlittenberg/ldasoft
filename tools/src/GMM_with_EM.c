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
 *  Library for Gaussian Mixture Model using the Expectation-Maximization Algorithm
 *
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

#include <util.h>

#include "GMM_with_EM.h"

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

void read_gmm_binary(struct GMM *gmm, char filename[])
{
    FILE *fptr = NULL;
    /* Read GMM results to binary for pick up by other processes */
    if( (fptr = fopen(filename,"rb"))!=NULL)
    {
        ufread(&gmm->NMODE, sizeof gmm->NMODE, 1, fptr);
        
        gmm->modes = malloc(gmm->NMODE*sizeof(struct MVG*));
        for(size_t n=0; n<gmm->NMODE; n++)
        {
            gmm->modes[n] = malloc(sizeof(struct MVG));
            alloc_MVG(gmm->modes[n],gmm->NP);
        }
        
        for(size_t n=0; n<gmm->NMODE; n++) read_MVG(gmm->modes[n],fptr);
        fclose(fptr);
    }
    else
    {
        fprintf(stderr,"Error reading %s\n",filename);
        fprintf(stderr,"Exiting to system\n");
        exit(1);
    }
}

void alloc_MVG(struct MVG *mode, size_t N)
{
    mode->size = N;
    mode->mu = gsl_vector_alloc(mode->size);
    mode->C = gsl_matrix_calloc(mode->size,mode->size);
    mode->L = gsl_matrix_calloc(mode->size,mode->size);
    mode->Cinv = gsl_matrix_calloc(mode->size,mode->size);
    mode->evectors = gsl_matrix_calloc(mode->size,mode->size);
    mode->evalues = gsl_vector_calloc(mode->size);
    mode->minmax = gsl_matrix_calloc(mode->size,2);
}

void free_MVG(struct MVG *mode)
{
    gsl_vector_free(mode->mu);
    gsl_matrix_free(mode->C);
    gsl_matrix_free(mode->L);
    gsl_matrix_free(mode->Cinv);
    gsl_matrix_free(mode->evectors);
    gsl_vector_free(mode->evalues);
    gsl_matrix_free(mode->minmax);
    free(mode);
}

void write_MVG(struct MVG *mode, FILE *fptr)
{
    //pack detC,p,and Neff into a vector to make life easier
    gsl_vector *temp = gsl_vector_alloc(3);
    gsl_vector_set(temp,0,mode->detC);
    gsl_vector_set(temp,1,mode->p);
    gsl_vector_set(temp,2,mode->Neff);

    //write 'em!
    gsl_vector_fwrite(fptr,mode->mu);
    gsl_matrix_fwrite(fptr,mode->C);
    gsl_matrix_fwrite(fptr,mode->L);
    gsl_matrix_fwrite(fptr,mode->Cinv);
    gsl_matrix_fwrite(fptr,mode->evectors);
    gsl_vector_fwrite(fptr,mode->evalues);
    gsl_vector_fwrite(fptr,temp);
    gsl_matrix_fwrite(fptr,mode->minmax);
    gsl_vector_free(temp);
}

void read_MVG(struct MVG *mode, FILE *fptr)
{
    //vector for holding packed detC,p,and Neff
    gsl_vector *temp = gsl_vector_alloc(3);

    //read 'em!
    gsl_vector_fread(fptr,mode->mu);
    gsl_matrix_fread(fptr,mode->C);
    gsl_matrix_fread(fptr,mode->L);
    gsl_matrix_fread(fptr,mode->Cinv);
    gsl_matrix_fread(fptr,mode->evectors);
    gsl_vector_fread(fptr,mode->evalues);
    gsl_vector_fread(fptr,temp);
    gsl_matrix_fread(fptr,mode->minmax);

    //unpack 'em!
    mode->detC = gsl_vector_get(temp,0);
    mode->p = gsl_vector_get(temp,1);
    mode->Neff = gsl_vector_get(temp,2);

    gsl_vector_free(temp);
}

void copy_MVG(struct MVG *origin, struct MVG *copy)
{
    copy->size = origin->size;
    gsl_vector_memcpy(copy->mu, origin->mu);
    gsl_matrix_memcpy(copy->C, origin->C);
    gsl_matrix_memcpy(copy->L, origin->L);
    gsl_matrix_memcpy(copy->Cinv, origin->Cinv);
    gsl_matrix_memcpy(copy->evectors, origin->evectors);
    gsl_matrix_memcpy(copy->minmax, origin->minmax);
    gsl_vector_memcpy(copy->evalues, origin->evalues);
    copy->detC = origin->detC;
    copy->p = origin->p;
    copy->Neff = origin->Neff;
}

double multivariate_gaussian(gsl_vector *x, struct MVG *mvg)
{
    
    size_t N = x->size;
    gsl_vector *dx = gsl_vector_alloc(N);
    gsl_vector *mu = mvg->mu;
    gsl_matrix *Cinv = mvg->Cinv;
    double detC = mvg->detC;
    
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
    
    return exp(-0.5*chi2)/sqrt(pow(2.0*M_PI,N)*detC);
}

void invert_gsl_matrix(gsl_matrix *A, gsl_matrix *Ainv, gsl_matrix *L, double *detA, double *R)
{
    int i,j;
    
    gsl_set_error_handler_off();
    
    //get size of matrix (assumed to be NxN)
    size_t N = A->size1;
    
    //some workspace
    gsl_permutation * permutation = gsl_permutation_alloc(N);
    
    //copy A into Ainv because LU decomposition destroys the matrix
    gsl_matrix_memcpy(Ainv,A);
    
    //cholesky decomposition
    gsl_linalg_cholesky_decomp(Ainv);
    
    //get condition number
    gsl_vector *work = gsl_vector_alloc(3*N);
    gsl_linalg_cholesky_rcond(Ainv, R, work);

    //inverse of A
    gsl_linalg_cholesky_invert(Ainv);
    
    //get deteriminant, need LU decomposition
    gsl_matrix_memcpy(L,A);
    gsl_linalg_LU_decomp(L,permutation,&i);
    *detA = gsl_linalg_LU_det(L,i);
    
    //recompute and save L
    gsl_matrix_memcpy(L,A);
    gsl_linalg_cholesky_decomp(L);
    for(i=0; i<N; i++) for(j=i+1; j<N; j++) gsl_matrix_set(L,i,j,0.0);

    //clean up
    gsl_vector_free (work);
    gsl_permutation_free (permutation);
}

void decompose_matrix(gsl_matrix *A, gsl_matrix *evec, gsl_vector *eval)
{
    //get size of matrix (assumed to be NxN)
    size_t N = A->size1;
    
    //get deteriminant, need LU decomposition
    gsl_matrix *Atemp = gsl_matrix_calloc(N,N);
    gsl_eigen_symmv_workspace * workspace = gsl_eigen_symmv_alloc (N);
    gsl_permutation * permutation = gsl_permutation_alloc(N);
    
    //copy A into Atemp because eigen_symmv destroys the matrix
    gsl_matrix_memcpy(Atemp,A);
    
    //the reason we're all here...
    gsl_eigen_symmv (Atemp, eval, evec, workspace);
    
    gsl_matrix_free (Atemp);
    gsl_eigen_symmv_free (workspace);
    gsl_permutation_free (permutation);
    
}

double log_likelihood(struct MVG **modes, struct Sample **samples, int NMCMC, int NMODE)
{
    
    double logL = 0.0;
    for(size_t i=0; i<NMCMC; i++)
    {
        double P = PMIN;
        for(size_t k=0; k<NMODE; k++)
        {
            P += modes[k]->p*gsl_vector_get(samples[i]->p,k);
        }
        logL += log(P);
    }
    return logL;
}

void print_1D_pdfs(struct MVG **modes, struct Sample **samples, size_t NMCMC, char root[], size_t ix)
{
    char filename[128];
    sprintf(filename,"%s_%i.dat",root,(int)ix);
    FILE *fptr = fopen(filename,"w");
    
    size_t NMODE = samples[0]->p->size;
    
    //get original parameter boundaries
    double pmin = gsl_matrix_get(modes[0]->minmax,ix,0);
    double pmax = gsl_matrix_get(modes[0]->minmax,ix,1);
    
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
        x = x0 + (double)n*dx;
        for(size_t k=0; k<NMODE; k++)
        {
            double mean = gsl_vector_get(modes[k]->mu,ix);
            double var  = gsl_matrix_get(modes[k]->C,ix,ix);
            
            /*
             Probability density p is weight * normal / Jacobian
             */
            p += modes[k]->p * exp( -0.5*(x-mean)*(x-mean)/var )/sqrt(2*M_PI*var) / ((pmax-pmin)*exp(-x)/pow(1. + exp(-x),2));
        }
        
        
        fprintf(fptr,"%.16g %.16g\n",sigmoid(x,pmin,pmax),p);
    }
    
    fclose(fptr);
    
    
    free(xvec);
}

void print_2D_pdfs(struct MVG **modes, struct Sample **samples, size_t NMCMC, char root[], size_t ix, size_t iy)
{
    char filename[128];
    sprintf(filename,"%s_%i_%i.dat",root,(int)ix,(int)iy);
    FILE *fptr = fopen(filename,"w");
    
    size_t NMODE = samples[0]->p->size;
    
    //get original parameter boundaries
    double pmin_x = gsl_matrix_get(modes[0]->minmax,ix,0);
    double pmax_x = gsl_matrix_get(modes[0]->minmax,ix,1);
    double pmin_y = gsl_matrix_get(modes[0]->minmax,iy,0);
    double pmax_y = gsl_matrix_get(modes[0]->minmax,iy,1);

    double *xvec = malloc(NMCMC*sizeof(double));
    double *yvec = malloc(NMCMC*sizeof(double));
    double xmin,xmax;
    double ymin,ymax;
    double x0 =  1e60;
    double xf = -1e60;
    double y0 =  1e60;
    double yf = -1e60;
    for(size_t k=0; k<NMODE; k++)
    {
        for(size_t i=0; i<NMCMC; i++) xvec[i] = gsl_vector_get(samples[i]->x,ix);
        gsl_stats_minmax(&xmin,&xmax,xvec,1,NMCMC);
        if(xmin<x0) x0 = xmin;
        if(xmax>xf) xf = xmax;
        
        for(size_t i=0; i<NMCMC; i++) yvec[i] = gsl_vector_get(samples[i]->x,iy);
        gsl_stats_minmax(&ymin,&ymax,yvec,1,NMCMC);
        if(ymin<y0) y0 = ymin;
        if(ymax>yf) yf = ymax;
    }
    
    double p;
    double x,y;
    double dx = (pmax_x-pmin_x)/100.;
    double dy = (pmax_y-pmin_y)/100.;
    for(int n=1; n<100; n++)
    {
        x = pmin_x + (double)n*dx;
        x = logit(x,pmin_x,pmax_x);
        
        for(int m=1; m<100; m++)
        {
            p = 0.0;

            y = pmin_y + (double)m*dy;
            y = logit(y,pmin_y,pmax_y);
            
            for(size_t k=0; k<NMODE; k++)
            {
                double mean_x = gsl_vector_get(modes[k]->mu,ix);
                double mean_y = gsl_vector_get(modes[k]->mu,iy);
                if(n==0 && m==0 && ix==2 && iy==5)
                {
                    printf("means: %g %g -> %g %g\n",mean_x,mean_y, sigmoid(mean_x,pmin_x,pmax_x),sigmoid(mean_y,pmin_y,pmax_y));
                }
                double C_xx  = gsl_matrix_get(modes[k]->C,ix,ix);
                double C_xy  = gsl_matrix_get(modes[k]->C,ix,iy);
                double C_yy  = gsl_matrix_get(modes[k]->C,iy,iy);
                double detC = C_xx*C_yy - C_xy*C_xy;
                double iC_xx = C_yy/detC;
                double iC_yy = C_xx/detC;
                double iC_xy = -C_xy/detC;

                /*
                 Probability density p is weight * normal / Jacobian
                 */
                p += modes[k]->p * exp( -0.5*( (x-mean_x)*(x-mean_x)*iC_xx + (y-mean_y)*(y-mean_y)*iC_yy + 2*(x-mean_x)*(y-mean_y)*iC_xy) )  * 1./(2.*M_PI*sqrt(detC)) / ((pmax_x-pmin_x)*exp(-x)/pow(1. + exp(-x),2)) / ((pmax_y-pmin_y)*exp(-y)/pow(1. + exp(-y),2)) ;
            }
            
            
            fprintf(fptr,"%.16g %.16g %.16g\n",sigmoid(x,pmin_x,pmax_x),sigmoid(y,pmin_y,pmax_y),p);
        }
    }
    
    fclose(fptr);

    
    free(xvec);
    free(yvec);
}

void print_2D_contours(struct MVG **modes, size_t NMODE, char root[], size_t x1, size_t x2)
{
    char filename[128];
    char filename2[128];
    FILE *fptr = NULL;
    FILE *fptr2 = NULL;

    struct MVG **submodes = malloc(NMODE*sizeof(struct MVG*));
    for(size_t k=0; k<NMODE; k++)
    {
        submodes[k] = malloc(sizeof(struct MVG));
        alloc_MVG(submodes[k], 2);
    }
    
    //get original parameter boundaries
    double pmin[2] = {gsl_matrix_get(modes[0]->minmax,x1,0),gsl_matrix_get(modes[0]->minmax,x2,0)};
    double pmax[2] = {gsl_matrix_get(modes[0]->minmax,x1,1),gsl_matrix_get(modes[0]->minmax,x2,1)};

    
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
        sprintf(filename2,"remapped_%s_%i_%i_%i.dat",root,(int)x1,(int)x2,(int)k);
        fptr2=fopen(filename2,"w");
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
            fprintf(fptr2,"%.16lg %.16lg ",x,y);
            x = sigmoid(x,pmin[0],pmax[0]);
            y = sigmoid(y,pmin[1],pmax[1]);
            fprintf(fptr,"%.16lg %.16lg ",x,y);
            
            x = 2.*( Rx*cos(angle)*cos(theta) + Ry*sin(angle)*sin(theta) ) + Cx;
            y = 2.*(-Rx*cos(angle)*sin(theta) + Ry*sin(angle)*cos(theta) ) + Cy;
            fprintf(fptr2,"%.16lg %.16lg ",x,y);
            x = sigmoid(x,pmin[0],pmax[0]);
            y = sigmoid(y,pmin[1],pmax[1]);
            fprintf(fptr,"%.16lg %.16lg ",x,y);
            
            x = 3.*( Rx*cos(angle)*cos(theta) + Ry*sin(angle)*sin(theta) ) + Cx;
            y = 3.*(-Rx*cos(angle)*sin(theta) + Ry*sin(angle)*cos(theta) ) + Cy;
            fprintf(fptr2,"%.16lg %.16lg ",x,y);
            x = sigmoid(x,pmin[0],pmax[0]);
            y = sigmoid(y,pmin[1],pmax[1]);
            fprintf(fptr,"%.16lg %.16lg ",x,y);
            
            fprintf(fptr,"\n");
            fprintf(fptr2,"\n");
        }
        fclose(fptr);
        fclose(fptr2);
    }
    
    for(size_t k=0; k<NMODE; k++) free_MVG(submodes[k]);
    free(submodes);
}

void print_model(struct MVG **modes, struct Sample **samples, size_t NMCMC, double logL, double BIC, size_t step)
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
                
                sprintf(filename,"2Dpdf_%i",(int)step);
                print_2D_pdfs(modes, samples, NMCMC, filename, m, n);

            }
        }
    }
}


int expectation_maximization(struct Sample **samples, struct MVG **modes, size_t NMCMC, double *logL, double *BIC)
{
    size_t NP    = samples[0]->x->size;
    size_t NMODE = samples[0]->p->size;
    
    // aliases to structures
    struct MVG *M = NULL;
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
            p = multivariate_gaussian(s->x, M);
            p = p > PMIN ? p : PMIN;

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
        if(modes[k]->Neff < 1.0 || modes[k]->Neff != modes[k]->Neff) return 1;
        modes[k]->p = modes[k]->Neff/(double)NMCMC;
    }
    
    /* check convergence with log likelihood & BIC */
    *logL = log_likelihood(modes, samples, NMCMC, NMODE);
    *BIC = -2.*(*logL) + (double)NMODE*((double)NP*((double)NP+3.)/2. + 1)*log((double)NMCMC);
    //printf(" logL = %g,  BIC = %g     ",*logL, *BIC);
    
    
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
        invert_gsl_matrix(modes[k]->C, modes[k]->Cinv, modes[k]->L, &modes[k]->detC, &R);
    }
    
    return 0;
}

int GMM_with_EM(struct MVG **modes, struct Sample **samples, size_t NMCMC, size_t NSTEP, gsl_rng *r, double *logL, double *BIC)
{
    size_t NP    = samples[0]->x->size;
    size_t NMODE = samples[0]->p->size;

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
        
        //get inverse, determinant, etc.
        invert_gsl_matrix(modes[k]->C, modes[k]->Cinv, modes[k]->L, &modes[k]->detC, &R);
    }
    
    /* EM Algorithm for Gaussian Mixture Models */
    size_t step=0;
    double BICmin = 1e60;
    while(step<NSTEP)
    {
        //printProgress((double)(step+1)/NSTEP);
        if(expectation_maximization(samples, modes, NMCMC, logL, BIC)) return 1;
        else
        {
            if(floor(*BIC) < floor(BICmin))
            {
                BICmin = *BIC;
                step=0;
            }
            step++;
        }
    }
    //printf("\n");
    return 0;
}

double logit(double x,double xmin,double xmax)
{
    return x;//log( (x-xmin)/(xmax-x) );
}

double sigmoid(double x,double xmin,double xmax)
{
    return x;//xmin + (1./(1. + exp(-x)))*(xmax - xmin);
}


double dsigmoid(double x, double xmin, double xmax)
{
    //double expyn = exp(-x);
    return 1.0;//(xmax-xmin) * expyn/((1.+expyn)*(1.+expyn));
}


void logit_mapping(gsl_vector *x_vec, gsl_vector *y_vec, double xmin, double xmax)
{
    size_t N = x_vec->size;

    for(size_t n=0; n<N; n++)
    {
        double x = gsl_vector_get(x_vec,n);
        if(x<xmin || x>xmax)
        {
            /*
             This is a safety check incase there are small differences
             between the prior setup by gb_catalog and the one used by
             the analysis.  This was happening for the fdot parameter
             which is a strong function of fmin, and inconsistencies
             in fmin between global_fit and gb_catalog runs.
             */
            if(x>xmax) x -= 2.*(x-xmax);
            if(x<xmin) x += 2.*(xmin-x);
        }
        double y = logit(x,xmin,xmax);//log( (x - xmin)/(xmax - x) );
        gsl_vector_set(y_vec,n,y);
    }
}

/* sigmoid function */
void sigmoid_mapping(gsl_vector *x_vec, gsl_vector *y_vec, double xmin, double xmax)
{
    size_t N = x_vec->size;

    for(size_t n=0; n<N; n++)
    {
        double y = gsl_vector_get(y_vec,n);
        double x = sigmoid(y,xmin,xmax);//xmin + (1./(1. + exp(-y)))*(xmax - xmin);
        gsl_vector_set(x_vec,n,x);
    }
}
