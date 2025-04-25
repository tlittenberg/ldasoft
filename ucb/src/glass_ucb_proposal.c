/*
 *  Copyright (C) 2019 Tyson B. Littenberg (MSFC-ST12), Neil J. Cornish, Kristen Lackeos
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

#include <glass_utils.h>

#include "glass_ucb_catalog.h"
#include "glass_ucb_model.h"
#include "glass_ucb_prior.h"
#include "glass_ucb_io.h"
#include "glass_ucb_waveform.h"
#include "glass_ucb_fstatistic.h"
#include "glass_ucb_proposal.h"


#define FIXME 0
#define SNRCAP 100.0 /* SNR cap on logL */

static void write_Fstat_animation(double fmin, double T, struct Proposal *proposal, char runDir[], double minp)
{
    char filename[MAXSTRINGSIZE];
    sprintf(filename,"%s/fstat/Fstat.gpi",runDir);
    FILE *fptr = fopen(filename,"w");
    
    fprintf(fptr,"!rm Fstat.mp4\n");
    fprintf(fptr,"set terminal pngcairo size 1024,512 enhanced font 'Verdana,12'\n");
    
    fprintf(fptr,"\n");
    
    fprintf(fptr,"# perula palette\n");
    fprintf(fptr,"set palette defined (0 '#352a87', 1 '#0363e1', 2 '#1485d4', 3 '#06a7c6', 4 '#38b99e', 5 '#92bf73', 6 '#d9ba56', 7 '#fcce2e', 8 '#f9fb0e')\n");
    
    fprintf(fptr,"\n");
    
    fprintf(fptr,"\n");
    
    fprintf(fptr,"unset colorbox\n");
    fprintf(fptr,"unset key\n");
    fprintf(fptr,"unset label\n");
    fprintf(fptr,"set tics scale 2 font ',8'\n");
    fprintf(fptr,"set xlabel '{/Symbol f}\n");
    fprintf(fptr,"set ylabel 'cos({/Symbol q})' offset 1\n");
    
    fprintf(fptr,"\n");
    
    fprintf(fptr,"\n");
    
    fprintf(fptr,"df0 = %g\n",proposal->matrix[0][1]/T);
    fprintf(fptr,"dph = 0.5*2.*pi/%i.\n",(int)proposal->matrix[1][0]);
    fprintf(fptr,"dth = 0.5*2./%i.\n"   ,(int)proposal->matrix[2][0]);
    
    fprintf(fptr,"\n");
    
    fprintf(fptr,"set pm3d map\n");
    
    fprintf(fptr,"\n");
    
    fprintf(fptr,"do for [ii=0:%i-1:+1]{\n",(int)proposal->matrix[0][0]);
    fprintf(fptr,"  set output sprintf('fstat_frame%%05.0f.png',ii)\n");
    fprintf(fptr,"  set size 1,1\n");
    //  fprintf(fptr,"  set multiplot title sprintf('fstat-frame%%05.0f.png',ii)\n");
    fprintf(fptr,"  set multiplot title sprintf('f = %%f ',(ii*df0+%g)*1000)\n",fmin);
    fprintf(fptr,"  set size 0.5,1\n");
    
    fprintf(fptr,"\n");
    
    fprintf(fptr,"  set origin 0,0\n");
    fprintf(fptr,"  set cbrange [%lg:%lg]\n",minp,proposal->maxp);
    fprintf(fptr,"  input = sprintf('skymap_%%05.0f.dat',ii)\n");
    fprintf(fptr,"  splot [0:2*pi] [-1:1] input u ($2+dph):($1+dth):3\n");
    
    fprintf(fptr,"\n");
    
    fprintf(fptr,"  set origin 0.5,0\n");
    fprintf(fptr,"  set cbrange [%lg:%lg]\n",log(minp),log(proposal->maxp));
    fprintf(fptr,"  input = sprintf('skymap_%%05.0f.dat',ii)\n");
    fprintf(fptr,"  splot [0:2*pi] [-1:1] input u ($2+dph):($1+dth):(log($3))\n");
    
    fprintf(fptr,"  unset multiplot\n");
    
    fprintf(fptr,"}\n");
    
    fprintf(fptr,"\n");
    
    fprintf(fptr,"# Create animation\n");
    fprintf(fptr,"CMD = 'ffmpeg -r 256 -f image2 -s 3840x1080 -start_number 0 -i fstat_frame%%05d.png -vframes %i -vcodec libx264 -crf 25  -pix_fmt yuv420p Fstat.mp4'\n",(int)proposal->matrix[0][0]);
    fprintf(fptr,"system(CMD)\n");
    
    fprintf(fptr,"\n");
    
    fprintf(fptr,"!rm *.png\n");
    fclose(fptr);
}

void setup_frequency_proposal(struct Data *data, struct Flags *flags)
{
    int BW = 20;
    double *power = data->p;
    double total  = 0.0;
    char filename[MAXSTRINGSIZE];
    sprintf(filename,"%s/data/frequency_proposal.dat",flags->runDir);
    FILE *temp = fopen(filename,"w");
    for(int i=0; i<data->NFFT-BW; i++)
    {
        power[i]=0.0;
        for(int n=i; n<i+BW; n++)
        {
            double SnA = data->noise->C[0][0][n];
            double SnE = data->noise->C[1][1][n];
            
            double AA = data->tdi->A[2*n]*data->tdi->A[2*n]+data->tdi->A[2*n+1]*data->tdi->A[2*n+1];
            double EE = data->tdi->E[2*n]*data->tdi->E[2*n]+data->tdi->E[2*n+1]*data->tdi->E[2*n+1];
            
            power[i] += AA/SnA + EE/SnE;
            total += power[i];
        }
    }
    for(int i=data->NFFT-BW; i<data->NFFT; i++)
    {
        power[i] = power[data->NFFT-BW-1];
        total += power[i];
    }
    
    data->pmax = 0.0;
    for(int i=0; i<data->NFFT; i++)
    {
        fprintf(temp,"%i %lg\n",i,power[i]);
        if(power[i]>data->pmax) data->pmax = power[i];
    }
    fclose(temp);
    
    
    //also get SNR^2 of data
    total = 0.0;
    for(int n=0; n<data->NFFT; n++)
    {
        double SnA = data->noise->C[0][0][n];
        double SnE = data->noise->C[1][1][n];

        double AA = data->tdi->A[2*n]*data->tdi->A[2*n]+data->tdi->A[2*n+1]*data->tdi->A[2*n+1];
        double EE = data->tdi->E[2*n]*data->tdi->E[2*n]+data->tdi->E[2*n+1]*data->tdi->E[2*n+1];
        
        total += AA/SnA + EE/SnE;
    }
    
    data->SNR2 = total - data->NFFT;
    printf("total=%g, N=%i\n",total, data->NFFT);
    if(data->SNR2<0.0)data->SNR2=0.0;
    printf("data-based SNR^2:  %g (%g)\n", data->SNR2, sqrt(data->SNR2));
    
}

void print_acceptance_rates(struct Proposal **proposal, int NProp, int ic, FILE *fptr)
{
    fprintf(fptr,"Acceptance rates for chain %i:\n", ic);
    fprintf(fptr," MCMC\n");
    for(int n=0; n<NProp; n++)
    {
        if(proposal[n]->weight > 0) fprintf(fptr,"   %.1e  [%s]\n", (double)proposal[n]->accept[ic]/(double)proposal[n]->trial[ic],proposal[n]->name);
    }
    fprintf(fptr," RJMCMC\n");
    for(int n=0; n<NProp; n++)
    {
        if(proposal[n]->rjweight > 0) fprintf(fptr,"   %.1e  [%s]\n", (double)proposal[n]->accept[ic]/(double)proposal[n]->trial[ic],proposal[n]->name);
    }
}

double draw_from_spectrum(struct Data *data, struct Model *model, struct Source *source, UNUSED struct Proposal *proposal, double *params, unsigned int *seed)
{
    //TODO: Work in amplitude
    
    //rejections ample for f
    double alpha;
    int q;
    do
    {
        params[0] = (double)model->prior[0][0] + rand_r_U_0_1(seed)*(double)(model->prior[0][1]-model->prior[0][0]);
        alpha     = rand_r_U_0_1(seed)*data->pmax;
        q = (int)(params[0]-data->qmin);
    }while(data->p[q]<alpha);

    //random draws for other parameters
    for(int n=1; n<UCB_MODEL_NP; n++) params[n] = model->prior[n][0] + rand_r_U_0_1(seed)*(model->prior[n][1]-model->prior[n][0]);
    
    return 0;
}

double draw_from_prior(struct Data *data, struct Model *model, struct Source *source, struct Proposal *proposal, double *params, unsigned int *seed)
{
    return draw_from_uniform_prior(data,model,source,proposal,params,seed);
}

double draw_from_gmm_prior(struct Data *data, struct Model *model, struct Source *source, struct Proposal *proposal, double *params, unsigned int *seed)
{
    double ran_no[UCB_MODEL_NP];
    
    //choose which entry
    int ngmm = (int)floor(rand_r_U_0_1(seed)*proposal->Ngmm);
    int NMODES = (int)proposal->gmm[ngmm]->NMODE;
    
    struct MVG **modes = proposal->gmm[ngmm]->modes;
    struct MVG *mode = NULL;
    
    //pick which mode
    int k;
    double p = 1.;
    do {
        k = (int)floor(rand_r_U_0_1(seed)*NMODES);
        mode = modes[k];
        p = rand_r_U_0_1(seed);
    } while (p>mode->p);
    

    //get vector of gaussian draws n;  y_i = x_mean_i + sum_j Lij^-1 * n_j
    for(int n=0; n<UCB_MODEL_NP; n++)
    {
        ran_no[n] = rand_r_N_0_1(seed);
    }
    
    double x[UCB_MODEL_NP];
    //the matrix multiplication...
    for(int n=0; n<UCB_MODEL_NP; n++)
    {
        //start at mean
        x[n] = mode->mu[n];
        
        //add contribution from each row of invC
        for(int m=0; m<UCB_MODEL_NP; m++) x[n] += ran_no[m]*mode->L[n][m];
        
        //map params from R back to interval
        double min = mode->minmax[n][0];
        double max = mode->minmax[n][1];
        
        x[n] = sigmoid(x[n],min,max);
    }
    
    //map to parameters
    source->f0 = x[0];
    source->costheta = x[1];
    source->phi = x[2];
    source->amp = exp(x[3]);
    source->cosi = x[4];
    source->psi = x[5];
    source->phi0 = x[6];
    if(UCB_MODEL_NP>7) source->dfdt = x[7];
    if(UCB_MODEL_NP>8) source->d2fdt2 = x[8];
    map_params_to_array(source, params, data->T);
    
    return gmm_prior_density(data, model, source, proposal, params);
}

double gmm_prior_density(struct Data *data, struct Model *model, struct Source *source, struct Proposal *proposal, double *params)
{
    double p = 0;
    
    /* sum over modes */
    for(int n=0; n<proposal->Ngmm; n++)
    {
        p += exp(evaluate_gmm_prior(data, proposal->gmm[n], params));
    }
        
    return log(p/(double)proposal->Ngmm);
}

double draw_from_uniform_prior(UNUSED struct Data *data, struct Model *model, UNUSED struct Source *source, UNUSED struct Proposal *proposal, double *params, unsigned int *seed)
{
    
    double logQ = 0.0;
    int n;
    
    //frequency
    n = 0;
    params[n] = model->prior[n][0] + rand_r_U_0_1(seed)*(model->prior[n][1]-model->prior[n][0]);
    logQ -= model->logPriorVolume[n];
    
    //sky location
    n = 1;
    params[n] = model->prior[n][0] + rand_r_U_0_1(seed)*(model->prior[n][1]-model->prior[n][0]);
    logQ -= model->logPriorVolume[n];
    
    n = 2;
    params[n] = model->prior[n][0] + rand_r_U_0_1(seed)*(model->prior[n][1]-model->prior[n][0]);
    logQ -= model->logPriorVolume[n];
    
    //amplitude
    n = 3;
    params[n] = model->prior[n][0] + rand_r_U_0_1(seed)*(model->prior[n][1]-model->prior[n][0]);
    logQ -= model->logPriorVolume[n];

    //inclination
    n = 4;
    params[n] = model->prior[n][0] + rand_r_U_0_1(seed)*(model->prior[n][1]-model->prior[n][0]);
    logQ -= model->logPriorVolume[n];
    
    //polarization
    n = 5;
    params[n] = model->prior[n][0] + rand_r_U_0_1(seed)*(model->prior[n][1]-model->prior[n][0]);
    logQ -= model->logPriorVolume[n];
    
    //phase
    n = 6;
    params[n] = model->prior[n][0] + rand_r_U_0_1(seed)*(model->prior[n][1]-model->prior[n][0]);
    logQ -= model->logPriorVolume[n];
    
    //fdot
    if(UCB_MODEL_NP>7)
    {
        n = 7;
        params[n] = model->prior[n][0] + rand_r_U_0_1(seed)*(model->prior[n][1]-model->prior[n][0]);
        logQ -= model->logPriorVolume[n];
    }
    
    //f-double-dot
    if(UCB_MODEL_NP>8)
    {
        n = 8;
        params[n] = model->prior[n][0] + rand_r_U_0_1(seed)*(model->prior[n][1]-model->prior[n][0]);
        logQ -= model->logPriorVolume[n];
    }
    
    
    return logQ;
}

double uniform_prior_density(struct Data *data, struct Model *model, UNUSED struct Source *source, UNUSED struct Proposal *proposal, double *params)
{
    
    double logQ = 0.0;
    int n;
    
    //frequency
    n = 0;
    logQ -= model->logPriorVolume[n];
    
    //sky location
    n = 1;
    logQ -= model->logPriorVolume[n];
    
    n = 2;
    logQ -= model->logPriorVolume[n];
    
    //amplitude
    n = 3;
    logQ -= model->logPriorVolume[n];
    
    //inclination
    n = 4;
    logQ -= model->logPriorVolume[n];
    
    //polarization
    n = 5;
    logQ -= model->logPriorVolume[n];
    
    //phase
    n = 6;
    logQ -= model->logPriorVolume[n];
    
    //fdot
    if(UCB_MODEL_NP>7)
    {
        n = 7;
        logQ -= model->logPriorVolume[n];
    }
    
    //f-double-dot
    if(UCB_MODEL_NP>8)
    {
        n = 8;
        logQ -= model->logPriorVolume[n];
    }
    
    
    return logQ;
}


double draw_from_extrinsic_prior(UNUSED struct Data *data, struct Model *model, UNUSED struct Source *source, UNUSED struct Proposal *proposal, double *params, unsigned int *seed)
{
    double logP = 0.0;
    
    for(int n=1; n<3; n++)
    {
        params[n] = model->prior[n][0] + rand_r_U_0_1(seed)*(model->prior[n][1]-model->prior[n][0]);
        logP -= model->logPriorVolume[n];
    }
    for(int n=4; n<7; n++)
    {
        params[n] = model->prior[n][0] + rand_r_U_0_1(seed)*(model->prior[n][1]-model->prior[n][0]);
        logP -= model->logPriorVolume[n];
    }
    
    return logP;
}

double draw_from_galaxy_prior(struct Model *model, struct Prior *prior, double *params, unsigned int *seed)
{
    double **uniform_prior = model->prior;
    double logP=-INFINITY;
    double alpha=prior->skymaxp;
    
    do
    {
        //sky location
        params[1] = model->prior[1][0] + rand_r_U_0_1(seed)*(model->prior[1][1]-model->prior[1][0]);
        params[2] = model->prior[2][0] + rand_r_U_0_1(seed)*(model->prior[2][1]-model->prior[2][0]);
        
        //map costheta and phi to index of skyhist array
        int i = (int)floor((params[1]-uniform_prior[1][0])/prior->dcostheta);
        int j = (int)floor((params[2]-uniform_prior[2][0])/prior->dphi);
        
        int k = i*prior->nphi + j;
        
        logP = prior->skyhist[k];
        alpha = log(rand_r_U_0_1(seed)*exp(prior->skymaxp));
    }while(alpha>logP);
    return logP;
}

double draw_calibration_parameters(struct Data *data, struct Model *model, unsigned int *seed)
{
    double logP = 0.0;
    
    /*apply calibration error to full signal model
    double dA,dphi;
    switch(data->Nchannel)
    {
        case 1:
            
            //amplitude
            dA = gsl_ran_gaussian(seed,CAL_SIGMA_AMP);
            model->calibration->dampX = dA;
            logP += log(gsl_ran_gaussian_pdf(dA,CAL_SIGMA_AMP));
            
            //phase
            dphi = 0.0;//gsl_ran_gaussian(seed,CAL_SIGMA_PHASE);
            model->calibration->dphiX = dphi;
            logP += 0.0;//log(gsl_ran_gaussian_pdf(dphi,CAL_SIGMA_PHASE));
            
            break;
        case 2:
            
            //amplitude
            dA = gsl_ran_gaussian(seed,CAL_SIGMA_AMP);
            model->calibration->dampA = dA;
            logP += log(gsl_ran_gaussian_pdf(dA,CAL_SIGMA_AMP));
            
            //phase
            dphi = 0.0;//gsl_ran_gaussian(seed,CAL_SIGMA_PHASE);
            model->calibration->dphiA = dphi;
            logP += 0.0;//log(gsl_ran_gaussian_pdf(dphi,CAL_SIGMA_PHASE));
            
            //amplitude
            dA = gsl_ran_gaussian(seed,CAL_SIGMA_AMP);
            model->calibration->dampE = dA;
            logP += log(gsl_ran_gaussian_pdf(dA,CAL_SIGMA_AMP));
            
            //phase
            dphi = 0.0;//gsl_ran_gaussian(seed,CAL_SIGMA_PHASE);
            model->calibration->dphiE = dphi;
            logP += 0.0;//log(gsl_ran_gaussian_pdf(dphi,CAL_SIGMA_PHASE));
            
            break;
        default:
            break;
    }end switch*/

    return logP;
}

double draw_from_fisher(UNUSED struct Data *data, struct Model *model, struct Source *source, UNUSED struct Proposal *proposal, double *params, unsigned int *seed)
{
    int i,j;
    //double sqNP = sqrt((double)source->NP);
    double Amps[UCB_MODEL_NP];
    double jump[UCB_MODEL_NP];
    
    //draw the eigen-jump amplitudes from N[0,1] scaled by evalue & dimension
    for(i=0; i<UCB_MODEL_NP; i++)
    {
        Amps[i] = rand_r_N_0_1(seed)/sqrt(source->fisher_evalue[i]);
        jump[i] = 0.0;
    }
    
    //choose one eigenvector to jump along
    i = (int)(rand_r_U_0_1(seed)*(double)UCB_MODEL_NP);
    for (j=0; j<UCB_MODEL_NP; j++) jump[j] += Amps[i]*source->fisher_evectr[j][i];
    
    //check jump value, set to small value if singular
    for(i=0; i<UCB_MODEL_NP; i++) if(jump[i]!=jump[i]) jump[i] = 0.01*source->params[i];
    
    //jump from current position
    for(i=0; i<UCB_MODEL_NP; i++) params[i] = source->params[i] + jump[i];
    
    //safety check for cos(latitude) parameters
    //cosine co-latitude
    if(params[1] >=  1.) params[1] = source->params[1] - jump[1];
    if(params[1] <= -1.) params[1] = source->params[1] - jump[1];
    //cosine inclination
    if(params[4] >=  1.) params[4] = source->params[4] - jump[4];
    if(params[4] <= -1.) params[4] = source->params[4] - jump[4];

    for(int j=0; j<UCB_MODEL_NP; j++)
    {
        if(params[j]!=params[j])
        {
            fprintf(stderr,"draw_from_fisher: params[%i]=%g, N[%g,%g]\n",j,params[j],source->params[j],jump[j]);
            return -INFINITY;
        }
    }
    
    //not updating Fisher between moves, proposal is symmetric
    return 0.0;
}

double draw_from_cdf(UNUSED struct Data *data, struct Model *model, struct Source *source, struct Proposal *proposal, double *params, unsigned int *seed)
{
    int N = proposal->size;
    double **cdf = proposal->matrix;
    
    for(int n=0; n<UCB_MODEL_NP; n++)
    {
        
        //draw n from U[0,N]
        double c_prime = rand_r_U_0_1(seed)*(double)(N-1);
        
        //map to nearest sample
        int c_minus = (int)floor(c_prime);
        int c_plus  = c_minus+1;
        
        //get value of pdf at cprime
        double p = 1./(cdf[n][c_plus] - cdf[n][c_minus]);
        
        //interpolate to get new parameter
        params[n] = (c_prime - c_minus)*(1./p) + cdf[n][c_minus];
        
        if(params[n]<model->prior[n][0] || params[n]>=model->prior[n][1]) return -INFINITY;
        
    }
    
    return cdf_density(data, model, source, proposal, params);
}

double cdf_density(UNUSED struct Data *data, struct Model *model, struct Source * source, struct Proposal *proposal, double *params)
{
    int N = proposal->size;
    double **cdf = proposal->matrix;
    double logP=0.0;
    
    int i,j;
    
    for(int n=0; n<UCB_MODEL_NP; n++)
    {
        if(params[n]<model->prior[n][0] || params[n]>=model->prior[n][1]) return -INFINITY;
        
        //find samples either end of p-value
        if(params[n]<cdf[n][0] || params[n]>= cdf[n][N-1]) return -INFINITY;
        
        i=binary_search(cdf[n],0,N,params[n]);
        j=i+1;
        while(cdf[n][j]==cdf[n][i]) j++;
        logP += log(  ((double)(j-i)/N) /  (cdf[n][j]-cdf[n][i])  );
    }
    return logP;
}

double draw_from_cov(UNUSED struct Data *data, struct Model *model, struct Source *source, struct Proposal *proposal, double *params, unsigned int *seed)
{
    double ran_no[UCB_MODEL_NP];

    int ns=(int)floor(rand_r_U_0_1(seed)*proposal->size);

    //pick which mode to propose to
    int mode;
    if (rand_r_U_0_1(seed) > (1.0-proposal->vector[4*ns])) mode = 2*ns;
    else mode = 2*ns+1;
    
    //define some helper pointers for ease of reading
    double *mean  = proposal->matrix[mode];
    double **Lij  = proposal->tensor[mode];
    
    //scale NP-dimensional jump to 1-sigma of the joint distribution
    double scale = 1.;
    
    //get vector of gaussian draws n;  y_i = x_mean_i + sum_j Lij^-1 * n_j
    for(int n=0; n<UCB_MODEL_NP; n++)
    {
        ran_no[n] = rand_r_N_0_1(seed)*scale;
    }
    
    //the matrix multiplication...
    for(int n=0; n<UCB_MODEL_NP; n++)
    {
        //start at mean
        params[n] = mean[n];
        
        //add contribution from each row of invC
        for(int k=0; k<UCB_MODEL_NP; k++) params[n] += ran_no[k]*Lij[n][k];
        
    }
        
    return cov_density(data, model, source, proposal, params);
}

/* when wrapping azimuthal parameters the sign matters! */
static double phase_wrapper(double x, double cycle)
{
    double a = fabs(x);
    double b = fabs(x + cycle);
    double c = fabs(x - cycle);
    double temp,min;
    
    temp = (a < b)    ? a : b;
    min =  (c < temp) ? c : temp;
    
    if(min==a)
        return x;
    else if(min==b)
        return x + cycle;
    else
        return x - cycle;
}
double cov_density(UNUSED struct Data *data, struct Model *model, struct Source * source, struct Proposal *proposal, double *params)
{
        
    /*
     Check that parameters are inside the prior range
     */
    
    //map angles over periodic boundary conditions
    //longitude
    if(params[2]>PI2)  params[2]-=PI2;
    if(params[2]<0.0)  params[2]+=PI2;
    
    //phase
    if(params[6]>PI2)  params[6]-=PI2;
    if(params[6]<0.0)  params[6]+=PI2;
    
    //psi
    if(params[5]>M_PI) params[5]-=M_PI;
    if(params[5]<0.0)  params[5]+=M_PI;
    
    
    //rejection sample everything else
    for(int n=0; n<UCB_MODEL_NP; n++)
    {
        if(params[n]<model->prior[n][0] || params[n]>model->prior[n][1])  return -INFINITY;
    }
    
    /*
     if everything is in bounds, evaluate the probability (density)
     */
    int Nmodes = proposal->size*2;
    double delta_x[UCB_MODEL_NP];
    
    //helper pointers
    double *mean    = NULL;
    double **invCij = NULL;
    
    double weight; //contribution of mode
    double norm;   //normalization of multivariate Gaussian
    double arg;    //argument of exponential
    
    long double p = 0.0;
    
    for(int i=0; i<Nmodes; i++)
    {
        //unpack proposal structure
        weight = proposal->vector[i*2];
        norm   = proposal->vector[i*2+1];
        mean   = proposal->matrix[i];
        invCij = proposal->tensor[Nmodes+i];
        
        //compute distances between mode and current params
        for(int n=0; n<UCB_MODEL_NP; n++)
        {
            delta_x[n] = params[n]-mean[n];
            
            //map parameters periodic on U[0,2pi]
            if(n==2) delta_x[n] = phase_wrapper(delta_x[n],PI2); //longitude
            if(n==6) delta_x[n] = phase_wrapper(delta_x[n],PI2); //phase
            
            //map parameters periodic on U[0,pi]
            if(n==5) delta_x[n] = phase_wrapper(delta_x[n],M_PI); //psi
            
        }
        
        //assemble argument of exponential
        arg = 0.0;
        for(int n=0; n<UCB_MODEL_NP; n++)
            for(int m=0; m<UCB_MODEL_NP; m++)
                arg += delta_x[n] * invCij[n][m] * delta_x[m];
        
        p += weight * norm * exp(-0.5*arg);
    }
    
    return log(p);
}

double fm_shift(struct Data *data, struct Model *model, struct Source *source, struct Proposal *proposal, double *params, unsigned int *seed)
{
    //doppler modulation frequency (in bins)
    double fm = data->T/YEAR;
    
    //perturb frequency by 1 fm
    int scale;
    do scale = (int)floor(-3. + 7.*rand_r_U_0_1(seed));
    while (scale == 0);
    
    params[0] += scale*fm;
    params[2] = params[2] + -0.2 + rand_r_U_0_1(seed)*0.4;
    params[7] -= 2.*scale*fm*fm;
    
    //perturb all parameters by Fisher Matrix (in liu of Jacaboian)
    draw_from_fisher(data, model, source, proposal, params, seed);
    
    //fm shift is symmetric
    return 0.0;
}

double psi_phi_jump(UNUSED struct Data *data, UNUSED struct Model *model, struct Source *source, UNUSED struct Proposal *proposal, double *params, unsigned int *seed)
{
    //proposal taking advantage of degeneracy between {psi,phi_0} -> {psi+/-pi/2,phi_0+/-pi}
    double scale = rand_r_U_0_1(seed)*PI2;
    double d_phi = scale;
    double d_psi = 0.5*scale;
    
    if(rand_r_U_0_1(seed) < 0.5 ) d_phi *= -1.;
    if(rand_r_U_0_1(seed) < 0.5 ) d_psi *= -1.;
    
    //jump from current position
    for(int i=0; i<UCB_MODEL_NP; i++) params[i] = source->params[i];
    
    params[5] += d_psi;
    params[6] += d_phi;
    
    //psi-phi jump is symmetric
    return 0.0;
    
}

void initialize_proposal(struct Orbit *orbit, struct Data *data, struct Prior *prior, struct Chain *chain, struct Flags *flags, struct Catalog *catalog, struct Proposal **proposal, int NMAX)
{
    int NC = chain->NC;
    double check  =0.0;
    double rjcheck=0.0;
    
    for(int i=0; i<UCB_PROPOSAL_NPROP; i++)
    {
        proposal[i] = malloc(sizeof(struct Proposal));

        proposal[i]->trial  = malloc(NC*sizeof(int));
        proposal[i]->accept = malloc(NC*sizeof(int));
        
        for(int ic=0; ic<NC; ic++)
        {
            proposal[i]->trial[ic]  = 1;
            proposal[i]->accept[ic] = 0;
        }
        
        switch(i)
        {
            case 0:
                sprintf(proposal[i]->name,"prior");
                proposal[i]->function = &draw_from_uniform_prior;
                proposal[i]->density  = &uniform_prior_density;
                proposal[i]->weight   = 0.1;
                proposal[i]->rjweight = 0.2;
                setup_prior_proposal(flags, prior, proposal[i]);
                check   += proposal[i]->weight;
                rjcheck += proposal[i]->rjweight;
                break;
            case 1:
                sprintf(proposal[i]->name,"fstat draw");
                setup_fstatistic_proposal(orbit, data, flags, proposal[i]);
                proposal[i]->function = &draw_from_fstatistic;
                proposal[i]->density  = &evaluate_fstatistic_proposal;
                proposal[i]->weight   = 0.2;
                proposal[i]->rjweight = 1.0; //that's a 1 all right.  don't panic
                check += proposal[i]->weight;
                break;
            case 2:
                sprintf(proposal[i]->name,"fstat jump");
                
                //re-use setup of fstat proposal from case 0
                proposal[i]->size   = proposal[1]->size;
                proposal[i]->norm   = proposal[1]->norm;
                proposal[i]->maxp   = proposal[1]->maxp;
                proposal[i]->vector = proposal[1]->vector;
                proposal[i]->matrix = proposal[1]->matrix;
                proposal[i]->tensor = proposal[1]->tensor;
                
                proposal[i]->function = &jump_from_fstatistic;
                proposal[i]->density  = &evaluate_fstatistic_proposal;
                proposal[i]->weight   = 0.2;
                proposal[i]->rjweight = 0.0;
                check   += proposal[i]->weight;
                rjcheck += proposal[i]->rjweight;
                break;
            case 3:
                sprintf(proposal[i]->name,"extrinsic prior");
                proposal[i]->function = &draw_from_extrinsic_prior;
                proposal[i]->density  = &prior_density;
                proposal[i]->weight   = 0.0;
                proposal[i]->rjweight = 0.0;
                check   += proposal[i]->weight;
                rjcheck += proposal[i]->rjweight;
                break;
            case 4:
                sprintf(proposal[i]->name,"fisher");
                proposal[i]->function = &draw_from_fisher;
                proposal[i]->density  = &symmetric_density;
                proposal[i]->weight = 1.0; //that's a 1 all right.  don't panic
                proposal[i]->rjweight = 0.0;
                //check   += proposal[i]->weight;
                rjcheck += proposal[i]->rjweight;
                break;
            case 5:
                sprintf(proposal[i]->name,"fm shift");
                proposal[i]->function = &fm_shift;
                proposal[i]->density  = &symmetric_density;
                proposal[i]->weight   = 0.1;
                proposal[i]->rjweight = 0.0;
                check   += proposal[i]->weight;
                rjcheck += proposal[i]->rjweight;
                break;
            case 6:
                sprintf(proposal[i]->name,"psi-phi jump");
                proposal[i]->function = &psi_phi_jump;
                proposal[i]->density  = &symmetric_density;
                proposal[i]->weight   = 0.1;
                proposal[i]->rjweight = 0.0;
                check   += proposal[i]->weight;
                rjcheck += proposal[i]->rjweight;
                break;
            case 7:
                sprintf(proposal[i]->name,"gmm draw");
                proposal[i]->function = &draw_from_gmm_prior;
                proposal[i]->density = &gmm_prior_density;
                proposal[i]->weight   = 0.0;
                proposal[i]->rjweight = 0.0;
                if(flags->catalog)
                {
                    setup_gmm_proposal(data, catalog, proposal[i]);
                    if(proposal[i]->Ngmm>0)
                    {
                        proposal[i]->weight   = 0.2;
                        proposal[i]->rjweight = 0.2;
                    }
                }
                check   += proposal[i]->weight;
                rjcheck += proposal[i]->rjweight;
                break;
            case 8:
                sprintf(proposal[i]->name,"cov draw");
                proposal[i]->function = &draw_from_cov;
                proposal[i]->density  = &cov_density;
                proposal[i]->weight  = 0.0;
                proposal[i]->rjweight = 0.0;
                if(flags->updateCov)
                {
                    setup_covariance_proposal(data, flags, proposal[i]);
                    proposal[i]->weight   = 0.2;
                    proposal[i]->rjweight = 0.2;
                }
                check   += proposal[i]->weight;
                rjcheck += proposal[i]->rjweight;
                break;
            default:
                break;
        }
    }
    //Fisher proposal fills in the cracks for fixed D moves
    proposal[4]->weight -= check;
    
    //Fstat proposal fills in the cracks for trans D moves
    proposal[1]->rjweight -= rjcheck;
    
    if(proposal[4]->weight<0.0 || proposal[1]->rjweight < 0.0)
    {
        fprintf(stderr,"Proposal weights not normalized (line %d of file %s)\n",__LINE__,__FILE__);
        exit(1);
    }
    
    if(!flags->quiet)
    {
        fprintf(stdout,"\n============ UCB Proposal Cocktail ============\n");
        fprintf(stdout,"   MCMC proposals:\n");
        for(int i=0; i<UCB_PROPOSAL_NPROP; i++)
        {
            if(proposal[i]->weight>0.0)fprintf(stdout,"     %i) %s %lg\n",i,proposal[i]->name,proposal[i]->weight);
        }
        fprintf(stdout,"   RJMCMC proposals:\n");
        for(int i=0; i<UCB_PROPOSAL_NPROP; i++)
        {
            if(proposal[i]->rjweight)fprintf(stdout,"     %i) %s %lg\n",i,proposal[i]->name,proposal[i]->rjweight);
        }
        fprintf(stdout,"===============================================\n");
    }
}

void initialize_vb_proposal(struct Orbit *orbit, struct Data *data, struct Prior *prior, struct Chain *chain, struct Flags *flags, struct Proposal **proposal, int NMAX)
{
    int NC = chain->NC;
    double check  =0.0;
    double rjcheck=0.0;
    
    for(int i=0; i<UCB_PROPOSAL_NPROP; i++)
    {
        proposal[i] = malloc(sizeof(struct Proposal));

        proposal[i]->trial  = malloc(NC*sizeof(int));
        proposal[i]->accept = malloc(NC*sizeof(int));
        
        for(int ic=0; ic<NC; ic++)
        {
            proposal[i]->trial[ic]  = 1;
            proposal[i]->accept[ic] = 0;
        }
        
        switch(i)
        {
            case 0:
                sprintf(proposal[i]->name,"prior");
                proposal[i]->function = &draw_from_uniform_prior;
                proposal[i]->density  = &prior_density;
                proposal[i]->weight   = 0.2;
                proposal[i]->rjweight = 0.0;
                setup_prior_proposal(flags, prior, proposal[i]);
                check   += proposal[i]->weight;
                rjcheck += proposal[i]->rjweight;
                break;
            case 1:
                sprintf(proposal[i]->name,"fstat draw");
                proposal[i]->function = &draw_from_fstatistic;
                proposal[i]->density  = &evaluate_fstatistic_proposal;
                proposal[i]->weight   = 0.0;
                proposal[i]->rjweight = 0.0; //that's a 1 all right.  don't panic
                check += proposal[i]->weight;
                break;
            case 2:
                sprintf(proposal[i]->name,"fstat jump");
                proposal[i]->function = &jump_from_fstatistic;
                proposal[i]->density  = &evaluate_fstatistic_proposal;
                proposal[i]->weight   = 0.0;
                proposal[i]->rjweight = 0.0;
                check   += proposal[i]->weight;
                rjcheck += proposal[i]->rjweight;
                break;
            case 3:
                sprintf(proposal[i]->name,"extrinsic prior");
                proposal[i]->function = &draw_from_extrinsic_prior;
                proposal[i]->density  = &prior_density;
                proposal[i]->weight   = 0.0;
                proposal[i]->rjweight = 0.0;
                check   += proposal[i]->weight;
                rjcheck += proposal[i]->rjweight;
                break;
            case 4:
                sprintf(proposal[i]->name,"fisher");
                proposal[i]->function = &draw_from_fisher;
                proposal[i]->density  = &symmetric_density;
                proposal[i]->weight = 1.0; //that's a 1 all right.  don't panic
                proposal[i]->rjweight = 0.0;
                //check   += proposal[i]->weight;
                rjcheck += proposal[i]->rjweight;
                break;
            case 5:
                sprintf(proposal[i]->name,"fm shift");
                proposal[i]->function = &fm_shift;
                proposal[i]->density  = &symmetric_density;
                proposal[i]->weight   = 0.0;
                proposal[i]->rjweight = 0.0;
                check   += proposal[i]->weight;
                rjcheck += proposal[i]->rjweight;
                break;
            case 6:
                sprintf(proposal[i]->name,"psi-phi jump");
                proposal[i]->function = &psi_phi_jump;
                proposal[i]->density  = &symmetric_density;
                proposal[i]->weight   = 0.2;
                proposal[i]->rjweight = 0.0;
                check   += proposal[i]->weight;
                rjcheck += proposal[i]->rjweight;
                break;
            case 7:
                sprintf(proposal[i]->name,"gmm draw");
                proposal[i]->function = &draw_from_gmm_prior;
                proposal[i]->density = &gmm_prior_density;
                proposal[i]->weight   = 0.0;
                proposal[i]->rjweight = 0.0;
                check   += proposal[i]->weight;
                rjcheck += proposal[i]->rjweight;
                break;
            case 8:
                sprintf(proposal[i]->name,"cov draw");
                proposal[i]->function = &draw_from_cov;
                proposal[i]->density  = &cov_density;
                proposal[i]->weight  = 0.0;
                proposal[i]->rjweight = 0.0;
                if(flags->updateCov)
                {
                    setup_covariance_proposal(data, flags, proposal[i]);
                    proposal[i]->weight   = 0.0;
                    proposal[i]->rjweight = 0.0;
                }
                check   += proposal[i]->weight;
                rjcheck += proposal[i]->rjweight;
                break;
            default:
                break;
        }
    }
    //Fisher proposal fills in the cracks for fixed D moves
    proposal[4]->weight -= check;
    
    //Fstat proposal fills in the cracks for trans D moves
    proposal[1]->rjweight -= rjcheck;
    
    if(proposal[4]->weight<0.0 || proposal[1]->rjweight < 0.0)
    {
        fprintf(stderr,"Proposal weights not normalized (line %d of file %s)\n",__LINE__,__FILE__);
        exit(1);
    }
    
    if(!flags->quiet)
    {
        fprintf(stdout,"\n============ VGB Proposal Cocktail ============\n");
        fprintf(stdout,"   MCMC proposals:\n");
        for(int i=0; i<UCB_PROPOSAL_NPROP; i++)
        {
            if(proposal[i]->weight>0.0)fprintf(stdout,"     %i) %s %lg\n",i,proposal[i]->name,proposal[i]->weight);
        }
        fprintf(stdout,"   RJMCMC proposals:\n");
        for(int i=0; i<UCB_PROPOSAL_NPROP; i++)
        {
            if(proposal[i]->rjweight)fprintf(stdout,"     %i) %s %lg\n",i,proposal[i]->name,proposal[i]->rjweight);
        }
        fprintf(stdout,"===============================================\n");
    }
}

void setup_fstatistic_proposal(struct Orbit *orbit, struct Data *data, struct Flags *flags, struct Proposal *proposal)
{
    /*
     -Step through 3D grid in frequency and sky location.
     - frequency resolution hard-coded to 1/4 of a bin
     - sky location resolution hard-coded to 30x30 bins
     -Compute F-statistic in each cell of the grid
     - cap F-statistic at SNRmax=20
     - normalize to make it a proper proposal (this part is a pain to get right...)
     */
    

    //grid sizes
    int n_f     = data->N/2;
    int n_theta = 20;
    int n_phi   = 20;
    
    if(flags->debug)
    {
        n_f/=4;
        n_theta/=3;
        n_phi/=3;
    }
    
    double d_f = (double)data->NFFT/(double)n_f/data->T;

    double d_theta = 2./(double)n_theta;
    double d_phi   = PI2/(double)n_phi;
    
    if(!flags->quiet)
    {
        fprintf(stdout,"\n============ F-statistic proposal ============\n");
        fprintf(stdout,"   n_f     = %i\n",n_f);
        fprintf(stdout,"   n_theta = %i\n",n_theta);
        fprintf(stdout,"   n_phi   = %i\n",n_phi);
        fprintf(stdout,"   cap     = %g\n",SNRCAP);
    }
            
    //allocate memory in proposal structure and pack up metadata
    /*
     proposal->matrix is 3x2 matrix.
     -rows are parameters {f,theta,phi}
     -columns are bin number and width {n,d}
     */
    proposal->matrix = malloc(3*sizeof(double *));
    for(int i=0; i<3; i++)
        proposal->matrix[i] = calloc(2,sizeof(double));
    
    proposal->matrix[0][0] = (double)n_f;
    proposal->matrix[0][1] = d_f;
    
    proposal->matrix[1][0] = (double)n_theta;
    proposal->matrix[1][1] = d_theta;
    
    proposal->matrix[2][0] = (double)n_phi;
    proposal->matrix[2][1] = d_phi;
    
    /*
     proposal->tensor holds the proposal density
     - n_f x n_theta x n_phi "tensor"
     */
    proposal->tensor = malloc(n_f*sizeof(double **));
    for(int i=0; i<n_f; i++)
    {
        proposal->tensor[i] = malloc(n_theta*sizeof(double *));
        for(int j=0; j<n_theta; j++)
        {
            proposal->tensor[i][j] = malloc(n_phi*sizeof(double));
            for(int k=0; k<n_phi; k++)
            {
                proposal->tensor[i][j][k] = 1.0;
            }
        }
    }
    
    double minp = +1e60;
    proposal->maxp = -1e60;
    
    /* compute or restore fisher-based proposal */
    char filename[MAXSTRINGSIZE];
    sprintf(filename,"%s/checkpoint/fstat_prop.bin",flags->runDir);
    FILE *propFile=NULL;
    int check=0;
    
    //first check if file exists
    if( (propFile=fopen(filename,"r")) )
    {
        fclose(propFile);
        check=1;
    }
    
    //if we are checkpointing and the files exist, use them
    if(flags->resume && check)
    {
        propFile=fopen(filename,"rb");
        
        fread(&proposal->norm, sizeof proposal->norm, 1, propFile);
        fread(&proposal->maxp, sizeof proposal->norm, 1, propFile);
        for(int i=0; i<n_f; i++)
        {
            for(int j=0; j<n_theta; j++)
            {
                for(int k=0; k<n_phi; k++)
                {
                    fread(&proposal->tensor[i][j][k],sizeof proposal->tensor[i][j][k], 1, propFile);
                }
            }
        }

        fclose(propFile);
    }
    
    //if no resume flag or no checkpoint file build the proposal
    else
    {
        build_fstatistic_proposal(orbit, data, flags, proposal);
        
        //store checkpointing files so we can skip this expensive step on resume
        propFile=fopen(filename,"wb");
        fwrite(&proposal->norm, sizeof proposal->norm, 1, propFile);
        fwrite(&proposal->maxp, sizeof proposal->norm, 1, propFile);
        for(int i=0; i<n_f; i++)
            for(int j=0; j<n_theta; j++)
                for(int k=0; k<n_phi; k++)
                    fwrite(&proposal->tensor[i][j][k],sizeof proposal->tensor[i][j][k], 1, propFile);
        fclose(propFile);


    }//end resume if/else

    //print diagnostics
    if(flags->verbose)
    {
        //get min
        for(int i=0; i<n_f; i++)
            for(int j=0; j<n_theta; j++)
                for(int k=0; k<n_phi; k++)
                    if(proposal->tensor[i][j][k]<minp) minp = proposal->tensor[i][j][k];

        char dirname[MAXSTRINGSIZE];
        sprintf(dirname,"%s/fstat",flags->runDir);
        mkdir(dirname,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        char filename[MAXSTRINGSIZE];
        
        for(int i=0; i<n_f; i++)
        {
            sprintf(filename,"%s/fstat/skymap_%05d.dat",flags->runDir,i);
            FILE *fptr=fopen(filename,"w");
            for (int j=0; j<n_theta; j++)
            {
                double theta = acos(-1. + (double)j*d_theta);
                
                for(int k=0; k<n_phi; k++)
                {
                    double phi = (double)k*d_phi;
                    
                    fprintf(fptr,"%.12g %.12g %.12g\n", cos(theta), phi, proposal->tensor[i][j][k]);
                }
                fprintf(fptr,"\n");
            }
            fclose(fptr);
        }
        write_Fstat_animation(data->fmin, data->T,proposal,flags->runDir, minp);
    }
    
    if(!flags->quiet)fprintf(stdout,"\n==============================================\n\n");
    fflush(stdout);
}

void build_fstatistic_proposal(struct Orbit *orbit, struct Data *data, struct Flags *flags, struct Proposal *proposal)
{
    //matrix to hold maximized extrinsic parameters
    double *Fparams = calloc(4,sizeof(double));

    int n_f     = (int)proposal->matrix[0][0];
    int n_theta = (int)proposal->matrix[1][0];
    int n_phi   = (int)proposal->matrix[2][0];
    
    double d_f     = proposal->matrix[0][1];
    double d_theta = proposal->matrix[1][1];
    double d_phi   = proposal->matrix[2][1];
    
    int n_fdot  = 10;
    double d_fdot = (1.1 - 0.1)/n_fdot; //scan over chirp mass range [0.1:1.1] Msolar

    double norm = 0.0;

    for(int i=0; i<n_f; i++)
    {
        if(n_f<100)
        {
            if(!flags->quiet) printProgress((double)i/(double)n_f);
        }
        else
        {
            if(i%(n_f/100)==0 && !flags->quiet)printProgress((double)i/(double)n_f);
        }
        
        double f;
        
        f = data->fmin + i*d_f;

        f += d_f/2; //evaluate logL in center of cell
        
        //loop over colatitude bins
        #pragma omp parallel for num_threads(flags->threads)
        for (int j=0; j<n_theta; j++)
        {
            //F-statistic for TDI variabls
            double logL_X,logL_AE, logLmax;

            //loop over longitude bins
            for(int k=0; k<n_phi; k++)
            {
                double theta = acos((-1. + (double)j*d_theta + d_theta/2)); //evaluate logL in center of cell
                double phi = (double)k*d_phi + d_phi/2; //evaluate logL in center of cell
                

                get_Fstat_logL(orbit, data, f, 0, theta, phi, &logL_X, &logL_AE, Fparams);

                logLmax = logL_AE;
                    
                if(logL_AE>1)
                {
                    for(int n=0; n<n_fdot; n++)
                    {
                        double fdot = ucb_fdot(0.1+n*d_fdot, f); //get fdot by scanning over chirpmass
                        
                        get_Fstat_logL(orbit, data, f, fdot, theta, phi, &logL_X, &logL_AE, Fparams);
                        
                        if(logL_AE>logLmax) logLmax = logL_AE;

                        if(logL_AE<1) break;
                    }
                }
                proposal->tensor[i][j][k] = logLmax*logLmax;

            }//end loop over longitude bins
        }//end loop over colatitude bins
    }//end loop over sub-bins

    //get normalization
    for(int i=0; i<n_f; i++)
        for(int j=0; j<n_theta; j++)
            for(int k=0; k<n_phi; k++)
                norm += proposal->tensor[i][j][k];
    
    //include correction for cell volume
    proposal->norm = (n_f*n_theta*n_phi)/norm;
    
    //normalize
    for(int i=0; i<n_f; i++)
        for(int j=0; j<n_theta; j++)
            for(int k=0; k<n_phi; k++)
                proposal->tensor[i][j][k] *= proposal->norm;

    //get max
    proposal->maxp = -1.e60;
    for(int i=0; i<n_f; i++)
        for(int j=0; j<n_theta; j++)
            for(int k=0; k<n_phi; k++)
                if(proposal->tensor[i][j][k]>proposal->maxp) proposal->maxp = proposal->tensor[i][j][k];
    
    free(Fparams);
}

void rebuild_fstatistic_proposal(struct Orbit *orbit, struct Data *data, struct Model *model, struct Flags *flags, struct Proposal *proposal)
{
    //store data
    double *Xsave = malloc(data->N*sizeof(double));
    double *Ysave = malloc(data->N*sizeof(double));
    double *Zsave = malloc(data->N*sizeof(double));
    double *Asave = malloc(data->N*sizeof(double));
    double *Esave = malloc(data->N*sizeof(double));

    memcpy(Xsave, data->tdi->X, data->N*sizeof(double));
    memcpy(Ysave, data->tdi->Y, data->N*sizeof(double));
    memcpy(Zsave, data->tdi->Z, data->N*sizeof(double));
    memcpy(Asave, data->tdi->A, data->N*sizeof(double));
    memcpy(Esave, data->tdi->E, data->N*sizeof(double));

    // send residual to fstat builder
    for(int n=0; n<data->N; n++)
    {
        if(data->Nchannel==2)
        {
            data->tdi->A[n] -= model->tdi->X[n];
            data->tdi->E[n] -= model->tdi->Y[n];
        }
        
        if(data->Nchannel==3)
        {
            data->tdi->X[n] -= model->tdi->X[n];
            data->tdi->Y[n] -= model->tdi->Y[n];
            data->tdi->Z[n] -= model->tdi->Z[n];
            
            //DFT F-stat uses A&E channels only
            XYZ2AE(data->tdi->X[n], data->tdi->Y[n], data->tdi->Z[n], &data->tdi->A[n], &data->tdi->E[n]);
        }
    }
    
    // convert DWT to DFT
    if(!strcmp(data->basis,"wavelet")) wavelet_layer_to_fourier_transform(data);
    
    build_fstatistic_proposal(orbit, data, flags, proposal);
    
    //restore data
    memcpy(data->tdi->X, Xsave, data->N*sizeof(double));
    memcpy(data->tdi->Y, Ysave, data->N*sizeof(double));
    memcpy(data->tdi->Z, Zsave, data->N*sizeof(double));
    memcpy(data->tdi->A, Asave, data->N*sizeof(double));
    memcpy(data->tdi->E, Esave, data->N*sizeof(double));

    free(Xsave);
    free(Ysave);
    free(Zsave);
    free(Asave);
    free(Esave);

}


void setup_gmm_proposal(struct Data *data, struct Catalog *catalog, struct Proposal *proposal)
{
    /* size of joint gmm */
    proposal->Ngmm = catalog->N;
    
    /* allocate space for gmm */
    proposal->gmm = malloc(proposal->Ngmm*sizeof(struct GMM*));
    
    /* point to different entries in catalog */
    for(int n=0; n<proposal->Ngmm; n++) proposal->gmm[n] = catalog->entry[n]->gmm;
}

void setup_prior_proposal(struct Flags *flags, struct Prior *prior, struct Proposal *proposal)
{
    /*
     awkwardly pack needed contents of prior structure
     into proposal for intput into standardized
     proposal density protocol
     */
    proposal->size         = flags->galaxyPrior;
    proposal->matrix       = malloc(2*sizeof(double *));
    proposal->matrix[0]    = calloc(3,sizeof(double));
    proposal->matrix[0][0] = prior->dcostheta;
    proposal->matrix[0][1] = prior->dphi;
    proposal->matrix[0][2] = (double)prior->nphi;
    proposal->matrix[1]    = prior->skyhist;
}


void setup_cdf_proposal(struct Data *data, struct Flags *flags, struct Proposal *proposal, int NMAX)
{
    if(!flags->quiet)fprintf(stdout,"\n============== Chain CDF proposal ==============\n\n");
    
    /*
     Use posterior samples from previous run to setup
     8x1D proposals from the marginalized posteriors
     */
    
    double junk;
    FILE *fptr=NULL;
    
    //parse chain file
    if(!flags->quiet)fprintf(stdout,"  reading chain file %s...\n",flags->cdfFile);
    fptr = fopen(flags->cdfFile,"r");
    proposal->size=0;
    while(!feof(fptr))
    {
        for(int j=0; j<UCB_MODEL_NP; j++)
        {
            int check = fscanf(fptr,"%lg",&junk);
            if(!check)
            {
                fprintf(stderr,"Error reading %s\n",flags->cdfFile);
                exit(1);
            }
        }
        proposal->size++;
    }
    rewind(fptr);
    proposal->size--;
    
    if(!flags->quiet)fprintf(stdout, "  samples in chain: %i\n",proposal->size);
    
    proposal->vector = calloc(proposal->size , sizeof(double));
    proposal->matrix = malloc(UCB_MODEL_NP * sizeof(double*));
    for(int j=0; j<UCB_MODEL_NP; j++) proposal->matrix[j] = calloc(proposal->size , sizeof(double));
    
    struct Model *temp = malloc(sizeof(struct Model));
    alloc_model(data,temp,NMAX);
    
    for(int n=0; n<proposal->size; n++)
    {
        scan_source_params(data, temp->source[0], fptr);
        for(int j=0; j<UCB_MODEL_NP; j++) proposal->matrix[j][n] = temp->source[0]->params[j];
    }
    free_model(temp);
    
    //now sort each row of the matrix
    for(int j=0; j<UCB_MODEL_NP; j++)
    {
        //fill up proposal vector
        for(int n=0; n<proposal->size; n++)
            proposal->vector[n] = proposal->matrix[j][n];
        
        //sort it
        double_sort(proposal->vector, proposal->size);
        
        //replace that row of the matrix
        for(int n=0; n<proposal->size; n++)
            proposal->matrix[j][n] = proposal->vector[n];
    }
    
    free(proposal->vector);
    
    fclose(fptr);
    
    if(!flags->quiet)fprintf(stdout,"\n================================================\n");
}

void test_covariance_proposal(struct Data *data, struct Flags *flags, struct Model *model, struct Prior *prior, struct Proposal *proposal, unsigned int *seed)
{
    fprintf(stdout,"\n======== Test Covariance matrix proposal =======\n");

    double logP;
    int N = 100000;

    //correction factor to rescale & renormalize covariance matrix
    double gamma = 0;// = (double)i/(double)N;
    double gamma2;// = gamma*gamma;
    
    int Nmodes = proposal->size*2;
    
    double **invCij = NULL;
    double **Lij  = NULL;
    double *mean = NULL;
    double ran_no[UCB_MODEL_NP];
    
    //renormalization asymptotically converges.  Do a few laps
    for(int j = 0; j<3; j++)
    {
        printf(" Iteration %i:{ ",j+1);
        for(int i=0; i<Nmodes; i++)
        {
            gamma = 0.0;
            
            //define some helper pointers for ease of reading
            mean = proposal->matrix[i];
            Lij  = proposal->tensor[i];
            
            //scale NP-dimensional jump to 1-sigma of the joint distribution
            double scale = 1.;
            
            for(int n=0; n<N; n++)
            {
                
                //get vector of gaussian draws n;  y_i = x_mean_i + sum_j Lij^-1 * n_j
                for(int m=0; m<UCB_MODEL_NP; m++)
                {
                    ran_no[m] = rand_r_N_0_1(seed)*scale;
                }
                
                //the matrix multiplication...
                for(int n=0; n<UCB_MODEL_NP; n++)
                {
                    //start at mean
                    model->source[0]->params[n] = mean[n];
                    
                    //add contribution from each row of invC
                    for(int k=0; k<UCB_MODEL_NP; k++) model->source[0]->params[n] += ran_no[k]*Lij[n][k];
                }
                
                
                logP = evaluate_prior(flags, data, model, prior, model->source[0]->params);
                if(logP>-INFINITY) gamma+=1.;
            }
            
            
            gamma  = gamma/(double)N;
            gamma2 = gamma*gamma;
            
            printf("%.2f ",gamma);
            
            proposal->vector[i*2+1]/=gamma;
            invCij = proposal->tensor[Nmodes+i];
            Lij = proposal->tensor[i];
            
            for(int n=0; n<UCB_MODEL_NP; n++)
            {
                for(int m=0; m<UCB_MODEL_NP; m++)
                {
                    invCij[n][m] /= gamma2;
                    Lij[n][m]    *= gamma;
                }
            }
        }
        printf("}\n");
    }
    fprintf(stdout,"================================================\n\n");
}

void setup_covariance_proposal(struct Data *data, struct Flags *flags, struct Proposal *proposal)
{
    if(!flags->quiet)fprintf(stdout,"\n========== Covariance matrix proposal ==========\n\n");
    
    
    if(!flags->quiet)fprintf(stdout,"   reading covariance matrix %s...\n",flags->covFile);
    FILE *fptr;
    int check;
    fptr = fopen(flags->covFile,"r");
    int no_sources;
    double junk;
    check = fscanf(fptr, "%d%lg%lg%lg%lg%lg%lg%lg", &no_sources, &junk, &junk, &junk, &junk, &junk, &junk , &junk);
    if(!check)
    {
        fprintf(stderr,"Error reading %s\n",flags->covFile);
        exit(1);
    }
    
    proposal->size = no_sources;
    if(!flags->quiet)fprintf(stdout,"\nproposal->size=%d\n",proposal->size);
    int Ncov=proposal->size*2;
    double alpha, detCij;
    
    
    //book keeping (weights and determinants
    proposal->vector = calloc(Ncov*2,sizeof(double));
    
    //centroids
    proposal->matrix = malloc(Ncov*sizeof(double*));
    for(int j=0; j<Ncov; j++) proposal->matrix[j] = calloc(UCB_MODEL_NP , sizeof(double));
    
    //covariance matrix and friends
    proposal->tensor = malloc((Ncov*2)*sizeof(double**));
    for(int j=0; j<Ncov*2; j++)
    {
        proposal->tensor[j] = malloc(UCB_MODEL_NP * sizeof(double*));
        for(int k=0; k<UCB_MODEL_NP; k++)
        {
            proposal->tensor[j][k] = calloc(UCB_MODEL_NP , sizeof(double));
        }
    }
    
    //create some aliases to make the code more readable
    double *mean    = NULL; //centroids of multivariate
    double **Cij    = NULL; //covariance matrix
    double **Lij    = NULL; //lower triangle of cholesky decomp. of Cij
    
    
    //repeat for every covariance matrix in file
    for(int n=0; n<Ncov; n++)
    {
        
        //second row has the weight, determinant, and then a bunch of zeroes
        check = fscanf(fptr, "%lg%lg%lg%lg%lg%lg%lg%lg", &alpha, &detCij, &junk, &junk, &junk, &junk, &junk, &junk);
        if(!check)
        {
            fprintf(stderr,"Error reading %s\n",flags->covFile);
            exit(1);
        }
        //relative weighting of mode
        proposal->vector[n*2]=alpha;
        
        //absorb normalization constants into stored determinant detCij
        //    proposal->vector[n*Ncov+1]=1./(pow(PI2,0.5*NP)*sqrt(detCij));
        proposal->vector[n*2+1]=1./(pow(PI2,0.5*UCB_MODEL_NP)*sqrt(detCij));
        
        //second row has the centroids
        mean = proposal->matrix[n];
        for(int i=0; i<UCB_MODEL_NP; i++)
        {
            check = fscanf(fptr, "%lg", &mean[i]);
            if(!check)
            {
                fprintf(stderr,"Error reading %s\n",flags->covFile);
                exit(1);
            }
        }
        
        //next NP rows have the covariance matrix
        Cij = proposal->tensor[Ncov+n];
        check=0.0;
        for(int i=0; i<UCB_MODEL_NP; i++) for(int j=0; j<UCB_MODEL_NP; j++) check+=fscanf(fptr, "%lg", &Cij[i][j]);
        if(!check)
        {
            fprintf(stderr,"Error reading %s\n",flags->covFile);
            exit(1);
        }
        
        //use gsl cholesky decomposition, preserving Cij
        /*
         this gives identical results to what is in the cov files
         */
        //Lij = proposal->tensor[n];
        //cholesky_decomp(Cij, Lij, NP);
        
        //next NP rows are the lower half of the cholesky decomp.
        Lij = proposal->tensor[n];
        check=0;
        for(int i=0; i<UCB_MODEL_NP; i++) for(int j=0; j<UCB_MODEL_NP; j++) check+=fscanf(fptr, "%lg", &Lij[i][j]);
        if(!check)
        {
            fprintf(stderr,"Error reading %s\n",flags->covFile);
            exit(1);
        }

        //get inverse of Cij (invert_matrix writes over contents of input)
        invert_matrix(Cij,UCB_MODEL_NP);
        
    }
    
    fclose(fptr);
    
    if(!flags->quiet)fprintf(stdout,"\n================================================\n");
}

double draw_from_fstatistic(struct Data *data, UNUSED struct Model *model, UNUSED struct Source *source, struct Proposal *proposal, double *params, unsigned int *seed)
{
    double logP = 0.0;
    
    int n_f     = (int)proposal->matrix[0][0];
    int n_theta = (int)proposal->matrix[1][0];
    int n_phi   = (int)proposal->matrix[2][0];
    
    double d_f     = proposal->matrix[0][1];
    double d_theta = proposal->matrix[1][1];
    double d_phi   = proposal->matrix[2][1];
    
    double q,costheta,phi,p,alpha;
    
    double i,j,k;
    
    //first draw from prior
    draw_from_uniform_prior(data, model, source, proposal, params, seed);
    
    //now rejection sample on f,theta,phi
    do
    {
        i = rand_r_U_0_1(seed)*n_f;
        j = rand_r_U_0_1(seed)*n_theta;
        k = rand_r_U_0_1(seed)*n_phi;
        
        q = (double)(data->qmin) + i*d_f*data->T;

        costheta = -1. + j*d_theta;
        phi      = k*d_phi;
                
        p = proposal->tensor[(int)i][(int)j][(int)k];
        alpha = rand_r_U_0_1(seed)*proposal->maxp;
                
    }while(p<alpha);
    
    params[0] = q;
    params[1] = costheta;
    params[2] = phi;
    
    logP = evaluate_fstatistic_proposal(data, model, source, proposal, params);
    
    return logP;
}

double jump_from_fstatistic(struct Data *data, struct Model *model, struct Source *source, struct Proposal *proposal, double *params, unsigned int *seed)
{
    double logP = 0.0;
    
    int n_f     = (int)proposal->matrix[0][0];
    int n_theta = (int)proposal->matrix[1][0];
    int n_phi   = (int)proposal->matrix[2][0];
    
    double d_f     = proposal->matrix[0][1];
    double d_theta = proposal->matrix[1][1];
    double d_phi   = proposal->matrix[2][1];
    
    double q,costheta,phi,p,alpha;
    
    double i=0,j,k;
    
    /* half the time do an fm shift, half the time completely reboot frequency */
    int fmFlag = 0;
    if(rand_r_U_0_1(seed)<-0.5) fmFlag=1;
    
    if(fmFlag)
    {
        fm_shift(data, model, source, proposal, params, seed);
        
        q = params[0];
        i = (int)floor((q - data->qmin)/(d_f*data->T));
        
        if(i<0.0 || i>n_f-1) return -INFINITY;
    }
    
    //now rejection sample on f,theta,phi
    do
    {
        if(!fmFlag)i = rand_r_U_0_1(seed)*n_f;
        j = rand_r_U_0_1(seed)*n_theta;
        k = rand_r_U_0_1(seed)*n_phi;
        
        q = (double)(data->qmin) + i*d_f*data->T;

        costheta = -1. + j*d_theta;
        phi      = k*d_phi;
        
        if(
           ((int)i < 0 || (int)i > n_f-1)     ||
           ((int)j < 0 || (int)j > n_theta-1) ||
           ((int)k < 0 || (int)k > n_phi-1)
           )
        {
            fprintf(stdout,"%i (%i) %i (%i) %i (%i)\n",(int)i,n_f,(int)j,n_theta,(int)k,n_phi);
            fflush(stdout);
            return -INFINITY;
        }
        else
            p = proposal->tensor[(int)i][(int)j][(int)k];
        alpha = rand_r_U_0_1(seed)*proposal->maxp;
                
    }while(p<alpha);
    
    params[0] = q;
    params[1] = costheta;
    params[2] = phi;
    
    logP = evaluate_fstatistic_proposal(data, model, source, proposal, params);
    
    return logP;
}


//double evaluate_fstatistic_proposal(struct Data *data, struct Proposal *proposal, double *params)
double evaluate_fstatistic_proposal(struct Data *data, UNUSED struct Model *model, UNUSED struct Source * source, struct Proposal *proposal, double *params)
{
    double logP = 0.0;
    
    //TODO: sky location proposal does not support galaxy prior!
    double *skyhist = NULL; //dummy pointer for sky location prior
    
    //sky location prior
    logP += evalaute_sky_location_prior(params, model->prior, model->logPriorVolume, 0, skyhist, 1, 1, 1);
    
    //everything else uses simple uniform priors
    logP += evaluate_uniform_priors(params, model->prior, model->logPriorVolume);
    
    double d_f     = proposal->matrix[0][1];
    double d_theta = proposal->matrix[1][1];
    double d_phi   = proposal->matrix[2][1];
    
    int n_f     = (int)proposal->matrix[0][0];
    int n_theta = (int)proposal->matrix[1][0];
    int n_phi   = (int)proposal->matrix[2][0];
    
    int i;
    i = (int)floor((params[0] - data->qmin)/(d_f*data->T));

    int j = (int)floor((params[1] - -1)/d_theta);
    int k = (int)floor((params[2])/d_phi);
    
    if      (i<0 || i>=n_f    ) return -INFINITY;
    else if (j<0 || j>=n_theta) return -INFINITY;
    else if (k<0 || k>=n_phi  ) return -INFINITY;
    else logP += log(proposal->tensor[i][j][k]);
    
    return logP;
}

double prior_density(struct Data *data, struct Model *model, UNUSED struct Source *source, struct Proposal *proposal, double *params)
{
    double logP = 0.0;
    double **uniform_prior = model->prior;
    
    //guard against nan's, but do so loudly
    for(int i=0; i<UCB_MODEL_NP; i++)
    {
        if(params[i]!=params[i])
        {
            fprintf(stderr,"WARNING: parameter %i not a number (line %d of file %s)\n",i,__LINE__,__FILE__);
            return -INFINITY;
        }
    }
    
    //unpack (obtuse) proposal structure
    int galaxyPriorFlag =      proposal->size;
    double dcostheta    =      proposal->matrix[0][0];
    double dphi         =      proposal->matrix[0][1];
    int nphi            = (int)proposal->matrix[0][2];
    double *skyhist     =      proposal->matrix[1];
    
    //sky location prior
    logP += evalaute_sky_location_prior(params, uniform_prior, model->logPriorVolume, galaxyPriorFlag , skyhist, dcostheta, dphi, nphi);
        
    //everything else uses simple uniform priors
    logP += evaluate_uniform_priors(params, uniform_prior, model->logPriorVolume);
    
    
    return logP;
    
}

double symmetric_density(UNUSED struct Data *data, UNUSED struct Model *model, UNUSED struct Source *source, UNUSED struct Proposal *proposal, UNUSED double *params)
{
    return 0.0;
}

