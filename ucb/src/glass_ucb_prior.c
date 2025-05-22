/*
 * Copyright 2019 Tyson B. Littenberg & Neil J. Cornish
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <glass_utils.h>

#include "glass_ucb_model.h"
#include "glass_ucb_io.h"
#include "glass_ucb_waveform.h"
#include "glass_ucb_catalog.h"
#include "glass_ucb_prior.h"

static double loglike(double *x, int D)
{
    return log(galaxy_distribution(x,GALAXY_A,GALAXY_Rb, GALAXY_Rd, GALAXY_Zd));    
}

void set_galaxy_prior(struct Flags *flags, struct Prior *prior)
{
    if(!flags->quiet)
    {
        fprintf(stdout,"\n============ Galaxy model sky prior ============\n");
        fprintf(stdout,"Monte carlo over galaxy model\n");
        fprintf(stdout,"   Distance to GC = %g kpc\n",GALAXY_RGC);
        fprintf(stdout,"   Disk Radius    = %g kpc\n",GALAXY_Rd);
        fprintf(stdout,"   Disk Height    = %g kpc\n",GALAXY_Zd);
        fprintf(stdout,"   Bulge Radius   = %g kpc\n",GALAXY_Rb);
        fprintf(stdout,"   Bulge Fraction = %g\n",    GALAXY_A);
    }
    double *x, *y;  // current and proposed parameters
    int D = 3;  // number of parameters
    int Nth = 200;  // bins in cos theta
    int Nph = 200;  // bins in phi
    int MCMC=100000000;
    int j;
    int ith, iph, cnt;
    double H, dOmega;
    double logLx, logLy;
    double alpha, beta, xx, yy, zz;
    double *xe, *xg;
    double r_ec, theta, phi;
    int mc;
    //  FILE *chain = NULL;
    
    if(flags->debug)
    {
        Nth /= 10;
        Nph /= 10;
        MCMC/=10;
    }
    
    unsigned int r = 150914;
    
    x =  (double*)calloc(D,sizeof(double));
    xe = (double*)calloc(D,sizeof(double));
    xg = (double*)calloc(D,sizeof(double));
    y =  (double*)calloc(D,sizeof(double));
    
    prior->skyhist = (double*)calloc((Nth*Nph),sizeof(double));
    
    prior->dcostheta = 2./(double)Nth;
    prior->dphi      = 2.*M_PI/(double)Nph;
    
    prior->ncostheta = Nth;
    prior->nphi      = Nph;
    
    prior->skymaxp = 0.0;
    
    for(ith=0; ith< Nth; ith++)
    {
        for(iph=0; iph< Nph; iph++)
        {
            prior->skyhist[ith*Nph+iph] = 0.0;
        }
    }
    
    // starting values for parameters
    x[0] = 0.5;
    x[1] = 0.4;
    x[2] = 0.1;
    
    logLx = loglike(x, D);
    
    cnt = 0;
    
    for(mc=0; mc<MCMC; mc++)
    {
        if(mc%(MCMC/100)==0)printProgress((double)mc/(double)MCMC);
        
        alpha = rand_r_U_0_1(&r);
        
        if(alpha > 0.7)  // uniform draw from a big box
        {
            
            y[0] = 20.0*GALAXY_Rd*(-1.0+2.0*rand_r_U_0_1(&r));
            y[1] = 20.0*GALAXY_Rd*(-1.0+2.0*rand_r_U_0_1(&r));
            y[2] = 40.0*GALAXY_Zd*(-1.0+2.0*rand_r_U_0_1(&r));
            
        }
        else
        {
            
            beta = rand_r_U_0_1(&r);
            
            if(beta > 0.8)
            {
                xx = 0.1;
            }
            else if (beta > 0.4)
            {
                xx = 0.01;
            }
            else
            {
                xx = 0.001;
            }
            for(j=0; j< 2; j++) y[j] = x[j] + rand_r_N_0_1(&r)*xx;
            y[2] = x[2] + 0.1*rand_r_N_0_1(&r)*xx;
            
        }
        
        
        logLy = loglike(y, D);
        
        H = logLy - logLx;
        beta = rand_r_U_0_1(&r);
        beta = log(beta);
        
        if(H > beta)
        {
            for(j=0; j< D; j++) x[j] = y[j];
            logLx = logLy;
        }
        
        if(mc%100 == 0)
        {
            
            /* rotate from galactic to ecliptic coordinates */
            
            xg[0] = x[0] - GALAXY_RGC;   // solar barycenter is offset from galactic center along x-axis (by convention)
            xg[1] = x[1];
            xg[2] = x[2];
            
            rotate_galtoeclip(xg, xe);
            
            r_ec = sqrt(xe[0]*xe[0]+xe[1]*xe[1]+xe[2]*xe[2]);
            
            theta = M_PI/2.0-acos(xe[2]/r_ec);
            
            phi = atan2(xe[1],xe[0]);
            
            if(phi<0.0) phi += 2.0*M_PI;
            
            ith = (int)(0.5*(1.0+sin(theta))*(double)(Nth));
            iph = (int)((2*M_PI-phi)/(2.0*M_PI)*(double)(Nph));
            
            cnt++;
            
            if(ith < 0 || ith > Nth -1) printf("%d %d\n", ith, iph);
            if(iph < 0 || iph > Nph -1) printf("%d %d\n", ith, iph);
            
            prior->skyhist[ith*Nph+iph] += 1.0;
            
        }
        
    }
    
    dOmega = 4.0*M_PI/(double)(Nth*Nph);
    
    double uni = 0.1;
    //fprintf(stderr,"\n   HACK:  setup_galaxy_prior() uni=%g\n",uni);
    yy = (1.0-uni)/(double)(cnt);
    zz = uni/(double)(Nth*Nph);
    
    yy /= dOmega;
    zz /= dOmega;
    
    for(ith=0; ith< Nth; ith++)
    {
        for(iph=0; iph< Nph; iph++)
        {
            xx = yy*prior->skyhist[ith*Nph+iph];
            prior->skyhist[ith*Nph+iph] = log(xx + zz);
            if(prior->skyhist[ith*Nph+iph]>prior->skymaxp) prior->skymaxp = prior->skyhist[ith*Nph+iph];
        }
    }
    
    if(flags->verbose)
    {
        FILE *fptr = fopen("skyprior.dat", "w");
        for(ith=0; ith< Nth; ith++)
        {
            xx = -1.0+2.0*((double)(ith)+0.5)/(double)(Nth);
            for(iph=0; iph< Nph; iph++)
            {
                yy = 2.0*M_PI*((double)(iph)+0.5)/(double)(Nph);
                fprintf(fptr,"%e %e %e\n", yy, xx, prior->skyhist[ith*Nph+iph]);
            }
            fprintf(fptr,"\n");
        }
        fclose(fptr);
    }
    
    free(x);
    free(y);
    free(xe);
    free(xg);
    if(!flags->quiet)fprintf(stdout,"\n================================================\n\n");
    fflush(stdout);
    
}

void set_uniform_prior(struct Flags *flags, struct Model *model, struct Data *data, int verbose)
{
    /*
     params[0] = source->f0*T;
     params[1] = source->costheta;
     params[2] = source->phi;
     params[3] = log(source->amp);
     params[4] = source->cosi;
     params[5] = source->psi;
     params[6] = source->phi0;
     params[7] = source->dfdt*T*T;
     */
    
    
    model->t0 = data->t0; //TODO move this to model setup
    
    //TODO: assign priors by parameter name, use mapper to get into vector (more robust to changes)
    
    //frequency bin
    model->prior[0][0] = data->fmin*data->T;
    model->prior[0][1] = data->fmax*data->T;
    
    //colatitude
    model->prior[1][0] = -1.0;
    model->prior[1][1] =  1.0;
    
    //longitude
    model->prior[2][0] = 0.0;
    model->prior[2][1] = PI2;
    
    //log amplitude
    model->prior[3][0] = -54.0;
    model->prior[3][1] = -46.0;
    
    //cos inclination
    model->prior[4][0] = -1.0;
    model->prior[4][1] =  1.0;
    
    //polarization
    model->prior[5][0] = 0.0;
    model->prior[5][1] = M_PI;
    
    //phase
    model->prior[6][0] = 0.0;
    model->prior[6][1] = PI2;
    
    //fdot (bins/Tobs)
    
    /* frequency derivative priors are a little trickier...*/
    double fmin = model->prior[0][0]/data->T;
    double fmax = model->prior[0][1]/data->T;
    
    /* emprical envelope functions from Gijs' MLDC catalog */
    double fdotmin = -0.000005*pow(fmin,(13./3.));
    double fdotmax = 0.0000008*pow(fmax,(11./3.));
        
    /* use prior on chirp mass to convert to priors on frequency evolution */
    if(flags->detached)
    {
        double Mcmin = 0.15;
        double Mcmax = 1.00;
        
        fdotmin = ucb_fdot(Mcmin, fmin);
        fdotmax = ucb_fdot(Mcmax, fmax);
    }
    
    double fddotmin = 11.0/3.0*fdotmin*fdotmin/fmax;
    double fddotmax = 11.0/3.0*fdotmax*fdotmax/fmin;
    
    if(!flags->detached)
    {
        fddotmin = -fddotmax;
    }
    if(verbose && !flags->quiet)
    {
        fprintf(stdout,"\n============== PRIORS ==============\n");
        if(flags->detached)fprintf(stdout,"  Assuming detached binary, Mchirp = [0.15,1]\n");
        fprintf(stdout,"  p(f)     = U[%g,%g]\n",fmin,fmax);
        fprintf(stdout,"  p(fdot)  = U[%g,%g]\n",fdotmin,fdotmax);
        fprintf(stdout,"  p(fddot) = U[%g,%g]\n",fddotmin,fddotmax);
        fprintf(stdout,"  p(lnA)   = U[%g,%g]\n",model->prior[3][0],model->prior[3][1]);
        fprintf(stdout,"====================================\n\n");
    }
    
    if(UCB_MODEL_NP>7)
    {
        model->prior[7][0] = fdotmin*data->T*data->T;
        model->prior[7][1] = fdotmax*data->T*data->T;
    }
    if(UCB_MODEL_NP>8)
    {
        model->prior[8][0] = fddotmin*data->T*data->T*data->T;
        model->prior[8][1] = fddotmax*data->T*data->T*data->T;
    }
    
    /*
     
     Learn to parse prior files with flexible format and
     reset uniform priors accordingly (doing the naive thing
     of setting a uniform distribution over the ~90% credible
     interval for example.  Will need to expand to at least
     gaussian priors.
     
     What to do about intervals so small that they are
     effectively delta functions?  Might anger some proposals...
     
     GAIA distance accuracy < 20%
     
     */
    if(flags->emPrior)
    {
        FILE *priorFile = fopen(flags->pdfFile,"r");
        char name[16];
        double min,max;
        
        while(fscanf(priorFile,"%s %lg %lg",name,&min,&max) != EOF)
        {
            
            if(strcmp("f0",name) == 0)
            {
                model->prior[0][0] = min*data->T;
                model->prior[0][1] = max*data->T;
            }
            
            else if(strcmp("dfdt",name) == 0)
            {
                model->prior[7][0] = min*data->T*data->T;
                model->prior[7][1] = max*data->T*data->T;
            }
            
            else if(strcmp("amp",name) == 0)
            {
                model->prior[3][0] = log(min);
                model->prior[3][1] = log(max);
            }
            
            else if(strcmp("phi",name) == 0)
            {
                model->prior[2][0] = min;
                model->prior[2][1] = max;
            }
            
            else if(strcmp("costheta",name) == 0)
            {
                model->prior[1][0] = min;
                model->prior[1][1] = max;
            }
            
            else if(strcmp("cosi",name) == 0)
            {
                model->prior[4][0] = min;
                model->prior[4][1] = max;
            }
            
            else if(strcmp("psi",name) == 0)
            {
                model->prior[5][0] = min;
                model->prior[5][1] = max;
            }
            
            else if(strcmp("phi0",name) == 0)
            {
                model->prior[6][0] = min;
                model->prior[6][1] = max;
            }
            
            else
            {
                fprintf(stdout, "unrecognized parameter in prior file: %s\n",name);
            }      
        }
        
        fclose(priorFile);
    }
    
    //set prior volume
    for(int n=0; n<UCB_MODEL_NP; n++) model->logPriorVolume[n] = log(model->prior[n][1]-model->prior[n][0]);
    
}

int check_range(double *params, double **uniform_prior)
{
    //nan check
    for(int n=0; n<UCB_MODEL_NP; n++) if(params[n]!=params[n]) return 1;
    
    //frequency bin (uniform)
    if(params[0]<uniform_prior[0][0] || params[0]>uniform_prior[0][1]) return 1;
    
    //cosine co-latitude
    if(params[1]<uniform_prior[1][0] || params[1]>uniform_prior[1][1]) return 1;

    //longitude
    if(params[2]<uniform_prior[2][0] || params[2]>=uniform_prior[2][1])
    {
        params[2] = atan2(sin(params[2]),cos(params[2]));
        if(params[2] < 0.0) params[2] += PI2;
    }

    //amplitude
    if(params[3]<uniform_prior[3][0] || params[3]>uniform_prior[3][1]) return 1;

    //cosine inclination
    if(params[4]<uniform_prior[4][0] || params[4]>uniform_prior[4][1]) return 1;
    
    //polarization
    while(params[5]<uniform_prior[5][0]) params[5] += M_PI;
    while(params[5]>uniform_prior[5][1]) params[5] -= M_PI;

    //phase
    while(params[6]<uniform_prior[6][0]) params[6] += PI2;
    while(params[6]>uniform_prior[6][1]) params[6] -= PI2;
    
    //fdot (bins/Tobs)
    if(UCB_MODEL_NP>7) if(params[7]<uniform_prior[7][0] || params[7]>uniform_prior[7][1]) return 1;
    
    //fddot
    if(UCB_MODEL_NP>8) if(params[8]<uniform_prior[8][0] || params[8]>uniform_prior[8][1]) return 1;

    return 0;
}

void set_gmm_prior(struct Flags *flags, struct Data *data, struct Prior *prior, struct Catalog *catalog)
{
    //get size of full catalog
    int N = catalog->N;
    
    //allocate gmm to include the full catalog
    prior->gmm = malloc(sizeof(struct GMM));
    prior->gmm->NMODE = 0;
    for(size_t n=0; n<N; n++) prior->gmm->NMODE += catalog->entry[n]->gmm->NMODE;
    prior->gmm->modes = malloc(prior->gmm->NMODE*sizeof(struct MVG *));
    
    for(size_t n=0; n<prior->gmm->NMODE; n++)
    {
        prior->gmm->modes[n] = malloc(sizeof(struct MVG));
        
        alloc_MVG(prior->gmm->modes[n], UCB_MODEL_NP);
    }

    //combine modes into one GMM
    size_t m=0;
    for(size_t n=0; n<N; n++)
    {
        for(size_t i=0; i<catalog->entry[n]->gmm->NMODE; i++)
        {
            copy_MVG(catalog->entry[n]->gmm->modes[i],prior->gmm->modes[m]);

            //(clumsily) renormalize modes
            prior->gmm->modes[m]->p /= (double)N;
            
            m++;
        }
    }
    
    //prior->gmm = catalog->entry[0]->gmm;
}

double evaluate_gmm_prior(struct Data *data, struct GMM *gmm, double *params)
{
    double *x = double_vector(UCB_MODEL_NP);
    
    /* pointers to GMM contents */
    struct MVG **modes = gmm->modes;
    size_t NMODES = gmm->NMODE;
    

    //pack parameters into source with correct units
    struct Source *source = malloc(sizeof(struct Source));
    alloc_source(source, data->N, data->Nchannel);
    
    map_array_to_params(source, params, data->T);
    x[0] = source->f0;
    x[1] = source->costheta;
    x[2] = source->phi;
    x[3] = log(source->amp);
    x[4] = source->cosi;
    x[5] = source->psi;
    x[6] = source->phi0;
    if(UCB_MODEL_NP>7)
        x[7] = source->dfdt;
    if(UCB_MODEL_NP>8)
        x[8] = source->d2fdt2;

    //map parameters to R
    double xmin,xmax,xn,yn, logJ = 0;
    for(size_t n=0; n<UCB_MODEL_NP; n++)
    {
        xmin = modes[0]->minmax[n][0];
        xmax = modes[0]->minmax[n][1];
        xn = x[n];
        if(xn < xmin || xn >= xmax)
        {            
            //clean up
            free_double_vector(x);
            free_source(source);
            return -INFINITY;
        }
        yn = logit(xn,xmin,xmax);
        x[n] = yn;
        
        //Jacobian
        //logJ -= log(dsigmoid(yn, xmin, xmax));
    }
    
    //sum over modes
    double P=0.0;
    for(size_t k=0; k<NMODES; k++)
        P += modes[k]->p*multivariate_gaussian(x,modes[k],UCB_MODEL_NP);
    
    //clean up
    free_double_vector(x);
    free_source(source);
    
    return log(P) + logJ;
}

double evaluate_prior(struct Flags *flags, struct Data *data, struct Model *model, struct Prior *prior, double *params)
{
    double logP=0.0;
    double **uniform_prior = model->prior;
    
    //guard against nan's, but do so loudly
    if(check_range(params, uniform_prior)) return -INFINITY;

    //update from existing runs prior
    if(flags->update) logP = evaluate_gmm_prior(data, prior->gmm, params);
    
    //blind search prior
    else
    {
        //sky location prior
        logP += evalaute_sky_location_prior(params, uniform_prior, model->logPriorVolume, flags->galaxyPrior, prior->skyhist, prior->dcostheta, prior->dphi, prior->nphi);
                
        //everything else uses simple uniform priors
        logP += evaluate_uniform_priors(params, uniform_prior, model->logPriorVolume);
    }
    
    return logP;
}

double evaluate_uniform_priors(double *params, double **uniform_prior, double *logPriorVolume)
{
    if(check_range(params, uniform_prior)) return -INFINITY;
    
    double logP = 0.0;
    //frequency bin (uniform)
    //TODO: is frequency logPriorVolume up to date?
    logP -= log(uniform_prior[0][1]-uniform_prior[0][0]);

    //amplitude
    logP -= logPriorVolume[3];

    //cosine inclination
    logP -= logPriorVolume[4];
    
    //polarization
    logP -= logPriorVolume[5];
    
    //phase
    logP -= logPriorVolume[6];
    
    //fdot (bins/Tobs)
    if(UCB_MODEL_NP>7) logP -= logPriorVolume[7];
    
    //fddot
    if(UCB_MODEL_NP>8) logP -= logPriorVolume[8];
    
    return logP;
}

double evalaute_sky_location_prior(double *params, double **uniform_prior, double *logPriorVolume, int galaxyFlag, double *skyhist, double dcostheta, double dphi, int nphi)
{
    
    double logP = 0.0;
    if(galaxyFlag)
    {
        if(params[1]<uniform_prior[1][0] || params[1]>uniform_prior[1][1]) return -INFINITY;
        
        if(uniform_prior[2][0] > 0.0 || uniform_prior[2][1] < PI2)
        {
            //rejection sample on reduced prior range
            if(params[2]<uniform_prior[2][0] || params[2]>uniform_prior[2][1]) return -INFINITY;
        }
        else
        {
            //periodic boundary conditions for full range
            while(params[2] < 0  ) params[2] += PI2;
            while(params[2] > PI2) params[2] -= PI2;
        }
        //map costheta and phi to index of skyhist array
        int i = (int)floor((params[1]-uniform_prior[1][0])/dcostheta);
        int j = (int)floor((params[2]-uniform_prior[2][0])/dphi);
        
        int k = i*nphi + j;
        
        logP += skyhist[k];
        
        //    FILE *fptr = fopen("prior.dat","a");
        //    fprintf(fptr,"%i %i %i %g\n",i,j,k,prior->skyhist[k]);
        //    fclose(fptr);
    }
    else
    {
        //colatitude (reflective)
        if(params[1]<uniform_prior[1][0] || params[1]>uniform_prior[1][1]) return -INFINITY;
        else logP -= logPriorVolume[1];
        
        //longitude (periodic)
        if(uniform_prior[2][0] > 0.0 || uniform_prior[2][1] < PI2)
        {
            //rejection sample on reduced prior range
            if(params[2]<uniform_prior[2][0] || params[2]>uniform_prior[2][1]) return -INFINITY;
        }
        else
        {
            if(params[2]<uniform_prior[2][0] || params[2]>=uniform_prior[2][1])
            {
                params[2] = atan2(sin(params[2]),cos(params[2]));
                if(params[2] < 0.0) params[2] += PI2;
            }
        }
        logP -= logPriorVolume[2];
    }
    return logP;
}


double evaluate_calibration_prior(struct Data *data, struct Model *model)
{
    
    double logP = 0.0;
    
    /* apply calibration error to full signal model */
    double dA,dphi;
    switch(data->Nchannel)
    {
        case 1:
            
            //amplitude
            dA   = model->calibration->dampX;
            logP += log(gaussian_pdf(dA,0,CAL_SIGMA_AMP));
            
            //phase
            dphi = model->calibration->dphiX;
            logP += log(gaussian_pdf(dphi,0,CAL_SIGMA_PHASE));
            
            break;
        case 2:
            
            //amplitude
            dA   = model->calibration->dampA;
            logP += log(gaussian_pdf(dA,0,CAL_SIGMA_AMP));
            
            //phase
            dphi = model->calibration->dphiA;
            logP += log(gaussian_pdf(dphi,0,CAL_SIGMA_PHASE));
            
            //amplitude
            dA   = model->calibration->dampE;
            logP += log(gaussian_pdf(dA,0,CAL_SIGMA_AMP));
            
            //phase
            dphi = model->calibration->dphiE;
            logP += log(gaussian_pdf(dphi,0,CAL_SIGMA_PHASE));
            break;
        case 3:
            
            //amplitude
            dA   = model->calibration->dampX;
            logP += log(gaussian_pdf(dA,0,CAL_SIGMA_AMP));
            
            //phase
            dphi = model->calibration->dphiX;
            logP += log(gaussian_pdf(dphi,0,CAL_SIGMA_PHASE));
            
            //amplitude
            dA   = model->calibration->dampY;
            logP += log(gaussian_pdf(dA,0,CAL_SIGMA_AMP));
            
            //phase
            dphi = model->calibration->dphiY;
            logP += log(gaussian_pdf(dphi,0,CAL_SIGMA_PHASE));

            //amplitude
            dA   = model->calibration->dampZ;
            logP += log(gaussian_pdf(dA,0,CAL_SIGMA_AMP));
            
            //phase
            dphi = model->calibration->dphiZ;
            logP += log(gaussian_pdf(dphi,0,CAL_SIGMA_PHASE));
            break;

        default:
            break;
    }//end switch

    return logP;
}



