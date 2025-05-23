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
#include <glass_noise.h>

#include "glass_ucb_model.h"
#include "glass_ucb_io.h"
#include "glass_ucb_catalog.h"
#include "glass_ucb_waveform.h"
#include "glass_ucb_fstatistic.h"
#include "glass_ucb_data.h"


void UCBInjectVerificationSet(struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Source *inj)
{
    //TODO: Combine this function w/ UCBInjectVerificationSource()

    //Book-keeping of injection time-frequency volume
    ucb_alignment(orbit, data, inj);

    ucb_waveform(orbit, data->format, data->T, data->t0, inj->params, UCB_MODEL_NP, inj->tdi->X, inj->tdi->Y, inj->tdi->Z, inj->tdi->A, inj->tdi->E, inj->BW, data->Nchannel);
    
    //Add waveform to data TDI channels
    for(int j=0; j<inj->BW; j++)
    {
        int i = j+inj->imin;
        if(i>0 && i<data->NFFT)
        {
            data->tdi->X[2*i]   += inj->tdi->X[2*j];
            data->tdi->X[2*i+1] += inj->tdi->X[2*j+1];

            data->tdi->Y[2*i]   += inj->tdi->Y[2*j];
            data->tdi->Y[2*i+1] += inj->tdi->Y[2*j+1];

            data->tdi->Z[2*i]   += inj->tdi->Z[2*j];
            data->tdi->Z[2*i+1] += inj->tdi->Z[2*j+1];

            data->tdi->A[2*i]   += inj->tdi->A[2*j];
            data->tdi->A[2*i+1] += inj->tdi->A[2*j+1];
            
            data->tdi->E[2*i]   += inj->tdi->E[2*j];
            data->tdi->E[2*i+1] += inj->tdi->E[2*j+1];
        }
    }

    //Get noise spectrum for data segment
    GetNoiseModel(data,orbit,flags);

    //Add Gaussian noise to injection
    if(flags->simNoise) AddNoise(data,data->tdi);

}

void UCBInjectVerificationSource(struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Source *inj)
{
    //TODO: support Michelson-only injection
    if(!flags->quiet) fprintf(stdout,"\n==== UCBInjectVerificationSource ====\n");
    
    FILE *fptr;
    
    /* structure for holding injection source parameters and waveforms */
    alloc_source(inj, data->N, data->Nchannel);
    
    /* Get injection parameters */
    double f0,dfdt,costheta,phi,m1,m2,D; //read from injection file
    double cosi,phi0,psi;                //drawn from prior
    double Mc,amp;                       //calculated
    
    FILE *injectionFile;
    FILE *paramFile;
    char filename[MAXSTRINGSIZE];
    char header[MAXSTRINGSIZE];
    
    struct TDI *tdi = data->tdi;

    for(int ii = 0; ii<flags->NINJ; ii++)
    {
        
        injectionFile = fopen(flags->injFile[ii],"r");
        if(!injectionFile)
            fprintf(stderr,"Missing injection file %s\n",flags->injFile[ii]);
        else
            fprintf(stdout,"Injecting verification binary %s  (%i/%i)\n",flags->injFile[ii],ii+1, flags->NINJ);
        
        //strip off header
        char *line = fgets(header, MAXSTRINGSIZE, injectionFile);
        if(line==NULL)
        {
            fprintf(stderr,"Error reading %s\n",flags->injFile[ii]);
            exit(1);
        }
        
        //parse injection parameters
        int check = fscanf(injectionFile,"%lg %lg %lg %lg %lg %lg %lg %lg",&f0,&dfdt,&costheta,&phi,&m1,&m2,&cosi,&D);
        if(!check)
        {
            fprintf(stderr,"Error reading %s\n",flags->injFile[ii]);
            exit(1);
        }
        
        
        //incoming distance in kpc, function expects pc
        D *= 1000.0;
        
        
        //draw extrinsic parameters
        
        //set RNG for injection
        unsigned int r = data->iseed;
        
        //TODO: support for verification binary priors
        phi0 = rand_r_U_0_1(&r)*M_PI*2.;
        psi  = rand_r_U_0_1(&r)*M_PI/4.;
        
        
        //compute derived parameters
        Mc  = chirpmass(m1,m2);
        amp = amplitude(Mc, f0, D);
                
        //set bandwidth of data segment centered on injection
        data->fmin = f0 - (data->NFFT/2)/data->T;
        data->fmax = f0 + (data->NFFT/2)/data->T;
        data->qmin = (int)(data->fmin*data->T);
        data->qmax = data->qmin+data->NFFT;
        
        //recompute fmin and fmax so they align with a bin
        data->fmin = data->qmin/data->T;
        data->fmax = data->qmax/data->T;
        
        if(!flags->quiet)fprintf(stdout,"Frequency bins for segment [%i,%i]\n",data->qmin,data->qmax);
        if(!flags->quiet) fprintf(stdout,"   ...start time  %g\n",data->t0);
        
        for(int n=0; n<data->N; n++)
        {
            inj->tdi->A[n] = 0.0;
            inj->tdi->E[n] = 0.0;
            inj->tdi->X[n] = 0.0;
            inj->tdi->Y[n] = 0.0;
            inj->tdi->Z[n] = 0.0;
        }
        
        //map parameters to vector
        inj->f0       = f0;
        inj->dfdt     = dfdt;
        inj->costheta = costheta;
        inj->phi      = phi;
        inj->amp      = amp;
        inj->cosi     = cosi;
        inj->phi0     = phi0;
        inj->psi      = psi;
        map_params_to_array(inj, inj->params, data->T);
        
        //save parameters to file
        sprintf(filename,"%s/injection_parameters_%i.dat",flags->runDir,ii);
        paramFile=fopen(filename,"w");
        fprintf(paramFile,"%lg ",data->t0);
        print_source_params(data, inj, paramFile);
        fprintf(paramFile,"\n");
        fclose(paramFile);
        
        //Book-keeping of injection time-frequency volume
        ucb_alignment(orbit, data, inj);
        
        //Simulate gravitational wave signal
        ucb_waveform(orbit, data->format, data->T, data->t0, inj->params, 8, inj->tdi->X, inj->tdi->Y, inj->tdi->Z, inj->tdi->A, inj->tdi->E, inj->BW, 2);
        
        //Add waveform to data TDI channels
        for(int n=0; n<inj->BW; n++)
        {
            int i = n+inj->imin;
            
            tdi->X[2*i]   = inj->tdi->X[2*n];
            tdi->X[2*i+1] = inj->tdi->X[2*n+1];

            tdi->Y[2*i]   = inj->tdi->Y[2*n];
            tdi->Y[2*i+1] = inj->tdi->Y[2*n+1];

            tdi->Z[2*i]   = inj->tdi->Z[2*n];
            tdi->Z[2*i+1] = inj->tdi->Z[2*n+1];

            tdi->A[2*i]   = inj->tdi->A[2*n];
            tdi->A[2*i+1] = inj->tdi->A[2*n+1];
            
            tdi->E[2*i]   = inj->tdi->E[2*n];
            tdi->E[2*i+1] = inj->tdi->E[2*n+1];
        }
        
        sprintf(filename,"%s/data/waveform_injection_%i.dat",flags->runDir,ii);
        fptr=fopen(filename,"w");
        for(int i=0; i<data->NFFT; i++)
        {
            double f = (double)(i+data->qmin)/data->T;
            switch(data->Nchannel)
            {
                case 1:
                    fprintf(fptr,"%lg %lg %lg\n", f, tdi->X[2*i],tdi->X[2*i+1]);
                    break;
                case 2:
                    fprintf(fptr,"%lg %lg %lg %lg %lg\n", f, tdi->A[2*i],tdi->A[2*i+1], tdi->E[2*i],tdi->E[2*i+1]);
                    break;
                case 3:
                    fprintf(fptr,"%lg %lg %lg %lg %lg %lg %lg\n", f, tdi->X[2*i],tdi->X[2*i+1], tdi->Y[2*i],tdi->Y[2*i+1], tdi->Z[2*i],tdi->Z[2*i+1]);
                    break;
            }
        }
        fclose(fptr);
        
        sprintf(filename,"%s/data/power_injection_%i.dat",flags->runDir,ii);
        fptr=fopen(filename,"w");
        for(int i=0; i<data->NFFT; i++)
        {
            double f = (double)(i+data->qmin)/data->T;
            switch(data->Nchannel)
            {
                case 1:
                    fprintf(fptr,"%.12g %lg\n", f, tdi->X[2*i]*tdi->X[2*i]+tdi->X[2*i+1]*tdi->X[2*i+1]);
                    break;
                case 2:
                    fprintf(fptr,"%.12g %lg %lg\n", f, tdi->A[2*i]*tdi->A[2*i]+tdi->A[2*i+1]*tdi->A[2*i+1], tdi->E[2*i]*tdi->E[2*i]+tdi->E[2*i+1]*tdi->E[2*i+1]);
                    break;
                case 3:
                    fprintf(fptr,"%.12g %lg %lg %lg\n", f, tdi->X[2*i]*tdi->X[2*i]+tdi->X[2*i+1]*tdi->X[2*i+1], tdi->Y[2*i]*tdi->Y[2*i]+tdi->Y[2*i+1]*tdi->Y[2*i+1],tdi->Z[2*i]*tdi->Z[2*i]+tdi->Z[2*i+1]*tdi->Z[2*i+1]);
                    break;
            }
        }
        fclose(fptr);
        
        //Get noise spectrum for data segment
        GetNoiseModel(data,orbit,flags);
        
        
        //Get injected SNR
        if(!flags->quiet) fprintf(stdout,"   ...injected SNR=%g\n",snr(inj, data->noise));
        
        //Add Gaussian noise to injection
        if(flags->simNoise)
            AddNoise(data,tdi);
        
        //Compute fisher information matrix of injection
        if(!flags->quiet) fprintf(stdout,"   ...computing Fisher Information Matrix of injection\n");
        
        ucb_fisher(orbit, data, inj, data->noise);
        
        fclose(injectionFile);
    }
    
    print_data(data,flags);
    
    if(!flags->quiet)fprintf(stdout,"================================================\n\n");
}
void UCBInjectSimulatedSource(struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Source **inj_vec)
{
    int k;
    FILE *fptr;
    
    /* Get injection parameters */
    double f0,dfdt,theta,phi,amp,iota,phi0,psi;//read from injection file
    
    FILE *injectionFile;
    FILE *paramFile;
    char filename[1024];
        
    struct TDI *tdi = data->tdi;
    struct Source *inj = NULL;
    
    int n_inj = 0; // keep track of number of injections
    for(int ii = 0; ii<flags->NINJ; ii++)
    {
        injectionFile = fopen(flags->injFile[ii],"r");
        if(!injectionFile)
            fprintf(stderr,"Missing injection file %s\n",flags->injFile[ii]);
        else
            fprintf(stdout,"Injecting simulated source %s  (%i/%i)\n",flags->injFile[ii],ii+1, flags->NINJ);
        
        
        //count sources in file
        int N=0;
        while(!feof(injectionFile))
        {
            int check = fscanf(injectionFile,"%lg %lg %lg %lg %lg %lg %lg %lg",&f0,&dfdt,&theta,&phi,&amp,&iota,&psi,&phi0);
            if(!check)
            {
                fprintf(stderr,"Error reading %s\n",flags->injFile[ii]);
                exit(1);
            }
            N++;
        }
        rewind(injectionFile);
        N--;
        
        for(int nn=0; nn<N; nn++)
        {
            if(n_inj<flags->DMAX)
            {
                alloc_source(inj_vec[n_inj], data->N, data->Nchannel);
                inj = inj_vec[n_inj];
            }
            else
            {
                printf("WARNING: number of injections (%i) exceeds size of model (%i)\n",n_inj,flags->DMAX);
                inj = inj_vec[flags->DMAX-1];
            }
            

            int check = fscanf(injectionFile,"%lg %lg %lg %lg %lg %lg %lg %lg",&f0,&dfdt,&theta,&phi,&amp,&iota,&psi,&phi0);
            if(!check)
            {
                fprintf(stderr,"Error reading %s\n",flags->injFile[ii]);
                exit(1);
            }
            
            //set bandwidth of data segment centered on injection
            if(nn==0 && !flags->strainData)
            {
                if(!strcmp("fourier",data->basis))
                {
                    data->fmin = f0 - (data->NFFT/2)/data->T;
                    data->fmax = f0 + (data->NFFT/2)/data->T;
                    data->qmin = (int)(data->fmin*data->T);
                    data->qmax = data->qmin+data->NFFT;

                    //recompute fmin and fmax so they align with a bin
                    data->fmin = data->qmin/data->T;
                    data->fmax = data->qmax/data->T;

                    if(!flags->quiet)fprintf(stdout,"Frequency bins for segment [%i,%i]\n",data->qmin,data->qmax);
                }
                if(!strcmp("wavelet",data->basis))
                {
                    
                    int q0 = (int)floor(f0/WAVELET_BANDWIDTH);

                    //pad by one layer above and below
                    data->lmin = q0-((data->Nlayer-1)/2);   //mimimum frequency layer
                    data->lmax = data->lmin + data->Nlayer; //maximum frequency layer
                    
                    //reset wavelet basis max and min ranges
                    wavelet_pixel_to_index(data->wdm,0,data->lmin,&data->wdm->kmin);
                    wavelet_pixel_to_index(data->wdm,0,data->lmax,&data->wdm->kmax);
                    
                    //recompute fmin and fmax so they align with a bin
                    data->fmin = data->lmin*WAVELET_BANDWIDTH;
                    data->fmax = data->lmax*WAVELET_BANDWIDTH;
                    data->qmin = (int)(data->fmin*data->T);
                    data->qmax = data->qmin+data->NFFT;

                    if(!flags->quiet)fprintf(stdout,"  Frequency layers [%i,%i]\n",data->lmin,data->lmax);
                    
                    printf("  fmin=%lg, fmax=%lg\n",data->fmin,data->fmax);
                }
                
                
                if(!flags->quiet)fprintf(stdout,"   ...start time: %g\n",data->t0);
            }
            
            for(int n=0; n<data->N; n++)
            {
                inj->tdi->A[n] = 0.0;
                inj->tdi->E[n] = 0.0;
                inj->tdi->X[n] = 0.0;
                inj->tdi->Y[n] = 0.0;
                inj->tdi->Z[n] = 0.0;
            }
            
            //map polarization angle into [0:pi], preserving relation to phi0
            if(psi>M_PI) psi  -= M_PI;
            if(phi0>PI2) phi0 -= PI2;
            
            //map parameters to vector
            inj->f0       = f0;
            inj->dfdt     = dfdt;
            inj->costheta = cos(M_PI/2. - theta);
            inj->phi      = phi;
            inj->amp      = amp;
            inj->cosi     = cos(iota);
            inj->phi0     = phi0;
            inj->psi      = psi;
            if(UCB_MODEL_NP>8)
                inj->d2fdt2 = 11.0/3.0*dfdt*dfdt/f0;
            //inj->d2fdt2 = fddot;
            
            map_params_to_array(inj, inj->params, data->T);
            
            //save parameters to file
            sprintf(filename,"%s/injection_parameters_%i.dat",flags->runDir,ii);
            if(nn==0)paramFile=fopen(filename,"w");
            else     paramFile=fopen(filename,"a");
            fprintf(paramFile,"%lg ",data->t0);
            print_source_params(data, inj, paramFile);
            fprintf(paramFile,"\n");
            fclose(paramFile);
            
            //Book-keeping of injection time-frequency volume
            if(!strcmp("fourier",data->basis))
            {
                ucb_alignment(orbit, data, inj);
                if(inj->qmax < data->qmin || inj->qmin > data->qmax)
                {
                    fprintf(stdout,"Injection %i is outside of the requested frequency segment\n",nn);
                    continue;
                }
                printf("   ...bandwidth : %i\n",inj->BW);
                printf("   ...fdot      : %g\n",inj->dfdt*data->T*data->T);
                printf("   ...fddot     : %g\n",inj->d2fdt2*data->T*data->T*data->T);
                printf("   ...t0        : %g\n",data->t0);
                
            }

            //Simulate gravitational wave signal
            if(!strcmp("fourier",data->basis))     
                ucb_waveform(orbit, data->format, data->T, data->t0, inj->params, UCB_MODEL_NP, inj->tdi->X, inj->tdi->Y, inj->tdi->Z, inj->tdi->A, inj->tdi->E, inj->BW, data->Nchannel);
            if(!strcmp("wavelet",data->basis)) 
                ucb_waveform_wavelet(orbit, data->wdm, data->T, data->t0, inj->params, inj->list, &inj->Nlist, inj->tdi->X, inj->tdi->Y, inj->tdi->Z);
            
            
            //Add waveform to data TDI channels
            if(!strcmp("fourier",data->basis))
            {
                printf("add fourier waveform to data\n");
                for(int n=0; n<inj->BW; n++)
                {
                    int i = n+inj->imin;
                    if(i>0 && i<data->NFFT)
                    {
                        tdi->X[2*i]   += inj->tdi->X[2*n];
                        tdi->X[2*i+1] += inj->tdi->X[2*n+1];

                        tdi->Y[2*i]   += inj->tdi->Y[2*n];
                        tdi->Y[2*i+1] += inj->tdi->Y[2*n+1];

                        tdi->Z[2*i]   += inj->tdi->Z[2*n];
                        tdi->Z[2*i+1] += inj->tdi->Z[2*n+1];

                        tdi->A[2*i]   += inj->tdi->A[2*n];
                        tdi->A[2*i+1] += inj->tdi->A[2*n+1];
                        
                        tdi->E[2*i]   += inj->tdi->E[2*n];
                        tdi->E[2*i+1] += inj->tdi->E[2*n+1];
                    }
                }
            }
            if(!strcmp("wavelet",data->basis))
            {
                for(int n=0; n<inj->Nlist; n++)
                {
                    k = inj->list[n];
                    tdi->X[k] += inj->tdi->X[k];
                    tdi->Y[k] += inj->tdi->Y[k];
                    tdi->Z[k] += inj->tdi->Z[k];
                }
            }
            
            sprintf(filename,"%s/data/waveform_injection_%i.dat",flags->runDir,n_inj);
            fptr=fopen(filename,"w");
            if(!strcmp("fourier",data->basis))
            {
                for(int i=0; i<data->NFFT; i++)
                {
                    double f = (double)(i+data->qmin)/data->T;
                    switch(data->Nchannel)
                    {
                        case 1:
                            fprintf(fptr,"%.14e %.14e %.14e\n", f, inj->tdi->X[2*i],tdi->X[2*i+1]);
                            break;
                        case 2:
                            fprintf(fptr,"%.14e %.14e %.14e %.14e %.14e\n", f, inj->tdi->A[2*i],tdi->A[2*i+1], tdi->E[2*i],tdi->E[2*i+1]);
                            break;
                        case 3:
                            fprintf(fptr,"%.14e %.14e %.14e %.14e %.14e %.14e %.14e\n", f, inj->tdi->X[2*i],inj->tdi->X[2*i+1], inj->tdi->Y[2*i],inj->tdi->Y[2*i+1],inj->tdi->Z[2*i],inj->tdi->Z[2*i+1]);
                            break;
                    }
                }
            }
            if(!strcmp("wavelet",data->basis))
            {
                for(int j=data->lmin; j<data->lmax; j++)
                {
                    double f = j*data->wdm->df;
                    for(int i=0; i<data->wdm->NT; i++)
                    {
                        double t = i*data->wdm->dt;
                        wavelet_pixel_to_index(data->wdm,i,j,&k);
                        k-=data->wdm->kmin;
                        fprintf(fptr,"%.14e %.14e %.14e %.14e %.14e\n", t, f, inj->tdi->X[k], inj->tdi->Y[k], inj->tdi->Z[k]);
                    }
                    fprintf(fptr,"\n");
                }
            }
            fclose(fptr);

            if(!strcmp("wavelet",data->basis))
            {
                sprintf(filename,"%s/waveform_injection_fourier_%i.dat",data->dataDir,n_inj);
                print_wavelet_fourier_spectra(data, inj->tdi, filename);
            }

            
            sprintf(filename,"%s/data/power_injection_%i.dat",flags->runDir,n_inj);
            fptr=fopen(filename,"w");
            if(!strcmp("fourier",data->basis))
            {
                for(int i=0; i<inj->BW; i++)
                {
                    double f = (double)(i+inj->qmin)/data->T;
                    switch(data->Nchannel)
                    {
                        case 1:
                            fprintf(fptr,"%.14e %.14e\n", f, inj->tdi->X[2*i]*inj->tdi->X[2*i]+inj->tdi->X[2*i+1]*inj->tdi->X[2*i+1]);
                            break;
                        case 2:
                            fprintf(fptr,"%.14e %.14e %.14e\n", f, inj->tdi->A[2*i]*inj->tdi->A[2*i]+inj->tdi->A[2*i+1]*inj->tdi->A[2*i+1], inj->tdi->E[2*i]*inj->tdi->E[2*i]+inj->tdi->E[2*i+1]*inj->tdi->E[2*i+1]);
                            break;
                        case 3:
                            fprintf(fptr,"%.14e %.14e %.14e %.14e\n", f, 
                                inj->tdi->X[2*i]*inj->tdi->X[2*i]+inj->tdi->X[2*i+1]*inj->tdi->X[2*i+1], 
                                inj->tdi->Y[2*i]*inj->tdi->Y[2*i]+inj->tdi->Y[2*i+1]*inj->tdi->Y[2*i+1],
                                inj->tdi->Z[2*i]*inj->tdi->Z[2*i]+inj->tdi->Z[2*i+1]*inj->tdi->Z[2*i+1]);
                            break;
                    }
                }
            }
            if(!strcmp("wavelet",data->basis))
            {
                for(int j=data->lmin; j<data->lmax; j++)
                {
                    double f = j*data->wdm->df;
                    for(int i=0; i<data->wdm->NT; i++)
                    {
                        double t = i*data->wdm->dt;
                        wavelet_pixel_to_index(data->wdm,i,j,&k);
                        k-=data->wdm->kmin;
                        fprintf(fptr,"%.14e %.14e %.14e %.14e %.14e\n", t, f, inj->tdi->X[k]*inj->tdi->X[k], inj->tdi->Y[k]*inj->tdi->Y[k], inj->tdi->Z[k]*inj->tdi->Z[k]);
                    }
                    fprintf(fptr,"\n");
                }
            }
            fclose(fptr);
            
            //Get noise spectrum for data segment
            if(!strcmp("fourier",data->basis)) GetNoiseModel(data,orbit,flags);
            if(!strcmp("wavelet",data->basis)) GetDynamicNoiseModel(data, orbit, flags);

            //Get injected SNR
            if(!flags->quiet)
            {
                if(!strcmp("fourier",data->basis))fprintf(stdout,"   ...injected SNR=%g\n",snr(inj, data->noise));
                if(!strcmp("wavelet",data->basis))fprintf(stdout,"   ...injected SNR=%g\n",snr_wavelet(inj,data->noise));
            }
            
            /*Add Gaussian noise to injection
            if(flags->simNoise && nn==0)
            {
                if(!strcmp("fourier",data->basis)) AddNoise(data,tdi);
                if(!strcmp("wavelet",data->basis)) AddNoiseWavelet(data,tdi);
            }*/

            //Compute fisher information matrix of injection
            if(!flags->quiet)fprintf(stdout,"   ...computing Fisher Information Matrix of injection\n");
            
            if(!strcmp("fourier",data->basis)) ucb_fisher(orbit, data, inj, data->noise);
            if(!strcmp("wavelet",data->basis)) ucb_fisher_wavelet(orbit, data, inj, data->noise);
            
            if(!flags->quiet)
            {
                fprintf(stdout,"\n Fisher Matrix:\n");
                for(int i=0; i<UCB_MODEL_NP; i++)
                {
                    fprintf(stdout," ");
                    for(int j=0; j<UCB_MODEL_NP; j++)
                    {
                        fprintf(stdout,"%+.2e ", inj->fisher_matrix[i][j]);
                    }
                    fprintf(stdout,"\n");
                }
                
                
                printf("\n Fisher std. errors:\n");
                for(int j=0; j<UCB_MODEL_NP; j++)  fprintf(stdout," %.2e\n", sqrt(inj->fisher_matrix[j][j]));
            }
            
            n_inj++;
            
        }//end nn loop over sources in file
        
        fclose(injectionFile);        
    }//end ii loop over injection files
    
    //print_data(data,flags);

    /* compute overlaps between injections */
    if(!strcmp(data->basis,"wavelet"))
    {
        printf("\n Match Matrix:\n");
        for(int i=0; i<n_inj; i++)
        {
            for(int j=0; j<n_inj; j++)
            {
                if(i==j) printf(" %+.2e",snr_wavelet(inj_vec[i],data->noise));
                else printf(" %+.2e",waveform_match_wavelet(inj_vec[i], inj_vec[j], data->noise));
            }
            printf("\n");
        }
        
        sprintf(filename,"%s/data/waveform_injection_fourier.dat",flags->runDir);
        print_wavelet_fourier_spectra(data, inj_vec[0]->tdi, filename);
    }

    if(!flags->quiet)fprintf(stdout,"================================================\n\n");
}

void GetVerificationBinary(struct Data *data, struct Flags *flags, struct Source *inj, FILE *vbFile)
{
    /* Get injection parameters */
    double f0,dfdt,costheta,phi,m1,m2,D; //read from injection file
    double cosi,phi0,psi;                //drawn from prior
    double Mc,amp;                       //calculated
    
    int check = fscanf(vbFile,"%lg %lg %lg %lg %lg %lg %lg %lg",&f0,&dfdt,&costheta,&phi,&m1,&m2,&cosi,&D);
    if(!check)
    {
        fprintf(stderr,"Error reading %s\n",flags->vbFile);
        exit(1);
    }
    
    //incoming distance in kpc, function expects pc
    D *= 1000.0;
    
    //compute derived parameters
    Mc  = chirpmass(m1,m2);
    amp = amplitude(Mc, f0, D);
    
    //initialize extrinsic parameters
    phi0 = 0.0;
    psi  = 0.0;
    
    //set bandwidth of data segment centered on injection
    data->fmin = f0 - (data->NFFT)/data->T;
    data->fmax = f0 + (data->NFFT)/data->T;
    data->qmin = (int)(data->fmin*data->T);
    data->qmax = data->qmin+data->NFFT;
    
    //recompute fmin and fmax so they align with a bin
    data->fmin = data->qmin/data->T;
    data->fmax = data->qmax/data->T;
    
    //map parameters to vector
    inj->f0       = f0;
    inj->dfdt     = dfdt;
    inj->costheta = costheta;
    inj->phi      = phi;
    inj->amp      = amp;
    inj->cosi     = cosi;
    inj->phi0     = phi0;
    inj->psi      = psi;
    map_params_to_array(inj, inj->params, data->T);
}



void UCBLoadCatalogCache(struct Data *data, struct Flags *flags, struct Catalog *catalog)
{
    /* check that file exists */
    FILE *catalog_file = NULL;
    if((catalog_file = fopen(flags->catalogFile, "r")) == NULL)
    {
        fprintf(stderr,"Error opening %s\n", flags->catalogFile);
        exit(1);
    }
    
    /* count sources in catalog cache file */
    char* line;
    char lineBuffer[MAXSTRINGSIZE];
    int Ncache=0;
    while((line = fgets(lineBuffer, MAXSTRINGSIZE, catalog_file)) != NULL) Ncache++;
    rewind(catalog_file);
    
    
    /* store entire cache file */
    char **cache; //!<contents of cache file
    cache = malloc(Ncache*sizeof(char *));
    for(int n=0; n<Ncache; n++)
    {
        cache[n] = (char *)malloc(MAXSTRINGSIZE);
        
        line = fgets(lineBuffer, MAXSTRINGSIZE, catalog_file);
        strcpy(cache[n],line);
        strtok(cache[n],"\n");
    }
    
    fclose(catalog_file);
    
    /* allocate enough space for the entire catalog */
    catalog->entry = malloc(Ncache*sizeof(struct Entry *));

    catalog->N = 0;
    
    /* only store catalog sources in current segment */
    double fmin = data->fmin + data->qpad/data->T;
    double fmax = data->fmax - data->qpad/data->T;

    /* parse catalog cache file */
    char *token[6];
    double f0;
    
    for(int n=0; n<Ncache; n++)
    {
        //split long string in cache file into tokens
        int i=0;
        token[i] = strtok(cache[n]," ");
        while(token[i]!=NULL)
        {
            i++;
            token[i] = strtok(NULL," \n");
        }

        f0 = atof(token[1]);
        
        if(f0>fmin && f0<fmax)
        {
            /* allocate  memory for catalog entry */
            create_empty_source(catalog, data->N, data->Nchannel);
            
            /* assign contents of cache file to entry */
            struct Entry *entry = catalog->entry[catalog->N-1];
            strcpy(entry->name,token[0]);
            entry->SNR = atof(token[2]);
            entry->evidence = atof(token[3]);
            strcpy(entry->path,token[4]);
        }
    }
    
    /* load catalog from cache file */
    for(int n=0; n<catalog->N; n++)
    {
        char filename[MAXSTRINGSIZE];
        
        struct Entry *entry = catalog->entry[n];
        struct Source *source = catalog->entry[n]->source[0];
        struct GMM *gmm = catalog->entry[n]->gmm;
        
        /* gaussian mixture model */
        gmm->NParams = (size_t)UCB_MODEL_NP;
        sprintf(filename,"%s%s_gmm.bin",entry->path,entry->name);
        read_gmm_binary(gmm, filename);
        
        /* source parameters */
        sprintf(filename,"%s%s_params.dat",entry->path,entry->name);
        FILE *fptr = NULL;
        if((fptr = fopen(filename,"r")) == NULL)
        {
            fprintf(stderr,"Error opening %s\n", filename);
            exit(1);
        }
        scan_source_params(data, source, fptr);
        fclose(fptr);

    }
}

