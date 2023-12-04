/*
 *  Copyright (C) 2019 Tyson B. Littenberg (MSFC-ST12), Neil J. Cornish
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

#include "glass_ucb_model.h"
#include "glass_ucb_io.h"
#include "glass_ucb_catalog.h"
#include "glass_ucb_waveform.h"
#include "glass_ucb_data.h"


void UCBInjectVerificationSet(struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Source *inj)
{
    //TODO: Combine this function w/ UCBInjectVerificationSource()

    //Book-keeping of injection time-frequency volume
    galactic_binary_alignment(orbit, data, inj);

    galactic_binary(orbit, data->format, data->T, data->t0, inj->params, UCB_MODEL_NP, inj->tdi->X, inj->tdi->Y, inj->tdi->Z, inj->tdi->A, inj->tdi->E, inj->BW, data->Nchannel);
    
    //Add waveform to data TDI channels
    for(int j=0; j<inj->BW; j++)
    {
        int i = j+inj->imin;
        if(i>0 && i<data->N)
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
        const gsl_rng_type *T = gsl_rng_default;
        gsl_rng *r = gsl_rng_alloc(T);
        gsl_rng_env_setup();
        gsl_rng_set (r, data->iseed);
        
        //TODO: support for verification binary priors
        phi0 = gsl_rng_uniform(r)*M_PI*2.;
        psi  = gsl_rng_uniform(r)*M_PI/4.;
        
        
        //compute derived parameters
        Mc  = chirpmass(m1,m2);
        amp = galactic_binary_Amp(Mc, f0, D);
                
        //set bandwidth of data segment centered on injection
        data->fmin = f0 - (data->N/2)/data->T;
        data->fmax = f0 + (data->N/2)/data->T;
        data->qmin = (int)(data->fmin*data->T);
        data->qmax = data->qmin+data->N;
        
        //recompute fmin and fmax so they align with a bin
        data->fmin = data->qmin/data->T;
        data->fmax = data->qmax/data->T;
        
        if(!flags->quiet)fprintf(stdout,"Frequency bins for segment [%i,%i]\n",data->qmin,data->qmax);
        if(!flags->quiet) fprintf(stdout,"   ...start time  %g\n",data->t0);
        
        for(int n=0; n<2*data->N; n++)
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
        galactic_binary_alignment(orbit, data, inj);
        
        //Simulate gravitational wave signal
        galactic_binary(orbit, data->format, data->T, data->t0, inj->params, 8, inj->tdi->X, inj->tdi->Y, inj->tdi->Z, inj->tdi->A, inj->tdi->E, inj->BW, 2);
        
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
        for(int i=0; i<data->N; i++)
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
        for(int i=0; i<data->N; i++)
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
        
        galactic_binary_fisher(orbit, data, inj, data->noise);
        
        fclose(injectionFile);
                
        gsl_rng_free(r);
    }
    
    print_data(data,tdi,flags);
    
    if(!flags->quiet)fprintf(stdout,"================================================\n\n");
}
void UCBInjectSimulatedSource(struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Source *inj)
{
    FILE *fptr;
    alloc_source(inj, data->N, data->Nchannel);
    
    /* Get injection parameters */
    double f0,dfdt,theta,phi,amp,iota,phi0,psi;//read from injection file
    
    FILE *injectionFile;
    FILE *paramFile;
    char filename[1024];
        
    struct TDI *tdi = data->tdi;

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
        
        //set RNG for injection
        const gsl_rng_type *T = gsl_rng_default;
        gsl_rng *r = gsl_rng_alloc(T);
        gsl_rng_env_setup();
        gsl_rng_set (r, data->iseed);
        
        for(int nn=0; nn<N; nn++)
        {
            int check = fscanf(injectionFile,"%lg %lg %lg %lg %lg %lg %lg %lg",&f0,&dfdt,&theta,&phi,&amp,&iota,&psi,&phi0);
            if(!check)
            {
                fprintf(stderr,"Error reading %s\n",flags->injFile[ii]);
                exit(1);
            }
            
            //set bandwidth of data segment centered on injection
            if(nn==0 && !flags->strainData)
            {
                data->fmin = f0 - (data->N/2)/data->T;
                data->fmax = f0 + (data->N/2)/data->T;
                data->qmin = (int)(data->fmin*data->T);
                data->qmax = data->qmin+data->N;
                
                //recompute fmin and fmax so they align with a bin
                data->fmin = data->qmin/data->T;
                data->fmax = data->qmax/data->T;
                
                if(!flags->quiet)fprintf(stdout,"Frequency bins for segment [%i,%i]\n",data->qmin,data->qmax);
                if(!flags->quiet)fprintf(stdout,"   ...start time: %g\n",data->t0);
            }
            
            for(int n=0; n<2*data->N; n++)
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
            galactic_binary_alignment(orbit, data, inj);
            
            if(inj->qmax < data->qmin || inj->qmin > data->qmax)
            {
                fprintf(stdout,"Injection %i is outside of the requested frequency segment\n",nn);
                continue;
            }
            
            printf("   ...bandwidth : %i\n",inj->BW);
            printf("   ...fdot      : %g\n",inj->dfdt*data->T*data->T);
            printf("   ...fddot     : %g\n",inj->d2fdt2*data->T*data->T*data->T);
            
            //Simulate gravitational wave signal
            printf("   ...t0        : %g\n",data->t0);
            galactic_binary(orbit, data->format, data->T, data->t0, inj->params, UCB_MODEL_NP, inj->tdi->X, inj->tdi->Y, inj->tdi->Z, inj->tdi->A, inj->tdi->E, inj->BW, 2);
            
            //Add waveform to data TDI channels
            for(int n=0; n<inj->BW; n++)
            {
                int i = n+inj->imin;
                if(i>0 && i<data->N)
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
            
            sprintf(filename,"%s/data/waveform_injection_%i.dat",flags->runDir,ii);
            fptr=fopen(filename,"w");
            for(int i=0; i<data->N; i++)
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
            for(int i=0; i<data->N; i++)
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
            if(!flags->quiet)fprintf(stdout,"   ...injected SNR=%g\n",snr(inj, data->noise));
            
            //Add Gaussian noise to injection
            if(flags->simNoise && nn==0)
                AddNoise(data,tdi);
            
            //Compute fisher information matrix of injection
            if(!flags->quiet)fprintf(stdout,"   ...computing Fisher Information Matrix of injection\n");
            
            galactic_binary_fisher(orbit, data, inj, data->noise);
            
            
            if(!flags->quiet)
            {
                fprintf(stdout,"\n Fisher Matrix:\n");
                for(int i=0; i<UCB_MODEL_NP; i++)
                {
                    fprintf(stdout," ");
                    for(int j=0; j<UCB_MODEL_NP; j++)
                    {
                        if(inj->fisher_matrix[i][j]<0)fprintf(stdout,"%.2e ", inj->fisher_matrix[i][j]);
                        else                          fprintf(stdout,"+%.2e ",inj->fisher_matrix[i][j]);
                    }
                    fprintf(stdout,"\n");
                }
                
                
                printf("\n Fisher std. errors:\n");
                for(int j=0; j<UCB_MODEL_NP; j++)  fprintf(stdout," %.4e\n", sqrt(inj->fisher_matrix[j][j]));
            }
            
        }//end nn loop over sources in file
        fclose(injectionFile);
        gsl_rng_free(r);
    }
    
    print_data(data,tdi,flags);

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
    amp = galactic_binary_Amp(Mc, f0, D);
    
    //initialize extrinsic parameters
    phi0 = 0.0;
    psi  = 0.0;
    
    //set bandwidth of data segment centered on injection
    data->fmin = f0 - (data->N/2)/data->T;
    data->fmax = f0 + (data->N/2)/data->T;
    data->qmin = (int)(data->fmin*data->T);
    data->qmax = data->qmin+data->N;
    
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

