//
//  GalacticBinaryData.c
//
//
//  Created by Littenberg, Tyson B. (MSFC-ZP12) on 2/3/17.
//
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "LISA.h"
#include "GalacticBinary.h"
#include "GalacticBinaryIO.h"
#include "GalacticBinaryData.h"
#include "GalacticBinaryMath.h"
#include "GalacticBinaryModel.h"
#include "GalacticBinaryWaveform.h"

void GalacticBinaryReadData(struct Data *data)
{
  
}
void GalacticBinarySimulateData(struct Data *data)
{
  
}

void GalacticBinaryInjectVerificationSource(struct Data **data_vec, struct Orbit *orbit, struct Flags *flags)
{
  //TODO: support Michelson-only injection
  fprintf(stdout,"\n==== GalacticBinaryInjectVerificationSource ====\n");
  
  FILE *fptr;
  
  /* Get injection parameters */
  double f0,dfdt,costheta,phi,m1,m2,D; //read from injection file
  double cosi,phi0,psi;                //drawn from prior
  double Mc,amp;                       //calculated
  
  FILE *injectionFile;
  FILE *paramFile;
  char filename[1024];
  
  for(int ii = 0; ii<flags->NF; ii++)
  {
    
    injectionFile = fopen(flags->injFile[ii],"r");
    if(!injectionFile)
      fprintf(stderr,"Missing injection file %s\n",flags->injFile[ii]);
    else
      fprintf(stdout,"Injecting verification binary %s  (%i/%i)\n",flags->injFile[ii],ii+1, flags->NF);
    
    fscanf(injectionFile,"%lg %lg %lg %lg %lg %lg %lg",&f0,&dfdt,&costheta,&phi,&m1,&m2,&D);
    
    //draw extrinsic parameters
    
    //set RNG for injection
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc(T);
    gsl_rng_env_setup();
    gsl_rng_set (r, data_vec[ii]->iseed);
    
    //TODO: support for verification binary priors
    cosi = -1.0 + gsl_rng_uniform(r)*2.0;
    phi0 = gsl_rng_uniform(r)*M_PI*2.;
    psi  = gsl_rng_uniform(r)*M_PI/4.;
    
    //compute derived parameters
    Mc  = chirpmass(m1,m2);
    amp = galactic_binary_Amp(Mc, f0, D, data_vec[0]->T);
    
    for(int jj=0; jj<flags->NT; jj++)
    {
      
      struct Data *data  = data_vec[ii];
      struct TDI *tdi = data->tdi[jj];
      
      //set bandwidth of data segment centered on injection
      data->fmin = f0 - (data->N/2)/data->T;
      data->fmax = f0 + (data->N/2)/data->T;
      data->qmin = (int)(data->fmin*data->T);
      data->qmax = data->qmin+data->N;
      
      //recompute fmin and fmax so they align with a bin
      data->fmin = data->qmin/data->T;
      data->fmax = data->qmax/data->T;
      
      if(jj==0)fprintf(stdout,"Frequency bins for segment [%i,%i]\n",data->qmin,data->qmax);
      fprintf(stdout,"   ...start time  %g\n",data->t0[jj]);
      
      
      struct Source *inj = data->inj;
      
      for(int n=0; n<2*data->N; n++)
      {
        inj->tdi->A[n] = 0.0;
        inj->tdi->E[n] = 0.0;
        inj->tdi->X[n] = 0.0;
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
      sprintf(filename,"injection_parameters_%i_%i.dat",ii,jj);
      paramFile=fopen(filename,"w");
      fprintf(paramFile,"%lg ",data->t0[jj]);
      print_source_params(data, inj, paramFile);
      fprintf(paramFile,"\n");
      fclose(paramFile);
      
      //Book-keeping of injection time-frequency volume
      galactic_binary_alignment(orbit, data, inj);
      
      //Simulate gravitational wave signal
      //double t0 = data->t0 + jj*(data->T + data->tgap);
      galactic_binary(orbit, data->T, data->t0[jj], inj->params, 8, inj->tdi->X, inj->tdi->A, inj->tdi->E, inj->BW, 2);
      
      //Add waveform to data TDI channels
      for(int n=0; n<inj->BW; n++)
      {
        int i = n+inj->imin;
        
        tdi->X[2*i]   = inj->tdi->X[2*n];
        tdi->X[2*i+1] = inj->tdi->X[2*n+1];
        
        tdi->A[2*i]   = inj->tdi->A[2*n];
        tdi->A[2*i+1] = inj->tdi->A[2*n+1];
        
        tdi->E[2*i]   = inj->tdi->E[2*n];
        tdi->E[2*i+1] = inj->tdi->E[2*n+1];
      }
      
      sprintf(filename,"waveform_injection_%i_%i.dat",ii,jj);
      fptr=fopen(filename,"w");
      for(int i=0; i<data->N; i++)
      {
        double f = (double)(i+data->qmin)/data->T;
        fprintf(fptr,"%lg %lg %lg %lg %lg",
                f,
                tdi->A[2*i],tdi->A[2*i+1],
                tdi->E[2*i],tdi->E[2*i+1]);
        fprintf(fptr,"\n");
      }
      fclose(fptr);
      
      sprintf(filename,"power_injection_%i_%i.dat",ii,jj);
      fptr=fopen(filename,"w");
      for(int i=0; i<data->N; i++)
      {
        double f = (double)(i+data->qmin)/data->T;
        fprintf(fptr,"%.12g %lg %lg ",
                f,
                tdi->A[2*i]*tdi->A[2*i]+tdi->A[2*i+1]*tdi->A[2*i+1],
                tdi->E[2*i]*tdi->E[2*i]+tdi->E[2*i+1]*tdi->E[2*i+1]);
        fprintf(fptr,"\n");
      }
      fclose(fptr);
      
      //Get noise spectrum for data segment
      for(int n=0; n<data->N; n++)
      {
        double f = data->fmin + (double)(n)/data->T;
        data->noise[jj]->SnA[n] = AEnoise(orbit->L, orbit->fstar, f);
        data->noise[jj]->SnE[n] = AEnoise(orbit->L, orbit->fstar, f);
      }
      
      //Get injected SNR
      fprintf(stdout,"   ...injected SNR=%g\n",snr(inj, data->noise[jj]));
      
      //Add Gaussian noise to injection
      gsl_rng_set (r, data->nseed+jj);
      
      if(!flags->zeroNoise)
      {
        printf("   ...adding Gaussian noise realization\n");
        
        for(int n=0; n<data->N; n++)
        {
          tdi->A[2*n]   += gsl_ran_gaussian (r, 1)*sqrt(data->noise[jj]->SnA[n])/2.;
          tdi->A[2*n+1] += gsl_ran_gaussian (r, 1)*sqrt(data->noise[jj]->SnA[n])/2.;
          
          tdi->E[2*n]   += gsl_ran_gaussian (r, sqrt(data->noise[jj]->SnE[n])/2.);
          tdi->E[2*n+1] += gsl_ran_gaussian (r, sqrt(data->noise[jj]->SnE[n])/2.);
        }
      }
      
      //Compute fisher information matrix of injection
      printf("   ...computing Fisher Information Matrix of injection\n");
      
      galactic_binary_fisher(orbit, data, inj, data->noise[jj]);
      
      /*
       printf("\n Fisher Matrix:\n");
       for(int i=0; i<8; i++)
       {
       fprintf(stdout," ");
       for(int j=0; j<8; j++)
       {
       if(inj->fisher_matrix[i][j]<0)fprintf(stdout,"%.2e ", inj->fisher_matrix[i][j]);
       else                          fprintf(stdout,"+%.2e ",inj->fisher_matrix[i][j]);
       }
       fprintf(stdout,"\n");
       }
       
       printf("\n Fisher std. errors:\n");
       for(int j=0; j<8; j++)  fprintf(stdout," %.4e\n", 1./sqrt(inj->fisher_evalue[j]));
       */
      
      
      sprintf(filename,"power_data_%i_%i.dat",ii,jj);
      fptr=fopen(filename,"w");
      
      for(int i=0; i<data->N; i++)
      {
        double f = (double)(i+data->qmin)/data->T;
        fprintf(fptr,"%.12g %lg %lg ",
                f,
                tdi->A[2*i]*tdi->A[2*i]+tdi->A[2*i+1]*tdi->A[2*i+1],
                tdi->E[2*i]*tdi->E[2*i]+tdi->E[2*i+1]*tdi->E[2*i+1]);
        fprintf(fptr,"\n");
      }
      fclose(fptr);
      fclose(injectionFile);
    }//end jj loop over time segments
    gsl_rng_free(r);
  }
  
  fprintf(stdout,"================================================\n\n");
}
void GalacticBinaryInjectSimulatedSource(struct Data **data_vec, struct Orbit *orbit, struct Flags *flags)
{
  //TODO: support Michelson-only injection
  fprintf(stdout,"\n==== GalacticBinaryInjectSimulatedSource ====\n");
  
  FILE *fptr;
  
  /* Get injection parameters */
  double f0,dfdt,theta,phi,amp,iota,phi0,psi;//read from injection file
  
  FILE *injectionFile;
  FILE *paramFile;
  char filename[1024];
  
  for(int ii = 0; ii<flags->NF; ii++)
  {
    
    injectionFile = fopen(flags->injFile[ii],"r");
    if(!injectionFile)
      fprintf(stderr,"Missing injection file %s\n",flags->injFile[ii]);
    else
      fprintf(stdout,"Injecting simulated source %s  (%i/%i)\n",flags->injFile[ii],ii+1, flags->NF);
    
    
    //count sources in file
    int N=0;
    while(!feof(injectionFile))
    {
      fscanf(injectionFile,"%lg %lg %lg %lg %lg %lg %lg %lg",&f0,&dfdt,&theta,&phi,&amp,&iota,&psi,&phi0);
      N++;
    }
    rewind(injectionFile);
    N--;
    
    //set RNG for injection
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc(T);
    gsl_rng_env_setup();
    gsl_rng_set (r, data_vec[ii]->iseed);
    
    for(int nn=0; nn<N; nn++)
    {
      fscanf(injectionFile,"%lg %lg %lg %lg %lg %lg %lg %lg",&f0,&dfdt,&theta,&phi,&amp,&iota,&psi,&phi0);
      
      
      for(int jj=0; jj<flags->NT; jj++)
      {
        
        struct Data *data  = data_vec[ii];
        struct TDI *tdi = data->tdi[jj];
        
        
        //set bandwidth of data segment centered on injection
        if(nn==0)
        {
          data->fmin = f0 - (data->N/2)/data->T;
          data->fmax = f0 + (data->N/2)/data->T;
          data->qmin = (int)(data->fmin*data->T);
          data->qmax = data->qmin+data->N;
          
          //recompute fmin and fmax so they align with a bin
          data->fmin = data->qmin/data->T;
          data->fmax = data->qmax/data->T;
          
          if(jj==0)fprintf(stdout,"Frequency bins for segment [%i,%i]\n",data->qmin,data->qmax);
          fprintf(stdout,"   ...start time  %g\n",data->t0[jj]);
        }
        
        
        struct Source *inj = data->inj;
        
        for(int n=0; n<2*data->N; n++)
        {
          inj->tdi->A[n] = 0.0;
          inj->tdi->E[n] = 0.0;
          inj->tdi->X[n] = 0.0;
        }
        
        //map parameters to vector
        inj->f0       = f0;
        inj->dfdt     = dfdt;
        inj->costheta = cos(M_PI/2. - theta);
        inj->phi      = phi;
        inj->amp      = amp;
        inj->cosi     = cos(iota);
        inj->phi0     = phi0;
        inj->psi      = psi;
        if(data->NP>8)
          inj->d2fdt2 = 11.0/3.0*dfdt*dfdt/f0;
        
        map_params_to_array(inj, inj->params, data->T);
        
        //save parameters to file
        sprintf(filename,"injection_parameters_%i_%i.dat",ii,jj);
        if(nn==0)paramFile=fopen(filename,"w");
        else     paramFile=fopen(filename,"a");
        fprintf(paramFile,"%lg ",data->t0[jj]);
        print_source_params(data, inj, paramFile);
        fprintf(paramFile,"\n");
        fclose(paramFile);
        
        //Book-keeping of injection time-frequency volume
        galactic_binary_alignment(orbit, data, inj);
        
        printf("   ...bandwidth : %i\n",inj->BW);
        printf("   ...fdot      : %g\n",inj->dfdt*data->T*data->T);
        printf("   ...fddot     : %g\n",inj->d2fdt2*data->T*data->T*data->T);
        
        //Simulate gravitational wave signal
        //double t0 = data->t0 + jj*(data->T + data->tgap);
        galactic_binary(orbit, data->T, data->t0[jj], inj->params, data->NP, inj->tdi->X, inj->tdi->A, inj->tdi->E, inj->BW, 2);
        
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
        
        sprintf(filename,"waveform_injection_%i_%i.dat",ii,jj);
        fptr=fopen(filename,"w");
        for(int i=0; i<data->N; i++)
        {
          double f = (double)(i+data->qmin)/data->T;
          fprintf(fptr,"%lg %lg %lg %lg %lg",
                  f,
                  tdi->A[2*i],tdi->A[2*i+1],
                  tdi->E[2*i],tdi->E[2*i+1]);
          fprintf(fptr,"\n");
        }
        fclose(fptr);
        
        sprintf(filename,"power_injection_%i_%i.dat",ii,jj);
        fptr=fopen(filename,"w");
        for(int i=0; i<data->N; i++)
        {
          double f = (double)(i+data->qmin)/data->T;
          fprintf(fptr,"%.12g %lg %lg ",
                  f,
                  tdi->A[2*i]*tdi->A[2*i]+tdi->A[2*i+1]*tdi->A[2*i+1],
                  tdi->E[2*i]*tdi->E[2*i]+tdi->E[2*i+1]*tdi->E[2*i+1]);
          fprintf(fptr,"\n");
        }
        fclose(fptr);
        
        //Get noise spectrum for data segment
        for(int n=0; n<data->N; n++)
        {
          double f = data->fmin + (double)(n)/data->T;
          data->noise[jj]->SnA[n] = AEnoise(orbit->L, orbit->fstar, f);
          data->noise[jj]->SnE[n] = AEnoise(orbit->L, orbit->fstar, f);
        }
        
        //Get injected SNR
        fprintf(stdout,"   ...injected SNR=%g\n",snr(inj, data->noise[jj]));
        
        //Add Gaussian noise to injection
        gsl_rng_set (r, data->nseed+jj);
        
        if(!flags->zeroNoise && nn==0)
        {
          printf("   ...adding Gaussian noise realization\n");
          
          for(int n=0; n<data->N; n++)
          {
            tdi->A[2*n]   += gsl_ran_gaussian (r, 1)*sqrt(data->noise[jj]->SnA[n])/2.;
            tdi->A[2*n+1] += gsl_ran_gaussian (r, 1)*sqrt(data->noise[jj]->SnA[n])/2.;
            
            tdi->E[2*n]   += gsl_ran_gaussian (r, sqrt(data->noise[jj]->SnE[n])/2.);
            tdi->E[2*n+1] += gsl_ran_gaussian (r, sqrt(data->noise[jj]->SnE[n])/2.);
          }
        }
        
        //Compute fisher information matrix of injection
        printf("   ...computing Fisher Information Matrix of injection\n");
        
        galactic_binary_fisher(orbit, data, inj, data->noise[jj]);
        
        
        printf("\n Fisher Matrix:\n");
        for(int i=0; i<data->NP; i++)
        {
          fprintf(stdout," ");
          for(int j=0; j<data->NP; j++)
          {
            if(inj->fisher_matrix[i][j]<0)fprintf(stdout,"%.2e ", inj->fisher_matrix[i][j]);
            else                          fprintf(stdout,"+%.2e ",inj->fisher_matrix[i][j]);
          }
          fprintf(stdout,"\n");
        }
        
        printf("\n Fisher std. errors:\n");
        for(int j=0; j<data->NP; j++)  fprintf(stdout," %.4e\n", sqrt(inj->fisher_matrix[j][j]));
        
        
        
        sprintf(filename,"power_data_%i_%i.dat",ii,jj);
        fptr=fopen(filename,"w");
        
        for(int i=0; i<data->N; i++)
        {
          double f = (double)(i+data->qmin)/data->T;
          fprintf(fptr,"%.12g %lg %lg ",
                  f,
                  tdi->A[2*i]*tdi->A[2*i]+tdi->A[2*i+1]*tdi->A[2*i+1],
                  tdi->E[2*i]*tdi->E[2*i]+tdi->E[2*i+1]*tdi->E[2*i+1]);
          fprintf(fptr,"\n");
        }
        fclose(fptr);
      }//end jj loop over segments
    }//end nn loop over sources in file
    fclose(injectionFile);
    gsl_rng_free(r);
  }
  
  fprintf(stdout,"================================================\n\n");
}

void GalacticBinaryCatalogSNR(struct Data *data, struct Orbit *orbit, struct Flags *flags)
{
  fprintf(stdout,"\n==== GalacticBinaryInjectSimulatedSource ====\n");
  
  /* Get injection parameters */
  double f0,dfdt,theta,phi,amp,iota,psi,phi0; //read from injection file
  
  FILE *injectionFile = fopen(flags->injFile[0],"r");
  if(!injectionFile)
    fprintf(stderr,"Missing catalog file %s\n",flags->injFile[0]);
  else
    fprintf(stdout,"Simulateing binary catalog %s\n",flags->injFile[0]);
  
  //count sources in file
  int N=0;
  while(!feof(injectionFile))
  {
    fscanf(injectionFile,"%lg %lg %lg %lg %lg %lg %lg %lg",&f0,&dfdt,&theta,&phi,&amp,&iota,&psi,&phi0);
    N++;
  }
  rewind(injectionFile);
  N--;
  
  fprintf(stdout,"Found %i sources in %s\n",N,flags->injFile[0]);
  
  FILE *outfile = fopen("snr.dat","w");
  for(int n=0; n<N; n++)
  {
    
    fscanf(injectionFile,"%lg %lg %lg %lg %lg %lg %lg %lg",&f0,&dfdt,&theta,&phi,&amp,&iota,&psi,&phi0);
    
    //set bandwidth of data segment centered on injection
    data->fmin = f0 - (data->N/2)/data->T;
    data->fmax = f0 + (data->N/2)/data->T;
    data->qmin = (int)(data->fmin*data->T);
    data->qmax = data->qmin+data->N;
    
    struct Source *inj = data->inj;
    
    for(int n=0; n<2*data->N; n++)
    {
      inj->tdi->A[n] = 0.0;
      inj->tdi->E[n] = 0.0;
      inj->tdi->X[n] = 0.0;
    }
    
    //map parameters to vector
    inj->f0       = f0;
    inj->dfdt     = dfdt;
    inj->costheta = cos(M_PI/2. - theta);
    inj->phi      = phi;
    inj->amp      = amp;
    inj->cosi     = cos(iota);
    inj->phi0     = phi0;
    inj->psi      = psi;
    
    map_params_to_array(inj, inj->params, data->T);
    
    //Book-keeping of injection time-frequency volume
    galactic_binary_alignment(orbit, data, inj);
    
    //Simulate gravitational wave signal
    double t0 = data->t0[0];
    galactic_binary(orbit, data->T, t0, inj->params, 8, inj->tdi->X, inj->tdi->A, inj->tdi->E, inj->BW, 2);
    
    //Get noise spectrum for data segment
    for(int n=0; n<data->N; n++)
    {
      double f = data->fmin + (double)(n)/data->T;
      data->noise[0]->SnA[n] = AEnoise(orbit->L, orbit->fstar, f);
      data->noise[0]->SnE[n] = AEnoise(orbit->L, orbit->fstar, f);
    }
    
    //Get injected SNR
    double SNR = snr(inj, data->noise[0]);
    double Mc  = galactic_binary_Mc(f0, dfdt, data->T);
    double dL  = galactic_binary_dL(f0, dfdt, amp, data->T);
    
    fprintf(outfile,"%g %g %g %g %g %g %g %g %g\n",f0,dfdt,amp,cos(iota),Mc,dL,cos(M_PI/2 - theta),phi,SNR);
  }
  
  fclose(injectionFile);
  
  fprintf(stdout,"================================================\n\n");
}
