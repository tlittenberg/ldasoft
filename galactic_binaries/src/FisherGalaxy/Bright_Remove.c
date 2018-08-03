/*********************************************************/
/*                                                       */
/*        Bright_Remove.c, Version 2.3, 4/28/2011        */
/*      Written by Neil Cornish & Tyson Littenberg       */
/*                                                       */
/* gcc -O2 -o Bright_Remove Bright_Remove.c arrays.c -lm */
/*                                                       */
/*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "arrays.h"
#include "Constants.h"
#include "Detector.h"
#include "Subroutines.h"


int main(int argc,char **argv)
{
  
  double f, fdot, theta, phi, A, iota, psi, phase;
  char Gfile[50];
  double *params;
  double *XfLS, *AALS, *EELS;
  double *XLS, *AA, *EE;
  double fonfs, Sm, Acut;
  long M, N, q;
  long i, k, cnt, cc1, mult, imax, imin;
  double SAE, SXYZ, sqT;
  double *Xnoise, *Xconf;
  double *AEnoise, *AEconf;
  double *XP, *AEP;
  double SNX, SNAE, SNRAE, SNRX;
  double SNRthres;
  
  FILE* Infile;
  FILE* Outfile;
  FILE* Xbright;
  FILE* Abright;
  
  if(argc !=5) KILL("Bright_Remove XAE.dat Noise.dat Bright.dat Orbit.dat\n");
  
  Xbright = fopen("BrightX.dat","w");
  Abright = fopen("BrightAE.dat","w");
  
  //Data structure for interpolating orbits from file
  struct lisa_orbit *LISAorbit;
  LISAorbit = &orbit;
  
  //Set up orbit structure (allocate memory, read file, cubic spline)
  sprintf(Gfile,"%s",argv[4]);
  initialize_orbit(Gfile, LISAorbit);
  double L     = LISAorbit->L;
  double fstar = LISAorbit->fstar;
  
  params = dvector(0,9);
  
  SNRthres = 7.0;
  
  if((TOBS/year) <= 8.0) mult = 8;
  if((TOBS/year) <= 4.0) mult = 4;
  if((TOBS/year) <= 2.0) mult = 2;
  if((TOBS/year) <= 1.0) mult = 1;
  
  XfLS = dvector(0,NFFT-1);  AALS = dvector(0,NFFT-1);  EELS = dvector(0,NFFT-1);
  
  for(i=0; i<NFFT; i++)
  {
    XfLS[i] = 0.0;
    AALS[i] = 0.0;
    EELS[i] = 0.0;
  }
  
  imax = (long)ceil(4.0e-2*TOBS);
  imin = (long)floor(1.0e-4*TOBS);
  sqT = sqrt(TOBS);
  
  Infile = fopen(argv[1],"r");
  for(i=1; i< imax; i++)
  {
    fscanf(Infile,"%lf%lf%lf%lf%lf%lf%lf\n", &f, &XfLS[2*i], &XfLS[2*i+1],
           &AALS[2*i], &AALS[2*i+1], &EELS[2*i], &EELS[2*i+1]);
  }
  fclose(Infile);
  
  XP = dvector(imin,imax);  AEP = dvector(imin,imax);
  Xnoise = dvector(imin,imax);  Xconf = dvector(imin,imax);
  AEnoise = dvector(imin,imax);  AEconf = dvector(imin,imax);
  
  Outfile = fopen("Power_0.dat","w");
  for(i=imin; i< imax; i++)
  {
    f = (double)(i)/TOBS;
    instrument_noise(f, fstar, L ,&SAE, &SXYZ);
    XP[i] = (2.0*(XfLS[2*i]*XfLS[2*i] + XfLS[2*i+1]*XfLS[2*i+1]));
    AEP[i] = (2.0*(AALS[2*i]*AALS[2*i]+AALS[2*i+1]*AALS[2*i+1]));
    fprintf(Outfile,"%e %e %e %e %e\n", f, XP[i], SXYZ, AEP[i], SAE);
  }
  fclose(Outfile);
  
  Outfile = fopen(argv[2],"r");
  for(i=imin; i<= imax; i++)
  {
    fscanf(Outfile,"%lf%lf%lf%lf%lf\n", &f, &Xnoise[i], &Xconf[i], &AEnoise[i], &AEconf[i]);
  }
  fclose(Outfile);
  
  Infile = fopen(argv[3],"r");
  
  printf("Starting Removal\n");
  
  cnt = 0;
  cc1 = 0;
  
  while ( !feof(Infile) )
  {
    fscanf(Infile, "%lf%lf%lf%lf%lf%lf%lf%lf\n", &f, &fdot, &theta, &phi, &A, &iota, &psi, &phase);
    
    params[0] = f;
    params[1] = 0.5*pi-theta;
    params[2] = phi;
    params[3] = A;
    params[4] = iota;
    params[5] = psi;
    params[6] = phase;
    params[7] = fdot;
    params[8] = 11.0/3.0*fdot*fdot/f;
    
    if(f > 4.0e-2) cc1++;
    
    if(f < 4.0e-2)
    {
      
      N = 32*mult;
      if(f > 0.001) N = 64*mult;
      if(f > 0.01) N = 256*mult;
      if(f > 0.03) N = 512*mult;
      if(f > 0.1) N = 1024*mult;
      
      
      fonfs = f/fstar;
      
      q = (long)(f*TOBS);
      
      instrument_noise(f, fstar, L, &SAE, &SXYZ);
      
      /*  calculate michelson noise  */
						
      /*inst2*/
      if(noiseFlag==1) Sm = SXYZ/(4.0*sin(f/fstar)*sin(f/fstar));
      if(noiseFlag==2) Sm = SXYZ/(4.0);//*sin(f/fstar)*sin(f/fstar));
      
      Acut = A*sqrt(TOBS/Sm);
      
      M = (long)(pow(2.0,(rint(log(Acut)/log(2.0))+1.0)));
      
      if(M < N) M = N;
      if(N < M) N = M;
      if(M > 8192) M = 8192;
      
      N = M;
      //N = M = 2*8192;
      
      XLS = dvector(1,2*M);
      AA = dvector(1,2*M);   EE = dvector(1,2*M);
      
      FAST_LISA(LISAorbit, params, N, M, XLS, AA, EE);
      
      /*inst2*/
      if(noiseFlag==1)
      {
        SNX = (SXYZ+Xconf[q]);//*sin(f/fstar)*sin(f/fstar);
        SNAE = (SAE+AEconf[q]);//*sin(f/fstar)*sin(f/fstar);
      }
      if(noiseFlag==2)
      {
        SNX = (SXYZ+Xconf[q])*sin(f/fstar)*sin(f/fstar);
        SNAE = (SAE+AEconf[q])*sin(f/fstar)*sin(f/fstar);
      }
      SNRAE = 0.0;
      SNRX = 0.0;
      for(i=1; i<=M; i++)
      {
        SNRX += 4.0*(XLS[2*i-1]*XLS[2*i-1]+XLS[2*i]*XLS[2*i]);
        SNRAE += 4.0*(AA[2*i-1]*AA[2*i-1]+AA[2*i]*AA[2*i]+EE[2*i-1]*EE[2*i-1]+EE[2*i]*EE[2*i]);
      }
      SNRAE *= TOBS/SNAE;
      SNRX *= TOBS/SNX;
      SNRAE = sqrt(SNRAE);
      SNRX = sqrt(SNRX);
      
      if(SNRX > SNRthres)
      {
        fprintf(Xbright, "%.16f %.10e %f %f %e %f %f %f\n", f, fdot, theta, phi, A, iota, psi, phase);
        for(i=1; i<=M; i++)
        {
          k = (q + i -1 - M/2);
          if(k>0)
          {
            if(noiseFlag==1)
            {
              XfLS[2*k] -= sqT*XLS[2*i-1];
              XfLS[2*k+1] -= sqT*XLS[2*i];
            }
            if(noiseFlag==2)
            {
              XfLS[2*k] -= sqT*XLS[2*i-1]/sin(f/fstar);
              XfLS[2*k+1] -= sqT*XLS[2*i]/sin(f/fstar);
            }
          }
        }
      }
      
      if(SNRAE > SNRthres)
      {
        fprintf(Abright, "%.16f %.10e %f %f %e %f %f %f\n", f, fdot, theta, phi, A, iota, psi, phase);
        for(i=1; i<=M; i++)
        {
          k = (q + i -1 - M/2);
          if(k>0)
          {
            if(noiseFlag==1)
            {
              AALS[2*k] -= sqT*AA[2*i-1];
              AALS[2*k+1] -= sqT*AA[2*i];
              EELS[2*k] -= sqT*EE[2*i-1];
              EELS[2*k+1] -= sqT*EE[2*i];
            }
            if(noiseFlag==2)
            {
              AALS[2*k] -= sqT*AA[2*i-1]/sin(f/fstar);
              AALS[2*k+1] -= sqT*AA[2*i]/sin(f/fstar);
              EELS[2*k] -= sqT*EE[2*i-1]/sin(f/fstar);
              EELS[2*k+1] -= sqT*EE[2*i]/sin(f/fstar);
            }
          }
          
        }
      }
      
      
      free_dvector(XLS,1,2*M);  free_dvector(AA,1,2*M);  free_dvector(EE,1,2*M);
      
    }
    else
    {
      // all the really high f sources will be detectable
      fprintf(Abright, "%.16f %.10e %f %f %e %f %f %f\n", f, fdot, theta, phi, A, iota, psi, phase);
      fprintf(Xbright, "%.16f %.10e %f %f %e %f %f %f\n", f, fdot, theta, phi, A, iota, psi, phase);
    }
    
    
    
  }
  
  printf("Removal Finished\n");
  
  printf("Number above 1e-2 Hz = %ld\n", cc1);
  
  
  Outfile = fopen("Galaxy_XAE_R1.dat","w");
  for(i=1; i< imax; i++)
  {
    f = (double)(i)/TOBS;
    fprintf(Outfile,"%e %e %e %e %e %e %e\n", f, XfLS[2*i], XfLS[2*i+1],
            AALS[2*i], AALS[2*i+1], EELS[2*i], EELS[2*i+1]);
  }
  fclose(Outfile);
  
  
  for(i=imin; i< imax; i++)
  {
    XP[i] = (2.0*(XfLS[2*i]*XfLS[2*i] + XfLS[2*i+1]*XfLS[2*i+1]));
    AEP[i] = (2.0*(AALS[2*i]*AALS[2*i]+AALS[2*i+1]*AALS[2*i+1]));
  }
  
  Outfile = fopen("Power_1.dat","w");
  for(i=imin; i< imax; i++)
  {
    f = (double)(i)/TOBS;
    instrument_noise(f, fstar, L, &SAE, &SXYZ);
    XP[i] = (2.0*(XfLS[2*i]*XfLS[2*i] + XfLS[2*i+1]*XfLS[2*i+1]));
    AEP[i] = (2.0*(AALS[2*i]*AALS[2*i]+AALS[2*i+1]*AALS[2*i+1]));
    fprintf(Outfile,"%e %e %e %e %e\n", f, XP[i], SXYZ, AEP[i], SAE);
  }
  fclose(Outfile);
  
  medianX(imin, imax, fstar, L, XP, Xnoise, Xconf);
  medianAE(imin, imax, fstar, L, AEP, AEnoise, AEconf);
  
  Outfile = fopen("Confusion_XAE_1.dat","w");
  for(i=imin; i<= imax; i++)
  {
    f = (double)(i)/TOBS;
    fprintf(Outfile,"%e %e %e %e %e\n", f, Xnoise[i], Xconf[i], AEnoise[i], AEconf[i]);
  }
  fclose(Outfile);
  
  Outfile = fopen("Confusion_XAE_DS.dat","w");
  for(i=imin; i<= imax; i++)
  {
    if(i%100==0)
    {
      f = (double)(i)/TOBS;
      fprintf(Outfile,"%e %e %e %e %e\n", f, Xnoise[i], Xconf[i], AEnoise[i], AEconf[i]);
    }
  }
  fclose(Outfile);
  
  return 0;
  
}


