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

#include <LISA.h>
#include <GalacticBinary.h>
#include <GalacticBinaryWaveform.h>

#include "arrays.h"
#include "Constants.h"
#include "Detector.h"
#include "Subroutines.h"

void readGalaxyFile(const char *filename, int imax, double *XfLS, double *AALS, double *EELS);


int main(int argc,char **argv)
{
  
  double f, fdot, theta, phi, A, iota, psi, phase;
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
  
  printf("***********************************************************************\n");
  printf("*\n");
  printf("* FisherGalaxy: Source Subtraction Tool\n");
  printf("*   Simulated Data:      %s\n",argv[1]);
  printf("*   Confusion Noise Fit: %s\n",argv[2]);
  printf("*   Candidate Sources:   %s\n",argv[3]);
  printf("*   Orbit File:          %s\n",argv[4]);
  
  /* Figure out TOBS and NFFT */
  Infile = fopen(argv[1],"r");
  double junk;
  double f1,f2;
  fscanf(Infile,"%lf%lf%lf%lf%lf%lf%lf\n", &f1, &junk, &junk, &junk, &junk, &junk, &junk);
  fscanf(Infile,"%lf%lf%lf%lf%lf%lf%lf\n", &f2, &junk, &junk, &junk, &junk, &junk, &junk);
  double TOBS = 1./(f2-f1);
  int    NFFT = (int)floor(TOBS*DT);
  fclose(Infile);
  /*****************************/

  
  printf("*   Observing Time:      %.1f year (%f s)\n",TOBS/year,TOBS);
  printf("*\n");
  printf("***********************************************************************\n");
  
  Xbright = fopen("BrightX.dat","w");
  Abright = fopen("BrightAE.dat","w");
  
  //Data structure for interpolating orbits from file
  struct Orbit *LISAorbit = malloc(sizeof(struct Orbit));
  
  
  //Set up orbit structure (allocate memory, read file, cubic spline)
  sprintf(LISAorbit->OrbitFileName,"%s",argv[4]);
  initialize_numeric_orbit(LISAorbit);
  double L     = LISAorbit->L;
  double fstar = LISAorbit->fstar;
  
  params = double_vector(9);
  
  SNRthres = 7.0;
  
  if((TOBS/year) <= 8.0) mult = 8;
  if((TOBS/year) <= 4.0) mult = 4;
  if((TOBS/year) <= 2.0) mult = 2;
  if((TOBS/year) <= 1.0) mult = 1;
  
  XfLS = double_vector(NFFT-1);
  AALS = double_vector(NFFT-1);
  EELS = double_vector(NFFT-1);
  
  for(i=0; i<NFFT; i++)
  {
    XfLS[i] = 0.0;
    AALS[i] = 0.0;
    EELS[i] = 0.0;
  }
  
  imax = (long)ceil(4.0e-2*TOBS);
  imin = (long)floor(1.0e-4*TOBS);
  sqT = sqrt(TOBS);
  
  readGalaxyFile(argv[1],imax,XfLS,AALS,EELS);
  
  XP = double_vector(NFFT/2);  AEP = double_vector(NFFT/2);
  Xnoise = double_vector(NFFT/2);  Xconf = double_vector(NFFT/2);
  AEnoise = double_vector(NFFT/2);  AEconf = double_vector(NFFT/2);
  
  for(i=0; i< NFFT/2; i++)
  {
    XP[i]  = (2.0*(XfLS[2*i]*XfLS[2*i] + XfLS[2*i+1]*XfLS[2*i+1]));
    AEP[i] = (2.0*(AALS[2*i]*AALS[2*i]+AALS[2*i+1]*AALS[2*i+1]));
  }
  
  printf("Reading Confusion Noise File\n");
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
  
  //count lines in file
  int NSIM = 0;
  while ( !feof(Infile) )
  {
    fscanf(Infile, "%lf%lf%lf%lf%lf%lf%lf%lf\n", &f, &fdot, &theta, &phi, &A, &iota, &psi, &phase);
    NSIM++;
  }
  rewind(Infile);
  NSIM--;
  
  for(int n=0; n<NSIM; n++)
  {
    if(n%(NSIM/100)==0)printProgress((double)n/(double)NSIM);
    
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
    Sm = SXYZ/(4.0*sin(f/fstar)*sin(f/fstar));
    
    Acut = A*sqrt(TOBS/Sm);
    
    M = galactic_binary_bandwidth(LISAorbit->L, LISAorbit->fstar, f, fdot, cos(params[1]), params[3], TOBS, N);
    
    XLS = double_vector(2*M);
    AA  = double_vector(2*M);
    EE  = double_vector(2*M);
    
    galactic_binary(LISAorbit, "phase", TOBS, 0, params, 9, XLS, AA, EE, M, 2);
    
    /*inst2*/
    SNX = (SXYZ+Xconf[q]);//*sin(f/fstar)*sin(f/fstar);
    SNAE = (SAE+AEconf[q]);//*sin(f/fstar)*sin(f/fstar);
    
    SNRAE = 0.0;
    SNRX = 0.0;
    for(i=1; i<=M; i++)
    {
      SNRX += 4.0*(XLS[2*i-1]*XLS[2*i-1]+XLS[2*i]*XLS[2*i]);
      SNRAE += 4.0*(AA[2*i-1]*AA[2*i-1]+AA[2*i]*AA[2*i]+EE[2*i-1]*EE[2*i-1]+EE[2*i]*EE[2*i]);
    }
    SNRAE *= TOBS/SNAE;
    SNRX  *= TOBS/SNX;
    SNRAE  = sqrt(SNRAE);
    SNRX   = sqrt(SNRX);
    
    if(SNRX > SNRthres)
    {
      fprintf(Xbright, "%.16f %.10e %f %f %e %f %f %f\n", f, fdot, theta, phi, A, iota, psi, phase);
      for(i=1; i<=M; i++)
      {
        k = (q + i -1 - M/2);
        if(k>0)
        {
          XfLS[2*k]   -= sqT*XLS[2*i-1];
          XfLS[2*k+1] -= sqT*XLS[2*i];
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
          AALS[2*k]   -= sqT*AA[2*i-1];
          AALS[2*k+1] -= sqT*AA[2*i];
          EELS[2*k]   -= sqT*EE[2*i-1];
          EELS[2*k+1] -= sqT*EE[2*i];
        }
      }
    }
    free_double_vector(XLS);
    free_double_vector(AA);
    free_double_vector(EE);
  }
  printProgress(1.0);

  printf("\nRemoval Finished\n");
  
  
  Outfile = fopen("Galaxy_XAE_R1.dat","w");
  for(i=1; i< imax; i++)
  {
    f = (double)(i)/TOBS;
    fprintf(Outfile,"%.12g %e %e %e %e %e %e\n", f, XfLS[2*i], XfLS[2*i+1],
            AALS[2*i], AALS[2*i+1], EELS[2*i], EELS[2*i+1]);
  }
  fclose(Outfile);
  
  
  for(i=1; i< NFFT/2; i++)
  {
    XP[i] = (2.0*(XfLS[2*i]*XfLS[2*i] + XfLS[2*i+1]*XfLS[2*i+1]));
    AEP[i] = (2.0*(AALS[2*i]*AALS[2*i]+AALS[2*i+1]*AALS[2*i+1]));
  }
  
  for(i=1; i< NFFT/2; i++)
  {
    XP[i]  = (2.0*(XfLS[2*i]*XfLS[2*i] + XfLS[2*i+1]*XfLS[2*i+1]));
    AEP[i] = (2.0*(AALS[2*i]*AALS[2*i]+AALS[2*i+1]*AALS[2*i+1]));
    instrument_noise((double)i/TOBS, LISAorbit->fstar, LISAorbit->L, &AEnoise[i], &Xnoise[i]);
  }
  
  printf("Estimating Confusion Noise\n");
  
  int divs = 100;  // must be even - used to compute median
  
  if(divs/2+1 > imin) imin = divs/2+1;
  if(imax > NFFT/2-divs/2-1) imax =  NFFT/2-divs/2-1;

  //spline_fit(0, divs, imin, imax, XP, Xnoise, Xconf, TOBS, fstar, L);
  //spline_fit(1, divs, imin, imax, AEP, AEnoise, AEconf, TOBS, fstar, L);
  //confusion_mcmc(AEP, AEnoise, AEconf, imin, imax, TOBS);
  confusion_mcmc(AEP, AEnoise, AEconf, (int)floor(0.0001*TOBS), (int)floor(0.006*TOBS), TOBS);

  
  medianX(imin, imax, fstar, L, XP, Xnoise, Xconf, TOBS);
  //medianAE(imin, imax, fstar, L, AEP, AEnoise, AEconf, TOBS);
  
  
  Outfile = fopen("Confusion_XAE_1.dat","w");
  for(i=imin; i<= imax; i++)
  {
    f = (double)(i)/TOBS;
    fprintf(Outfile,"%.12g %e %e %e %e\n", f, Xnoise[i], Xconf[i], AEnoise[i], AEconf[i]);
  }
  fclose(Outfile);
  
  Outfile = fopen("Confusion_XAE_DS.dat","w");
  for(i=imin; i<= imax; i++)
  {
    if(i%100==0)
    {
      f = (double)(i)/TOBS;
      fprintf(Outfile,"%.12g %e %e %e %e\n", f, Xnoise[i], Xconf[i], AEnoise[i], AEconf[i]);
    }
  }
  fclose(Outfile);
  
  return 0;
  
}


void readGalaxyFile(const char *filename, int imax, double *XfLS, double *AALS, double *EELS)
{
  printf("Reading Galaxy File\n");
  FILE *Infile = fopen(filename,"r");
  
  double f;
  for(int i=1; i< imax; i++)
  {
    if(i%(imax/100)==0)printProgress((double)i/(double)imax);
    fscanf(Infile,"%lf%lf%lf%lf%lf%lf%lf\n", &f, &XfLS[2*i], &XfLS[2*i+1],
           &AALS[2*i], &AALS[2*i+1], &EELS[2*i], &EELS[2*i+1]);
  }
  printProgress(1.0);
  printf("\n");
  fclose(Infile);
}
