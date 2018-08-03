/*********************************************************/
/*                                                       */
/*           Galaxy.c, Version 2.3, 4/28/2011            */
/*      Written by Neil Cornish & Tyson Littenberg       */
/*                                                       */
/*        gcc -O2 -o Galaxy Galaxy.c arrays.c -lm        */
/*                                                       */
/*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "arrays.h"
#include "Constants.h"
#include "Detector.h"
#include "Subroutines.h"

#define C 299792458.
#define TSUN  4.9169e-6
#define PC 3.0856775807e16

double M_fdot(double f, double fdot)
{
  return pow( fdot*pow(f,-11./3.)*(5./96.)*pow(M_PI,-8./3.)  ,  3./5.)/TSUN;
}

double galactic_binary_dL(double f0, double dfdt, double A)
{
  double f    = f0;//T;
  double fd = dfdt;//(T*T);
  double amp   = A;
  return ((5./48.)*(fd/(M_PI*M_PI*f*f*f*amp))*C/PC); //seconds  !check notes on 02/28!
}

int main(int argc,char **argv)
{
  
  double f, fdot, theta, phi, A, iota, psi, phase;
  char Gfile[50];
  double *params;
  double *XfLS, *AALS, *EELS;
  double *XLS, *AA, *EE;
  double fonfs, Sm, Acut;
  long M, N, q;
  long i, k, cnt, cnt2, mult, imax;
  long rseed;
  double SAE, SXYZ, sqT;
  double XR, XI, AR, AI, ER, EI;
  double alpha;
  
  FILE* Infile;
  FILE* Outfile;
  
  if(argc != 3) KILL("Galaxy Galaxy.dat Orbits.dat\n");
  
  params = dvector(0,9);
  
  if((TOBS/year) <= 8.0) mult = 8;
  if((TOBS/year) <= 4.0) mult = 4;
  if((TOBS/year) <= 2.0) mult = 2;
  if((TOBS/year) <= 1.0) mult = 1;
  
  Infile  = fopen(argv[1],"r");
  Outfile = fopen("Bright.dat","w");
  
  //Data structure for interpolating orbits from file
  struct lisa_orbit *LISAorbit;
  LISAorbit = &orbit;
  
  //Set up orbit structure (allocate memory, read file, cubic spline)
  sprintf(Gfile,"%s",argv[2]);
  initialize_orbit(Gfile, LISAorbit);
		
  XfLS = dvector(0,NFFT-1);  AALS = dvector(0,NFFT-1);  EELS = dvector(0,NFFT-1);
  
  
  for(i=0; i<NFFT; i++)
  {
    XfLS[i] = 0.0;
    AALS[i] = 0.0;
    EELS[i] = 0.0;
  }
  
  printf("Starting Simulation\n");
  
  cnt = 0;
  cnt2 = 0;
  
  rseed = -924514346;
  
  double Amax=0;
  
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
    
    alpha = ran2(&rseed);
    
    if((fdot > 0.0) || (alpha < 0.1))  // keep all detached and 10% of AMCVn
    {
      
      cnt++;
      
      if (cnt%100000 == 0) printf("%ld sources simulated, %ld bright\n", cnt, cnt2);
      
      N = 32*mult;
      if(f > 0.001) N = 64*mult;
      if(f > 0.01) N = 256*mult;
      if(f > 0.03) N = 512*mult;
      if(f > 0.1) N = 1024*mult;
      
      
      fonfs = f/LISAorbit->fstar;
      
      q = (long)(f*TOBS);
      
      instrument_noise(f, LISAorbit->fstar, LISAorbit->L, &SAE, &SXYZ);
      
      /*  calculate michelson noise  */
      
      /*inst2*/
      if(noiseFlag==1)Sm = SXYZ/(4.0*sin(fonfs)*sin(fonfs));
      if(noiseFlag==2)Sm = SXYZ/(4.0);//*sin(fonfs)*sin(fonfs));
      
      Acut = A*sqrt(TOBS/Sm);
      if(f>0.008)
      {
        if(Acut>Amax)
        {
          Amax=Acut;
          fprintf(stdout,"New max: f=%g, Acut=%g, D=%g, M=%g\n",f,Acut,galactic_binary_dL(f, fdot, A),M_fdot(f,fdot));
        }
      }
      if(Acut > 2.0)
      {
        cnt2++;
        fprintf(Outfile, "%.16f %.10e %f %f %e %f %f %f\n", f, fdot, theta, phi, A, iota, psi, phase);
      }
      
      if(f < 4.0e-2)
      {
        
        M = (long)(pow(2.0,(rint(log(Acut)/log(2.0))+1.0)));
        
        if(M < N) M = N;
        if(N < M) N = M;
        if(M > 8192) M = 8192;
        //N=M=2*8192;
        N=M;
        
        XLS = dvector(1,2*M);
        AA = dvector(1,2*M);   EE = dvector(1,2*M);
        
        //added orbit structure to FAST_LISA function call -- TBL
        FAST_LISA(LISAorbit, params, N, M, XLS, AA, EE);
        
        for(i=1; i<=M; i++)
        {
          
          k = (q + i -1 - M/2);
          
          //if(2*k+1 > NFFT) printf("bad\n");
          if(k>0)
          {
            XfLS[2*k] += XLS[2*i-1];
            XfLS[2*k+1] += XLS[2*i];
            AALS[2*k] += AA[2*i-1];
            AALS[2*k+1] += AA[2*i];
            EELS[2*k] += EE[2*i-1];
            EELS[2*k+1] += EE[2*i];
          }
          
        }
        
        
        free_dvector(XLS,1,2*M);  free_dvector(AA,1,2*M);  free_dvector(EE,1,2*M);
        
      }
      
    }
    
  }
  
  printf("Simulation Finished\n");
  
  imax = (long)ceil(4.0e-2*TOBS);
  sqT = sqrt(TOBS);
  
  rseed = -7584529636;
  
  Outfile = fopen("Galaxy_XAE.dat","w");
  for(i=1; i< imax; i++)
  {
    f = (double)(i)/TOBS;
    fonfs = f/LISAorbit->fstar;
    instrument_noise(f, LISAorbit->fstar, LISAorbit->L, &SAE, &SXYZ);
    XR = 0.5*sqrt(SXYZ) * gasdev2(&rseed);
    XI = 0.5*sqrt(SXYZ) * gasdev2(&rseed);
    AR = 0.5*sqrt(SAE) * gasdev2(&rseed);
    AI = 0.5*sqrt(SAE) * gasdev2(&rseed);
    ER = 0.5*sqrt(SAE) * gasdev2(&rseed);
    EI = 0.5*sqrt(SAE) * gasdev2(&rseed);
    if(noiseFlag==1)fprintf(Outfile,"%e %e %e %e %e %e %e\n", f, sqT*XfLS[2*i]+XR, sqT*XfLS[2*i+1]+XI,
                            sqT*AALS[2*i]+AR, sqT*AALS[2*i+1]+AI, sqT*EELS[2*i]+ER, sqT*EELS[2*i+1]+EI);
    if(noiseFlag==2)fprintf(Outfile,"%e %e %e %e %e %e %e\n", f, (sqT/sin(fonfs))*XfLS[2*i]+XR, (sqT/sin(fonfs))*XfLS[2*i+1]+XI,
                            (sqT/sin(fonfs))*AALS[2*i]+AR, (sqT/sin(fonfs))*AALS[2*i+1]+AI, (sqT/sin(fonfs))*EELS[2*i]+ER, (sqT/sin(fonfs))*EELS[2*i+1]+EI);
  }
  fclose(Outfile);
  
  return 0;
  
}




