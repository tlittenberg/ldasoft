/***********************************************************************/
/*                                                                     */
/*                  Galaxy.c, Version 3.0, 8/03/2018                   */
/*             Written by Neil Cornish & Tyson Littenberg              */
/*                                                                     */
/*        gcc -O2 -o Galaxy Galaxy.c Subroutines.c arrays.c -lm        */
/*                                                                     */
/***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

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
  long i, k, count, mult, imax;
  double SAE, SXYZ, sqT;
  double XR, XI, AR, AI, ER, EI;
  
  FILE* Infile;
  FILE* Outfile;
  
  if(argc != 4) KILL("Galaxy Galaxy.dat Orbits.dat TOBS\n");
 
  printf("***********************************************************************\n");
  printf("*\n");
  printf("* FisherGalaxy: Galaxy Simulation Tool\n");
  printf("*   Galaxy Simulation: %s\n",argv[1]);
  printf("*   Orbit File:        %s\n",argv[2]);
  printf("*   Observing Time:    %.1f year\n",atof(argv[3])/year);
  printf("*\n");
  printf("***********************************************************************\n");
  
  double TOBS = (double)atof(argv[3]);
  int    NFFT = (int)floor(TOBS*DT);
  
  //set RNG for noise
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *r = gsl_rng_alloc(T);
  gsl_rng_env_setup();
  gsl_rng_set (r, -924514346);

  
  params = malloc(sizeof(double)*9);
  
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
		
  XfLS = malloc(sizeof(double)*NFFT);
  AALS = malloc(sizeof(double)*NFFT);
  EELS = malloc(sizeof(double)*NFFT);
  
  
  for(i=0; i<NFFT; i++)
  {
    XfLS[i] = 0.0;
    AALS[i] = 0.0;
    EELS[i] = 0.0;
  }
  
  printf("Starting Simulation\n");
  
  //count lines in file
  int NSIM = 0;
  int decade = 1;
  while ( !feof(Infile) )
  {
    fscanf(Infile, "%lf%lf%lf%lf%lf%lf%lf%lf\n", &f, &fdot, &theta, &phi, &A, &iota, &psi, &phase);
    NSIM++;
    if(NSIM%decade==0)
    {
      decade*=10;
      printf("read %i lines of file\n");
    }
  }
  rewind(Infile);
  NSIM--;
  
  count=0;
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
    
    
    fonfs = f/LISAorbit->fstar;
    
    q = (long)(f*TOBS);
    
    instrument_noise(f, LISAorbit->fstar, LISAorbit->L, &SAE, &SXYZ);
    
    /*  calculate michelson noise  */
    Sm = SXYZ/(4.0*sin(fonfs)*sin(fonfs));

    /* rough guess at SNR */
    Acut = A*sqrt(TOBS/Sm);
    if(Acut > 2.0)
    {
      count++;
      fprintf(Outfile, "%.16f %.10e %f %f %e %f %f %f\n", f, fdot, theta, phi, A, iota, psi, phase);
    }
    
    M = galactic_binary_bandwidth(LISAorbit->L, LISAorbit->fstar, f, fdot, cos(params[1]), params[3], TOBS, N);
    
    XLS = dvector(1,2*M);
    AA  = dvector(1,2*M);
    EE  = dvector(1,2*M);
    
    FAST_LISA(LISAorbit, TOBS, params, N, M, XLS, AA, EE);
    
    for(i=1; i<=M; i++)
    {
      
      k = (q + i - 1 - M/2);
      
      if(k>0 && k<NFFT/2)
      {
        XfLS[2*k]   += XLS[2*i-1];
        XfLS[2*k+1] += XLS[2*i];
        AALS[2*k]   += AA[2*i-1];
        AALS[2*k+1] += AA[2*i];
        EELS[2*k]   += EE[2*i-1];
        EELS[2*k+1] += EE[2*i];
      }
      
    }
    
    
    free_dvector(XLS,1,2*M);
    free_dvector(AA,1,2*M);
    free_dvector(EE,1,2*M);
  }
  fclose(Infile);
  printProgress(1.0);

  printf("\nSimulation Finished\n");
  
  imax = (long)ceil(4.0e-2*TOBS);
  sqT = sqrt(TOBS);
  
  Outfile = fopen("Galaxy_XAE.dat","w");
  for(i=1; i< imax; i++)
  {
    f = (double)(i)/TOBS;
    fonfs = f/LISAorbit->fstar;
    instrument_noise(f, LISAorbit->fstar, LISAorbit->L, &SAE, &SXYZ);
    XR = 0.5 * sqrt(SXYZ) * gsl_ran_ugaussian(r);
    XI = 0.5 * sqrt(SXYZ) * gsl_ran_ugaussian(r);
    AR = 0.5 * sqrt(SAE)  * gsl_ran_ugaussian(r);
    AI = 0.5 * sqrt(SAE)  * gsl_ran_ugaussian(r);
    ER = 0.5 * sqrt(SAE)  * gsl_ran_ugaussian(r);
    EI = 0.5 * sqrt(SAE)  * gsl_ran_ugaussian(r);
    fprintf(Outfile,"%.12g %.12g %.12g %.12g %.12g %.12g %.12g\n", f, sqT*XfLS[2*i]+XR, sqT*XfLS[2*i+1]+XI, sqT*AALS[2*i]+AR, sqT*AALS[2*i+1]+AI, sqT*EELS[2*i]+ER, sqT*EELS[2*i+1]+EI);
  }
  fclose(Outfile);
  
  
  
  
  free(XfLS);
  free(AALS);
  free(EELS);

  return 0;
  
}




