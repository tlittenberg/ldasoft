/*
*  Copyright (C) 2019 Neil J. Cornish, Tyson B. Littenberg (MSFC-ST12)
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <LISA.h>
#include <GalacticBinary.h>
#include <GalacticBinaryIO.h>
#include <GalacticBinaryWaveform.h>

#include "arrays.h"
#include "Constants.h"
#include "Detector.h"
#include "Subroutines.h"

int main(int argc,char **argv)
{
  
  double f, fdot, theta, phi, A, iota, psi, phase;
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
  struct Orbit *LISAorbit = malloc(sizeof(struct Orbit));

  //Set up orbit structure (allocate memory, read file, cubic spline)
  sprintf(LISAorbit->OrbitFileName,"%s",argv[2]);
  initialize_numeric_orbit(LISAorbit);

  XfLS = malloc(sizeof(double)*NFFT);
  AALS = malloc(sizeof(double)*NFFT);
  EELS = malloc(sizeof(double)*NFFT);
  
  
  for(i=0; i<NFFT; i++)
  {
    XfLS[i] = 0.0;
    AALS[i] = 0.0;
    EELS[i] = 0.0;
  }
  
  time_t rawtime;
  struct tm * timeinfo;
  
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  printf ( "Starting Simulation at: %s", asctime (timeinfo) );

  //count lines in file
  int NSIM = 0;
  int decade = 1;
  long time0 = time(0);
  while ( !feof(Infile) )
  {
    if(fscanf(Infile, "%lf%lf%lf%lf%lf%lf%lf%lf\n", &f, &fdot, &theta, &phi, &A, &iota, &psi, &phase)!=8) break;
    NSIM++;
    if(NSIM%decade==0)
    {
      decade*=10;
      fprintf(stdout,"\r read %i lines of file at t=%li s",NSIM, time(0)-time0);
      fflush(stdout);
    }
  }
  rewind(Infile);
  NSIM--;
  
  count=0;
  for(int n=0; n<NSIM; n++)
  {
    if(n%(NSIM/100)==0)printProgress((double)n/(double)NSIM);
    
    fscanf(Infile, "%lf%lf%lf%lf%lf%lf%lf%lf\n", &f, &fdot, &theta, &phi, &A, &iota, &psi, &phase);
    
    // hack for astrid simulation
    //theta-=0.5*M_PI;
    
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
    
    SAE  = AEnoise(LISAorbit->L,LISAorbit->fstar,f);
    SXYZ = XYZnoise(LISAorbit->L,LISAorbit->fstar,f);

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
    
    XLS = double_vector(2*M);
    AA  = double_vector(2*M);
    EE  = double_vector(2*M);
    
    galactic_binary(LISAorbit, "phase", TOBS, 0, params, 9, XLS, AA, EE, M, 2);
    
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
    
    
    free_double_vector(XLS);
    free_double_vector(AA);
    free_double_vector(EE);
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
    SAE  = AEnoise(LISAorbit->L,LISAorbit->fstar,f);
    SXYZ = XYZnoise(LISAorbit->L,LISAorbit->fstar,f);
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




