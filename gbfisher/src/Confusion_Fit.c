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

#include "arrays.h"
#include "Constants.h"
#include "Detector.h"
#include "Subroutines.h"


int main(int argc,char **argv)
{
  
  double f;
  double *XfLS, *AALS, *EELS;
  double *XP, *AEP;
  long i, imax, imin;
  long rseed;
  double *Xnoise, *Xconf;
  double *AEnoise, *AEconf;
  
  FILE* Infile;
  FILE* Outfile;
  
  if(argc !=3) KILL("Confusion_Fit Galaxy.dat Orbit.dat\n");
  
  printf("***********************************************************************\n");
  printf("*\n");
  printf("* FisherGalaxy: Confusion Noise Fit Tool\n");
  printf("*   Simulated Data: %s\n",argv[1]);
  printf("*   Orbit File:     %s\n",argv[2]);
  
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

  printf("*   Observing Time: %.1f year (%f s)\n",TOBS/YEAR,TOBS);
  printf("*\n");
  printf("***********************************************************************\n");

  XfLS = double_vector(NFFT-1);  AALS = double_vector(NFFT-1);  EELS = double_vector(NFFT-1);
  
  
  imax = (long)ceil(4.0e-2*TOBS);
  imin = (long)floor(1.0e-4*TOBS);
  
  XfLS = double_vector(NFFT-1);  AALS = double_vector(NFFT-1); EELS = double_vector(NFFT-1);
  
  
  //Data structure for interpolating orbits from file
  struct Orbit *LISAorbit = malloc(sizeof(struct Orbit));

  //Set up orbit structure (allocate memory, read file, cubic spline)
  sprintf(LISAorbit->OrbitFileName,"%s",argv[2]);
  initialize_numeric_orbit(LISAorbit);

  double L     = LISAorbit->L;
  double fstar = LISAorbit->fstar;
  
  rseed = -7584529636;
  
  printf("Reading Data File\n");
  
  Infile = fopen(argv[1],"r");
  for(i=1; i< imax; i++)
  {
    fscanf(Infile,"%lf%lf%lf%lf%lf%lf%lf\n", &f, &XfLS[2*i], &XfLS[2*i+1],
           &AALS[2*i], &AALS[2*i+1], &EELS[2*i], &EELS[2*i+1]);
  }
  fclose(Infile);
  
  printf("Estimating Confusion Noise\n");

  XP = double_vector(NFFT/2);  AEP = double_vector(NFFT/2);
  Xnoise = double_vector(NFFT/2);  Xconf = double_vector(NFFT/2);
  AEnoise = double_vector(NFFT/2);  AEconf = double_vector(NFFT/2);
  
  rseed = -7584529636;
  
  for(i=0; i< NFFT/2; i++)
  {
    XP[i]  = (2.0*(XfLS[2*i]*XfLS[2*i] + XfLS[2*i+1]*XfLS[2*i+1]));
    AEP[i] = (2.0*(AALS[2*i]*AALS[2*i]+AALS[2*i+1]*AALS[2*i+1]));
  }
  
  int divs = 100;  // must be even - used to compute median
  
  if(divs/2+1 > imin) imin = divs/2+1;
  if(imax > NFFT/2-divs/2-1) imax =  NFFT/2-divs/2-1;
  
  //spline_fit(0, divs, imin, imax, XP, Xnoise, Xconf, TOBS, fstar, L);
  spline_fit(1, divs, imin, imax, AEP, AEnoise, AEconf, TOBS, fstar, L);

  
  medianX(imin, imax, fstar, L, XP, Xnoise, Xconf, TOBS);
  //medianAE(imin, imax, fstar, L, AEP, AEnoise, AEconf, TOBS);
  
  Outfile = fopen("Confusion_XAE_0.dat","w");
  for(i=imin; i<= imax; i++)
  {
    f = (double)(i)/TOBS;
    fprintf(Outfile,"%.12g %e %e %e %e\n", f, Xnoise[i], Xconf[i], AEnoise[i], AEconf[i]);
  }
  fclose(Outfile);
  
  return 0;
  
}
