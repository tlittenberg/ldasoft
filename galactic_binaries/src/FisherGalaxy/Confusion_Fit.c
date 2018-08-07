/*********************************************************/
/*                                                       */
/*       Confusion_Fit.c, Version 2.3, 4/28/2011         */
/*      Written by Neil Cornish & Tyson Littenberg       */
/*                                                       */
/* gcc -O2 -o Confusion_Fit Confusion_Fit.c arrays.c -lm */
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
  
  double f;
  char Gfile[50];
  double *XfLS, *AALS, *EELS;
  double *XP, *AEP;
  long i, imax, imin;
  long rseed;
  double SAE, SXYZ;
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

  printf("*   Observing Time: %.1f year (%f s)\n",TOBS/year,TOBS);
  printf("*\n");
  printf("***********************************************************************\n");

  XfLS = dvector(0,NFFT-1);  AALS = dvector(0,NFFT-1);  EELS = dvector(0,NFFT-1);
  
  
  imax = (long)ceil(4.0e-2*TOBS);
  imin = (long)floor(1.0e-4*TOBS);
  
  XfLS = dvector(0,NFFT-1);  AALS = dvector(0,NFFT-1); EELS = dvector(0,NFFT-1);
  
  
  //Data structure for interpolating orbits from file
  struct lisa_orbit *LISAorbit;
  LISAorbit = &orbit;
  
  //Set up orbit structure (allocate memory, read file, cubic spline)
  sprintf(Gfile,"%s",argv[2]);
  initialize_orbit(Gfile, LISAorbit);
  
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

  XP = dvector(imin,imax);  AEP = dvector(imin,imax);
  Xnoise = dvector(imin,imax);  Xconf = dvector(imin,imax);
  AEnoise = dvector(imin,imax);  AEconf = dvector(imin,imax);
  
  rseed = -7584529636;
  
  for(i=imin; i< imax; i++)
  {
    f = (double)(i)/TOBS;
    instrument_noise(f, fstar, L, &SAE, &SXYZ);
    XP[i] = (2.0*(XfLS[2*i]*XfLS[2*i] + XfLS[2*i+1]*XfLS[2*i+1]));
    //XP[i] = SXYZ*0.5*(pow(gasdev2(&rseed), 2.0)+pow(gasdev2(&rseed), 2.0));
    AEP[i] = (2.0*(AALS[2*i]*AALS[2*i]+AALS[2*i+1]*AALS[2*i+1]));
  }
  
  Outfile = fopen("Galaxy_XAE_Pow.dat","w");
  for(i=imin; i< imax; i++)
  {
    f = (double)(i)/TOBS;
    instrument_noise(f, fstar, L, &SAE, &SXYZ);
    fprintf(Outfile,"%.12g %e %e %e %e\n", f, XP[i], AEP[i], SXYZ, SAE);
  }
  fclose(Outfile);
  
  medianX(imin, imax, fstar, L, XP, Xnoise, Xconf, TOBS);
  medianAE(imin, imax, fstar, L, AEP, AEnoise, AEconf, TOBS);
  
  Outfile = fopen("Confusion_XAE_0.dat","w");
  for(i=imin; i<= imax; i++)
  {
    f = (double)(i)/TOBS;
    fprintf(Outfile,"%.12g %e %e %e %e\n", f, Xnoise[i], Xconf[i], AEnoise[i], AEconf[i]);
  }
  fclose(Outfile);
  
  Outfile = fopen("Noise_Pow.dat","w");
  for(i=imin; i< imax; i++)
  {
    f = (double)(i)/TOBS;
    instrument_noise(f, fstar, L, &SAE, &SXYZ);
    fprintf(Outfile,"%.12g %e %e %e %e %e %e\n", f, Xnoise[i], Xconf[i], SXYZ, AEnoise[i], AEconf[i], SAE);
  }
  fclose(Outfile);
  
  return 0;
  
}
