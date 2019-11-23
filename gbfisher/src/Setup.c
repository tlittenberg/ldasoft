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

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stddef.h>
#include "Constants.h"
#include "Subroutines.h"

int main(int argc,char **argv)
{
  
  FILE* Output;
  double kappa, lambda, L, Sacc, Sps, dt, T, ec, fstar;
  long NFFT;
  
  Output = fopen("Detector.h","w");
  
  if(argc !=8) KILL("./Setup kappa, lambda, Larm (m), Sacc (m^2 s^-4 Hz^-1), Spos (m^2 Hz^-1), Cadence (s), Nsamples\n");
  
  kappa = atof(argv[1]);
  lambda = atof(argv[2]);
  L = atof(argv[3]);
  Sacc = atof(argv[4]);
  Sps = atof(argv[5]);
  dt = atof(argv[6]);
  NFFT = atoi(argv[7]);
  
  T = dt*(double)(NFFT);
  fstar = clight/(2.0*pi*L);
  ec = (L/(2.0*sqrt(3.0)*Rgc));
  
  
  fprintf(Output, "         /* ----------------  DETECTOR CONSTANTS  ---------------- */               \n\n");
  
  fprintf(Output, " /* Observation time (seconds) */ \n");
  fprintf(Output, "#define T %.12f\n\n", T);
  
  fprintf(Output, " /* Number of data points */ \n");
  fprintf(Output, "#define NFFT %ld\n\n", NFFT);
  
  fprintf(Output, "#define dt %f\n\n", dt);
  
  fprintf(Output, " /* Initial azimuthal position of the guiding center */ \n");
  fprintf(Output, "#define kappa %f\n\n", kappa);
  
  fprintf(Output, " /* Initial orientation of the LISA constellation */ \n");
  fprintf(Output, "#define lambda %f\n\n", lambda);
  
  fprintf(Output, " /* Mean arm length of the LISA detector (meters) */ \n");
  fprintf(Output, "#define LARM %e\n\n", L);
  
  fprintf(Output, " /* Photon shot noise power */ \n");
  fprintf(Output, "#define Sps %e\n\n", Sps);
  
  fprintf(Output, " /* Acceleration noise power */ \n");
  fprintf(Output, "#define Sacc %e\n\n", Sacc);
  
  fprintf(Output, " /* Transfer frequency */ \n");
  fprintf(Output, "#define FSTAR %.10f\n\n", fstar);
  
  fprintf(Output, " /* LISA orbital eccentricity */ \n");
  fprintf(Output, "#define ec %.10f\n\n", ec);
  
  fprintf(Output, " /* details about noise model */ \n");
  fprintf(Output, "#define noiseFlag 2\n");
  fprintf(Output, "#define redden 1\n");
  fprintf(Output, "#define fred 0.0001\n\n");
  
  fclose(Output);
  
  return 0;
}

