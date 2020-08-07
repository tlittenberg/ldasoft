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
#include "Constants.h"
#include "Detector.h"
#include "arrays.h"
#include <stdio.h>
#include <stddef.h>

void KILL(char*);

int main(int argc,char **argv)
{
  
  FILE* Output;
  FILE* Input;
  double f, fdot, Amp, theta, phi, iota, psi, phase, fddot;
  double fmin, fmax;
  double bw1, bw2, max, maxf;
  double *den;
  int bins, i;
  
  if(argc !=2) KILL("./Denisty detections.dat");
  
  Input = fopen(argv[1],"r");
  
  fmin = 1.0e-4;
  fmax = 5.0e-3;
  
  bins = 100;
  
  den = dvector(0,bins);
  
  for(i=0; i <= bins; i++)
  {
    den[i] = 0.0;
  }
  
  bw1 = (fmax-fmin)/(double)(bins);
  bw2 = 1.0/(bw1*T);
  
  while ( !feof(Input) )
  {
    int check = fscanf(Input, "%lf%lf%lf%lf%lf%lf%lf%lf\n", &f, &fdot, &theta, &phi, &Amp, &iota, &psi, &phase);
    if(!check)
    {
        fprintf(stderr,"Error reading %s\n",argv[1]);
        exit(1);
    }
    i = (int)((f-fmin)/bw1);
    
    if(i <= bins) den[i] += bw2;
    
  }
  
  Output = fopen("Density.dat","w");
  
  max = 0.0;
  maxf = fmin;
  for(i=0; i <= bins; i++)
  {
    fprintf(Output,"%e %e\n", fmin+(double)(i)*bw1, den[i]);
    if(den[i] > max)
	   {
       max = den[i];
       maxf = fmin+(double)(i)*bw1;
     }
  }
  
  //printf("f:\n %e\n  max den:\n %e\n", maxf, max);
  printf("%g\n",max);
  fclose(Output);
  
  
  return 0;
}

void KILL(char* Message)
{
  printf("\a\n");
  printf("%s",Message);
  printf("Terminating the program.\n\n\n");
  exit(1);
  
  
  return;
}
