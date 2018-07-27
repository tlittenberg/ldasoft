/*********************************************************/
/*                                                       */
/*          Density.c, Version 2.3, 4/28/2011            */                                                            
/*      Written by Neil Cornish & Tyson Littenberg       */                                          
/*                                                       */
/*       gcc -o Density Density.c arrays.c -lm           */
/*                                                       */
/*********************************************************/


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
    fscanf(Input, "%lf%lf%lf%lf%lf%lf%lf%lf\n", &f, &fdot, &theta, &phi, &Amp, &iota, &psi, &phase);

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
