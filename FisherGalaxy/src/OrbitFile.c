/* gcc -O2 -o OrbitFile OrbitFile.c arrays.c -lm */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "arrays.h"
#include "Constants.h"
#include "Detector.h"
#include "Subroutines.h"
#include "../../galactic_binaries/src/LISA.h"

/*************************************************************************/
/*        Rigid approximation position of each LISA spacecraft           */
/*************************************************************************/
//void analytic_orbits(double t, double *x, double *y, double *z)
//{
//  double alpha;
//  double beta1, beta2, beta3;
//  double sa, sb, ca, cb;
//  
//  alpha = 2.*pi*fm*t + KAPPA;
//  
//  beta1 = 0. + LAMBDA;
//  beta2 = 2.*pi/3. + LAMBDA;
//  beta3 = 4.*pi/3. + LAMBDA;
//  
//  sa = sin(alpha);
//  ca = cos(alpha);
//  
//  
//  sb = sin(beta1);
//  cb = cos(beta1);
//  x[1] = AU*ca + AU*ECC*(sa*ca*sb - (1. + sa*sa)*cb);
//  y[1] = AU*sa + AU*ECC*(sa*ca*cb - (1. + ca*ca)*sb);
//  z[1] = -sq3*AU*ECC*(ca*cb + sa*sb);
//  
//  
//  sb = sin(beta2);
//  cb = cos(beta2);
//  x[2] = AU*ca + AU*ECC*(sa*ca*sb - (1. + sa*sa)*cb);
//  y[2] = AU*sa + AU*ECC*(sa*ca*cb - (1. + ca*ca)*sb);
//  z[2] = -sq3*AU*ECC*(ca*cb + sa*sb);
//  
//  sb = sin(beta3);
//  cb = cos(beta3);
//  x[3] = AU*ca + AU*ECC*(sa*ca*sb - (1. + sa*sa)*cb);
//  y[3] = AU*sa + AU*ECC*(sa*ca*cb - (1. + ca*ca)*sb);
//  z[3] = -sq3*AU*ECC*(ca*cb + sa*sb);
//  
//}
/*************************************************************************/


int main(int argc,char **argv)
{
  //-----------------------------------------------------------
  //
  // Eccentric Inclined Orbits a la MLDC
  //
  //-----------------------------------------------------------
  
  FILE *Outfile = fopen("EccentricInclined.txt","w");
  
  int i,j;

  struct Orbit *orbit = malloc(sizeof(struct Orbit));
  initialize_analytic_orbit(orbit);
  
  double t;
  double *x = malloc(sizeof(double)*4);//dvector(1,3);
  double *y = malloc(sizeof(double)*4);//dvector(1,3);
  double *z = malloc(sizeof(double)*4);//dvector(1,3);

  for(i=0; i<366*10; i++)
  {
    t = (double)i*24.*60.*60.;
    analytic_orbits(orbit,t,x,y,z);
    fprintf(Outfile,"%.12g ",t);
    for(j=1; j<=3; j++)
    {
      fprintf(Outfile,"%.12g %.12g %.12g ",x[j],y[j],z[j]);
    }
    fprintf(Outfile,"\n");
  }

  fclose(Outfile);
  
  return 0;
  
}


