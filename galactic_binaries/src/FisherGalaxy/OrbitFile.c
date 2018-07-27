/* gcc -O2 -o OrbitFile OrbitFile.c arrays.c -lm */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "arrays.h"
#include "Constants.h"
#include "Detector.h"

void spacecraft(double t,  double *x, double *y, double *z);

int main(int argc,char **argv)
{
  //-----------------------------------------------------------
  //
  // Eccentric Inclined Orbits a la MLDC
  //
  //-----------------------------------------------------------

  FILE *Outfile = fopen("EccentricInclined.txt","w");

  int n,i,j;
  double t,*x,*y,*z;
  x=dvector(1,3);
  y=dvector(1,3);
  z=dvector(1,3);
  for(i=0; i<366*5; i++)
  {
    t = (double)i*24.*60.*60.;
    spacecraft(t,x,y,z);
    fprintf(Outfile,"%.12g ",t);
    for(j=1; j<=3; j++)
    {
      fprintf(Outfile,"%.12g %.12g %.12g ",x[j],y[j],z[j]);
    }
    fprintf(Outfile,"\n");
  }
  free_dvector(x,1,3);
  free_dvector(y,1,3);
  free_dvector(z,1,3);
  fclose(Outfile);

  return 0;

}


/*************************************************************************/
/*        Rigid approximation position of each LISA spacecraft           */
/*************************************************************************/
void spacecraft(double t, double *x, double *y, double *z)
{
  double alpha;
  double beta1, beta2, beta3;
  double sa, sb, ca, cb;

  alpha = 2.*pi*fm*t + kappa;

  beta1 = 0. + lambda;
  beta2 = 2.*pi/3. + lambda;
  beta3 = 4.*pi/3. + lambda;

  sa = sin(alpha);
  ca = cos(alpha);


  sb = sin(beta1);
  cb = cos(beta1);
  x[1] = AU*ca + AU*ec*(sa*ca*sb - (1. + sa*sa)*cb);
  y[1] = AU*sa + AU*ec*(sa*ca*cb - (1. + ca*ca)*sb);
  z[1] = -sq3*AU*ec*(ca*cb + sa*sb);


  sb = sin(beta2);
  cb = cos(beta2);
  x[2] = AU*ca + AU*ec*(sa*ca*sb - (1. + sa*sa)*cb);
  y[2] = AU*sa + AU*ec*(sa*ca*cb - (1. + ca*ca)*sb);
  z[2] = -sq3*AU*ec*(ca*cb + sa*sb);

  sb = sin(beta3);
  cb = cos(beta3);
  x[3] = AU*ca + AU*ec*(sa*ca*sb - (1. + sa*sa)*cb);
  y[3] = AU*sa + AU*ec*(sa*ca*cb - (1. + ca*ca)*sb);
  z[3] = -sq3*AU*ec*(ca*cb + sa*sb);
  
}
/*************************************************************************/

