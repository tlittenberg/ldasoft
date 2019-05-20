/*
 */

/***************************  REQUIRED LIBRARIES  ***************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

/* Square root of 3 */
#define SQ3   1.73205080757

/* Astronomical unit (meters) */
#define AU 1.49597870660e11

/* Speed of light (m/s) */
#define CLIGHT 299792458.

/* Mean arm length of constellation (m) */
#define Larm 2.5e9

static void analytic_orbits(double t, double *x, double *y, double *z)
{
  
  double L     = Larm;
  double fstar = CLIGHT/(2.0*M_PI*L);
  double ecc   = Larm/(2.0*SQ3*AU);
  double R     = AU*ecc;

  //double alpha = PI2*t/YEAR;
  double alpha = 2*M_PI*t*3.168753575e-8;
  
  /*
   double beta1 = 0.;
   double beta2 = 2.0943951023932; //2.*pi/3.;
   double beta3 = 4.18879020478639;//4.*pi/3.;
   */
  
  double sa = sin(alpha);
  double ca = cos(alpha);
  
  double sa2  = sa*sa;
  double ca2  = ca*ca;
  double saca = sa*ca;
  double AUca = AU*ca;
  double AUsa = AU*sa;
  
  double sb,cb;
  
  sb = 0.0;//sin(beta1);
  cb = 1.0;//cos(beta1);
  x[0] = AUca + R*(saca*sb - (1. + sa2)*cb);
  y[0] = AUsa + R*(saca*cb - (1. + ca2)*sb);
  z[0] = -SQ3*R*(ca*cb + sa*sb);
  
  sb = 0.866025403784439;//sin(beta2);
  cb = -0.5;//cos(beta2);
  x[1] = AUca + R*(saca*sb - (1. + sa2)*cb);
  y[1] = AUsa + R*(saca*cb - (1. + ca2)*sb);
  z[1] = -SQ3*R*(ca*cb + sa*sb);
  
  sb = -0.866025403784438;//sin(beta3);
  cb = -0.5;//cos(beta3);
  x[2] = AUca + R*(saca*sb - (1. + sa2)*cb);
  y[2] = AUsa + R*(saca*cb - (1. + ca2)*sb);
  z[2] = -SQ3*R*(ca*cb + sa*sb);
  
}

int main(int argc, char* argv[])
{
  double seconds_per_day = 86400.0*(365);
  double sampling_rate = 1.0/86400.0; //Hz
  double dt = 1./sampling_rate; //time step (s)
  
  //number of samples
  int N = (int)floor(seconds_per_day*sampling_rate);
  
  
  //cartesian location of spacecraft
  double *x   = malloc(3*sizeof(double));
  double *y   = malloc(3*sizeof(double));
  double *z   = malloc(3*sizeof(double));
  double *x_0 = malloc(3*sizeof(double));
  double *y_0 = malloc(3*sizeof(double));
  double *z_0 = malloc(3*sizeof(double));
  double *vx  = malloc(3*sizeof(double));
  double *vy  = malloc(3*sizeof(double));
  double *vz  = malloc(3*sizeof(double));

  FILE *orbit_file = fopen("orbit_file.dat","w");
  FILE *velocity_file = fopen("velocity_file.dat","w");

  
  double t = 0.0;
  analytic_orbits(t,x_0,y_0,z_0);
  for(int i=0; i<N; i++)
  {
    t += dt;
    fprintf(orbit_file,"%.12g ",t);
    fprintf(velocity_file,"%.12g ",t);
    analytic_orbits(t,x,y,z);
    for(int j=0; j<3; j++)
    {
      vx[j] = (x[j]-x_0[j])/dt;
      vy[j] = (y[j]-y_0[j])/dt;
      vz[j] = (z[j]-z_0[j])/dt;
    }

    
    for(int j=0; j<3; j++)
    {
      fprintf(orbit_file,"%.16g %.16g %.16g ",x[j],y[j],z[j]);
      fprintf(velocity_file,"%.16g %.16g %.16g ",vx[j],vy[j],vz[j]);
    }
    fprintf(orbit_file,"\n");
    fprintf(velocity_file,"\n");
    for(int j=0; j<3; j++)
    {
      x_0[j] = x[j];
      y_0[j] = y[j];
      z_0[j] = z[j];
    }

  }
  
}


