/*
 * gcc histograms.c -lm -o histograms
 *   -makes marginalized pdfs from MCMC chains
 *	-ignores first two columns as iteration and log-likelihood
 */

/***************************  REQUIRED LIBRARIES  ***************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <fftw3.h>
#include <gsl/gsl_sf_gamma.h>


/*************  PROTOTYPE DECLARATIONS FOR INTERNAL FUNCTIONS  **************/
double binomial(double x, double y);
double sinc(double n, double D);
double func(double x);

double blackman (int n, int N, double D);
double blackman3(int n, int N, double D);
double lagrange (int n, int N, double D);
double truncate (int n, int N, double D);

void fft(double *x, int N);

void print_usage(FILE *fptr);

/* ============================  MAIN PROGRAM  ============================ */

int main(int argc, char* argv[])
{
  
  
  int N = 16384*4;
  int M = 1000;
  double x;
  double *h = malloc(N*sizeof(double));
  double *w = malloc(N*sizeof(double));
  double *H = malloc(N*sizeof(double));
  
  double fs = 10;
  double D  = 0.5;
  double T = (double)N/fs;
  
  FILE *out;
  
  
  for(int i=0; i<N; i++)
  {
    h[i] = 0.0;
    w[i] = 0.0;
    H[i] = 0.0;
  }
  
  
  /* create filters */
  for(int i=0; i<M; i++)
  {
    x = (double)(i-M/2);
    int j = N/2+(i-M/2);
    h[j] = sinc(x,D);
    w[j] = blackman(x,M/2,D);
    H[j] = h[j]*w[j];
  }
  
  out = fopen("sinc.dat","w");
  for(int i=0; i<N; i++)fprintf(out,"%lg %lg %lg %lg\n",(double)(i-N/2),h[i],w[i],H[i]);
  fclose(out);

  
  /* FFT filters */
  out = fopen("dft.dat","w");
  fft(H,N);
  for(int i=0; i<N/2; i++)
  {
    x = (double)(i/T);
    fprintf(out,"%lg %lg %lg\n",x,H[2*i],H[2*i+1]);
  }
  fclose(out);
  
  
  return 0;
}

void fft(double *x, int N)
{
  double in[N];
  fftw_complex out[N/2+1];
  unsigned flags;
  fftw_plan plan = fftw_plan_dft_r2c_1d(N,in,out,flags);
  
  //load data int FFTW data types
  for(int i=0; i<N; i++)
  {
    in[i] = x[i];
  }
  
  fftw_execute(plan);
  
  //export data
  for(int i=0; i<N/2; i++)
  {
    x[2*i]   = out[i][0];
    x[2*i+1] = out[i][1];
  }
  fftw_destroy_plan(plan);
}

double func(double x)
{
  double f = 0.001;
  return 1.+sin(2.0*M_PI*f*x);
}

double sinc(double n, double D)
{
  if(fabs(n-D)<1e-8) return 1.0;
  else return sin(M_PI*(n-D))/(M_PI*(n-D));
}


double blackman3(int n, int N, double D)
{
  return blackman(n,N,D)*blackman(n,N,D)*blackman(n,N,D);
}

double blackman(int n, int N, double D)
{
  return 0.42 + 0.5*cos(M_PI*(double)(n-D)/(double)(N-1)) + 0.08*cos(2.0*M_PI*(double)(n-D)/(double)(N-1));
}

double truncate(int n, int N, double D)
{
  return 1.0;
}

double lagrange(int n, int N, double D)
{
  double td = 0.5*(double)(N)+D;
  
  return -M_PI * (double)N/sin(M_PI*td) * binomial(td, (double)N) * binomial((double)N,(double)(n+(double)(N)/2.));
}

double binomial(double x, double y)
{
  return gsl_sf_gamma(x+1.)/gsl_sf_gamma(y+1.)/gsl_sf_gamma(x-y+1.);
}

void print_usage(FILE *fptr)
{
  fprintf(fptr,"Usage: fdi D N window\n");
  fprintf(fptr,"Window options:\n");
  fprintf(fptr," truncate\n");
  fprintf(fptr," blackman\n");
}


