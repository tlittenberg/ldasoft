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
  
  
  int N = 16384;
  int M = atoi(argv[1]);
  double x;
  double *h = malloc(N*sizeof(double));
  double *w = malloc(N*sizeof(double));
  double *H = malloc(N*sizeof(double));
  double *H_blackman  = malloc(N*sizeof(double));
  double *H_blackman3 = malloc(N*sizeof(double));
  double *H_lagrange  = malloc(N*sizeof(double));
  double *H_truncate  = malloc(N*sizeof(double));
  
  double fs = 10;
  double D  = 0.5;
  double T = (double)N/fs;
  
  FILE *out;
  
  
  for(int i=0; i<N; i++)
  {
    h[i] = 0.0;
    w[i] = 0.0;
    H[i] = 0.0;
    H_blackman[i]  = 0.0;
    H_blackman3[i] = 0.0;
    //H_lagrange[i]  = 0.0;
    H_truncate[i]  = 0.0;
  }
  
  
  /* create filters */
  for(int i=0; i<M; i++)
  {
    x = (double)(i-M/2);
    int j = N/2+(i-M/2);
    h[j] = sinc(x,D);
    w[j] = blackman(x,M/2,D);
    H_blackman[j]  = h[j]*blackman(x,M/2,D);
    H_blackman3[j] = h[j]*blackman3(x,M/2,D);
    //H_lagrange[j]  = h[j]*lagrange(x,M/2,D);
    H_truncate[j]  = h[j]*truncate(x,M/2,D);
  }
  
  out = fopen("sinc.dat","w");
  for(int i=0; i<N; i++)fprintf(out,"%lg %lg %lg %lg %lg\n",(double)(i-N/2),h[i],H_blackman[i],H_blackman3[i],H_truncate[i]);
  fclose(out);

  
  /* FFT filters */
  char filename[128];
  sprintf(filename,"dft_%i.dat",M);
  out = fopen("dft.dat","w");
  fft(H,N);
  fft(H_blackman,N);
  fft(H_blackman3,N);
  //fft(H_lagrange,N);
  fft(H_truncate,N);
  
  
  double max_blackman = -1e60;
  double max_blackman3= -1e60;
  //double max_lagrange= -1e60;
  double max_truncate= -1e60;
  
  double Re_diff;
  double Im_diff;
  double mag;
  
  
  for(int i=0; i<N/2; i++)
  {
    x = (double)(i/T);
    fprintf(out,"%lg %lg %lg\n",x,H[2*i],H[2*i+1]);
    
    
    Re_diff = fabs(H_blackman[2*i])   - cos(2.0*M_PI*x*0.5/fs);
    Im_diff = fabs(H_blackman[2*i+1]) - sin(2.0*M_PI*x*0.5/fs);
    mag = sqrt(Re_diff*Re_diff + Im_diff*Im_diff);
    if(x < 1.0 && mag > max_blackman) max_blackman=mag;

    Re_diff = fabs(H_blackman3[2*i])   - cos(2.0*M_PI*x*0.5/fs);
    Im_diff = fabs(H_blackman3[2*i+1]) - sin(2.0*M_PI*x*0.5/fs);
    mag = sqrt(Re_diff*Re_diff + Im_diff*Im_diff);
    if(x < 1.0 && mag > max_blackman3) max_blackman3=mag;

    
//    Re_diff = fabs(H_lagrange[2*i])   - cos(2.0*M_PI*x*0.5/fs);
//    Im_diff = fabs(H_lagrange[2*i+1]) - sin(2.0*M_PI*x*0.5/fs);
//    mag = sqrt(Re_diff*Re_diff + Im_diff*Im_diff);
//    if(x < 1.0 && mag > max_lagrange) max_lagrange=mag;
//    
    
    Re_diff = fabs(H_truncate[2*i])   - cos(2.0*M_PI*x*0.5/fs);
    Im_diff = fabs(H_truncate[2*i+1]) - sin(2.0*M_PI*x*0.5/fs);
    mag = sqrt(Re_diff*Re_diff + Im_diff*Im_diff);
    if(x < 1.0 && mag > max_truncate) max_truncate=mag;

    
  }
  printf("%i ",M);
  printf("%g ",max_truncate);
  printf("%g ",max_blackman);
  printf("%g ",max_blackman3);
  //printf("%g ",max_lagrange);
  printf("\n");
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
  
  printf("Lagrange: %i %g\n",n,-M_PI * (double)N/sin(M_PI*td) * binomial(td, (double)N) * binomial((double)N,(double)(n+(double)(N)/2.)));
  
  return -M_PI * (double)N/sin(M_PI*td) * binomial(td, (double)N) * binomial((double)N,(double)(n+(double)(N)/2.));
}

double binomial(double x, double y)
{
  //int gsl_sf_lngamma_sgn_e (double x, gsl_sf_result * result_lg, double * sgn)
  gsl_sf_result result_lg;
  double sgn;
  
  if(x-y+1 <= 0) return 1;
  
  printf("x=%g,y=%g, gamma=%g\n",x,y,gsl_sf_gamma(x+1.)/gsl_sf_gamma(y+1.)/gsl_sf_gamma(x-y+1.));
  
//  int t1 = gsl_sf_lngamma_sgn_e(x+1.,&result_lg,&sgn);
//  double g1 = sgn * exp(result_lg.val);
//
//  int t2 = gsl_sf_lngamma_sgn_e(y+1.,&result_lg,&sgn);
//  double g2 = sgn * exp(result_lg.val);
//
//  int t3 = gsl_sf_lngamma_sgn_e(x-y+1.,&result_lg,&sgn);
//  double g3 = sgn * exp(result_lg.val);

//  printf("%g %g %g %g %g\n",x,y,g1,g2,g3);
  return gsl_sf_gamma(x+1.)/gsl_sf_gamma(y+1.)/gsl_sf_gamma(x-y+1.);
//return g1/g2/g3;
}

void print_usage(FILE *fptr)
{
  fprintf(fptr,"Usage: fdi D N window\n");
  fprintf(fptr,"Window options:\n");
  fprintf(fptr," truncate\n");
  fprintf(fptr," blackman\n");
}


