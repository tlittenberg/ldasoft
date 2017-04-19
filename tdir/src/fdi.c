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

#include <gsl/gsl_sf_gamma.h>


/*************  PROTOTYPE DECLARATIONS FOR INTERNAL FUNCTIONS  **************/
double binomial(double x, double y);
double sinc(double x, double fs);
double func(double x);

double blackman (int n, int N, double D);
double blackman3(int n, int N, double D);
double lagrange (int n, int N, double D);
double truncate (int n, int N, double D);

void convolve(double *d, double *h, double *hg, int NMAX, int N);

void print_usage(FILE *fptr);

/* ============================  MAIN PROGRAM  ============================ */

int main(int argc, char* argv[])
{
  if(argc!=4)
  {
    print_usage(stderr);
    exit(1);
  }
  int N=atoi(argv[2]);
  int NMAX=8192;
  double *t  = malloc(sizeof(double)*NMAX);
  double *d  = malloc(sizeof(double)*NMAX);
  double *h  = malloc(sizeof(double)*N);
  double *hg = malloc(sizeof(double)*NMAX);
  double fs = 1.0;
  double D = atof(argv[1]);
  double T = N;
  
  double (*w)(int,int,double); //window
  
  FILE *fptr;
  
  fptr = fopen("data.dat","w");
  
  for(int i=0; i<NMAX; i++)
  {
    t[i] = (double)(i);
    d[i] = func(t[i]);
    fprintf(fptr,"%lg %lg %lg\n",t[i],h[i],d[i]);
  }
  fclose(fptr);
  
  
  int windowFlag=-1;
  
  if(!strcmp(argv[3],"truncate"))       windowFlag = 0;
  else if(!strcmp(argv[3],"blackman"))  windowFlag = 1;
  else if(!strcmp(argv[3],"blackman3")) windowFlag = 2;
  else if(!strcmp(argv[3],"lagrange"))  windowFlag = 3;
  
  printf("window flag for function %s = %i\n",argv[3],windowFlag);
  switch(windowFlag)
  {
    case 0:
      w = &truncate;
      break;
    case 1:
      w = &blackman;
      break;
    case 2:
      w = &blackman3;
      break;
    case 3:
      w = &lagrange;
      break;
    default:
      fprintf(stderr,"Unknown window function %s\n", argv[3]);
      print_usage(stderr);
      exit(1);
      break;
  }
  
  char filename[128];
  sprintf(filename,"filter_%g_%i_%s.dat",D,N,argv[3]);
  fptr = fopen(filename,"w");
  
  for(int i=0; i<N; i++)
  {
    h[i] = (*w)(i-N/2,N,D)*sinc((double)(i-N/2)-D,fs);
    fprintf(fptr,"%i %g %g %g\n",i-N/2,(*w)(i-N/2,N,D),sinc((double)(i-N/2)-D,fs),h[i]);
  }
  
  fclose(fptr);

  convolve(d, h, hg, NMAX, N);
  
  sprintf(filename,"result_%g_%i_%s.dat",D,N,argv[3]);
  fptr = fopen(filename,"w");
  for(int i=0; i<NMAX; i++)
    fprintf(fptr,"%lg %lg %lg\n",t[i],hg[i],func(t[i]+D));
  fclose(fptr);

  
  
  
  return 0;
}

double func(double x)
{
  double f = 0.001;
  return 1.+sin(2.0*M_PI*f*x);
}

double sinc(double t, double fs)
{
  if(fabs(t)<1e-8) return 1.0;
  else return sin(M_PI*fs*t)/(M_PI*fs*t);
}

void convolve(double *d, double *h, double *hg, int NMAX, int N)
{
  for(int n=0; n<NMAX; n++)
  {
    hg[n]=0.0;
    for(int m=0; m<N; m++)
    {
      int k = n+(m-N/2);
//      if(k>=N)k=k-N;
//      if(k<0) k=k+N;
      if(k>=0 && k<NMAX)hg[n] += h[m]*d[k];
    }
  }
}

double blackman3(int n, int N, double D)
{
  return blackman(n,N,D)*blackman(n,N,D)*blackman(n,N,D);
}

double blackman(int n, int N, double D)
{
  return 0.42 + 0.5*cos(M_PI*(double)n/(double)(N-1)) + 0.08*cos(2.0*M_PI*(double)n/(double)(N-1));
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


