//
//  GalacticBinaryMath.c
//  
//
//  Created by Littenberg, Tyson B. (MSFC-ZP12) on 1/15/17.
//
//
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf.h>

#include "LISA.h"
#include "GalacticBinary.h"
#include "GalacticBinaryMath.h"


double chirpmass(double m1, double m2)
{
  return pow(m1*m2,3./5.)/pow(m1+m2,1./5.);
}


double fourier_nwip(double *a, double *b, double *Sn, int n)
{
  int i, j, k;
  double arg, product;
  double ReA, ReB, ImA, ImB;
  
  arg = 0.0;
  for(i=0; i<n; i++)
  {
    j = i * 2;
    k = j + 1;
    ReA = a[j]; ImA = a[k];
    ReB = b[j]; ImB = b[k];
    product = ReA*ReB + ImA*ImB;
    arg += product/Sn[i];
    //if(arg!=arg){printf("WTF?  %i,ReA=%g,ImA=%g,ReB=%g,ImB=%g,Sn=%g\n",i,ReA,ImA,ReB,ImB,Sn[i]);}
  }
  
  return(4.0*arg);
}

double snr(struct Source *source, struct Noise *noise)
{
  double snr2=0.0;
  switch(source->tdi->Nchannel)
  {
    case 1: //Michelson
      snr2 += fourier_nwip(source->tdi->X,source->tdi->X,noise->SnX,source->tdi->N);
      break;
    case 2: //A&E
      snr2 += fourier_nwip(source->tdi->A,source->tdi->A,noise->SnA,source->tdi->N);
      snr2 += fourier_nwip(source->tdi->E,source->tdi->E,noise->SnE,source->tdi->N);
      break;
  }
  
  return(sqrt(snr2));
}

void matrix_eigenstuff(double **matrix, double **evector, double *evalue, int N)
{
  int i,j;
  
  // Don't let errors kill the program (yikes)
  gsl_set_error_handler_off ();
  int err=0;
  
  // Find eigenvectors and eigenvalues
  gsl_matrix *GSLfisher = gsl_matrix_alloc(N,N);
  gsl_matrix *GSLcovari = gsl_matrix_alloc(N,N);
  gsl_matrix *GSLevectr = gsl_matrix_alloc(N,N);
  gsl_vector *GSLevalue = gsl_vector_alloc(N);
  
  for(i=0; i<N; i++)
  {
    for(j=0; j<N; j++)
    {
      if(matrix[i][j]!=matrix[i][j])fprintf(stderr,"GalacticBinaryMath.c:83: WARNING: nan matrix element, now what?\n");
      gsl_matrix_set(GSLfisher,i,j,matrix[i][j]);
    }
  }
  
  // sort and put them into evec
  gsl_eigen_symmv_workspace * workspace = gsl_eigen_symmv_alloc (N);
  gsl_permutation * permutation = gsl_permutation_alloc(N);
  err += gsl_eigen_symmv (GSLfisher, GSLevalue, GSLevectr, workspace);
  err += gsl_eigen_symmv_sort (GSLevalue, GSLevectr, GSL_EIGEN_SORT_ABS_ASC);

  // eigenvalues destroy matrix
  for(i=0; i<N; i++) for(j=0; j<N; j++) gsl_matrix_set(GSLfisher,i,j,matrix[i][j]);

  err += gsl_linalg_LU_decomp(GSLfisher, permutation, &i);
  err += gsl_linalg_LU_invert(GSLfisher, permutation, GSLcovari);
  
  if(err>0)
  {
    fprintf(stderr,"GalacticBinaryMath.c:98: WARNING: singluar matrix, treating matrix as diagonal\n");
    fflush(stderr);
    for(i=0; i<N; i++)for(j=0; j<N; j++)
    {
      evector[i][j] = 0.0;
      if(i==j)
      {
        evector[i][j]=1.0;
        evalue[i]=1./matrix[i][j];
      }
    }
    
  }
  else
  {
    
    //unpack arrays from gsl inversion
    for(i=0; i<N; i++)
    {
      evalue[i] = gsl_vector_get(GSLevalue,i);
      for(j=0; j<N; j++)
      {
        evector[i][j] = gsl_matrix_get(GSLevectr,i,j);
        if(evector[i][j] != evector[i][j]) evector[i][j] = 0.;
      }
    }
    
    //for(i=0;i<N-1;i++)for(j=i+1;j<N;j++) gsl_matrix_set(GSLcovari,j,i, gsl_matrix_get(GSLcovari,i,j) );

    //copy covariance matrix back into Fisher
    for(i=0; i<N; i++)
    {
      for(j=0; j<N; j++)
      {
        matrix[i][j] = gsl_matrix_get(GSLcovari,i,j);
      }
    }

    //cap minimum size eigenvalues
    for(i=0; i<N; i++)
    {
      if(evalue[i] != evalue[i] || evalue[i] <= 10.0) evalue[i] = 10.;
    }
  }
  
  gsl_vector_free (GSLevalue);
  gsl_matrix_free (GSLfisher);
  gsl_matrix_free (GSLcovari);
  gsl_matrix_free (GSLevectr);
  gsl_eigen_symmv_free (workspace);
  gsl_permutation_free (permutation);
}

void invert_matrix(double **matrix, int N)
{
  int i,j;
  
  // Don't let errors kill the program (yikes)
  gsl_set_error_handler_off ();
  int err=0;
  
  // Find eigenvectors and eigenvalues
  gsl_matrix *GSLmatrix = gsl_matrix_alloc(N,N);
  gsl_matrix *GSLinvrse = gsl_matrix_alloc(N,N);
  
  for(i=0; i<N; i++)
  {
    for(j=0; j<N; j++)
    {
      if(matrix[i][j]!=matrix[i][j])fprintf(stderr,"GalacticBinaryMath.c:172: WARNING: nan matrix element, now what?\n");
      gsl_matrix_set(GSLmatrix,i,j,matrix[i][j]);
    }
  }
  
  gsl_permutation * permutation = gsl_permutation_alloc(N);
  
  err += gsl_linalg_LU_decomp(GSLmatrix, permutation, &i);
  err += gsl_linalg_LU_invert(GSLmatrix, permutation, GSLinvrse);
  
  if(err>0)
  {
    fprintf(stderr,"GalacticBinaryMath.c:184: WARNING: singluar matrix\n");
    fflush(stderr);
  }
  else
  {
    //copy covariance matrix back into Fisher
    for(i=0; i<N; i++)
    {
      for(j=0; j<N; j++)
      {
        matrix[i][j] = gsl_matrix_get(GSLinvrse,i,j);
      }
    }
  }
  
  gsl_matrix_free (GSLmatrix);
  gsl_matrix_free (GSLinvrse);
  gsl_permutation_free (permutation);
}


/* ********************************************************************************** */
/*																					                                          */
/*                                   Fourier Tools                                    */
/*																					                                          */
/* ********************************************************************************** */

void dfour1(double data[], unsigned long nn, int isign)
{
  unsigned long n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  double tempr,tempi;
  
  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) {
    if (j > i) {
      SWAP(data[j],data[i]);
      SWAP(data[j+1],data[i+1]);
    }
    m=n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (n > mmax) {
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax);
    wtemp=gsl_sf_sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=gsl_sf_sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
        j=i+mmax;
        tempr=wr*data[j]-wi*data[j+1];
        tempi=wr*data[j+1]+wi*data[j];
        data[j]=data[i]-tempr;
        data[j+1]=data[i+1]-tempi;
        data[i] += tempr;
        data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}

void drealft(double data[], unsigned long n, int isign)
{
  void dfour1(double data[], unsigned long nn, int isign);
  unsigned long i,i1,i2,i3,i4,np3;
  double c1=0.5,c2,h1r,h1i,h2r,h2i;
  double wr,wi,wpr,wpi,wtemp,theta;
  
  theta=3.141592653589793/(double) (n>>1);
  if (isign == 1) {
    c2 = -0.5;
    dfour1(data,n>>1,1);
  } else {
    c2=0.5;
    theta = -theta;
  }
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  wr=1.0+wpr;
  wi=wpi;
  np3=n+3;
  for (i=2;i<=(n>>2);i++) {
    i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
    h1r=c1*(data[i1]+data[i3]);
    h1i=c1*(data[i2]-data[i4]);
    h2r = -c2*(data[i2]+data[i4]);
    h2i=c2*(data[i1]-data[i3]);
    data[i1]=h1r+wr*h2r-wi*h2i;
    data[i2]=h1i+wr*h2i+wi*h2r;
    data[i3]=h1r-wr*h2r+wi*h2i;
    data[i4] = -h1i+wr*h2i+wi*h2r;
    wr=(wtemp=wr)*wpr-wi*wpi+wr;
    wi=wi*wpr+wtemp*wpi+wi;
  }
  if (isign == 1) {
    data[1] = (h1r=data[1])+data[2];
    data[2] = h1r-data[2];
  } else {
    data[1]=c1*((h1r=data[1])+data[2]);
    data[2]=c1*(h1r-data[2]);
    dfour1(data,n>>1,-1);
  }
}

double power_spectrum(double *data, int n)
{
  int i,j;
  double Re, Im;
  
  i = 2*n;
  j = i+1;
  Re = data[i];
  Im = data[j];
  
  return (Re*Re + Im*Im);
}
