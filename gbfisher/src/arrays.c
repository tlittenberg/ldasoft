/* CAUTION: This is the ANSI C (only) version of the Numerical Recipes
 utility file nrutil.c.  Do not confuse this file with the same-named
 file nrutil.c that is supplied in the 'misc' subdirectory.
 *That* file is the one from the book, and contains both ANSI and
 traditional K&R versions, along with #ifdef macros to select the
 correct version.  *This* file contains only ANSI C.               */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>


int *int_vector(int N)
{
   return malloc( (N+1) * sizeof(int) );
}

void free_int_vector(int *v)
{
   free(v);
}

double *double_vector(int N)
{
  return malloc( (N+1) * sizeof(double) );
}

void free_double_vector(double *v)
{
   free(v);
}

double **double_matrix(int N, int M)
{
  int i;
  double **m = malloc( (N+1) * sizeof(double *));
  
  for(i=0; i<N+1; i++)
  {
    m[i] = malloc( (M+1) * sizeof(double));
  }

  return m;
}

void free_double_matrix(double **m, int N)
{
  int i;
  for(i=0; i<N+1; i++) free_double_vector(m[i]);
  free(m);
}

double ***double_tensor(int N, int M, int L)
{
  int i,j;
  
  double ***t = malloc( (N+1) * sizeof(double **));
  for(i=0; i<N+1; i++)
  {
    t[i] = malloc( (M+1) * sizeof(double *));
    for(j=0; j<M+1; j++)
    {
      t[i][j] = malloc( (L+1) * sizeof(double));
    }
  }
  
  return t;
}

void free_double_tensor(double ***t, int N, int M)
{
  int i;
  
  for(i=0; i<N+1; i++) free_double_matrix(t[i],M);
  
  free(t);
}
