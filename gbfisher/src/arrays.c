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
