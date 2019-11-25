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

int *int_vector(int N);
void free_int_vector(int *v);

double *double_vector(int N);
void free_double_vector(double *v);

double **double_matrix(int N, int M);
void free_double_matrix(double **m, int N);

double ***double_tensor(int N, int M, int L);
void free_double_tensor(double ***t, int N, int M);
