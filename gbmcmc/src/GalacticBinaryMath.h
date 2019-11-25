/*
*  Copyright (C) 2019 Tyson B. Littenberg (MSFC-ST12), Neil J. Cornish
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



#ifndef GalacticBinaryMath_h
#define GalacticBinaryMath_h

#include <stdio.h>
#include <stdlib.h>

double chirpmass(double m1, double m2);

double ipow(double x, int n);

double fourier_nwip(double *a, double *b, double *Sn, int n);
double snr(struct Source *source, struct Noise *noise);
double waveform_match(struct Source *a, struct Source *b, struct Noise *noise);

int binary_search(double *array, int nmin, int nmax, double x);

void matrix_eigenstuff(double **matrix, double **evector, double *evalue, int N);
void invert_matrix(double **matrix, int N);
void matrix_multiply(double **A, double **B, double **AB, int N);
void cholesky_decomp(double **A, double **L, int N);

/* ********************************************************************************** */
/*                                                                                    */
/*                                    Fourier Tools                                   */
/*                                                                                    */
/* ********************************************************************************** */

void fftw_wrapper(double *data, int N, int flag);

double power_spectrum(double *data, int n);


#endif /* GalacticBinaryMath_h */
