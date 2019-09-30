//
//  GalacticBinaryMath.h
//
//
//  Created by Littenberg, Tyson B. (MSFC-ZP12) on 1/15/17.
//
//

#ifndef GalacticBinaryMath_h
#define GalacticBinaryMath_h

#include <stdio.h>
#include <stdlib.h>

double swap, tempr;
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

double chirpmass(double m1, double m2);

double ipow(double x, int n);

double fourier_nwip(double *a, double *b, double *Sn, int n);
double snr(struct Source *source, struct Noise *noise);

int binary_search(double *array, int nmin, int nmax, double x);

void matrix_eigenstuff(double **matrix, double **evector, double *evalue, int N);
void invert_matrix(double **matrix, int N);
void matrix_multiply(double **A, double **B, double **AB, int N);
void cholesky_decomp(double matrix[8][8], double **tensor_out, int N);

/* ********************************************************************************** */
/*                                                                                    */
/*                                    Fourier Tools                                   */
/*                                                                                    */
/* ********************************************************************************** */

void dfour1(double data[], unsigned long nn, int isign);
void drealft(double data[], unsigned long n, int isign);

double power_spectrum(double *data, int n);


#endif /* GalacticBinaryMath_h */
