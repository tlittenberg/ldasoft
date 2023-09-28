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

/**
@file GalacticBinaryMath.h
\brief Frequently encountered math functions for `gbmcmc`.
*/


#ifndef GalacticBinaryMath_h
#define GalacticBinaryMath_h

#include <stdio.h>
#include <stdlib.h>

/**
\brief Compute chirp mass from component masses
  
 @param m1 mass of primary
 @param m2 mass of secondary
 @return \f$ \mathcal{M} \equiv (m_1 m_2)^{3/5} (m_1+m_2)^{-1/5} \f$
 */
double chirpmass(double m1, double m2);

/**
\brief Compute integer powers of by brute force multiplication
 
 Faster than pow() for small integers
  
 @param x variable
 @param n exponent
 @return \f$x^n =  x\times x \times\ ...\times x\f$
 */
double ipow(double x, int n);

/**
\brief Compute power of complex amplitude in single element of data
 
 Assumes the `data` array has alternating real and imaginary terms.
 
 @param data complex amplitude array
 @param n desired sample
 @return \f$P =  d^*_n d_n \f$
 */
double power_spectrum(double *data, int n);

/**
\brief Fourier-domain noise weighted inner product
  
 @param a complex amplitude array
 @param b complex amplitude array
 @param invC inverse covariance matrix
 @param n number of frequency bins in sum
 @return \f$(a|b) =  4 \sum_n a^*_n b_n C^-1_n \f$
 */
double fourier_nwip(double *a, double *b, double *invC, int n);

/**
\brief Signal to noise ratio
  
 Calls nwip() to compute inner product of waveform with itself summed over all frequencies and data channels `I`. If 4-link data, use single TDI channel `I` = `X`. If 6-link, use `I` = {`A`, `E`}
 
 @param source structure containing source parameters and waveform \f$h\f$
 @param noise structure containing noise model \f$S_n\f$
 @return \f$\rho =  \sqrt{\sum_I (h_I|h_I)} \f$
 */
double snr(struct Source *source, struct Noise *noise);

/**
\brief Analytic approximation to SNR
  
 Not exactly what is in the paper. Expression has been calibrated against snr().
 
 @param A gravitational wave amplitude \f$\mathcal{A}\f$
 @param Sn noise power spectral density
 @param Sf correction based on sources proximity to transfer frequency
 @param sqT \f$\sqrt{T}\f$
 @return \f$ \rho \approx  \frac{1}{2} \mathcal{A} \sqrt{T} S_f / S_n \f$
 */
double analytic_snr(double A, double Sn, double Sf, double sqT);

/**
\brief Compute prior on SNR
  
 Shouldn't this be in GalacticBinaryPrior.h?
 Peak of distribution `SNRPEAK` \f$\rho_*\f$ defined in Constants.h
 
 @param SNR signal to noise ratio \f$\rho\f$ of source
 @return \f$ p(\rho) = \frac{3\rho}{4 \rho_*^2 \left(1 + \frac{\rho}{4\rho_*}
 \right)^{5}} \f$
 */
double snr_prior(double SNR);

/**
\brief Compute match between waveforms
   
 @param a waveform \f$h_a\f$
 @param b waveform \f$h_b\f$
 @return \f$  M = \frac{(h_a|h_b)}{\sqrt{(h_a|h_a)+(h_b|h_b)}} \f$
 */
double waveform_match(struct Source *a, struct Source *b, struct Noise *noise);

/**
\brief Compute distance between waveforms
   
 @param a waveform \f$h_a\f$
 @param b waveform \f$h_b\f$
 @return \f$  D = (h_a-h_b | h_a-h_b) \f$
 */
double waveform_distance(struct Source *a, struct Source *b, struct Noise *noise);

/**
\brief Our implementation of the recursive binary search algorithm
   
 @param array list of values to search over
 @param nmin starting index of search
 @param nmax stopping index of search
 @param x value to search for in array
 @return index in array of x
 */
int binary_search(double *array, int nmin, int nmax, double x);

/**
\brief Computes eigenvectors and eigenvalues of matrix
   
 @param[in] N size of matrix
 @param[in] matrix the input \f$N\times N\f$ matrix
 @param[out] evector 2D array of eigenvector components
 @param[out] evalue 1D array of eigenvalues, corresponding with columns of evector
 */
void matrix_eigenstuff(double **matrix, double **evector, double *evalue, int N);
/**
\brief Matrix inversion with replacement
   
 @param[in] N size of matrix
 @param[in,out] matrix the input \f$N\times N\f$ matrix, replaced with inverse
 */
void invert_matrix(double **matrix, int N);

/**
\brief Matrix multiplication
   
 @param[in] N size of matrix
 @param[in] A \f$N\times N\f$ matrix \f$A\f$
 @param[in] B \f$N\times N\f$ matrix \f$B\f$
 @param[out] AB the matrix product \f$AB\f$
 */
void matrix_multiply(double **A, double **B, double **AB, int N);

/**
\brief Wrapper to GSL Cholesky decomposition routine
   
 Copies matrices into `gsl_matrix` structures and then fills the lower half of \f$L\f$ with the result.  GSL returns the original matrix in the upper half of \f$L\f$. We replace the upper half of \f$L\f$ with all 0s. Contents of \f$A\f$ are unaltered.
 
 @param[in] N size of matrix
 @param[in] A \f$N\times N\f$ matrix \f$A\f$
 @param[out] L Cholesky decomposition of \f$A = LL^{-1}\f$
 */
void cholesky_decomp(double **A, double **L, int N);



#endif /* GalacticBinaryMath_h */
