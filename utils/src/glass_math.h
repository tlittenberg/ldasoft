//
//  glass_math.h
//  utils
//
//  Created by Tyson Littenberg on 11/29/23.
//

/**
@file glass_math.h
\brief Common math functions used throughout GLASS library
 */

#ifndef math_h
#define math_h

/**
 \brief Analytic in-place inversion of noise covariance matrix
 
 Substitutes contents of noise->Cij[n] with inverse, and sets deteriminent noise->detC[n]
 */
void invert_noise_covariance_matrix(struct Noise *noise);


/**
\brief Compute chirp mass from component masses
  
 @param m1 mass of primary
 @param m2 mass of secondary
 @return \f$ \mathcal{M} \equiv (m_1 m_2)^{3/5} (m_1+m_2)^{-1/5} \f$
 */
double chirpmass(double m1, double m2);

/**
 \brief Compute GW amplitude from intrinsic parameters
 
 @param Mc chirp mass: \f$\mathcal{M}\ [{\rm M}_\odot]\f$
 @param f0 initial GW frequency: \f$ f_0\ [{\rm Hz}]\f$
 @param D luminosity distance: \f$ D_L [{\rm pc}]\f$
 @return \f$ \mathcal{A} = 2 \frac{ \mathcal{M}^{5/3} (\pi f_0)^{2/3} }{D_L} \f$
 */
double amplitude(double Mc, double f0, double D);

/**
\brief Compute integer powers of by brute force multiplication
 
 Faster than pow() for small integers
  
 @param x variable
 @param n exponent
 @return \f$x^n =  x\times x \times\ ...\times x\f$
 */
double ipow(double x, int n);

/**
\brief Tukey window time series data
   
 @param data[in,out] time series to be windowed
 @param alpha[in] size of window filter in [s]
 @param N[in] size of time series
 */
void tukey(double *data, double alpha, int N);

/**
\brief Rearrange output of GSL RFT 
 
 Real and Imaginary components from GSL RFT functions are not ordered the way the rest of the package expects.
  
 @param x[out] array to be filled with ordered fourier coefficients
 @param x_gsl[in] input array with gsl-formatted courier coefficients
 @param N[in] size of arrays
 */
void unpack_gsl_rft_output(double *x, double *x_gsl, int N);

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
 @param N number of frequency bins in sum
 @return \f$(a|b) =  4 \sum_n a^*_n b_n C^-1_n \f$
 */
double fourier_nwip(double *a, double *b, double *invC, int N);
double wavelet_nwip(double *a, double *b, double *invC, int *list, int N);

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

/**
 \brief Wrapper to `GSL` cubic spline interpolation routines.

 @param[in] N number of spline points
 @param[in] x vector of independent-variable spline points
 @param[in] y vector of dependent-variable spline points
 @param[in] Nint number of interpolated points
 @param[out] xint vector of interpolated independent-variable points
 @param[out] yint vector of interpolated dependent-variable points
 */
void CubicSplineGSL(int N, double *x, double *y, int Nint, double *xint, double *yint);

/**
\brief GLASS implementation of DBSCAN clustering algorithm
 
 Density based clustering algorithm implemented here for 1D
 data with simple Euclidean distance measure
  
 @param[in] X set of data points to cluster
 @param[in] eps maximum distance between two samples to be considered neighbors
 @param[in] min minimum number of samples in a neighborhood to be considered a cluster
 @param[out] C cluster assignments mapping C[n] = M means X[n] is assigned to cluster M
 @param[out] K total number of clusters found
 */
void dbscan(gsl_vector *X, double eps, int min, int *C, int *K);

/** 
\brief Transform periodic to linear phase

 @param[in] N size of phase array
 @param[in,out] phase input from [0,2pi] replaced with unwrapped version
*/
void unwrap_phase(int N, double *phase);


double simpson_integration_3(double f0, double f1, double f2, double h);
double simpson_integration_5(double f0, double f1, double f2, double f3, double f4, double h);

/**
\brief wrapper for qsort() specific to integer arrays
*/
void integer_sort(int *x, int N);

void list_union(int *A, int *B, int NA, int NB, int *AUB, int *NAUB);

#endif /* math_h */
