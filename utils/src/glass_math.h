//
//  math.h
//  utils
//
//  Created by Tyson Littenberg on 11/29/23.
//

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
\brief Compute integer powers of by brute force multiplication
 
 Faster than pow() for small integers
  
 @param x variable
 @param n exponent
 @return \f$x^n =  x\times x \times\ ...\times x\f$
 */
double ipow(double x, int n);

void tukey(double *data, double alpha, int N);
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
 @param n number of frequency bins in sum
 @return \f$(a|b) =  4 \sum_n a^*_n b_n C^-1_n \f$
 */
double fourier_nwip(double *a, double *b, double *invC, int n);

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
  
 @param X set of data points to cluster
 @param eps maximum distance between two samples to be considered neighbors
 @param min minimum number of samples in a neighborhood to be considered a cluster
 @param C cluster assignments mapping C[n] = M means X[n] is assigned to cluster M
 @param K total number of clusters found
 */
void dbscan(gsl_vector *X, double eps, int min, int *C, int *K);


#endif /* math_h */
