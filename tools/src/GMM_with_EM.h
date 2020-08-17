/*
 *  Author: Tyson B. Littenberg (MSFC-ST12)
 *  Created: 07.27.2020
 *
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
 *
 *  Headers for Gaussian Mixture Model using the Expectation-Maximization Algorithm
 *
 */

/**
@file GMM_with_EM.h
\brief Library of functions for Gaussian Mixture Model fit to input data samples using Expectation-Maximization Algorithm.
 */

#ifndef GMM_with_EM_h
#define GMM_with_EM_h

/** @name Flags for command line parser */
///@{
#define no_arg 0   //!< no argument
#define req_arg 1  //!< required argument
#define opt_arg 2  //!< optional argument
///@}

/** @name Boundaries for variables */
///@{
#define BUFFER_SIZE 1024 //!< max size of `char` buffers
#define PMIN 1.e-16 //!< floor on probabilitiy densities (avoid singularities at \f$p=0\f$
///@}

/** @name Progress bar settings */
///@{
#define PBSTR "||||||||||||||||||||||||||||||||||||||||" //!< displayed progress
#define PBWIDTH 40 //!<size of progress bar
///@}



/**
 * \brief Show usage
 */
void printUsage(const char *program);

/**
 * \brief Data structure for individual samples,
 * along with their relative weights and probabilites in each mode.
 */
struct Sample
{
    gsl_vector *x; //!< location in parameter space
    gsl_vector *p; //!< p(x) for each mode
    gsl_vector *w; //!< weight sample for mode
};

/**
 * \brief Data structure for each multvariate Gaussian,
 * including parameters and covariance matrix products used in calculations
 */
struct MVG
{
    size_t size; //!< dimension of mvg
    gsl_vector *mu; //!< means
    gsl_matrix *C; //!< covariance matrix
    gsl_matrix *L; //!< LU decomposition of C: \f$ C = LL^{-1}\f$
    gsl_matrix *Cinv; //!< inverse covariance matrix
    gsl_matrix *evectors; //!< eigenvectors
    gsl_vector *evalues; //!< eigenvalues
    double detC; //!< determinant of covariance matrix
    double p; //!< prior for Mode (i.e. weighting)
    double Neff; //!< effective number of samples in mode
};

void alloc_MVG(struct MVG *mode, size_t N);

void free_MVG(struct MVG *mode);

/**
 * \brief Evaluates the probability density of a multviariate Gaussian with input mean \f$\mu\f$
 * and covariance matrix \f$ C \f$.
 * \param[in] x vector of location to evaluate multivariate gaussian
 * \param[in] mvg data structure for multivariate gaussian
 * \param[out] probability \f$ \frac{\exp -\frac{1}{2}\left(x - \mu\right)^T C^{-1} \left(x - \mu\right)}{\sqrt{(2\pi)^N \det C}} \f$
 */
double multivariate_gaussian(gsl_vector *x, struct MVG *mvg);

/**
 * \brief Wrapper for computing matrix inversion using
 * `GSL` functions.
 * \param[in] A input symmetric matrix
 * \param[out] Ainv inverse of `A`
 * \param[out] L LU decomposition of `A`
 * \param[out] detA determinant of `A`
 * \param[out] R recipricol condition number [0,1]
 */
void invert_gsl_matrix(gsl_matrix *A, gsl_matrix *Ainv, gsl_matrix *L, double *detA, double *R);

/**
 * \brief Wrapper for computing matrix eigenvalue decomposition using
 * `GSL` functions.
 * \param[in] A input symmetric matrix
 * \param[out] evec matrix where each row has components of eigenvector
 * \param[out] eval vector where each element is the eigenvalue associated with the same row of `evec`.
 */
void decompose_matrix(gsl_matrix *A, gsl_matrix *evec, gsl_vector *eval);

/**
 * \brief Log-likelihood of Gaussian Mixture Model
 * \param[in] modes parameters of each Gaussian including relative weights \f$\alpha_{k}\f$
 * \param[in] samples data points \f$x_i\f$ and probabilities \f$p(x_i|k)\f$ for each Gaussian \f$k\f$
 * \param[in] NMCMC number of samples
 * \param[in] NMODE number of modes
 * \param[out] log-likelihood \f$\log L = \sum_k^{\rm NMCMC} \log \sum_i^{\rm NMODE} \alpha_k p(x_i | k)\f$
 *
 */
double log_likelihood(struct MVG **modes, struct Sample **samples, int NMCMC, int NMODE);

/**
 * \brief Print joint 1D distributions for each parameter to file
 */
void print_1D_pdfs(struct MVG **modes, struct Sample **samples, size_t NMCMC, char root[], size_t ix);

/**
 * \brief Print 1,2, and 3\f$\sigma\f$ contours of each individual
 *  Gaussian for in the model for each parameter pair
 */
void print_2D_contours(struct MVG **modes, size_t NMODE, char root[], size_t x1, size_t x2);

/**
 * \brief Wrapper for print_1D_pdfs() and print_2D_contours()
 */
void print_model(struct MVG **modes, struct Sample **samples, size_t NMCMC, double logL, double BIC, size_t step);


/**
 * \brief The Expectation-Maximization (EM) Algorithm
 *
 * Execulte one iteration of the EM algorithm to update fit to multivariate gaussian model (could be generalized)
 *
 * *E-step*: Compute probability for each sample to belong to each mode
 *
 * *M-step*: Recompute mean, covariance, and relative weight of newly weighted samples
 *
 * \param[in] samples data points \f$x_i\f$ and probabilities \f$p(x_i|k)\f$ for each Gaussian \f$k\f$
 * \param[in] NMCMC number of samples
 * \param[out] logL log-likelihood of input model
 * \param[out] BIC Bayesian Information Criteria (BIC) for input model
 * \param[in,out] modes parameters of each Gaussian including relative weights \f$\alpha_{k}\f$, updated by M-step.
 */
void expectation_maximization(struct Sample **samples, struct MVG **modes, size_t NMCMC, double *logL, double *BIC);

/**
 * \brief Gaussian Mixture Model fit with Expectation-Maximization (EM) Algorithm
 *
 * Full Gaussian Mixture model fit to input `samples` using `EM` model
 *
 * \param[in] modes data structure of multivariate Gaussians in GMM
 * \param[in] samples data structure with all chain samples to be fit
 * \param[in] NMCMC number of chain samples
 * \param[in] NSTEP number of iterations of EM algorithm
 * \param[in] r `GSL` random number generator seed
 * \param[out] logL log-likelihood of input model
 * \param[out] BIC Bayesian Information Criteria (BIC) for input model
*/
void GMM_with_EM(struct MVG **modes, struct Sample **samples, size_t NMCMC, size_t NSTEP, gsl_rng *r, double *logL, double *BIC);

#endif /* GMM_with_EM_h */
