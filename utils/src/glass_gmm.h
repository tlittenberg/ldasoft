/*
 * Copyright 2020 Tyson B. Littenberg
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
@file glass_gmm.h
\brief Library of functions for Gaussian Mixture Model fit to input data samples using Expectation-Maximization Algorithm.
 */

#ifndef gmm_h
#define gmm_h

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


/**
 * \brief Gaussian Mixture Model structure with metadata and pointer to mulitvariate Gaussian
 */
struct GMM
{
    size_t NParams;
    size_t NMODE;
    struct MVG **modes;
};

/**
 * \brief Parse GMM binary file and populate structre
 */
void read_gmm_binary(struct GMM *gmm, char filename[]);


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
    double *x; //!< location in parameter space
    double *p; //!< p(x) for each mode
    double *w; //!< weight sample for mode
};

/**
 * \brief Data structure for each multvariate Gaussian,
 * including parameters and covariance matrix products used in calculations
 */
struct MVG
{
    size_t size; //!< dimension of mvg
    double *mu; //!< means
    double **C; //!< covariance matrix
    double **L; //!< LU decomposition of C: \f$ C = LL^{-1}\f$
    double **Cinv; //!< inverse covariance matrix
    double **evectors; //!< eigenvectors
    double *evalues; //!< eigenvalues
    double **minmax; //!< min and max range for samples
    double detC; //!< determinant of covariance matrix
    double p; //!< prior for Mode (i.e. weighting)
    double Neff; //!< effective number of samples in mode
};

void alloc_MVG(struct MVG *mode, size_t N);

void free_MVG(struct MVG *mode);

/**
 * \brief This function writes the contents of MVG structure `mode` to the stream `fptr` in binary format.
 */
void write_MVG(struct MVG *mode, FILE *fptr);


/**
 * \brief This function reads the  contents of MVG structure `mode` from the stream `fptr` in binary format.
 */
void read_MVG(struct MVG *mode, FILE *fptr);

/**
 * \brief Deep copy of MVG structure
 */
void copy_MVG(struct MVG *origin, struct MVG *copy);

/**
 * \brief Evaluates the probability density of a multviariate Gaussian with input mean \f$\mu\f$
 * and covariance matrix \f$ C \f$.
 * \param[in] x vector of location to evaluate multivariate gaussian
 * \param[in] mvg data structure for multivariate gaussian
 * \param[out] probability \f$ \frac{\exp -\frac{1}{2}\left(x - \mu\right)^T C^{-1} \left(x - \mu\right)}{\sqrt{(2\pi)^N \det C}} \f$
 */
double multivariate_gaussian(double *x, struct MVG *mvg, int N);
double multivariate_gaussian_no_min(double *x, struct MVG *mvg);


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
void print_1D_pdfs(struct MVG **modes, struct Sample **samples, size_t NMODE, size_t NMCMC, char root[], size_t ix);

/**
 * \brief Print joint 2D distributions for each parameter to file
 */
void print_2D_pdfs(struct MVG **modes, struct Sample **samples, size_t NMODE, size_t NMCMC, char root[], size_t ix, size_t iy);

/**
 * \brief Print 1,2, and 3\f$\sigma\f$ contours of each individual
 *  Gaussian for in the model for each parameter pair
 */
void print_2D_contours(struct MVG **modes, size_t NMODE, char root[], size_t x1, size_t x2);

/**
 * \brief Wrapper for print_1D_pdfs() and print_2D_contours()
 */
void print_model(struct MVG **modes, struct Sample **samples, size_t NP, size_t NMODE, size_t NMCMC, double logL, double BIC, size_t step);


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
 * \param[in] NP size of parameter space
 * \param[in] NMODE number of modes
 * \param[in] NMCMC number of samples
 * \param[out] logL log-likelihood of input model
 * \param[out] BIC Bayesian Information Criteria (BIC) for input model
 * \param[in,out] modes parameters of each Gaussian including relative weights \f$\alpha_{k}\f$, updated by M-step.
 * \return `0` if successful, `1` if singular due to zero weight in one of the modes.
 */
int expectation_maximization(struct Sample **samples, struct MVG **modes, size_t NP, size_t NMODE, size_t NMCMC, double *logL, double *BIC);

/**
 * \brief Gaussian Mixture Model fit with Expectation-Maximization (EM) Algorithm
 *
 * Full Gaussian Mixture model fit to input `samples` using `EM` model
 *
 * \param[in] modes data structure of multivariate Gaussians in GMM
 * \param[in] samples data structure with all chain samples to be fit
 * \param[in] NP size of parameter space
 * \param[in] NMODE number of modes
 * \param[in] NMCMC number of chain samples
 * \param[in] NSTEP number of iterations of EM algorithm
 * \param[in] r random number generator seed
 * \param[out] logL log-likelihood of input model
 * \param[out] BIC Bayesian Information Criteria (BIC) for input model
 * \return `0` if successfule, `1` if singular due to zero weight in one of the modes
*/
int GMM_with_EM(struct MVG **modes, struct Sample **samples, size_t NP, size_t NMODE, size_t NMCMC, size_t NSTEP, unsigned int *r, double *logL, double *BIC);


double logit(double x,double xmin,double xmax);
double sigmoid(double x,double xmin,double xmax);
double dsigmoid(double x, double xmin, double xmax);
void logit_mapping(double *x_vec, double *y_vec, double xmin, double xmax, int N);
void sigmoid_mapping(double *x_vec, double *y_vec, double xmin, double xmax, int N);

#endif /* gmm_h */
