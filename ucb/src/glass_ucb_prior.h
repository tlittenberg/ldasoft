/*
 * Copyright 2019 Tyson B. Littenberg & Neil J. Cornish
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
 @file glass_ucb_prior.h
 \brief Functions supporting prior distributions.
 
 Includes functions for createing and evaluating model priors.
 */


#ifndef ucb_prior_h
#define ucb_prior_h


///@name Calibration prior
///@{
#define CAL_SIGMA_PHASE 0.35 //!< 1-\f$\sigma\f$ phase error (rad) \f${\sim}30^\circ\f$
#define CAL_SIGMA_AMP 0.20 //!< 1-\f$\sigma\f$ fractional amplitude error
///@}

/*!
 \brief Prototype structure for prior distributions.
 
 Generic data structure for holding all information needed by prior distributions.
 Structure contains parameters for different supported priors and various book-keeping scalars, vectors, and matrices to
 hold needed metadata.

*/
struct Prior
{    
    ///@name Uniform prior
    ///@{
    double **prior; //!<upper and lower bounds for uniform priors \f$ [\theta_{\rm min},\theta_{\rm max}]\f$
    double logPriorVolume; //!<prior volume \f$ -\sum \log(\theta_{\rm max}-\theta_{\rm min})\f$
    ///@}

    ///@name Uniform prior
    ///@{
    double *skyhist; //!<2D histogram of prior density on sky
    double dcostheta; //!<size of `skyhist` bins in \f$\cos\theta\f$ direction
    double dphi; //!<size of `skyhist` bins in \f$\phi\f$ direction
    double skymaxp; //!<max prior density of `skyhist`
    int ncostheta; //!<number of `skyhist` bins in \f$\cos\theta\f$ direction
    int nphi; //!<number of `skyhist` bins in \f$\phi\f$ direction
    ///@}
    
    ///@name workspace
    ///@{
    double *vector;  //!<utility 1D array for prior metadata
    double **matrix; //!<utility 2D array for prior metadata
    double ***tensor;//!<utility 3D array for prior metadata
    ///@}

    /// Gaussian Mixture Model prior
    struct GMM *gmm;
};
/**
 \brief Checks that parameters \f$\vec x\f$ are within prior volume \f$V\f$.
 
 @returns 0 if \f$\vec x \in V\f$
 @returns 1 else
 */
int check_range(double *params, double **uniform_prior);

/**
\brief Set up sky location prior assuming galactic distribution
 
 Uses axially symmetric disk and bulge model of galaxy,
 randomly distributes points within the assumed distribution,
 and bins the points based on their sky location from Earth
 to use as a prior.
 */
void set_galaxy_prior(struct Flags *flags, struct Prior *prior);

/**
\brief Sets Gaussian Mixture Model prior for source model
 */
void set_gmm_prior(struct Flags *flags, struct Data *data, struct Prior *prior, struct Catalog *catalog);

/**
 \brief Sets Uniform prior for source model
 */
void set_uniform_prior(struct Flags *flags, struct Model *model, struct Data *data, int verbose);

/**
 \brief Computes joint prior for input parameters `params`
 
 @param UCB parameters `params` \f$ \vec x\f$
 @returns \f$ \log p(\vec x)\f$
 */
double evaluate_prior(struct Flags *flags, struct Data *data, struct Model *model, struct Prior *prior, double *params);

/**
 \brief Computes prior for sky location parameters \f$\vec\Omega\f$
 
 Depending on Flags::galaxyFlag, evaluates either uniform or galaxy prior
 for sky location parameters
 
 @param UCB parameters `params` \f$ \vec x\f$
 @returns \f$ \log p({\vec\Omega})\f$
 */
double evalaute_sky_location_prior(double *params, double **uniform_prior, double *logPriorVolume, int galaxyFlag, double *skyhist, double dcostheta, double dphi, int nphi);

/**
 \brief Computes uniform prior for parameters \f$\vec x\f$

 \f$ p(\vec x) = \prod_i \frac{1}{\Delta x_i} \f$
 
 @param UCB parameters `params` \f$ \vec x\f$
 @returns \f$ \log p({\vec x})\f$
 */
double evaluate_uniform_priors(double *params, double **uniform_prior, double *logPriorVolume);
double evaluate_gmm_prior(struct Data *data, struct GMM *gmm, double *params);



#endif /* ucb_prior_h */
