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
 @file GalacticBinaryPrior.h
 \brief Functions supporting prior distributions.
 
 Includes functions for createing and evaluating model priors.
 */


#ifndef GalacticBinaryPrior_h
#define GalacticBinaryPrior_h

///@name Galaxy Prior
///@{
#define GALAXY_RGC 7.2 //!<distance from solar BC to GC (kpc)
#define GALAXY_A  0.25 //!<bulge fraction
#define GALAXY_Rb 0.8  //!< bulge radius (kpc)
#define GALAXY_Rd 2.5  //!< disk radius (kpc)
#define GALAXY_Zd 0.4  //!< disk height (kpc)
///@}

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

    ///@name Gaussian Mixture Model prior
    ///@{
    size_t NP; //!<number of parameters in mixture model
    size_t NMODE; //!<number of modes in mixture model
    struct MVG **modes; //!<data structure for multivariate Gaussian
    ///@}
    
    ///@name workspace
    ///@{
    double *vector;  //!<utility 1D array for prior metadata
    double **matrix; //!<utility 2D array for prior metadata
    double ***tensor;//!<utility 3D array for prior metadata
    ///@}
};

int check_range(double *params, double **uniform_prior, int NP);
void set_galaxy_prior(struct Flags *flags, struct Prior *prior);
void set_gmm_prior(struct Flags *flags, struct Data *data, struct Prior *prior);
void set_uniform_prior(struct Flags *flags, struct Model *model, struct Data *data, int verbose);
double evaluate_prior(struct Flags *flags, struct Data *data, struct Model *model, struct Prior *prior, double *params);
double evaluate_snr_prior(struct Data *data, struct Model *model, double *params);
double evalaute_sky_location_prior(double *params, double **uniform_prior, double *logPriorVolume, int galaxyFlag, double *skyhist, double dcostheta, double dphi, int nphi);
double evaluate_uniform_priors(double *params, double **uniform_prior, double *logPriorVolume, int NP);
double evaluate_gmm_prior(struct Data *data, struct MVG **modes, int NMODES, double *params);



#endif /* GalacticBinaryPrior_h */
