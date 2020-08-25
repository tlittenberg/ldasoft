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


#ifndef GalacticBinaryPrior_h
#define GalacticBinaryPrior_h

/* ----------------  MILKY WAY MODEL  ------------------- */

/* distance from solar BC to GC (kpc) */
#define GALAXY_RGC 7.2

/* bulge fraction */
#define GALAXY_A  0.25

/* bulge radius (kpc) */
#define GALAXY_Rb 0.8

/* disk radius (kpc) */
#define GALAXY_Rd 2.5

/* disk height (kpc) */
#define GALAXY_Zd 0.4

/* ----------------  CALIBRATION MODEL  ------------------- */

/* stddev in phase error (radians) */
//#define CAL_SIGMA_PHASE 0.175 // ~10^deg
#define CAL_SIGMA_PHASE 0.35 // ~30^deg

/* stddev in fractional amplitude error */
//#define CAL_SIGMA_AMP 0.1 // 10%
#define CAL_SIGMA_AMP 0.20 // 20%

/* ----------------  MISC  ------------------- */

struct Prior
{
    //uniform prior
    double **prior;
    double logPriorVolume;
    
    //galaxy prior
    double *skyhist;
    double dcostheta;
    double dphi;
    double skymaxp;
    int ncostheta;
    int nphi;
        
    //gaussian mixture model prior
    size_t NP;
    size_t NMODE;
    struct MVG **modes; //!<data structure for multivariate Gaussian

    //workspace
    double *vector;  //!<utility 1D array for prior metadata
    double **matrix; //!<utility 2D array for prior metadata
    double ***tensor;//!<utility 3D array for prior metadata

};

int check_range(double *params, double **uniform_prior, int NP);
void set_galaxy_prior(struct Flags *flags, struct Prior *prior);
void set_gmm_prior(struct Flags *flags, struct Data *data, struct Prior *prior);
void set_uniform_prior(struct Flags *flags, struct Model *model, struct Data *data, int verbose);
double evaluate_prior(struct Flags *flags, struct Data *data, struct Model *model, struct Prior *prior, double *params);
double evaluate_snr_prior(struct Data *data, struct Model *model, double *params);
double evalaute_sky_location_prior(double *params, double **uniform_prior, double *logPriorVolume, int galaxyFlag, double *skyhist, double dcostheta, double dphi, int nphi);
double evaluate_uniform_priors(double *params, double **uniform_prior, double *logPriorVolume, int NP);
double evaluate_gmm_prior(struct Data *data, struct Prior *prior, double *params);



#endif /* GalacticBinaryPrior_h */
