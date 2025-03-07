/*
 *  Copyright (C) 2023 Tyson B. Littenberg (MSFC-ST12)
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
 @file glass_noise_model.h
 \brief Functions and structures defining noise models.

 */
#ifndef noise_model_h
#define noise_model_h

#define MIN_SPLINE_STENCIL 5
#define MIN_SPLINE_SPACING 256.

/*!
 * \brief Data strictire for spline model used in noise model.
 *
 * The SplineModel structure stores metadata and state of spline model,
 * including size of model, current state of spline control points
 * and likelihood, and interpolated model.
 */
struct SplineModel
{
    int Nmin; //!< Minimum number of spline control points
    int Nmax; //!< Maximum number of spline control points
    int Nchannel; //!< Number of TDI channels
    double logL; //!< Log likelhood of spline model
    struct Noise *spline; //!< Spline control points
    struct Noise *psd; //!< Reconstructed noise model
};

struct InstrumentModel
{
    int Nlink;         //!< number of interferometer links (6 for 3-spacecraft constellation)
    double logL;       //!< log Likelihood of model
    double *soms;      //!< optical metrology system noise parameters
    double *sacc;      //!< acceleration noise parameters
    struct Noise *psd; //!< power and cross spectral densities
    
    /** @name Link level noise parameters */
     ///@{
    double sacc12;
    double sacc13;
    double sacc21;
    double sacc23;
    double sacc31;
    double sacc32;
    double soms12;
    double soms13;
    double soms21;
    double soms23;
    double soms31;
    double soms32;
    ///@}
};

struct ForegroundModel
{
    int Nparams;       //!< number of foreground parameters (7)
    double Tobs;       //!< observation time (in seconds)
    double *sgal;      //!< galactic foreground parameters
    double logL;       //!< log Likelihood of model
    struct Noise *psd; //!< power and cross spectral densities

    /** @name galactic foreground parameters from Digman & Cornish 10.3847/1538-4357/ac9139*/
     ///@{
    double Amp;   //!< overall amplitude  (~2.13e-37)
    double f1;
    double alpha; //!< spectral index of low frequency part (~1.58)
    double fk;
    double f2;    //!< some other frequency parameter? (~0.53 mHz)
     ///@}

    struct GalaxyModulation *modulation; //!< time-series of modulation
};

/**
 \brief Converts physical noise parameters to array expected by InstrumentModel
 */
void map_noise_params_to_array(struct InstrumentModel *model);

/**
 \brief Converts array expected by InstrumentModel to
 physical noise parameters
 */
void map_array_to_noise_params(struct InstrumentModel *model);

/**
 \brief Converts phenomenological foreground parameters to array expected by ForegroundModel
 */
void map_foreground_params_to_array(struct ForegroundModel *model);

/**
 \brief Converts array expected by ForegroundModel to
 phenomenological foreground parameters
 */
void map_array_to_foreground_params(struct ForegroundModel *model);

/**
 \brief Allocates spline model structure and contents.
 */
void alloc_spline_model(struct SplineModel *model, int Ndata, int Nlayer, int Nchannel, int Nspline);

/**
 \brief Allocates instrument model structure and contents.
 */
void alloc_instrument_model(struct InstrumentModel *model, int Ndata, int Nlayer, int Nchannel);

/**
 \brief Allocates galactic foreground model structure and contents.
 */
void alloc_foreground_model(struct ForegroundModel *model, int Ndata, int Nlayer, int Nchannel);

/**
 \brief Free allocated spline model.
 */
void free_spline_model(struct SplineModel *model);

/**
 \brief Free allocated spline model.
 */
void free_instrument_model(struct InstrumentModel *model);

/**
 \brief Free allocated spline model.
 */
void free_foreground_model(struct ForegroundModel *model);

/**
 \brief Deep copy of SplineModel structure from `origin` into `copy`
 */
void copy_spline_model(struct SplineModel *origin, struct SplineModel *copy);

/**
 \brief Deep copy of InstrumentModel structure from `origin` into `copy`
 */
void copy_instrument_model(struct InstrumentModel *origin, struct InstrumentModel *copy);

/**
 \brief Deep copy of ForegroundModel structure from `origin` into `copy`
 */
void copy_foreground_model(struct ForegroundModel *origin, struct ForegroundModel *copy);

/**
 \brief Wrapper to `CubicSplineGSL` functions for generating PSD model based on current state of `model`
 */
void generate_spline_noise_model(struct SplineModel *model);

/**
 \brief Compute instrument model contribution to noise covariance matrix based on current state of `model`
 */
void generate_instrument_noise_model(struct Orbit *orbit, struct InstrumentModel *model);
void generate_instrument_noise_model_wavelet(struct Wavelets *wdm, struct Orbit *orbit, struct InstrumentModel *model);

/**
 \brief Compute galactic foreground contribution to covariance matrix based on current state of `model`
 */
void generate_galactic_foreground_model(struct ForegroundModel *model);
void generate_galactic_foreground_model_wavelet(struct Wavelets *wdm, struct ForegroundModel *model);

/**
 \brief Add components to `full` noise covariance matrix `C`
 */
void generate_full_covariance_matrix(struct Noise *full, struct Noise *component, int Nchannel);
void generate_full_dynamic_covariance_matrix(struct Wavelets *wdm, struct InstrumentModel *inst, struct ForegroundModel *conf, struct Noise *full);

/**
\brief Compute spline model only where interpolant changes
 
 Interpolates spline points only in the vicinity of `new_knot`
 */
void update_spline_noise_model(struct SplineModel *model, int new_knot, int min_knot, int max_knot);

/**
 \brief Log likelihood for noise model.
 
 @param data Data structure
 @param model SplineModel structure containing current state
 @return \f$  \ln p({\rm data}|{\rm spline}) \f$
 */
double noise_log_likelihood(struct Data *data, struct Noise *noise);

/**
 \brief Change in log likelihood for noise model.
 
 @param data Data structure
 @param model_x Current SplineModel structure containing current state
 @param model_y Current SplineModel structure containing current state
 @return \f$  \ln p({\rm data}|{\rm spline}_y)-\ln p({\rm data}|{\rm spline}_x) \f$
 */
double noise_delta_log_likelihood(struct Data *data, struct SplineModel *model_x, struct SplineModel *model_y, double fmin, double fmax, int ic);

/**
 \brief Set initial state of spline `model`
 */
void initialize_spline_model(struct Orbit *orbit, struct Data *data, struct SplineModel *model, int Nspline);

/**
 \brief Set initial state of instrument noise `model`
 */
void initialize_instrument_model(struct Orbit *orbit, struct Data *data, struct InstrumentModel *model);
void initialize_instrument_model_wavelet(struct Orbit *orbit, struct Data *data, struct InstrumentModel *model);
/**
 \brief Set initial state of instrument noise `model`
 */
void initialize_foreground_model(struct Orbit *orbit, struct Data *data, struct ForegroundModel *model);
void initialize_foreground_model_wavelet(struct Orbit *orbit, struct Data *data, struct ForegroundModel *model);

void GetDynamicNoiseModel(struct Data *data, struct Orbit *orbit, struct Flags *flags);
void GetStationaryNoiseModel(struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Noise *noise);

#endif /* noise_model_h */
