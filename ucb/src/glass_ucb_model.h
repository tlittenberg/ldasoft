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
 @file glass_ucb_model.h
 \brief Functions for defining, manipulating, and evaluating UCB model
 */


#ifndef ucb_model_h
#define ucb_model_h

#include <stdio.h>

#define UCB_MODEL_NP 8 ///< Number of source parameters for UCB model

/**
\brief Hierarchical structure of UCB model
 */
struct Model
{
    ///@name Source parameters
    ///@{
    int Nmax;   //!<maximum number of signals in model
    int Neff;   //!<effective maximum number of signals during burn in
    int Nlive;  //!<current number of signals in model
    struct Source **source; //!<source structures for each signal in the model
    ///@}
    
    /// Noise parameters
    struct Noise *noise;
    
    /// Calibration parameters
    struct Calibration *calibration;
    
    ///@name TDI structures
    ///@{
    struct TDI *tdi; //!<joint signal model
    struct TDI *residual; //!<joint residual
    ///@}
    
    ///@name Segment start time
    ///@{
    double t0; //!<start time
    ///@}
    
    ///@name Source parameter priors
    ///@{
    double **prior; //!<upper and lower bounds for uniform priors \f$ [\theta_{\rm min},\theta_{\rm max}]\f$
    double *logPriorVolume; //!<prior volume \f$ -\Sum \log(\theta_{\rm max}-\theta_{\rm min})\f$
    ///@}
    
    ///@name Model likelihood
    ///@{
    double logL; //!<unnormalized log likelihood \f$ -(d-h|d-h)/2 \f$
    double logLnorm; //!<normalization of log likelihood \f$ \propto -\log \det C \f$
    ///@}

    ///@name Wavelet bookkeeping
    ///@{
    int *list; //!<list of active wavelet pixels
    int Nlist; //!<number of active wavelet pixels
    ///@}
};

/**
\brief Structure containing parameters and meta data for a single galactic binary.
*/
struct Source
{
    /// Array containing parameters to be passed to ucb_waveform()
    double *params;

    /// Instrument response to signal with Source::params \f$ h(\vec\theta) \f$
    struct TDI *tdi;
    

    ///@name Intrinsic Parameters
    ///@{
    double f0;    //!< Gravitational wave frequency at start of observations \f$f_0 = 2/P\f$
    double dfdt;  //!< First time derivitive of GW frequemcy \f$ \frac{df}{dt}\f$, see galactic_binary_fdot()
    double d2fdt2;//!< Second time derivitive of GW frequemcy \f$ \frac{d^2f}{dt^2}\f$
    double amp;   //!< Gravitational wave amplitude \f$ \mathcal{A} = 2 \frac{\mathcal{M}^{5/3} (\pi f)^{2/3}}{D_L} \f$, see galactic_binary_Amp()
    ///@}
    
    ///@name Extrinsic Parameters
    ///@{
    double psi;      //!< Polarization angle
    double cosi;     //!< Cosine of orbital inclination angle
    double phi0;     //!< Gravitational wave phase at start of observations
    double phi;      //!< Ecliptic longitude
    double costheta; //!< Cosine of ecliptic co-latitude
    ///@}
    
    ///@name Derived Parameters
    ///See data.c for functions to convert between intrinsic and derived parameters
    ///@{
    double m1;  //!< Primary component mass \f$m_1 \geq m_2\f$
    double m2;  //!< Secondary component mass \f$m_2 \leq m_1\f$
    double Mc;  //!< Chirp mass \f$\mathcal{M}\f$, see galactic_binary_Mc()
    double D;   //!< Luminosity Distance \f$D_L\f$, see galactic_binary_dL()
    ///@}
    
    ///@name Template waveform alignment
    ///See galactic_binary_alignment()
    ///@{
    int BW;   //!< Signal bandwidth in frequency bins, see galactic_binary_bandwidth()
    int qmin; //!< Minimum frequency bin of signal (relative to 0 Hz)
    int qmax; //!< Maximum frequency bin of signal (relative to 0 Hz)
    int imin; //!< Minimum frequency bin of signal relative to segment start
    int imax; //!< Maximum frequency bin of signal relative to segment start
    ///@}

    ///@name Fisher Information Matrix
    ///See galactic_binary_fisher()
    ///@{
    double **fisher_matrix; //!<Fisher approximation to inverse covariance matrix
    double **fisher_evectr; //!<Eigenvectors of covariance matrix
    double *fisher_evalue;  //!<Eigenvalues of covariance matrix
    int fisher_update_flag; //!<1 if fisher needs update, 0 if not
    ///@}

    ///@name Wavelet bookkeeping
    ///@{
    int *list; //!<list of active wavelet pixels
    int Nlist; //!<number of active wavelet pixels
    ///@}

};

/**
\brief Create galactic binary model waveform

 Computes LISA response to signal with parameters indexed by `source_id`
 in the Model::Source, and creates meta template of all sources in the model.
 */
void generate_signal_model(struct Orbit *orbit, struct Data *data, struct Model *model, int source_id);
void generate_signal_model_wavelet(struct Orbit *orbit, struct Data *data, struct Model *model, int source_id);

/**
\brief Modify galactic binary model waveform

 Computes LISA response to signal with parameters indexed by `source_id`
 in the Model::Source, and updates meta template of all sources in the model.
 */
void update_signal_model(struct Orbit *orbit, struct Data *data, struct Model *model_x, struct Model *model_y, int source_id);
void update_signal_model_wavelet(struct Orbit *orbit, struct Data *data, struct Model *model_x, struct Model *model_y, int source_id);

/**
\brief F-statistic maximization of galactic binary parameters

 Wrapper of UCBFstatistic.c functions for maximizing waveform
 over \f$(\mathcal{A},\cos\iota,\psi,\varphi_0)\f$ by filtering on
 original data.
 
 \todo test filtering on residuals
 */
void maximize_signal_model(struct Orbit *orbit, struct Data *data, struct Model *model, int source_id);

/**
 \brief Create LISA instrument noise model
 
 Computes \f$S_n(f)\f$ from Model::noise.
 */
void generate_noise_model(struct Data *data, struct Model *model);
void generate_noise_model_wavelet(struct Data *data, struct Model *model);

/**
 \brief Create LISA instrument calibration model
 
 Computes phase and amplitude corrections from Model::calibration parameters.
 */
void generate_calibration_model(struct Data *data, struct Model *model);

/**
\brief Apply amplitude and phase corrections.

Computes new LISA instrument response Model::tdi after applying
 amplitude and phase corrections from calibration parameters.
 */
void apply_calibration_model(struct Data *data, struct Model *model);

/**
 \brief Compute argument of Gaussian likelihood
 
 Computes residual of data and meta-template from Model and noise weighted inner product.
 @return \f$ -\frac{1}{2}(d-h|d-h) \f$
 */
double gaussian_log_likelihood(struct Data *data, struct Model *model);

/**
 \brief Compute normalization of Gaussian likelihood for constant noise level
 
 For noise models that are (approximately) a constant over the band
 of \f$N\f$ Fourier bins, parameterized by a multiplyer \f$\eta\f$
 
 @return \f$ \sum_{\rm TDI} N\log\eta_{\rm TDI} \f$
 */
double gaussian_log_likelihood_constant_norm(struct Data *data, struct Model *model);

/**
 \brief Compute normalization of Gaussian likelihood for arbitrary noise level
 
 For noise models that are free to vary over the band
 of \f$N\f$ Fourier bins
 
 @return \f$ \sum_{\rm TDI} \sum_f \log S_{n,{\rm TDI}}(f) \f$
 */
double gaussian_log_likelihood_model_norm(struct Data *data, struct Model *model);

/**
 \brief Compute difference in log Likelihood from changing parameters of one source
 
 Updates residuals from changing single source `source_id` and computes change in likelihood sum only over frequency range where waveforms were non-zero.
 @return \f$ -\frac{1}{2}\left[(d-h_{\rm new}|d-h__{\rm new}) - (d-h_{\rm old}|d-h_{\rm old})\right] \f$
 */
double delta_log_likelihood(struct Data *data, struct Model *model_x, struct Model *model_y, int source_id);

double gaussian_log_likelhood_wavelet(struct Data *data, struct Model *model);

/**
 \brief Check for increase in maximum log likelihood
 */
int update_max_log_likelihood(struct Model **model, struct Chain *chain, struct Flags *flags);

/**
 \brief Converts physical UCB parameters to array expected by ucb_waveform.c
 */
void map_params_to_array(struct Source *source, double *params, double T);

/**
 \brief Converts array expected by ucb_waveform.c to
 physical UCB parameters
 */
void map_array_to_params(struct Source *source, double *params, double T);



/** @name Allocate memory for structures */
///@{
void alloc_model(struct Data *data, struct Model *model, int Nmax);
void alloc_source(struct Source *source, int N, int Nchannel);
///@}

/**
 \brief Shallow copy of Data structure
 */
void copy_data(struct Data *origin, struct Data *copy);

/** @name Deep copy structure contents */
///@{
void copy_source(struct Source *origin, struct Source *copy);
void copy_model(struct Model *origin, struct Model *copy);
void copy_model_lite(struct Model *origin, struct Model *copy);
///@}

/** @name Free memory for structures */
///@{
void free_model(struct Model *model);
void free_source(struct Source *source);
///@}

/**
 \brief Deep comparison of Model contents
 
 @return 0 if models are identical
 @return 1 if models are different
 */
int compare_model(struct Model *a, struct Model *b);

void remove_signal_model(struct Data *data, struct Model *model, struct Source *source);
void remove_signal_model_wavelet(struct Data *data, struct Model *model, struct Source *source);

void add_signal_model(struct Data *data, struct Model *model, struct Source *source);
void add_signal_model_wavelet(struct Data *data, struct Model *model, struct Source *source);

void waveform_check(struct Orbit *orbit, struct Data *data, struct Model *model, struct Source *source);

#endif /* ucb_model_h */
