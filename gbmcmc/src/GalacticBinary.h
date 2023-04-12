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


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <util.h>

/**
@file GalacticBinary.h
\brief Definition of data structures.
*/

#define MAXSTRINGSIZE 1024 //!<maximum number of characters for a line in a data file.

/*!
 * \brief Analaysis segment and meta data about size of segment, location in full data stream, and LISA observation parameters.
 *
 * The Data structure stores information about the data being analyzed,
 * including size of data, where it is located in the full LISA band,
 * information about gaps, etc.
 *
 * The structure stores the Fourier-domain TDI data channels, and a
 * fiducial model for the instrument noise. It also has memory allocated
 * for storing the reconstructed waveforms, residuals, and noise models.
 */
struct Data
{
    /** @name Size of Data and Model */
     ///@{
    int N;        //!<number of frequency bins
    int Nmax;     //!<max size of frequency segment (for global fit)
    int NT;       //!<number of time segments
    int Nchannel; //!<number of data channels
    int DMAX;     //!<max dimension of signal model
    double logN;  //!<log total number of data points \f$ \log ( 2 \times N \times N_{\rm channel} \times N_{\rm T} )\f$
    ///@}

    /** @name Random Number Generator Seeds */
     ///@{
    long cseed; //!<seed for MCMC set by `[--chainseed=INT; default=150914]`
    long nseed; //!<seed for noise realization set by `[--noiseseed=INT; default=151226]`
    long iseed; //!<seed for injection parameters set by `[--injseed=INT; default=151012]`
    ///@}

    /** @name Time of Observations */
     ///@{
    double T; //!<observation time
    double sqT; //!<\f$\sqrt{T}\f$
    double *t0;   //!<start times of segments
    double *tgap; //!<time between segments
    ///@}

    
    /** @name Analysis Frequency Segment */
     ///@{
    int qmin; //!<minimum frequency bin of segment
    int qmax; //!<maximum frequency bin of segment
    int qpad; //!<number of frequency bins padding ends of segment
    
    double fmin; //!<minimum frequency of segment
    double fmax; //!<maximum frequency of segment
    double sine_f_on_fstar; //!<\f$sin(f * 2\pi L/c)\f$

    //some manipulations of f,fmin for likelihood calculation
    double sum_log_f; //!<\f$\sum \log(f)\f$ appears in some normalizations
    double logfmin; //!<\f$\log(f_{\rm min})\f$ appears in some normalizations
    double logfmax; //!<\f$\log(f_{\rm max})\f$ appears spline setup for bayesline
    ///@}

    /** @name TDI Data and Noise */
     ///@{
    struct TDI **tdi; //!<TDI data channels as seen by sampler
    struct TDI **raw; //!<TDI data channels unaltered from input
    struct Noise **noise; //!<Reference noise model
    /**
     \brief Convention for data format
     
     format = "phase" for phase difference (distance)
     format = "frequency" for fractional frequency (velocity) **Use for matching LDC Radler**
     format = "sangria" for fractional frequency w/ LDC Sangria-era TDI & phase conventions
     */
    char format[16];
    char dataDir[PATH_BUFSIZE]; //!<Directory for storing data files

    //Spectrum proposal
    double *p; //!<power spectral density of data
    double pmax; //!<maximum power spectrial density
    double SNR2; //!<estimated \f${\rm SNR}^2\f$ of data
    ///@}

    /** @name Model Reconstructions */
     ///@{
    int Nwave; //!<Number of samples for computing posterior reconstructions
    int downsample; //!<Downsample factor for getting the desired number of samples

    double ****h_rec; //!<Store waveform reconstruction samples \f$ 2N \times N_\rm{channel} \times NT \times NMCMC \f$
    double ****h_res; //!<Store data residual samples \f$ 2N \times N_\rm{channel} \times NT \times NMCMC \f$
    double ****r_pow; //!<Store residual power samples \f$ N \times N_\rm{channel} \times NT \times NMCMC \f$
    double ****h_pow; //!<Store waveform power samples \f$ N \times N_\rm{channel} \times NT \times NMCMC \f$
    double ****S_pow; //!<Store noise power samples \f$ N \times N_\rm{channel} \times NT \times NMCMC \f$
    char fileName[PATH_BUFSIZE]; //!<place holder for filnames
    ///@}

    /** @name Signal Injections */
     ///@{
    int NP; //!<number of parameters of injection
    struct Source *inj; //!<injected source structure
    ///@}

    
    /** @name Already known sources */
    ///@{
    int Ncache; //!<number of sources in the cache file
    char **cache; //!<contents of cache file
    struct Catalog *catalog; //!< data and metadata for known sources
    ///@}
};

/*!
 *\brief Run flags set or changed from default values when parsing command line for `gb_mcmc` **and** `gb_catalog`.
 *
 *Descriptions of each flag includes default settings and command line arguments to change/set flags.
 *Boolean flags are defined as integers with `0==FALSE` and `1==TRUE`.
 *
 *To see how the defaults are set and then adjusted according to the command line, see parse().
 */
struct Flags
{
    /** @name Run Flags  */
     ///@{
    int verbose;    //!<`[--verbose; default=FALSE]`: increases file output (all chains, less downsampling, etc.)
    int quiet;      //!<`[--quiet; default=FALSE]`: decreases file output (less diagnostic data, only what is needed for post processing/recovery)
    int NMCMC;      //!<`[--steps=INT; default=100000]`: number of MCMC samples
    int NBURN;      //!<number of burn in samples, equals Flags::NMCMC.
    int NINJ;       //!<`[--inj=FILENAME]`: number of injections = number of `--inj` instances in command line
    int NDATA;      //!<`[default=1]`: number of frequency segments, equal to Flags::NINJ.
    int NT;         //!<`[--segments=INT; default=1]`: number of time segments
    int NVB;        //!<number of known binaries for `vb_mcmc`
    int DMAX;       //!<`[--sources=INT; default=10]`: max number of sources
    int simNoise;   //!<`[--sim-noise; default=FALSE]`: simulate random noise realization and add to data
    int fixSky;     //!<`[--fix-sky; default=FALSE]`: hold sky location fixed to injection parameters.  Set to `TRUE` if Flags::knownSource=`TRUE`.
    int fixFdot;
    int fixFreq;    //!<`[--fix-freq; default=FALSE]`: hold GW frequency fixed to injection parameters
    int galaxyPrior;//!<`[--galaxy-prior; default=FALSE]`: use model of galaxy for sky location prior
    int snrPrior;   //!<`[--snr-prior; default=FALSE]`: use SNR prior for amplitude
    int emPrior;    //!<`[--em-prior=FILENAME]`: use input data file with EM-derived parameters for priors.
    int knownSource;//!<`[--known-source; default=FALSE]`: injection is known binary, will need polarization and phase to be internally generated. Sets Flags::fixSky = `TRUE`.
    int detached;   //!<`[--detached; default=FALSE]`: assume binary is detached, fdot prior becomes \f$U[\dot{f}(\mathcal{M}_c=0.15),\dot{f}(\mathcal{M}_c=1.00)]\f$
    int strainData; //!<`[--data=FILENAME; default=FALSE]`: read data from ASCII file instead of simulate internally.
    int hdf5Data;   //!<'[--hdf5Data=FILENAME; default=FALSE]`: read data from LDC HDF5 file (compatible w/ Sangria dataset).
    int orbit;      //!<`[--orbit=FILENAME; default=FALSE]`: use numerical spacecraft ephemerides supplied in `FILENAME`. `--orbit` argument sets flag to `TRUE`.
    int prior;      //!<`[--prior; default=FALSE]`: set log-likelihood to constant for testing detailed balance.
    int debug;      //!<`[--debug; default=FALSE]`: coarser settings for proposals and verbose output for debugging
    int cheat;      //!<start sampler at injection values
    int burnin;     //!<`[--no-burnin; default=TRUE]`: chain is in the burn in phase
    int maximize;   //!<maximize over extrinsic parameter during burn in phase.
    int update;     //!<`[--update=FILENAME; default=FALSE]`: use Gaussian Mixture Model approximation to previous posterior as current prior.
    int updateCov;  //!<`[--update-cov=FILENAME; default=FALSE]`: updating fit from covariance matrix files built from chain samples, passed as `FILENAME`, used in draw_from_cov().
    int match;      //!<[--match=FLOAT; default=0.8]`: match threshold for chain sample clustering in post processing.
    int rj;         //!<--no-rj; default=TRUE]`: flag for determining if trans dimensional MCMC moves (RJMCMC) are enabled.
    int gap;        //!<`[--fit-gap; default=FALSE]`: flag for determining if model includes fit to time-gap in the data.
    int calibration;//!<`[--calibration; default=FALSE]`: flag for determining if model is marginalizing over calibration  uncertainty.
    int confNoise;  //!<`[--conf-noise; default=FALSE]`: include model of confusion noise in \f$S_n(f)\f$, either for simulating noise or as starting value for parameterized noise model.
    int resume;     //!<`[--resume; default=FALSE]`: restart sampler from run state saved during checkpointing. Starts from scratch if no checkpointing files are found.
    int catalog;    //!<`[--catalog=FILENAME; default=FALSE]`: use list of previously detected sources supplied in `FILENAME` to clean bandwidth padding (`gb_mcmc`) or for building family tree (`gb_catalog`).
    int threads;    //!<number of openMP threads for parallel tempering
    int psd;        //!<`[--psd=FILENAME; default=FALSE]`: use PSD input as ASCII file from command line
    ///@}

    
    /** @name Input File Names
     */
     ///@{
    char runDir[PATH_BUFSIZE];       //!<store `DIRECTORY` to serve as top level directory for output files.
    char vbFile[PATH_BUFSIZE];       //!<store `FILENAME` of list of known binaries `vb_mcmc`
    char **injFile;                   //!<`[--inj=FILENAME]`: list of injection files. Can support up to `NINJ=10` separate injections.
    char noiseFile[PATH_BUFSIZE];    //!<file containing reconstructed noise model for `gb_catalog` to compute SNRs against.
    char cdfFile[PATH_BUFSIZE];      //!<store `FILENAME` of input chain file from Flags::update.
    char gmmFile[PATH_BUFSIZE];      //!<store `FILENAME` of input gmm file from Flags::update.
    char covFile[PATH_BUFSIZE];      //!<store `FILENAME` of input covariance matrix file from Flags::updateCov.
    char matchInfile1[PATH_BUFSIZE]; //!<input waveform \f$A\f$ for computing match \f$(h_A|h_B)\f$
    char matchInfile2[PATH_BUFSIZE]; //!<input waveform \f$B\f$ for computing match \f$(h_A|h_B)\f$
    char pdfFile[PATH_BUFSIZE];      //!<store `FILENAME` of input priors for Flags:knownSource.
    char psdFile[PATH_BUFSIZE];      //!<store `FILENAME` of input psd file from Flags::psd.
    char catalogFile[PATH_BUFSIZE];  //!<store `FILENAME` containing previously identified detections from Flags::catalog for cleaning padding regions
     ///@}
};

/**
 \brief Structure containing settings and housekeeping data for each of the parallel chains.
 */
struct Chain
{
    /// Number of chains
    int NC;
    
    /// Number of proposals being used by chain
    int NP;

    /// Array containing current order of chains in temperature ladder
    int *index;
    
    /// Size of model being sampled by chain.  Depricated?
    int **dimension;
    
    /// Array tracking acceptance rate of parallel chain exchanges.
    double *acceptance;
    
    /// Array of chain temperatures
    double *temperature;
    
    /// Array storing \f$\langle \log p(d|\vec\theta)^{1/T} \rangle\f$ for thermodynamic integration.
    double *avgLogL;
    
    /// Annealing temperature for all chains used during burnin **DEPRICATED**
    double annealing;
    
    /// Store the maximum value of \f$\log p(d|\vec\theta)\f$ encountered by the chain.  Used for determining when burn-in is complete.
    double logLmax;
    
    /** @name Random Number Generator (RNG) Data Types
     Thread-safe random number generator data types used by `GSL`
     */
     ///@{
    /// Needed for initializing RNG
    const gsl_rng_type **T;
    
    /// Seed for RNGs for each parallel chain.
    gsl_rng **r;
    ///@}
    
    /** @name Chain File Pointers
     By default only the cold chain `M=0` is saved.  When Flags::verbose = `TRUE` files for each of the parallel chain are written.
     */
    ///@{
    
    /**
     \brief Noise parameter chain file: `chains/noise_chain.dat.M`
     
     Columns:  `step | logL | logL_norm | etaA | eta E`
     */
    FILE **noiseFile;
    
    /**
     \brief Markov chain state (iterations, likelihoods, etc.): `chains/model_chain.dat.M`
     
     Columns: `step | N_live | logL | logL_norm | t_0`
     */
    FILE **chainFile;
    
    /**
     \brief Calibration parameter chain file: `chains/calibration_chain.dat.M`
     
     Columns:
     */
    FILE **calibrationFile;
    
    /**
     
     \brief Signal parameter chain files for discrete models, only for cold chains: `chains/dimension_chain.dat.D`
     
     Columns: `f | fdot | A | cos_colat | long | cos_inc | psi | phi`
     */
    FILE **dimensionFile;
    
    /**
     \brief Full signal parameter chain files: `chains/dimension_chain.dat.M`

     Columns: `f | fdot | A | cos_colat | long | cos_inc | psi | phi`
     */
    FILE **parameterFile;
    
    /**
     \brief Log-likelhood values for each parallel chain.
     
     Columns: `step | logL_0 | logL_1 | ... | logL_NC-1`
     */
    FILE *likelihoodFile;
    
    /**
     \brief Temperature of each parallel chain to monitor adaptive temperature spacing.
     
     Columns: `step | 1/T_0 | 1/T_1 | ... | 1/T_NC-1`
     */
    FILE *temperatureFile;
    ///@}
    
    char chainDir[PATH_BUFSIZE]; //!<store chain directory.
    char chkptDir[PATH_BUFSIZE]; //!<store checkpoint directory.

};

/**
\brief Structure containing parameters and meta data for a single galactic binary.
*/
struct Source
{
    /// Number of parameters in signal model
    int NP;

    /// Array containing parameters to be passed to galactic_binary()
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
    ///See GalacticBinary.c for functions to convert between intrinsic and derived parameters
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
};

/**
\brief Structure containing parameters and meta data for noise model in narrow-band analysis segment.
*/
struct Noise
{
    /// Number of data samples fit by noise model
    int N;
    
    ///@name Constant Noise Parameters
    ///Each \f$\eta\f$ is a multiplier to the assumed noise level \f$S_n\f$ stored in Data structure. One per channel (`X` for 4-link, `A`,`E` for 6-link)
    ///@{
    double etaA;
    double etaE;
    double etaX;
    ///@}
    
    ///@name Noise Model
    ///Composite noise model to use over the analysis window \f$\eta_I \times Sn_I\f$
    ///@{
    double *f;
    double *SnA;
    double *SnE;
    double *SnX;
    double *transfer;
    ///@}
    
    ///@name UNDER CONSTRUCTION! noise parameters for power-law fit
    ///Parameters are noise levels at reference frequency \f$S_{n,0}\f$ and power law index \f$\alpha\f$.
    ///@{
    double SnA_0;
    double SnE_0;
    double SnX_0;
    double alpha_A;
    double alpha_E;
    double alpha_X;
    ///@}

};

/**
\brief Structure containing calibration parameters
 */
struct Calibration
{
    ///@name Amplitude parameters for each TDI channel
    ///@{
    double dampA;
    double dampE;
    double dampX;
    ///@}
    
    ///@name Overall phase parameters for each TDI channel
    ///@{
    double dphiA;
    double dphiE;
    double dphiX;
    ///@}

    ///@name Phase correction to Re and Im part of TDI channels
    ///@{
    double real_dphiA;
    double real_dphiE;
    double real_dphiX;
    double imag_dphiA;
    double imag_dphiE;
    double imag_dphiX;
    ///@}

};

/**
\brief Hierarchical structure of GBMCMC model
 */
struct Model
{
    ///@name Source parameters
    ///@{
    int NT;     //!<number of time segments
    int NP;     //!<maximum number of signal parameters
    int Nmax;   //!<maximum number of signals in model
    int Nlive;  //!<current number of signals in model
    struct Source **source; //!<source structures for each signal in the model
    ///@}
    
    /// Noise parameters
    struct Noise **noise;
    
    /// Calibration parameters
    struct Calibration **calibration;
    
    ///@name TDI structures
    ///@{
    struct TDI **tdi; //!<joint signal model
    struct TDI **residual; //!<joint residual
    ///@}
    
    ///@name Segment start time
    ///@{
    double *t0; //!<start time
    double *t0_min; //!<lower prior bound on start time
    double *t0_max; //!<upper prior bound on start time
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
};


