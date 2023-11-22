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
@file data.h
\brief Definition of GLASS data structures.
*/

#ifndef data_h
#define data_h

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define MAXSTRINGSIZE 1024 //!<maximum number of characters for `path+filename` strings

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
    double t0;   //!<start times of segments
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
    struct TDI *tdi; //!<TDI data channels as seen by sampler
    struct TDI *raw; //!<TDI data channels unaltered from input
    struct Noise *noise; //!<Reference noise model
    /**
     \brief Convention for data format
     
     format = "phase" for phase difference (distance)
     format = "frequency" for fractional frequency (velocity) **Use for matching LDC Radler**
     format = "sangria" for fractional frequency w/ LDC Sangria-era TDI & phase conventions
     */
    char format[16];
    char dataDir[MAXSTRINGSIZE]; //!<Directory for storing data files

    //Spectrum proposal
    double *p; //!<power spectral density of data
    double pmax; //!<maximum power spectrial density
    double SNR2; //!<estimated \f${\rm SNR}^2\f$ of data
    ///@}

    /** @name Model Reconstructions */
     ///@{
    int Nwave; //!<Number of samples for computing posterior reconstructions
    int downsample; //!<Downsample factor for getting the desired number of samples

    double ***h_rec; //!<Store waveform reconstruction samples \f$ 2N \times N_\rm{channel} \times NMCMC \f$
    double ***h_res; //!<Store data residual samples \f$ 2N \times N_\rm{channel} \times NMCMC \f$
    double ***r_pow; //!<Store residual power samples \f$ N \times N_\rm{channel} \times NMCMC \f$
    double ***h_pow; //!<Store waveform power samples \f$ N \times N_\rm{channel} \times NMCMC \f$
    double ***S_pow; //!<Store noise power samples \f$ N \times N_\rm{channel} \times NMCMC \f$
    char fileName[MAXSTRINGSIZE]; //!<place holder for filnames
    ///@}

    /** @name Signal Injections */
     ///@{
    //int NP; //!<number of parameters of injection
    //struct Source *inj; //!<injected source structure
    ///@}

    
    /** @name Already known sources */
    ///@{
    //int Ncache; //!<number of sources in the cache file
    //char **cache; //!<contents of cache file
    //struct Catalog *catalog; //!< data and metadata for known sources
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
    char runDir[MAXSTRINGSIZE];       //!<store `DIRECTORY` to serve as top level directory for output files.
    char vbFile[MAXSTRINGSIZE];       //!<store `FILENAME` of list of known binaries `vb_mcmc`
    char **injFile;                   //!<`[--inj=FILENAME]`: list of injection files. Can support up to `NINJ=10` separate injections.
    char noiseFile[MAXSTRINGSIZE];    //!<file containing reconstructed noise model for `gb_catalog` to compute SNRs against.
    char cdfFile[MAXSTRINGSIZE];      //!<store `FILENAME` of input chain file from Flags::update.
    char gmmFile[MAXSTRINGSIZE];      //!<store `FILENAME` of input gmm file from Flags::update.
    char covFile[MAXSTRINGSIZE];      //!<store `FILENAME` of input covariance matrix file from Flags::updateCov.
    char matchInfile1[MAXSTRINGSIZE]; //!<input waveform \f$A\f$ for computing match \f$(h_A|h_B)\f$
    char matchInfile2[MAXSTRINGSIZE]; //!<input waveform \f$B\f$ for computing match \f$(h_A|h_B)\f$
    char pdfFile[MAXSTRINGSIZE];      //!<store `FILENAME` of input priors for Flags:knownSource.
    char psdFile[MAXSTRINGSIZE];      //!<store `FILENAME` of input psd file from Flags::psd.
    char catalogFile[MAXSTRINGSIZE];  //!<store `FILENAME` containing previously identified detections from Flags::catalog for cleaning padding regions
     ///@}
};

/**
 \brief Structure containing settings and housekeeping data for each of the parallel chains.
 */
struct Chain
{
    /// Number of chains
    int NC;
    
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
    FILE **foregroundFile;
    
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
    
    char chainDir[MAXSTRINGSIZE]; //!<store chain directory.
    char chkptDir[MAXSTRINGSIZE]; //!<store checkpoint directory.

};



/**
\brief Structure containing parameters and meta data for noise model in narrow-band analysis segment.
*/
struct Noise
{
    /// Number of data samples fit by noise model
    int N;
    
    /// Number of TDI channels
    int Nchannel;
    
    ///@name Constant Noise Parameters
    ///Each \f$\eta\f$ is a multiplier to the assumed noise level \f$S_n\f$ stored in Data structure. One per channel (`X` for 4-link, `A`,`E` or `X`,`Y`,`Z` for 6-link)
    ///@{
    double *eta;
    ///@}
    
    ///@name Noise Model
    ///Composite noise model to use over the analysis window \f$\eta_I \times Sn_I\f$
    ///@{
    double *f;
    
    double ***C;    //!<Covariance matrix
    double ***invC; //!<Inverse covariance matrix>
    double *detC;   //!<Determinent of covariance matrix
    
    double *transfer;
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
    double dampY;
    double dampZ;
    ///@}
    
    ///@name Overall phase parameters for each TDI channel
    ///@{
    double dphiA;
    double dphiE;
    double dphiX;
    double dphiY;
    double dphiZ;
    ///@}

    ///@name Phase correction to Re and Im part of TDI channels
    ///@{
    double real_dphiA;
    double real_dphiE;
    double real_dphiX;
    double real_dphiY;
    double real_dphiZ;
    double imag_dphiA;
    double imag_dphiE;
    double imag_dphiX;
    double imag_dphiY;
    double imag_dphiZ;
    ///@}

};

/**
 \brief Show progress bar
 */
void printProgress (double percentage);

/**
 \brief Print git hash of running version
 */
void print_version(FILE *fptr);

/**
  \brief Name and create output directories
 */
void setup_run_directories(struct Flags *flags, struct Data *data, struct Chain *chain);

/**
 \brief Loads spacecraft orbits, either analytically or from tabulated ephemerides.
 */
void initialize_orbit(struct Data *data, struct Orbit *orbit, struct Flags *flags);

/**
 \brief Allocates and initializes Chain structure and prepares output files.
 */
void initialize_chain(struct Chain *chain, struct Flags *flags, long *seed, const char *mode);

/** @name Allocate memory for structures */
///@{
void alloc_data(struct Data *data, struct Flags *flags);
void alloc_noise(struct Noise *noise, int NFFT, int Nchannel);
void alloc_calibration(struct Calibration *calibration);
///@}

/**
 \brief Shallow copy of Data structure
 */
void copy_data(struct Data *origin, struct Data *copy);

/** @name Deep copy structure contents */
///@{
void copy_noise(struct Noise *origin, struct Noise *copy);
void copy_Cij(double ***origin, double ***copy, int M, int N);
void copy_calibration(struct Calibration *origin, struct Calibration *copy);
///@}

/** @name Free memory for structures */
///@{
void free_noise(struct Noise *noise);
void free_chain(struct Chain *chain, struct Flags *flags);
void free_calibration(struct Calibration *calibration);
///@}

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
 \brief Reads data from external source input using `--data` flag
 */
void ReadData(struct Data *data, struct Orbit *orbit, struct Flags *flags);

/**
 \brief Reads LDC-formatted HDF5 data using `--h5-data` flag
 */
void ReadHDF5(struct Data *data, struct TDI *tdi);

/**
 \brief Reads ASCII data using `--data` flag
 */
void ReadASCII(struct Data *data, struct TDI *tdi);

/**
 \brief Get theoretical noise PSDs for TDI channels
 */
void GetNoiseModel(struct Data *data, struct Orbit *orbit, struct Flags *flags);

/**
 \brief Add simulated Gaussian noise realization to data
 */
void AddNoise(struct Data *data, struct TDI *tdi);

/**
 \brief Generate noise-only simulated data
 */
void SimulateData(struct Data *data, struct Orbit *orbit, struct Flags *flags);

/**
 \brief Wrapper function that calls data print functions
 */
void print_data(struct Data *data, struct TDI *tdi, struct Flags *flags);

/**
 \brief Print list of command line arguments and defaults
 */
void print_usage();

/**
 \brief Parse command line
 */
void parse(int argc, char **argv, struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Chain *chain, int Nmax);

/**
 \brief Check if file exists and handle missing files cleanly
 */
int checkfile(char filename[]);

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

#endif /* data_h */
