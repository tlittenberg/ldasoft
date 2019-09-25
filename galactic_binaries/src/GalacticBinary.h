#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


struct Data
{
  int N;        //number of frequency bins
  int NT;       //number of time segments
  int Nchannel; //number of data channels
  int DMAX;     //max dimension of signal model

  long cseed; //seed for MCMC
  long nseed; //seed for noise realization
  long iseed; //seed for injection parameters
  
  double T;
  
  int qmin;
  int qmax;
  
  double fmin;
  double fmax;
  
  double *t0;   //start times of segments
  double *tgap; //time between segments
  
  //Response
  struct TDI **tdi;
  struct Noise **noise;
  
  //Reconstructed signal
  int Nwave;
  int downsample;
  double ****h_rec; // N x Nchannel x NT x NMCMC
  double ****h_res; // N x Nchannel x NT x NMCMC
  double ****r_pow; // N x Nchannel x NT x NMCMC
  double ****h_pow; // N x Nchannel x NT x NMCMC
  double ****S_pow; // N x Nchannel x NT x NMCMC
  
  //Injection
  int NP; //number of parameters of injection
  struct Source *inj;
  
  //Spectrum proposal
  double *p;
  double pmax;
  double SNR2;
  
  //
  char fileName[128];
  
  /*
   Data format string 
   'phase'     ==> LISA Simulator esque
   'frequency' ==> Synthetic LISA esque
  */
  char format[16];
  
};

struct Flags
{
  int verbose;
  int NMCMC; //number of MCMC steps
  int NBURN; //number of Burn-in steps
  int NINJ; //number of frequency segments;
  int NDATA;  //number of frequency segments;
  int NT;    //number of time segments
  int DMAX;  //max number of sources
  int simNoise;
  int fixSky;
  int fixFreq;
  int galaxyPrior;
  int snrPrior;
  int emPrior;
  int knownSource;
  int detached;
  int strainData;
  int orbit;
  int prior;
  int debug;
  int cheat;
  int burnin;
  int update;
  int updateCov;
  int match;
  int rj;
  int gap; //are we fitting for a time-gap in the data?
  int calibration; //are we marginalizing over calibration  uncertainty?
  int confNoise; //include model of confusion noise in Sn(f)
  
  char **injFile;
  char cdfFile[128];
  char covFile[128];
  char matchInfile1[128];
  char matchInfile2[128];
  char pdfFile[128];
};

struct Chain
{
  //Number of chains
  int NC;
  int NP;
  int NS;
  int *index;
  int **dimension;
  double *acceptance;
  double *temperature;
  double *avgLogL;
  double annealing;
  double logLmax;
  
  //thread-safe RNG
  const gsl_rng_type **T;
  gsl_rng **r;
  
  //chain files
  FILE **noiseFile;
  FILE **chainFile;
  FILE **calibrationFile;
  FILE **dimensionFile;
  FILE **parameterFile;
  FILE *likelihoodFile;
  FILE *temperatureFile;
};

struct Source
{
  //Intrinsic
  double m1;
  double m2;
  double f0;

  //Extrinisic
  double psi;
  double cosi;
  double phi0;

  double D;
  double phi;
  double costheta;

  //Derived
  double amp;
  double dfdt;
  double d2fdt2;
  double Mc;
  
  //Book-keeping
  int BW;
  int qmin;
  int qmax;
  int imin;
  int imax;
  
  //Response
  struct TDI *tdi;
  
  //Fisher matrix
  double **fisher_matrix;
  double **fisher_evectr;
  double *fisher_evalue;

  //Package parameters for waveform generator
  int NP;
  double *params;

};

struct Noise
{
  int N;
  
  double etaA;
  double etaE;
  double etaX;
  
  double *SnA;
  double *SnE;
  double *SnX;
};

struct Calibration
{
  double dampA;
  double dampE;
  double dampX;
  double dphiA;
  double dphiE;
  double dphiX;
  double real_dphiA;
  double real_dphiE;
  double real_dphiX;
  double imag_dphiA;
  double imag_dphiE;
  double imag_dphiX;
};

struct Model
{
  //Source parameters
  int NT;     //number of time segments
  int NP;     //maximum number of signal parameters
  int Nmax;   //maximum number of signals in model
  int Nlive;  //current number of signals in model
  struct Source **source;
  
  //Noise parameters
  struct Noise **noise;
  
  //Calibration parameters
  struct Calibration **calibration;
  
  //TDI
  struct TDI **tdi;
  struct TDI **residual;
  
  //Start time for segment for model
  double *t0;
  double *t0_min;
  double *t0_max;
  
  //Source parameter priors
  double **prior;
  double *logPriorVolume;
  
  //Model likelihood
  double logL;
  double logLnorm;
};


