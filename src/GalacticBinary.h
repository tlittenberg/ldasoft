#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


struct Data
{
  int N;        //number of frequency bins
  int Nchannel; //number of data channels
  
  long cseed; //seed for MCMC
  long nseed; //seed for noise realization
  long iseed; //seed for injection parameters
  
  double T;
  
  int qmin;
  int qmax;
  
  double fmin;
  double fmax;
  
  double t0;
  
  //Response
  struct TDI *tdi;
  struct Noise *noise;
  
  //Reconstructed signal
  int Nwave;
  int downsample;
  double ***h_rec; // N x Nchannel x NMCMC 
  double ***h_res; // N x Nchannel x NMCMC
  double ***h_pow; // N x Nchannel x NMCMC
  double ***S_pow; // N x Nchannel x NMCMC
  
  char injFile[1024];
};

struct Flags
{
  int verbose;
  int injection;
  int zeroNoise;
};

struct Chain
{
  //Number of chains
  int NC;
  int *index;
  double *acceptance;
  double *temperature;
  double *avgLogL;
  
  //thread-safe RNG
  const gsl_rng_type **T;
  gsl_rng **r;
  
  //chain files
  FILE **noiseFile;
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
  double Mc;
  
  //Package parameters for waveform generator
  double *params;
  
  //Response
  struct TDI *tdi;

  //Book-keeping
  int BW;
  int qmin;
  int qmax;
  int imin;
  int imax;
  
  double t0;
  
  //Fisher matrix
  double **fisher_matrix;
  double **fisher_evectr;
  double *fisher_evalue;

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

struct Model
{
  //Source parameters
  int Nmax;   //maximum number of signals in model
  int Nlive;  //current number of signals in model
  struct Source **source;
  
  //Noise parameters
  struct Noise *noise;
  
  //TDI
  struct TDI *tdi;
  
  //Source parameter priors
  double **prior;
  double logPriorVolume;
  
  //Model likelihood
  double logL;
  double logLnorm;
};
