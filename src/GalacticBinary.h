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
  int index;
  const gsl_rng_type *T;
  gsl_rng *r;
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

};

struct Noise
{
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
  int Nmax;
  struct Source **source;
  
  //Noise parameters
  struct Noise *noise;
  
  //TDI
  struct TDI *tdi;
  
};
