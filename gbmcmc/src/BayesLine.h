#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>

#include <gsl/gsl_sort.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_statistics.h>

//#define SAmin  1.0e-48
//#define SAmax  1.0e-30
//#define LQmin  1.0e2  // allowing Q's below ~ 1e2 starts a fight with the spline model as lines get too fat
//#define LQmax 1.0e8
//#define LAmin  1.0e-44
//#define LAmax  1.0e-30
#define lAwidth 2.3    // 2.3 is one decade
#define zeta 1.0
#define kappa_BL 0.8
#define FSTEP 10.0   // the stencil separation in Hz for the spline model. Going below 2 Hz is dangerous - will fight with line model
// as we are getting down to the line width of the broadest lines


typedef struct
{
  int n;
  int size;

  int *larray;

  double *Q;
  double *A;
  double *f;

}lorentzianParams;

typedef struct
{
  int tmax;
  int ncut;
  int nmin;
  int tfull;
  int sgmts;

  double df;
  double fny;
  double Tobs;
  double fmin;
  double fmax;
  double flow;
  double fgrid;
  double fstep;
  double fhigh;
  double cadence;

}dataParams;

typedef struct
{
  int n;
  double *points;
  double *data;

}splineParams;

typedef struct
{
  double SAmin;
  double SAmax;
  double LQmin;
  double LQmax;
  double LAmin;
  double LAmax;

  //double *invsigma; //variances for each frequency bin
  double *sigma; //variances for each frequency bin
  double *upper; //variances for each frequency bin
  double *lower; //variances for each frequency bin
  double *mean;     //means for each frequency bin

}BayesLinePriors;

struct BayesLineParams
{
  dataParams *data;
  splineParams *spline;
  splineParams *spline_x;
  lorentzianParams *lines_x;
  lorentzianParams *lines_full;
  BayesLinePriors *priors;

  double *Snf;
  double *Sna;
  double *fa;
  double *freq;
  double *power;
  double *spow;
  double *sfreq;
  double *Sbase;
  double *Sline;

  int constantLogLFlag;

  double TwoDeltaT;
  gsl_rng *r;
  
  FILE *splineChainFile;
  FILE *lineChainFile;
};

void BayesLineFree(struct BayesLineParams *bptr);
void BayesLineSetup(struct BayesLineParams *bptr, double *freqData, double fmin, double fmax, double deltaT, double Tobs);

void BayesLineRJMCMC(struct BayesLineParams *bayesline, double *freqData, double *psd, double *invpsd, double *splinePSD, int N, int cycle, double beta, int priorFlag);
void BayesLineSearch(struct BayesLineParams *bptr, double *freqData, double fmin, double fmax, double deltaT, double Tobs);

void BayesLineNonMarkovianFit          (struct BayesLineParams *bayesline, int *nj);
void BayesLineLorentzSplineMCMC        (struct BayesLineParams *bayesline, double heat, int steps, int focus, int priorFlag, double *dan);
void BayesLineMarkovianSplineOnly      (struct BayesLineParams *bayesline, int nspline, int jj);
void BayesLineMarkovianFocusedAnalysis (struct BayesLineParams *bayesline);

double loglike_fit_spline(double *respow, double *Snf, int ncut);

double loglike_pm        (double *respow, double *Sn, double *Snx, int ilow, int ihigh);
double loglike_single    (double *respow, double *Sn, double *Snx, int ilowx, int ihighx, int ilowy, int ihighy);

double sample(double *fprop, double pmax, dataParams *data, gsl_rng *r);
double lprop(double f, double *fprop, dataParams *data);

void full_spectrum_single(double *Sn, double *Snx, double *Sbasex, double *sfreq, dataParams *data, lorentzianParams *line_x, lorentzianParams *line_y, int ii, int *ilowx, int *ihighx, int *ilowy, int *ihighy);
void full_spectrum_add_or_subtract(double *Snew, double *Sold, double *Sbase, double *sfreq, dataParams *data, lorentzianParams *lines, int ii, int *ilow, int *ihigh, int flag);
void full_spectrum_spline(double *Sline, double *Sbase, double *sfreq, dataParams *data, lorentzianParams *lines);

void spectrum_spline(double *Sn, double *Sbase, double *sfreq, dataParams *data, lorentzianParams *lines, splineParams *spline);

void SpecFitSpline    (BayesLinePriors *priors, int zeroLogL, int steps, double *freq, double *power, splineParams *spline, double *Snf, int ncut, gsl_rng *r);
void LorentzSplineFit (BayesLinePriors *priors, int zeroLogL, int steps, dataParams *data, lorentzianParams *lines, splineParams *spline, double *sfreq, double *spow, gsl_rng *r);

void CubicSplineGSL(int N, double *x, double *y, int Nint, double *xint, double *yint);

void create_dataParams(dataParams *data, double *f, int n);

void create_lorentzianParams(lorentzianParams *lines, int size);
void copy_lorentzianParams(lorentzianParams *origin, lorentzianParams *copy);
void destroy_lorentzianParams(lorentzianParams *lines);

void create_splineParams(splineParams *spline, int size);
void copy_splineParams(splineParams *origin, splineParams *copy);
void destroy_splineParams(splineParams *spline);

void copy_bayesline_params(struct BayesLineParams *origin, struct BayesLineParams *copy);
void print_line_model(FILE *fptr, struct BayesLineParams *bayesline);
void print_spline_model(FILE *fptr, struct BayesLineParams *bayesline);
void parse_line_model(FILE *fptr, struct BayesLineParams *bayesline);
void parse_spline_model(FILE *fptr, struct BayesLineParams *bayesline);



