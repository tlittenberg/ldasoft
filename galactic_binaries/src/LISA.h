//
//  LISA.h
//  
//
//  Created by Littenberg, Tyson B. (MSFC-ZP12) on 1/15/17.
//
//

#ifndef LISA_h
#define LISA_h

#include <stdio.h>

/* Photon shot noise power */
#define Sps 8.321000e-23

/* Acceleration noise power */
#define Sacc 9.000000e-30

/* Mean arm length of constellation (m) */
#define Larm 2.5e9


struct Orbit
{
  char OrbitFileName[1024];
  
  int Norb;
  
  double L;
  double fstar;
  double ecc;
  double R;
  
  double *t;
  double **x;
  double **y;
  double **z;
  double **dx;
  double **dy;
  double **dz;
  
  void (*orbit_function)(struct Orbit*,double,double*,double*,double*);
};

struct TDI
{
  //Michelson
  double *X;
  double *Y;
  double *Z;
  
  //Noise-orthogonal
  double *A;
  double *E;
  double *T;
  
  //Number of data channels
  int Nchannel;
  
  //Number of frequency bins
  int N;
};


void interpolate_orbits(struct Orbit *orbit, double t, double *x, double *y, double *z);
void analytic_orbits(struct Orbit *orbit, double t, double *x, double *y, double *z);

void initialize_analytic_orbit(struct Orbit *orbit);
void initialize_numeric_orbit(struct Orbit *orbit);
void free_orbit(struct Orbit *orbit);

void LISA_spline(double *x, double *y, int n, double yp1, double ypn, double *y2);
void LISA_splint(double *xa, double *ya, double *y2a, int n, double x, double *y);
void LISA_tdi(double L, double fstar, double T, double ***d, double f0, long q, double *M, double *A, double *E, int BW, int NI);
double AEnoise(double L, double fstar, double f);
double GBnoise(double T, double f);

/* Fractional frequency versions of TDI & Sn(f) codes */
void LISA_tdi_FF(double L, double fstar, double T, double ***d, double f0, long q, double *M, double *A, double *E, int BW, int NI);
double AEnoise_FF(double L, double fstar, double f);
double GBnoise_FF(double T, double fstar, double f);


#endif /* LISA_h */
