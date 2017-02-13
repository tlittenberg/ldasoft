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
#define Sps 1.21e-22

/* Acceleration noise power */
#define Sacc 21.00e-30

struct Orbit
{
  char OrbitFileName[1024];
  
  int Norb;
  
  double L;
  double fstar;
  
  double *t;
  double **x;
  double **y;
  double **z;
  double **dx;
  double **dy;
  double **dz;
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


void spacecraft(struct Orbit *orbit, double tint, double *xint, double *yint, double *zint);

void initialize_orbit(struct Orbit *orbit);
void free_orbit(struct Orbit *orbit);

void LISA_spline(double *x, double *y, int n, double yp1, double ypn, double *y2);
void LISA_splint(double *xa, double *ya, double *y2a, int n, double x, double *y);
void LISA_tdi(double L, double fstar, double T, double ***d, double f0, long q, double *M, double *A, double *E, int BW, int NI);
double AEnoise(double L, double fstar, double f);

#endif /* LISA_h */
