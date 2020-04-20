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

#ifndef LISA_h
#define LISA_h

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

/* Photon shot noise power */
#define SPS 8.321000e-23

/* Acceleration noise power */
#define SACC 9.000000e-30

/* Position noise? */
#define SLOC 2.89e-24

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
    
    gsl_spline **dx;
    gsl_spline **dy;
    gsl_spline **dz;
    gsl_interp_accel *acc;
    
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


void print_LISA_ASCII_art(FILE *fptr);

void interpolate_orbits(struct Orbit *orbit, double t, double *x, double *y, double *z);
void analytic_orbits(struct Orbit *orbit, double t, double *x, double *y, double *z);

void initialize_analytic_orbit(struct Orbit *orbit);
void initialize_numeric_orbit(struct Orbit *orbit);
void free_orbit(struct Orbit *orbit);

void LISA_tdi(double L, double fstar, double T, double ***d, double f0, long q, double *M, double *A, double *E, int BW, int NI);
double AEnoise(double L, double fstar, double f);
double XYZnoise(double L, double fstar, double f);
double GBnoise(double T, double f);

/* Fractional frequency versions of TDI & Sn(f) codes */
void LISA_tdi_FF(double L, double fstar, double T, double ***d, double f0, long q, double *M, double *A, double *E, int BW, int NI);
double AEnoise_FF(double L, double fstar, double f);
double GBnoise_FF(double T, double fstar, double f);


#endif /* LISA_h */
