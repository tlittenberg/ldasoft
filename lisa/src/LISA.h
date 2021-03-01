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
@file LISA.h
\brief Codes defining LISA orbits, noise assumptions, and TDI method.
*/

#ifndef LISA_h
#define LISA_h

#include <stdio.h>
#include <hdf5.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "Constants.h"

/// Mean arm length of constellation (m) for baseline LISA configuration
#define Larm 2.5e9

/** @name Component Noise Levels For Phase Data */
///@{

/// Photon shot noise power
#define SPS 8.321000e-23

/// Acceleration noise power
#define SACC 9.000000e-30

/// Position noise power when using phase data
#define SLOC 2.89e-24
///@}


/**
 * \brief Ephemerides of individual spacecraft and metadata for using orbits in waveform modeling.
 *
 * If numerical orbit files are provided, they are interpolated to the sample rate of the data using `GSL` cubic spline functions.
 *
 * If not, the eccentric inclined analytic model is computed once at the data sampling rate and stored.
 */
struct Orbit
{
    /// Filename input from `--orbit` command line argument when using numerical orbits
    char OrbitFileName[1024];
    
    /// Size of orbit arrays
    int Norb;
    
    /** @name Constellation Parameters */
    ///@{
    double L;     //!< Average armlength of constellation
    double fstar; //!< Transfer frequency \f$f_* = 1/(L/c)\f$.
    double ecc;   //!< Eccentricity of spacecraft orbits
    double R;     //!< Distance to constellation guiding center from Sun (1 AU)
    ///@}
    
    double *t;  //!<time step relative to start of mission (seconds)

    /** @name Spacecraft Ephemerides
     Cartesian ecliptic coordinates, in meters, where the x-y plane is the ecliptic.
     */
    ///@{
    double **x;
    double **y;
    double **z;
    ///@}
    
    /** @name Derivatives of Orbits for Cubic Spline Interpolation
     Stores the derivatives for the cubic spline for the location of each spacecraft (\f$dx,dy,dz\f$) and some internal `GSL` workspace `acc`.
     */
    ///@{
    gsl_spline **dx;
    gsl_spline **dy;
    gsl_spline **dz;
    gsl_interp_accel *acc;
    ///@}
    
    
    /**
     \brief Function prototyp for retreiving spacecraft locations.
     \param[in] time (double)
     \param[in] orbit data (Orbit*)
     \param[out] eccliptic carteisian location of each spacecraft
     
     If an orbit file is input this points to interpolate_orbits() which uses `GSL` cubic splines to interpolate the ephemerides at the needed time steps
     
     Otherwise this points to analytic_orbits() which is passed an arbitrary time \f$t\f$ and returns the spacecraft location.
     */
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
    
    //Number of data bins
    int N;
    
    //Data cadence
    double delta;
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
double Tnoise_FF(double L, double fstar, double f);
double GBnoise_FF(double T, double fstar, double f);
double XYZnoise_FF(double L, double fstar, double f);

void test_noise_model(struct Orbit *orbit);

void alloc_tdi(struct TDI *tdi, int NFFT, int Nchannel);
void copy_tdi(struct TDI *origin, struct TDI *copy);
void free_tdi(struct TDI *tdi);

/* LDC HDF5 */
void LISA_Read_HDF5_LDC_TDI(struct TDI *tdi, char *fileName);

#endif /* LISA_h */
