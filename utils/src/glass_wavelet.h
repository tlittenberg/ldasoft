/*
 *  Copyright (C) 2024 Neil J. Cornish & Tyson B. Littenberg (MSFC-ST12)
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

#ifndef glass_wavelet_h
#define glass_wavelet_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_spline.h>

#include "glass_constants.h"

#define WAVELET_FILTER_CONSTANT 4  // filter steepness in frequency
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

struct Wavelets
{
    /** @name define time-frequency-fdot grid */
     ///@{
    int NF;    //!<frequency layers
    int NT;    //!<time layers
    double dt; //!<time cadence
    double df; //!<pixel bandwidth>
    
    int frequency_steps; //!<frequency steps
    int fdot_steps;      //!<fdot steps
    double d_fdot;       //!<fractional fdot increment
    ///@}

    int oversample;   //!<oversampling factor for basis functions
    
    double Omega;
    double dOmega;
    double domega;
    double inv_root_dOmega;
    double B;
    double A;
    double BW;
    double deltaf;
    double *fdot;
    double *window;
    double norm;
    
    int N;
    double T;

    gsl_vector **table; //!< lookup table of wavelet coefficients
};

void initialize_wavelet(struct Wavelets *wdm, int N, int NF, double dt);
void wavelet_transform(struct Wavelets *wdm, double *data);
void wavelet_index_to_pixel(struct Wavelets *wdm, int *i, int *j, int k);
void wavelet_pixel_to_index(struct Wavelets *wdm, int i, int j, int *k);
void wavelet_transfrom_from_table(struct Wavelets *wdm, double BW, double *Phase, double *freq, double *Amp, int *list, double *data);
void wavelet_lookup_table(struct Wavelets *wdm);

#endif /* glass_wavelet_h */
