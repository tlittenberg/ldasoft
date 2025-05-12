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

#include "glass_constants.h"

#define WAVELET_FILTER_CONSTANT 4  // filter steepness in frequency
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

struct Wavelets
{
    /** @name define time-frequency-fdot grid */
     ///@{
    int NF;    //!<total number of frequency layers
    int NT;    //!<total number of time layers
    double cadence; //!<data sample cadence
    double dt; //!<wavelet pixel duration
    double df; //!<wavelet pixel bandwidth
    
    int frequency_steps; //!<frequency steps for wavelet filters
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

    int kmin; //!<minimum wavelet index
    int kmax; //!<maximum wavelet index
    
    int *n_table; //!< number of terms in the lookup table at each frequency
    double **table; //!< lookup table of wavelet coefficients
};

void initialize_wavelet(struct Wavelets *wdm, double T);
void wavelet_window_frequency(struct Wavelets *wdm, double *window, int Nlayers);
void wavelet_transform(struct Wavelets *wdm, double *data);
void wavelet_transform_inverse(struct Wavelets *wdm, double *data);
void wavelet_to_fourier_transform(struct Wavelets *wdm, double *data);
void wavelet_index_to_pixel(struct Wavelets *wdm, int *i, int *j, int k);
void wavelet_pixel_to_index(struct Wavelets *wdm, int i, int j, int *k);
void wavelet_transform_by_layers(struct Wavelets *wdm, int jmin, int Nlayers, double *window, double *data);
void wavelet_transform_from_table(struct Wavelets *wdm, double *phase, double *freq, double *freqd, double *amp, int *jmin, int *jmax, double *wave, int *list, int *rlist, int Nmax);
void active_wavelet_list(struct Wavelets *wdm, double *freqX, double *freqY, double *freqZ, double *fdotX, double *fdotY, double *fdotZ, int *wavelet_list, int *reverse_list, int *Nwavelet, int *jmin, int *jmax);

#endif /* glass_wavelet_h */
