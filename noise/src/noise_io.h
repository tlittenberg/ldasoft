/*
 *  Copyright (C) 2023 Tyson B. Littenberg (MSFC-ST12)
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

#ifndef noise_io_h
#define noise_io_h

#include <stdio.h>
#include "noise.h"



/**
 \brief Print current state of spline model to ASCII
 */
void print_spline_state(struct SplineModel *model, FILE *fptr, int step);

/**
 \brief Print current state of instrument model to ASCII
 */
void print_instrument_state(struct InstrumentModel *model, FILE *fptr);

/**
 \brief Print current state of instrument model to ASCII
 */
void print_foreground_state(struct ForegroundModel *model, FILE *fptr);

/**
 \brief Print full PSD model to file named `filename`
 */
void print_noise_model(struct Noise *noise, char filename[]);

/**
 \brief Print data whitened by modeled variance to file named `filename`
 */
void print_whitened_data(struct Data *data, struct Noise *noise, char filename[]);

/**
 \brief Print PSD 5, 25, 50, 75, and 95 quantiles of noise model posteriors
 */
void print_noise_reconstruction(struct Data *data, struct Flags *flags);

#endif /* noise_io_h */
