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

#ifndef GalacticBinaryWaveform_h
#define GalacticBinaryWaveform_h

#include <stdio.h>

double galactic_binary_Amp(double Mc, double f0, double D, double T);

double galactic_binary_fdot(double Mc, double f0, double T);

double galactic_binary_Mc(double f0, double dfdt, double T);

double galactic_binary_dL(double f0, double dfdt, double A);

void galactic_binary_fisher(struct Orbit *orbit, struct Data *data, struct Source *source, struct Noise *noise);

void galactic_binary_alignment(struct Orbit *orbit, struct Data *data, struct Source *source);

int galactic_binary_bandwidth(double L, double fstar, double f, double fdot, double costheta, double A, double T, int N);

void galactic_binary(struct Orbit *orbit, char *format, double T, double t0, double params[], int NP, double *X, double *A, double *E, int BW, int NI);

#endif /* GalacticBinaryWaveform_h */
