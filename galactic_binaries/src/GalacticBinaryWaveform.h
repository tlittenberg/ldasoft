//
//  GalacticBinaryWaveform.h
//  
//
//  Created by Littenberg, Tyson B. (MSFC-ZP12) on 1/15/17.
//
//

#ifndef GalacticBinaryWaveform_h
#define GalacticBinaryWaveform_h

#include <stdio.h>

double galactic_binary_Amp(double Mc, double f0, double D, double T);

double galactic_binary_fdot(double Mc, double f0, double T);

double galactic_binary_Mc(double f0, double dfdt, double T);

double galactic_binary_dL(double f0, double dfdt, double A, double T);

void galactic_binary_fisher(struct Orbit *orbit, struct Data *data, struct Source *source, struct Noise *noise);

void galactic_binary_alignment(struct Orbit *orbit, struct Data *data, struct Source *source);

int galactic_binary_bandwidth(double L, double fstar, double f, double A, double T, int N);

void galactic_binary(struct Orbit *orbit, double T, double t0, double params[], int NP, double *X, double *A, double *E, int BW, int NI);

#endif /* GalacticBinaryWaveform_h */
