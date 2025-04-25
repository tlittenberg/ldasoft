//
//  glass_utils.h
//  glass
//
//  Created by Tyson Littenberg on 11/29/23.
//

/** 
 @file glass_utils.h
 \brief Libary for GLASS package
 
 Including
 - GSL dependencies
 - physical constants
 - LISA constellation 
 - common math functions
 */

#ifndef utils_h
#define utils_h

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <time.h>
#include <hdf5.h>
#include <omp.h>
#include <lapacke.h>

#include <sys/stat.h>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_complex.h>

#include "glass_constants.h"
#include "glass_lisa.h"
#include "glass_wavelet.h"
#include "glass_data.h"
#include "glass_math.h"
#include "glass_gmm.h"
#include "glass_healpix.h"
#include "glass_galaxy.h"

/**
\brief wrapper for rand\_r() -- threadsafe RNG for U[0,1]
*/
double rand_r_U_0_1(unsigned int *seed);

/**
\brief Threadsafe RNG for N[0,1] using Marsaglia polar method
*/
double rand_r_N_0_1(unsigned int *seed);

int *int_vector(int N);
void free_int_vector(int *v);

int **int_matrix(int N, int M);
void free_int_matrix(int **m, int N);

double *double_vector(int N);
void free_double_vector(double *v);

double **double_matrix(int N, int M);
void free_double_matrix(double **m, int N);

double ***double_tensor(int N, int M, int L);
void free_double_tensor(double ***t, int N, int M);

#endif /* utils_h */
