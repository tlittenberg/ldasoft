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

#include <sys/stat.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_statistics_double.h>

#include "glass_constants.h"
#include "glass_lisa.h"
#include "glass_wavelet.h"
#include "glass_data.h"
#include "glass_math.h"
#include "glass_gmm.h"
#include "glass_healpix.h"
#include "glass_galaxy.h"

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
