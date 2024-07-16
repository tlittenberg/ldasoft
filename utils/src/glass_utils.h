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
#include "glass_data.h"
#include "glass_math.h"
#include "glass_gmm.h"

#endif /* utils_h */
