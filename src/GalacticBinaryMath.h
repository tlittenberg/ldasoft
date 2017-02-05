//
//  GalacticBinaryMath.h
//  
//
//  Created by Littenberg, Tyson B. (MSFC-ZP12) on 1/15/17.
//
//

#ifndef GalacticBinaryMath_h
#define GalacticBinaryMath_h

#include <stdio.h>

double swap, tempr;
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

double chirpmass(double m1, double m2);

double ipow(double x, int n);

/* ********************************************************************************** */
/*																					  */
/*                                    Fourier Tools                                   */
/*																					  */
/* ********************************************************************************** */

void dfour1(double data[], unsigned long nn, int isign);
void drealft(double data[], unsigned long n, int isign);

double power_spectrum(double *data, int n);

#endif /* GalacticBinaryMath_h */
