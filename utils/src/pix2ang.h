// This function is a wrapper around the astrometry.net functions to match healpix's conventions.
// from Robert Rosati, 2025, ST-12 / UAH
// based on astrometry.net and astropy-healpix
// https://github.com/astropy/astropy-healpix/tree/main

#ifndef _PIX2ANG_H
#define _PIX2ANG_H
void astropy_pix2ang_ring(long nside, long ipix, double *theta, double *phi);
#endif
