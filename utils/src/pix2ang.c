#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "astrometry/astrometry.h"
#include "pix2ang.h"

void astropy_pix2ang_ring(long nside, long ipix, double *theta, double *phi) {
    int64_t xy = healpixl_ring_to_xy(ipix,nside);
    healpixl_to_radec(xy,nside,0.5,0.5,phi,theta);
    *theta = M_PI/2.0 - *theta;
}
