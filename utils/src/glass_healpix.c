/* -----------------------------------------------------------------------------
 *
 *  Copyright (C) 1997-2016 Krzysztof M. Gorski, Eric Hivon, Martin Reinecke,
 *                          Benjamin D. Wandelt, Anthony J. Banday,
 *                          Matthias Bartelmann,
 *                          Reza Ansari & Kenneth M. Ganga
 * 
 *  Copyright (C) 2024 Neil J. Cornish & Tyson B. Littenberg (MSFC-ST12) 
 *
 *  This file is modified from the original source code part of HEALPix.
 *
 *  HEALPix is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  HEALPix is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with HEALPix; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix see http://healpix.sourceforge.net
 *
 *---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void pix2ang_ring( long nside, long ipix, double *theta, double *phi) {
    /*
     c=======================================================================
     c     gives theta and phi corresponding to pixel ipix (RING)
     c     for a parameter nside
     c=======================================================================
     */
    
    int nl2, nl4, npix, ncap, iring, iphi, ip, ipix1;
    double  fact1, fact2, fodd, hip, fihip;
    double PI=M_PI;
    //      PARAMETER (pi     = 3.1415926535897932384626434d0)
    //      parameter (ns_max = 8192) ! 2^13 : largest nside available
    
    int ns_max=8192;
    
    if( nside<1 || nside>ns_max ) {
        fprintf(stderr, "%s (%d): nside out of range: %ld\n", __FILE__, __LINE__, nside);
        exit(0);
    }
    npix = 12*nside*nside;      // ! total number of points
    if( ipix<0 || ipix>npix-1 ) {
        fprintf(stderr, "%s (%d): ipix out of range: %ld\n", __FILE__, __LINE__, ipix);
        exit(0);
    }
    
    ipix1 = ipix + 1; // in {1, npix}
    nl2 = 2*nside;
    nl4 = 4*nside;
    ncap = 2*nside*(nside-1);// ! points in each polar cap, =0 for nside =1
    fact1 = 1.5*nside;
    fact2 = 3.0*nside*nside;
    
    if( ipix1 <= ncap ) {  //! North Polar cap -------------
        
        hip   = ipix1/2.;
        fihip = floor(hip);
        iring = (int)floor( sqrt( hip - sqrt(fihip) ) ) + 1;// ! counted from North pole
        iphi  = ipix1 - 2*iring*(iring - 1);
        
        *theta = acos( 1. - iring*iring / fact2 );
        *phi   = (1.*iphi - 0.5) * PI/(2.*iring);
    }
    else if( ipix1 <= nl2*(5*nside+1) ) {//then ! Equatorial region ------
        
        ip    = ipix1 - ncap - 1;
        iring = (int)floor( ip / nl4 ) + nside;// ! counted from North pole
        iphi  = (int)fmod(ip,nl4) + 1;
        
        fodd  = 0.5 * (1 + fmod((double)(iring+nside),2));//  ! 1 if iring+nside is odd, 1/2 otherwise
        *theta = acos( (nl2 - iring) / fact1 );
        *phi   = (1.*iphi - fodd) * PI /(2.*nside);
    }
    else {//! South Polar cap -----------------------------------
        
        ip    = npix - ipix1 + 1;
        hip   = ip/2.;
        /* bug corrige floor instead of 1.* */
        fihip = floor(hip);
        iring = (int)floor( sqrt( hip - sqrt(fihip) ) ) + 1;//     ! counted from South pole
        iphi  = (int)(4.*iring + 1 - (ip - 2.*iring*(iring-1)));
        
        *theta = acos( -1. + iring*iring / fact2 );
        *phi   = (1.*iphi - 0.5) * PI/(2.*iring);
    }
}
