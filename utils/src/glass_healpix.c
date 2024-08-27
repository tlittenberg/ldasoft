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

void ang2pix_ring( const long nside, double theta, double phi, long *ipix) {
    /*
     c=======================================================================
     c     gives the pixel number ipix (RING)
     c     corresponding to angles theta and phi
     c=======================================================================
     */
    
    int nl2, nl4, ncap, npix, jp, jm, ipix1;
    double  z, za, tt, tp, tmp;
    int ir, ip, kshift;
    
    double piover2 = 0.5*M_PI;
    double PI=M_PI;
    double twopi=2.0*M_PI;
    double z0=2.0/3.0;
    long ns_max=8192;
    
    if( nside<1 || nside>ns_max ) {
        fprintf(stderr, "%s (%d): nside out of range: %ld\n", __FILE__, __LINE__, nside);
        exit(0);
    }
    
    if( theta<0. || theta>PI) {
        fprintf(stderr, "%s (%d): theta out of range: %f\n", __FILE__, __LINE__, theta);
        exit(0);
    }
    
    z = cos(theta);
    za = fabs(z);
    if( phi >= twopi)  phi = phi - twopi;
    if (phi < 0.)     phi = phi + twopi;
    tt = phi / piover2;//  ! in [0,4)
    
    nl2 = 2*nside;
    nl4 = 4*nside;
    ncap  = nl2*(nside-1);// ! number of pixels in the north polar cap
    npix  = 12*nside*nside;
    
    if( za <= z0 ) {
        
        jp = (int)floor(nside*(0.5 + tt - z*0.75)); /*index of ascending edge line*/
        jm = (int)floor(nside*(0.5 + tt + z*0.75)); /*index of descending edge line*/
        
        ir = nside + 1 + jp - jm;// ! in {1,2n+1} (ring number counted from z=2/3)
        kshift = 0;
        if (fmod(ir,2)==0.) kshift = 1;// ! kshift=1 if ir even, 0 otherwise
        
        ip = (int)floor( ( jp+jm - nside + kshift + 1 ) / 2 ) + 1;// ! in {1,4n}
        if( ip>nl4 ) ip = ip - nl4;
        
        ipix1 = ncap + nl4*(ir-1) + ip ;
    }
    else {
        
        tp = tt - floor(tt);//      !MOD(tt,1.d0)
        tmp = sqrt( 3.*(1. - za) );
        
        jp = (int)floor( nside * tp * tmp );// ! increasing edge line index
        jm = (int)floor( nside * (1. - tp) * tmp );// ! decreasing edge line index
        
        ir = jp + jm + 1;//        ! ring number counted from the closest pole
        ip = (int)floor( tt * ir ) + 1;// ! in {1,4*ir}
        if( ip>4*ir ) ip = ip - 4*ir;
        
        ipix1 = 2*ir*(ir-1) + ip;
        if( z<=0. ) {
            ipix1 = npix - 2*ir*(ir+1) + ip;
        }
    }
    *ipix = ipix1 - 1;// ! in {0, npix-1}
    
}

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