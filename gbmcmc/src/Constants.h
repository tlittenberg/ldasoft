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

/**
 @file Constants.h
 \brief Mathematical and physical constants
 */


#ifndef Constants_h
#define Constants_h

#define FIXME 0
#define SNRPEAK 10 ///< Peak of SNR prior

/* --------------  MATHEMATICAL CONSTANTS  -------------- */
/* Some square roots */
#define SQ2 1.4142135623731  ///< \f$\sqrt{2}\f$
#define SQ3 1.73205080757    ///< \f$\sqrt{3}\f$
#define SQ8 2.82842712474619 ///< \f$\sqrt{8}\f$

/* Pi's and frinds */
//use math.h (M_PI) for PI
#define PI2   6.283185307179586      ///< \f$2\pi\f$
#define PIon2 1.57079632679          ///< \f$\pi/2\f$
#define PIon4 0.78539816339          ///< \f$\pi/4\f$
#define RT2PI 2.5066282746310005024  ///< \f$\sqrt{2\pi}\f$

/* Natural log of 2 */
#define LN2 0.693147180559945 ///< \f$\ln{2}\f$


/* ----------------  NATURAL CONSTANTS  ----------------- */

//! Gravitational Constant \f$G\f$ [m\f$^3\f$ kg\f$^{-1}\f$ s\f$^{-2}\f$]
#define GNEWTON 6.67408e-11

//! Speed of light \f$c\f$ [m/s]
#define CLIGHT 299792458.

//! Mass of the Sun \f$M_\odot\f$ [kg]
#define MSUN 1.9889e30

//! Mass of the Sun \f$M_\odot G c^{-3}\f$ [s]
#define TSUN  4.9169e-6

//! Radius of the Sun \f$_\odot c^{-1}\f$ [s]
#define RSUN 2.32060541029

//! Parsec [m]
#define PC 3.0856775807e16

//! Year [s]
#define YEAR 31457280.0

//! Astronomical unit [m]
#define AU 1.49597870660e11

/**
 \brief Orbital radius of the LISA guiding center [m]
 
 \todo \c RGC belongs in \c LISA.h
*/
#define RGC (1.0*AU)

#endif /* Constants_h */
