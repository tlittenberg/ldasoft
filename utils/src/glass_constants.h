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


#ifndef constants_h
#define constants_h

#define FIXME 0
#define SNRPEAK 5 ///< Peak of SNR prior

/* --------------  MATHEMATICAL CONSTANTS  -------------- */
/* Some square roots */
#define SQ3 1.73205080757    ///< \f$\sqrt{3}\f$
#define SQ8 2.82842712474619 ///< \f$\sqrt{8}\f$

/* Pi's and frinds
 use math.h where possible
 M_PI for PI,
 M_PI_2 for PI/2
 M_PI_4 for PI/4
 */
#define PI2   6.283185307179586      ///< \f$2\pi\f$
#define RT2PI 2.5066282746310005024  ///< \f$\sqrt{2\pi}\f$

/* Convert between angles*/
#define RAD2DEG 0.01745329251 ///< 1 deg [rad]
#define DEG2RAD 57.2957795131 ///< 1 rad [deg]

/* ----------------  NATURAL CONSTANTS  ----------------- */

//! Gravitational Constant \f$G\f$ [m\f$^3\f$ kg\f$^{-1}\f$ s\f$^{-2}\f$]
#define GNEWTON 6.67408e-11

//! Speed of light \f$c\f$ [m/s]
#define CLIGHT 299792458.

//! Mass of the Sun \f$M_\odot\f$ [kg]
#define MSUN 1.98848e30

//! Mass of the Sun \f$M_\odot G c^{-3}\f$ [s]
#define TSUN  4.9255e-6

//! Radius of the Sun \f$_\odot c^{-1}\f$ [s]
#define RSUN 2.32060541029

//! Parsec [m]
#define PC 3.0856775815e16

//! Year [s]
#define YEAR 3.15581497632e7

//! Day [s]
#define DAY 86400.0

//! Astronomical unit [m]
#define AU 1.49597870700e11

/**
 \brief Orbital radius of the LISA guiding center [m]
 
 \todo \c RGC belongs in \c LISA.h
*/
#define RGC (1.0*AU)

#endif /* constants_h */
