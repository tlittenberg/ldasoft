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



#ifndef Constants_h
#define Constants_h

#define FIXME 0
#define SNRPEAK 5

/* --------------  MATHEMATICAL CONSTANTS  -------------- */
/* Square root of 3 */
#define SQ3   1.73205080757

/* Pi's and frinds */
//use math.h (M_PI) for PI
#define PI2   6.283185307179586
#define PIon2 1.57079632679
#define PIon4 0.78539816339
#define RT2PI 2.5066282746310005024  // sqrt(2 Pi)

/* Natural log of 2 */
#define LN2 0.693147180559945


/* ----------------  NATURAL CONSTANTS  ----------------- */

/* Speed of light (m/s) */
#define CLIGHT 299792458.

/* Mass of the Sun (s) */
#define TSUN  4.9169e-6

/* Radius of the Sun (s) */
#define RSUN 2.32060541029

/* Number of meters in a parsec */
#define PC 3.0856775807e16

/* Number of seconds in a year */
#define YEAR 31457280.0

/* Astronomical unit (meters) */
#define AU 1.49597870660e11

/* Gravitational Constant (m^3 kg^-1 s^-2) */
#define G 6.67408e-11


#endif /* Constants_h */
