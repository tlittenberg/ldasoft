/* 
 :tlittenb:20110129 
 - Put constants that need to be changed at top of file
 - Made dependencies explicit
 :tlittenb:20010217
 - Remove BW (Dynamically calculated a la MLDC)
*/
	/* ------------------   RUN CONSTANTS  ------------------ */
#define verbose 1

#define N          256//2048   //number of fourier bins
#define Nblock      16         //number of blocks in window
#define NW          32         //size of noise wings (usually BW/2)

//#define TMAX        800        //maximum temperature during S.A.
#define ANNEALING  69000        //duration of annealing phase

#define CHARACTER_BURNIN 0    //burnin for characterization runs
#define INTRINSIC_BURNIN 69000 //searching over f,dfdt,sky
#define EXTRINSIC_BURNIN 1000 //holding f,dfdt,sky constant

	/* ------------------  DATA CONSTANTS  ------------------ */

/* MLDC_4 2 year observation time */
#define T   62914560.0 
#define Tsq 3.9582418599936e15 
#define sqT 7.93186989303279e3

/* Model dimension */
#define NP 8         //number of parameters
#define NI 2         //number of interferometer channels
#define NN (N/NW)    //number of noise parameters per channel

/* Data dimension */
#define N2   (2*N) //number of time samples 
#define Non2 (N/2) //half number of time samples

/* Size of noise features */
#define NB (N/NN)  //width of noise blocks

/* Differential evolution history */
#define NH 1000  //number of samples stored in history

          /* --------------  MATHEMATICAL CONSTANTS  -------------- */

 /* Set the values of pi */
#define PI 3.141592653589793

#define pi2 6.283185307179586

#define PI2 6.283185307179586

#define pi4 12.566370614359172

#define PIon2 1.570796326794897

#define PIon4 0.785398163397448

#define log2pi 1.837877066409345

/* Natural log of 2 */
#define ln2 0.693147180559945
/* ----------------  NATURAL CONSTANTS  ----------------- */

 /* Speed of light (m/s) */
#define c 299792458.

          /* ----------------  DETECTOR CONSTANTS  ---------------- */

 /* Initial azimuthal position of the guiding center */
#define kappa 0.0

 /* Initial orientation of the LISA constellation */
#define lambda 0.0

 /* Orbital radius of the guiding center */
#define Rgc AU

 /* Mean arm length of the LISA detector (meters) */
#define L 3.0e9

 /* Eccentricity of the LISA spacecraft orbits */
#define ecc (L/(2.0*sqrt(3.0)*Rgc))

 /* Photon shot noise power */
#define Sps 4.0e-22
 
 /* Acceleration noise power */
#define Sacc 9.0e-30

 /* Transfer frequency */ 
#define fstar 0.00954269032

 /* LISA orbital eccentricity */
//#define ec 0.009648370435

 /* LISA modulation frequency */
//#define fm 3.168753575e-8
#define fm 3.178914439e-8 

/* Light travel time between spacecraft */
#define LonC 16.6782047599076

			/* --------------  OPTIMIZATION CONSTANTS  -------------- */

 /* LISA modulation orbital frequency */
#define omega_m 1.99098659045e-07 //2pi*fm

 /* LISA spacecraft maximum height above orbital plane */
#define Ae3 AU*ecc*sqrt(3.0)// = L/2
