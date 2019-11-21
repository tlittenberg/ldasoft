/* ----------------  DETECTOR CONSTANTS  ---------------- */

/* Observation time (seconds) */
//#define TOBS 125829120.0
//#define TOBS 62914560
//#define TOBS 31457280

/* Number of data points */
//#define NFFT 125829120
//#define NFFT 62914560
//#define NFFT 31457280

#define DT 1.000000

/* Initial azimuthal position of the guiding center */
#define KAPPA 0.000000

/* Initial orientation of the LISA constellation */
#define LAMBDA 0.000000

/* Mean arm length of the LISA detector (meters) */
#define LARM 2.5e+09 //L3 LISA
//#define LARM 5.0e9 //Classic LISA

/* Photon shot noise power */
#define SPS 8.321000e-23 // L3 LISA
//#define SPS 3.24e-22 // Classic LISA

/* Acceleration noise power */
#define SACC 9.000000e-30 // L3 LISA
//#define SACC 9.000000e-30 // Classic LISA

/* Transfer frequency = c/2piL */
#define FSTAR 0.0190853806 // L3 LISA
//#define FSTAR 0.00954269031847388 // Classic LISA

/* LISA orbital eccentricity = L/(2 sqrt(3) R)*/
#define ECC 0.0048241852 // L3 LISA
//#define ECC 0.0096483704387378 // Classic LISA
