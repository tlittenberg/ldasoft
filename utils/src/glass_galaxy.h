/*
 * Copyright 2024 Tyson B. Littenberg & Neil J. Cornish
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


///@name Galaxy Parameters
///@{
#define GALAXY_RGC 7.2 //!< distance from SSB to GC (kpc)
#define GALAXY_A  0.25 //!< bulge fraction
#define GALAXY_Rb 0.8  //!< bulge radius (kpc)
#define GALAXY_Rd 2.5  //!< disk radius (kpc)
#define GALAXY_Zd 0.4  //!< disk height (kpc)
///@}

#define NSIDE 16 //!< Healpix resolution for galaxy modulation calculations
#define LMAX 4   //!< Maximum l for spherical harmonic decomposition of galaxy modulation
struct GalaxyModulation
{
    int N;
    double alpha_0;
    double alphamax;
    double ***Plm;
    double ***XXR;
    double ***XXI;
    double ***YYR;
    double ***YYI;
    double ***ZZR;
    double ***ZZI;
    double ***XYR;
    double ***XYI;
    double ***YZR;
    double ***YZI;
    double ***XZR;
    double ***XZI;
    double *t;

    struct CubicSpline *XX_spline;
    struct CubicSpline *YY_spline;
    struct CubicSpline *ZZ_spline;
    struct CubicSpline *XY_spline;
    struct CubicSpline *XZ_spline;
    struct CubicSpline *YZ_spline;

    long Npix;
    double *skytheta;
    double *skyphi;
};

double galaxy_distribution(double *x, double bulge_to_disk, double bulge_radius, double disk_radius, double disk_height);
double galaxy_foreground(double f, double A, double f1, double alpha, double fk, double f2);
void rotate_galtoeclip(double *xg, double *xe);
void rotate_ecliptogal(double *xg, double *xe);

void initialize_galaxy_modulation(struct GalaxyModulation *gm, struct Wavelets *wdm, struct Orbit *orbit, double Tobs, double t0);
void galaxy_modulation(struct GalaxyModulation *gm, double *params);
