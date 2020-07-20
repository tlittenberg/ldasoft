/*
 *  Copyright (C) 2019 Neil J. Cornish, Tyson B. Littenberg (MSFC-ST12)
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "arrays.h"
#include "Detector.h"
#include "Subroutines.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_randist.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf.h>

#define DETECTION_THRESHOLD 7.0 //SNR
#define SKY_LOCATION_THRESHOLD 1.0  //sq deg
#define DISTANCE_THRESHOLD 0.1 //fractional
#define FDOT_THRESHOLD 0.2 //fractional

void FISHER(struct Orbit *orbit, double TOBS, double *params, long N, long M, double SXYZ, double SAE, double *SigmaX, double *SigmaAE, double *DrawAE);

static void print_fisher_results(FILE *fptr, double *params, double *sigmas, double SNR, double TOBS)
{
    /**
     
     1. f (Hz)
     2. co-latitude (rad)
     3. longitude (rad)
     4. amplitude
     5. inclination (rad)
     6. polarization angle (rad)
     7. ref phase (rad)`
     8. df/dt (s^-2)
     9. d^2f/dt^2 (s^-3)
     10. sigma_f (Hz)
     11. sigma_co-latitude (radians)
     12. sigma_longitude (radians)
     13. sigma_amplitude (fractional)
     14. sigma_inclination (radians)
     15. sigma_polarization (radians)
     16. sigma_phase (radians)
     17. sigma_fdot (fractional)
     18. sigma_fddot (fractional)
     19. sigma_omega (sq. deg.)
     20. SNR
     
     */
    
    //print parameters
    for(int i=0; i<9; i++) fprintf(fptr, "%e ", params[i]); //then everybody else
    
    //print uncertainties
    for(int i=0; i<9; i++) fprintf(fptr, "%e ", sigmas[i]); //then everybody else
    
    //print 2D sky location error in square degrees
    fprintf(fptr, "%e ", sigmas[9]);
    
    //print SNR
    fprintf(fptr, "%f\n", SNR);
}

static void get_3d_location(double *params, double *x, double *y, double *z)
{
    double f = params[0];
    double fdot = params[7];
    double amp = params[3];
    double colatitude = params[1];
    double longitude = params[2];
    
    double r  = galactic_binary_dL(f, fdot, amp); //pc
    
    double theta = colatitude;
    double phi = longitude;
    
    *x = r*sin(theta)*cos(phi)/1000.;// convert to kpc;
    *y = r*sin(theta)*sin(phi)/1000.;
    *z = r*cos(theta)/1000.;
}


int main(int argc,char **argv)
{
    double *params, *DrawAE;
    double *X, *A, *E;
    double *SigmaX, *SigmaAE;
    double SNRX, SNR;
    double SAE, SXYZ;
    long M, N, q;
    long i, j, k, mult;
    
    double f, fdot, amp, theta, phi, iota, psi, phase;//, fddot;
    
    int N_X, N_2D_X, N_3D_X, N_FDOT_X, N_FDDOT_X;
    int N_AE, N_2D_AE, N_3D_AE, N_FDOT_AE, N_FDDOT_AE;
    int N_SINGULAR, N_SNR;
    
    //alias some of the more relevant uncertainties
    double dOmega;
    double dfdot;
    double dfddot;
    double dA;
    
    //cartesian ecliptic coordinates for location of source
    double x, y, z;
    
    int imin, imax;
    double *Xc, *AEc, *Xnoise, *Anoise;

    if(argc != 7) KILL("Fisher_Galaxy detections.dat confusion.dat sigmasX.dat sigmasAE.dat DrawAE.dat Orbit.dat\n");

    //parse command line and set up file streams
    FILE* sourceFile    = fopen(argv[1],"r");
    FILE* noiseFile     = fopen(argv[2],"r");
    FILE* resultsFileX  = fopen(argv[3], "w");
    FILE* resultsFileAE = fopen(argv[4], "w");
    FILE* drawFile      = fopen(argv[5], "w");
    FILE* skyFile       = fopen("sky_3d.dat","w");
    FILE* allSkyFile    = fopen("sky_all.dat","w");

    
    printf("***********************************************************************\n");
    printf("*\n");
    printf("* Fisher Parameter Estimation Tool\n");
    printf("*   Candidate Sources:   %s\n",argv[1]);
    printf("*   Confusion Noise Fit: %s\n",argv[2]);
    //printf("*   X-Channel Results:   %s\n",argv[3]);
    printf("*   AE-Channel Results:  %s\n",argv[4]);
    printf("*   Perturbed Params:    %s\n",argv[5]);
    printf("*   Orbit File:          %s\n",argv[6]);
    
    /* Figure out TOBS and NFFT from confusion.dat*/
    double junk;
    double f1,f2;
    fscanf(noiseFile,"%lf%lf%lf%lf%lf", &f1, &junk, &junk, &junk, &junk);
    fscanf(noiseFile,"%lf%lf%lf%lf%lf", &f2, &junk, &junk, &junk, &junk);
    double TOBS = 1./(f2-f1);
    rewind(noiseFile);
    
    printf("*   Observing Time:      %.1f year (%f s)\n",TOBS/YEAR,TOBS);
    printf("*\n");
    printf("***********************************************************************\n");
    
    

    //Data structure for interpolating orbits from file
    struct Orbit *LISAorbit = malloc(sizeof(struct Orbit));
    
    //Set up orbit structure (allocate memory, read file, cubic spline)
    sprintf(LISAorbit->OrbitFileName,"%s",argv[6]);
    initialize_numeric_orbit(LISAorbit);
    
    double L    = LISAorbit->L;
    double fstar= LISAorbit->fstar;
    
    //Read in confusion noise
    imin = (int)floor(FISHERGALAXY_FMIN*TOBS);
    imax = (int)ceil(FISHERGALAXY_FMAX*TOBS);
    
    Xc = double_vector(imax);
    AEc = double_vector(imax);
    Xnoise = double_vector(imax);
    Anoise = double_vector(imax);
    
    noiseFile = fopen(argv[2],"r");
    for(i=imin; i<=imax; i++) fscanf(noiseFile,"%lf%lf%lf%lf%lf", &f, &Xnoise[i], &Xc[i], &Anoise[i], &AEc[i]);
    fclose(noiseFile);
    
    
    //get scale factor for bandwidth set by frequency resolution
    if((TOBS/YEAR) <= 8.0) mult = 8;
    if((TOBS/YEAR) <= 4.0) mult = 4;
    if((TOBS/YEAR) <= 2.0) mult = 2;
    if((TOBS/YEAR) <= 1.0) mult = 1;
    
    //initialize counters for various detection types
    N_X = N_2D_X = N_3D_X = N_FDOT_X = N_FDDOT_X = 0;
    N_AE = N_2D_AE = N_3D_AE = N_FDOT_AE = N_FDDOT_AE = 0;
    N_SINGULAR = N_SNR = 0;
    
    //get number of candidate sources
    int COUNT = 0;
    int loudSNRcount=0;
    while ( !feof(sourceFile) )
    {
        fscanf(sourceFile, "%lf%lf%lf%lf%lf%lf%lf%lf\n", &f, &fdot, &theta, &phi, &amp, &iota, &psi, &phase);
        COUNT++;
    }
    double **loudSNR = double_matrix(19,COUNT-1);
    for(i=0; i<COUNT; i++)for(j=0; j<17; j++) loudSNR[j][i] = 0.0;
    
    rewind(sourceFile);
    
    
    //arrays for parameters and their errors
    params = double_vector(8);
    SigmaX = double_vector(8+1);
    SigmaAE = double_vector(8+1);
    DrawAE = double_vector(8+1);
    
    //Loop through candidate file and comput Fisher for each source
    for(int n=0; n<COUNT; n++)
    {
        if(n%(COUNT/100)==0)printProgress((double)n/(double)COUNT);
        
        //parse input file
        fscanf(sourceFile, "%lf%lf%lf%lf%lf%lf%lf%lf", &f, &fdot, &theta, &phi, &amp, &iota, &psi, &phase);
        
        //map parameters to array for waveform generator
        params[0] = f;
        params[1] = theta;
        params[2] = phi;
        params[3] = amp;
        params[4] = iota;
        params[5] = psi;
        params[6] = phase;
        params[7] = fdot;
        params[8] = 11.0/3.0*fdot*fdot/f;
        
        //Get signal bandwidth
        N = 64*mult;
        if(f > 0.001) N = 128*mult;
        if(f > 0.01)  N = 512*mult;
        if(f > 0.03)  N = 1024*mult;
        if(f > 0.1)   N = 2048*mult;
        M = 2*galactic_binary_bandwidth(LISAorbit->L, LISAorbit->fstar, f, fdot, cos(params[1]), params[3], TOBS, N);
        
        //get noise levels at carrier frequency
        /// \todo use noise from confusion file instaed of compute on the fly
        SAE  = AEnoise_FF(L,fstar,f);
        SXYZ = XYZnoise_FF(L,fstar,f);
        
        //add confusion noise
        q = (long)(f*TOBS);
        if(q>=imin && q<imax)
        {
            SAE  += AEc[q];
            SXYZ += Xc[q];
        }
        
        //get waveforms
        X  = double_vector(2*M);
        A  = double_vector(2*M);
        E  = double_vector(2*M);
        
        FAST_LISA(LISAorbit, TOBS, params, M, X, A, E);
        
        //get SNRs
        SNRX = sqrt(Sum(X,X,M,SXYZ,TOBS));
        SNR  = sqrt(Sum(A,A,M,SAE,TOBS)+Sum(E,E,M,SAE,TOBS));
        
        //free waveform memory for next iteration
        free_double_vector(X);
        free_double_vector(A);
        free_double_vector(E);
        
        //select events above SNR detection threshold
        if(SNR > DETECTION_THRESHOLD)
        {
            //update counter
            N_AE++;
            
            //get 3d location for detached systems
            if(fdot>0.0)
            {
                get_3d_location(params, &x, &y, &z);
                fprintf(allSkyFile,"%lg %lg %lg\n",x,y,z);
                fflush(allSkyFile);
            }
            
            //get Fisher estimates for errors
            FISHER(LISAorbit, TOBS, params, N, M, SXYZ, SAE, SigmaX, SigmaAE, DrawAE);
            
            //nan check
            int err=0;
            for(i=0; i<10; i++) if(isnan(SigmaAE[i])) err=1;
            if(err)
            {
                N_SINGULAR++;
                continue;
            }
            
            //crazy snr check
            ///\todo really need to fix this...
            if(SNR>1e7)
            {
                N_SNR++;
                continue;
            }


            for(i=0; i<9; i++)
            {
                loudSNR[i][loudSNRcount] = params[i];
                loudSNR[i+9][loudSNRcount] = SigmaAE[i];
            }
            loudSNR[18][loudSNRcount] = SigmaAE[9];
            loudSNR[19][loudSNRcount] = SNR;
            loudSNRcount++;
            
            
            //print main output file with parameters & errors
            print_fisher_results(resultsFileAE,params,SigmaAE,SNR,TOBS);
            
            //alias some of the more relevant uncertainties
            dOmega = SigmaAE[9];
            dfdot  = SigmaAE[7];
            dfddot = SigmaAE[8];
            dA     = SigmaAE[3];
            
            //check for well-localized sources
            if(dOmega < SKY_LOCATION_THRESHOLD)
            {
                N_2D_AE++;   // 2D mapping
                
                //check for 3D mapping
                if( (fdot > 0.0) && (dfdot < DISTANCE_THRESHOLD) && (dA < DISTANCE_THRESHOLD))
                {
                    N_3D_AE++;  // 3D mapping
                    get_3d_location(params, &x, &y, &z);
                    fprintf(skyFile,"%lg %lg %lg\n",x,y,z);
                    fflush(skyFile);
                }
            }
            
            if(dfdot  < FDOT_THRESHOLD) N_FDOT_AE++;   // measure fdot to 20%
            if(dfddot < FDOT_THRESHOLD) N_FDDOT_AE++;  // measure fddot to 20%
            
            
            //Cut on 4-link detections
            if(SNRX > DETECTION_THRESHOLD)
            {
                N_X++;  /* SNR > 7 with just Michelson */
                
                print_fisher_results(resultsFileX,params,SigmaX,SNRX,TOBS);
                
                //alias some of the more relevant uncertainties
                dOmega = SigmaX[9];
                dfdot  = SigmaX[7];
                dfddot = SigmaX[8];
                dA     = SigmaX[3];
                
                if(dOmega < SKY_LOCATION_THRESHOLD)
                {
                    N_2D_X++;   // 2D mapping
                    if((fdot > 0.0) && (dfdot < DISTANCE_THRESHOLD) && (dA < DISTANCE_THRESHOLD))
                        N_3D_X++;  // 3D mapping
                }
                
                if(dfdot  < FDOT_THRESHOLD) N_FDOT_X++;
                if(dfddot < FDOT_THRESHOLD) N_FDDOT_X++;
                
            }
            
            
            //print fair draw from Fisher
            f     = DrawAE[0];
            fdot  = DrawAE[7];
            theta = DrawAE[1];
            phi   = DrawAE[2];
            amp   = DrawAE[3];
            iota  = DrawAE[4];
            psi   = DrawAE[5];
            phase = DrawAE[6];
            fprintf(drawFile, "%.16f %.10e %f %f %e %f %f %f\n", f, fdot, theta, phi, amp, iota, psi, phase);
            
        } //AE SNR>7
    }
    printf("\n");
    printf("*************** AE *****************\n");
    printf("SNR > 7  = %d\n", N_AE);
    printf("2D sky mapping = %d\n", N_2D_AE);
    printf("3D sky mapping = %d\n", N_3D_AE);
    printf("fdot measured to 20 percent = %d\n", N_FDOT_AE);
    printf("FIXME: fddot measured to 20 percent = %d\n", N_FDDOT_AE);
    printf("Sources with singular matrices = %d\n", N_SINGULAR);
    printf("Sources with absurd SNRs = %d\n", N_SNR);

    //  printf("*************** X *****************\n");
    //  printf("SNR > 7  = %d\n", N_X);
    //  printf("2D sky mapping = %d\n", N_2D_X);
    //  printf("3D sky mapping = %d\n", N_3D_X);
    //  printf("fdot measured to 20 percent = %d\n", N_FDOT_X);
    //  printf("fddot measured to 20 percent = %d\n", N_FDDOT_X);
    
    
    
    fclose(resultsFileAE);
    fclose(resultsFileX);
    fclose(skyFile);
    fclose(allSkyFile);
    
    
    //Print 100 loudest detections & errors
    
    //arrays to hold sorted indicies
    size_t *sort_ascending  = malloc(loudSNRcount*(sizeof(size_t)));
    size_t *sort_descending = malloc(loudSNRcount*(sizeof(size_t)));
    
    gsl_sort_index(sort_ascending, loudSNR[19], 1, loudSNRcount);
    
    for(i=0; i<loudSNRcount; i++)
    {
        j = loudSNRcount - 1 - i;
        sort_descending[i] = sort_ascending[j];
    }
    
    FILE *outfile = fopen("Brightest.dat","w");
    
    int MAX = 100;
    if(MAX>loudSNRcount) MAX=loudSNRcount;
    for(i=0; i<MAX; i++)
    {
        j = sort_descending[i];
        for(k=0; k<20; k++) fprintf(outfile,"%g ",loudSNR[k][j]);
        fprintf(outfile,"\n");
    }
    
    return 0;
}


void FISHER(struct Orbit *orbit, double TOBS, double *Params, long N, long M, double SXYZ, double SAE, double *SigmaX, double *SigmaAE, double *DrawAE)
{
    
    int i, j;
    double *ParamsP, *ParamsM;
    double *templateP1, *templateP2;
    double *templateM1, *templateM2;
    double *XP, *XM;
    double **XD, **temD1, **temD2;
    double **Fisher1, **Fisher2;
    double **Cov1,**Cov2;
    double epsilon;
    int d;
    
    d = 8;
    
    epsilon = 1.0e-6;
    
    // Plus and minus parameters:
    
    ParamsP = double_vector(d-1);
    ParamsM = double_vector(d-1);
    
    // Plus and minus templates for each channel:
    
    XP = double_vector(2*M);
    templateP1 = double_vector(2*M);
    templateP2 = double_vector(2*M);
    
    XM = double_vector(2*M);
    templateM1 = double_vector(2*M);
    templateM2 = double_vector(2*M);
    
    // Differences for each channel
    
    XD = double_matrix(d-1, 2*M);
    temD1 = double_matrix(d-1, 2*M);
    temD2 = double_matrix(d-1, 2*M);
    
    // The individual Fisher matrices:
    
    Fisher1 = double_matrix(d-1,d-1);
    Fisher2 = double_matrix(d-1,d-1);
    Cov1 = double_matrix(d-1,d-1);
    Cov2 = double_matrix(d-1,d-1);
    
    // Random draw to perturb params by Fisher
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc(T);
    gsl_rng_env_setup();
    gsl_rng_set (r, time(0));
    
    double *u = malloc(d*sizeof(double));
    for (i=0; i<d; i++) u[i] = gsl_ran_ugaussian(r);
    
    for (i = 0 ; i < d ; i++)
    {
        for (j = 0 ; j < d ; j++)
        {
            ParamsP[j] = Params[j];
            ParamsM[j] = Params[j];
        }
        
        if(i==0)
        {
            ParamsP[0] *= TOBS;
            ParamsM[0] *= TOBS;
        }
        if(i==3)
        {
            ParamsP[3] = log(ParamsP[3]);
            ParamsM[3] = log(ParamsM[3]);
        }
        if(i==7)
        {
            ParamsP[7] *= TOBS*TOBS;
            ParamsM[7] *= TOBS*TOBS;
        }
        if(i==8)
        {
            ParamsP[8] *= TOBS*TOBS*TOBS;
            ParamsM[8] *= TOBS*TOBS*TOBS;
        }
        
        
        ParamsP[i] += epsilon;
        ParamsM[i] -= epsilon;
        
        
        if(i==0)
        {
            ParamsP[0] /= TOBS;
            ParamsM[0] /= TOBS;
        }
        if(i==3)
        {
            ParamsP[3] = exp(ParamsP[3]);
            ParamsM[3] = exp(ParamsM[3]);
        }
        if(i==7)
        {
            ParamsP[7] /= TOBS*TOBS;
            ParamsM[7] /= TOBS*TOBS;
        }
        if(i==8)
        {
            ParamsP[8] /= TOBS*TOBS*TOBS;
            ParamsM[8] /= TOBS*TOBS*TOBS;
        }
        
        FAST_LISA(orbit, TOBS, ParamsP, M, XP, templateP1, templateP2);
        FAST_LISA(orbit, TOBS, ParamsM, M, XM, templateM1, templateM2);
        
        for (j = 0 ; j < 2*M ; j++)
        {
            XD[i][j] = XP[j] - XM[j];
            temD1[i][j] = templateP1[j] - templateM1[j];
            temD2[i][j] = templateP2[j] - templateM2[j];
        }
        
    }
    
    free_double_vector(XP);
    free_double_vector(XM);
    free_double_vector(templateP1);
    free_double_vector(templateP2);
    free_double_vector(templateM1);
    free_double_vector(templateM2);
    
    // Calculate Fisher1 (using just one channel)
    
    for (i = 0 ; i < d ; i++)
    {
        
        for (j = i ; j < d ; j++)
        {
            Fisher1[i][j] =  Sum(XD[i],XD[j],M,SXYZ,TOBS);
        }
    }
    
    free_double_matrix(XD,d-1);
    
    
    for (i = 0; i < d; i++)
    {
        for (j = i+1; j < d; j++)
            Fisher1[j][i] = Fisher1[i][j];
    }
    
    for (i = 0; i < d; i++)
    {
        for (j = 0; j < d; j++)
            Fisher1[i][j] /= (4.0*epsilon*epsilon);
    }
    
    // Calculate Fisher2: (using A + E channels)
    
    for (i = 0 ; i < d ; i++)
    {
        for (j = i ; j < d ; j++)
        {
            Fisher2[i][j] =  Sum(temD1[i],temD1[j],M,SAE,TOBS)+Sum(temD2[i],temD2[j],M,SAE,TOBS);
        }
    }
    
    free_double_matrix(temD1,d-1);
    free_double_matrix(temD2,d-1);
    
    for (i = 0; i < d; i++)
    {
        for (j = i+1; j < d; j++)
            Fisher2[j][i] = Fisher2[i][j];
    }
    
    for (i = 0; i < d; i++)
    {
        for (j = 0; j < d; j++)
        {
            Fisher2[i][j] /= (4.0*epsilon*epsilon);
        }
    }
    
    double **evector = double_matrix(d-1,d-1);
    double *evalue = double_vector(d-1);
    matrix_eigenstuff(Fisher2,evector,evalue,d);
    for (i = 0; i < d; i++)for (j = 0; j < d; j++) Cov2[i][j] = Fisher2[i][j];
    
    //get perturbed params (ignoring fddot)
    
    for (i = 0; i < 8; i++)
    {
        DrawAE[i] = 0.0;
        for (j = 0; j < 8; j++)
        {
            /**
             \todo debug and turn back on Draw
             */
            DrawAE[i] += 0.0;//evector[i][j]*u[j]/sqrt(evalue[j])/sqrt(8.);
        }
        if(i==0)DrawAE[i] = DrawAE[i]/TOBS + Params[i];
        else if(i==3)DrawAE[i] = exp(DrawAE[i] + log(Params[i]));
        else if(i==7)DrawAE[i] = DrawAE[i]/TOBS/TOBS + Params[i];
        else DrawAE[i] = DrawAE[i] + Params[i];
        
    }
    
    free_double_vector(evalue);
    free_double_matrix(evector,d-1);
    
    
    
    //get Fisher estimated errors
    for (i = 0; i < d; i++) SigmaAE[i] = sqrt(Cov2[i][i]);
    
    SigmaAE[9] = 2.0*M_PI*sin(Params[1])*sqrt(Cov2[1][1]*Cov2[2][2] - Cov2[2][1]*Cov2[2][1])/RAD2DEG/RAD2DEG;
    
    d = 8;
    
    
    evector = double_matrix(d-1,d-1);
    evalue = double_vector(d-1);
    matrix_eigenstuff(Fisher1,evector,evalue,d);
    for (i = 0; i < d; i++)for (j = 0; j < d; j++) Cov1[i][j] = Fisher1[i][j];
    free_double_vector(evalue);
    free_double_matrix(evector,d-1);
    
    for (i = 0; i < d; i++) SigmaX[i] = sqrt(Cov1[i][i]);
    
    SigmaX[9] = 2.0*M_PI*sin(Params[1])*sqrt(Cov1[1][1]*Cov1[2][2] - Cov1[2][1]*Cov1[2][1])/RAD2DEG/RAD2DEG;
    
    
    for (j = 0; j < 2; j++)
    {
        if(SigmaX[d-1] > 0.5)
        {
            SigmaX[d-1] = 1.0e10;
            d -= 1;
            
            evector = double_matrix(d-1,d-1);
            evalue = double_vector(d-1);
            matrix_eigenstuff(Fisher1,evector,evalue,d);
            for (i = 0; i < d; i++)for (j = 0; j < d; j++) Cov1[i][j] = Fisher1[i][j];
            free_double_vector(evalue);
            free_double_matrix(evector,d-1);
            
            for (i = 0; i < d; i++) SigmaX[i] = sqrt(Cov2[i][i]);
            
        }
    }
    
    
    //convert back to desired units
    SigmaAE[0]/=TOBS; //absolute uncertainty in frequency (Hz)
    SigmaAE[7]/=fabs(Params[7]*TOBS*TOBS); //fractional uncertainty in fdot
    SigmaAE[8]/=fabs(Params[8]*TOBS*TOBS*TOBS); //fractional uncertainty in fddot
    
    SigmaX[0]/=TOBS; //absolute uncertainty in frequency (Hz)
    SigmaX[7]/=fabs(Params[7]*TOBS*TOBS); //fractional uncertainty in fdot
    SigmaX[8]/=fabs(Params[8]*TOBS*TOBS*TOBS); //fractional uncertainty in fddot
    
    
    d = 8;
    
    free_double_vector(ParamsP);
    free_double_vector(ParamsM);
    
    free_double_matrix(Fisher1,d-1);
    free_double_matrix(Fisher2,d-1);
    
    free_double_matrix(Cov1, d-1);
    free_double_matrix(Cov2, d-1);
    
}


