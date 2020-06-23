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


void FISHER(struct Orbit *orbit, double TOBS, double *params, long N, long M, double SXYZ, double SAE, double *SigmaX, double *SigmaAE, double *DrawAE);

int main(int argc,char **argv)
{
  double *params, *DrawAE;
  double *XLS;
  double *AA, *EE;
  double *SigmaX, *SigmaAE;
  double fonfs, Sm, Acut;
  double SNRX, SNR;
  double SAE, SXYZ;
  long M, N, q;
  long i, j, k, mult;
  double fix, Amp;
  double x, y ,z, r;
  double f, fdot, A=0.0, theta, phi, iota, psi, phase, fddot;
  int cntx, cnt1, cnt2, cnt3, cnt4, cnt5, cnt6, cnt7, cnt8, cnt9, cnt10;
  int imin, imax;
  double *Xc, *AEc, *Xnoise, *Anoise;
  FILE* vfile;
  FILE* afile;
  FILE* pfile;
  //FILE* xfile;
  FILE* sky;
  FILE* sky2;
  FILE* amps;
  
  
  
  if(argc != 7) KILL("Fisher_Galaxy detections.dat confusion.dat sigmasX.dat sigmasAE.dat DrawAE.dat Orbit.dat\n");

  printf("***********************************************************************\n");
  printf("*\n");
  printf("* Fisher Parameter Estimation Tool\n");
  printf("*   Candidate Sources:   %s\n",argv[1]);
  printf("*   Confusion Noise Fit: %s\n",argv[2]);
  //printf("*   X-Channel Results:   %s\n",argv[3]);
  printf("*   AE-Channel Results:  %s\n",argv[4]);
  printf("*   Perturbed Params:    %s\n",argv[5]);
  printf("*   Orbit File:          %s\n",argv[6]);
  
  /* Figure out TOBS and NFFT */
  FILE *Infile = fopen(argv[2],"r");
  double junk;
  double f1,f2;
  fscanf(Infile,"%lf%lf%lf%lf%lf\n", &f1, &junk, &junk, &junk, &junk);
  fscanf(Infile,"%lf%lf%lf%lf%lf\n", &f2, &junk, &junk, &junk, &junk);
  double TOBS = 1./(f2-f1);
  fclose(Infile);
  /*****************************/

  printf("*   Observing Time:      %.1f year (%f s)\n",TOBS/YEAR,TOBS);
  printf("*\n");
  printf("***********************************************************************\n");
  
  vfile = fopen(argv[1],"r");
  
  //Data structure for interpolating orbits from file
  struct Orbit *LISAorbit = malloc(sizeof(struct Orbit));

  //Set up orbit structure (allocate memory, read file, cubic spline)
  sprintf(LISAorbit->OrbitFileName,"%s",argv[6]);
  initialize_numeric_orbit(LISAorbit);

  double L    = LISAorbit->L;
  double fstar= LISAorbit->fstar;
  
  imax = (int)ceil(4.0e-2*TOBS);
  imin = (int)floor(1.0e-4*TOBS);
  
  Xc = double_vector(imax);  AEc = double_vector(imax);
  Xnoise = double_vector(imax);  Anoise = double_vector(imax);
  
  amps = fopen(argv[2],"r");
  for(i=imin; i<= imax; i++)
  {
    fscanf(amps,"%lf%lf%lf%lf%lf\n", &f, &Xnoise[i], &Xc[i], &Anoise[i], &AEc[i]);
    SAE  = AEnoise_FF(L,fstar,f);
    SXYZ = XYZnoise_FF(L,fstar,f);

    // printf("%e %e %e\n", f, SXYZ+Xc[i], SAE+AEc[i]);
  }
  fclose(amps);
  
  params = double_vector(8);
  SigmaX = double_vector(8+1);
  SigmaAE = double_vector(8+1);
  DrawAE = double_vector(8+1);
  
  
  fix = M_PI/180.0;  /* one degree in radians */
  
  
  if((TOBS/YEAR) <= 8.0) mult = 8;
  if((TOBS/YEAR) <= 4.0) mult = 4;
  if((TOBS/YEAR) <= 2.0) mult = 2;
  if((TOBS/YEAR) <= 1.0) mult = 1;
  
  
  //xfile = fopen(argv[3], "w");
  afile = fopen(argv[4], "w");
  pfile = fopen(argv[5], "w");
  sky = fopen("sky_3d.dat","w");
  sky2 = fopen("sky_all.dat","w");
  
  
  cntx = cnt1 = cnt2 = cnt3 = cnt4 = cnt5 = cnt6 = cnt7 = cnt8 = cnt9 = cnt10 = 0;
  while ( !feof(vfile) )
  {
    
    fscanf(vfile, "%lf%lf%lf%lf%lf%lf%lf%lf", &f, &fdot, &theta, &phi, &Amp, &iota, &psi, &phase);
    cntx++;
  }
  
  int COUNT = cntx--;
  int loudSNRcount=0;
  double **loudSNR = double_matrix(16,COUNT-1);
  for(i=0; i<COUNT; i++)for(j=0; j<17; j++) loudSNR[j][i] = 0.0;
  
  rewind(vfile);
  cntx=0;
  
  
  
  double th;// = M_PI - params[1];
  double ph;// = params[2];
  //  double r;//  = galactic_binary_dL(params[0], params[7], params[3]);
  //  double x;// = r*sin(th)*cos(ph);// + GC;
  //  double y;// = r*sin(th)*sin(ph);
  //  double z;// = r*cos(th);
  
  
  //count lines in file
  int NSIM = 0;
  while ( !feof(vfile) )
  {
    fscanf(vfile, "%lf%lf%lf%lf%lf%lf%lf%lf", &f, &fdot, &theta, &phi, &A, &iota, &psi, &phase);
    NSIM++;
  }
  rewind(vfile);
  NSIM--;
  
  for(int n=0; n<NSIM; n++)
  {
    if(n%(NSIM/100)==0)printProgress((double)n/(double)NSIM);

    fscanf(vfile, "%lf%lf%lf%lf%lf%lf%lf%lf", &f, &fdot, &theta, &phi, &Amp, &iota, &psi, &phase);
    
    params[0] = f;
    params[1] = theta;
    params[2] = phi;
    params[3] = Amp;
    params[4] = iota;
    params[5] = psi;
    params[6] = phase;
    params[7] = fdot;
    params[8] = 11.0/3.0*fdot*fdot/f;
    
    /*
     fscanf(vfile, "%lf%lf%lf%lf%lf%lf%lf%lf\n", &Amp, &iota, &psi, &phase, &f, &fdot, &theta, &phi);
     f = f/T;
     fdot = fdot/T/T;
     Amp = exp(Amp);
     iota = cos(iota);
     theta = acos(theta);
     params[0] = f;
     params[1] = theta;//0.5*M_PI - theta;
     params[2] = phi;
     params[3] = Amp;
     params[4] = iota;
     params[5] = psi;
     params[6] = phase;
     params[7] = fdot;
     params[8] = 11.0/3.0*params[7]*params[7]/params[0];
     */
    cntx++;
    
    N = 64*mult;
    if(f > 0.001) N = 128*mult;
    if(f > 0.01) N = 512*mult;
    if(f > 0.03) N = 1024*mult;
    if(f > 0.1) N = 2048*mult;
    
    fonfs = f/fstar;
    
    q = (long)(f*TOBS);
    
    SAE  = AEnoise_FF(L,fstar,f);
    SXYZ = XYZnoise_FF(L,fstar,f);

    if(q <= imax)
    {
      SAE  += AEc[q];
      SXYZ += Xc[q];
    }
    
    
    /*  calculate michelson noise  */
    Sm = SXYZ/(4.0*sin(f/fstar)*sin(f/fstar));
    
    Acut = A*sqrt(TOBS/Sm);
    
    M = 2*galactic_binary_bandwidth(LISAorbit->L, LISAorbit->fstar, f, fdot, cos(params[1]), params[3], TOBS, N);
    if(M>N)
    {
      printf("bandwidth=%li, max size=%li\n",M,N);
      exit(0);
    }
    
    XLS = double_vector(2*M);
    AA  = double_vector(2*M);
    EE  = double_vector(2*M);
    
    FAST_LISA(LISAorbit, TOBS, params, M, XLS, AA, EE);

    SNRX = sqrt(Sum(XLS,XLS,M,SXYZ,TOBS));
    SNR  = sqrt(Sum(AA,AA,M,SAE,TOBS)+Sum(EE,EE,M,SAE,TOBS));
    
    free_double_vector(XLS);
    free_double_vector(AA);
    free_double_vector(EE);
    
    if(SNR > 7.0)
    {
      
      cnt1++;
      
      th = M_PI - params[1];
      ph = params[2];
      r  = galactic_binary_dL(params[0], params[7], params[3]);
      x = r*sin(th)*cos(ph);// + GC;
      y = r*sin(th)*sin(ph);
      z = r*cos(th);
      fprintf(sky2,"%lg %lg %lg\n",x/1000.,y/1000.,z/1000.);
      fflush(sky2);
      
      FISHER(LISAorbit, TOBS, params, N, M, SXYZ, SAE, SigmaX, SigmaAE, DrawAE);
      
      loudSNR[0][loudSNRcount] = log10(SNR);        //1:2
      loudSNR[1][loudSNRcount] = params[0];        //3:4
      loudSNR[2][loudSNRcount] = fabs(params[7]*TOBS*TOBS);      //5:6
      loudSNR[3][loudSNRcount] = theta;  //7:8
      loudSNR[4][loudSNRcount] = phi;        //9:10
      loudSNR[5][loudSNRcount] = Amp;        //11:12
      loudSNR[6][loudSNRcount] = iota;      //13:14
      loudSNR[7][loudSNRcount] = psi;        //15:16
      loudSNR[8][loudSNRcount] = phase;      //17:18
      loudSNR[9][loudSNRcount]  = log10(SigmaAE[0]/TOBS/params[0]);  //19:20
      loudSNR[10][loudSNRcount] = log10(SigmaAE[7]/fabs(params[7]*TOBS*TOBS)); //21:22
      loudSNR[11][loudSNRcount] = log10(SigmaAE[1]);    //23:24
      loudSNR[12][loudSNRcount] = log10(SigmaAE[2]);    //25:26
      loudSNR[13][loudSNRcount] = log10(SigmaAE[3]);  //27:28
      loudSNR[14][loudSNRcount] = SigmaAE[4];    //29:30
      loudSNR[15][loudSNRcount] = SigmaAE[5];    //31:32
      loudSNR[16][loudSNRcount] = SigmaAE[6];    //33:34
      loudSNRcount++;
      
      f     = DrawAE[0];
      fdot  = DrawAE[7];
      theta = DrawAE[1];
      phi   = DrawAE[2];
      A     = DrawAE[3];
      iota  = DrawAE[4];
      psi   = DrawAE[5];
      phase = DrawAE[6];
      fprintf(pfile, "%.16f %.10e %f %f %e %f %f %f\n", f, fdot, theta, phi, A, iota, psi, phase);

//      for(i = 0; i < 8; i++)
//      {
//        fprintf(pfile, "%.16g ", DrawAE[i]);
//      }
//      fprintf(pfile,"\n");
      
      fprintf(afile, "%e %e %e ", params[0], params[7], params[8]);
      for(i = 1; i < 7; i++)
      {
        fprintf(afile, "%e ", params[i]);
      }
      
      fdot = SigmaAE[7]/(fabs(params[7]*TOBS*TOBS));  // fraction uncertainity in fdot
      fddot = SigmaAE[8]/(fabs(params[8]*TOBS*TOBS*TOBS));  // fraction uncertainity in fddot
      fprintf(afile, "%e ", SigmaAE[0]/TOBS);  // convert back to Hz
      fprintf(afile, "%e %e ", fdot, fddot);

      for(i = 1; i < 7; i++)
      {
        fprintf(afile, "%e ", SigmaAE[i]);
      }
      fprintf(afile, "%e ", SigmaAE[9]/fix/fix);
      fprintf(afile, "%f\n", SNR);
      
      if((SigmaAE[1] < fix) && (SigmaAE[2] < fix))
      {
        cnt2++;   // 2D mapping
        //        if(fdot>0)
        //        {
        //        double th = M_PI - params[1];
        //        double ph = params[2];
        //        double r  = galactic_binary_dL(params[0], params[7], params[3]);
        //        double x = r*sin(th)*cos(ph);// + GC;
        //        double y = r*sin(th)*sin(ph);
        //        double z = r*cos(th);
        //        fprintf(sky,"%lg %lg %lg\n",x/1000.,y/1000.,z/1000.);
        //        }
        if((fdot < 0.1) && (SigmaAE[3] < 0.1))
        {
          cnt3++;  // 3D mapping
          th = M_PI - params[1];
          ph = params[2];
          r  = galactic_binary_dL(params[0], params[7], params[3]);
          x = r*sin(th)*cos(ph);// + GC;
          y = r*sin(th)*sin(ph);
          z = r*cos(th);
          fprintf(sky,"%lg %lg %lg\n",x/1000.,y/1000.,z/1000.);
          fflush(sky);
        }
      } //less than 1 sq degree
      
      if(fdot  < 0.2) cnt4++;  // measure fdot to 20%
      if(fddot < 0.2) cnt5++;  // measure fddot to 20%
      
//      if(SNRX > 7.0)
//      {
//        cnt6++;  /* SNR > 7 with just Michelson */
//
//        loudSNR[0][loudSNRcount] = log10(SNR);        //1:2
//        loudSNR[1][loudSNRcount] = params[0];        //3:4
//        loudSNR[2][loudSNRcount] = fabs(params[7]*TOBS*TOBS);      //5:6
//        loudSNR[3][loudSNRcount] = 0.5*M_PI-theta;  //7:8
//        loudSNR[4][loudSNRcount] = phi;        //9:10
//        loudSNR[5][loudSNRcount] = Amp;        //11:12
//        loudSNR[6][loudSNRcount] = iota;      //13:14
//        loudSNR[7][loudSNRcount] = psi;        //15:16
//        loudSNR[8][loudSNRcount] = phase;      //17:18
//        loudSNR[9][loudSNRcount]  = log10(SigmaX[0]/TOBS/params[0]);  //19:20
//        loudSNR[10][loudSNRcount] = log10(SigmaX[7]/fabs(params[7]*TOBS*TOBS)); //21:22
//        loudSNR[11][loudSNRcount] = log10(SigmaX[1]);    //23:24
//        loudSNR[12][loudSNRcount] = log10(SigmaX[2]);    //25:26
//        loudSNR[13][loudSNRcount] = log10(SigmaX[3]);  //27:28
//        loudSNR[14][loudSNRcount] = SigmaX[4];    //29:30
//        loudSNR[15][loudSNRcount] = SigmaX[5];    //31:32
//        loudSNR[16][loudSNRcount] = SigmaX[6];    //33:34
//        loudSNRcount++;
//
//
//
//        fprintf(xfile, "%e %e %e ", params[0], params[7], params[8]);
//        for(i = 1; i < 7; i++)
//        {
//          fprintf(xfile, "%e ", params[i]);
//        }
//
//        fdot = SigmaX[7]/(fabs(params[7]*TOBS*TOBS));  // fraction uncertainity in fdot
//        fddot = SigmaX[8]/(fabs(params[8]*TOBS*TOBS*TOBS));  // fraction uncertainity in fddot
//        fprintf(xfile, "%e ", SigmaX[0]/TOBS);  // convert back to Hz
//        fprintf(xfile, "%e %e ", fdot, fddot);
//        //SigmaX[3] *= Amp;
//
//        for(i = 1; i < 7; i++)
//        {
//          fprintf(xfile, "%e ", SigmaX[i]);
//        }
//        fprintf(xfile, "%e ", SigmaX[9]/fix/fix);
//        fprintf(xfile, "%f\n", SNRX);
//
//        if((SigmaX[1] < fix) && (SigmaX[2] < fix))
//        {
//          cnt7++;   // 2D mapping
//          if((fdot < 0.1) && (SigmaAE[3] < 0.1))
//          {
//            cnt8++;  // 3D mapping
//          }
//        }
//
//        if(fdot < 0.2)
//        {
//          cnt9++;  // measure fdot to 20%
//        }
//
//        if(fddot < 0.2)
//        {
//          cnt10++;  // measure fddot to 20%
//          // printf("%e %f\n", f, SNRX);
//        }
//      } // Michelson SNR > 7
    } //AE SNR>7
  }
  printf("\n");
  printf("*************** AE *****************\n");
  printf("SNR > 7  = %d\n", cnt1);
  printf("2D sky mapping = %d\n", cnt2);
  printf("3D sky mapping = %d\n", cnt3);
  printf("fdot measured to 20 percent = %d\n", cnt4);
  printf("FIXME: fddot measured to 20 percent = %d\n", cnt5);
  
//  printf("*************** X *****************\n");
//  printf("SNR > 7  = %d\n", cnt6);
//  printf("2D sky mapping = %d\n", cnt7);
//  printf("3D sky mapping = %d\n", cnt8);
//  printf("fdot measured to 20 percent = %d\n", cnt9);
//  printf("fddot measured to 20 percent = %d\n", cnt10);
  
  
  
  fclose(afile);
  //fclose(xfile);
  
  
  //sort cell->occupancy so I search from most occupied to least occupied cell
  size_t *sort_ascending  = malloc(COUNT*(sizeof(size_t)));
  size_t *sort_descending = malloc(COUNT*(sizeof(size_t)));
  
  gsl_sort_index(sort_ascending, loudSNR[0], 1, COUNT);
  
  for(i=0; i<COUNT; i++)
  {
    j = COUNT - 1 - i;
    sort_descending[i] = sort_ascending[j];
  }
  
  FILE *outfile = fopen("Brightest.dat","w");
    
    int MAX = 100;
    if(MAX>COUNT)MAX=COUNT;
  for(i=0; i<100; i++)
  {
    j = sort_descending[i];
    for(k=0; k<17; k++) fprintf(outfile,"%g ",loudSNR[k][j]);
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
  
  epsilon = 1.0e-5;
  
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
  
  SigmaAE[9] = 2.0*M_PI*sin(Params[1])*sqrt(Cov2[1][1]*Cov2[2][2] - Cov2[2][1]*Cov2[2][1]);
  
  d = 8;

  
  evector = double_matrix(d-1,d-1);
  evalue = double_vector(d-1);
  matrix_eigenstuff(Fisher1,evector,evalue,d);
  for (i = 0; i < d; i++)for (j = 0; j < d; j++) Cov1[i][j] = Fisher1[i][j];
  free_double_vector(evalue);
  free_double_matrix(evector,d-1);
  
  for (i = 0; i < d; i++) SigmaX[i] = sqrt(Cov1[i][i]);
  
  SigmaX[9] = 2.0*M_PI*sin(Params[1])*sqrt(Cov1[1][1]*Cov1[2][2] - Cov1[2][1]*Cov1[2][1]);
  
  
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
  
  d = 8;
  
  free_double_vector(ParamsP);
  free_double_vector(ParamsM);
  
  free_double_matrix(Fisher1,d-1);
  free_double_matrix(Fisher2,d-1);
  
  free_double_matrix(Cov1, d-1);
  free_double_matrix(Cov2, d-1);
  
}


