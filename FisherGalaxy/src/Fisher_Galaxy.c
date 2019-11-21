/***************************************************************************/
/*                                                                         */
/*              Fisher_Galaxy.c, Version 2.3, 4/28/2011                    */
/*             Written by Neil Cornish & Tyson Littenberg                  */
/*                                                                         */
/* gcc -O2 -o Fisher_Galaxy Fisher_Galaxy.c arrays.c -lcblas -lclapack -lm */
/*                                                                         */
/***************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "arrays.h"
#include "Constants.h"
#include "Detector.h"
#include "Subroutines.h"

#include <LISA.h>
#include <GalacticBinary.h>
#include <GalacticBinaryMath.h>
#include <GalacticBinaryWaveform.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf.h>


void FISHER(struct Orbit *orbit, double TOBS, double *params, long N, long M, double SXYZ, double SAE, double *SigmaX, double *SigmaAE, double *DrawAE);
void indexx(unsigned long n, double *arr, unsigned long *indx);


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
  double *Xc, *AEc, *Xnoise, *AEnoise;
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

  printf("*   Observing Time:      %.1f year (%f s)\n",TOBS/year,TOBS);
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
  
  Xc = dvector(imin,imax);  AEc = dvector(imin,imax);
  Xnoise = dvector(imin,imax);  AEnoise = dvector(imin,imax);
  
  amps = fopen(argv[2],"r");
  for(i=imin; i<= imax; i++)
  {
    fscanf(amps,"%lf%lf%lf%lf%lf\n", &f, &Xnoise[i], &Xc[i], &AEnoise[i], &AEc[i]);
    instrument_noise(f, fstar, L, &SAE, &SXYZ);
    // printf("%e %e %e\n", f, SXYZ+Xc[i], SAE+AEc[i]);
  }
  fclose(amps);
  
  params = dvector(0,8);
  SigmaX = dvector(0,8+1);
  SigmaAE = dvector(0,8+1);
  DrawAE = dvector(0,8+1);
  
  
  fix = pi/180.0;  /* one degree in radians */
  
  
  if((TOBS/year) <= 8.0) mult = 8;
  if((TOBS/year) <= 4.0) mult = 4;
  if((TOBS/year) <= 2.0) mult = 2;
  if((TOBS/year) <= 1.0) mult = 1;
  
  
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
  double **loudSNR = dmatrix(0,16,0,COUNT-1);
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
    params[1] = 0.5*pi-theta;
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
     params[1] = theta;//0.5*pi - theta;
     params[2] = phi;
     params[3] = Amp;
     params[4] = iota;
     params[5] = psi;
     params[6] = phase;
     params[7] = fdot;
     params[8] = 11.0/3.0*params[7]*params[7]/params[0];
     */
    cntx++;
    
    N = 32*mult;
    if(f > 0.001) N = 64*mult;
    if(f > 0.01) N = 256*mult;
    if(f > 0.03) N = 512*mult;
    if(f > 0.1) N = 1024*mult;
    
    fonfs = f/fstar;
    
    q = (long)(f*TOBS);
    
    instrument_noise(f, fstar, L, &SAE, &SXYZ);
    
    if(q <= imax)
    {
      SAE  += AEc[q];
      SXYZ += Xc[q];
    }
    
    
    /*  calculate michelson noise  */
    Sm = SXYZ/(4.0*sin(f/fstar)*sin(f/fstar));
    
    Acut = A*sqrt(TOBS/Sm);
    
    M = galactic_binary_bandwidth(LISAorbit->L, LISAorbit->fstar, f, fdot, cos(params[1]), params[3], TOBS, N);
    if(M>N)
    {
      printf("bandwidth=%li, max size=%li\n",M,N);
      exit(0);
    }
    
    XLS = dvector(1,2*M);
    AA  = dvector(1,2*M);
    EE  = dvector(1,2*M);
    
    //FAST_LISA(LISAorbit, TOBS, params, N, M, XLS, AA, EE);
    galactic_binary(LISAorbit, "phase", TOBS, 0, params, 9, XLS, AA, EE, M, 2);

    SNRX = sqrt(Sum(XLS,XLS,M,SXYZ,TOBS));
    SNR  = sqrt(Sum(AA,AA,M,SAE,TOBS)+Sum(EE,EE,M,SAE,TOBS));
    
    free_dvector(XLS,1,2*M);
    free_dvector(AA,1,2*M);
    free_dvector(EE,1,2*M);
    
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
      loudSNR[3][loudSNRcount] = 0.5*pi-theta;  //7:8
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
      theta = 0.5*M_PI - DrawAE[1];
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
//        loudSNR[3][loudSNRcount] = 0.5*pi-theta;  //7:8
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
  printf("fddot measured to 20 percent = %d\n", cnt5);
  
//  printf("*************** X *****************\n");
//  printf("SNR > 7  = %d\n", cnt6);
//  printf("2D sky mapping = %d\n", cnt7);
//  printf("3D sky mapping = %d\n", cnt8);
//  printf("fdot measured to 20 percent = %d\n", cnt9);
//  printf("fddot measured to 20 percent = %d\n", cnt10);
  
  
  
  fclose(afile);
  //fclose(xfile);
  
  
  //sort cell->occupancy so I search from most occupied to least occupied cell
  unsigned long *dumb_sort = lvector(1,COUNT);
  unsigned long *sort = lvector(0,COUNT-1);
  
  indexx(COUNT, loudSNR[0]-1, dumb_sort);
  for(i=0; i<COUNT; i++)
  {
    j = COUNT - i;
    sort[i] = dumb_sort[j]-1;
  }
  
  FILE *outfile = fopen("Brightest.dat","w");
  for(i=0; i<COUNT; i++)
  {
    j = sort[i];
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
  
  d = 9;
  
  epsilon = 1.0e-5;
  
  // Plus and minus parameters:
  
  ParamsP = dvector(0, d-1);
  ParamsM = dvector(0, d-1);
  
  // Plus and minus templates for each channel:
  
  XP = dvector(1, 2*M);
  templateP1 = dvector(1, 2*M);
  templateP2 = dvector(1, 2*M);
  
  XM = dvector(1, 2*M);
  templateM1 = dvector(1, 2*M);
  templateM2 = dvector(1, 2*M);
  
  // Differences for each channel
  
  XD = dmatrix(0, d-1, 1, 2*M);
  temD1 = dmatrix(0, d-1, 1, 2*M);
  temD2 = dmatrix(0, d-1, 1, 2*M);
  
  // The individual Fisher matrices:
  
  Fisher1 = dmatrix(0, d-1, 0, d-1);
  Fisher2 = dmatrix(0, d-1, 0, d-1);
  Cov1 = dmatrix(0, d-1, 0, d-1);
  Cov2 = dmatrix(0, d-1, 0, d-1);
  
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
    
    //       printf("Constructing plus template");
    
    //FAST_LISA(orbit, TOBS, ParamsP, N, M, XP, templateP1, templateP2);
    galactic_binary(orbit, "phase", TOBS, 0, ParamsP, 9, XP, templateP1, templateP2, M, 2);

    //       printf("Constructing minus template");
    
    //FAST_LISA(orbit, TOBS, ParamsM, N, M, XM, templateM1, templateM2);
    galactic_binary(orbit, "phase", TOBS, 0, ParamsM, 9, XP, templateM1, templateM2, M, 2);

    for (j = 1 ; j <= 2*M ; j++)
    {
      XD[i][j] = XP[j] - XM[j];
      temD1[i][j] = templateP1[j] - templateM1[j];
      temD2[i][j] = templateP2[j] - templateM2[j];
      //if(i==7)printf("%i %g %g %g\n",j,temD1[i][j] , templateP1[j] , templateM1[j]);
    }
    
  }
  
  free_dvector(XP, 1, 2*M);
  free_dvector(XM, 1, 2*M);
  free_dvector(templateP1, 1, 2*M);
  free_dvector(templateP2, 1, 2*M);
  free_dvector(templateM1, 1, 2*M);
  free_dvector(templateM2, 1, 2*M);
  
  // Calculate Fisher1 (using just one channel)
  
  for (i = 0 ; i < d ; i++)
  {
    
    for (j = i ; j < d ; j++)
    {
      Fisher1[i][j] =  Sum(XD[i],XD[j],M,SXYZ,TOBS);
    }
  }
  
  free_dmatrix(XD, 0, d-1, 1, 2*M);
  
  
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
  
  free_dmatrix(temD1, 0, d-1, 1, 2*M);
  free_dmatrix(temD2, 0, d-1, 1, 2*M);
  
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
  
  /*
   nA = d;
   lwork = 4*d;
   IPIV = ivector(0,(d-1));
   work = dvector(0,lwork-1);
   A = dvector(0,(d*d-1));
   W = dvector(0,(d-1));
   
   for (i = 0; i < d; i++)
   {
   for (j = 0; j < d; j++)
   {
   A[i*d+j] = Fisher2[i][j];
   }
   }
   
   
   nA = d; lwork = 4*d;
   
   dpotrf_("U", &nA, A, &nA, &info);
   dpotri_("U", &nA, A, &nA, &info);
   
   
   for (i = 0; i < d; i++)
   {
   for (j = 0; j < d; j++)
   {
   Cov2[i][j] = A[i*d+j];
   }
   }
   */
  double **evector = dmatrix(0,d-1,0,d-1);
  double *evalue = dvector(0,d-1);
  matrix_eigenstuff(Fisher2,evector,evalue,d);
  for (i = 0; i < d; i++)for (j = 0; j < d; j++) Cov2[i][j] = Fisher2[i][j];
  
  //get perturbed params (ignoring fddot)
  
  //are evectors normalized?
//  double df = 0.0;
//  for (i = 0; i < 8; i++)
//  {
//    double norm = 0.0;
//    for (j = 0; j < 8; j++)
//    {
//      norm += evector[i][j]*evector[i][j];
//    }
////    printf("f-component[%i]: %g * %g / %g\n",i, u[i], evector[0][i], sqrt(evalue[i]) );
//    df+=evector[0][i]/sqrt(evalue[i]);
//    //printf("|evector[%i]| = %g\n",i,sqrt(norm));
//  }
//  printf("df = %g\n",df);
  
  for (i = 0; i < 8; i++)
  {
    DrawAE[i] = 0.0;
    for (j = 0; j < 8; j++)
    {
      DrawAE[i] += evector[i][j]*u[j]/sqrt(evalue[j])/sqrt(8.);
    }
    if(i==0)DrawAE[i] = DrawAE[i]/TOBS + Params[i];
    else if(i==3)DrawAE[i] = exp(DrawAE[i] + log(Params[i]));
    else if(i==7)DrawAE[i] = DrawAE[i]/TOBS/TOBS + Params[i];
    else DrawAE[i] = DrawAE[i] + Params[i];
    
  }
  
  free_dvector(evalue,0,d-1);
  free_dmatrix(evector,0,d-1,0,d-1);
  
  
  
  //get Fisher estimated errors
  for (i = 0; i < d; i++) SigmaAE[i] = sqrt(Cov2[i][i]);
  
  SigmaAE[9] = 2.0*pi*sin(Params[1])*sqrt(Cov2[1][1]*Cov2[2][2] - Cov2[2][1]*Cov2[2][1]);
  
  
//  for (j = 0; j < 2; j++)
//  {
//    if(SigmaAE[d-1] > 0.5)
//    {
//      SigmaAE[d-1] = 1.0e10;
//      d -= 1;
//      /*
//       for (i = 0; i < d; i++)
//       {
//       for (j = 0; j < d; j++)
//       {
//       Cov2[i][j] = Fisher2[i][j];
//       }
//       }
//       
//       nA = d; lwork = 4*d;
//       
//       dpotrf_("U", &nA, A, &nA, &info);
//       dpotri_("U", &nA, A, &nA, &info);
//       
//       for (i = 0; i < d; i++) SigmaAE[i] = sqrt(A[i*d+i]);
//       */
//      evector = dmatrix(0,d-1,0,d-1);
//      evalue = dvector(0,d-1);
//      matrix_eigenstuff(Fisher2,evector,evalue,d);
//      for (i = 0; i < d; i++)for (j = 0; j < d; j++) Cov2[i][j] = Fisher2[i][j];
//      free_dvector(evalue,0,d-1);
//      free_dmatrix(evector,0,d-1,0,d-1);
//      
//      for (i = 0; i < d; i++) SigmaAE[i] = sqrt(Cov2[i][i]);
//      
//    }
//  }
  
  // printf("%d %f %f\n", d, SigmaAE[7], SigmaAE[8]);
  
  d = 9;
  
  /*
   for (i = 0; i < d; i++)
   {
   for (j = 0; j < d; j++)
   {
   A[i*d+j] = Fisher1[i][j];
   }
   }
   
   nA = d;
   lwork = 4*d;
   
   dpotrf_("U", &nA, A, &nA, &info);
   dpotri_("U", &nA, A, &nA, &info);
   
   for (i = 0; i < d; i++)
   {
   for (j = 0; j < d; j++)
   {
   Cov1[i][j] = A[i*d+j];
   }
   }
   
   for (i = 0; i < d; i++) SigmaX[i] = sqrt(A[i*d+i]);
   SigmaX[9] = 2.0*pi*sin(Params[1])*sqrt(Cov1[1][1]*Cov1[2][2] - Cov1[2][1]*Cov1[2][1]);
   
   for (j = 0; j < 2; j++)
   {
   if(SigmaX[d-1] > 0.5)
   {
   SigmaX[d-1] = 1.0e10;
   d -= 1;
   for (i = 0; i < d; i++)
   {
   for (j = 0; j < d; j++)
   {
   A[i*d+j] = Fisher1[i][j];
   }
   }
   
   nA = d; lwork = 4*d;
   
   dpotrf_("U", &nA, A, &nA, &info);
   dpotri_("U", &nA, A, &nA, &info);
   
   for (i = 0; i < d; i++) SigmaX[i] = sqrt(A[i*d+i]);
   
   }
   }*/
  
  evector = dmatrix(0,d-1,0,d-1);
  evalue = dvector(0,d-1);
  matrix_eigenstuff(Fisher1,evector,evalue,d);
  for (i = 0; i < d; i++)for (j = 0; j < d; j++) Cov1[i][j] = Fisher1[i][j];
  free_dvector(evalue,0,d-1);
  free_dmatrix(evector,0,d-1,0,d-1);
  
  for (i = 0; i < d; i++) SigmaX[i] = sqrt(Cov1[i][i]);
  
  SigmaX[9] = 2.0*pi*sin(Params[1])*sqrt(Cov1[1][1]*Cov1[2][2] - Cov1[2][1]*Cov1[2][1]);
  
  
  for (j = 0; j < 2; j++)
  {
    if(SigmaX[d-1] > 0.5)
    {
      SigmaX[d-1] = 1.0e10;
      d -= 1;
      /*
       for (i = 0; i < d; i++)
       {
       for (j = 0; j < d; j++)
       {
       Cov2[i][j] = Fisher2[i][j];
       }
       }
       
       nA = d; lwork = 4*d;
       
       dpotrf_("U", &nA, A, &nA, &info);
       dpotri_("U", &nA, A, &nA, &info);
       
       for (i = 0; i < d; i++) SigmaAE[i] = sqrt(A[i*d+i]);
       */
      evector = dmatrix(0,d-1,0,d-1);
      evalue = dvector(0,d-1);
      matrix_eigenstuff(Fisher1,evector,evalue,d);
      for (i = 0; i < d; i++)for (j = 0; j < d; j++) Cov1[i][j] = Fisher1[i][j];
      free_dvector(evalue,0,d-1);
      free_dmatrix(evector,0,d-1,0,d-1);
      
      for (i = 0; i < d; i++) SigmaX[i] = sqrt(Cov2[i][i]);
      
    }
  }
  
  
  //  free_dvector(work, 0, lwork-1);
  //  free_dvector(A, 0,(d*d-1));
  //  free_dvector(W, 0,(d-1));
  //  free_ivector(IPIV, 0,(d-1));
  
  d = 9;
  
  free_dvector(ParamsP, 0, d-1);
  free_dvector(ParamsM, 0, d-1);
  
  free_dmatrix(Fisher1, 0, d-1, 0, d-1);
  free_dmatrix(Fisher2, 0, d-1, 0, d-1);
  
  free_dmatrix(Cov1, 0, d-1, 0, d-1);
  free_dmatrix(Cov2, 0, d-1, 0, d-1);
  
}



/*************************************************************************************/
/*                                                                                   */
/*                                    Subroutines                                    */
/*                                                                                   */
/*************************************************************************************/



/*************************************************************************/
/*        Rigid approximation position of each LISA spacecraft           */
/*************************************************************************/
//void spacecraft(double t, double *x, double *y, double *z)
//{
//
//	double alpha;
//	double beta1, beta2, beta3;
//	double sa, sb, ca, cb;
//
//	alpha = 2.*pi*fm*t + kappa;
//
//	beta1 = 0. + lambda;
//	beta2 = 2.*pi/3. + lambda;
//	beta3 = 4.*pi/3. + lambda;
//
//	sa = sin(alpha);
//	ca = cos(alpha);
//
//
//	sb = sin(beta1);
//	cb = cos(beta1);
//	x[1] = AU*ca + AU*ec*(sa*ca*sb - (1. + sa*sa)*cb);
//	y[1] = AU*sa + AU*ec*(sa*ca*cb - (1. + ca*ca)*sb);
//	z[1] = -sq3*AU*ec*(ca*cb + sa*sb);
//
//
//	sb = sin(beta2);
//	cb = cos(beta2);
//	x[2] = AU*ca + AU*ec*(sa*ca*sb - (1. + sa*sa)*cb);
//	y[2] = AU*sa + AU*ec*(sa*ca*cb - (1. + ca*ca)*sb);
//	z[2] = -sq3*AU*ec*(ca*cb + sa*sb);
//
//	sb = sin(beta3);
//	cb = cos(beta3);
//	x[3] = AU*ca + AU*ec*(sa*ca*sb - (1. + sa*sa)*cb);
//	y[3] = AU*sa + AU*ec*(sa*ca*cb - (1. + ca*ca)*sb);
//	z[3] = -sq3*AU*ec*(ca*cb + sa*sb);
//
//}



void indexx(unsigned long n, double *arr, unsigned long *indx)
{
  unsigned long i,indxt,ir=n,j,k,l=1;
  int jstack=0,*istack,M=7;
  double a;
  double tempr;
  
  istack=ivector(1,NSTACK);
  for (j=1;j<=n;j++) indx[j]=j;
  for (;;) {
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
        indxt=indx[j];
        a=arr[indxt];
        for (i=j-1;i>=l;i--) {
          if (arr[indx[i]] <= a) break;
          indx[i+1]=indx[i];
        }
        indx[i+1]=indxt;
      }
      if (jstack == 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
    } else {
      k=(l+ir) >> 1;
      //SWAP(indx[k],indx[l+1]);
      tempr=indx[k];indx[k]=indx[l+1];indx[l+1]=tempr;
      if (arr[indx[l]] > arr[indx[ir]]) {
        //SWAP(indx[l],indx[ir]);
        tempr=indx[l],indx[l]=indx[ir];indx[ir]=tempr;
      }
      if (arr[indx[l+1]] > arr[indx[ir]]) {
        //SWAP(indx[l+1],indx[ir]);
        tempr=indx[l+1],indx[l+1]=indx[ir];indx[ir]=tempr;
        
      }
      if (arr[indx[l]] > arr[indx[l+1]]) {
        //SWAP(indx[l],indx[l+1]);
        tempr=indx[l],indx[l]=indx[l+1];indx[l+1]=tempr;
      }
      i=l+1;
      j=ir;
      indxt=indx[l+1];
      a=arr[indxt];
      for (;;) {
        do i++; while (arr[indx[i]] < a);
        do j--; while (arr[indx[j]] > a);
        if (j < i) break;
        //SWAP(indx[i],indx[j]);
        tempr=indx[i];indx[i]=indx[j];indx[j]=tempr;
      }
      indx[l+1]=indx[j];
      indx[j]=indxt;
      jstack += 2;
      //if (jstack > NSTACK) nrerror("NSTACK too small in indexx.");
      if (ir-i+1 >= j-l) {
        istack[jstack]=ir;
        istack[jstack-1]=i;
        ir=j-1;
      } else {
        istack[jstack]=j-1;
        istack[jstack-1]=l;
        l=i;
      }
    }
  }
  free_ivector(istack,1,NSTACK);
}



