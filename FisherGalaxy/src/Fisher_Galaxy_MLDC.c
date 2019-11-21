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
#include "arrays.h"
#include "Constants.h"
#include "Detector.h"

struct lisa_orbit{
  int N;
  
  double L;
  double fstar;
  
  double *t;
  double **x;
  double **y;
  double **z;
  double **dx;
  double **dy;
  double **dz;
};
struct lisa_orbit orbit;

//void spacecraft(double t, double *x, double *y, double *z);
void spacecraft(struct lisa_orbit *orbit, double tint, double *xint, double *yint, double *zint);
void initialize_orbit(char OrbitFile[], struct lisa_orbit *orbit);
void LISA_spline(double *x, double *y, int n, double yp1, double ypn, double *y2);
void LISA_splint(double *xa, double *ya, double *y2a, int n, double x, double *y);

void convolve(long N, double *a, long M, double *b, double *cn);
void XYZ(double L, double fstar, double ***d, double f0, long q, long M, double *XLS, double *ALS, double *ELS);
void FAST_LISA(struct lisa_orbit *orbit, double *params, long N, long M, double *XLS, double *ALS, double *ELS);
void instrument_noise(double f, double fstar, double L, double *SAE, double *SXYZ);
double AE_confusion_noise(double f, double *carray, double *farray, int Nfit);
double X_confusion_noise(double f, double *carrayX, double *farrayX, int NfitX);
double Sum(double *AA, double *EE, long M, double SN);
void FISHER(struct lisa_orbit *orbit, double *params, long N, long M, double SXYZ, double SAE, double *SigmaX, double *SigmaAE);
void dfour1(double data[], unsigned long nn, int isign);
void KILL(char*);

int main(int argc,char **argv)
{
  char Gfile[50];
  double *params;
  double *XfLS, *YfLS, *ZfLS;
  double *XfSL, *YfSL, *ZfSL;
  double *XLS, *YLS, *ZLS;
  double *AA, *EE;
  double *carray, *farray;
  double *SigmaX, *SigmaAE;
  double Aconf;
  int Nfit;
  double *carrayX, *farrayX;
  double Xconf;
  int NfitX;
  double fonfs, Sn, Sm, Acut;
  double SNRX, SNRA, SNRE, SNR, SNRin;
  double SAE, SXYZ;
  double SNR_av, SNRX_av;
  double SNR_max, SNRX_max, SNR_min, SNRX_min;
  long M, N, q;
  long rSeed;
  long i, j, k, cnt, mult, ii, Nav;
  double m1, m2, Mc, Porb, Porb_dot, l, b, DL;
  double fix, Mfac, Afac, Amp;
  double x, y ,z, r;
  double x_ec, y_ec, z_ec;
  double f, fdot, A, theta, phi, iota, psi, phase, fddot;
  int cntx, cnt1, cnt2, cnt3, cnt4, cnt5, cnt6, cnt7, cnt8, cnt9, cnt10;
  int imin, imax;
  double *Xc, *AEc, *Xnoise, *AEnoise;
  FILE* vfile;
  FILE* afile;
  FILE* xfile;
  FILE* sky;
  FILE* amps;
  
  
  if(argc != 6) KILL("Fisher_Galaxy detections.dat confusion.dat sigmasX.dat sigmasAE.dat Orbit.dat\n");
  
  vfile = fopen(argv[1],"r");
  
  //Data structure for interpolating orbits from file
  struct lisa_orbit *LISAorbit;
  LISAorbit = &orbit;
  
  //Set up orbit structure (allocate memory, read file, cubic spline)
  sprintf(Gfile,"%s",argv[5]);
  initialize_orbit(Gfile, LISAorbit);
  
  double L    = LISAorbit->L;
  double fstar= LISAorbit->fstar;
  
  imax = (int)ceil(1.0e-2*T);
  imin = (int)floor(1.0e-4*T);
  
  Xc = dvector(imin,imax);  AEc = dvector(imin,imax);
  Xnoise = dvector(imin,imax);  AEnoise = dvector(imin,imax);
  
  amps = fopen(argv[2],"r");
  for(i=imin; i<= imax; i++)
  {
    fscanf(amps,"%lg%lg%lg%lg%lg\n", &f, &Xnoise[i], &Xc[i], &AEnoise[i], &AEc[i]);
    instrument_noise(f, fstar, L, &SAE, &SXYZ);
    // printf("%e %e %e\n", f, SXYZ+Xc[i], SAE+AEc[i]);
  }
  fclose(amps);
  
  params = dvector(0,8);
  SigmaX = dvector(0,8);
  SigmaAE = dvector(0,8);
  
  
  fix = pi/180.0;  /* one degree in radians */
  
  
  if((T/year) <= 8.0) mult = 8;
  if((T/year) <= 4.0) mult = 4;
  if((T/year) <= 2.0) mult = 2;
  if((T/year) <= 1.0) mult = 1;
  
  xfile = fopen(argv[3], "w");
  afile = fopen(argv[4], "w");
  
  
  cntx = cnt1 = cnt2 = cnt3 = cnt4 = cnt5 = cnt6 = cnt7 = cnt8 = cnt9 = cnt10 = 0;
  
  while ( !feof(vfile) )
  {
    
    fscanf(vfile, "%lg%lg%lg%lg%lg%lg%lg%lg\n", &f, &fdot, &theta, &phi, &Amp, &iota, &psi, &phase);
    
    params[0] = f/T;
    params[1] = 0.5*pi-theta;
    params[2] = phi;
    params[3] = Amp;
    params[4] = iota;
    params[5] = psi;
    params[6] = phase;
    params[7] = fdot;
    params[8] = 11.0/3.0*fdot*fdot/f;
    
    cntx++;
    
    if(cntx%1000 == 0) printf("%d processed\n", cntx);
    
    N = 32*mult;
    if(f > 0.001) N = 64*mult;
    if(f > 0.01) N = 256*mult;
    if(f > 0.03) N = 512*mult;
    if(f > 0.1) N = 1024*mult;
    
    fonfs = f/fstar;
    
    q = (long)(f*T);
    
    instrument_noise(f, fstar, L, &SAE, &SXYZ);
    
    if(q <= imax)
    {
      SAE += AEc[q];
      SXYZ += Xc[q];
    }
    
    
    /*  calculate michelson noise  */
    Sm = SXYZ/(4.0*sin(f/fstar)*sin(f/fstar));
    
    Acut = A*sqrt(T/Sm);
    
    M = (long)(pow(2.0,(rint(log(Acut)/log(2.0))+1.0)));
    
    if(M < N) M = N;
    if(M > 8192) M = 8192;
    N=M;
    
    XLS = dvector(1,2*M);
    AA = dvector(1,2*M);   EE = dvector(1,2*M);
    
    FAST_LISA(LISAorbit, params, N, M, XLS, AA, EE);
    
    SNRX = sqrt(Sum(XLS,XLS,M,SXYZ));
    SNR = sqrt(Sum(AA,AA,M,SAE)+Sum(EE,EE,M,SAE));
    
    // printf("SNR_X = %f  SNR_in = %f SNR_AE = %f\n", SNRX, SNRin, SNR);
    
    free_dvector(XLS,1,2*M);
    free_dvector(AA,1,2*M);  free_dvector(EE,1,2*M);
    
    if(SNR > 7.0)
    {
      
      cnt1++;
      
      FISHER(LISAorbit, params, N, M, SXYZ, SAE, SigmaX, SigmaAE);
      
      fprintf(afile, "%e %e %e ", params[0], params[7], params[8]);
      for(i = 1; i < 7; i++)
      {
        fprintf(afile, "%e ", params[i]);
      }
      
      fdot = SigmaAE[7]/(fabs(params[7]*T*T));  // fraction uncertainity in fdot
      fddot = SigmaAE[8]/(fabs(params[8]*T*T*T));  // fraction uncertainity in fddot
      fprintf(afile, "%e ", SigmaAE[0]/T);  // convert back to Hz
      fprintf(afile, "%e %e ", fdot, fddot);
      for(i = 1; i < 7; i++)
      {
        fprintf(afile, "%e ", SigmaAE[i]);
      }
      fprintf(afile, "%f\n", SNR);
      
      if((SigmaAE[1] < fix) && (SigmaAE[2] < fix))
      {
        cnt2++;   // 2D mapping
        if((fdot < 0.1) && (SigmaAE[3] < 0.1))
        {
          cnt3++;  // 3D mapping
        }
      }
      
      if(fdot < 0.2)
      {
        cnt4++;  // measure fdot to 20%
      }
      
      if(fddot < 0.2)
      {
        cnt5++;  // measure fddot to 20%
        printf("%e %f\n", f, SNR);
      }
      
      
      
      if(SNRX > 7.0)
      {
        cnt6++;  /* SNR > 7 with just Michelson */
        
        fprintf(xfile, "%e %e %e ", params[0], params[7], params[8]);
        for(i = 1; i < 7; i++)
        {
          fprintf(xfile, "%e ", params[i]);
        }
        
        fdot = SigmaX[7]/(fabs(params[7]*T*T));  // fraction uncertainity in fdot
        fddot = SigmaX[8]/(fabs(params[8]*T*T*T));  // fraction uncertainity in fddot
        fprintf(xfile, "%e ", SigmaX[0]/T);  // convert back to Hz
        fprintf(xfile, "%e %e ", fdot, fddot);
        for(i = 1; i < 7; i++)
        {
          fprintf(xfile, "%e ", SigmaX[i]);
        }
        fprintf(xfile, "%f\n", SNRX);
        
        if((SigmaX[1] < fix) && (SigmaX[2] < fix))
        {
          cnt7++;   // 2D mapping
          if((fdot < 0.1) && (SigmaAE[3] < 0.1))
          {
            cnt8++;  // 3D mapping
          }
        }
        
        if(fdot < 0.2)
        {
          cnt9++;  // measure fdot to 20%
        }
        
        if(fddot < 0.2)
        {
          cnt10++;  // measure fddot to 20%
          // printf("%e %f\n", f, SNRX);
        }
        
      }
      
      
      
      
    }
    
    
  }
  
  printf("*************** AE *****************\n");
  printf("SNR > 7  = %d\n", cnt1);
  printf("2D sky mapping = %d\n", cnt2);
  printf("3D sky mapping = %d\n", cnt3);
  printf("fdot measured to 20 percent = %d\n", cnt4);
  printf("fddot measured to 20 percent = %d\n", cnt5);
  
  printf("*************** X *****************\n");
  printf("SNR > 7  = %d\n", cnt6);
  printf("2D sky mapping = %d\n", cnt7);
  printf("3D sky mapping = %d\n", cnt8);
  printf("fdot measured to 20 percent = %d\n", cnt9);
  printf("fddot measured to 20 percent = %d\n", cnt10);
  
  fclose(afile);
  fclose(xfile);
  
  return;
  
}


double Sum(double *AA, double *EE, long M, double SN)
{
  long i;
  double sm;
  
  sm = 0.0;
  
  for(i=1; i<=M; i++)
  {
    sm += (AA[2*i-1]*EE[2*i-1]+AA[2*i]*EE[2*i]);
  }
  
  sm *= 4.0*T/SN;
  
  return sm;
  
}

void FISHER(struct lisa_orbit *orbit, double *Params, long N, long M, double SXYZ, double SAE, double *SigmaX, double *SigmaAE)
{
  
  int i, j, k, flag;
  int ii, jj, ibad;
  double *ParamsP, *ParamsM;
  double *templateP1, *templateP2;
  double *templateM1, *templateM2;
  double *XP, *XM;
  double **XD, **temD1, **temD2;
  double *templateA, *templateB;
  double **Fisher1, **Fisher2;
  double epsilon, max;
  double fdot, fddot;
  int d;
  
  double *A, *B;
  double *W;
  int *IPIV;
  int lwork;
  int nA;
  double *work;
  int info;
  
  d = 9;
  
  epsilon = 1.0e-6;
  
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
  
  for (i = 0 ; i < d ; i++)
  {
    for (j = 0 ; j < d ; j++)
    {
      ParamsP[j] = Params[j];
      ParamsM[j] = Params[j];
    }
    
    if(i==0)
    {
      ParamsP[0] *= T;
      ParamsM[0] *= T;
    }
    if(i==3)
    {
      ParamsP[3] = log(ParamsP[3]);
      ParamsM[3] = log(ParamsM[3]);
    }
    if(i==7)
    {
      ParamsP[7] *= T*T;
      ParamsM[7] *= T*T;
    }
    if(i==8)
    {
      ParamsP[8] *= T*T*T;
      ParamsM[8] *= T*T*T;
    }
    
    
    ParamsP[i] += epsilon;
    ParamsM[i] -= epsilon;
    
    if(i==0)
    {
      ParamsP[0] /= T;
      ParamsM[0] /= T;
    }
    if(i==3)
    {
      ParamsP[3] = exp(ParamsP[3]);
      ParamsM[3] = exp(ParamsM[3]);
    }
    if(i==7)
    {
      ParamsP[7] /= T*T;
      ParamsM[7] /= T*T;
    }
    if(i==8)
    {
      ParamsP[8] /= T*T*T;
      ParamsM[8] /= T*T*T;
    }
    
    
    //       printf("Constructing plus template");
    
    FAST_LISA(orbit, ParamsP, N, M, XP, templateP1, templateP2);
    
    //       printf("Constructing minus template");
    
    FAST_LISA(orbit, ParamsM, N, M, XM, templateM1, templateM2);
    
    for (j = 1 ; j <= 2*M ; j++)
    {
      XD[i][j] = XP[j] - XM[j];
      temD1[i][j] = templateP1[j] - templateM1[j];
      temD2[i][j] = templateP2[j] - templateM2[j];
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
      Fisher1[i][j] =  Sum(XD[i],XD[j],M,SXYZ);
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
      Fisher2[i][j] =  Sum(temD1[i],temD1[j],M,SAE)+Sum(temD2[i],temD2[j],M,SAE);
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
      Fisher2[i][j] /= (4.0*epsilon*epsilon);
  }
  
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
  
  for (i = 0; i < d; i++) SigmaAE[i] = sqrt(A[i*d+i]);
  
  for (j = 0; j < 2; j++)
  {
    if(SigmaAE[d-1] > 0.5)
    {
      SigmaAE[d-1] = 1.0e10;
      d -= 1;
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
      
      for (i = 0; i < d; i++) SigmaAE[i] = sqrt(A[i*d+i]);
      
    }
  }
  
  // printf("%d %f %f\n", d, SigmaAE[7], SigmaAE[8]);
  
  d = 9;
  
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
  
  for (i = 0; i < d; i++) SigmaX[i] = sqrt(A[i*d+i]);
  
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
  }
  
  
  
  free_dvector(work, 0, lwork-1);
  free_dvector(A, 0,(d*d-1));
  free_dvector(W, 0,(d-1));
  free_ivector(IPIV, 0,(d-1));
  
  d = 9;
  
  free_dvector(ParamsP, 0, d-1);
  free_dvector(ParamsM, 0, d-1);
  
  free_dmatrix(Fisher1, 0, d-1, 0, d-1);
  free_dmatrix(Fisher2, 0, d-1, 0, d-1);
  
}


void FAST_LISA(struct lisa_orbit *orbit, double *params, long N, long M, double *XLS, double *ALS, double *ELS)
{
  
  /*   Indicies   */
  int i,j,n,m;
  /*   Carrier frequency bin  */
  long q;
  /*   Gravitational Wave basis vectors   */
  double *u,*v,*k;
  /*   Polarization basis tensors   */
  double **eplus, **ecross;
  /*   Spacecraft position and separation vector   */
  double *x, *y, *z;
  double *r12, *r13, *r21, *r23, *r31, *r32;
  /*   Dot products   */
  double *kdotx, **kdotr;
  /*   Convenient quantities   */
  double **dplus, **dcross;
  /*   GW Source data   */
  double Mc, theta, phi, psi, D, iota, Amp, Aplus, Across, f0, fdot, fddot, phio;
  double costh, sinth, cosph, sinph, cosi, cosps, sinps;
  /*   Time and distance variables   */
  double *xi, t;
  /*   Gravitational wave frequency & ratio of f and transfer frequency f*  */
  double L, fstar, *f, *fonfs;
  /*   LISA response to slow terms (Real & Imaginary pieces)   */
  //Static quantities (Re and Im)
  double DPr, DPi, DCr, DCi;
  //Time varrying quantities (Re & Im) broken up into convenient segments
  double **TR, **TI, tran1r, tran1i, tran2r, tran2i;
  //Miscellaneous constants used to speed up calculations
  double arg1, arg2, sinc, df;
  /*   Fourier coefficients before FFT and after convolution  */
  //Time series of slowly evolving terms at each vertex
  double *data12, *data13, *data21, *data23, *data31, *data32;
  //Fourier coefficients of slowly evolving terms (numerical)
  double *a12, *a13, *a21, *a23, *a31, *a32;
  //Fourier coefficients of rapidly evolving terms (analytical)
  double *b;
  //Fourier coefficients of entire response (convolution)
  double *c12, *c13, *c21, *c23, *c31, *c32;
  //Package cij's into proper form for TDI subroutines
  double ***d;
  /*   Miscellaneous  */
  double xm, fstep, power, aevol;
  
  
  /*   Allocating Arrays   */
  
  u = dvector(1,3); v = dvector(1,3); k = dvector(1,3);
  
  kdotx = dvector(1,3); kdotr = dmatrix(1,3,1,3);
  
  xi = dvector(1,3);
  
  f = dvector(1,3);
  
  fonfs = dvector(1,3);
  
  eplus  = dmatrix(1,3,1,3); ecross = dmatrix(1,3,1,3);
  
  dplus  = dmatrix(1,3,1,3); dcross = dmatrix(1,3,1,3);
  
  x = dvector(1,3); y = dvector(1,3); z = dvector(1,3);
  
  r12 = dvector(1,3); r21 = dvector(1,3); r31 = dvector(1,3);
  r13 = dvector(1,3); r23 = dvector(1,3); r32 = dvector(1,3);
  
  TR = dmatrix(1,3,1,3); TI = dmatrix(1,3,1,3);
  
  data12 = dvector(1,2*N); data21 = dvector(1,2*N); data31 = dvector(1,2*N);
  data13 = dvector(1,2*N); data23 = dvector(1,2*N); data32 = dvector(1,2*N);
  
  a12 = dvector(1,2*N+2); a21 = dvector(1,2*N+2); a31 = dvector(1,2*N+2);
  a13 = dvector(1,2*N+2); a23 = dvector(1,2*N+2); a32 = dvector(1,2*N+2);
  
  d = d3tensor(1,3,1,3,1,2*M);
  
  /*   Gravitational Wave source parameters   */
  
  f0 = params[0];
  costh =cos(params[1]);
  phi = params[2];
  Amp = params[3];
  cosi = cos(params[4]);
  psi = params[5];
  phio = params[6];
  fdot = params[7];
  fddot = params[8];
  
  //Calculate cos and sin of sky position, inclination, polarization
  sinth = sqrt(1.0-costh*costh);
  cosph = cos(phi);     sinph = sin(phi);
  cosps = cos(2.*psi);  sinps = sin(2.*psi);
  
  //Calculate GW polarization amplitudes
  Aplus = Amp*(1.+cosi*cosi);
  Across = -2.0*Amp*cosi;
  
  //Calculate carrier frequency bin
  q = (long)(f0*T);
  
  df = 2.0*pi*(((double)q)/T);
  
  //Calculate constant pieces of transfer functions
  DPr = Aplus*cosps;
  DPi = -Across*sinps;
  DCr = -Aplus*sinps;
  DCi = -Across*cosps;
  
  
  /*   Tensor construction for buildingslowly evolving LISA response   */
  //Gravitational Wave source basis vectors
  u[1] =  costh*cosph;  u[2] =  costh*sinph;  u[3] = -sinth;
  v[1] =  sinph;        v[2] = -cosph;        v[3] =  0.;
  k[1] = -sinth*cosph;  k[2] = -sinth*sinph;  k[3] = -costh;
  
  //GW polarization basis tensors
  for(i=1;i<=3;i++)
  {
    for(j=1;j<=3;j++)
    {
      eplus[i][j]  = u[i]*u[j] - v[i]*v[j];
      ecross[i][j] = u[i]*v[j] + v[i]*u[j];
      
    }
  }
  
  
  /*****************************   Main Loop   **********************************/
  L = orbit->L;
  fstar = orbit->fstar;
  for(n=1; n<=N; n++)
  {
    //First time sample must be at t=0 for phasing
    t = T*(double)(n-1)/(double)N;
    
    //Calculate position of each spacecraft at time t
    spacecraft(orbit, t, x, y, z);
    
    for(i=1; i<=3; i++)
    {
      kdotx[i] = (x[i]*k[1]+y[i]*k[2]+z[i]*k[3])/clight;
      //Wave arrival time at spacecraft i
      xi[i]    = t - kdotx[i];
      //First order approximation to frequency at spacecraft i
      f[i]     = f0 + fdot*xi[i]+0.5*fddot*xi[i]*xi[i];
      //Ratio of true frequencto transfer frequency
      fonfs[i] = f[i]/orbit->fstar;
    }
    
    
    //Unit separation vector from spacecrafts i to j
    r12[1] = (x[2] - x[1])/L;   r13[1] = (x[3] - x[1])/L;   r23[1] = (x[3] - x[2])/L;
    r12[2] = (y[2] - y[1])/L;   r13[2] = (y[3] - y[1])/L;   r23[2] = (y[3] - y[2])/L;
    r12[3] = (z[2] - z[1])/L;   r13[3] = (z[3] - z[1])/L;   r23[3] = (z[3] - z[2])/L;
    //Make use of symmetry
    for(i=1; i<=3; i++)
    {
      r21[i] = -r12[i];
      r31[i] = -r13[i];
      r32[i] = -r23[i];
    }
    
    //Zero arrays to be summed
    dplus[1][2] = dplus[1][3] = dplus[2][1] = dplus[2][3] = dplus[3][1] = dplus[3][2] = 0.;
    dcross[1][2] = dcross[1][3] = dcross[2][1] = dcross[2][3] = dcross[3][1] = dcross[3][2] = 0.;
    //Convenient quantities d+ & dx
    for(i=1; i<=3; i++)
    {
      for(j=1; j<=3; j++)
      {
        dplus[1][2]  += r12[i]*r12[j]*eplus[i][j];   dcross[1][2] += r12[i]*r12[j]*ecross[i][j];
        dplus[2][3]  += r23[i]*r23[j]*eplus[i][j];   dcross[2][3] += r23[i]*r23[j]*ecross[i][j];
        dplus[1][3]  += r13[i]*r13[j]*eplus[i][j];   dcross[1][3] += r13[i]*r13[j]*ecross[i][j];
      }
    }
    //Makng use of symmetry
    dplus[2][1] = dplus[1][2];  dcross[2][1] = dcross[1][2];
    dplus[3][2] = dplus[2][3];  dcross[3][2] = dcross[2][3];
    dplus[3][1] = dplus[1][3];  dcross[3][1] = dcross[1][3];
    
    //Zero arrays to be summed
    kdotr[1][2] = kdotr[1][3] = kdotr[2][1] = kdotr[2][3] = kdotr[3][1] = kdotr[3][2] = 0.;
    for(i=1; i<=3; i++)
    {
      kdotr[1][2] += k[i]*r12[i];   kdotr[1][3] += k[i]*r13[i];   kdotr[2][3] += k[i]*r23[i];
    }
    //Making use of antisymmetry
    kdotr[2][1] = -kdotr[1][2];
    kdotr[3][1] = -kdotr[1][3];
    kdotr[3][2] = -kdotr[2][3];
    
    //Calculating Transfer function
    for(i=1; i<=3; i++)
    {
      for(j=1; j<=3; j++)
      {
        if(i!=j)
        {
          //Argument of transfer function
          arg1 = 0.5*fonfs[i]*(1 - kdotr[i][j]);
          //Argument of complex exponentials
          arg2 = 2.0*pi*f0*xi[i]+pi*fdot*xi[i]*xi[i]+pi*fddot*xi[i]*xi[i]*xi[i]/3.0 +phio - df*t;
          //Transfer function
          sinc = 0.25*sin(arg1)/arg1;
          //Evolution of amplitude
          aevol = 1.0 + 0.66666666666666666666*fdot/f0*xi[i];
          ///Real and imaginary pieces of time series (no complex exponential)
          tran1r = aevol*(dplus[i][j]*DPr + dcross[i][j]*DCr);
          tran1i = aevol*(dplus[i][j]*DPi + dcross[i][j]*DCi);
          //Real and imaginry components of complex exponential
          tran2r = cos(arg1 + arg2);
          tran2i = sin(arg1 + arg2);
          //Real & Imaginary part of the slowly evolving signal
          TR[i][j] = sinc*(tran1r*tran2r - tran1i*tran2i);
          TI[i][j] = sinc*(tran1r*tran2i + tran1i*tran2r);
        }
      }
    }
    
    //Fill  time series data arrays with slowly evolving signal.
    //dataij corresponds to fractional arm length difference yij
    data12[2*n-1] = TR[1][2];   data21[2*n-1] = TR[2][1];   data31[2*n-1] = TR[3][1];
    data12[2*n]   = TI[1][2];   data21[2*n]   = TI[2][1];   data31[2*n]   = TI[3][1];
    data13[2*n-1] = TR[1][3];   data23[2*n-1] = TR[2][3];   data32[2*n-1] = TR[3][2];
    data13[2*n]   = TI[1][3];   data23[2*n]   = TI[2][3];   data32[2*n]   = TI[3][2];
  }
  
  /*   Numerical Fourier transform of slowly evolving signal   */
  dfour1(data12, N, -1);  dfour1(data21, N, -1);  dfour1(data31, N, -1);
  dfour1(data13, N, -1);  dfour1(data23, N, -1);  dfour1(data32, N, -1);
  //Unpack arrays from dfour1.c and normalize
  for(i=1; i<=N; i++)
  {
    a12[i] = data12[N+i]/(double)N;  a21[i] = data21[N+i]/(double)N;  a31[i] = data31[N+i]/(double)N;
    a12[i+N] = data12[i]/(double)N;  a21[i+N] = data21[i]/(double)N;  a31[i+N] = data31[i]/(double)N;
    a13[i] = data13[N+i]/(double)N;  a23[i] = data23[N+i]/(double)N;  a32[i] = data32[N+i]/(double)N;
    a13[i+N] = data13[i]/(double)N;  a23[i+N] = data23[i]/(double)N;  a32[i+N] = data32[i]/(double)N;
  }
  
  a12[2*N+1] = data12[N+1]/(double)N;  a21[2*N+1] = data21[N+1]/(double)N;  a31[2*N+1] = data31[N+1]/(double)N;
  a12[2*N+2] = data12[N+2]/(double)N;  a21[2*N+2] = data21[N+2]/(double)N;  a31[2*N+2] = data31[N+2]/(double)N;
  a13[2*N+1] = data13[N+1]/(double)N;  a23[2*N+1] = data23[N+1]/(double)N;  a32[2*N+1] = data32[N+1]/(double)N;
  a13[2*N+2] = data13[N+2]/(double)N;  a23[2*N+2] = data23[N+2]/(double)N;  a32[2*N+2] = data32[N+2]/(double)N;
  
  /*   Renormalize so that the resulting time series is real   */
  for(i=1; i<=2*M; i++)
  {
    d[1][2][i] = 0.5*a12[i];  d[2][1][i] = 0.5*a21[i];  d[3][1][i] = 0.5*a31[i];
    d[1][3][i] = 0.5*a13[i];  d[2][3][i] = 0.5*a23[i];  d[3][2][i] = 0.5*a32[i];
  }
  
  /*   Call subroutines for synthesizing different TDI data channels  */
  /*   X Y Z-Channel   */
  XYZ(L, fstar, d, f0, q, M, XLS, ALS, ELS);
  
  
  /*   Deallocate Arrays   */
  
  free_dvector(u,1,3); free_dvector(v,1,3); free_dvector(k,1,3);
  
  free_dvector(kdotx,1,3); free_dmatrix(kdotr,1,3,1,3);
  
  free_dvector(xi,1,3);
  
  free_dvector(f,1,3);
  
  free_dvector(fonfs,1,3);
  
  free_dmatrix(eplus,1,3,1,3); free_dmatrix(ecross,1,3,1,3);
  
  free_dmatrix(dplus,1,3,1,3); free_dmatrix(dcross,1,3,1,3);
  
  free_dvector(x,1,3); free_dvector(y,1,3); free_dvector(z,1,3);
  
  free_dvector(r12,1,3); free_dvector(r21,1,3); free_dvector(r31,1,3);
  free_dvector(r13,1,3); free_dvector(r23,1,3); free_dvector(r32,1,3);
  
  free_dmatrix(TR,1,3,1,3); free_dmatrix(TI,1,3,1,3);
  
  free_dvector(data12,1,2*N); free_dvector(data21,1,2*N); free_dvector(data31,1,2*N);
  free_dvector(data13,1,2*N); free_dvector(data23,1,2*N); free_dvector(data32,1,2*N);
  
  free_dvector(a12,1,2*N+2); free_dvector(a21,1,2*N+2); free_dvector(a31,1,2*N+2);
  free_dvector(a13,1,2*N+2); free_dvector(a23,1,2*N+2); free_dvector(a32,1,2*N+2);
  
  free_d3tensor(d,1,3,1,3,1,2*M);
  
  
  return;
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

void spacecraft(struct lisa_orbit *orbit, double tint, double *xint, double *yint, double *zint)
{
  int i;
  for(i=0; i<3; i++)
  {
    LISA_splint(orbit->t, orbit->x[i], orbit->dx[i], orbit->N, tint, &xint[i+1]);
    LISA_splint(orbit->t, orbit->y[i], orbit->dy[i], orbit->N, tint, &yint[i+1]);
    LISA_splint(orbit->t, orbit->z[i], orbit->dz[i], orbit->N, tint, &zint[i+1]);
  }
}

void LISA_spline(double *x, double *y, int n, double yp1, double ypn, double *y2)
// Unlike NR version, assumes zero-offset arrays.  CHECK THAT THIS IS CORRECT.
{
  int i, k;
  double p, qn, sig, un, *u;
  u = dvector(0, n-2);
  // Boundary conditions: Check which is best.
  if (yp1 > 0.99e30)
    y2[0] = u[0] = 0.0;
  else {
    y2[0] = -0.5;
    u[0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }
  
  
  for(i = 1; i < n-1; i++) {
    sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p = sig*y2[i-1] + 2.0;
    y2[i] = (sig-1.0)/p;
    u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  // More boundary conditions.
  if (ypn > 0.99e30)
    qn = un = 0.0;
  else {
    qn = 0.5;
    un = (3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }
  y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  for (k = n-2; k >= 0; k--)
    y2[k] = y2[k]*y2[k+1]+u[k];
  free_dvector(u, 0, n-2);
}

void LISA_splint(double *xa, double *ya, double *y2a, int n, double x, double *y)
{
  // Unlike NR version, assumes zero-offset arrays.  CHECK THAT THIS IS CORRECT.
  int klo,khi,k;
  double h,b,a;
  
  klo=0;
  khi=n-1;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (xa[k] > x) khi=k;
    else klo=k;
  }
  h=xa[khi]-xa[klo];
  
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

void initialize_orbit(char OrbitFile[], struct lisa_orbit *orbit)
{
  int n,i;
  double junk;
  
  FILE *infile = fopen(OrbitFile,"r");
  
  //how big is the file
  n=0;
  while(!feof(infile))
  {
    /*columns of orbit file:
     t sc1x sc1y sc1z sc2x sc2y sc2z sc3x sc3y sc3z
     */
    n++;
    fscanf(infile,"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",&junk,&junk,&junk,&junk,&junk,&junk,&junk,&junk,&junk,&junk);
  }
  n--;
  rewind(infile);
  
  //allocate memory for local workspace
  int N = n;
  double *t  = dvector(0,N-1);
  double **x  = dmatrix(0,2,0,N-1);
  double **y  = dmatrix(0,2,0,N-1);
  double **z  = dmatrix(0,2,0,N-1);
  double **dx = dmatrix(0,2,0,N-1);
  double **dy = dmatrix(0,2,0,N-1);
  double **dz = dmatrix(0,2,0,N-1);
  
  //allocate memory for orbit structure
  orbit->N = n;
  
  orbit->t  = dvector(0,orbit->N-1);
  orbit->x  = dmatrix(0,2,0,orbit->N-1);
  orbit->y  = dmatrix(0,2,0,orbit->N-1);
  orbit->z  = dmatrix(0,2,0,orbit->N-1);
  orbit->dx = dmatrix(0,2,0,orbit->N-1);
  orbit->dy = dmatrix(0,2,0,orbit->N-1);
  orbit->dz = dmatrix(0,2,0,orbit->N-1);
  
  //read in orbits
  for(n=0; n<orbit->N; n++)
  {
    //First time sample must be at t=0 for phasing
    fscanf(infile,"%lg",&t[n]);
    for(i=0; i<3; i++) fscanf(infile,"%lg %lg %lg",&x[i][n],&y[i][n],&z[i][n]);
    
    orbit->t[n] = t[n];
    
    //Repackage orbit positions into arrays for interpolation
    for(i=0; i<3; i++)
    {
      orbit->x[i][n] = x[i][n];
      orbit->y[i][n] = y[i][n];
      orbit->z[i][n] = z[i][n];
    }
  }
  fclose(infile);
  
  //calculate derivatives for cubic spline
  for(i=0; i<3; i++)
  {
    LISA_spline(t, orbit->x[i], N, 1.e30, 1.e30, orbit->dx[i]);
    LISA_spline(t, orbit->y[i], N, 1.e30, 1.e30, orbit->dy[i]);
    LISA_spline(t, orbit->z[i], N, 1.e30, 1.e30, orbit->dz[i]);
  }
  
  //calculate average arm length
  printf("estimating average armlengths -- assumes evenly sampled orbits\n");
  double L12=0.0;
  double L23=0.0;
  double L31=0.0;
  double x1,x2,x3,y1,y2,y3,z1,z2,z3;
  for(n=0; n<orbit->N; n++)
  {
    x1 = orbit->x[0][n];
    x2 = orbit->x[1][n];
    x3 = orbit->x[2][n];
    
    y1 = orbit->y[0][n];
    y2 = orbit->y[1][n];
    y3 = orbit->y[2][n];
    
    z1 = orbit->z[0][n];
    z2 = orbit->z[1][n];
    z3 = orbit->z[2][n];
    
    
    //L12
    L12 += sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) );
    
    //L23
    L23 += sqrt( (x3-x2)*(x3-x2) + (y3-y2)*(y3-y2) + (z3-z2)*(z3-z2) );
    
    //L31
    L31 += sqrt( (x1-x3)*(x1-x3) + (y1-y3)*(y1-y3) + (z1-z3)*(z1-z3) );
  }
  L12 /= (double)orbit->N;
  L31 /= (double)orbit->N;
  L23 /= (double)orbit->N;
  
  printf("Average arm lengths for the constellation:\n");
  printf("  L12 = %g\n",L12);
  printf("  L31 = %g\n",L31);
  printf("  L23 = %g\n",L23);
  
  //are the armlenghts consistent?
  double L = (L12+L31+L23)/3.;
  printf("Fractional deviation from average armlength for each side:\n");
  printf("  L12 = %g\n",fabs(L12-L)/L);
  printf("  L31 = %g\n",fabs(L31-L)/L);
  printf("  L23 = %g\n",fabs(L23-L)/L);
  
  //store armlenght & transfer frequency in orbit structure.
  orbit->L     = L;
  orbit->fstar = clight/(2.0*pi*L);
  
  //free local memory
  free_dvector(t,0,N-1);
  free_dmatrix(x ,0,2,0,N-1);
  free_dmatrix(y ,0,2,0,N-1);
  free_dmatrix(z ,0,2,0,N-1);
  free_dmatrix(dx,0,2,0,N-1);
  free_dmatrix(dy,0,2,0,N-1);
  free_dmatrix(dz,0,2,0,N-1);
}
/*************************************************************************/







void XYZ(double L, double fstar, double ***d, double f0, long q, long M, double *XLS, double *ALS, double *ELS)
{
  int i;
  double fonfs, sqT;
  double ReX, ImX, ReY, ImY, ReZ, ImZ;
  double c3, s3, c2, s2, c1, s1;
  double f,power;
  double *X, *Y, *Z;
  double *YLS, *ZLS;
  double phiLS, phiSL, cLS, sLS, cSL, sSL;
  
  X = dvector(1,2*M);  Y = dvector(1,2*M);  Z = dvector(1,2*M);
  YLS = dvector(1,2*M);  ZLS = dvector(1,2*M);
  
  phiLS = 2.0*pi*f0*(dt/2.0-L/clight);
  /* phiSL = pi/2.0-2.0*pi*f0*(L/clight); */
  cLS = cos(phiLS);
  sLS = sin(phiLS);
  /*cSL = cos(phiSL);
   sSL = sin(phiSL); */
  
  sqT = sqrt(T);
  for(i=1; i<=M; i++)
  {
    
    f = ((double)(q + i-1 - M/2))/T;
    fonfs = f/fstar;
    c3 = cos(3.*fonfs);  c2 = cos(2.*fonfs);  c1 = cos(1.*fonfs);
    s3 = sin(3.*fonfs);  s2 = sin(2.*fonfs);  s1 = sin(1.*fonfs);
    
    
    X[2*i-1] = (d[1][2][2*i-1]-d[1][3][2*i-1])*c3 + (d[1][2][2*i]-d[1][3][2*i])*s3 +
    (d[2][1][2*i-1]-d[3][1][2*i-1])*c2 + (d[2][1][2*i]-d[3][1][2*i])*s2 +
    (d[1][3][2*i-1]-d[1][2][2*i-1])*c1 + (d[1][3][2*i]-d[1][2][2*i])*s1 +
    (d[3][1][2*i-1]-d[2][1][2*i-1]);
    X[2*i]   = (d[1][2][2*i]-d[1][3][2*i])*c3 - (d[1][2][2*i-1]-d[1][3][2*i-1])*s3 +
    (d[2][1][2*i]-d[3][1][2*i])*c2 - (d[2][1][2*i-1]-d[3][1][2*i-1])*s2 +
    (d[1][3][2*i]-d[1][2][2*i])*c1 - (d[1][3][2*i-1]-d[1][2][2*i-1])*s1 +
    (d[3][1][2*i]-d[2][1][2*i]);
    
    Y[2*i-1] = (d[2][3][2*i-1]-d[2][1][2*i-1])*c3 + (d[2][3][2*i]-d[2][1][2*i])*s3 +
    (d[3][2][2*i-1]-d[1][2][2*i-1])*c2 + (d[3][2][2*i]-d[1][2][2*i])*s2+
    (d[2][1][2*i-1]-d[2][3][2*i-1])*c1 + (d[2][1][2*i]-d[2][3][2*i])*s1+
    (d[1][2][2*i-1]-d[3][2][2*i-1]);
    Y[2*i]   = (d[2][3][2*i]-d[2][1][2*i])*c3 - (d[2][3][2*i-1]-d[2][1][2*i-1])*s3+
    (d[3][2][2*i]-d[1][2][2*i])*c2 - (d[3][2][2*i-1]-d[1][2][2*i-1])*s2+
    (d[2][1][2*i]-d[2][3][2*i])*c1 - (d[2][1][2*i-1]-d[2][3][2*i-1])*s1+
    (d[1][2][2*i]-d[3][2][2*i]);
    
    Z[2*i-1] = (d[3][1][2*i-1]-d[3][2][2*i-1])*c3 + (d[3][1][2*i]-d[3][2][2*i])*s3+
    (d[1][3][2*i-1]-d[2][3][2*i-1])*c2 + (d[1][3][2*i]-d[2][3][2*i])*s2+
    (d[3][2][2*i-1]-d[3][1][2*i-1])*c1 + (d[3][2][2*i]-d[3][1][2*i])*s1+
    (d[2][3][2*i-1]-d[1][3][2*i-1]);
    
    Z[2*i]   = (d[3][1][2*i]-d[3][2][2*i])*c3 - (d[3][1][2*i-1]-d[3][2][2*i-1])*s3+
    (d[1][3][2*i]-d[2][3][2*i])*c2 - (d[1][3][2*i-1]-d[2][3][2*i-1])*s2+
    (d[3][2][2*i]-d[3][1][2*i])*c1 - (d[3][2][2*i-1]-d[3][1][2*i-1])*s1+
    (d[2][3][2*i]-d[1][3][2*i]);
    
    
    /* XSL[2*i-1] = 2.0*fonfs*(X[2*i-1]*cSL-X[2*i]*sSL);
     XSL[2*i] = -2.0*fonfs*(X[2*i-1]*sSL+X[2*i]*cSL);
     YSL[2*i-1] = 2.0*fonfs*(Y[2*i-1]*cSL-Y[2*i]*sSL);
     YSL[2*i] = -2.0*fonfs*(Y[2*i-1]*sSL+Y[2*i]*cSL);
     ZSL[2*i-1] = 2.0*fonfs*(Z[2*i-1]*cSL-Z[2*i]*sSL);
     ZSL[2*i] = -2.0*fonfs*(Z[2*i-1]*sSL+Z[2*i]*cSL); */
    
    XLS[2*i-1] = (X[2*i-1]*cLS-X[2*i]*sLS);
    XLS[2*i] = -(X[2*i-1]*sLS+X[2*i]*cLS);
    YLS[2*i-1] = (Y[2*i-1]*cLS-Y[2*i]*sLS);
    YLS[2*i] = -(Y[2*i-1]*sLS+Y[2*i]*cLS);
    ZLS[2*i-1] = (Z[2*i-1]*cLS-Z[2*i]*sLS);
    ZLS[2*i] = -(Z[2*i-1]*sLS+Z[2*i]*cLS);
    
    
  }
  
  
  for(i=1; i<=2*M; i++)
  {
    // A channel
    ALS[i] = (2.0*XLS[i] - YLS[i] - ZLS[i])/3.0;
    // E channel
    ELS[i] = (ZLS[i]-YLS[i])/sq3;
  }
  
  free_dvector(X,1,2*M);  free_dvector(Y,1,2*M);  free_dvector(Z,1,2*M);
  free_dvector(YLS,1,2*M);  free_dvector(ZLS,1,2*M);
  
}

void instrument_noise(double f, double fstar, double L, double *SAE, double *SXYZ)
{
  //Power spectral density of the detector noise and transfer frequency
  double Sn, red, confusion_noise;
  double f1, f2;
  double A1, A2, slope;
  FILE *outfile;
  
  red = (2.0*pi*1.0e-4)*(2.0*pi*1.0e-4);
  
  
  // Calculate the power spectral density of the detector noise at the given frequency
  
  //NGO reddening
  //*SAE = 16.0/3.0*pow(sin(f/fstar),2.0)*( ( (2.0+cos(f/fstar))*Sps + 2.0*(3.0+2.0*cos(f/fstar)+cos(2.0*f/fstar))*Sacc*(1.+(0.0001/f))*(1.0/pow(2.0*pi*f,4)+ red/pow(2.0*pi*f,6))) / pow(2.0*L,2.0));
  //*SXYZ = 4.0*pow(sin(f/fstar),2.0)*( ( 4.0*Sps + 8.0*(1.0+pow(cos(f/fstar),2.0))*Sacc*(1.+(0.0001/f))*(1.0/pow(2.0*pi*f,4)+ red/pow(2.0*pi*f,6))) / pow(2.0*L,2.0));
  
  *SAE = 16.0/3.0*pow(sin(f/fstar),2.0)*( ( (2.0+cos(f/fstar))*Sps + 2.0*(3.0+2.0*cos(f/fstar)+cos(2.0*f/fstar))*Sacc*(1.0/pow(2.0*pi*f,4)+ red/pow(2.0*pi*f,6))) / pow(2.0*L,2.0));
  *SXYZ = 4.0*pow(sin(f/fstar),2.0)*( ( 4.0*Sps + 8.0*(1.0+pow(cos(f/fstar),2.0))*Sacc*(1.0/pow(2.0*pi*f,4)+ red/pow(2.0*pi*f,6))) / pow(2.0*L,2.0));
  
}

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void dfour1(double data[], unsigned long nn, int isign)
{
  unsigned long n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  double tempr,tempi;
  
  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) {
    if (j > i) {
      SWAP(data[j],data[i]);
      SWAP(data[j+1],data[i+1]);
    }
    m=n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (n > mmax) {
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
        j=i+mmax;
        tempr=wr*data[j]-wi*data[j+1];
        tempi=wr*data[j+1]+wi*data[j];
        data[j]=data[i]-tempr;
        data[j+1]=data[i+1]-tempi;
        data[i] += tempr;
        data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}
#undef SWAP

void KILL(char* Message)
{
  printf("\a\n");
  printf("%s",Message);
  printf("Terminating the program.\n\n\n");
  exit(1);
  
  
  return;
}

