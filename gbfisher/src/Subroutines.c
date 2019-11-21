#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "arrays.h"
#include "Constants.h"
#include "Detector.h"
#include "Subroutines.h"

#include <LISA.h>

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress (double percentage)
{
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush (stdout);
}

void XYZ(double L, double fstar, double TOBS, double ***d, double f0, long q, long M, double *XLS, double *ALS, double *ELS)
{
  int i;
  double fonfs, sqT;
  double c3, s3, c2, s2, c1, s1;
  double f;
  double *X, *Y, *Z;
  double *YLS, *ZLS;
  double phiLS, cLS, sLS;
  
  X = dvector(1,2*M);  Y = dvector(1,2*M);  Z = dvector(1,2*M);
  YLS = dvector(1,2*M);  ZLS = dvector(1,2*M);
  
  phiLS = 2.0*pi*f0*(DT/2.0-L/clight);
  /* phiSL = pi/2.0-2.0*pi*f0*(L/clight); */
  cLS = cos(phiLS);
  sLS = sin(phiLS);
  /*cSL = cos(phiSL);
   sSL = sin(phiSL); */
  
  sqT = sqrt(TOBS);
  for(i=1; i<=M; i++)
  {
    
    f = ((double)(q + i-1 - M/2))/TOBS;
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

double M_fdot(double f, double fdot)
{
  return pow( fdot*pow(f,-11./3.)*(5./96.)*pow(M_PI,-8./3.)  ,  3./5.)/TSUN;
}

void instrument_noise(double f, double fstar, double L, double *SAE, double *SXYZ)
{
  
  double SLOC = 2.89e-24; //what is this?

  double fonfstar = f/fstar;
  double trans = pow(sin(fonfstar),2.0);
  double red = 1.0 + 16.0*(pow((2.0e-5/f), 10.0) + pow(1.0e-4/f,2));
  
  double cosfstar = cos(fonfstar);
  
  double L4 = 4.0*L*L;
  double f4 = 1./(pow(PI2*f,4));
  
  *SAE  = (16.0/3.0)*trans*( (2.0+cosfstar)*(SPS + SLOC) + 2.0*( 3.0 + 2.0*cosfstar + cos(2.0*fonfstar) ) * ( SLOC/2.0 + SACC*f4*red ) ) / L4;
  *SXYZ =      (4.0)*trans*(          (4.0)*(SPS + SLOC) + 8.0*( 1.0 +                pow(cosfstar,2.0) ) * ( SLOC/2.0 + SACC*f4*red ) ) / L4;
}

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

//void LISA_spline(double *x, double *y, int n, double yp1, double ypn, double *y2)
//// Unlike NR version, assumes zero-offset arrays.  CHECK THAT THIS IS CORRECT.
//{
//  int i, k;
//  double p, qn, sig, un, *u;
//  u = dvector(0, n-2);
//  // Boundary conditions: Check which is best.
//  if (yp1 > 0.99e30)
//    y2[0] = u[0] = 0.0;
//  else {
//    y2[0] = -0.5;
//    u[0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
//  }
//  
//  
//  for(i = 1; i < n-1; i++) {
//    sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
//    p = sig*y2[i-1] + 2.0;
//    y2[i] = (sig-1.0)/p;
//    u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
//    u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
//  }
//  // More boundary conditions.
//  if (ypn > 0.99e30)
//    qn = un = 0.0;
//  else {
//    qn = 0.5;
//    un = (3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
//  }
//  y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.0);
//  for (k = n-2; k >= 0; k--)
//    y2[k] = y2[k]*y2[k+1]+u[k];
//  free_dvector(u, 0, n-2);
//}
//
//void LISA_splint(double *xa, double *ya, double *y2a, int n, double x, double *y)
//{
//  // Unlike NR version, assumes zero-offset arrays.  CHECK THAT THIS IS CORRECT.
//  int klo,khi,k;
//  double h,b,a;
//  
//  klo=0;
//  khi=n-1;
//  while (khi-klo > 1) {
//    k=(khi+klo) >> 1;
//    if (xa[k] > x) khi=k;
//    else klo=k;
//  }
//  h=xa[khi]-xa[klo];
//  
//  a=(xa[khi]-x)/h;
//  b=(x-xa[klo])/h;
//  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
//}

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
  //printf("estimating average armlengths -- assumes evenly sampled orbits\n");
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
  
  /*
  printf("Average arm lengths for the constellation:\n");
  printf("  L12 = %g\n",L12);
  printf("  L31 = %g\n",L31);
  printf("  L23 = %g\n",L23);
  */
  //are the armlenghts consistent?
  double L = (L12+L31)/2.;//+L23)/3.;
  /*
  printf("Fractional deviation from average armlength for each side:\n");
  printf("  L12 = %g\n",fabs(L12-L)/L);
  printf("  L31 = %g\n",fabs(L31-L)/L);
  printf("  L23 = %g\n",fabs(L23-L)/L);
   */

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


double Sum(double *AA, double *EE, long M, double SN, double TOBS)
{
  long i;
  double sm;
  
  sm = 0.0;
  
  for(i=1; i<=M; i++)
  {
    sm += (AA[2*i-1]*EE[2*i-1]+AA[2*i]*EE[2*i]);
  }
  
  sm *= 4.0*TOBS/SN;
  
  return sm;
  
}


/*****************************************************/
/*                                                   */
/*        Median-based Confusion Noise Fitting       */
/*                                                   */
/*****************************************************/



void medianX(long imin, long imax, double fstar, double L, double *XP, double *Xnoise, double *Xconf, double TOBS)
{
  printf(" Median fit to X-channel confusion noise\n");

  double f;
  double SAE, SXYZ;
  double chi;
  long i, j, k;
  long segs;
  long rseed;
  int Npoly;
  double *XX;
  double *fdata, *mdata, *pcx, *pcy, *inst;
  double chix, chiy, fit, alpha, beta, mul, conf;
  double lfmin, lfmax, dlf, lf, ln4;
  FILE *Xfile;
  
  XX = dvector(0,100);
  
  rseed = -546214;
  
  segs = (int)((double)(imax-imin)/101.0);
  
  lfmin = log((double)(imin-101)/TOBS);
  lfmax = log((double)(imin+101*(segs))/TOBS);
  
  
  Npoly = 30;
  
  dlf = (lfmax-lfmin)/(double)(Npoly);
  ln4 = log(1.0e-4);
  
  fdata = dvector(0,segs-1);
  mdata = dvector(0,segs-1);
  inst = dvector(0,segs-1);
  pcx = dvector(0,Npoly);
  pcy = dvector(0,Npoly);
  
  for(i=0; i < segs; i++)
  {
    for(j=0; j<=100; j++) XX[j] = XP[imin+101*i+j];
    f = (double)(imin+101*i-50)/TOBS;
    instrument_noise(f, fstar, L, &SAE, &SXYZ);
    inst[i] = log(SXYZ*1.0e40);
    chi=quickselect(XX, 101, 51);
    //printf("%e %e\n", f, chi/0.72);
    fdata[i] = log(f);
    mdata[i] = log(chi/0.72*1.0e40);
  }
  
  // initialize fit
  for(i=1; i < Npoly; i++)
  {
    f = exp(lfmin+(double)(i)*dlf);
    j = (long)((f*TOBS-(double)(imin-50))/101.0);
    //printf("%ld %ld\n", i, j);
    pcx[i] = mdata[j];
    //printf("%e %e\n", f, exp(pcx[i])*1.0e-40);
  }
  pcx[0] = pcx[1];
  pcx[Npoly] = pcx[Npoly-1];
  //printf("%ld\n", segs);
  
  
  chix = 0.0;
  for(i=0; i < segs; i++)
  {
    lf = log((double)(imin+101*i-50)/TOBS);
    j = (long)floor((lf-lfmin)/dlf);
    fit = pcx[j]+((pcx[j+1]-pcx[j])/dlf)*(lf-(lfmin+(double)(j)*dlf));
    f = exp(-1.0*((lfmin+((double)(j)+0.5)*dlf)-ln4));
    chix += 500.0*(mdata[i] - fit)*(mdata[i] - fit)*f;
    //printf("%ld %e\n", i, (lf-(lfmin+(double)(j)*dlf))/dlf);
    //printf("%e %e %e\n", exp(lf), exp(mdata[i])*1.0e-40, exp(fit)*1.0e-40);
  }
  
  for(k=0; k < 10000; k++)
  {
    if(k%(10000/100)==0)printProgress((double)k/10000.);
    beta = pow(10.0,-0.0001*(10000.0-(double)(k))/10000.0);
    
    for(j=0; j <= Npoly; j++) pcy[j] = pcx[j];
    j = (long)((double)(Npoly+1)*ran2(&rseed));
    mul = 1.0;
    alpha = ran2(&rseed);
    if(alpha > 0.3) mul = 10.0;
    if(alpha > 0.7) mul = 0.01;
    pcy[j] = pcx[j]+mul*0.01*gasdev2(&rseed)/sqrt(beta);
    
    chiy = 0.0;
    for(i=0; i < segs; i++)
    {
      lf = log((double)(imin+101*i-50)/TOBS);
      j = (long)floor((lf-lfmin)/dlf);
      fit = pcy[j]+((pcy[j+1]-pcy[j])/dlf)*(lf-(lfmin+(double)(j)*dlf));
      f = exp(-1.0*((lfmin+((double)(j)+0.5)*dlf)-ln4));
      chiy += 500.0*(mdata[i] - fit)*(mdata[i] - fit)*f;
    }
    
    alpha = log(ran2(&rseed));
    
    if(beta*(chix-chiy) > alpha)
    {
      chix = chiy;
      for(j=0; j <= Npoly; j++) pcx[j] = pcy[j];
    }
    //if(k%100==0) printf("%ld %.10e %.10e\n", k, chix, chiy);
  }
  printProgress(1.0);
  
  printf("\n Store X-channel results\n");
  
  Xfile = fopen("Xfit.dat","w");
  for(i=0; i < segs; i++)
  {
    lf = log((double)(imin+101*i-50)/TOBS);
    j = (long)floor((lf-lfmin)/dlf);
    fit = pcx[j]+((pcx[j+1]-pcx[j])/dlf)*(lf-(lfmin+(double)(j)*dlf));
    fprintf(Xfile, "%e %e %e\n", exp(lf), exp(fit)*1.0e-40, exp(mdata[i])*1.0e-40);
  }
  fclose(Xfile);
  
  Xfile = fopen("Xf.dat","w");
  for(i=0; i <= Npoly; i++)
  {
    lf = lfmin+(double)(i)*dlf;
    fprintf(Xfile, "%e %e\n", exp(lf), exp(pcx[i])*1.0e-40);
  }
  fclose(Xfile);
  
  
  
  for(i=imin; i <= imax; i++)
  {
    f = (double)(i)/TOBS;
    lf = log(f);
    j = (long)floor((lf-lfmin)/dlf);
    fit = pcx[j]+((pcx[j+1]-pcx[j])/dlf)*(lf-(lfmin+(double)(j)*dlf));
    instrument_noise(f, fstar, L, &SAE, &SXYZ);
    alpha = exp(fit)*1.0e-40;
    conf = alpha -SXYZ;
    if(conf < SXYZ/30.0) conf = 1.0e-46;
    Xnoise[i] = alpha;
    Xconf[i] = conf;
    
    
    
    //Xnoise[i]*=pow(sin(f/fstar),2.0);
    //Xconf[i]*=pow(sin(f/fstar),2.0);
  }
  
  free_dvector(XX, 0,100);
  free_dvector(fdata, 0,segs-1);
  free_dvector(mdata, 0,segs-1);
  free_dvector(pcx, 0,Npoly);
  free_dvector(pcy, 0,Npoly);
  
  return;
}

void medianAE(long imin, long imax, double fstar, double L, double *AEP, double *AEnoise, double *AEconf, double TOBS)
{
  printf(" Median fit to AE-channel confusion noise\n");
  double f;
  double SAE, SXYZ;
  double chi;
  long i, j, k;
  long segs;
  long rseed;
  int Npoly;
  double *XX;
  double *fdata, *mdata, *pcx, *pcy, *inst;
  double chix, chiy, fit, alpha, beta, mul, conf;
  double lfmin, lfmax, dlf, lf, ln4;
  FILE *Xfile;
  
  XX = dvector(0,100);
  
  rseed = -546214;
  
  segs = (int)((double)(imax-imin)/101.0);
  
  lfmin = log((double)(imin-101)/TOBS);
  lfmax = log((double)(imin+101*(segs))/TOBS);
  
  
  Npoly = 30;
  
  dlf = (lfmax-lfmin)/(double)(Npoly);
  ln4 = log(1.0e-4);
  
  fdata = dvector(0,segs-1);
  mdata = dvector(0,segs-1);
  inst = dvector(0,segs-1);
  pcx = dvector(0,Npoly);
  pcy = dvector(0,Npoly);
  
  for(i=0; i < segs; i++)
  {
    for(j=0; j<=100; j++) XX[j] = AEP[imin+101*i+j];
    f = (double)(imin+101*i-50)/TOBS;
    instrument_noise(f, fstar, L, &SAE, &SXYZ);
    inst[i] = log(SAE*1.0e40);
    chi=quickselect(XX, 101, 51);
    //printf("%e %e\n", f, chi/0.72);
    fdata[i] = log(f);
    mdata[i] = log(chi/0.72*1.0e40);
  }
  
  // initialize fit
  for(i=1; i < Npoly; i++)
  {
    f = exp(lfmin+(double)(i)*dlf);
    j = (long)((f*TOBS-(double)(imin-50))/101.0);
    //printf("%ld %ld\n", i, j);
    pcx[i] = mdata[j];
    //printf("%e %e\n", f, exp(pcx[i])*1.0e-40);
  }
  pcx[0] = pcx[1];
  pcx[Npoly] = pcx[Npoly-1];
  //printf("%ld\n", segs);
  
  
  chix = 0.0;
  for(i=0; i < segs; i++)
  {
    lf = log((double)(imin+101*i-50)/TOBS);
    j = (long)floor((lf-lfmin)/dlf);
    fit = pcx[j]+((pcx[j+1]-pcx[j])/dlf)*(lf-(lfmin+(double)(j)*dlf));
    f = exp(-1.0*((lfmin+((double)(j)+0.5)*dlf)-ln4));
    chix += 500.0*(mdata[i] - fit)*(mdata[i] - fit)*f;
    //printf("%ld %e\n", i, (lf-(lfmin+(double)(j)*dlf))/dlf);
    //printf("%e %e %e\n", exp(lf), exp(mdata[i])*1.0e-40, exp(fit)*1.0e-40);
  }
  
  for(k=0; k < 10000; k++)
  {
    if(k%(10000/100)==0)printProgress((double)k/10000.);
    beta = pow(10.0,-0.0001*(10000.0-(double)(k))/10000.0);
    
    for(j=0; j <= Npoly; j++) pcy[j] = pcx[j];
    j = (long)((double)(Npoly+1)*ran2(&rseed));
    mul = 1.0;
    alpha = ran2(&rseed);
    if(alpha > 0.3) mul = 10.0;
    if(alpha > 0.7) mul = 0.01;
    pcy[j] = pcx[j]+mul*0.01*gasdev2(&rseed)/sqrt(beta);
    
    chiy = 0.0;
    for(i=0; i < segs; i++)
    {
      lf = log((double)(imin+101*i-50)/TOBS);
      j = (long)floor((lf-lfmin)/dlf);
      fit = pcy[j]+((pcy[j+1]-pcy[j])/dlf)*(lf-(lfmin+(double)(j)*dlf));
      f = exp(-1.0*((lfmin+((double)(j)+0.5)*dlf)-ln4));
      chiy += 500.0*(mdata[i] - fit)*(mdata[i] - fit)*f;
    }
    
    alpha = log(ran2(&rseed));
    
    if(beta*(chix-chiy) > alpha)
    {
      chix = chiy;
      for(j=0; j <= Npoly; j++) pcx[j] = pcy[j];
    }
    //if(k%100==0) printf("%ld %.10e %.10e\n", k, chix, chiy);
  }
  printProgress(1.0);
  
  printf("\n Store AE-channel results\n");
  
  Xfile = fopen("Afit.dat","w");
  for(i=0; i < segs; i++)
  {
    lf = log((double)(imin+101*i-50)/TOBS);
    j = (long)floor((lf-lfmin)/dlf);
    fit = pcx[j]+((pcx[j+1]-pcx[j])/dlf)*(lf-(lfmin+(double)(j)*dlf));
    fprintf(Xfile, "%e %e %e\n", exp(lf), exp(fit)*1.0e-40, exp(mdata[i])*1.0e-40);
  }
  fclose(Xfile);
  
  for(i=imin; i <= imax; i++)
  {
    f = (double)(i)/TOBS;
    lf = log(f);
    j = (long)floor((lf-lfmin)/dlf);
    fit = pcx[j]+((pcx[j+1]-pcx[j])/dlf)*(lf-(lfmin+(double)(j)*dlf));
    instrument_noise(f, fstar, L, &SAE, &SXYZ);
    alpha = exp(fit)*1.0e-40;
    conf = alpha -SAE;
    if(conf < SAE/30.0) conf = 1.0e-46;
    AEnoise[i] = alpha;
    AEconf[i] = conf;
    
    
    //AEnoise[i]*=pow(sin(f/fstar),2.0);
    //AEconf[i]*=pow(sin(f/fstar),2.0);
    
  }
  
  free_dvector(XX, 0,100);
  free_dvector(fdata, 0,segs-1);
  free_dvector(mdata, 0,segs-1);
  free_dvector(pcx, 0,Npoly);
  free_dvector(pcy, 0,Npoly);
  
  return;
}

void KILL(char* Message)
{
  printf("\a\n");
  printf("%s",Message);
  printf("Terminating the program.\n\n\n");
  exit(1);
}

double gasdev2(long *idum)
{
  double ran2(long *idum);
  static int iset=0;
  static double gset;
  double fac, rsq, v1, v2;
  
  if(*idum < 0) iset = 0;
  if(iset == 0){
    do{
      v1 = 2.0 * ran2(idum)-1.0;
      v2 = 2.0 * ran2(idum)-1.0;
      rsq = v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset = v1 * fac;
    iset = 1;
    return(v2*fac);
  } else {
    iset = 0;
    return (gset);
  }
}

double ran2(long *idum)
{
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  double temp;
  
  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) *idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

// for n odd, median given by k = (n+1)/2. for n even, median is average of k=n/2, k=n/2+1 elements

double quickselect(double *arr, int n, int k)
{
  unsigned long i,ir,j,l,mid;
  double a,temp;
  
  l=0;
  ir=n-1;
  for(;;) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
        SWAP(arr[l],arr[ir]);
      }
      return arr[k];
    }
    else {
      mid=(l+ir) >> 1;
      SWAP(arr[mid],arr[l+1]);
      if (arr[l] > arr[ir]) {
        SWAP(arr[l],arr[ir]);
      }
      if (arr[l+1] > arr[ir]) {
        SWAP(arr[l+1],arr[ir]);
      }
      if (arr[l] > arr[l+1]) {
        SWAP(arr[l],arr[l+1]);
      }
      i=l+1;
      j=ir;
      a=arr[l+1];
      for (;;) {
        do i++; while (arr[i] < a);
        do j--; while (arr[j] > a);
        if (j < i) break;
        SWAP(arr[i],arr[j]);
      }
      arr[l+1]=arr[j];
      arr[j]=a;
      if (j >= k) ir=j-1;
      if (j <= k) l=i;
    }
  }
}

/*****************************************************/
/*                                                   */
/*        Spline-based Confusion Noise Fitting       */
/*                                                   */
/*****************************************************/

void spline_fit(int flag, int divs, long imin, long imax, double *XP, double *Xnoise, double *Xconf, double T, double fstar, double L)
{
  double f;
  double SAE, SXYZ;
  double chi;
  long i, j;
  long segs;
  long rseed;
  int Npoly;
  double *XX;
  double *fdata, *mdata, *pcx, *pcy, *inst;
  double *adata, *sdata;
  double conf;
  double lfmin, lfmax, x, dlf, ln4;
  int divsp, divh;
  double mean, var;
  
  divsp = divs+1;
  divh = divs/2;
  
  XX = dvector(0,divs);
  
  rseed = -546214;
  
  segs = (int)((double)(imax-imin)/(double)(divsp));
  
  lfmin = log((double)(imin)/T);
  lfmax = log((double)(imin+divs*(segs))/T);
  
  Npoly = 20;
  
  dlf = (lfmax-lfmin)/(double)(Npoly);
  ln4 = log(1.0e-4);
  
  fdata = dvector(0,segs-1);
  mdata = dvector(0,segs-1);
  adata = dvector(0,segs-1);
  sdata = dvector(0,segs-1);
  inst = dvector(0,segs-1);
  pcx = dvector(0,Npoly);
  pcy = dvector(0,Npoly);
  
  if(flag==0)   printf(" Spline fit to X-channel confusion noise\n");
  if(flag==1)   printf(" Spline fit to AE-channel confusion noise\n");

  for(i=0; i < segs; i++)
  {
    for(j=0; j<=divs; j++)
    {
      XX[j] = XP[imin-divh+divsp*i+j]*1.0e40;
    }

    mean = 0.0;
    var = 0.0;
    for(j=0; j<=divs; j++)
    {
      x = log(XX[j]);
      mean += x;
      var += x*x;
    }
    mean /= (double)(divs+1);
    var /= (double)(divs+1);
    var = sqrt(var-mean*mean);
    var /= sqrt((double)(divs+1));  // deviation of mean
    mean += 0.57721566490153286060;  // have to add Euler's constant since taking average of log
    adata[i] = mean;
    sdata[i] = var;
    f = (double)(imin+divsp*i)/T;
    instrument_noise(f, fstar, L, &SAE, &SXYZ);
    if(flag == 0) inst[i] = log(SXYZ*1.0e40);
    if(flag == 1) inst[i] = log(SAE*1.0e40);
    
    chi=quickselect(XX, divsp, (divh+1));
    fdata[i] = log(f);
    mdata[i] = log(chi/0.72);
  }

  splineMCMC(imin, imax, segs, fdata, mdata, sdata, Xnoise, T);

  for(i=imin; i <= imax; i++)
  {
    f = (double)(i)/T;
    instrument_noise(f, fstar, L, &SAE, &SXYZ);
    if(flag == 0)
    {
      conf = Xnoise[i] -SXYZ;
      if(conf < SXYZ/30.0) conf = 1.0e-46;
    }
    if(flag == 1)
    {
      conf = Xnoise[i] -SAE;
      if(conf < SAE/30.0) conf = 1.0e-46;
    }
    
    Xconf[i] = conf;
  }
  
  
  return;
}

void splineMCMC(int imin, int imax, int ND, double *datax, double *datay, double *sigma, double *Xnoise, double T)
{
  
  
  int N, Nf, Nx, Ny;
  int mc, i, j, ii, test;
  long seed;
  int ltest;
  int fixedD;
  double logLx, logLy;
  double logpx, logpy;
  double *ref;
  int *activex, *activey;
  double alpha, beta, H, av;
  double model;
  double sh;
  double max, min;
  double fmin, fmax;
  double rmin, rmax;
  double maxy, miny;
  double x, y;
  double q;
  double f;
  double *mdl;
  double lmax, lmin;
  double lambdax, lambday;
  int acc=0;
  
  N = 100000;   // number of MCMC steps
  
  Nf = 20; // maximum number of spline control points
  
  seed = -946673524;
  
  // log likelihood fixed test. set this flag = 1 to fix likelihood
  ltest = 0;
  
  
  maxy = -1.0e60;
  miny = 1.0e60;
  
  lmax = 100.0;
  lmin = -100.0;
  
  lambdax = lambday = 1.0;
  
  
  max = datax[ND-1];
  min = datax[0];
  
  // printf("min %e max %e\n", min, max);
  
  double *sdatay, *sderiv;
  double *spoints, *sdatax, *tdata, *tpoints;
  
  mdl = dvector(0,ND);
  ref = dvector(1,Nf);
  spoints = dvector(1,Nf);
  tdata = dvector(1,Nf);
  tpoints = dvector(1,Nf);
  sdatax = dvector(1,Nf);
  sdatay = dvector(1,Nf);
  sderiv = dvector(1,Nf);
  
  activex = ivector(1,Nf);
  activey = ivector(1,Nf);
  
  // only start with 3 active points
  Nx = Ny = Nf;
  
  for(i=1; i<= Nf; i++)
  {
    activex[i] = 1;
    activey[i] = 1;
  }
  
  activex[1] = 1;
  activey[1] = 1;
  activex[Nf] = 1;
  activey[Nf] = 1;
  activex[Nf/4] = 1;
  activey[Nf/4] = 1;
  activex[Nf/2] = 1;
  activey[Nf/2] = 1;
  activex[3*Nf/4] = 1;
  activey[3*Nf/4] = 1;
  
  rmin = 1.0e10;
  rmax = -1.0e10;
  for(j=0; j< ND; j++)
  {
    if(datay[j] > rmax) rmax = datay[j];
    if(datay[j] < rmin) rmin = datay[j];
  }
  
  x = 0.2*(rmax-rmin);
  
  rmin -= x;
  rmax += x;
  
  
  
  
  // inititate fit
  for(i=1; i<= Nf; i++)
  {
    spoints[i]= min + (max-min)/(double)(Nf-1)*(double)(i-1);
    j = -1;
    do
    {
      j++;
      //printf("fucking do while loops: %d %d %f %f\n", i, j, spoints[i], datax[j]);
    }while(spoints[i] > datax[j]);
    
    sdatax[i] = datay[j];
    if(j > 0)
    {
      sdatax[i] = datay[j-1] +(datay[j]-datay[j-1])/(datax[j]-datax[j-1])*(spoints[i]-datax[j-1]);
    }
  }
  
  
  fmin = exp(spoints[1]);
  fmax = exp(spoints[Nf]);
  
  
  i = 0;
  for(ii=1; ii<= Nf; ii++)
  {
    if(activey[ii] == 1)  // only use active points
    {
      i++;
      tpoints[i] = spoints[ii];
      tdata[i] = sdatax[ii];
    }
  }
  
  spline(tpoints, tdata, Nx, 1.0e31, 1.0e31, sderiv);
  
  av = 0.0;
  for(j=0; j< ND; j++)
  {
    splint(tpoints, tdata, sderiv, Nx, datax[j], &model);
    mdl[j] = model;
    y = (datay[j]-model)/sigma[j];
    av -= y*y;
  }
  logLx = av/2.0;
  
  beta = exp(lambdax);
  y = 0.0;
  for(i=2; i< ND; i++)
  {
    x = ((mdl[i]-mdl[i-1])/(datax[i]-datax[i-1])- (mdl[i-1]-mdl[i-2])/(datax[i-1]-datax[i-2]))/(datax[i]-datax[i-2]);
    x /= beta;
    y += x*x;
  }
  // 1/(sqrt(2Pi) beta) exp(-x*x/(2 beta^2))
  logpx = -y/2.0-(double)(ND)*lambdax;
  
  
  // set this flag to 1 if you want to do a fixed dimension run
  fixedD = 1;
  
  double logLmax=-1e60;
  if(ltest == 1)
  {
    logLx = 0.0;
    logpx = 0.0;
  }
  
  // start the RJMCMC
  for(mc=0; mc< N; mc++)
  {
    
    sh = 1.0/sqrt((double)(Nx));
    
    lambday = lambdax;
    Ny = Nx;
    
    test = 0;
    
    for(i=1; i<= Nf; i++)
    {
      sdatay[i] = sdatax[i];
      activey[i] = activex[i];
    }
    
    alpha = ran2(&seed);
    
    q = 0.5;
    if(fixedD == 1) q = 10.0;
    
    if(alpha > q)   // propose a dimension change
    {
      // Note that the spline points at the ends are never added or subtracted
      
      
      alpha = ran2(&seed);
      
      if(alpha < 0.5)
      {
        Ny = Nx + 1;
      }
      else
      {
        Ny = Nx - 1;
      }
      
      if(Ny < Nx) // propose a kill
      {
        if(Ny > 1 && Ny <= Nf)
        {
          
          do
          {
            i = 2 + (int)(ran2(&seed)*(double)(Nf-2)); // pick one to kill
          } while(activex[i] == 0);  // can't kill it if already dead
          activey[i] = 0;
        }
        else
        {
          test = 1;
        }
      }
      else
      {
        if(Ny >= 1 && Ny < Nf)
        {
          
          do
          {
            i = 2 + (int)(ran2(&seed)*(double)(Nf-2)); // pick one to add
          } while(activex[i] == 1);  // can't add if already in use
          activey[i] = 1;
          
          sdatay[i] = rmin + (rmax-rmin)*ran2(&seed);  // draw from prior
          
        }
        else
        {
          test = 1;
        }
        
        
      }
      
      
    }
    else     // within dimension update
    {
      
      Ny = Nx;
      
      alpha = ran2(&seed);
      
      if(alpha > 0.6)  // update all points
      {
        
        for(ii=1; ii<= Nf; ii++)
        {
          // variety of jump sizes
          if(alpha > 0.8)
          {
            sdatay[ii] += sh*1.0e-1*gasdev2(&seed);
          }
          else if (alpha > 0.5)
          {
            sdatay[ii] += sh*1.0e-2*gasdev2(&seed);
          }
          else if (alpha > 0.3)
          {
            sdatay[ii] += sh*1.0e-3*gasdev2(&seed);
          }
          else
          {
            sdatay[ii] += sh*1.0e-4*gasdev2(&seed);
          }
          
        }
        
      }
      else  if(alpha > 0.1) // just update one
      {
        
        do
        {
          ii = (int)(ran2(&seed)*(double)(Nf));
        }while(activey[ii] == 0);
        
        sdatay[ii] += sh*1.0e-1*gasdev2(&seed);
        
      }
      else  // birth/death
      {
        do
        {
          i = 2 + (int)(ran2(&seed)*(double)(Nf-2)); // pick one to kill
        } while(activex[i] == 0);  // can't kill it if already dead
        activey[i] = 0;
        
        do
        {
          i = 2 + (int)(ran2(&seed)*(double)(Nf-2)); // pick one to add
        } while(activey[i] == 1);  // can't add if already in use
        activey[i] = 1;
        
        sdatay[i] = rmin + (rmax-rmin)*ran2(&seed);  // draw from prior
        
        
      }
      
      
    }
    
    
    alpha = ran2(&seed);
    if(alpha > 0.9)
    {
      lambday = lmin+(lmax-lmin)*ran2(&seed);  // uniform draw from the prior
    }
    else if (alpha > 0.7)
    {
      lambday = lambdax + 0.1*gasdev2(&seed);
    }
    else if (alpha > 0.3)
    {
      lambday = lambdax + 0.01*gasdev2(&seed);
    }
    else
    {
      lambday = lambdax + 0.001*gasdev2(&seed);
    }
    
    
    // check that proposed values are within the prior range
    if(Ny < 2 || Ny > Nf) test = 1;
    if(lambday < lmin || lambday > lmax) test = 1;
    
    
    
    logLy = 0.0;
    
    if(test == 0)
    {
      
      if(ltest == 0)
      {
        i = 0;
        for(ii=1; ii<= Nf; ii++)
        {
          if(activey[ii] == 1)  // only use active points
          {
            i++;
            tpoints[i] = spoints[ii];
            tdata[i] = sdatay[ii];
          }
        }
        
        
        spline(tpoints, tdata, Ny, 1.0e31, 1.0e31, sderiv);
        av = 0.0;
        for(j=0; j< ND; j++)
        {
          splint(tpoints, tdata, sderiv, Ny, datax[j], &model);
          mdl[j] = model;
          y = (datay[j]-model)/sigma[j];
          av -= y*y;
        }
        logLy = av/2.0;
        
        beta = exp(lambday);
        y = 0.0;
        for(i=2; i< ND; i++)
        {
          x = ((mdl[i]-mdl[i-1])/(datax[i]-datax[i-1])- (mdl[i-1]-mdl[i-2])/(datax[i-1]-datax[i-2]))/(datax[i]-datax[i-2]);
          x /= beta;
          y += x*x;
        }
        // 1/(sqrt(2Pi) beta) exp(-x*x/(2 beta^2))
        logpy = -y/2.0-(double)(ND)*lambday;
        
      }
      else
      {
        logLy = 0.0;
        logpy = 0.0;
      }
      
    }
    
    
    
    H = (logLy-logLx) +logpy  - logpx;
    
    alpha = log(ran2(&seed));
    
    if((H > alpha) && (test==0))
    {
      acc++;
      logLx = logLy;
      logpx = logpy;
      lambdax = lambday;
      Nx = Ny;
      for(i=1; i<= Nf; i++)
      {
        sdatax[i] = sdatay[i];
        activex[i] = activey[i];
      }
      
      if(fixedD == 1)
      {
        if(logLx > logLmax)
        {
          logLmax=logLx;
          i = 0;
          for(ii=1; ii<= Nf; ii++)
          {
            if(activex[ii] == 1)  // only use active points
            {
              i++;
              tpoints[i] = spoints[ii];
              tdata[i] = sdatax[ii];
            }
          }        }
      }
    }
    if(mc%(N/100)==0)printProgress((double)mc/(double)N);
  }
  printProgress(1);
  printf("\n");
  
  if(fixedD == 0)
  {
  i = 0;
  for(ii=1; ii<= Nf; ii++)
  {
    if(activex[ii] == 1)  // only use active points
    {
      i++;
      tpoints[i] = spoints[ii];
      tdata[i] = sdatax[ii];
    }
  }
  }
  
  spline(tpoints, tdata, Nx, 1.0e31, 1.0e31, sderiv);
  
  for(i=imin; i <= imax; i++)
  {
    f = (double)(i)/T;
    splint(tpoints, tdata, sderiv, Nx, log(f), &model);
    Xnoise[i] = exp(model)*1.0e-40;
  }
  
  
  free_dvector(mdl,0,ND);
  free_dvector(ref,1,Nf);
  free_dvector(spoints,1,Nf);
  free_dvector(tdata,1,Nf);
  free_dvector(tpoints,1,Nf);
  free_dvector(sdatax,1,Nf);
  free_dvector(sdatay,1,Nf);
  free_dvector(sderiv,1,Nf);
  free_ivector(activex,1,Nf);
  free_ivector(activey, 1,Nf);
  
  
}

void splint(double *xa, double *ya, double *y2a, int n, double x, double *y)
{
  int klo,khi,k;
  double h,b,a;
  
  klo=1;
  khi=n;
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

#define NRANSI
void spline(double *x, double *y, int n, double yp1, double ypn, double *y2)
{
  int i,k;
  double p,qn,sig,un;
  double *u;
  
  u=dvector(1,n-1);
  if (yp1 > 0.99e30)
    y2[1]=u[1]=0.0;
  else {
    y2[1] = -0.5;
    u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
  }
  for (i=2;i<=n-1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)   qn=un=0.0;
  
  else {
    qn=0.5;
    un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
  }
  y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
  
  for (k=n-1;k>=1;k--) y2[k]=y2[k]*y2[k+1]+u[k];
  
  free_dvector(u,1,n-1);
}
#undef NRANSI






/*****************************************************/
/*                                                   */
/*         tanh-based Confusion Noise Fitting        */
/*                                                   */
/*****************************************************/

double confusion_fit(double f, double logA, double alpha, double beta, double kappa, double gamma, double fk)
{
  return pow(f,-5./3.) * exp( logA - pow(f,alpha) + beta*f*sin(kappa*f)) * (1.0 + tanh(gamma*(fk - f)));
}

void confusion_mcmc(double *data, double *noise, double *conf, int imin, int imax, double T)
{
  printf(" Parameterized fit to AE-channel confusion noise\n");
  double logA = log(1e-45);
  double alpha = 0.138;
  double beta = -221;
  double kappa = 521;
  double gamma = 1680;
  double fk = 0.00113;

  double logA_y = log(1e-45);
  double alpha_y = 0.1;
  double beta_y = 200;
  double kappa_y = 500;
  double gamma_y = 1000;
  double fk_y = 0.001;

  
  double logL,logL_y,logLmax;
  
  double *Sc   = malloc(imax*sizeof(double));
  double *Sc_y = malloc(imax*sizeof(double));

  
  long seed = -123456;
  
  //gasdev2(&seed)
  //ran2(&seed)
  
  
  logL = 0.0;
  for(int i=imin; i<imax; i++) Sc[i] = confusion_fit((double)i/T, logA, alpha, beta, kappa, gamma, fk);
  for(int i=imin; i<imax; i++) logL += -0.5*data[i]/(noise[i]+Sc[i]) - 0.5*log(Sc[i]+noise[i]);
  logLmax = logL;
  
  
//  FILE *testfile = fopen("testchain.dat","w");
//  FILE *testfile2;

  double fk_min = (double)imin/T;
  double fk_max = (double)imax/T;
  
  double gamma_min = 900;
  double gamma_max = 2000;
  
  double beta_min = -400;
  double beta_max = 400;
  
  int N = 50000;
  for(int n=0; n<N; n++)
  {
    if(n%(N/100)==0)printProgress((double)n/(double)N);

//    fprintf(testfile,"%lg %lg %lg %lg %lg %lg %lg\n",logL-(double)imax,logA, alpha,beta,kappa,gamma,fk);
//    fflush(testfile);
//
//    if(n%100==0)
//    {
//      testfile2 = fopen("testfit.dat","w");
//
//      for(int i=imin; i<imax; i++)
//      {
//        fprintf(testfile2,"%lg %lg %lg %lg\n",(double)i/T, data[i], noise[i], Sc[i]);
//      }
//
//      fclose(testfile2);
//    }


    double scale = 1.0;
    double draw = ran2(&seed);
    if(draw<0.3) scale = 10.0;
    else if (draw<0.6) scale = 1.0;
    else scale = 0.1;
    
    logA_y  = logA  + gasdev2(&seed)*1.0*scale;
    alpha_y = alpha + gasdev2(&seed)*0.01*scale;
    beta_y  = beta  + gasdev2(&seed)*100*scale;
    kappa_y = kappa + gasdev2(&seed)*10*scale;
    gamma_y = gamma + gasdev2(&seed)*200*scale;
    fk_y    = fk    + gasdev2(&seed)*.0001*scale;
    
    //check priors
    if(fk_y < fk_min || fk_y > fk_max) continue;
    if(gamma_y < gamma_min || gamma_y > gamma_max) continue;
    if(beta_y < beta_min || beta_y > beta_max) continue;

    
    //if still in the loop, calculate the likelihood
    logL_y = 0.0;
    for(int i=imin; i<imax; i++) Sc_y[i] = confusion_fit((double)i/T, logA_y, alpha_y, beta_y, kappa_y, gamma_y, fk_y);
    for(int i=imin; i<imax; i++)
    {
      logL_y += -0.5*data[i]/(noise[i]+Sc_y[i]) - 0.5*log(Sc_y[i]+noise[i]);
    }
    
    
    if(logL_y - logL > log(ran2(&seed)))
    {
      logL  = logL_y;
      logA  = logA_y;
      alpha = alpha_y;
      beta  = beta_y;
      kappa = kappa_y;
      gamma = gamma_y;
      fk    = fk_y;
      
      for(int i=imin; i<imax; i++) Sc[i] = Sc_y[i];
    }
    
    //store maxL
    if(logL>logLmax)
    {
      logLmax = logL;
      for(int i=imin; i<imax; i++) conf[i] = Sc[i];
    }
    
    
  }

  
  for(int i=imin; i<imax; i++) noise[i] += Sc[i];

  printProgress(1);
  printf("\n");

  
  free(Sc);
  free(Sc_y);
  
  
}



