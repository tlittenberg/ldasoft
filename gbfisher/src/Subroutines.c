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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_randist.h>

#include "arrays.h"
#include "Detector.h"
#include "Subroutines.h"

#include <LISA.h>
#include <Constants.h>
#include <GalacticBinary.h>
#include <GalacticBinaryIO.h>

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60



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

  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *seed = gsl_rng_alloc(T);
  gsl_rng_env_setup();

  
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
  
  XX = double_vector(100);
  
  rseed = -546214;
  
  segs = (int)((double)(imax-imin)/101.0);
  
  lfmin = log((double)(imin-101)/TOBS);
  lfmax = log((double)(imin+101*(segs))/TOBS);
  
  
  Npoly = 30;
  
  dlf = (lfmax-lfmin)/(double)(Npoly);
  ln4 = log(1.0e-4);
  
  fdata = double_vector(segs-1);
  mdata = double_vector(segs-1);
  inst = double_vector(segs-1);
  pcx = double_vector(Npoly);
  pcy = double_vector(Npoly);
  
  for(i=0; i < segs; i++)
  {
    for(j=0; j<=100; j++) XX[j] = XP[imin+101*i+j];
    f = (double)(imin+101*i-50)/TOBS;
    SAE  = AEnoise(L,fstar,f);
    SXYZ = XYZnoise(L,fstar,f);
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
    j = (long)((double)(Npoly+1)*gsl_rng_uniform(seed));
    mul = 1.0;
    alpha = gsl_rng_uniform(seed);
    if(alpha > 0.3) mul = 10.0;
    if(alpha > 0.7) mul = 0.01;
    pcy[j] = pcx[j]+mul*0.01*gsl_ran_gaussian(seed,1)/sqrt(beta);
    
    chiy = 0.0;
    for(i=0; i < segs; i++)
    {
      lf = log((double)(imin+101*i-50)/TOBS);
      j = (long)floor((lf-lfmin)/dlf);
      fit = pcy[j]+((pcy[j+1]-pcy[j])/dlf)*(lf-(lfmin+(double)(j)*dlf));
      f = exp(-1.0*((lfmin+((double)(j)+0.5)*dlf)-ln4));
      chiy += 500.0*(mdata[i] - fit)*(mdata[i] - fit)*f;
    }
    
    alpha = log(gsl_rng_uniform(seed));
    
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
    SAE  = AEnoise(L,fstar,f);
    SXYZ = XYZnoise(L,fstar,f);
    alpha = exp(fit)*1.0e-40;
    conf = alpha -SXYZ;
    if(conf < SXYZ/30.0) conf = 1.0e-46;
    Xnoise[i] = alpha;
    Xconf[i] = conf;
    
    
    
    //Xnoise[i]*=pow(sin(f/fstar),2.0);
    //Xconf[i]*=pow(sin(f/fstar),2.0);
  }
  
  free_double_vector(XX);
  free_double_vector(fdata);
  free_double_vector(mdata);
  free_double_vector(pcx);
  free_double_vector(pcy);
  
  return;
}

void medianAE(long imin, long imax, double fstar, double L, double *AEP, double *AEinst, double *AEconf, double TOBS)
{
  printf(" Median fit to AE-channel confusion noise\n");
  
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *seed = gsl_rng_alloc(T);
  gsl_rng_env_setup();

  double f;
  double SAE, SXYZ;
  double chi;
  long i, j, k;
  long segs;
  int Npoly;
  double *XX;
  double *fdata, *mdata, *pcx, *pcy, *inst;
  double chix, chiy, fit, alpha, beta, mul, conf;
  double lfmin, lfmax, dlf, lf, ln4;
  FILE *Xfile;
  
  XX = double_vector(100);
    
  segs = (int)((double)(imax-imin)/101.0);
  
  lfmin = log((double)(imin-101)/TOBS);
  lfmax = log((double)(imin+101*(segs))/TOBS);
  
  
  Npoly = 30;
  
  dlf = (lfmax-lfmin)/(double)(Npoly);
  ln4 = log(1.0e-4);
  
  fdata = double_vector(segs-1);
  mdata = double_vector(segs-1);
  inst = double_vector(segs-1);
  pcx = double_vector(Npoly);
  pcy = double_vector(Npoly);
  
  for(i=0; i < segs; i++)
  {
    for(j=0; j<=100; j++) XX[j] = AEP[imin+101*i+j];
    f = (double)(imin+101*i-50)/TOBS;
    SAE  = AEnoise(L,fstar,f);
    SXYZ = XYZnoise(L,fstar,f);
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
    j = (long)((double)(Npoly+1)*gsl_rng_uniform(seed));
    mul = 1.0;
    alpha = gsl_rng_uniform(seed);
    if(alpha > 0.3) mul = 10.0;
    if(alpha > 0.7) mul = 0.01;
    pcy[j] = pcx[j]+mul*0.01*gsl_ran_gaussian(seed,1)/sqrt(beta);
    
    chiy = 0.0;
    for(i=0; i < segs; i++)
    {
      lf = log((double)(imin+101*i-50)/TOBS);
      j = (long)floor((lf-lfmin)/dlf);
      fit = pcy[j]+((pcy[j+1]-pcy[j])/dlf)*(lf-(lfmin+(double)(j)*dlf));
      f = exp(-1.0*((lfmin+((double)(j)+0.5)*dlf)-ln4));
      chiy += 500.0*(mdata[i] - fit)*(mdata[i] - fit)*f;
    }
    
    alpha = log(gsl_rng_uniform(seed));
    
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
    SAE  = AEnoise(L,fstar,f);
    SXYZ = XYZnoise(L,fstar,f);
    alpha = exp(fit)*1.0e-40;
    conf = alpha -SAE;
    if(conf < SAE/30.0) conf = 1.0e-46;
    AEinst[i] = alpha;
    AEconf[i] = conf;
    
  }
  
  free_double_vector(XX);
  free_double_vector(fdata);
  free_double_vector(mdata);
  free_double_vector(pcx);
  free_double_vector(pcy);
  
  return;
}

void KILL(char* Message)
{
  printf("\a\n");
  printf("%s",Message);
  printf("Terminating the program.\n\n\n");
  exit(1);
}

double quickselect(double *arr, int n, int k)
{
  double *arr_sorted = malloc(n*sizeof(double));
  for(int i=0; i<n; i++) arr_sorted[i] = arr[i];
  
  double val = gsl_stats_select(arr_sorted, 1, n, k);
  
  free(arr_sorted);
  return val;
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
  
  XX = double_vector(divs);

  segs = (int)((double)(imax-imin)/(double)(divsp));
  
  lfmin = log((double)(imin)/T);
  lfmax = log((double)(imin+divs*(segs))/T);
  
  Npoly = 20;
  
  dlf = (lfmax-lfmin)/(double)(Npoly);
  ln4 = log(1.0e-4);
  
  fdata = double_vector(segs-1);
  mdata = double_vector(segs-1);
  adata = double_vector(segs-1);
  sdata = double_vector(segs-1);
  inst = double_vector(segs-1);
  pcx = double_vector(Npoly);
  pcy = double_vector(Npoly);
  
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
    SAE  = AEnoise(L,fstar,f);
    SXYZ = XYZnoise(L,fstar,f);
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
    SAE  = AEnoise(L,fstar,f);
    SXYZ = XYZnoise(L,fstar,f);
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
  
  const gsl_rng_type *RNGT = gsl_rng_default;
  gsl_rng *seed = gsl_rng_alloc(RNGT);
  gsl_rng_env_setup();
  
  gsl_spline *cspline;
  gsl_interp_accel *splineacc = gsl_interp_accel_alloc();


  N = 100000;   // number of MCMC steps
  
  Nf = 20; // maximum number of spline control points
    
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
  
  mdl = double_vector(ND);
  ref = double_vector(Nf);
  spoints = double_vector(Nf);
  tdata = double_vector(Nf);
  tpoints = double_vector(Nf);
  sdatax = double_vector(Nf);
  sdatay = double_vector(Nf);
  sderiv = double_vector(Nf);
  
  activex = int_vector(Nf);
  activey = int_vector(Nf);
  
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
  
  cspline = gsl_spline_alloc(gsl_interp_cspline, Nx);
  gsl_spline_init(cspline,tpoints,tdata,Nx);
  
  av = 0.0;
  for(j=0; j< ND; j++)
  {
    model = gsl_spline_eval(cspline,datax[j],splineacc);
    mdl[j] = model;
    y = (datay[j]-model)/sigma[j];
    av -= y*y;
  }
  logLx = av/2.0;
  
  gsl_spline_free (cspline);

  
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
    
    alpha = gsl_rng_uniform(seed);
    
    q = 0.5;
    if(fixedD == 1) q = 10.0;
    
    if(alpha > q)   // propose a dimension change
    {
      // Note that the spline points at the ends are never added or subtracted
      
      
      alpha = gsl_rng_uniform(seed);
      
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
            i = 2 + (int)(gsl_rng_uniform(seed)*(double)(Nf-2)); // pick one to kill
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
            i = 2 + (int)(gsl_rng_uniform(seed)*(double)(Nf-2)); // pick one to add
          } while(activex[i] == 1);  // can't add if already in use
          activey[i] = 1;
          
          sdatay[i] = rmin + (rmax-rmin)*gsl_rng_uniform(seed);  // draw from prior
          
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
      
      alpha = gsl_rng_uniform(seed);
      
      if(alpha > 0.6)  // update all points
      {
        
        for(ii=1; ii<= Nf; ii++)
        {
          // variety of jump sizes
          if(alpha > 0.8)
          {
            sdatay[ii] += sh*1.0e-1*gsl_ran_gaussian(seed,1);
          }
          else if (alpha > 0.5)
          {
            sdatay[ii] += sh*1.0e-2*gsl_ran_gaussian(seed,1);
          }
          else if (alpha > 0.3)
          {
            sdatay[ii] += sh*1.0e-3*gsl_ran_gaussian(seed,1);
          }
          else
          {
            sdatay[ii] += sh*1.0e-4*gsl_ran_gaussian(seed,1);
          }
          
        }
        
      }
      else  if(alpha > 0.1) // just update one
      {
        
        do
        {
          ii = (int)(gsl_rng_uniform(seed)*(double)(Nf));
        }while(activey[ii] == 0);
        
        sdatay[ii] += sh*1.0e-1*gsl_ran_gaussian(seed,1);
        
      }
      else  // birth/death
      {
        do
        {
          i = 2 + (int)(gsl_rng_uniform(seed)*(double)(Nf-2)); // pick one to kill
        } while(activex[i] == 0);  // can't kill it if already dead
        activey[i] = 0;
        
        do
        {
          i = 2 + (int)(gsl_rng_uniform(seed)*(double)(Nf-2)); // pick one to add
        } while(activey[i] == 1);  // can't add if already in use
        activey[i] = 1;
        
        sdatay[i] = rmin + (rmax-rmin)*gsl_rng_uniform(seed);  // draw from prior
        
        
      }
      
      
    }
    
    
    alpha = gsl_rng_uniform(seed);
    if(alpha > 0.9)
    {
      lambday = lmin+(lmax-lmin)*gsl_rng_uniform(seed);  // uniform draw from the prior
    }
    else if (alpha > 0.7)
    {
      lambday = lambdax + 0.1*gsl_ran_gaussian(seed,1);
    }
    else if (alpha > 0.3)
    {
      lambday = lambdax + 0.01*gsl_ran_gaussian(seed,1);
    }
    else
    {
      lambday = lambdax + 0.001*gsl_ran_gaussian(seed,1);
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
        
        cspline = gsl_spline_alloc(gsl_interp_cspline, Ny);
        gsl_spline_init(cspline,tpoints,tdata,Ny);
        
        av = 0.0;
        for(j=0; j< ND; j++)
        {
          model = gsl_spline_eval(cspline,datax[j],splineacc);

          mdl[j] = model;
          y = (datay[j]-model)/sigma[j];
          av -= y*y;
        }
        logLy = av/2.0;
        gsl_spline_free (cspline);

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
    
    alpha = log(gsl_rng_uniform(seed));
    
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
  
  
  cspline = gsl_spline_alloc(gsl_interp_cspline, Nx);
  gsl_spline_init(cspline,tpoints,tdata,Nx);
  
  for(i=imin; i <= imax; i++)
  {
    f = (double)(i)/T;
    model = gsl_spline_eval(cspline,log(f),splineacc);

    Xnoise[i] = exp(model)*1.0e-40;
  }
  gsl_spline_free (cspline);

  
  free_double_vector(mdl);
  free_double_vector(ref);
  free_double_vector(spoints);
  free_double_vector(tdata);
  free_double_vector(tpoints);
  free_double_vector(sdatax);
  free_double_vector(sdatay);
  free_double_vector(sderiv);
  free_int_vector(activex);
  free_int_vector(activey);
  
  gsl_interp_accel_free (splineacc);

}


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

      const gsl_rng_type *RNGT = gsl_rng_default;
      gsl_rng *seed = gsl_rng_alloc(RNGT);
      gsl_rng_env_setup();

  
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
    double draw = gsl_rng_uniform(seed);
    if(draw<0.3) scale = 10.0;
    else if (draw<0.6) scale = 1.0;
    else scale = 0.1;
    
    logA_y  = logA  + gsl_ran_gaussian(seed,1)*1.0*scale;
    alpha_y = alpha + gsl_ran_gaussian(seed,1)*0.01*scale;
    beta_y  = beta  + gsl_ran_gaussian(seed,1)*100*scale;
    kappa_y = kappa + gsl_ran_gaussian(seed,1)*10*scale;
    gamma_y = gamma + gsl_ran_gaussian(seed,1)*200*scale;
    fk_y    = fk    + gsl_ran_gaussian(seed,1)*.0001*scale;
    
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
    
    
    if(logL_y - logL > log(gsl_rng_uniform(seed)))
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



