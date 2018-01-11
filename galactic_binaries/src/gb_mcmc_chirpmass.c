/*
 gcc gb_mcmc_chirpmass.c -lm -lgsl -o gb_mcmc_chirpmass
 */

/***************************  REQUIRED LIBRARIES  ***************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "Constants.h"

#define RAD2ARCMIN 3437.75
#define RAD2DEGREE 57.295833313961

/* ============================  MAIN PROGRAM  ============================ */

double M_fdot(double f, double fdot)
{
  return pow( fdot*pow(f,-11./3.)*(5./96.)*pow(M_PI,-8./3.)  ,  3./5.)/TSUN;
}

double M_fddot(double f, double fddot)
{
  return pow( fddot*pow(f,-19./3.)*(25./33792.)*pow(M_PI,-16./3.)  ,  3./10.)/TSUN;
}

double M_fdot_fddot(double f, double fdot, double fddot)
{
  //Eq 174 in SVN:LISA-WDWD-Test/Mathematica/with-Tyson.nb.pdf
  return pow(5.,3./5.)*pow((-(fdot*sqrt(fddot))/(pow(f,19./6.)*(30.*sqrt(f*fddot) - 11.*sqrt(33.)*fdot))),3./5.)/8./pow(M_PI,8./5.)/TSUN;
}

double beta(double f, double fdot, double fddot)
{
  //Eq 175 in SVN:LISA-WDWD-Test/Mathematica/with-Tyson.nb.pdf
  return (11.*(-3*sqrt(f*fddot)+sqrt(33.)*fdot)*pow((sqrt(fddot)*fdot)/(-30.*sqrt(f*fddot) + 11.*sqrt(33.)*fdot),2./5.)*pow(5./M_PI,2./5.))/(12.*pow(f,11./10.)*sqrt(fddot));
}

double galactic_binary_dL(double f0, double dfdt, double A)
{
  double f    = f0;//T;
  double fd = dfdt;//(T*T);
  double amp   = A;
  return ((5./48.)*(fd/(M_PI*M_PI*f*f*f*amp))*C/PC); //seconds  !check notes on 02/28!
}

int main(int argc, char* argv[])
{

  if(argc!=3)
  {
    fprintf(stdout,"Usage: gb_mcmc_chirpmass parameter_chain.dat.0 9\n");
    return 1;
  }
  
  FILE *ifile = fopen(argv[1],"r");
  
  int NP = atoi(argv[2]);
  
  char filename[128];
  sprintf(filename,"mchirp_%s",argv[1]);
  FILE *ofile = fopen(filename,"w");
  
  sprintf(filename,"frequencies_%s",argv[1]);
  FILE *ofile2 = fopen(filename,"w");

  sprintf(filename,"localization_%s",argv[1]);
  FILE *ofile3 = fopen(filename,"w");

  //parse file
  
  /*
   0)  index
   1)  logL
   2)  t0
   3)  tgap
   4)  f
   5)  fdot
   6)  A
   7)  phi
   8)  theta
   9)  cosi
   10) psi
   11) phase
   12) fddot
   */
  
  int i,n;
  double logL,t0,f,fdot,A,phi,theta,cosi,psi,phase,fddot;
  double Mc,B;
  
  while(!feof(ifile))
  {
    if(NP==9)fscanf(ifile,"%lg%lg%lg%lg%lg%lg%lg%lg%lg",&f,&fdot,&A,&phi,&theta,&cosi,&psi,&phase,&fddot);
    if(NP==8)
    {
      fscanf(ifile,"%lg%lg%lg%lg%lg%lg%lg%lg",&f,&fdot,&A,&phi,&theta,&cosi,&psi,&phase);
      fddot = 11.0/3.0*fdot*fdot/f;
    }
    if(fdot>0.0&&fddot>0.0)
    {
      Mc = M_fdot_fddot(f,fdot,fddot);
      B = beta(f,fdot,fddot);
      fprintf(ofile3,"%.12g %.12g %.12g %.12g %.12g\n",phi*RAD2ARCMIN,theta*RAD2ARCMIN,galactic_binary_dL(f,fdot,A),acos(cosi)*RAD2DEGREE,M_fdot(f,fdot));
      fprintf(ofile2,"%.12g %.12g %.12g\n",f,fdot,fddot);
      if(Mc==Mc && B==B) fprintf(ofile,"%.12g %.12g %.12g %.12g\n",M_fdot(f,fdot),M_fddot(f,fddot),M_fdot_fddot(f,fdot,fddot),beta(f,fdot,fddot));
    }
  }
  
  return 0;
}

