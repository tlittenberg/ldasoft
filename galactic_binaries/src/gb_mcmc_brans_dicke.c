/*
 gcc gb_mcmc_chirpmass.c -lm -lgsl -o gb_mcmc_chirpmass
 */

/***************************  REQUIRED LIBRARIES  ***************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

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

double reduced_mass_ratio(double m1, double m2)
{
  return m1*m2/pow(m1+m2,2);
}

void component_masses(double Mc, double eta, double *m1, double *m2)
{
  double M = Mc*pow(eta,-3./5.);
  double mass1 = 0.5*(M + sqrt(M*M - 4.*M*M*eta));
  double mass2 = M - mass1;
  *m1 = mass1;
  *m2 = mass2;
}

double sensitivity(double m1, double m2, double r1, double r2)
{
  double s1 = m1*TSUN / r1*RSUN;
  double s2 = m2*TSUN / r2*RSUN;
  return pow(s1-s2,2);
}

double white_dwarf_radius(double M)
{
  double mu_e = 2.00;
  double M_ch = 1.44;
  
  return 0.0126*pow(2./mu_e,5./3.)*pow(M,-1./3.)*pow(1.-pow(M/M_ch,4./3.),1./2.);
}

int main(int argc, char* argv[])
{

  if(argc!=1)
  {
    fprintf(stdout,"Usage: gb_mcmc_chirpmass\n");
    return 1;
  }
  
  FILE *ifile;
  FILE *ofile;
  
  //set RNG for noise
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *r = gsl_rng_alloc(T);
  gsl_rng_env_setup();
  gsl_rng_set (r, atoi(argv[2]));
  
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
  
  int count_pos  = 0;
  int count_neg  = 0;
  int count_real = 0;
  int count_imag = 0;
  
  double m1,sm1;
  double m2,sm2;
  double r1,sr1;
  double r2,sr2;
  
  double mass_1;
  double mass_2;
  double radius_1;
  double radius_2;
  double eta;
  double S2;
  
  double omega_BD;

  //get injection parameters
  ifile = fopen("injection_parameters_0_0.dat","r");
  fscanf(ifile,"%i%lg%lg%lg%lg%lg%lg%lg%lg%lg",&i,&f,&fdot,&A,&phi,&theta,&cosi,&psi,&phase,&fddot);
  Mc  = M_fdot(f, fdot);
  eta = 0.24;
  component_masses(Mc,eta,&m1,&m2);
  r1 = white_dwarf_radius(m1);
  r2 = white_dwarf_radius(m2);
  
  //fix uncertainty on mass and radius
  sm1 = sm2 = 0.01;
  sr1 = sr2 = 0.01;
  
  ifile = fopen("chains/parameter_chain.dat.0","r");
  ofile = fopen("brans_dicke.dat","w");
  
  while(!feof(ifile))
  {
    fscanf(ifile,"%lg%lg%lg%lg%lg%lg%lg%lg%lg",&f,&fdot,&A,&phi,&theta,&cosi,&psi,&phase,&fddot);
    if(fdot>0.0&&fddot>0.0)
    {
      count_pos++;
      Mc = M_fdot_fddot(f,fdot,fddot);
      B = beta(f,fdot,fddot);
      if(Mc==Mc && B==B)
      {
        count_real++;
        
        for(int k=0; k<100; k++)
        {
          //draw M1
          mass_1 = m1 + gsl_ran_gaussian(r,sm1);
          
          //draw R1
          radius_1 = r1 + gsl_ran_gaussian(r,sr1);
          
          //draw M2
          mass_2 = m2 + gsl_ran_gaussian(r,sm2);
          
          //draw R2
          radius_2 = r2 + gsl_ran_gaussian(r,sr1);
          
          //get omega_BD
          eta = reduced_mass_ratio(mass_1, mass_2);
          S2  = sensitivity(mass_1, mass_2, radius_2, radius_2);
          
          omega_BD = (5./48.) * S2 * pow(eta,2./5.)/B;
          
          fprintf(ofile,"%.12g %.12g %.12g %.12g %.12g %.12g\n",Mc, B, eta, S2, omega_BD, 1./omega_BD);
        }
      }
      else
      {
        count_imag++;
      }
    }
    count_neg++;
  }
  printf("Drawn masses & radii: M1 = %g +/- %g, R1 = %g +/- %g\n",m1,sm1,r1,sr1);
  printf("                      M2 = %g +/- %g, R2 = %g +/- %g\n",m2,sm2,r2,sr2);
  printf("Bayes factor for positive frequency evolution: %g\n",(double)count_pos/(double)count_neg);
  
  return 0;
}

