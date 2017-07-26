//
//  GalacticBinaryPrior.c
//  
//
//  Created by Littenberg, Tyson B. (MSFC-ZP12) on 2/6/17.
//
//
#include <math.h>

#include "Constants.h"
#include "GalacticBinary.h"
#include "GalacticBinaryPrior.h"

void set_uniform_prior(struct Model *model, struct Data *data)
{
  /*
  params[0] = source->f0*T;
  params[1] = source->costheta;
  params[2] = source->phi;
  params[3] = log(source->amp);
  params[4] = source->cosi;
  params[5] = source->psi;
  params[6] = source->phi0;
  params[7] = source->dfdt*T*T;
  */
  
  
  //TODO:  make t0 a parameter
  for(int i=0; i<model->NT; i++)
  {
    model->t0[i] = data->t0[i];
    model->t0_min[i] = data->t0[i] - 20.0;
    model->t0_max[i] = data->t0[i] + 20.0;
  }
  
  //TODO: assign priors by parameter name, use mapper to get into vector (more robust to changes)
  
  //frequency bin
  //double qpad = 10;
  model->prior[0][0] = data->qmin;//+qpad;
  model->prior[0][1] = data->qmax;//-qpad;
  
  //colatitude
  model->prior[1][0] = -1.0;
  model->prior[1][1] =  1.0;
  
  //longitude
  model->prior[2][0] = 0.0;
  model->prior[2][1] = PI2;
  
  //log amplitude
  model->prior[3][0] = -51.0;//-54
  model->prior[3][1] = -53.0;
  
  //cos inclination
  model->prior[4][0] = -1.0;
  model->prior[4][1] =  1.0;

  //polarization
  model->prior[5][0] = 0.0;
  model->prior[5][1] = PI2;

  //phase
  model->prior[6][0] = 0.0;
  model->prior[6][1] = PI2;
  
  //fdot (bins/Tobs)
  
  /* frequency derivative priors are a little trickier...*/
  double fmin = model->prior[0][0]/data->T;
  double fmax = model->prior[0][1]/data->T;
  
  /* emprical envolope functions from Gijs' MLDC catalog */
  double fdotmin = -0.000005*pow(fmin,(13./3.));
  double fdotmax = 0.0000008*pow(fmax,(11./3.));
  
  double fddotmin = 11.0/3.0*fdotmin*fdotmin/fmax;
  double fddotmax = 11.0/3.0*fdotmax*fdotmax/fmin;
  
  if(data->NP>7)
  {
    model->prior[7][0] = fdotmin*data->T*data->T;
    model->prior[7][1] = fdotmax*data->T*data->T;
  }
  if(data->NP>8)
  {
    model->prior[8][0] = fddotmin*data->T*data->T*data->T;
    model->prior[8][1] = fddotmax*data->T*data->T*data->T;
  }

  //set prior volume
  model->logPriorVolume = 0.0;
  for(int n=0; n<data->NP; n++) model->logPriorVolume -= log(model->prior[n][1]-model->prior[n][0]);
  
}

double evaluate_uniform_prior(struct Model *model, double *params)
{
  double **prior = model->prior;
  
  //frequency bin (uniform)
  if(params[0]<prior[0][0] || params[0]>prior[0][1]) return -INFINITY;
  
  //colatitude (reflective)
  if(params[1]<prior[1][0] || params[1]>prior[1][1]) return -INFINITY;
//  while(params[1] < prior[1][0] || params[1] > prior[1][1])
//  {
//    if(params[1] < prior[1][0] ) params[1] = 2.0*prior[1][0] - params[1];
//    if(params[1] > prior[1][1] ) params[1] = 2.0*prior[1][1] - params[1];
//  }
  
  //longitude (periodic)
  while(params[2] < prior[2][0]) params[2] += prior[2][1]-prior[2][0];
  while(params[2] > prior[2][1]) params[2] -= prior[2][1]-prior[2][0];
  
  //log amplitude (step)
  if(params[3]<prior[3][0] || params[3]>prior[3][1]) return -INFINITY;
  
  //cosine inclination (reflective)
  if(params[4]<prior[4][0] || params[4]>prior[4][1]) return -INFINITY;
//  while(params[4] < prior[4][0] || params[4] > prior[4][1])
//  {
//    if(params[4] < prior[4][0] ) params[4] = 2.0*prior[4][0] - params[4];
//    if(params[4] > prior[4][1] ) params[4] = 2.0*prior[4][1] - params[4];
//  }
  
  //polarization
  while(params[5] < prior[5][0]) params[5] += prior[5][1]-prior[5][0];
  while(params[5] > prior[5][1]) params[5] -= prior[5][1]-prior[5][0];
  
  //phase
  while(params[6] < prior[6][0]) params[6] += prior[6][1]-prior[6][0];
  while(params[6] > prior[6][1]) params[6] -= prior[6][1]-prior[6][0];
  
  //fdot (bins/Tobs)
  if(params[7]<prior[7][0] || params[7]>prior[7][1]) return -INFINITY;

  return model->logPriorVolume;
}

static double loglike(double *x, int D)
{
  double u, rsq, z, s, ll;
  
  z = x[2];
  rsq = x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
  u = sqrt(x[0]*x[0]+x[1]*x[1]);
  
  s = 1.0/cosh(z/GALAXY_Zd);
  
  // Note that overall rho0 in density is irrelevant since we are working with ratios of likelihoods in the MCMC
  
  ll = log(GALAXY_A*exp(-rsq/(GALAXY_Rb*GALAXY_Rb))+(1.0-GALAXY_A)*exp(-u/GALAXY_Rd)*s*s);
  
  return(ll);
  
}


static void rotate_galtoeclip(double *xg, double *xe)
{
  xe[0] = -0.05487556043*xg[0] + 0.4941094278*xg[1] - 0.8676661492*xg[2];
  
  xe[1] = -0.99382137890*xg[0] - 0.1109907351*xg[1] - 0.00035159077*xg[2];
  
  xe[2] = -0.09647662818*xg[0] + 0.8622858751*xg[1] + 0.4971471918*xg[2];
}

void setup_galaxy_prior(struct Flags *flags, double *skyhist, int Nth, int Nph)
{
  double *x, *y;  // current and proposed parameters
  int D = 3;  // number of parameters
//  int Nth = 200;  // bins in cos theta
//  int Nph = 200;  // bins in phi
  int j;
  int ith, iph, cnt;
  double H, dOmega;
  double logLx, logLy;
  double alpha, beta, xx, yy, zz;
  double *xe, *xg;
  double r_ec, theta, phi;
//  double *skyhist;
  int mc;
  FILE *chain = NULL;
  
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  
  x = (double*)malloc(sizeof(double)* D);
  xe = (double*)malloc(sizeof(double)* D);
  xg = (double*)malloc(sizeof(double)* D);
  y = (double*)malloc(sizeof(double)* D);
  
  skyhist = (double*)malloc(sizeof(double)* (Nth*Nph));
  
  for(ith=0; ith< Nth; ith++)
  {
    for(iph=0; iph< Nph; iph++)
    {
      skyhist[ith*Nph+iph] = 0.0;
    }
  }
  
  // starting values for parameters
  x[0] = 0.5;
  x[1] = 0.4;
  x[2] = 0.1;
  
  logLx = loglike(x, D);
  
  cnt = 0;
  
  if(flags->verbose)chain = fopen("chain.dat", "w");
  
  for(mc=0; mc< 1000000000; mc++)
  {
    alpha = gsl_rng_uniform(r);
    
    if(alpha > 0.7)  // uniform draw from a big box
    {
      
      y[0] = 20.0*GALAXY_Rd*(-1.0+2.0*gsl_rng_uniform(r));
      y[1] = 20.0*GALAXY_Rd*(-1.0+2.0*gsl_rng_uniform(r));
      y[2] = 40.0*GALAXY_Zd*(-1.0+2.0*gsl_rng_uniform(r));
      
    }
    else
    {
      
      beta = gsl_rng_uniform(r);
      
      if(beta > 0.8)
      {
        xx = 0.1;
      }
      else if (beta > 0.4)
      {
        xx = 0.01;
      }
      else
      {
        xx = 0.001;
      }
      for(j=0; j< 2; j++) y[j] = x[j] + gsl_ran_gaussian(r,xx);
      y[2] = x[2] + 0.1*gsl_ran_gaussian(r,xx);
      
    }
    
    
    logLy = loglike(y, D);
    
    H = logLy - logLx;
    beta = gsl_rng_uniform(r);
    beta = log(beta);
    
    if(H > beta)
    {
      for(j=0; j< D; j++) x[j] = y[j];
      logLx = logLy;
    }
    
    if(mc%100 == 0)
    {
      
      /* rotate from galactic to ecliptic coordinates */
      
      xg[0] = x[0] - GALAXY_RGC;   // solar barycenter is offset from galactic center along x-axis (by convention)
      xg[1] = x[1];
      xg[2] = x[2];
      
      rotate_galtoeclip(xg, xe);
      
      r_ec = sqrt(xe[0]*xe[0]+xe[1]*xe[1]+xe[2]*xe[2]);
      
      theta = M_PI/2.0-acos(xe[2]/r_ec);
      
      phi = atan2(xe[1],xe[0]);
      
      if(phi<0.0) phi += 2.0*M_PI;
      
      if(mc%1000 == 0 && flags->verbose) fprintf(chain,"%d %e %e %e %e %e %e %e\n", mc/1000, logLx, x[0], x[1], x[2], theta, phi, r_ec);
      
      ith = (int)(0.5*(1.0+sin(theta))*(double)(Nth));
      iph = (int)(phi/(2.0*M_PI)*(double)(Nph));
      cnt++;
      
      if(ith < 0 || ith > Nth -1) printf("%d %d\n", ith, iph);
      if(iph < 0 || iph > Nph -1) printf("%d %d\n", ith, iph);
      
      skyhist[ith*Nph+iph] += 1.0;
      
    }
    
  }
  
  if(flags->verbose)fclose(chain);
  
  dOmega = 4.0*M_PI/(double)(Nth*Nph);
  
  double uni = 0.1;
  yy = (1.0-uni)/(double)(cnt);
  zz = uni/(double)(Nth*Nph);
  
  yy /= dOmega;
  zz /= dOmega;
  
  for(ith=0; ith< Nth; ith++)
  {
    for(iph=0; iph< Nph; iph++)
    {
      xx = yy*skyhist[ith*Nph+iph];
      //if(fabs(xx)  > 0.0) printf("%e %e %e\n", xx, zz, yy);
      skyhist[ith*Nph+iph] = xx + zz;
    }
  }
  
  if(flags->verbose)
  {
    chain = fopen("skyprior.dat", "w");
    for(ith=0; ith< Nth; ith++)
    {
      xx = -1.0+2.0*((double)(ith)+0.5)/(double)(Nth);
      xx *= -1.0;    // plotting flip?
      for(iph=0; iph< Nph; iph++)
      {
        yy = 2.0*M_PI*((double)(iph)+0.5)/(double)(Nph);
        yy = 2.0*M_PI - yy;  // plotting flip?
        fprintf(chain,"%e %e %e\n", yy, xx, skyhist[ith*Nph+iph]);
      }
      fprintf(chain,"\n");
    }
    fclose(chain);
  }
  
  free(x);
  free(y);
  free(xe);
  free(xg);
  gsl_rng_free (r);
  
}
