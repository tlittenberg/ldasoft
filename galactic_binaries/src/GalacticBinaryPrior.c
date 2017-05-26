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
  double qpad = 10;
  model->prior[0][0] = data->qmin+qpad;
  model->prior[0][1] = data->qmax-qpad;
  
  //colatitude
  model->prior[1][0] = -1.0;
  model->prior[1][1] =  1.0;
  
  //longitude
  model->prior[2][0] = 0.0;
  model->prior[2][1] = PI2;
  
  //log amplitude
  model->prior[3][0] = -55.0;
  model->prior[3][1] = -45.0;
  
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
  while(params[1] < prior[1][0] || params[1] > prior[1][1])
  {
    if(params[1] < prior[1][0] ) params[1] = 2.0*prior[1][0] - params[1];
    if(params[1] > prior[1][1] ) params[1] = 2.0*prior[1][1] - params[1];
  }
  
  //longitude (periodic)
  while(params[2] < prior[2][0]) params[2] += prior[2][1]-prior[2][0];
  while(params[2] > prior[2][1]) params[2] -= prior[2][1]-prior[2][0];
  
  //log amplitude (step)
  if(params[3]<prior[3][0] || params[3]>prior[3][1]) return -INFINITY;
  
  //cosine inclination (reflective)
  while(params[4] < prior[4][0] || params[4] > prior[4][1])
  {
    if(params[4] < prior[4][0] ) params[4] = 2.0*prior[4][0] - params[4];
    if(params[4] > prior[4][1] ) params[4] = 2.0*prior[4][1] - params[4];
  }
  
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
