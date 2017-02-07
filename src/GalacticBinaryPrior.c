//
//  GalacticBinaryPrior.c
//  
//
//  Created by Littenberg, Tyson B. (MSFC-ZP12) on 2/6/17.
//
//
#include <math.h>

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
  
  //TODO: assign priors by parameter name, use mapper to get into vector (more robust to changes)
  
  //frequency bin
  model->prior[0][0] = data->qmin;
  model->prior[0][1] = data->qmax;
  
  //colatitude
  model->prior[1][0] = -1.0;
  model->prior[1][1] =  1.0;
  
  //longitude
  model->prior[2][0] = 0.0;
  model->prior[2][1] = 2.0*M_PI;
  
  //log amplitude
  model->prior[3][0] = -55.0;
  model->prior[3][1] = -47.0;
  
  //cos inclination
  model->prior[4][0] = -1.0;
  model->prior[4][1] =  1.0;

  //polarization
  model->prior[5][0] = 0.0;
  model->prior[5][1] = M_PI/4.0;

  //phase
  model->prior[6][0] = 0.0;
  model->prior[6][1] = 2.0*M_PI;
  
  //fdot (bins/Tobs)
  model->prior[7][0] = -10.0;
  model->prior[7][1] = +10.0;
  
  //set prior volume
  model->logPriorVolume = 0.0;
  for(int n=0; n<8; n++) model->logPriorVolume -= log(model->prior[n][1]-model->prior[n][0]);
  
}

double evaluate_uniform_prior(struct Model *model, double *params)
{
  double **prior = model->prior;

  for(int i=0; i<8; i++)
  {
    if(params[i]<prior[i][0] || params[i]>prior[i][1]) return -1.0e-60;
  }
  
  return model->logPriorVolume;
}
