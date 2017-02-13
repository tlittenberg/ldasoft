
/***************************  REQUIRED LIBRARIES  ***************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "omp.h"

/*************  PROTOTYPE DECLARATIONS FOR INTERNAL FUNCTIONS  **************/

#include "LISA.h"
#include "GalacticBinary.h"
#include "GalacticBinaryIO.h"
#include "GalacticBinaryData.h"
#include "GalacticBinaryPrior.h"
#include "GalacticBinaryModel.h"
#include "GalacticBinaryProposal.h"
#include "GalacticBinaryWaveform.h"

void ptmcmc(struct Model **model, struct Chain *chain);

/* ============================  MAIN PROGRAM  ============================ */

int main(int argc, char *argv[])
{
  
  time_t start, stop;
  start = time(NULL);
  
  int NC = 12; //number of chains
  int Nmax=10; //max number of waveforms

  /* Allocate data structures */
  struct Flags *flags  = malloc(sizeof(struct Flags));
  struct Data  *data   = malloc(sizeof(struct Data));
  struct Orbit *orbit  = malloc(sizeof(struct Orbit));
  struct Chain *chain  = malloc(sizeof(struct Chain));
  struct Model **model = malloc(sizeof(struct Model*)*NC);
  struct Model **trial = malloc(sizeof(struct Model*)*NC);
  
  
  /* Parse command line and set defaults/flags */
  parse(argc,argv,data,orbit,flags);

  
  /* Load spacecraft ephemerides */
  initialize_orbit(orbit);

  
  /* Initialize data structures */
  data->tdi = malloc(sizeof(struct TDI));
  alloc_tdi(data->tdi, data->N, data->Nchannel);
  data->noise = malloc(sizeof(struct Noise));
  alloc_noise(data->noise, data->N);

  
  /* Inject gravitational wave signal */
  if(flags->injection) GalacticBinaryInjectVerificationSource(data,orbit,flags);

  /* Initialize parallel chain */
  initialize_chain(chain, &data->cseed, NC);
  
  /* Initialize data models */
  int ic;
#pragma omp parallel for private(ic) shared(model,chain,data,orbit,trial)
  for(ic=0; ic<NC; ic++)
  {
    model[ic] = malloc(sizeof(struct Model));
    alloc_model(model[ic],Nmax,data->N,data->Nchannel);

    trial[ic] = malloc(sizeof(struct Model));
    alloc_model(trial[ic],Nmax,data->N,data->Nchannel);

    set_uniform_prior(model[ic], data);
    
    for(int n=0; n<Nmax; n++)
    {
      draw_from_prior(model[ic], model[ic]->source[n]->params, chain->r[ic]);
      galactic_binary_fisher(orbit, data, model[ic]->source[n], data->noise);
      map_array_to_params(model[ic]->source[n], model[ic]->source[n]->params, data->T);
    }
    
    model[ic]->logL = gaussian_log_likelihood(orbit, data, model[ic]);
  }

  FILE *temp=fopen("chain.dat","w");
  for(int mcmc = 0; mcmc < 1000; mcmc++)
  {
#pragma omp parallel for private(ic) shared(model,chain,data,orbit,trial)
    for(ic=0; ic<NC; ic++)
    {
      for(int steps=0; steps < 100; steps++)
      {
        
        copy_model(model[ic],trial[ic]);
        
        if(gsl_rng_uniform(chain->r[ic]) < 0.1) draw_from_prior(model[ic], trial[ic]->source[0]->params, chain->r[ic]);
        
        else draw_from_fisher(model[ic], model[ic]->source[0], trial[ic]->source[0]->params, chain->r[ic]);
        
        map_array_to_params(trial[ic]->source[0], trial[ic]->source[0]->params, data->T);
        
        trial[ic]->logL = gaussian_log_likelihood(orbit, data, trial[ic]);
        
        double H = (trial[ic]->logL - model[ic]->logL)/chain->temperature[ic]; //hastings ratio
        double a = log(gsl_rng_uniform(chain->r[ic]));                         //transition probability
        
        if(H > a) copy_model(trial[ic],model[ic]);
      }
    }
    ptmcmc(model,chain);
    fprintf(temp, "%lg %lg %lg\n",model[chain->index[0]]->logL, model[chain->index[0]]->source[0]->params[0],model[chain->index[0]]->source[0]->params[3]);
  }
  
  //print total run time
  stop = time(NULL);
  
  if(flags->verbose) printf(" ELAPSED TIME = %g second\n",(double)(stop-start));
  
  
  //free memory and exit cleanly
  for(ic=0; ic<NC; ic++)
  {
    free_model(model[ic]);
    free_model(trial[ic]);
  }
  free_orbit(orbit);
  free_noise(data->noise);
  free_tdi(data->tdi);
  free(model);
  free(trial);
  free(chain);
  free(data);
  
  return 0;
}

void ptmcmc(struct Model **model, struct Chain *chain)
{
  int a, b;
  int olda, oldb;
  
  double heat1, heat2;
  double logL1, logL2;
  double dlogL;
  double H;
  double alpha;
  double beta;
  
  int NC = chain->NC;
  
  //b = (int)(ran2(seed)*((double)(chain->NC-1)));
  for(b=NC-1; b>0; b--)
  {
    a = b - 1;
    //chain[a]->A=0;
    
    olda = chain->index[a];
    oldb = chain->index[b];
    
    heat1 = chain->temperature[a];
    heat2 = chain->temperature[b];
    
    logL1 = model[olda]->logL;
    logL2 = model[oldb]->logL;
    
    //Hot chains jump more rarely
    if(gsl_rng_uniform(chain->r[a])<1.0)
    {
      dlogL = logL2 - logL1;
      H  = (heat2 - heat1)/(heat2*heat1);
      
      alpha = exp(dlogL*H);
      beta  = gsl_rng_uniform(chain->r[a]);
      
      //chain->ptprop[a]++;
      
      if(alpha >= beta)
      {
        //chain->ptacc[a]++;
        chain->index[a] = oldb;
        chain->index[b] = olda;
        //chain->A[a]=1;
      }
    }
  }
}


