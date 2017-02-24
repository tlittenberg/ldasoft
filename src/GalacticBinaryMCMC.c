
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
void adapt_temperature_ladder(struct Chain *chain, int mcmc);

void galactic_binary_mcmc(struct Orbit *orbit, struct Data *data, struct Model **model, struct Model **trial, struct Chain *chain, int ic);
void noise_model_mcmc(struct Orbit *orbit, struct Data *data, struct Model **model, struct Model **trial, struct Chain *chain, int ic);


/* ============================  MAIN PROGRAM  ============================ */

int main(int argc, char *argv[])
{
  
  time_t start, stop;
  start = time(NULL);
  
  int NC    = 8;    //number of chains
  int Nmax  = 10;   //max number of waveforms
  int NMCMC = 3000; //number of MCMC steps
  int NBURN = 3000; //number of burn-in steps

  
  /* Allocate data structures */
  struct Flags *flags  = malloc(sizeof(struct Flags));
  struct Orbit *orbit  = malloc(sizeof(struct Orbit));
  struct Chain *chain  = malloc(sizeof(struct Chain));
  struct Data  **data  = malloc(sizeof(struct Data*)*Nmax);
  struct Model **model = malloc(sizeof(struct Model*)*NC);
  struct Model **trial = malloc(sizeof(struct Model*)*NC);
  
  
  /* Parse command line and set defaults/flags */
  for(int i=0; i<Nmax; i++) data[i] = malloc(sizeof(struct Data));
  parse(argc,argv,data,orbit,flags,Nmax);
  
  
  /* Load spacecraft ephemerides */
  initialize_orbit(orbit);
  
  
  /* Initialize data structures */
  alloc_data(data, NMCMC, Nmax);

  
  /* Inject gravitational wave signal */
  GalacticBinaryInjectVerificationSource(data,orbit,flags);

  
  /* Initialize parallel chain */
  initialize_chain(chain, flags, &data[0]->cseed, NC);
  
  
  /* Initialize data models */
  int ic;
#pragma omp parallel for private(ic) shared(model,chain,data,orbit,trial)
  for(ic=0; ic<NC; ic++)
  {
    model[ic] = malloc(sizeof(struct Model));
    alloc_model(model[ic],Nmax,data[0]->N,data[0]->Nchannel);

    trial[ic] = malloc(sizeof(struct Model));
    alloc_model(trial[ic],Nmax,data[0]->N,data[0]->Nchannel);

    set_uniform_prior(model[ic], data[0]);
    
    //set noise model
    copy_noise(data[0]->noise, model[ic]->noise);

    //set signal model
    for(int n=0; n<Nmax; n++)
    {
      draw_from_prior(model[ic], model[ic]->source[n], model[ic]->source[n]->params, chain->r[ic]);
      galactic_binary_fisher(orbit, data[0], model[ic]->source[n], data[0]->noise);
      map_array_to_params(model[ic]->source[n], model[ic]->source[n]->params, data[0]->T);
    }
    
    // Form master model & compute likelihood of starting position
    generate_noise_model(data[0], model[ic]);
    generate_signal_model(orbit, data[0], model[ic]);

    model[ic]->logL     = gaussian_log_likelihood(orbit, data[0], model[ic]);
    model[ic]->logLnorm = gaussian_log_likelihood_constant_norm(data[0], model[ic]);

  }

  
  /* The MCMC loop */
  for(int mcmc = -NBURN; mcmc < NMCMC; mcmc++)
  {
#pragma omp parallel for private(ic) shared(model,chain,data,orbit,trial)
    for(ic=0; ic<NC; ic++)
    {
      for(int steps=0; steps < 100; steps++)
      {
        galactic_binary_mcmc(orbit, data[0], model, trial, chain, ic);
        noise_model_mcmc(orbit, data[0], model, trial, chain, ic);
      }
      
      //update fisher matrix for each chain
      for(int n=0; n<model[ic]->Nlive; n++)
      {
        galactic_binary_fisher(orbit, data[0], model[ic]->source[n], data[0]->noise);
      }
      
    }
    ptmcmc(model,chain);
    adapt_temperature_ladder(chain, mcmc+NBURN);

    print_chain_files(data[0], model, chain, flags, mcmc);
    
    //store reconstructed waveform
    if(mcmc>0 && mcmc%data[0]->downsample==0)
    {
      print_chain_state(data[0], chain, model[chain->index[0]], stdout, mcmc);
      save_waveforms(data[0], model[chain->index[0]], mcmc/data[0]->downsample);
      for(ic=0; ic<NC; ic++) chain->avgLogL[ic] += model[chain->index[ic]]->logL + model[chain->index[ic]]->logLnorm;
    }
  }
  
  //print aggregate run files/results
  print_reconstructed_waveforms(data[0]);
  
  FILE *chainFile = fopen("avg_log_likelihood.dat","w");
  for(ic=0; ic<NC; ic++) fprintf(chainFile,"%lg %lg\n",1./chain->temperature[ic],chain->avgLogL[ic]/(double)(NMCMC/data[0]->downsample));
  fclose(chainFile);
  
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
  free_noise(data[0]->noise);
  free_tdi(data[0]->tdi);
  free_chain(chain,flags);
  free(model);
  free(trial);
  free(data[0]);
  
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
  for(b=1; b<NC; b++)
  {
    a = b - 1;
    chain->acceptance[a]=0;
    
    olda = chain->index[a];
    oldb = chain->index[b];
    
    heat1 = chain->temperature[a];
    heat2 = chain->temperature[b];
    
    logL1 = model[olda]->logL + model[olda]->logLnorm;
    logL2 = model[oldb]->logL + model[oldb]->logLnorm;
    
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
        chain->acceptance[a]=1;
      }
    }
  }
}

void adapt_temperature_ladder(struct Chain *chain, int mcmc)
{
  int ic;
  
  int NC = chain->NC;
  
  double S[NC];
  double A[NC][2];
  
  double nu=1;
  //double t0=100;
  double t0=1000.;
  
  for(ic=1; ic<NC-1; ic++)
  {
    S[ic] = log(chain->temperature[ic] - chain->temperature[ic-1]);
    A[ic][0] = chain->acceptance[ic-1];
    A[ic][1] = chain->acceptance[ic];
  }
  
  ic=0;
  for(ic=1; ic<NC-1; ic++)
  {
    S[ic] += (A[ic][0] - A[ic][1])*(t0/((double)mcmc+t0))/nu;
    //S[ic] += (A[ic][0] - A[ic][1])/nu;
    
    chain->temperature[ic] = chain->temperature[ic-1] + exp(S[ic]);
    
    if(chain->temperature[ic]/chain->temperature[ic-1] < 1.1) chain->temperature[ic] = chain->temperature[ic-1]*1.1;
  }//end loop over ic
}//end adapt function

void noise_model_mcmc(struct Orbit *orbit, struct Data *data, struct Model **model, struct Model **trial, struct Chain *chain, int ic)
{
  double logH  = 0.0; //(log) Hastings ratio
  double loga  = 1.0; //(log) transition probability
  
  double logPx  = 0.0; //(log) prior density for model x (current state)
  double logPy  = 0.0; //(log) prior density for model y (proposed state)

  //get right model pointer for current chain
  int nc = chain->index[ic];
  
  //shorthand pointers
  struct Model *model_x = model[nc];
  struct Model *model_y = trial[nc];
  
  copy_model(model_x,model_y);
  
  //choose proposal distribution
  switch(data->Nchannel)
  {
    case 1:
      model_y->noise->etaX = model_x->noise->etaX + 0.1*gsl_ran_gaussian(chain->r[ic],1);
      break;
    case 2:
      model_y->noise->etaA = model_x->noise->etaA + 0.1*gsl_ran_gaussian(chain->r[ic],1);
      model_y->noise->etaE = model_x->noise->etaE + 0.1*gsl_ran_gaussian(chain->r[ic],1);
      break;
  }
  
  //get priors for x and y
  switch(data->Nchannel)
  {
    case 1:
      if(model_y->noise->etaX < 0.1 || model_y->noise->etaX>10) logPy=-INFINITY;
      break;
    case 2:
      if(model_y->noise->etaA < 0.1 || model_y->noise->etaA>10.) logPy=-INFINITY;
      if(model_y->noise->etaE < 0.1 || model_y->noise->etaE>10.) logPy=-INFINITY;
      break;
  }
  
  if(logPy > -INFINITY)
  {
    //  Form master template
    generate_noise_model(data, model_y);
    
    //get likelihood for y
    model_y->logL     = gaussian_log_likelihood(orbit, data, model_y);
    model_y->logLnorm = gaussian_log_likelihood_constant_norm(data, model_y);
    
    /*
     H = [p(d|y)/p(d|x)]/T x p(y)/p(x) x q(x|y)/q(y|x)
     */
    logH += ( (model_y->logL+model_y->logLnorm) - (model_x->logL+model_x->logLnorm) )/chain->temperature[ic]; //delta logL
    logH += logPy  - logPx;                                         //priors
    
    loga = log(gsl_rng_uniform(chain->r[ic]));
    if(logH > loga) copy_model(model_y,model_x);
  }

}

void galactic_binary_mcmc(struct Orbit *orbit, struct Data *data, struct Model **model, struct Model **trial, struct Chain *chain, int ic)
{
  double logH  = 0.0; //(log) Hastings ratio
  double loga  = 1.0; //(log) transition probability
  
  double logPx  = 0.0; //(log) prior density for model x (current state)
  double logPy  = 0.0; //(log) prior density for model y (proposed state)
  double logQyx = 0.0; //(log) proposal denstiy from x->y
  double logQxy = 0.0; //(log) proposal density from y->x

  //get right model pointer for current chain
  int nc = chain->index[ic];

  //shorthand pointers
  struct Model *model_x = model[nc];
  struct Model *model_y = trial[nc];
  
  copy_model(model_x,model_y);
  
  //pick a source to update
  int n  = (int)(gsl_rng_uniform(chain->r[ic])*(double)model_x->Nlive);
  
  //more shorthand pointers
  struct Source *source_x = model_x->source[n];
  struct Source *source_y = model_y->source[n];
  
  //choose proposal distribution
  double draw = gsl_rng_uniform(chain->r[ic]);
  if(draw < 0.1)
  /* no proposal density terms because proposal is symmetric */
    draw_from_prior(model_x, source_y, source_y->params, chain->r[ic]);
  else if(draw<0.9)
  /* no proposal density terms because proposal is symmetric */
    draw_from_fisher(model_x, source_x, source_y->params, chain->r[ic]);
  else
  /* no proposal density terms because proposal is symmetric */
    fm_shift(data, model_x, source_x, source_y->params, chain->r[ic]);

  map_array_to_params(source_y, source_y->params, data->T);

  //get priors for x and y
  logPx = evaluate_uniform_prior(model_x, source_x->params);
  logPy = evaluate_uniform_prior(model_y, source_y->params);

  if(logPy > -INFINITY)
  {
    //  Form master template
    generate_signal_model(orbit, data, model_y);
    
    //get likelihood for y
    model_y->logL = gaussian_log_likelihood(orbit, data, model_y);
    
    /*
     H = [p(d|y)/p(d|x)]/T x p(y)/p(x) x q(x|y)/q(y|x)
     */
    logH += (model_y->logL - model_x->logL)/chain->temperature[ic]; //delta logL
    logH += logPy  - logPx;                                         //priors
    logH += logQxy - logQyx;                                        //proposals
    
    loga = log(gsl_rng_uniform(chain->r[ic]));
    if(logH > loga) copy_model(model_y,model_x);
  }
}


