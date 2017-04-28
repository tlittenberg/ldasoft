
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

#define FIXME 0

void ptmcmc(struct Model ****model, struct Chain *chain, struct Flags *flags);
void adapt_temperature_ladder(struct Chain *chain, int mcmc);

void galactic_binary_mcmc(struct Orbit *orbit, struct Data **data, struct Model **model, struct Model **trial, struct Chain *chain, struct Flags *flags, int ic);
void galactic_binary_drmc(struct Orbit *orbit, struct Data **data, struct Model **model, struct Model **trial, struct Chain *chain, struct Flags *flags, int ic);
void data_mcmc(struct Orbit *orbit, struct Data ***data, struct Model ***model, struct Chain *chain, struct Flags *flags, int ic);
void noise_model_mcmc(struct Orbit *orbit, struct Data *data, struct Model *model, struct Model *trial, struct Chain *chain, struct Flags *flags, int ic, int Nseg);


/* ============================  MAIN PROGRAM  ============================ */

int main(int argc, char *argv[])
{
  
  time_t start, stop;
  start = time(NULL);
  
  int NMAX  = 12;   //max number of waveforms & time segments
  int NMCMC = 10000; //number of MCMC steps
  int NBURN = 10000; //number of burn-in steps
  
  
  /* Allocate data structures */
  struct Flags *flags   = malloc(sizeof(struct Flags));
  struct Orbit *orbit   = malloc(sizeof(struct Orbit));
  struct Chain *chain   = malloc(sizeof(struct Chain));
  struct Data  ***data  = malloc(sizeof(struct Data**)*NMAX); //data[source][segment]

  
  /* Parse command line and set defaults/flags */
  for(int i=0; i<NMAX; i++)
  {
    data[i] = malloc(sizeof(struct Data*)*NMAX);
    for(int j=0; j<NMAX; j++) data[i][j] = malloc(sizeof(struct Data));
  }
  parse(argc,argv,data,orbit,flags,chain,NMAX);
  int NC = chain->NC;

  /* Allocate model structures */
  struct Model ***trial = malloc(sizeof(struct Model**)*NC);  //trial[chain][segment]
  struct Model ****model= malloc(sizeof(struct Model***)*NC); //model[chain][source][segment]


  /* Load spacecraft ephemerides */
  switch(flags->orbit)
  {
    case 0:
      initialize_analytic_orbit(orbit);
      break;
    case 1:
      initialize_numeric_orbit(orbit);
      break;
    default:
      fprintf(stderr,"unsupported orbit type\n");
      return(1);
      break;
  }
  
  /* Initialize data structures */
  alloc_data(data, flags, NMCMC);
  
  
  /* Inject gravitational wave signal */
  if(flags->knownSource)
    GalacticBinaryInjectVerificationSource(data,orbit,flags);
  else
    GalacticBinaryInjectSimulatedSource(data,orbit,flags);

  /* Initialize data-dependent proposal */
  setup_frequency_proposal(data[0][0]);
  

  /* Initialize parallel chain */
  initialize_chain(chain, flags, &data[0][0]->cseed);
  
  
  /* Initialize data models */
  int ic;
#pragma omp parallel for private(ic) shared(model,chain,data,orbit,trial)
  for(ic=0; ic<NC; ic++)
  {
    trial[ic] = malloc(sizeof(struct Model * ) * flags->segment);
    for(int nseg = 0; nseg < flags->segment; nseg++)
    {
      trial[ic][nseg] = malloc(sizeof(struct Model));
      alloc_model(trial[ic][nseg],NMAX,data[0][nseg]->N,data[0][nseg]->Nchannel,data[0][nseg]->NP);
    }
    
    
    model[ic] = malloc(sizeof(struct Model **) * flags->injection);
    
    for(int i=0; i<flags->injection; i++)
    {
      model[ic][i] = malloc(sizeof(struct Model *) * flags->segment);
      for(int nseg = 0; nseg < flags->segment; nseg++)
      {
        model[ic][i][nseg] = malloc(sizeof(struct Model));
        
        struct Model *model_ptr = model[ic][i][nseg];
        struct Data  *data_ptr  = data[i][nseg];
        
        alloc_model(model_ptr,NMAX,data_ptr->N,data_ptr->Nchannel, data_ptr->NP);
        
        set_uniform_prior(model_ptr, data_ptr);
        
        //set noise model
        copy_noise(data_ptr->noise, model_ptr->noise);
        
        //set signal model
        for(int n=0; n<NMAX; n++)
        {
          if(flags->cheat)
          {
            struct Source *inj = data_ptr->inj;
            //map parameters to vector
            model_ptr->source[n]->NP       = inj->NP;
            model_ptr->source[n]->f0       = inj->f0;
            model_ptr->source[n]->dfdt     = inj->dfdt;
            model_ptr->source[n]->costheta = inj->costheta;
            model_ptr->source[n]->phi      = inj->phi;
            model_ptr->source[n]->amp      = inj->amp;
            model_ptr->source[n]->cosi     = inj->cosi;
            model_ptr->source[n]->phi0     = inj->phi0;
            model_ptr->source[n]->psi      = inj->psi;
            model_ptr->source[n]->d2fdt2   = inj->d2fdt2;
            map_params_to_array(model_ptr->source[n], model_ptr->source[n]->params, data_ptr->T);

          }
          else draw_from_prior(data_ptr, model_ptr, model_ptr->source[n], model_ptr->source[n]->params, chain->r[ic]);
          galactic_binary_fisher(orbit, data_ptr, model_ptr->source[n], data_ptr->noise);
          map_array_to_params(model_ptr->source[n], model_ptr->source[n]->params, data_ptr->T);
        }
        
        // Form master model & compute likelihood of starting position
        generate_noise_model(data_ptr, model_ptr);
        generate_signal_model(orbit, data_ptr, model_ptr);

        model_ptr->logL     = gaussian_log_likelihood(orbit, data_ptr, model_ptr);
        model_ptr->logLnorm = gaussian_log_likelihood_constant_norm(data_ptr, model_ptr);
        
        if(ic==0) chain->logLmax += model_ptr->logL + model_ptr->logLnorm;
        
      }//end loop over sources
    }//end loop over segments
  }//end loop over chains
    
  /* The MCMC loop */
  for(int mcmc = -NBURN; mcmc < NMCMC; mcmc++)
  {
    if(mcmc<0) flags->burnin=1;
    else       flags->burnin=0;
    
    //set annealinging tempurature during burnin
    if(flags->burnin)
    {
      chain->annealing = data[0][0]->SNR2*pow(data[0][0]->SNR2,-((double)mcmc+(double)NBURN)/((double)NBURN/(double)10));
      if(chain->annealing<1.0)chain->annealing=1.0;
      printf("T=%g; SNReff = %g\n",chain->annealing, sqrt(data[0][0]->SNR2/chain->annealing));
    }
#pragma omp parallel for private(ic) shared(model,chain,data,orbit,trial,flags)
    for(ic=0; ic<NC; ic++)
    {
      for(int i=0; i<flags->injection; i++)
      {
        struct Model **model_ptr = model[chain->index[ic]][i];
        struct Model **trial_ptr = trial[chain->index[ic]];
        struct Data  **data_ptr  = data[i];
        
        for(int steps=0; steps < 100; steps++)
        {
          galactic_binary_mcmc(orbit, data_ptr, model_ptr, trial_ptr, chain, flags, ic);
          
          /*for(int j=0; j<flags->segment; j++)
           noise_model_mcmc(orbit, data[i], model_ptr[j], trial_ptr[j], chain, flags, ic, j);*/
        }//loop over MCMC steps
        
        //delayed rejection mode-hopper
        if(mcmc<0)galactic_binary_drmc(orbit, data_ptr, model_ptr, trial_ptr, chain, flags, ic);
        
        //update fisher matrix for each chain
        if(mcmc%100==0)
        {
          for(int n=0; n<model_ptr[0]->Nlive; n++)
          {
            galactic_binary_fisher(orbit, data_ptr[0], model_ptr[0]->source[n], data_ptr[0]->noise);
          }//loop over sources in model
        }
      }//loop over frequency segments
      
      //update start time for data segments
      data_mcmc(orbit, data, model[chain->index[ic]], chain, flags, ic);
      
    }// (parallel) loop over chains
    
    ptmcmc(model,chain,flags);
    adapt_temperature_ladder(chain, mcmc+NBURN);
    
    print_chain_files(data[FIXME][FIXME], model, chain, flags, mcmc);
    
    //track maximum log Likelihood
    if(mcmc%100)update_max_log_likelihood(model, chain, flags);
    
    //store reconstructed waveform
    print_waveform_draw(data, model[chain->index[0]], flags);
    
    if(mcmc%data[FIXME][FIXME]->downsample==0)print_chain_state(data[FIXME][FIXME], chain, model[chain->index[0]][0], flags, stdout, mcmc);
    if(mcmc>0 && mcmc%data[FIXME][FIXME]->downsample==0)
    {
      for(int i=0; i<flags->injection; i++)save_waveforms(data[i][FIXME], model[chain->index[0]][i][FIXME], mcmc/data[i][FIXME]->downsample);
      for(ic=0; ic<NC; ic++)
      {
        for(int i=0; i<flags->injection; i++)
        {
          for(int j=0; j<flags->segment; j++)
          {
            chain->avgLogL[ic] += model[chain->index[ic]][i][j]->logL + model[chain->index[ic]][i][j]->logLnorm;
          }
        }
      }
    }
  }
  
  //print aggregate run files/results
  for(int i=0; i<flags->injection; i++)print_waveforms_reconstruction(data[FIXME][i],i);
  
  FILE *chainFile = fopen("avg_log_likelihood.dat","w");
  for(ic=0; ic<NC; ic++) fprintf(chainFile,"%lg %lg\n",1./chain->temperature[ic],chain->avgLogL[ic]/(double)(NMCMC/data[FIXME][FIXME]->downsample));
  fclose(chainFile);
  
  //print total run time
  stop = time(NULL);
  
  if(flags->verbose) printf(" ELAPSED TIME = %g second\n",(double)(stop-start));
  
  
  //free memory and exit cleanly
  for(ic=0; ic<NC; ic++)
  {
    free_model(model[ic][FIXME][FIXME]);
    free_model(trial[ic][FIXME]);
  }
  if(flags->orbit)free_orbit(orbit);
  //free_noise(data[0]->noise[FIXME]);
  //free_tdi(data[0]->tdi[FIXME]);
  free_chain(chain,flags);
  //free(model[FIXME][FIXME]);
  //free(trial[FIXME][FIXME]);
  //free(data[0]);
  
  return 0;
}

void ptmcmc(struct Model ****model, struct Chain *chain, struct Flags *flags)
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
    chain->acceptance[a]=0;
    
    olda = chain->index[a];
    oldb = chain->index[b];
    
    heat1 = chain->temperature[a];
    heat2 = chain->temperature[b];
    
    logL1 = 0.0;
    logL2 = 0.0;
    for(int i=0; i<flags->injection; i++)
    {
      for(int j=0; j<flags->segment; j++)
      {
        logL1 += model[olda][i][j]->logL + model[olda][i][j]->logLnorm;
        logL2 += model[oldb][i][j]->logL + model[oldb][i][j]->logLnorm;
      }
    }
    
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
  
  double nu=10;
  //double t0=100;
  double t0=10000.;
  
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

void noise_model_mcmc(struct Orbit *orbit, struct Data *data, struct Model *model, struct Model *trial, struct Chain *chain, struct Flags *flags, int ic, int Nseg)
{
  double logH  = 0.0; //(log) Hastings ratio
  double loga  = 1.0; //(log) transition probability
  
  double logPx  = 0.0; //(log) prior density for model x (current state)
  double logPy  = 0.0; //(log) prior density for model y (proposed state)
  
  //shorthand pointers
  struct Model *model_x = model;
  struct Model *model_y = trial;
  
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
    if(!flags->prior)
    {
      logH += ( (model_y->logL+model_y->logLnorm) - (model_x->logL+model_x->logLnorm) )/chain->temperature[ic]; //delta logL
      if(flags->burnin) logH /= chain->annealing;
    }
    logH += logPy  - logPx;                                         //priors
    
    loga = log(gsl_rng_uniform(chain->r[ic]));
    if(logH > loga) copy_model(model_y,model_x);
  }
  
}

void galactic_binary_mcmc(struct Orbit *orbit, struct Data **data, struct Model **model, struct Model **trial, struct Chain *chain, struct Flags *flags, int ic)
{
  double logH  = 0.0; //(log) Hastings ratio
  double loga  = 1.0; //(log) transition probability
  
  double logPx  = 0.0; //(log) prior density for model x (current state)
  double logPy  = 0.0; //(log) prior density for model y (proposed state)
  double logQyx = 0.0; //(log) proposal denstiy from x->y
  double logQxy = 0.0; //(log) proposal density from y->x
  
  //shorthand pointers
  struct Model **model_x = model;
  struct Model **model_y = trial;
  
  for(int i=0; i<flags->segment; i++)
  {
    copy_model(model_x[i],model_y[i]);
  }
  //pick a source to update
  int n = (int)(gsl_rng_uniform(chain->r[ic])*(double)model_x[0]->Nlive);
  
  //more shorthand pointers
  struct Source *source_x = model_x[0]->source[n];
  struct Source *source_y = model_y[0]->source[n];
  
  
  //choose proposal distribution
  double draw = gsl_rng_uniform(chain->r[ic]);
  if(draw < 0.1)
  /* no proposal density terms because proposal is symmetric */
    draw_from_prior(data[0], model_x[0], source_y, source_y->params, chain->r[ic]);
  else if(draw < 0.2)
  {
    draw_from_prior(data[0], model_x[0], source_y, source_y->params, chain->r[ic]);
    //TODO: Spectrum draw is not compatible with multiple sources
//    draw_from_spectrum(data[0], model_x[0], source_y, source_y->params, chain->r[ic]);
//    logQyx = log(data[0]->p[(int)(source_y->params[0]-data[0]->qmin)]);
//    logQxy = log(data[0]->p[(int)(source_x->params[0]-data[0]->qmin)]);
  }
  else if(draw<0.3)
    draw_from_extrinsic_prior(data[0], model_x[0], source_y, source_y->params, chain->r[ic]);
  else if(draw<0.8)
  /* no proposal density terms because proposal is symmetric */
    draw_from_fisher(data[0], model_x[0], source_x, source_y->params, chain->r[ic]);
  else
  {
  /* no proposal density terms because proposal is symmetric */
    fm_shift(data[0], model_x[0], source_x, source_y->params, chain->r[ic]);
  }
  
  map_array_to_params(source_y, source_y->params, data[0]->T);
  
  //hold sky position fixed to injected value
  if(flags->fixSky)
  {
    source_y->costheta = data[0]->inj->costheta;
    source_y->phi      = data[0]->inj->phi;
    map_params_to_array(source_y, source_y->params, data[0]->T);
  }
  
  //copy params for segment 0 into higher segments
  for(int i=1; i<flags->segment; i++)
  {
    copy_source(model_y[0]->source[n],model_y[i]->source[n]);
    map_params_to_array(model_y[i]->source[n], model_y[i]->source[n]->params, data[i]->T);
  }
  
  //get priors for x and y
  logPx = evaluate_uniform_prior(model_x[0], source_x->params);
  logPy = evaluate_uniform_prior(model_y[0], source_y->params);
  
  if(logPy > -INFINITY)
  {
    if(!flags->prior)
    {
      for(int i=0; i<flags->segment; i++)
      {
        //  Form master template
        generate_signal_model(orbit, data[i], model_y[i]);
        
        //get likelihood for y
        model_y[i]->logL = gaussian_log_likelihood(orbit, data[i], model_y[i]);
        
        /*
         H = [p(d|y)/p(d|x)]/T x p(y)/p(x) x q(x|y)/q(y|x)
         */
        logH += (model_y[i]->logL - model_x[i]->logL)/chain->temperature[ic]; //delta logL
        if(flags->burnin) logH /= chain->annealing;
      }
    }
    logH += logPy  - logPx;  //priors
    logH += logQxy - logQyx; //proposals
        
    loga = log(gsl_rng_uniform(chain->r[ic]));
    if(logH > loga) for(int i=0; i<flags->segment; i++) copy_model(model_y[i],model_x[i]);
    
  }
}

void galactic_binary_drmc(struct Orbit *orbit, struct Data **data, struct Model **model, struct Model **trial, struct Chain *chain, struct Flags *flags, int ic)
{
  double logH  = 0.0; //(log) Hastings ratio
  double loga  = 1.0; //(log) transition probability
  
  //shorthand pointers
  struct Data *data_ptr  = data[0];
  struct Model **model_x = model;
  struct Model **model_y = trial;
  struct Model **temp = malloc(sizeof(struct Model *) * flags->segment);
  for(int i = 0; i < flags->segment; i++)
  {
    temp[i] = malloc(sizeof(struct Model));
    alloc_model(temp[i],model_x[0]->Nmax,data_ptr->N,data_ptr->Nchannel, data_ptr->NP);
  }
  
  for(int i=0; i<flags->segment; i++)
  {
    copy_model(model_x[i],model_y[i]);
    copy_model(model_x[i],temp[i]);
  }
  
  //pick a source to update
  int n = (int)(gsl_rng_uniform(chain->r[ic])*(double)model_x[0]->Nlive);
  
  //more shorthand pointers
  struct Source *source_x = model_x[0]->source[n];
  
  // force fm-shift for temp model
  fm_shift(data[0], model_x[0], source_x, temp[0]->source[n]->params, chain->r[ic]);
  
  // do a bunch of MCMC steps to evolve from fm-shift
  for(int i=0; i<20; i++)galactic_binary_mcmc(orbit, data, temp, model_y, chain, flags, ic);
    
  // test current likelihood for model y to original likelihood of model x
  for(int i=0; i<flags->segment; i++) logH += (model_y[i]->logL - model_x[i]->logL)/chain->temperature[ic];
  if(flags->burnin) logH /= chain->annealing;
  loga = log(gsl_rng_uniform(chain->r[ic]));
  if(logH > loga) for(int i=0; i<flags->segment; i++) copy_model(model_y[i],model_x[i]);
  
  for(int i = 0; i < flags->segment; i++)
  {
    free_model(temp[i]);
  }
  free(temp);

}

void data_mcmc(struct Orbit *orbit, struct Data ***data, struct Model ***model, struct Chain *chain, struct Flags *flags, int ic)
{
  double logH  = 0.0; //(log) Hastings ratio
  double loga  = 1.0; //(log) transition probability
  double logQ  = 0.0;
  
  struct Model ***trial = malloc(sizeof(struct Model **) * flags->injection);
  
  for(int i=0; i<flags->injection; i++)
  {
    trial[i] = malloc(sizeof(struct Model *) * flags->segment);
    for(int j = 0; j < flags->segment; j++)
    {
      trial[i][j] = malloc(sizeof(struct Model));
      
      alloc_model(trial[i][j],model[i][j]->Nmax,data[i][j]->N,data[i][j]->Nchannel, data[i][j]->NP);
      
      set_uniform_prior(trial[i][j], data[i][j]);
      
      copy_model(model[i][j],trial[i][j]);
    }
  }
  
  logQ = 0.0;
  for(int i=1; i<flags->segment; i++)
  {
    
    logQ += t0_shift(data[0][i], trial[0][i], trial[0][i]->source[0], trial[0][i]->source[0]->params, chain->r[ic]);
    
    for(int j=1; j<flags->injection; j++)
    {
      trial[j][i]->t0 = trial[0][i]->t0;
    }
  }
  
  for(int j=0; j<flags->injection; j++)
  {
    for(int i=0; i<flags->segment; i++)
    {
      // Form master template
      generate_signal_model(orbit, data[j][i], trial[j][i]);
//      generate_signal_model(orbit, data[j][i], model[j][i]);
      
      // get likelihood for y
      trial[j][i]->logL = gaussian_log_likelihood(orbit, data[j][i], trial[j][i]);

      if(ic==0)
      {
        if(i==0)
        {
          if(trial[j][i]->logL!=model[j][i]->logL)
          {
            printf("%i %i %lg %lg %lg\n",i,j,trial[j][i]->t0-i*2621440, trial[j][i]->logL, model[j][i]->logL);
            printf("%lg %lg %i: ",data[j][i]->T, model[j][i]->t0, model[j][i]->source[0]->imin);
            print_source_params(data[j][i], model[j][i]->source[0], stdout);
            printf("\n");
            printf("%lg %lg %i: ",data[j][i]->T, trial[j][i]->t0, trial[j][i]->source[0]->imin);
            print_source_params(data[j][i], trial[j][i]->source[0], stdout);
            printf("\n");

          }
        }
      }
      
      /*
       H = [p(d|y)/p(d|x)]/T x p(y)/p(x) x q(x|y)/q(y|x)
       */
      if(!flags->prior)
      {
        logH += (trial[j][i]->logL - model[j][i]->logL)/chain->temperature[ic];
        if(flags->burnin) logH /= chain->annealing;
      }
      logH += logQ; //delta logL
    }
  }
  
  loga = log(gsl_rng_uniform(chain->r[ic]));
  if(logH > loga)
  {
    for(int j=0; j<flags->injection; j++)
    {
      for(int i=0; i<flags->segment; i++)
      {
        copy_model(trial[j][i],model[j][i]);
      }
    }
  }
  
  for(int i=0; i<flags->injection; i++)
  {
    for(int j = 0; j < flags->segment; j++)
    {
      free_model(trial[i][j]);
    }
    free(trial[i]);
  }
  free(trial);
  
}


