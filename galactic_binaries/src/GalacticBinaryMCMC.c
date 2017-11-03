
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
#include "Constants.h"
#include "GalacticBinary.h"
#include "GalacticBinaryIO.h"
#include "GalacticBinaryData.h"
#include "GalacticBinaryPrior.h"
#include "GalacticBinaryModel.h"
#include "GalacticBinaryProposal.h"
#include "GalacticBinaryWaveform.h"


void ptmcmc(struct Model ***model, struct Chain *chain, struct Flags *flags);
void adapt_temperature_ladder(struct Chain *chain, int mcmc);

void galactic_binary_mcmc(struct Orbit *orbit, struct Data *data, struct Model *model, struct Model *trial, struct Chain *chain, struct Flags *flags, struct Prior *prior, struct Proposal **proposal, int ic);
void galactic_binary_drmc(struct Orbit *orbit, struct Data *data, struct Model *model, struct Model *trial, struct Chain *chain, struct Flags *flags, struct Prior *prior, struct Proposal **proposal, int ic);
void galactic_binary_rjmcmc(struct Orbit *orbit, struct Data *data, struct Model *model, struct Model *trial, struct Chain *chain, struct Flags *flags, struct Prior *prior, struct Proposal **proposal, int ic);

void data_mcmc(struct Orbit *orbit, struct Data **data, struct Model **model, struct Chain *chain, struct Flags *flags, struct Proposal **proposal, int ic);
void noise_model_mcmc(struct Orbit *orbit, struct Data *data, struct Model *model, struct Model *trial, struct Chain *chain, struct Flags *flags, int ic);

/* ============================  MAIN PROGRAM  ============================ */

int main(int argc, char *argv[])
{
  
  int ic;
  
  time_t start, stop;
  start = time(NULL);
  
  int NMAX = 5;   //max number of frequency & time segments
  int DMAX = 5;   //100; //max number of GB waveforms
  
  /* Allocate data structures */
  struct Flags *flags = malloc(sizeof(struct Flags));
  struct Orbit *orbit = malloc(sizeof(struct Orbit));
  struct Chain *chain = malloc(sizeof(struct Chain));
  struct Data  **data = malloc(sizeof(struct Data*)*NMAX); //data[NF]
  
  
  /* Parse command line and set defaults/flags */
  for(int i=0; i<NMAX; i++)
  {
    data[i] = malloc(sizeof(struct Data));
    data[i]->t0   = malloc( NMAX * sizeof(double) );
    data[i]->tgap = malloc( NMAX * sizeof(double) );
  }
  parse(argc,argv,data,orbit,flags,chain,NMAX);
  int NC = chain->NC;
  
  
  /* Allocate model structures */
  struct Model **trial = malloc(sizeof(struct Model*)*NC);//trial[chain]
  struct Model ***model= malloc(sizeof(struct Model**)*NC); //model[chain][source][segment]
  
  
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
  alloc_data(data, flags);

  /* Inject strain data */
  if(flags->strainData)
  {
    GalacticBinaryReadData(data,orbit,flags);
  }
  else
  {
    /* Inject gravitational wave signal */
    if(flags->knownSource)
      GalacticBinaryInjectVerificationSource(data,orbit,flags);
    else
      GalacticBinaryInjectSimulatedSource(data,orbit,flags);
  }
  
  /* Initialize data-dependent proposal */
  setup_frequency_proposal(data[0]);
  
  /* Initialize parallel chain */
  initialize_chain(chain, flags, &data[0]->cseed);
  
  /* Initialize MCMC proposals */
  printf("chain->NP=%i\n",chain->NP);
  struct Proposal ***proposal = malloc(NMAX*sizeof(struct Proposal**));
  for(int j=0; j<NMAX; j++)
  {
    proposal[j] = malloc((chain->NP+1)*sizeof(struct Proposal*));
    for(int i=0; i<chain->NP+1; i++) proposal[j][i] = malloc(sizeof(struct Proposal));
  
  }
  for(int j=0; j<flags->NF; j++) initialize_proposal(orbit, data[j], chain, flags, proposal[j], DMAX);

  /* Initialize priors */
  struct Prior *prior = malloc(sizeof(struct Prior));
  if(flags->skyPrior) set_galaxy_prior(flags, prior);

  
  /* Initialize data models */
  for(ic=0; ic<NC; ic++)
  {
    //printf("initialize model\n");

    trial[ic] = malloc(sizeof(struct Model));
    alloc_model(trial[ic],DMAX,data[0]->N,data[0]->Nchannel,data[0]->NP, data[0]->NT);
    
    model[ic] = malloc(sizeof(struct Model *) * flags->NF);
    
    //loop over frequency segments
    for(int i=0; i<flags->NF; i++)
    {
      //printf("frequency segment %i\n",i);

      model[ic][i] = malloc(sizeof(struct Model));
      
      struct Model *model_ptr = model[ic][i];
      struct Data  *data_ptr  = data[i];
      
      alloc_model(model_ptr,DMAX,data_ptr->N,data_ptr->Nchannel, data_ptr->NP, flags->NT);
      
      if(ic==0)set_uniform_prior(flags, model_ptr, data_ptr, 1);
      else     set_uniform_prior(flags, model_ptr, data_ptr, 0);
      
      //set noise model
      for(int j=0; j<flags->NT; j++) copy_noise(data_ptr->noise[j], model_ptr->noise[j]);
      
      //set signal model
      for(int n=0; n<DMAX; n++)
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
        else
        {
          draw_from_prior(data_ptr, model_ptr, model_ptr->source[n], proposal[i][0], model_ptr->source[n]->params , chain->r[ic]);
        }
        map_array_to_params(model_ptr->source[n], model_ptr->source[n]->params, data_ptr->T);
        galactic_binary_fisher(orbit, data_ptr, model_ptr->source[n], data_ptr->noise[0]);
      }
      
      // Form master model & compute likelihood of starting position
      generate_noise_model(data_ptr, model_ptr);
      generate_signal_model(orbit, data_ptr, model_ptr);
      
      if(!flags->prior)
      {
        model_ptr->logL     = gaussian_log_likelihood(orbit, data_ptr, model_ptr);
        model_ptr->logLnorm = gaussian_log_likelihood_constant_norm(data_ptr, model_ptr);
      }
      else model_ptr->logL = model_ptr->logLnorm = 0.0;
      
      if(ic==0) chain->logLmax += model_ptr->logL + model_ptr->logLnorm;
      
    }//end loop over frequency segments
  }//end loop over chains
  
  
  /* The MCMC loop */
  for(int mcmc = -flags->NBURN; mcmc < flags->NMCMC; mcmc++)
  {
    if(mcmc<0) flags->burnin=1;
    else       flags->burnin=0;
    
    //set annealinging tempurature during burnin
//    if(flags->burnin)
//    {
//      chain->annealing = data[0]->SNR2*pow(data[0]->SNR2,-((double)mcmc+(double)flags->NBURN)/((double)flags->NBURN/(double)10))/40.;
//      if(chain->annealing<1.0)chain->annealing=1.0;
//      chain->annealing=1.0;
//      //printf("annealing=%g\n",chain->annealing);
//    }
    chain->annealing=1.0;
    
    // (parallel) loop over chains
    //#pragma omp parallel for private(ic) shared(flags,model,trial,chain,orbit,proposal)
    for(ic=0; ic<NC; ic++)
    {
      
      //loop over frequency segments
      for(int i=0; i<flags->NF; i++)
      {
        struct Model *model_ptr = model[chain->index[ic]][i];
        struct Model *trial_ptr = trial[chain->index[ic]];
        struct Data  *data_ptr  = data[i];
        
        
        for(int steps=0; steps < 100; steps++)
        {
          //for(int j=0; j<model_ptr->Nlive; j++)
            galactic_binary_mcmc(orbit, data_ptr, model_ptr, trial_ptr, chain, flags, prior, proposal[i], ic);

          noise_model_mcmc(orbit, data_ptr, model_ptr, trial_ptr, chain, flags, ic);
        }//loop over MCMC steps
        
        
        //reverse jump birth/death move
        if(flags->rj)galactic_binary_rjmcmc(orbit, data_ptr, model_ptr, trial_ptr, chain, flags, prior, proposal[i], ic);

        //delayed rejection mode-hopper
        //if(model_ptr->Nlive>0 && mcmc<0 && ic<NC/2)galactic_binary_drmc(orbit, data_ptr, model_ptr, trial_ptr, chain, flags, prior, proposal[i], ic);

        //update fisher matrix for each chain
        if(mcmc%100==0)
        {
          for(int n=0; n<model_ptr->Nlive; n++)
          {
            galactic_binary_fisher(orbit, data_ptr, model_ptr->source[n], data_ptr->noise[FIXME]);
          }
        }
      }//end loop over frequency segments
      
      //update start time for data segments
      //data_mcmc(orbit, data, model[chain->index[ic]], chain, flags, proposal[i], ic);
      
    }// end (parallel) loop over chains
    
    ptmcmc(model,chain,flags);
    adapt_temperature_ladder(chain, mcmc+flags->NBURN);
    
    print_chain_files(data[FIXME], model, chain, flags, mcmc);
    
    //track maximum log Likelihood
    if(mcmc%100)
    {
      if(update_max_log_likelihood(model, chain, flags)) mcmc = -flags->NBURN;
    }
    
    //store reconstructed waveform
    print_waveform_draw(data, model[chain->index[0]], flags);
    
    //update run status
    if(mcmc%data[FIXME]->downsample==0)
    {
      for(int i=0; i<flags->NF; i++)
      {
        print_chain_state(data[i], chain, model[chain->index[0]][i], flags, stdout, mcmc);
        fprintf(stdout,"Sources: %i\n",model[chain->index[0]][i]->Nlive);
        print_acceptance_rates(proposal[i], chain->NP, 0, stdout);
      }
    }
    
    //dump waveforms to file, update avgLogL for thermodynamic integration
    if(mcmc>0 && mcmc%data[FIXME]->downsample==0)
    {
      for(int i=0; i<flags->NF; i++)save_waveforms(data[i], model[chain->index[0]][i], mcmc/data[i]->downsample);
      for(ic=0; ic<NC; ic++)
      {
        chain->dimension[ic][model[chain->index[ic]][0]->Nlive]++;
        for(int i=0; i<flags->NF; i++)
          chain->avgLogL[ic] += model[chain->index[ic]][i]->logL + model[chain->index[ic]][i]->logLnorm;
      }
    }
    
  }// end MCMC loop
  
  //print aggregate run files/results
  for(int i=0; i<flags->NF; i++)print_waveforms_reconstruction(data[i],i);
  
  FILE *chainFile = fopen("avg_log_likelihood.dat","w");
  for(ic=0; ic<NC; ic++) fprintf(chainFile,"%lg %lg\n",1./chain->temperature[ic],chain->avgLogL[ic]/(double)(flags->NMCMC/data[FIXME]->downsample));
  fclose(chainFile);
  
  FILE *zFile = fopen("evidence.dat","w");
  for(int i=0; i<DMAX; i++) fprintf(zFile,"%i %i\n",i,chain->dimension[0][i]);
  fclose(zFile);
  
  //print total run time
  stop = time(NULL);
  
  if(flags->verbose) printf(" ELAPSED TIME = %g second\n",(double)(stop-start));
  
  
  //free memory and exit cleanly
  for(ic=0; ic<NC; ic++)
  {
    for(int i=0; i<flags->NF; i++) free_model(model[ic][i]);
    free_model(trial[ic]);
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

void ptmcmc(struct Model ***model, struct Chain *chain, struct Flags *flags)
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
    for(int i=0; i<flags->NF; i++)
    {
      logL1 += model[olda][i]->logL + model[olda][i]->logLnorm;
      logL2 += model[oldb][i]->logL + model[oldb][i]->logLnorm;
    }
    
    //Hot chains jump more rarely
    if(gsl_rng_uniform(chain->r[a])<1.0)
    {
      dlogL = logL2 - logL1;
      H  = (heat2 - heat1)/(heat2*heat1);
      
      alpha = exp(dlogL*H);
      beta  = gsl_rng_uniform(chain->r[a]);
      
      if(alpha >= beta)
      {
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
    
    if(chain->temperature[ic]/chain->temperature[ic-1] < 1.01) chain->temperature[ic] = chain->temperature[ic-1]*1.01;
  }//end loop over ic
}//end adapt function

void noise_model_mcmc(struct Orbit *orbit, struct Data *data, struct Model *model, struct Model *trial, struct Chain *chain, struct Flags *flags, int ic)
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
  for(int i=0; i<flags->NT; i++)
  {
    switch(data->Nchannel)
    {
      case 1:
        model_y->noise[i]->etaX = model_x->noise[i]->etaX + 0.1*gsl_ran_gaussian(chain->r[ic],1);
        break;
      case 2:
        model_y->noise[i]->etaA = model_x->noise[i]->etaA + 0.1*gsl_ran_gaussian(chain->r[ic],1);
        model_y->noise[i]->etaE = model_x->noise[i]->etaE + 0.1*gsl_ran_gaussian(chain->r[ic],1);
        break;
    }
    
    //get priors for x and y
    switch(data->Nchannel)
    {
      case 1:
        if(model_y->noise[i]->etaX < 0.1 || model_y->noise[i]->etaX>10) logPy=-INFINITY;
        break;
      case 2:
        if(model_y->noise[i]->etaA < 0.1 || model_y->noise[i]->etaA>10.) logPy=-INFINITY;
        if(model_y->noise[i]->etaE < 0.1 || model_y->noise[i]->etaE>10.) logPy=-INFINITY;
        break;
    }
  }
  
  if(logPy > -INFINITY)
  {
    
    if(!flags->prior)
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
      if(flags->burnin) logH /= chain->annealing;
    }
    logH += logPy  - logPx;                                         //priors
    
    loga = log(gsl_rng_uniform(chain->r[ic]));
    if(logH > loga) copy_model(model_y,model_x);
  }
  
}

void galactic_binary_mcmc(struct Orbit *orbit, struct Data *data, struct Model *model, struct Model *trial, struct Chain *chain, struct Flags *flags, struct Prior *prior, struct Proposal **proposal, int ic)
{
  double logH  = 0.0; //(log) Hastings ratio
  double loga  = 1.0; //(log) transition probability
  
  double logPx  = 0.0; //(log) prior density for model x (current state)
  double logPy  = 0.0; //(log) prior density for model y (proposed state)
  double logQyx = 0.0; //(log) proposal denstiy from x->y
  double logQxy = 0.0; //(log) proposal density from y->x
  
  //shorthand pointers
  struct Model *model_x = model;
  struct Model *model_y = trial;
  
  copy_model(model_x,model_y);
  
  //pick a source to update
  int n = (int)(gsl_rng_uniform(chain->r[ic])*(double)model_x->Nlive);
  
  //more shorthand pointers
  struct Source *source_x = model_x->source[n];
  struct Source *source_y = model_y->source[n];
  
  
  //choose proposal distribution
  int trial_n;
  double trial_w;
  int nprop=0;
  
  while(nprop<1)
  {
    trial_n = (int)floor((chain->NP+1)*gsl_rng_uniform(chain->r[ic]));
    trial_w = gsl_rng_uniform(chain->r[ic]);
    if(trial_w < proposal[trial_n]->weight) nprop = trial_n;
  }
  proposal[nprop]->trial[ic]++;
  
  (*proposal[nprop]->function)(data, model_x, source_y, proposal[nprop], source_y->params, chain->r[ic]);
  
  //TODO: Fix this
  if(!strcmp(proposal[nprop]->name,"fstat"))
  {
    logQyx = evaluate_fstatistic_proposal(data, proposal[nprop], source_y->params);
    logQxy = evaluate_fstatistic_proposal(data, proposal[nprop], source_x->params);
  }
  
  if(!strcmp(proposal[nprop]->name,"cdf draw"))
  {
    logQyx = cdf_density(model_x, source_y, proposal[nprop]);
    logQxy = cdf_density(model_x, source_x, proposal[nprop]);
  }
  
  map_array_to_params(source_y, source_y->params, data->T);
  
  //hold sky position fixed to injected value
  if(flags->fixSky)
  {
    source_y->costheta = data->inj->costheta;
    source_y->phi      = data->inj->phi;
    map_params_to_array(source_y, source_y->params, data->T);
  }
  
  //copy params for segment 0 into higher segments
  copy_source(model_y->source[n],model_y->source[n]);
  map_params_to_array(model_y->source[n], model_y->source[n]->params, data->T);
  
  //get priors for x and y
  logPx = evaluate_prior(flags, model_x, prior, source_x->params);
  logPy = evaluate_prior(flags, model_y, prior, source_y->params);
  
  if(logPy > -INFINITY)
  {
    if(!flags->prior)
    {
      //  Form master template
      generate_signal_model(orbit, data, model_y);
      
      //get likelihood for y
      model_y->logL = gaussian_log_likelihood(orbit, data, model_y);
      
      /*
       H = [p(d|y)/p(d|x)]/T x p(y)/p(x) x q(x|y)/q(y|x)
       */
      logH += (model_y->logL - model_x->logL)/chain->temperature[ic]; //delta logL
      if(flags->burnin) logH /= chain->annealing;
      //        if(!strcmp(proposal[nprop]->name,"cdf draw") && ic==0)
      //        {
      //          printf("cdf logH=%g, logLx=%g, logLy=%g\n",logH,model_x[i]->logL , model_y[i]->logL);
      //          printf("   dlogQ=%g, logQxy=%g, logQyx=%g\n",logQxy - logQyx,logQxy,logQyx);
      //        }
    }
    logH += logPy  - logPx;  //priors
    logH += logQxy - logQyx; //proposals
    
    loga = log(gsl_rng_uniform(chain->r[ic]));
    if(logH > loga)
    {
      proposal[nprop]->accept[ic]++;
      copy_model(model_y,model_x);
    }
  }
}

void galactic_binary_rjmcmc(struct Orbit *orbit, struct Data *data, struct Model *model, struct Model *trial, struct Chain *chain, struct Flags *flags, struct Prior *prior, struct Proposal **proposal, int ic)
{
  double logH  = 0.0; //(log) Hastings ratio
  double loga  = 1.0; //(log) transition probability
  
  double logPx  = 0.0; //(log) prior density for model x (current state)
  double logPy  = 0.0; //(log) prior density for model y (proposed state)
  double logQyx = 0.0; //(log) proposal denstiy from x->y
  double logQxy = 0.0; //(log) proposal density from y->x
  
  //shorthand pointers
  struct Model *model_x = model;
  struct Model *model_y = trial;
  
  copy_model(model_x,model_y);
  
  int freqflag=0;
  if(gsl_rng_uniform(chain->r[ic])<1.5) freqflag=1;
  
  //proposal[2]->trial[ic]++;

  /* pick birth or death move */
  if(gsl_rng_uniform(chain->r[ic])<0.5)/* birth move */
  {
    //ny=nx+1
    model_y->Nlive++;
    
    //slot new source in at end of live  source array
    int create = model_y->Nlive-1;
    
    if(model_y->Nlive<model_x->Nmax)
    {
      //draw new parameters
      if(freqflag) draw_from_fstatistic(data, model_y, model_y->source[create], proposal[2], model_y->source[create]->params, chain->r[ic]);
      else         draw_from_prior(data, model_y, model_y->source[create], proposal[0], model_y->source[create]->params, chain->r[ic]);

      
      map_array_to_params(model_y->source[create], model_y->source[create]->params, data->T);

      logQxy = 0;
      logQyx = model_x->logPriorVolume;
      if(freqflag)
      {
//        logQyx += log(model->prior[0][1]-model->prior[1][0]);                         //undo uniform frequency contribution to volume
//        logQyx += log(data->p[(int)(model_y->source[create]->params[0]-data->qmin)]); //add power-based frequency proposal
        logQyx += evaluate_fstatistic_proposal(data, proposal[2], model_y->source[create]->params);
      }
      
      //copy params for segment 0 into higher segments
//      copy_source(model_y->source[create],model_y->source[create]);
//      map_params_to_array(model_y->source[create], model_y->source[create]->params, data->T);
      
    }
    else logPy = -INFINITY;
  }
  else /* death move */
  {
    //ny=nx-1
    model_y->Nlive--;
    
    //pick source to kill
    int kill = (int)(gsl_rng_uniform(chain->r[ic])*(double)model_x->Nlive);
    
    if(model_y->Nlive>-1)
    {
      //copy params for segment 0 into higher segments
//      copy_source(model_y->source[kill],model_y->source[kill]);
//      map_params_to_array(model_y->source[kill], model_y->source[kill]->params, data->T);
      
      logQxy = model_x->logPriorVolume;
      logQyx = 0;
      if(freqflag)
      {
        //logQxy += log(model->prior[0][1]-model->prior[1][0]);
        //logQxy += log(data->p[(int)(model_y->source[kill]->params[0]-data->qmin)]);
        logQxy += evaluate_fstatistic_proposal(data, proposal[2], model_y->source[kill]->params);

      }
      
      //consolodiate parameter structure
      for(int j=kill; j<model_x->Nlive; j++)
      {
        copy_source(model_x->source[j+1],model_y->source[j]);
      }
    }
    else logPy = -INFINITY;
  }
  
  for(int n=0; n<model_x->Nlive; n++) logPx +=  evaluate_prior(flags, model_x, prior, model_x->source[n]->params);
  for(int n=0; n<model_y->Nlive; n++) logPy +=  evaluate_prior(flags, model_y, prior, model_y->source[n]->params);
  
  
  /* Hasting's ratio */
  if(logPy > -INFINITY && !flags->prior)
  {
    //  Form master template
    generate_signal_model(orbit, data, model_y);
    
    //get likelihood for y
    model_y->logL = gaussian_log_likelihood(orbit, data, model_y);
    
    /*
     H = [p(d|y)/p(d|x)]/T x p(y)/p(x) x q(x|y)/q(y|x)
     */
    logH += (model_y->logL - model_x->logL)/chain->temperature[ic]; //delta logL
    if(flags->burnin) logH /= chain->annealing;
  }
  
  
  logH += logPy  - logPx;  //priors
  logH += logQxy - logQyx; //proposals
  
//  if(model_y->Nlive > model_x->Nlive && ic==0)
//    if(ic==0)
//    {
//      FILE *fptr = fopen("proposal.dat","a");
//      //    fprintf(stdout,"%lg %lg %lg %lg %g ",model_y->logL+model_y->logLnorm, model_x->logL+model_x->logLnorm, logH, logPy  - logPx, logQxy - logQyx);
//      //    fprintf(stdout,"%i -> %i ",model_x->Nlive, model_y->Nlive);
//      print_source_params(data, model_y->source[model_y->Nlive-1], fptr);
//      fprintf(fptr,"\n");
//      fclose(fptr);
//    }
  
  loga = log(gsl_rng_uniform(chain->r[ic]));
  if(logH > loga)
  {
    //proposal[2]->accept[ic]++;
    copy_model(model_y,model_x);
  }
  
}

void galactic_binary_drmc(struct Orbit *orbit, struct Data *data, struct Model *model, struct Model *trial, struct Chain *chain, struct Flags *flags, struct Prior *prior, struct Proposal **proposal, int ic)
{
  double logH  = 0.0; //(log) Hastings ratio
  double loga  = 1.0; //(log) transition probability
  
  //shorthand pointers
  struct Data *data_ptr  = data;
  struct Model *model_x = model;
  struct Model *model_y = trial;
  struct Model *temp = malloc(sizeof(struct Model) );
  alloc_model(temp,model_x->Nmax,data_ptr->N,data_ptr->Nchannel, data_ptr->NP, flags->NT);
  
  //keep track of DR trials for acceptance rates
  proposal[0]->trial[ic]++;
  
  /* Initialize custom MCMC proposal for DR move */
  struct Proposal **proposal_temp = malloc((chain->NP+1)*sizeof(struct Proposal*));
  double check=0.0;
  for(int i=0; i<=chain->NP; i++)
  {
    proposal_temp[i] = malloc(sizeof(struct Proposal));
    proposal_temp[i]->trial  = malloc(chain->NC*sizeof(int));
    proposal_temp[i]->accept = malloc(chain->NC*sizeof(int));
    
    for(int j=0; j<chain->NC; j++)
    {
      proposal_temp[i]->trial[j]  = 1;
      proposal_temp[i]->accept[j] = 0;
    }
    
    switch(i)
    {
      case 1:
        sprintf(proposal_temp[i]->name,"fisher");
        proposal_temp[i]->function = &draw_from_fisher;
        proposal_temp[i]->weight = 0.1;
        break;
      case 2:
        sprintf(proposal_temp[i]->name,"fisher");
        proposal_temp[i]->function = &draw_from_fisher;
        proposal_temp[i]->weight = 0.1;
        break;
      case 3:
        sprintf(proposal_temp[i]->name,"fisher");
        proposal_temp[i]->function = &draw_from_fisher;
        proposal_temp[i]->weight = 0.1;
        break;
      case 4:
        sprintf(proposal_temp[i]->name,"fisher");
        proposal_temp[i]->function = &draw_from_fisher;
        proposal_temp[i]->weight = 0.5;
        break;
      case 5:
        sprintf(proposal_temp[i]->name,"fisher");
        proposal_temp[i]->function = &draw_from_fisher;
        proposal_temp[i]->weight = 0.2;
        break;
      default:
        sprintf(proposal_temp[0]->name,"fisher");
        proposal_temp[0]->function = &draw_from_fisher;
        proposal_temp[0]->weight = 0.0;
        break;
    }
    check+=proposal_temp[i]->weight;
  }
  if(check<0.99999)
  {
    fprintf(stderr,"Proposal weights not normalized (line %d of file %s)\n",__LINE__,__FILE__);
    exit(1);
  }
  
  copy_model(model_x,model_y);
  copy_model(model_x,temp);
  
  //pick a source to update
  int n = (int)(gsl_rng_uniform(chain->r[ic])*(double)model_x->Nlive);
  
  //more shorthand pointers
  struct Source *source_x = model_x->source[n];
  
  
  //  // force fm-shift for temp model
  //    fm_shift(data[0], model_x[0], source_x, temp[0]->source[n]->params, chain->r[ic]);
  //    double jump = temp[0]->source[n]->params[0] - model_x[0]->source[n]->params[0];
  //    temp[0]->source[n]->params[0] = model_x[0]->source[n]->params[0] + jump*sqrt(chain->annealing);
  //    temp[0]->source[n]->params[7] = model_x[0]->source[n]->params[7] - 2*jump*sqrt(chain->annealing);
  //  }
  //  else if(gsl_rng_uniform(chain->r[ic])<0.7)
  //  {
  //    // force extrinsic update
  //    draw_from_extrinsic_prior(data[0], model_x[0], source_x, temp[0]->source[n]->params, chain->r[ic]);
  //  }
  
  if(gsl_rng_uniform(chain->r[ic])<0.5)
    draw_from_fstatistic(data, model_x, source_x, proposal[2], temp->source[n]->params, chain->r[ic]);
  else
    fm_shift(data, model_x, source_x, proposal[4], temp->source[n]->params, chain->r[ic]);
  
  
  generate_signal_model(orbit, data, temp);
  temp->logL = gaussian_log_likelihood(orbit, data, temp);

  //if(ic==0)printf("delayed rejection hastings ratio: %g: ",temp->logL);
  
  
  // do a bunch of MCMC steps to evolve from forced jump
  for(int i=0; i<100; i++)galactic_binary_mcmc(orbit, data, temp, model_y, chain, flags, prior, proposal_temp, ic);
  
  // test current likelihood for model y to original likelihood of model x
  logH += (model_y->logL - model_x->logL)/chain->temperature[ic];
  //if(ic==0)printf("%g = %g - %g\n",logH,model_y[0]->logL,model_x[0]->logL);
  if(flags->burnin) logH /= chain->annealing;
  loga = log(gsl_rng_uniform(chain->r[ic]));
  if(logH > loga)
  {
    proposal[chain->NP]->accept[ic]++;
    copy_model(model_y,model_x);
  }
  
  free_model(temp);
  
  for(int i=0; i<=chain->NP; i++)
  {
    free(proposal_temp[i]->trial);
    free(proposal_temp[i]->accept);
    free(proposal_temp[i]);
  }
  free(proposal_temp);
  
}

void data_mcmc(struct Orbit *orbit, struct Data **data, struct Model **model, struct Chain *chain, struct Flags *flags, struct Proposal **proposal, int ic)
{
  double logH  = 0.0; //(log) Hastings ratio
  double loga  = 1.0; //(log) transition probability
  double logQ  = 0.0;
  
  struct Model **trial = malloc(sizeof(struct Model *) * flags->NF);
  
  for(int i=0; i<flags->NF; i++)
  {
    trial[i] = malloc(sizeof(struct Model));
    
    
    alloc_model(trial[i],model[i]->Nmax,data[i]->N,data[i]->Nchannel, data[i]->NP, flags->NT);
    
    set_uniform_prior(flags, trial[i], data[i], 0);
    
    copy_model(model[i],trial[i]);
  }
  
  logQ = 0.0;
  
  logQ += t0_shift(data[0], trial[0], trial[0]->source[0], proposal[0], trial[0]->source[0]->params, chain->r[ic]);
  
  for(int j=1; j<flags->NF; j++)
  {
    for(int i=0; i<flags->NT; i++)
    {
      trial[j]->t0[i] = trial[0]->t0[i];
    }
  }
  
  for(int j=0; j<flags->NF; j++)
  {
    // Form master template
    generate_signal_model(orbit, data[j], trial[j]);
    //      generate_signal_model(orbit, data[j], model[j]);
    
    // get likelihood for y
    trial[j]->logL = gaussian_log_likelihood(orbit, data[j], trial[j]);
    
    /*
     H = [p(d|y)/p(d|x)]/T x p(y)/p(x) x q(x|y)/q(y|x)
     */
    if(!flags->prior)
    {
      logH += (trial[j]->logL - model[j]->logL)/chain->temperature[ic];
      if(flags->burnin) logH /= chain->annealing;
    }
    logH += logQ; //delta logL
  }
  
  loga = log(gsl_rng_uniform(chain->r[ic]));
  
  if(logH > loga)
  {
    for(int j=0; j<flags->NF; j++)
    {
      copy_model(trial[j],model[j]);
    }
  }
  
  for(int i=0; i<flags->NF; i++)
  {
    free_model(trial[i]);
  }
  free(trial);
  
}




