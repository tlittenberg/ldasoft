//
//  GalacticBinaryProposal.c
//  
//
//  Created by Littenberg, Tyson B. (MSFC-ZP12) on 2/6/17.
//
//

#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "LISA.h"
#include "Constants.h"
#include "GalacticBinary.h"
#include "GalacticBinaryIO.h"
#include "GalacticBinaryPrior.h"
#include "GalacticBinaryModel.h"
#include "GalacticBinaryWaveform.h"
#include "GalacticBinaryProposal.h"

void setup_frequency_proposal(struct Data *data)
{
  int BW = 20;
  double *power = data->p;
  double total  = 0.0;
  FILE *temp = fopen("temp2.dat","w");
  for(int i=0; i<data->N-BW; i++)
  {
    power[i]=0.0;
    for(int n=i; n<i+BW; n++)
    {
      double SnA = data->noise->SnA[n];
      double SnE = data->noise->SnE[n];
      
      double AA = data->tdi->A[2*n]*data->tdi->A[2*n]+data->tdi->A[2*n+1]*data->tdi->A[2*n+1];
      double EE = data->tdi->E[2*n]*data->tdi->E[2*n]+data->tdi->E[2*n+1]*data->tdi->E[2*n+1];
      
      power[i] += AA/SnA + EE/SnE;
      total += power[i];
    }
  }
  for(int i=data->N-BW; i<data->N; i++)
  {
    power[i] = power[data->N-BW-1];
    total += power[i];
  }
  
  data->pmax = 0.0;
  for(int i=0; i<data->N; i++)
  {
    fprintf(temp,"%i %lg\n",i,power[i]);
    if(power[i]>data->pmax) data->pmax = power[i];
  }
  fclose(temp);
  
  
  //also get SNR^2 of data
  total = 0.0;
  for(int n=0; n<data->N; n++)
  {
      double SnA = data->noise->SnA[n];
      double SnE = data->noise->SnE[n];
      
      double AA = data->tdi->A[2*n]*data->tdi->A[2*n]+data->tdi->A[2*n+1]*data->tdi->A[2*n+1];
      double EE = data->tdi->E[2*n]*data->tdi->E[2*n]+data->tdi->E[2*n+1]*data->tdi->E[2*n+1];
      
      total += AA/SnA + EE/SnE;
  }

  data->SNR2 = total - data->N;
  printf("total=%g, N=%i\n",total, data->N);
  if(data->SNR2<0.0)data->SNR2=0.0;
  data->SNR2*=4.0;//why the factor of 4?
  printf("data-based SNR^2:  %g (%g)\n", data->SNR2, sqrt(data->SNR2));

}

void print_acceptance_rates(struct Proposal **proposal, int NP, int ic, FILE *fptr)
{
  fprintf(fptr,"Acceptance rates for chain %i:\n", ic);
  for(int n=0; n<NP+1; n++)
  {
    fprintf(fptr,"   %.1e  [%s]\n", (double)proposal[n]->accept[ic]/(double)proposal[n]->trial[ic],proposal[n]->name);
  }
}

double draw_from_spectrum(struct Data *data, struct Model *model, struct Source *source, UNUSED struct Proposal *proposal, double *params, gsl_rng *seed)
{
  //TODO: Work in amplitude
  
  //rejections ample for f
  int check = 1;
  double alpha;
  int q;
  int count=0;
  while(check)
  {
    params[0] = model->prior[0][0] + gsl_rng_uniform(seed)*(model->prior[0][1]-model->prior[0][0]);
    alpha     = gsl_rng_uniform(seed)*data->pmax;
    q = (int)(params[0]-data->qmin);
    if(alpha<data->p[q]) check = 0;
    count++;
  }
  
  //random draws for other parameters
  for(int n=1; n<source->NP; n++) params[n] = model->prior[n][0] + gsl_rng_uniform(seed)*(model->prior[n][1]-model->prior[n][0]);

  return 0;
}

double draw_from_prior(UNUSED struct Data *data, struct Model *model, UNUSED struct Source *source, UNUSED struct Proposal *proposal, double *params, gsl_rng *seed)
{
  for(int n=0; n<source->NP; n++) params[n] = model->prior[n][0] + gsl_rng_uniform(seed)*(model->prior[n][1]-model->prior[n][0]);
  
  for(int j=0; j<source->NP; j++)
  {
    if(params[j]!=params[j]) fprintf(stderr,"draw_from_prior: params[%i]=%g, U[%g,%g]\n",j,params[j],model->prior[j][0],model->prior[j][1]);
  }

  return model->logPriorVolume;
}

double draw_from_extrinsic_prior(UNUSED struct Data *data, struct Model *model, UNUSED struct Source *source, UNUSED struct Proposal *proposal, double *params, gsl_rng *seed)
{
  for(int n=1; n<3; n++) params[n] = model->prior[n][0] + gsl_rng_uniform(seed)*(model->prior[n][1]-model->prior[n][0]);
  for(int n=4; n<7; n++) params[n] = model->prior[n][0] + gsl_rng_uniform(seed)*(model->prior[n][1]-model->prior[n][0]);
  
  for(int j=0; j<source->NP; j++)
  {
    if(params[j]!=params[j]) fprintf(stderr,"draw_from_prior: params[%i]=%g, U[%g,%g]\n",j,params[j],model->prior[j][0],model->prior[j][1]);
  }
  
  return model->logPriorVolume;
}

double draw_from_fisher(UNUSED struct Data *data, struct Model *model, struct Source *source, UNUSED struct Proposal *proposal, double *params, gsl_rng *seed)
{
  int i,j;
  int NP=source->NP;
  double sqNP = sqrt((double)source->NP);
  double Amps[NP];
  double jump[NP];
  
  //draw the eigen-jump amplitudes from N[0,1] scaled by evalue & dimension
  for(i=0; i<NP; i++)
  {
    //Amps[i] = gsl_ran_gaussian(seed,1)/sqrt(source->fisher_evalue[i])/sqNP;
    Amps[i] = gsl_ran_gaussian(seed,1)/sqrt(source->fisher_evalue[i]);
    jump[i] = 0.0;
  }
  
  //decompose eigenjumps into paramter directions
  /*for(i=0; i<NP; i++) for (j=0; j<NP; j++)
  {
    jump[j] += Amps[i]*source->fisher_evectr[j][i];
    if(jump[j]!=jump[j])jump[j]=0.0;
  }*/
  
  //choose one eigenvector to jump along
  i = (int)(gsl_rng_uniform(seed)*(double)NP);
  for (j=0; j<NP; j++) jump[j] += Amps[i]*source->fisher_evectr[j][i];
  
  //jump from current position
  for(i=0; i<NP; i++) params[i] = source->params[i] + jump[i];
  
  for(int j=0; j<NP; j++)
  {
    if(params[j]!=params[j]) fprintf(stderr,"draw_from_fisher: params[%i]=%g, N[%g,%g]\n",j,params[j],source->params[j],jump[j]);
  }

  //not updating Fisher between moves, proposal is symmetric
  return 0.0;
}

double draw_from_cdf(UNUSED struct Data *data, struct Model *model, struct Source *source, struct Proposal *proposal, double *params, gsl_rng *seed)
{
  int N = proposal->size;
  int NP = source->NP;
  double **cdf = proposal->matrix;
  double logP=0.0;

  for(int n=0; n<NP; n++)
  {
    //draw p-value
    double p = gsl_rng_uniform(seed);

    //find samples either end of p-value
    int i = (int)floor(p*N);

    //linear interpolation betweein i and i+1
    //y = y0 + (x - x0)*((y1-y0)/(x1-x0))
    double x0 = (double)i/(double)N;
    double y0;
    if(i==0) y0 = model->prior[n][0];
    else     y0 = cdf[n][i];
    double x1 = (double)(i+1)/(double)N;
    double y1;
    if(i==N-1) y1 = model->prior[n][1];
    else       y1 = cdf[n][i+1];

    params[n] = y0 + (p-x0)*((y1-y0)/(x1-x0));
    logP += log(p);
  }

  return logP;
}

double cdf_density(struct Model *model, struct Source *source, struct Proposal *proposal)
{
  int N = proposal->size;
  int NP = source->NP;
  double **cdf = proposal->matrix;
  double logP=0.0;
  double *params = source->params;
  
  double x0,x1,y0,y1;
  
  int ii;
  for(int n=0; n<NP; n++)
  {
    //find samples either end of p-value
    int i = 0;
    while(params[n]>cdf[n][i]) i++;
    
    ii=i+1;
    while(cdf[n][i]==cdf[n][ii])ii++;
    
    //linear interpolation betweein i and i+1
    //y = y0 + (x - x0)*((y1-y0)/(x1-x0))
    if(i==0) x0 = model->prior[n][0];
    else     x0 = cdf[n][i];
    y0 = (double)i/(double)N;

    if(i>N-2)
    {
      x1 = model->prior[n][1];
      y1 = 1.0;
    }
    else
    {
      x1 = cdf[n][ii];
      y1 = (double)(ii)/(double)N;
    }
    
    double p = y0 + (params[n]-x0)*((y1-y0)/(x1-x0));
    
    logP += log(p);
  }
  
  return logP;
}


double fm_shift(struct Data *data, struct Model *model, struct Source *source, struct Proposal *proposal, double *params, gsl_rng *seed)
{
  //doppler modulation frequency (in bins)
  double fm = data->T/YEAR;
  
  //update all parameters with a draw from the fisher
  if(gsl_rng_uniform(seed)<0.5) draw_from_fisher(data, model, source, proposal, params, seed);
  else for(int n=1; n<source->NP; n++) params[n] = model->prior[n][0] + gsl_rng_uniform(seed)*(model->prior[n][1]-model->prior[n][0]);

  //perturb frequency by 1 fm
  double scale = floor(6*gsl_ran_gaussian(seed,1));
  
  params[0] += scale*fm;
  //params[7] += scale*fm*fm;
  
  //fm shift is symmetric
  return 0.0;
}

double t0_shift(UNUSED struct Data *data, struct Model *model, UNUSED struct Source *source, UNUSED struct Proposal *proposal, UNUSED double *params, gsl_rng *seed)
{
  //uniform draw
  if(gsl_rng_uniform(seed) < 0.5 )
    model->t0 = model->t0_min + gsl_rng_uniform(seed)*(model->t0_max - model->t0_min);
  
  //gaussian draw
  else
    model->t0 += 3.0*gsl_ran_gaussian(seed,1);
  
  
  //t0 shift is symmetric
  if(model->t0 < model->t0_min || model->t0 >= model->t0_max) return -INFINITY;
  else return 0.0;
}

void initialize_proposal(struct Data *data, struct Chain *chain, struct Flags *flags, struct Proposal **proposal, int NMAX)
{
  int NC = chain->NC;
  double check=0.0;
  for(int i=0; i<chain->NP+1; i++)
  {
    
    proposal[i]->trial  = malloc(NC*sizeof(int));
    proposal[i]->accept = malloc(NC*sizeof(int));
    
    for(int ic=0; ic<NC; ic++)
    {
      proposal[i]->trial[ic]  = 1;
      proposal[i]->accept[ic] = 0;
    }
    
    switch(i)
    {
      case 0:
        /*
         delayed rejection proposal does not fit in with others' protocal
         -must have zero weight
         -must be last in the list
         */
        sprintf(proposal[i]->name,"delayed rejection");
        proposal[i]->weight = 0.0;
        break;
        
      case 1:
        sprintf(proposal[i]->name,"prior");
        proposal[i]->function = &draw_from_prior;
        proposal[i]->weight = 0.1;
        check+=proposal[i]->weight;
        break;
      case 2:
        sprintf(proposal[i]->name,"spectrum");
        proposal[i]->function = &draw_from_spectrum;
        proposal[i]->weight = 0.0;
        check+=proposal[i]->weight;
        break;
      case 3:
        sprintf(proposal[i]->name,"extrinsic prior");
        proposal[i]->function = &draw_from_extrinsic_prior;
        proposal[i]->weight = 0.1;
        check+=proposal[i]->weight;
        break;
      case 4:
        sprintf(proposal[i]->name,"fisher");
        proposal[i]->function = &draw_from_fisher;
        proposal[i]->weight = 1.0; //that's a 1 all right.  don't panic
        break;
      case 5:
        sprintf(proposal[i]->name,"fm shift");
        proposal[i]->function = &fm_shift;
        proposal[i]->weight = 0.2;
        check+=proposal[i]->weight;
        break;
      case 6:
        sprintf(proposal[i]->name,"cdf draw");
        proposal[i]->function = &draw_from_cdf;
        proposal[i]->weight = 0.2;
        check+=proposal[i]->weight;
        //parse cdf file
        FILE *fptr = fopen(flags->cdfFile,"r");
        proposal[i]->size=0;
        double junk;
        while(!feof(fptr))
        {
          fscanf(fptr,"%lg",&junk);
          for(int j=0; j<data->NP; j++) fscanf(fptr,"%lg",&junk);
          proposal[i]->size++;
        }
        rewind(fptr);
        proposal[i]->size--;
        proposal[i]->matrix = malloc(data->NP * sizeof(double*));
        for(int j=0; j<data->NP; j++) proposal[i]->matrix[j] = malloc(proposal[i]->size * sizeof(double));
        
        struct Model *temp = malloc(sizeof(struct Model));
        alloc_model(temp,NMAX,data->N,data->Nchannel, data->NP);
        
        for(int n=0; n<proposal[i]->size; n++)
        {
          fscanf(fptr,"%lg",&junk);
          scan_source_params(data, temp->source[0], fptr);
          for(int j=0; j<data->NP; j++) proposal[i]->matrix[j][n] = temp->source[0]->params[j];
          
        }
        free_model(temp);
        free(temp);
        fclose(fptr);
        break;
        
      default:
        break;
    }
  }
  //Fisher proposal fills in the cracks
  proposal[4]->weight -= check;
  
  if(proposal[4]->weight<0.0)
  {
    fprintf(stderr,"Proposal weights not normalized (line %d of file %s)\n",__LINE__,__FILE__);
    exit(1);
  }
}
