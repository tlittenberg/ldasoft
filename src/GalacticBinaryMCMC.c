
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
#include "GalacticBinaryModel.h"
#include "GalacticBinaryWaveform.h"

/* ============================  MAIN PROGRAM  ============================ */

int main(int argc, char *argv[])
{
  
  time_t start, stop;
  start = time(NULL);
  
  int NC = 10;

  /* Allocate data structures */
  struct Flags *flags  = malloc(sizeof(struct Flags));
  struct Data  *data   = malloc(sizeof(struct Data));
  struct Orbit *orbit  = malloc(sizeof(struct Orbit));
  struct Chain **chain = malloc(sizeof(struct Chain*)*NC);
  struct Model **model = malloc(sizeof(struct Model*)*NC);
  
  
  /* Parse command line and set defaults/flags */
  parse(argc,argv,data,orbit,flags);

  
  /* Load spacecraft ephemerides */
  initialize_orbit(orbit);

  
  /* Initialize data structures */
  data->tdi = malloc(sizeof(struct TDI));
  initialize_tdi(data->tdi, data->N, data->Nchannel);
  data->noise = malloc(sizeof(struct Noise));
  initialize_noise(data->noise, data->N);

  
  /* Inject gravitational wave signal */
  if(flags->injection) GalacticBinaryInjectVerificationSource(data,orbit,flags);

  
  /* Initialize data models */
  int ic;
  for(ic=0; ic<NC; ic++)
  {
    chain[ic] = malloc(sizeof(struct Chain));
    initialize_chain(chain[ic],data->cseed, ic);
    
    model[ic] = malloc(sizeof(struct Model));
    initialize_model(model[ic],NC,data->N,data->Nchannel);
  }


  
  #pragma omp parallel for
  for(ic=0; ic<NC; ic++)
  {
    printf("Hello from process %i\n",ic);
  }

  //print total run time
  stop = time(NULL);
  
  if(flags->verbose) printf(" ELAPSED TIME = %g second\n",(double)(stop-start));
  
  return 0;
}

