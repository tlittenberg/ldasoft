
/***************************  REQUIRED LIBRARIES  ***************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//#include "omp.h"

/*************  PROTOTYPE DECLARATIONS FOR INTERNAL FUNCTIONS  **************/

#include "LISA.h"
#include "Constants.h"
#include "BayesLine.h"
#include "GalacticBinary.h"
#include "GalacticBinaryIO.h"
#include "GalacticBinaryData.h"
#include "GalacticBinaryMath.h"
#include "GalacticBinaryPrior.h"
#include "GalacticBinaryModel.h"
#include "GalacticBinaryProposal.h"
#include "GalacticBinaryWaveform.h"



void hallowelt(struct Flags *flags);
/* ============================  MAIN PROGRAM  ============================ */

int main(int argc, char *argv[])
{
    
    if(argc!=3)
    {
        fprintf(stdout,"Usage: gb_match chain_file match_file\n");
        return 0;
    }
    
//    int NMAX = 10;   //max number of frequency & time segments
    
    /* Allocate data structures */
    struct Flags *flags = malloc(sizeof(struct Flags));
//    struct Orbit *orbit = malloc(sizeof(struct Orbit));
//    struct Chain *chain = malloc(sizeof(struct Chain));
//    struct Data  **data = malloc(sizeof(struct Data*)*NMAX); //data[NF]

    FILE *chain_file = fopen(argv[1],"r");
    FILE *match_file = fopen(argv[2],"w");
    
    double f0,dfdt,theta,phi,amp,iota,phi0,psi;//read parameters from chain file
    int n;
    //count sources in file
    int N=0;
    while(!feof(chain_file))
    {
        fscanf(chain_file,"%lg %lg %lg %lg %lg %lg %lg %lg",&f0,&dfdt,&theta,&phi,&amp,&iota,&psi,&phi0);
        N++;
    }
    rewind(chain_file);
    
    
    
    
    for(n=0; n<N; n++)
    {
        fprintf(match_file, "%lg %lg %lg %lg %lg %lg %lg %lg\n",f0,dfdt,theta,phi,amp,iota,psi,phi0);
    }
    printf("\n");
    fclose(chain_file);
    fclose(match_file);
    
    
    
    hallowelt(flags);
    
//    if(flags->orbit)free_orbit(orbit);
//    free_chain(chain,flags);

    
    return 0;
}

void hallowelt(struct Flags *flags)
{
    fprintf(stdout,"Hallo welt!!\n");
}





