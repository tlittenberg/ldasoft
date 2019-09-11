
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
#include "GalacticBinaryPrior.h"
#include "GalacticBinaryModel.h"
#include "GalacticBinaryProposal.h"
#include "GalacticBinaryWaveform.h"


void hallowelt(struct Flags *flags);
/* ============================  MAIN PROGRAM  ============================ */

int main(int argc, char *argv[])
{
    
    /* Allocate data structures */
    struct Flags *flags = malloc(sizeof(struct Flags));

    hallowelt(flags);
    return 0;
}

void hallowelt(struct Flags *flags)
{
    fprintf(stdout,"Hallo welt!!\n");
}





