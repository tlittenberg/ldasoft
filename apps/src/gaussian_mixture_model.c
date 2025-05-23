/*
 * Copyright 2020 Tyson B. Littenberg
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 *
 * Usage: gaussian_mixture_model -h
 * Example: Fit 4 gaussians to the first two columns of `chain.dat`
 *          gaussian_mixture_model --file chain.dat --modes 4 --nparams 2
 *
 */

/**
@file gaussian_mixture_model.c
\brief Gaussian Mixture Model fit to input data samples using Expectation-Maximization Algorithm.
 
 \todo More generalization is possible.
*/


/*
 REQUIRED LIBRARIES
 */

#include <glass_utils.h>

/* map parameters from R to [xmin,xmax] with sigmoid function
for(size_t n=0; n<NP; n++)
{
    logit(params[n], y_vec, pxminmax);

    for(size_t i=0; i<NMCMC; i++)
    {
        y = y_vec[i];
        samples[i]->x[n] = y;
    }
}*/

/**
 * \brief Main function for data handling and calling GM_with_EM() algorithm
 *
 */
int main(int argc, char* argv[])
{
    // chain file
    FILE *chainFile=NULL;
    
    // dimension of model
    size_t NP = 0;
    
    // maximum number of modes
    size_t NMODE = 2;
    
    // number of samples
    size_t NMCMC = 0;
        
    // number of EM iterations
    size_t NSTEP = 500;

    // thinning rate of input chain
    size_t NTHIN = 1;
    
    // array for flagging parameters to log
    size_t *LFLAG = calloc(100,sizeof(size_t));
    
    // rng
    unsigned int r = 150914;
    
    
    /* parse command line */
    struct option cmd_opt[] =
    {
        { "help",     no_arg,  0, 'h' },
        { "file",     req_arg, 0, 'f' },
        { "log",      req_arg, 0, 'l' },
        { "modes",    req_arg, 0, 'm' },
        { "nparams",  req_arg, 0, 'n' },
        { "seed",     req_arg, 0, 's' },
        { "thin",     req_arg, 0, 't' },
        { 0, 0, 0, 0 }
    };
    
    char args[] = "hl:f:m:n:s:t:";
    char *program = argv[0];
    
    if(argc==1)
    {
        printUsage(program);
        exit(0);
    }
    while(1)
    {
        int opt_indx = 0;
        int c;
        
        c = getopt_long(argc, argv, args, cmd_opt, &opt_indx );
        if(c==-1) // end of options
            break;
        
        switch(c)
        {
            case 'f':
                chainFile = fopen(optarg,"r");
                break;
            case 'h': // help
                printUsage(program);
                exit(0);
            case 'l':
                LFLAG[(size_t)atoi(optarg)] = 1;
                break;
            case 'm':
                NMODE = (size_t)atoi(optarg);
                break;
            case 'n':
                NP = (size_t)atoi(optarg);
                break;
            case 's':
                r = (unsigned int)atoi(optarg);
                break;
            case 't':
                NTHIN = (size_t)atoi(optarg);
                break;
            default:
                fprintf(stderr,"unknown error while parsing options\n" );
                exit(1);
        }
    }
    
    /* count lines in file */
    char* line;
    char lineBuffer[BUFFER_SIZE];
    while((line = fgets(lineBuffer, BUFFER_SIZE, chainFile)) != NULL) NMCMC++;
    rewind(chainFile);
    
    //thin chain
    NMCMC /= NTHIN;
    
    
    /* allocate memory in data structures*/
    
    // raw chain samples
    double **params = double_matrix(NP,NMCMC);
    
    // mapped chain samples
    struct Sample **samples = malloc(NMCMC*sizeof(struct Sample*));
    for(size_t n=0; n<NMCMC; n++)
    {
        samples[n] = malloc(sizeof(struct Sample));
        samples[n]->x = double_vector(NP);
        samples[n]->p = double_vector(NMODE);
        samples[n]->w = double_vector(NMODE);
    }
    
    // covariance matrices for different modes
    struct MVG **modes = malloc(NMODE*sizeof(struct MVG*));
    for(size_t n=0; n<NMODE; n++)
    {
        modes[n] = malloc(sizeof(struct MVG));
        alloc_MVG(modes[n],NP);
    }

    // Logistic mapping of samples onto R
    double y;
    double pmin,pmax;
    double *y_vec = double_vector(NMCMC);

    /* parse chain file */
    double value;
    char *column;
    for(size_t i=0; i<NMCMC; i++)
    {
        for(size_t j=0; j<NTHIN; j++) line = fgets(lineBuffer, BUFFER_SIZE, chainFile);
        
        column=strtok(line," ");
        
        for(size_t n=0; n<NP; n++)
        {
            sscanf(column, "%lg", &value);
            if(LFLAG[n]) value = log(value);
            params[n][i] = value;
            column=strtok(NULL," ");
        }
    }
    
    
    /* Get max and min for each parameter */
    int *sorted_indicies = int_vector(NMCMC);
    for(size_t n=0; n<NP; n++)
    {
        index_sort(sorted_indicies,params[n],NMCMC);
        pmin = params[n][sorted_indicies[0]] - (params[n][sorted_indicies[9]]-params[n][sorted_indicies[0]]); //pad min to avoid infs in mapping
        pmax = params[n][sorted_indicies[NMCMC-1]] + (params[n][sorted_indicies[NMCMC-1]]-params[n][sorted_indicies[NMCMC-10]]); //pad min to avoid infs in mapping
        
        /* cpopy max and min into each MVG structure */
        for(size_t k=0; k<NMODE; k++)
        {
            modes[k]->minmax[n][0] = pmin;
            modes[k]->minmax[n][1] = pmax;
        }
    }
    free_int_vector(sorted_indicies);
    
    
    /* map params to R with logit function */
    for(size_t n=0; n<NP; n++)
    {
        double pmin = modes[0]->minmax[n][0];
        double pmax = modes[0]->minmax[n][1];
        logit_mapping(params[n], y_vec, pmin, pmax, NMCMC);

        for(size_t i=0; i<NMCMC; i++)
        {
            y = y_vec[i];
            samples[i]->x[n] = y;
        }
    }
    
    /* scatter plot of remapped parameters */
    FILE *scatter = fopen("remapped_samples.dat","w");
    for(size_t i=0; i<NMCMC; i++)
    {
        for(size_t n=0; n<NP; n++)
        {
            fprintf(scatter, "%lg ", samples[i]->x[n]);
        }
        fprintf(scatter,"\n");
    }
    fclose(scatter);

    /* The main Gaussian Mixture Model with Expectation Maximization function */
    double logL, BIC;
    if(GMM_with_EM(modes,samples, NP, NMODE, NMCMC,NSTEP,&r,&logL,&BIC)) return 1;
    

    /* Write GMM results to binary for pick up by other processes */
    char filename[BUFFER_SIZE];
    sprintf(filename,"gmm_%i.bin",(int)NMODE);
    FILE *fptr = fopen(filename,"wb");
    for(size_t n=0; n<NMODE; n++) write_MVG(modes[n],fptr);
    fclose(fptr);

    /* Read GMM results to binary for pick up by other processes */
    sprintf(filename,"gmm_%i.bin",(int)NMODE);
    fptr = fopen(filename,"rb");
    for(size_t n=0; n<NMODE; n++) read_MVG(modes[n],fptr);
    fclose(fptr);

    
    print_model(modes, samples, NP, NMODE, NMCMC, logL, BIC, NMODE);

    
    for(size_t n=0; n<NMODE; n++) fprintf(stdout,"%i %g\n",(int)n, modes[n]->p);

    
    /* clean up */
    for(size_t n=0; n<NMCMC; n++)
    {
        free_double_vector(samples[n]->x);
        free_double_vector(samples[n]->p);
        free_double_vector(samples[n]->w);
        free(samples[n]);
    }
    free(samples);
    
    // covariance matrices for different modes
    for(size_t n=0; n<NMODE; n++) free_MVG(modes[n]);
    free(modes);

    free(LFLAG);

    
    return 0;
}



