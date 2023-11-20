//
//  gb_mcmc_frequency_spacing.c
//  tools
//
//  Created by Littenberg, Tyson B. (MSFC-ZP12) on 11/2/21.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct freq_data
{
    int N;
    int Nmax;
    int qpad;
    int qmin;
    int qmax;
    double fmin;
    double fmax;
    double T;
};

int main(int argc, char* argv[])
{
    if(argc!=7)
    {
        printf("Usage: gbmcmc_frequency_spacing Nnode Nmin Nmax qpad fmax T\n");
        return 0;
    }
    
    struct freq_data *data = malloc(sizeof(struct freq_data));
    
    //how many ucb nodes of max size
    int N_node = atoi(argv[1]);
    
    //parse command line
    data->N = atoi(argv[2]);
    data->Nmax = atoi(argv[3]);
    data->qpad = atoi(argv[4]);
    double fstart;
    double fmax = (double)atof(argv[5]);
    data->T = (double)atof(argv[6]);

    int procID_min = 0;
    int procID_max = N_node-1;

    //how many sections?
    int Smin =  (int)round(log(data->N)/log(2.));
    int Smax =  (int)round(log(data->Nmax)/log(2.));
    int N_seg = Smax - Smin + 1;
    
    //number of nodes per section
    int narray[N_seg];
    int n;
    printf("segments = %i\n",N_seg);
    int N_node_temp = N_node;
    if(N_seg==1) narray[0] = N_node;
    else
    {
        for(int seg=N_seg-1; seg>=0; seg--)
        {
            N_node_temp/=2;
            printf("seg %i, size = %i\n",seg,N_node_temp);
        }
    }
    return(1);
        
    //remainder
    int n_extra = N_node - (int)(n*N_seg);
    
    //which section am I in?
    for(int procID=procID_min; procID<=procID_max; procID++)
    {
        int k = (procID - procID_min)/n;
        if(k>N_seg-1) k = N_seg-1;
        
        //size of nodes in my section
        data->N = (int)round(pow(2,Smin+k));
        
        //start bin of my node?
        int Nsum=0;
        for(int node=procID_min; node<procID; node++)
        {
            k = (node - procID_min)/n;
            if(k>N_seg-1) k = N_seg-1;
            Nsum += (int)round(pow(2,Smin+k));
        }
        data->fmin = fstart + Nsum/data->T;
        data->fmax = data->fmin + data->N/data->T;
        data->qmin = (int)(data->fmin*data->T);
        data->qmax = data->qmin+data->N;
        
        printf("Node %i: size = %i, f=[%g,%g]\n",procID,data->N,data->fmin,data->fmax);
        
        //add padding
        data->N += 2*data->qpad;
        data->qmin -= data->qpad;
        data->qmax += data->qpad;
        data->fmin = (double)data->qmin/data->T;
        data->fmax = (double)data->qmax/data->T;
        
    }

}
