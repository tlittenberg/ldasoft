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
        printf("Usage: gbmcmc_frequency_spacing Nnode Nmin Nmax qpad fmin T\n");
        return 0;
    }
    
    struct freq_data *data = malloc(sizeof(struct freq_data));
    
    //how many ucb nodes
    int N_node = atoi(argv[1]);
    
    //parse command line
    data->N = atoi(argv[2]);
    data->Nmax = atoi(argv[3]);
    data->qpad = atoi(argv[4]);
    double fstart = (double)atof(argv[5]);
    data->T = (double)atof(argv[6]);

    int procID_min = 0;
    int procID_max = N_node-1;

    //how many section sizes?
    int Smin =  (int)round(log(data->N)/log(2.));
    int Smax =  (int)round(log(data->Nmax)/log(2.));
    int N_seg = Smax - Smin + 1;
    
//    printf("N_seg=%i\n",N_seg);
    
    //integer part of nodes per section
    int n = (int)floor((double)N_node/(double)N_seg);
        
    int k[procID_max+1];
    int counter = 1;
    int counter2= 0;
    int div=4;
    int start = procID_min;
    for(int procID=procID_min; procID<=procID_max; procID++)
    {
        k[procID] = (procID - procID_min)/n;
        if(k[procID]>N_seg-1) k[procID] = N_seg-1;
        
        if( procID-start > (procID_max-start)/div )
        {
            start = procID;
            if(counter2==0)
            {
                counter = 0;
                div=2;
            }
//            else if(counter2==1)
//            {
//                counter = 0;
//                div=2;
//            }
            else counter++;
            //if(counter2>0)counter++;
            //if(counter2>0)div=2;
            counter2++;
        }
                
        k[procID]=counter;
        if(k[procID]>N_seg-1) k[procID] = N_seg-1;
        
        //printf("k[%i]=%i\n",procID,k[procID]);
    }
    //exit(1);
    
    //remainder
    int n_extra = N_node - (int)(n*N_seg);
    
    //which section am I in?
    for(int procID=procID_min; procID<=procID_max; procID++)
    {
        k[procID]=0;
        if(data->fmin <= 0.001) k[procID]=1;
        if(data->fmin > 0.001 && data->fmin <= .008) k[procID]=0;
        if(data->fmin > .008 && data->fmin <= .009) k[procID]=1;
        if(data->fmin > .009 && data->fmin <= 0.011) k[procID]=2;
        //if(data->fmin > 0.01 && data->fmin <= 0.02) k[procID]=3;
        if(data->fmin > 0.011)k[procID]=3;

        //size of nodes in my section
        data->N = (int)round(pow(2,Smin+k[procID]));
        
        //start bin of my node?
        int Nsum=0;
        for(int node=procID_min; node<procID; node++)
        {
            Nsum += (int)round(pow(2,Smin+k[node]));
        }
        data->fmin = fstart + Nsum/data->T;
        data->fmax = data->fmin + data->N/data->T;
        data->qmin = (int)(data->fmin*data->T);
        data->qmax = data->qmin+data->N;
        
        printf("%i %i %.16g %.16g\n",procID,data->N,data->fmin,data->fmax);
        
        //add padding
        data->N += 2*data->qpad;
        data->qmin -= data->qpad;
        data->qmax += data->qpad;
        data->fmin = (double)data->qmin/data->T;
        data->fmax = (double)data->qmax/data->T;
        
    }

}
