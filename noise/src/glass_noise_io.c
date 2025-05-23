/*
 * Copyright 2023 Tyson B. Littenberg
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


#include <glass_utils.h>

#include "glass_noise.h"


void print_spline_state(struct SplineModel *model, FILE *fptr, int step)
{
    fprintf(fptr,"%i %.12g %i\n",step,model->logL,model->spline->N);
}

void print_instrument_state(struct InstrumentModel *model, FILE *fptr)
{
    for(int i=0; i<model->Nlink; i++)fprintf(fptr,"%.12g ", model->sacc[i]);
    for(int i=0; i<model->Nlink; i++)fprintf(fptr,"%.12g ", model->soms[i]);
}

void print_foreground_state(struct ForegroundModel *model, FILE *fptr)
{
    for(int i=0; i<model->Nparams; i++)fprintf(fptr,"%.12g ", model->sgal[i]);
}

void print_noise_model(struct Noise *noise, char filename[])
{
    FILE *fptr = fopen(filename,"w");
    for(int i=0; i<noise->N; i++)
    {
        fprintf(fptr,"%lg ",noise->f[i]);
        for(int j=0; j<noise->Nchannel; j++)
            fprintf(fptr,"%lg ",noise->C[j][j][i]);
        fprintf(fptr,"%lg ",noise->C[0][1][i]);
        fprintf(fptr,"%lg ",noise->C[0][2][i]);
        fprintf(fptr,"%lg ",noise->C[1][2][i]);
        fprintf(fptr,"\n");
    }
    fclose(fptr);
    
}

void print_whitened_data(struct Data *data, struct Noise *noise, char filename[])
{
    FILE *fptr = fopen(filename,"w");
    for(int i=0; i<noise->N; i++)
    {
        fprintf(fptr,"%lg ",noise->f[i]);
        fprintf(fptr,"%lg %lg ",data->tdi->X[2*i]/sqrt(noise->C[0][0][i]),data->tdi->X[2*i+1]/sqrt(noise->C[0][0][i]));
        fprintf(fptr,"%lg %lg ",data->tdi->Y[2*i]/sqrt(noise->C[1][1][i]),data->tdi->Y[2*i+1]/sqrt(noise->C[1][1][i]));
        fprintf(fptr,"%lg %lg ",data->tdi->Z[2*i]/sqrt(noise->C[2][2][i]),data->tdi->Z[2*i+1]/sqrt(noise->C[2][2][i]));
        fprintf(fptr,"\n");
    }
        
    fclose(fptr);
}

void print_noise_reconstruction(struct Data *data, struct Flags *flags)
{
    FILE *fptr_Snf;
    char filename[MAXSTRINGSIZE];
    
    sprintf(filename,"%s/power_noise_reconstruction.dat",data->dataDir);
    fptr_Snf=fopen(filename,"w");
    
    for(int i=0; i<data->NFFT; i++)
    {
        double f = (double)(i+data->qmin)/data->T;
        fprintf(fptr_Snf,"%.12g ",f);
        
        for(int j=0; j<data->Nchannel; j++)
        {
            double_sort(data->S_pow[i][j],data->Nwave);

            double S_med   = get_quantile_from_sorted_data(data->S_pow[i][j], data->Nwave, 0.50);
            double S_lo_50 = get_quantile_from_sorted_data(data->S_pow[i][j], data->Nwave, 0.25);
            double S_hi_50 = get_quantile_from_sorted_data(data->S_pow[i][j], data->Nwave, 0.75);
            double S_lo_90 = get_quantile_from_sorted_data(data->S_pow[i][j], data->Nwave, 0.05);
            double S_hi_90 = get_quantile_from_sorted_data(data->S_pow[i][j], data->Nwave, 0.95);
            
            fprintf(fptr_Snf,"%lg ",S_med);
            fprintf(fptr_Snf,"%lg ",S_lo_50);
            fprintf(fptr_Snf,"%lg ",S_hi_50);
            fprintf(fptr_Snf,"%lg ",S_lo_90);
            fprintf(fptr_Snf,"%lg ",S_hi_90);
        }
        fprintf(fptr_Snf,"\n");
    }
    fclose(fptr_Snf);
    
}
