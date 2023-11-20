/*
 *  Copyright (C) 2023 Tyson B. Littenberg (MSFC-ST12)
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

#include <lisa.h>
#include <data.h>

#include "noise.h"

void print_spline_state(struct SplineModel *model, FILE *fptr, int step)
{
    fprintf(fptr,"%i %.12g %i\n",step,model->logL,model->spline->N);
}

void print_instrument_state(struct InstrumentModel *model, FILE *fptr, int step)
{
    fprintf(fptr,"%i %.12g ",step,model->logL);
    for(int i=0; i<model->Nlink; i++)fprintf(fptr,"%.12g ", model->sacc[i]);
    for(int i=0; i<model->Nlink; i++)fprintf(fptr,"%.12g ", model->soms[i]);
    fprintf(fptr,"\n");
}

void print_foreground_state(struct ForegroundModel *model, FILE *fptr, int step)
{
    fprintf(fptr,"%i %.12g ",step,model->logL);
    for(int i=0; i<model->Nparams; i++)fprintf(fptr,"%.12g ", model->sgal[i]);
    fprintf(fptr,"\n");
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
    
    for(int i=0; i<data->N; i++)
    {
        gsl_sort(data->S_pow[i][0],1,data->Nwave);
        gsl_sort(data->S_pow[i][1],1,data->Nwave);
        
        
        double f = (double)(i+data->qmin)/data->T;
        fprintf(fptr_Snf,"%.12g ",f);
        
        for(int j=0; j<data->Nchannel; j++)
        {
            double S_med   = gsl_stats_median_from_sorted_data   (data->S_pow[i][j], 1, data->Nwave);
            double S_lo_50 = gsl_stats_quantile_from_sorted_data (data->S_pow[i][j], 1, data->Nwave, 0.25);
            double S_hi_50 = gsl_stats_quantile_from_sorted_data (data->S_pow[i][j], 1, data->Nwave, 0.75);
            double S_lo_90 = gsl_stats_quantile_from_sorted_data (data->S_pow[i][j], 1, data->Nwave, 0.05);
            double S_hi_90 = gsl_stats_quantile_from_sorted_data (data->S_pow[i][j], 1, data->Nwave, 0.95);
            
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
