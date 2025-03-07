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

#include <glass_utils.h>

#include "glass_noise.h"


void map_array_to_noise_params(struct InstrumentModel *model)
{
    model->sacc12 = model->sacc[0];
    model->sacc21 = model->sacc[1];

    model->sacc23 = model->sacc[2];
    model->sacc32 = model->sacc[3];

    model->sacc13 = model->sacc[4];
    model->sacc31 = model->sacc[5];

    model->soms12 = model->soms[0];
    model->soms21 = model->soms[1];
        
    model->soms23 = model->soms[2];
    model->soms32 = model->soms[3];
    
    model->soms13 = model->soms[4];
    model->soms31 = model->soms[5];

}

void map_noise_params_to_array(struct InstrumentModel *model)
{
    model->sacc[0] = model->sacc12;
    model->sacc[1] = model->sacc21;
        
    model->sacc[2] = model->sacc23;
    model->sacc[3] = model->sacc32;
    
    model->sacc[4] = model->sacc13;
    model->sacc[5] = model->sacc31;

    model->soms[0] = model->soms12;
    model->soms[1] = model->soms21;
    
    model->soms[2] = model->soms23;
    model->soms[3] = model->soms32;

    model->soms[4] = model->soms13;
    model->soms[5] = model->soms31;
}

void map_array_to_foreground_params(struct ForegroundModel *model)
{
    model->Amp   = exp(model->sgal[0]);
    model->f1    = exp(model->sgal[1]);
    model->alpha = model->sgal[2];
    model->fk    = exp(model->sgal[3]);
    model->f2    = exp(model->sgal[4]);
}

void map_foreground_params_to_array(struct ForegroundModel *model)
{
    model->sgal[0] = log(model->Amp);
    model->sgal[1] = log(model->f1);
    model->sgal[2] = model->alpha;
    model->sgal[3] = log(model->fk);
    model->sgal[4] = log(model->f2);
}

void alloc_spline_model(struct SplineModel *model, int Ndata, int Nlayer, int Nchannel, int Nspline)
{
    model->Nchannel = Nchannel;
    
    model->psd = malloc(sizeof(struct Noise));
    model->spline=malloc(sizeof(struct Noise));
    
    alloc_noise(model->psd, Ndata, Nlayer, Nchannel);
    alloc_noise(model->spline, Nspline, Nlayer, Nchannel);
}

void alloc_instrument_model(struct InstrumentModel *model, int Ndata, int Nlayer, int Nchannel)
{
    model->Nlink=6;
    model->soms = malloc(model->Nlink*sizeof(double));
    model->sacc = malloc(model->Nlink*sizeof(double));
    model->psd = malloc(sizeof(struct Noise));
    alloc_noise(model->psd, Ndata, Nlayer, Nchannel);
}

void alloc_foreground_model(struct ForegroundModel *model, int Ndata, int Nlayer, int Nchannel)
{
    model->Nparams=5;
    model->sgal = malloc(model->Nparams*sizeof(double));
    model->psd = malloc(sizeof(struct Noise));
    alloc_noise(model->psd, Ndata, Nlayer, Nchannel);
}

void free_spline_model(struct SplineModel *model)
{
    free_noise(model->psd);
    free_noise(model->spline);
    free(model);
}

void free_instrument_model(struct InstrumentModel *model)
{
    free(model->soms);
    free(model->sacc);
    free_noise(model->psd);
    free(model);
}

void free_foreground_model(struct ForegroundModel *model)
{
    free_noise(model->psd);
    free(model->sgal);
    free(model);
}

void copy_spline_model(struct SplineModel *origin, struct SplineModel *copy)
{
    copy->Nchannel = origin->Nchannel;

    //Spline parameters
    copy_noise(origin->spline,copy->spline);
    
    //Noise model parameters
    copy_noise(origin->psd,copy->psd);
    
    //Model likelihood
    copy->logL = origin->logL;
    
    //Priors
    copy->Nmin = origin->Nmin;
    copy->Nmax = origin->Nmax;
}

void copy_instrument_model(struct InstrumentModel *origin, struct InstrumentModel *copy)
{
    //Noise model parameters
    copy_noise(origin->psd,copy->psd);

    //Model likelihood
    copy->logL = origin->logL;

    //Instrument model Parameters
    copy->Nlink = origin->Nlink;
    memcpy(copy->soms, origin->soms, origin->Nlink*sizeof(double));
    memcpy(copy->sacc, origin->sacc, origin->Nlink*sizeof(double));
}

void copy_foreground_model(struct ForegroundModel *origin, struct ForegroundModel *copy)
{
    //Noise model parameters
    copy_noise(origin->psd,copy->psd);
    
    //Model likelihood
    copy->logL = origin->logL;

    //Foreground model Parameters
    copy->Tobs = origin->Tobs;
    copy->Nparams = origin->Nparams;
    memcpy(copy->sgal, origin->sgal, origin->Nparams*sizeof(double));
    
}

void update_spline_noise_model(struct SplineModel *model, int new_knot, int min_knot, int max_knot)
{
    struct Noise *psd = model->psd;
    struct Noise *spline = model->spline;
    
    /* find location in data vector of knot */
    double T = 1./(psd->f[1] - psd->f[0]);

    gsl_spline **cspline = malloc(model->Nchannel*sizeof(gsl_spline *));
    gsl_interp_accel **acc = malloc(model->Nchannel*sizeof(gsl_interp_accel *));
    
    /* have to recompute the spline everywhere (derivatives on boundary) */
    for(int n=0; n<model->Nchannel; n++)
    {
        cspline[n] = gsl_spline_alloc(gsl_interp_akima, spline->N);
        acc[n] = gsl_interp_accel_alloc();
        gsl_spline_init(cspline[n],spline->f,spline->C[n][n],spline->N);
    }
    
    int imin = (int)((spline->f[min_knot]-psd->f[0])*T);
    int imax = (int)((spline->f[max_knot]-psd->f[0])*T);
    
    
    for(int i=imin; i<imax; i++)
    {
        for(int n=0; n<model->Nchannel; n++)
        {
            psd->C[n][n][i]=gsl_spline_eval(cspline[n],psd->f[i],acc[n]);
            
            /*
             apply transfer function
             -this catches the sharp features in the spectrum from f/fstar
             -without needing to interpolate
             */
            psd->C[n][n][i]+=psd->transfer[i];
        }
    }
    invert_noise_covariance_matrix(psd);

    for(int n=0; n<model->Nchannel; n++)
    {
        gsl_spline_free(cspline[n]);
        gsl_interp_accel_free(acc[n]);
    }
    free(cspline);
    free(acc);

}


void generate_spline_noise_model(struct SplineModel *model)
{
    struct Noise *psd = model->psd;
    struct Noise *spline = model->spline;
    
    for(int i=0; i<model->Nchannel; i++)
        CubicSplineGSL(spline->N, spline->f, spline->C[i][i], psd->N, psd->f, psd->C[i][i]);

    //set a floor on Sn so likelihood doesn't go crazy
    for(int n=0; n<psd->N; n++)
    {
        /*
         apply transfer function
         -this catches the sharp features in the spectrum from f/fstar
         -without needing to interpolate
         */
        for(int i=0; i<model->Nchannel; i++)
            psd->C[i][i][n]+=psd->transfer[n];
        
    }
    invert_noise_covariance_matrix(psd);
}

void generate_instrument_noise_model(struct Orbit *orbit, struct InstrumentModel *model)
{
    double f,f2;
    double x;
    double cosx;
    double oms_transfer_function;
    double acc_transfer_function;
    double tdi_transfer_function;
    double Sacc;
    double Soms;
    double Sacc12,Sacc21,Sacc13,Sacc31,Sacc23,Sacc32;
    double Soms12,Soms21,Soms13,Soms31,Soms23,Soms32;

    map_array_to_noise_params(model);

    double acc_units = 2.81837551648e-19;//1./(PI2*CLIGHT)/(PI2*CLIGHT)
    double oms_units = 4.39256635604e-16;//PI2*PI2/CLIGHT/CLIGHT
    double acc_fmax  = 0.01; //Hz (maximum frequency for computing acc contribution
    
    for(int n=0; n<model->psd->N; n++)
    {
        f = model->psd->f[n];
        x = f/orbit->fstar;
        f2= f*f;
        
        //at low frequency use linear approximation for trig functions
        if(x<0.1)
        {
            cosx = 1.0;
            tdi_transfer_function = x*x;
        }
        else
        {
            cosx = cos(x);
            tdi_transfer_function = noise_transfer_function(x);
        }

        cosx = cos(x);
        tdi_transfer_function = noise_transfer_function(x);

        if(f<acc_fmax) //at f<~.1mHz acceleration noise is < 1% of the total noise budget
            acc_transfer_function = acc_units / f2 * (1.0 + pow(0.4e-3/f,2)) * (1.0 + pow(f/8.0e-3,4));
        else
            acc_transfer_function = 0.0;

        oms_transfer_function = oms_units * f2 * (1.0 + pow(2.0e-3/f,4));

        //absorb factor of -8 into tdi transfer function to cut down on multiplications
        tdi_transfer_function *= -8.0;
        
        switch(model->psd->Nchannel)
        {
            case 1:
                Sacc = model->sacc12 * acc_transfer_function;
                Soms = model->soms12 * oms_transfer_function;
                model->psd->C[0][0][n] = XYZnoise_FF(orbit->L, orbit->fstar, f, Sacc, Soms);
                break;
            case 2:
                Sacc = model->sacc12 * acc_transfer_function;
                Soms = model->soms12 * oms_transfer_function;
                model->psd->C[0][0][n] = AEnoise_FF(orbit->L, orbit->fstar, f, Sacc, Soms);
                model->psd->C[1][1][n] = AEnoise_FF(orbit->L, orbit->fstar, f, Sacc, Soms);
                model->psd->C[0][1][n] = 0.0;
                break;
            case 3:
                
                if(f<acc_fmax) //at f<~.1mHz acceleration noise is < 1% of the total noise budget
                {
                    Sacc12 = model->sacc12 * acc_transfer_function;
                    Sacc21 = model->sacc21 * acc_transfer_function;
                    Sacc23 = model->sacc23 * acc_transfer_function;
                    Sacc32 = model->sacc32 * acc_transfer_function;
                    Sacc13 = model->sacc13 * acc_transfer_function;
                    Sacc31 = model->sacc31 * acc_transfer_function;
                }
                
                Soms12 = model->soms12 * oms_transfer_function;
                Soms21 = model->soms21 * oms_transfer_function;
                Soms23 = model->soms23 * oms_transfer_function;
                Soms32 = model->soms32 * oms_transfer_function;
                Soms13 = model->soms13 * oms_transfer_function;
                Soms31 = model->soms31 * oms_transfer_function;
            
                //Start from scratch
                model->psd->C[0][0][n] = 0.;
                model->psd->C[1][1][n] = 0.;
                model->psd->C[2][2][n] = 0.;

                model->psd->C[0][1][n] = 0.;
                model->psd->C[0][2][n] = 0.;
                model->psd->C[1][2][n] = 0.;

                //OMS noise
                model->psd->C[0][0][n] += (Soms12+Soms21+Soms13+Soms31)*.25;
                model->psd->C[1][1][n] += (Soms23+Soms32+Soms21+Soms12)*.25;
                model->psd->C[2][2][n] += (Soms31+Soms13+Soms32+Soms23)*.25;

                model->psd->C[0][1][n] += (Soms12 + Soms21)*.5;
                model->psd->C[0][2][n] += (Soms13 + Soms31)*.5;
                model->psd->C[1][2][n] += (Soms32 + Soms23)*.5;
                
                //Acceleration noise
                if(f<acc_fmax) //at f<~.1mHz acceleration noise is < 1% of the total noise budget
                {
                    //distribute through factors of 2
                    /*
                    model->psd->C[0][0][n] += 2. * ( (Sacc12 + Sacc13)*.5 + ((Sacc21 + Sacc31)*.5)*cosx*cosx );
                    model->psd->C[1][1][n] += 2. * ( (Sacc23 + Sacc21)*.5 + ((Sacc32 + Sacc12)*.5)*cosx*cosx );
                    model->psd->C[2][2][n] += 2. * ( (Sacc31 + Sacc32)*.5 + ((Sacc13 + Sacc23)*.5)*cosx*cosx );

                    model->psd->C[0][1][n] += 4. * ( (Sacc12 + Sacc21)*.5 );
                    model->psd->C[0][2][n] += 4. * ( (Sacc13 + Sacc31)*.5 );
                    model->psd->C[1][2][n] += 4. * ( (Sacc32 + Sacc23)*.5 );
                    */
                    model->psd->C[0][0][n] += ( (Sacc12 + Sacc13) + ((Sacc21 + Sacc31))*cosx*cosx );
                    model->psd->C[1][1][n] += ( (Sacc23 + Sacc21) + ((Sacc32 + Sacc12))*cosx*cosx );
                    model->psd->C[2][2][n] += ( (Sacc31 + Sacc32) + ((Sacc13 + Sacc23))*cosx*cosx );

                    model->psd->C[0][1][n] += 2.*(Sacc12 + Sacc21);
                    model->psd->C[0][2][n] += 2.*(Sacc13 + Sacc31);
                    model->psd->C[1][2][n] += 2.*(Sacc32 + Sacc23);
                }
                

                //TDI transfer functions (note tdi_transfer_function has abosrbed a factor of -8)
                model->psd->C[0][0][n] *= -2. * tdi_transfer_function;
                model->psd->C[1][1][n] *= -2. * tdi_transfer_function;
                model->psd->C[2][2][n] *= -2. * tdi_transfer_function;

                model->psd->C[0][1][n] *= tdi_transfer_function * cosx ;
                model->psd->C[0][2][n] *= tdi_transfer_function * cosx ;
                model->psd->C[1][2][n] *= tdi_transfer_function * cosx ;
                
                //Symmetry
                model->psd->C[1][0][n] = model->psd->C[0][1][n];
                model->psd->C[2][0][n] = model->psd->C[0][2][n];
                model->psd->C[2][1][n] = model->psd->C[1][2][n];
                
                break;
        }
        
        //Normalization
        for(int i=0; i<model->psd->Nchannel; i++)
            for(int j=0; j<model->psd->Nchannel; j++)
                model->psd->C[i][j][n] /= 2.0;
    }
}

void generate_instrument_noise_model_wavelet(struct Wavelets *wdm, struct Orbit *orbit, struct InstrumentModel *model)
{

    /* 
    oversampled frequency grid 
    */
    struct InstrumentModel *grid = malloc(sizeof(struct InstrumentModel));

    // active layers
    int imin = (int)round(model->psd->f[0]/wdm->df);
    int imax = (int)round(model->psd->f[model->psd->N-1]/wdm->df)+1;

    // initialize data models
    alloc_instrument_model(grid, 2*wdm->NF, imax-imin, 3);

    // set up psd frequency grid
    for(int n=0; n<grid->psd->N; n++)
        grid->psd->f[n] = wdm->df/2.0*(n+1);

    // initialize noise levels
    for(int i=0; i<grid->Nlink; i++)
    {
        grid->soms[i] = model->soms[i];
        grid->sacc[i] = model->sacc[i];
    }

    // get noise covariance matrix for initial parameters
    generate_instrument_noise_model(orbit,grid);
    
    /*
    integrate instrument noise over each frequency layer
    */
    double ***C     = model->psd->C;
    double ***Cgrid = grid->psd->C;
    for(int i=imin; i<imax; i++)
    {
        int j = 2*i-2;
        for(int n=0; n<3; n++)
            for(int m=n; m<3; m++)
                C[n][m][i-imin] = simpson_integration_3(Cgrid[n][m][j],Cgrid[n][m][j+1],Cgrid[n][m][j+2],1.0);
    }

    //NOTE: normalization fudge factor
    for(int i=0; i<model->psd->N; i++)
        for(int n=0; n<3; n++)
            for(int m=n; m<3; m++)
                C[n][m][i]/=8.;

    free_instrument_model(grid);

}

void generate_galactic_foreground_model(struct ForegroundModel *model)
{
    double f;
    double Sgal;
    
    map_array_to_foreground_params(model);
    
    for(int n=0; n<model->psd->N; n++)
    {
        f = model->psd->f[n];
        double fstar=CLIGHT/(PI2*LARM);
        double x = f/fstar;
        double t = 4.*x*x*sin(x)*sin(x); //transfer function

        //skip confusion noise calculation at high frequencies where contribution is negligible
        int high_f_check = 0;
        if(f>model->f1*10) high_f_check = 1;
        
        if(high_f_check) Sgal = 0.0;
        else Sgal = galaxy_foreground(f,model->Amp, model->f1, model->alpha, model->fk, model->f2);
        
        switch(model->psd->Nchannel)
        {
            case 1:
                model->psd->C[0][0][n] = Sgal;
                break;
            case 2:
                model->psd->C[0][0][n] = model->psd->C[1][1][n] = 1.5*t*Sgal;
                model->psd->C[0][1][n] = model->psd->C[1][0][n] = 0;
                break;
            case 3:
                model->psd->C[0][0][n] = model->psd->C[1][1][n] = model->psd->C[2][2][n] = t*Sgal;
                model->psd->C[0][1][n] = model->psd->C[0][2][n] = model->psd->C[1][2][n] = -0.5*t*Sgal;
                model->psd->C[1][0][n] = model->psd->C[2][0][n] = model->psd->C[2][1][n] = -0.5*t*Sgal;
                break;
        }
    }
}

void generate_galactic_foreground_model_wavelet(struct Wavelets *wdm, struct ForegroundModel *model)
{
    /* 
    oversampled frequency grid 
    */
    struct ForegroundModel *grid = malloc(sizeof(struct ForegroundModel));

    // active layers
    int imin = (int)round(model->psd->f[0]/wdm->df);
    int imax = (int)round(model->psd->f[model->psd->N-1]/wdm->df)+1;

    // initialize data models
    alloc_foreground_model(grid, 2*wdm->NF, imax-imin, 3);

    // set up psd frequency grid
    for(int n=0; n<grid->psd->N; n++)
        grid->psd->f[n] = wdm->df/2.0*(n+1);

    // initialize foreground parameters levels
    grid->Tobs  = model->Tobs;
    grid->Amp   = model->Amp;
    grid->f1    = model->f1;
    grid->alpha = model->alpha;
    grid->fk    = model->fk;
    grid->f2    = model->f2;
    map_foreground_params_to_array(grid);

    // get noise covariance matrix for initial parameters
    generate_galactic_foreground_model(grid);

    /*
    integrate foreground over each frequency layer
    */
    double ***C     = model->psd->C;
    double ***Cgrid = grid->psd->C;

    for(int i=imin; i<imax; i++)
    {
        int j = 2*i-2;
        for(int n=0; n<3; n++)
            for(int m=n; m<3; m++)
                C[n][m][i-imin] = simpson_integration_3(Cgrid[n][m][j],Cgrid[n][m][j+1],Cgrid[n][m][j+2],1.0);
    }

    //NOTE: undo isotropc -1/2 on covariance hardcoded in generate_galactic_foreground_model()
    for(int i=0; i<model->psd->N; i++)
        for(int n=0; n<3; n++)
            for(int m=n; m<3; m++)
                if(n!=m) C[n][m][i]*=-2.;


    //NOTE: normalization fudge factor
    for(int i=0; i<model->psd->N; i++)
        for(int n=0; n<3; n++)
            for(int m=n; m<3; m++)
                C[n][m][i]/=8.;



    free_foreground_model(grid);
}

void generate_full_covariance_matrix(struct Noise *full, struct Noise *component, int Nchannel)
{
    for(int n=0; n<full->N; n++)
    {
        for(int i=0; i<Nchannel; i++)
        {
            for(int j=0; j<Nchannel; j++)
            {
                full->C[i][j][n] += component->C[i][j][n];
            }
        }
    }
}

void generate_full_dynamic_covariance_matrix(struct Wavelets *wdm, struct InstrumentModel *inst, struct ForegroundModel *conf, struct Noise *full)
{
    int k;
    int jmin=(int)round(inst->psd->f[0]/wdm->df);
    int jmax=(int)round(inst->psd->f[inst->psd->N-1]/wdm->df)+1; 

    for(int i=0; i<wdm->NT; i++)
    {
        double t = i*wdm->dt;
        for(int j=jmin; j<jmax; j++)
        {
            wavelet_pixel_to_index(wdm,i,j,&k);

            k-=wdm->kmin;

            //stationary instrument noise 
            full->C[0][0][k] = inst->psd->C[0][0][j-jmin];
            full->C[1][1][k] = inst->psd->C[1][1][j-jmin];
            full->C[2][2][k] = inst->psd->C[2][2][j-jmin];
            full->C[0][1][k] = inst->psd->C[0][1][j-jmin];
            full->C[0][2][k] = inst->psd->C[0][2][j-jmin];
            full->C[1][2][k] = inst->psd->C[1][2][j-jmin];

            //modulated galactic foreground
            full->C[0][0][k] += conf->psd->C[0][0][j-jmin]*gsl_spline_eval(conf->modulation->XX_spline, t, conf->modulation->acc);
            full->C[1][1][k] += conf->psd->C[1][1][j-jmin]*gsl_spline_eval(conf->modulation->YY_spline, t, conf->modulation->acc);
            full->C[2][2][k] += conf->psd->C[2][2][j-jmin]*gsl_spline_eval(conf->modulation->ZZ_spline, t, conf->modulation->acc);
            full->C[0][1][k] += conf->psd->C[0][1][j-jmin]*gsl_spline_eval(conf->modulation->XY_spline, t, conf->modulation->acc); 
            full->C[0][2][k] += conf->psd->C[0][2][j-jmin]*gsl_spline_eval(conf->modulation->XZ_spline, t, conf->modulation->acc); 
            full->C[1][2][k] += conf->psd->C[1][2][j-jmin]*gsl_spline_eval(conf->modulation->YZ_spline, t, conf->modulation->acc); 

            //noise covariance matrix is symmetric
            full->C[1][0][k] = full->C[0][1][k]; 
            full->C[2][0][k] = full->C[0][2][k]; 
            full->C[2][1][k] = full->C[1][2][k]; 
        } //loop over frequency layers
    } //loop over time slices
}

static void generate_full_stationary_covariance_matrix(struct Wavelets *wdm, struct InstrumentModel *inst, struct ForegroundModel *conf, struct Noise *full)
{
    int k;
    int jmin=(int)round(inst->psd->f[0]/wdm->df);
    int jmax=(int)round(inst->psd->f[inst->psd->N-1]/wdm->df)+1;

    for(int i=0; i<wdm->NT; i++)
    {
        for(int j=jmin; j<jmax; j++)
        {
            wavelet_pixel_to_index(wdm,i,j,&k);

            k-=wdm->kmin;

            //stationary instrument noise
            full->C[0][0][k] = inst->psd->C[0][0][j-jmin];
            full->C[1][1][k] = inst->psd->C[1][1][j-jmin];
            full->C[2][2][k] = inst->psd->C[2][2][j-jmin];
            full->C[0][1][k] = inst->psd->C[0][1][j-jmin];
            full->C[0][2][k] = inst->psd->C[0][2][j-jmin];
            full->C[1][2][k] = inst->psd->C[1][2][j-jmin];

            //modulated galactic foreground
            full->C[0][0][k] += conf->psd->C[0][0][j-jmin];
            full->C[1][1][k] += conf->psd->C[1][1][j-jmin];
            full->C[2][2][k] += conf->psd->C[2][2][j-jmin];
            full->C[0][1][k] -= conf->psd->C[0][1][j-jmin]/2.;
            full->C[0][2][k] -= conf->psd->C[0][2][j-jmin]/2.;
            full->C[1][2][k] -= conf->psd->C[1][2][j-jmin]/2.;

            //noise covariance matrix is symmetric
            full->C[1][0][k] = full->C[0][1][k];
            full->C[2][0][k] = full->C[0][2][k];
            full->C[2][1][k] = full->C[1][2][k];
        }// loop over frequency layers
    }// loop over time slices
}

double noise_log_likelihood(struct Data *data, struct Noise *noise)
{
    double logL = 0.0;
    
    struct TDI *tdi = data->tdi;
    
    int N = data->NFFT;
    
    switch(data->Nchannel)
    {
        case 1:
            logL += -0.5*fourier_nwip(tdi->X, tdi->X, noise->invC[0][0], N);
            break;
        case 2:
            logL += -0.5*fourier_nwip(tdi->A, tdi->A, noise->invC[0][0], N);
            logL += -0.5*fourier_nwip(tdi->E, tdi->E, noise->invC[1][1], N);
            break;
        case 3:
            logL += -0.5*fourier_nwip(tdi->X, tdi->X, noise->invC[0][0], N);
            logL += -0.5*fourier_nwip(tdi->Y, tdi->Y, noise->invC[1][1], N);
            logL += -0.5*fourier_nwip(tdi->Z, tdi->Z, noise->invC[2][2], N);
            logL += -fourier_nwip(tdi->X, tdi->Y, noise->invC[0][1], N);
            logL += -fourier_nwip(tdi->X, tdi->Z, noise->invC[0][2], N);
            logL += -fourier_nwip(tdi->Y, tdi->Z, noise->invC[1][2], N);
            break;
    }
    for(int n=0; n<N; n++)
        logL -= log(noise->detC[n]);
    
    return logL;
}

double noise_delta_log_likelihood(struct Data *data, struct SplineModel *model_x, struct SplineModel *model_y, double fmin, double fmax,int ic)
{
    double dlogL = 0.0;
    
    struct TDI *tdi = data->tdi;
    struct Noise *psd_x = model_x->psd;
    struct Noise *psd_y = model_y->psd;

    int N = (int)floor((fmax-fmin)*data->T);
    int imin = (int)floor((fmin-psd_x->f[0])*data->T);
    if(imin<0)imin=0;
    
    /* remove contribution for current state x */
    dlogL -= -0.5*fourier_nwip(tdi->A+2*imin, tdi->A+2*imin, psd_x->invC[0][0]+imin, N);
    dlogL -= -0.5*fourier_nwip(tdi->E+2*imin, tdi->E+2*imin, psd_x->invC[1][1]+imin, N);
    for(int n=imin; n<imin+N; n++)
        dlogL += log(psd_x->detC[n]);

    /* add contribution for proposed state y */
    dlogL += -0.5*fourier_nwip(tdi->A+2*imin, tdi->A+2*imin, psd_y->invC[0][0]+imin, N);
    dlogL += -0.5*fourier_nwip(tdi->E+2*imin, tdi->E+2*imin, psd_y->invC[1][1]+imin, N);
    for(int n=imin; n<imin+N; n++)
        dlogL -= log(psd_y->detC[n]);
    
    return dlogL;
}

void initialize_spline_model(struct Orbit *orbit, struct Data *data, struct SplineModel *model, int Nspline)
{
    
    // Initialize data models
    alloc_spline_model(model, data->NFFT, data->Nlayer, data->Nchannel, Nspline);
    
    //set max and min spline points
    model->Nmin = MIN_SPLINE_STENCIL;
    model->Nmax = Nspline;
    model->Nchannel = data->Nchannel;
    
    //set up psd frequency grid
    for(int n=0; n<model->psd->N; n++)
    {
        double f = data->fmin + (double)n/data->T;
        double Spm, Sop;
        get_noise_levels("sangria", f, &Spm, &Sop);
        model->psd->f[n] = f;
        if(model->Nchannel==2) model->psd->transfer[n] = AEnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);//noise_transfer_function(f/orbit->fstar);
        else model->psd->transfer[n] = XYZnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);
    }
    
    //divide into Nspline control points
    double logdf = (log(data->fmax) - log(data->fmin))/(Nspline-1);
    for(int i=0; i<Nspline; i++)
    {
        double f = exp(log(data->fmin) + (double)i*logdf);
        model->spline->f[i] = f;
        
        for(int n=0; n<model->Nchannel; n++)
        {
            /* initialize model to theoretical level without transfer function applied */
            for(int m=0; m<model->Nchannel; m++)
                model->spline->C[n][m][i] = 0.0;
        }
    }
    //shift first spline control point by half a bin to avoid rounding problems
    model->spline->f[0] -= 0.5/data->T;
    
    generate_spline_noise_model(model);
    model->logL = noise_log_likelihood(data, model->psd);
}

void initialize_instrument_model(struct Orbit *orbit, struct Data *data, struct InstrumentModel *model)
{
    // initialize data models
    alloc_instrument_model(model, data->NFFT, data->Nlayer, data->Nchannel);
    
    // set up psd frequency grid
    for(int n=0; n<model->psd->N; n++)
        model->psd->f[n] = data->fmin + (double)n/data->T;

    // initialize noise levels
    for(int i=0; i<model->Nlink; i++)
    {
        model->soms[i] = 2.25e-22;
        model->sacc[i] = 9.00e-30;
    }
    
    // get noise covariance matrix for initial parameters
    generate_instrument_noise_model(orbit,model);
}

void initialize_instrument_model_wavelet(struct Orbit *orbit, struct Data *data, struct InstrumentModel *model)
{
    // wavelet basis
    struct Wavelets *wdm = data->wdm;

    // initialize data models
    alloc_instrument_model(model, data->qmax-data->qmin, data->Nlayer, data->Nchannel);
    
    // set up psd frequency grid
    for(int n=0; n<model->psd->N; n++)
        model->psd->f[n] = (data->qmin+n)*wdm->df;

    // initialize noise levels
    for(int i=0; i<model->Nlink; i++)
    {
        model->soms[i] = 1.28e-22;
        model->sacc[i] = 5.76e-30;
    }
    
    // get noise covariance matrix for initial parameters
    generate_instrument_noise_model_wavelet(wdm,orbit,model);
}

void initialize_foreground_model(struct Orbit *orbit, struct Data *data, struct ForegroundModel *model)
{
    // initialize data models
    alloc_foreground_model(model, data->NFFT, data->Nlayer, data->Nchannel);
    
    // set up psd frequency grid
    for(int n=0; n<model->psd->N; n++)
        model->psd->f[n] = data->fmin + (double)n/data->T;

    // initialize constant foreground parameters levels 
    model->Tobs  =  data->T;
    model->Amp   =  1.2826e-44;
    model->alpha =  1.629667;
    model->f2    =  4.810781e-4;

    // initialize time dependent foreground parameter elves
    double af1 = -2.235e-1;
    double bf1 = -2.7040844;
    double afk = -3.60976122e-1;
    double bfk = -2.37822436;

    model->f1  =  pow(10., af1*log10(model->Tobs/YEAR) + bf1);
    model->fk  =  pow(10., afk*log10(model->Tobs/YEAR) + bfk);

    map_foreground_params_to_array(model);
    
    // get noise covariance matrix for initial parameters
    generate_galactic_foreground_model(model);

}

void initialize_foreground_model_wavelet(struct Orbit *orbit, struct Data *data, struct ForegroundModel *model)
{
    //wavelet basis
    struct Wavelets *wdm = data->wdm;

    // initialize data models
    alloc_foreground_model(model, data->qmax-data->qmin, data->Nlayer, data->Nchannel);
    
    // set up psd frequency grid
    for(int n=0; n<model->psd->N; n++)
        model->psd->f[n] = (data->qmin+n)*wdm->df;

    // initialize constant foreground parameters levels 
    model->Tobs  =  data->T;
    model->Amp   =  1.2826e-44;
    model->alpha =  1.629667;
    model->f2    =  4.810781e-4;

    // initialize time dependent foreground parameter elves
    double af1 = -2.235e-1;
    double bf1 = -2.7040844;
    double afk = -3.60976122e-1;
    double bfk = -2.37822436;

    model->f1  =  pow(10., af1*log10(model->Tobs/YEAR) + bf1);
    model->fk  =  pow(10., afk*log10(model->Tobs/YEAR) + bfk);
    map_foreground_params_to_array(model);
    
    // get noise covariance matrix for initial parameters
    generate_galactic_foreground_model_wavelet(wdm,model);

        // get galaxy modulation
    model->modulation = malloc(sizeof(struct GalaxyModulation));
    initialize_galaxy_modulation(model->modulation, data->wdm, orbit, data->T, data->t0);
    
    /**************************************************
     * Compute galaxy modulation
    **************************************************/
    
    double *galaxy_params = double_vector(6); // defines galaxy shape
    galaxy_params[0] = 0.25; // A 0.25    bulge fraction
    galaxy_params[1] = 0.8;  // Rb 0.8    bulge radius (kpc)
    galaxy_params[2] = 2.5;  // Rd 2.5    disk radius (kpc)
    galaxy_params[3] = 0.4;  // Zd 0.4    disk height (kpc)
    galaxy_params[4] = 7.2;  // RGC 7.2   distance from solar BC to GC (kpc)
    galaxy_params[5] = 3.5;  // Rcut 3.5  radius out to which all sources are found (kpc)

    //computes the modulation of the confusion noise
    galaxy_modulation(model->modulation, galaxy_params);
    free_double_vector(galaxy_params);
}

void GetDynamicNoiseModel(struct Data *data, struct Orbit *orbit, struct Flags *flags)
{
    /**************************************************
     * Compute instrument noise levels
    **************************************************/
    struct InstrumentModel *inst_noise = malloc(sizeof(struct InstrumentModel));
    initialize_instrument_model_wavelet(orbit, data, inst_noise);

    /**************************************************
     * Compute galactic foreground noise levels
    **************************************************/
    struct ForegroundModel *conf_noise = malloc(sizeof(struct ForegroundModel));
    initialize_foreground_model_wavelet(orbit, data, conf_noise);

    /**************************************************
     * Combine noise components
    **************************************************/
    generate_full_dynamic_covariance_matrix(data->wdm, inst_noise, conf_noise, data->noise);
    invert_noise_covariance_matrix(data->noise);

    char filename[128];
    sprintf(filename,"%s/power_noise.dat",data->dataDir);
    FILE *fptr=fopen(filename,"w");
    int k;
    for(int j=data->qmin; j<data->qmax; j++)
    {
        double f = j*data->wdm->df;
        for(int i=0; i<data->wdm->NT; i++)
        {
            double t = i*data->wdm->dt;
            wavelet_pixel_to_index(data->wdm,i,j,&k);
            k-=data->wdm->kmin;
            fprintf(fptr,"%lg %lg %.14e\n", t, f, data->noise->C[0][0][k]);
        }
        fprintf(fptr,"\n");
    }
    fclose(fptr);

    free_instrument_model(inst_noise);
    free_foreground_model(conf_noise);
}

void GetStationaryNoiseModel(struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Noise *noise)
{
    /**************************************************
     * Compute instrument noise levels
    **************************************************/
    struct InstrumentModel *inst_noise = malloc(sizeof(struct InstrumentModel));
    initialize_instrument_model_wavelet(orbit, data, inst_noise);

    /**************************************************
     * Compute galactic foreground noise levels
    **************************************************/
    struct ForegroundModel *conf_noise = malloc(sizeof(struct ForegroundModel));
    initialize_foreground_model_wavelet(orbit, data, conf_noise);

    /**************************************************
     * Combine noise components
    **************************************************/
    generate_full_stationary_covariance_matrix(data->wdm, inst_noise, conf_noise, noise);
    invert_noise_covariance_matrix(noise);

    char filename[128];
    sprintf(filename,"%s/power_stationary_noise.dat",data->dataDir);
    FILE *fptr=fopen(filename,"w");
    int k;
    for(int j=data->qmin; j<data->qmax; j++)
    {
        double f = j*data->wdm->df;
        for(int i=0; i<data->wdm->NT; i++)
        {
            double t = i*data->wdm->dt;
            wavelet_pixel_to_index(data->wdm,i,j,&k);
            k-=data->wdm->kmin;
            fprintf(fptr,"%lg %lg %.14e\n", t, f,noise->C[0][0][k]);
        }
        fprintf(fptr,"\n");
    }
    fclose(fptr);

    free_instrument_model(inst_noise);
    free_foreground_model(conf_noise);
}
 
