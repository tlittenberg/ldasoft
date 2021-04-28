/*
 *  Copyright (C) 2019 Tyson B. Littenberg (MSFC-ST12), Neil J. Cornish
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

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <sys/stat.h>

#include <gsl/gsl_fft_complex.h>

#include <LISA.h>

#include "GalacticBinary.h"
#include "GalacticBinaryMath.h"
#include "GalacticBinaryModel.h"
#include "GalacticBinaryWaveform.h"


double galactic_binary_Amp(double Mc, double f0, double D)
{
    double f = f0;//T;
    double M = Mc*TSUN;
    double dL= D*PC/CLIGHT;
    
    return 2.*pow(pow(M,5)*pow(M_PI*f,2),1./3.)/dL;
}
double galactic_binary_fdot(double Mc, double f0)
{
    double f = f0;
    double M = Mc*TSUN;
    double Q = 19.2;//96./5.
    
    return Q*pow(pow(M_PI,8)*pow(M,5)*pow(f,11),1./3.);  
}
double galactic_binary_Mc(double f0, double dfdt)
{
    double f = f0;
    double fd = dfdt;
    double pi83 = 21.170591578193; //pow(pi,8./3.)

    return pow(fd/(96./5.)/pi83/pow(f,11./3.), 3./5.)/TSUN;
}

double galactic_binary_dL(double f0, double dfdt, double A)
{
    double f    = f0;
    double fd = dfdt;
    double amp   = A;
    return ((5./48.)*(fd/(M_PI*M_PI*f*f*f*amp))*CLIGHT/PC); //seconds  !check notes on 02/28!
}

void galactic_binary_fisher(struct Orbit *orbit, struct Data *data, struct Source *source, struct Noise *noise)
{
    //TODO:  galactic_binary_fisher should compute joint Fisher
    int i,j,n;
    
    int NP = source->NP;
    
    double epsilon    = 1.0e-7;
    double invepsilon2= 1./(2.*epsilon);
    double invstep;
    
    // Plus and minus parameters:
    double *params_p = calloc(NP,sizeof(double));
    double *params_m = calloc(NP,sizeof(double));
    
    // Plus and minus templates for each detector:
    struct Source *wave_p = malloc(sizeof(struct Source));
    struct Source *wave_m = malloc(sizeof(struct Source));
    alloc_source(wave_p, data->N, data->Nchannel, NP);
    alloc_source(wave_m, data->N, data->Nchannel, NP);
    
    // TDI variables to hold derivatives of h
    struct TDI **dhdx = malloc(NP*sizeof(struct TDI *));
    for(n=0; n<NP; n++)
    {
        dhdx[n] = malloc(sizeof(struct TDI));
        alloc_tdi(dhdx[n], data->N, data->Nchannel);
    }
    
    /* assumes all the parameters are log or angle */
    int N2 = data->N*2;
    for(i=0; i<NP; i++)
    {
        //step size for derivatives
        invstep = invepsilon2;
        
        // copy parameters
        for(j=0; j<NP; j++)
        {
            wave_p->params[j] = source->params[j];
            wave_m->params[j] = source->params[j];
        }
        
        // perturb parameters
        wave_p->params[i] += epsilon;
        wave_m->params[i] -= epsilon;
        
        // complete info in source structure
        map_array_to_params(wave_p, wave_p->params, data->T);
        map_array_to_params(wave_m, wave_m->params, data->T);
        
        // clean up TDI arrays, just in case
        for(j=0; j<N2; j++)
        {
            wave_p->tdi->X[j]=0.0;
            wave_p->tdi->A[j]=0.0;
            wave_p->tdi->E[j]=0.0;
            wave_m->tdi->X[j]=0.0;
            wave_m->tdi->A[j]=0.0;
            wave_m->tdi->E[j]=0.0;
        }
        
        // align perturbed waveforms in data array
        galactic_binary_alignment(orbit, data, wave_p);
        galactic_binary_alignment(orbit, data, wave_m);
        
        // compute perturbed waveforms
        galactic_binary(orbit, data->format, data->T, data->t0[0], wave_p->params, NP, wave_p->tdi->X, wave_p->tdi->A, wave_p->tdi->E, wave_p->BW, wave_p->tdi->Nchannel);
        galactic_binary(orbit, data->format, data->T, data->t0[0], wave_m->params, NP, wave_m->tdi->X, wave_m->tdi->A, wave_m->tdi->E, wave_m->BW, wave_m->tdi->Nchannel);
        
        // central differencing derivatives of waveforms w.r.t. parameters
        switch(source->tdi->Nchannel)
        {
            case 1:
                for(n=0; n<N2; n++)
                {
                    dhdx[i]->X[n] = (wave_p->tdi->X[n] - wave_m->tdi->X[n])*invstep;
                }
                break;
            case 2:
                for(n=0; n<N2; n++)
                {
                    dhdx[i]->A[n] = (wave_p->tdi->A[n] - wave_m->tdi->A[n])*invstep;
                    dhdx[i]->E[n] = (wave_p->tdi->E[n] - wave_m->tdi->E[n])*invstep;
                }
                break;
        }
    }
    
    // Calculate fisher matrix
    for(i=0; i<NP; i++)
    {
        for(j=i; j<NP; j++)
        {
            //source->fisher_matrix[i][j] = 10.0; //fisher gets a "DC" level to keep the inversion stable
            switch(source->tdi->Nchannel)
            {
                case 1:
                    source->fisher_matrix[i][j] += fourier_nwip(dhdx[i]->X, dhdx[j]->X, noise->SnX, data->N);
                    break;
                case 2:
                    source->fisher_matrix[i][j] += fourier_nwip(dhdx[i]->A, dhdx[j]->A, noise->SnA, data->N);
                    source->fisher_matrix[i][j] += fourier_nwip(dhdx[i]->E, dhdx[j]->E, noise->SnE, data->N);
                    break;
            }
            if(source->fisher_matrix[i][j]!=source->fisher_matrix[i][j])
            {
                fprintf(stderr,"WARNING: nan matrix element (line %d of file %s)\n",__LINE__,__FILE__);
                fprintf(stderr, "fisher_matrix[%i][%i], Snf=[%g,%g]\n",i,j,noise->SnA[data->N/2],noise->SnE[data->N/2]);
                for(int k=0; k<NP; k++)
                {
                    fprintf(stderr,"source->params[%i]=%g\n",k,source->params[k]);
                }
                source->fisher_matrix[i][j] = 10.0;
            }
            source->fisher_matrix[j][i] = source->fisher_matrix[i][j];
        }
    }
    
    // Calculate eigenvalues and eigenvectors of fisher matrix
    matrix_eigenstuff(source->fisher_matrix, source->fisher_evectr, source->fisher_evalue, NP);
    
    free(params_p);
    free(params_m);
    free_source(wave_p);
    free_source(wave_m);
    
    for(n=0; n<NP; n++) free_tdi(dhdx[n]);
    free(dhdx);
}


int galactic_binary_bandwidth(double L, double fstar, double f, double fdot, double costheta, double A, double T, int N)
{
    int Nmin = 16;
    int Nmax = (int)pow(2,(int)log2((double)(N/2)));
    
    double sqT=sqrt(T);
    
    double sf = sin(f/fstar); //sin(f/f*)
    double sn = AEnoise(L,fstar,f);
    
    //Doppler spreading
    double sintheta = sin(acos(costheta));
    double bw = 8*T*((4.+PI2*f*(AU/CLIGHT)*sintheta)/YEAR + fdot*T);
    int DS = (int)pow(2,(int)log2(bw-1)+1);
    if(DS > Nmax) DS = Nmax;
    if(DS < Nmin) DS = Nmin;
    
    
    //Sinc spreading
    double SNRm = analytic_snr(A, sn, sf, sqT);
    
    int SS = (int)pow(2,(int)log2(SNRm-1)+1);
    
    if(SS > Nmax) SS = Nmax;
    if(SS < Nmin) SS = Nmin;
    
    return (DS > SS) ? DS : SS; //return largest spread as bandwidth
}

void galactic_binary_alignment(struct Orbit *orbit, struct Data *data, struct Source *source)
{
    map_array_to_params(source, source->params, data->T);
    
    source->BW   = 2*galactic_binary_bandwidth(orbit->L, orbit->fstar, source->f0, source->dfdt, source->costheta, source->amp, data->T, data->N);
    source->qmin = (int)(source->f0*data->T) - source->BW/2;
    source->qmax = source->qmin+source->BW;
    source->imin = source->qmin - data->qmin;
    source->imax = source->imin + source->BW;  
}

void galactic_binary(struct Orbit *orbit, char *format, double T, double t0, double *params, int NP, double *X, double *A, double *E, int BW, int NI)
{
    /*   Indicies   */
    int i,j,n;
    /*   Carrier frequency bin  */
    long q;
    /*   Bandwidth      */
    int BW2   = BW*2;
    double invBW2 = 1./(double)BW2;
    
    /*   Gravitational Wave basis vectors   */
    double u[4],v[4],k[4];
    /*   Polarization basis tensors   */
    double eplus[4][4], ecross[4][4];
    /*   Spacecraft position and separation vector   */
    double *x, *y, *z;
    double r12[4]={0},r13[4]={0},r23[4]={0};
    /*   Dot products   */
    double kdotx[4]={0},kdotr[4][4];
    /*   Convenient quantities   */
    double dplus[4][4],dcross[4][4];
    /*   GW source parameters   */
    double phi, psi, amp, Aplus, Across, f0, dfdt, d2fdt2, phi0;
    double costh, sinth, cosph, sinph, cosi, cosps, sinps;
    /*   Time and distance variables   */
    double t, xi[4] = {0};
    /*   Gravitational wave frequency & ratio of f and transfer frequency f*  */
    double f[4] = {0},fonfs[4] = {0};
    /*   LISA response to slow terms (Real & Imaginary pieces)   */
    //Static quantities (Re and Im)
    double DPr, DPi, DCr, DCi;
    //Time varrying quantities (Re & Im) broken up into convenient segments
    double TR[4][4], TI[4][4];
    //Miscellaneous constants used to speed up calculations
    double df;
    /*   Fourier coefficients before FFT and after convolution  */
    //Time series of slowly evolving terms at each vertex
    double *data12, *data13, *data21, *data23, *data31, *data32;
    //Fourier coefficients of slowly evolving terms (numerical)
    double a12[BW2+3], a13[BW2+3], a21[BW2+3], a23[BW2+3], a31[BW2+3], a32[BW2+3];
    //Package cij's into proper form for TDI subroutines
    double ***d;
    
    /*   Allocating Arrays   */
    x = calloc(4,sizeof(double));
    y = calloc(4,sizeof(double));
    z = calloc(4,sizeof(double));
    
    data12 = calloc((BW2+1),sizeof(double));
    data21 = calloc((BW2+1),sizeof(double));
    data31 = calloc((BW2+1),sizeof(double));
    data13 = calloc((BW2+1),sizeof(double));
    data23 = calloc((BW2+1),sizeof(double));
    data32 = calloc((BW2+1),sizeof(double));
    
    d = malloc(sizeof(double**)*4);
    for(i=0; i<4; i++)
    {
        d[i] = malloc(sizeof(double*)*4);
        for(j=0; j<4; j++)
        {
            	d[i][j] = calloc((BW2+1),sizeof(double));
		dplus[i][j] = 0.0;
		dcross[i][j] = 0.0;
		eplus[i][j] = 0.0;
		ecross[i][j] = 0.0;
		kdotr[i][j] = 0.0;
		TR[i][j] = 0.0;
		TI[i][j] = 0.0;
        }
    }
    
    /*   Gravitational Wave source parameters   */
    
    f0     = params[0]/T;
    costh  = params[1];
    phi    = params[2];
    amp    = exp(params[3]);
    cosi   = params[4];
    psi    = params[5];
    phi0   = params[6];
    dfdt   = 0.0;
    d2fdt2 = 0.0;
    if(NP>7)
        dfdt   = params[7]/(T*T);
    if(NP>8)
        d2fdt2 = params[8]/(T*T*T);
    
    //Calculate carrier frequency bin
    q = (long)(f0*T);
    
    //Calculate cos and sin of sky position, inclination, polarization
    sinth	= sqrt(1.0 - costh*costh); //sin(theta) >= 0 (theta -> 0,pi)
    cosph	= cos(phi);
    sinph	= sin(phi);
    cosps	= cos(2.*psi);
    sinps	= sin(2.*psi);
    
    //Calculate GW polarization amplitudes
    Aplus  =  amp*(1.+cosi*cosi);
    Across = -amp*(2.0*cosi);
    
    df = PI2*(((double)q)/T);
    
    //Calculate constant pieces of transfer functions
    DPr =  Aplus*cosps;
    DPi = -Across*sinps;
    DCr = -Aplus*sinps;
    DCi = -Across*cosps;
    
    /*   Tensor construction for buildingslowly evolving LISA response   */
    //Gravitational Wave source basis vectors
    u[1] =  costh*cosph;  u[2] =  costh*sinph;  u[3] = -sinth;
    v[1] =  sinph;        v[2] = -cosph;        v[3] =  0.;
    k[1] = -sinth*cosph;  k[2] = -sinth*sinph;  k[3] = -costh;
    
    //GW polarization basis tensors
    /*
     * LDC convention:
     * https://gitlab.in2p3.fr/LISA/LDC/-/blob/develop/ldc/waveform/fastGB/GB.cc#L530
     */
    for(i=1;i<=3;i++)
    {
        for(j=1;j<=3;j++)
        {
            eplus[i][j]  = v[i]*v[j] - u[i]*u[j];
            ecross[i][j] = u[i]*v[j] + v[i]*u[j];
        }
    }
    
    
    /* Main loop over signal bandwidth */
    for(n=1; n<=BW; n++)
    {
        //First time sample must be at t=0 for phasing
        t = t0 + T*(double)(n-1)/(double)BW;
        
        //Calculate position of each spacecraft at time t
        (*orbit->orbit_function)(orbit, t, x, y, z);
        
        for(i=1; i<=3; i++)
        {
            kdotx[i] = (x[i]*k[1]+y[i]*k[2]+z[i]*k[3])/CLIGHT;
            
            //Wave arrival time at spacecraft i
            xi[i] = t - kdotx[i];
            
            //Zeroeth order approximation to frequency at spacecraft i
            f[i] = f0;
            
            //First order in frequency
            if(NP>7) f[i] += dfdt*xi[i];
            
            //Second order in frequency
            if(NP>8) f[i] += 0.5*d2fdt2*xi[i]*xi[i];
            
            //Ratio of true frequency to transfer frequency
            fonfs[i] = f[i]/orbit->fstar;
        }
        
        
        //Unit separation vector from spacecrafts i to j
        r12[1] = (x[2] - x[1])/orbit->L;   r13[1] = (x[3] - x[1])/orbit->L;   r23[1] = (x[3] - x[2])/orbit->L;
        r12[2] = (y[2] - y[1])/orbit->L;   r13[2] = (y[3] - y[1])/orbit->L;   r23[2] = (y[3] - y[2])/orbit->L;
        r12[3] = (z[2] - z[1])/orbit->L;   r13[3] = (z[3] - z[1])/orbit->L;   r23[3] = (z[3] - z[2])/orbit->L;
        
        
        //Zero arrays to be summed
        dplus[1][2]  = dplus[1][3]  = dplus[2][1]  = dplus[2][3]  = dplus[3][1]  = dplus[3][2]  = 0.;
        dcross[1][2] = dcross[1][3] = dcross[2][1] = dcross[2][3] = dcross[3][1] = dcross[3][2] = 0.;
        
        //Convenient quantities d+ & dx
        for(i=1; i<=3; i++)
        {
            for(j=1; j<=3; j++)
            {
                dplus[1][2]  += r12[i]*r12[j]*eplus[i][j];   dcross[1][2] += r12[i]*r12[j]*ecross[i][j];
                dplus[2][3]  += r23[i]*r23[j]*eplus[i][j];   dcross[2][3] += r23[i]*r23[j]*ecross[i][j];
                dplus[1][3]  += r13[i]*r13[j]*eplus[i][j];   dcross[1][3] += r13[i]*r13[j]*ecross[i][j];
            }
        }
        
        //Make use of symmetry
        dplus[2][1] = dplus[1][2];  dcross[2][1] = dcross[1][2];
        dplus[3][2] = dplus[2][3];  dcross[3][2] = dcross[2][3];
        dplus[3][1] = dplus[1][3];  dcross[3][1] = dcross[1][3];
        
        //Zero arrays to be summed
        kdotr[1][2] = kdotr[1][3] = kdotr[2][1] = kdotr[2][3] = kdotr[3][1] = kdotr[3][2] = 0.;
        for(i=1; i<=3; i++)
        {
            kdotr[1][2] += k[i]*r12[i];   kdotr[1][3] += k[i]*r13[i];   kdotr[2][3] += k[i]*r23[i];
        }
        
        //Make use of antisymmetry
        kdotr[2][1] = -kdotr[1][2];
        kdotr[3][1] = -kdotr[1][3];
        kdotr[3][2] = -kdotr[2][3];
        
        //Calculating Transfer function
        for(i=1; i<=3; i++)
        {
            //Argument of complex exponentials
            /*
             * LDC phase parameter in key files is
             * -phi0, hence the -phi0 in arg2
             */
            double arg2 = PI2*f0*xi[i] - phi0 - df*t;
            
            
            //First order frequency evolution
            if(NP>7) arg2 += M_PI*dfdt*xi[i]*xi[i];
            
            //Second order frequency evolution
            if(NP>8) arg2 += (M_PI/3.0)*d2fdt2*xi[i]*xi[i]*xi[i];
            
            //Evolution of amplitude
            double aevol = 1.0;
            
            //First order amplitude evolution
            if(NP>7) aevol += 0.66666666666666666666*dfdt/f0*xi[i];
            
            //Second order amplitude evolution
            //if(NP>8) aevol += const.*d2fdt2*xi[i]*xi[i]/f0;
            
            for(j=1; j<=3; j++)
            {
                if(i!=j)
                {
                    //Argument of transfer function
                    /*
                     * Set to match LDC convention
                     *
                     https://gitlab.in2p3.fr/LISA/LDC/-/blob/develop/ldc/waveform/fastGB/GB.cc#L601
                     */
                    double arg1 = 0.5*fonfs[i]*(1.0 + kdotr[i][j]);
                    
                    //Transfer function
                    double sinc = 0.25*sinf(arg1)/arg1;
                    
                    //Real and imaginary pieces of time series (no complex exponential)
                    double tran1r = aevol*(dplus[i][j]*DPr + dcross[i][j]*DCr);
                    double tran1i = aevol*(dplus[i][j]*DPi + dcross[i][j]*DCi);
                    
                    //Real and imaginry components of complex exponential
                    double tran2r = cos(arg1 + arg2);
                    double tran2i = sin(arg1 + arg2);
                    
                    //Real & Imaginary part of the slowly evolving signal
                    TR[i][j] = sinc*(tran1r*tran2r - tran1i*tran2i);
                    TI[i][j] = sinc*(tran1r*tran2i + tran1i*tran2r);
                }
            }
        }
        
        //Fill  time series data arrays with slowly evolving signal->
        //dataij corresponds to fractional arm length difference yij
        j = 2*n;
        i = j-1;
        data12[i] = TR[1][2];   data21[i] = TR[2][1];   data31[i] = TR[3][1];
        data12[j] = TI[1][2];   data21[j] = TI[2][1];   data31[j] = TI[3][1];
        data13[i] = TR[1][3];   data23[i] = TR[2][3];   data32[i] = TR[3][2];
        data13[j] = TI[1][3];   data23[j] = TI[2][3];   data32[j] = TI[3][2];
    }
    
    /*   Numerical Fourier transform of slowly evolving signal   */
    gsl_fft_complex_radix2_forward (data12+1, 1, BW);
    gsl_fft_complex_radix2_forward (data21+1, 1, BW);
    gsl_fft_complex_radix2_forward (data31+1, 1, BW);
    gsl_fft_complex_radix2_forward (data13+1, 1, BW);
    gsl_fft_complex_radix2_forward (data23+1, 1, BW);
    gsl_fft_complex_radix2_forward (data32+1, 1, BW);
    
    //Unpack arrays from fft and normalize
    for(i=1; i<=BW; i++)
    {
        j = i + BW;
        a12[i] = data12[j]*invBW2;  a21[i] = data21[j]*invBW2;  a31[i] = data31[j]*invBW2;
        a12[j] = data12[i]*invBW2;  a21[j] = data21[i]*invBW2;  a31[j] = data31[i]*invBW2;
        a13[i] = data13[j]*invBW2;  a23[i] = data23[j]*invBW2;  a32[i] = data32[j]*invBW2;
        a13[j] = data13[i]*invBW2;  a23[j] = data23[i]*invBW2;  a32[j] = data32[i]*invBW2;
    }
    
    /*   Renormalize so that the resulting time series is real   */
    for(i=1; i<=BW2; i++)
    {
        d[1][2][i] = a12[i];  d[2][1][i] = a21[i];  d[3][1][i] = a31[i];
        d[1][3][i] = a13[i];  d[2][3][i] = a23[i];  d[3][2][i] = a32[i];
    }
    
    /*   Call subroutines for synthesizing different TDI data channels  */
    if(strcmp("phase",format) == 0)
        LISA_tdi(orbit->L, orbit->fstar, T, d, f0, q, X-1, A-1, E-1, BW, NI);
    else if(strcmp("frequency",format) == 0)
        LISA_tdi_FF(orbit->L, orbit->fstar, T, d, f0, q, X-1, A-1, E-1, BW, NI);
    else
    {
        fprintf(stderr,"Unsupported data format %s",format);
        exit(1);
    }
    
    /*   Free Arrays   */
    free(x);
    free(y);
    free(z);
    
    free(data12);
    free(data21);
    free(data31);
    free(data13);
    free(data23);
    free(data32);
    
    for(i=0; i<4; i++)
    {
        for(j=0; j<4; j++)
        {
            free(d[i][j]);
        }
        free(d[i]);
    }
    free(d);
    
    
    return;
}
