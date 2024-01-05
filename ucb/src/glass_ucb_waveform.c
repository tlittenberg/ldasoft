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

#include <glass_utils.h>

#include "glass_ucb_model.h"
#include "glass_ucb_waveform.h"


double analytic_snr(double A, double Sn, double Sf, double sqT)
{
    return 0.5*A*sqT*Sf/sqrt(Sn); //not exactly what's in paper--calibrated against (h|h)
}

double snr(struct Source *source, struct Noise *noise)
{
    double snr2=0.0;
    switch(source->tdi->Nchannel)
    {
        case 1: //Michelson
            snr2 += fourier_nwip(source->tdi->X,source->tdi->X,noise->invC[0][0],source->tdi->N);
            break;
        case 2: //A&E
            snr2 += fourier_nwip(source->tdi->A,source->tdi->A,noise->invC[0][0],source->tdi->N);
            snr2 += fourier_nwip(source->tdi->E,source->tdi->E,noise->invC[1][1],source->tdi->N);
            break;
        case 3: //XYZ
            snr2 += fourier_nwip(source->tdi->X,source->tdi->X,noise->invC[0][0],source->tdi->N);
            snr2 += fourier_nwip(source->tdi->Y,source->tdi->Y,noise->invC[1][1],source->tdi->N);
            snr2 += fourier_nwip(source->tdi->Z,source->tdi->Z,noise->invC[2][2],source->tdi->N);
            snr2 += fourier_nwip(source->tdi->X,source->tdi->Y,noise->invC[0][1],source->tdi->N)*2.;
            snr2 += fourier_nwip(source->tdi->X,source->tdi->Z,noise->invC[0][2],source->tdi->N)*2.;
            snr2 += fourier_nwip(source->tdi->Y,source->tdi->Z,noise->invC[1][2],source->tdi->N)*2.;
            break;
    }
    
    return(sqrt(snr2));
}

double snr_prior(double SNR)
{
    //SNRPEAK defined in glass_constants.h
    double dfac  = 1.+SNR/(4.*SNRPEAK);
    double dfac5 = ipow(dfac,5);
    return (3.*SNR)/(4.*SNRPEAK*SNRPEAK*dfac5);
}

double waveform_match(struct Source *a, struct Source *b, struct Noise *noise)
{
    int N = a->tdi->N;
    int NFFT = 2*N;
    double match=0;
    
    double *a_A = calloc(NFFT,sizeof(double));
    double *a_E = calloc(NFFT,sizeof(double));
    double *b_A = calloc(NFFT,sizeof(double));
    double *b_E = calloc(NFFT,sizeof(double));
    
    int qmin = a->qmin - a->imin;
    
    
    //Align waveforms into arrays for summing
    for(int i=0; i<a->BW; i++)
    {
        int j = i+a->qmin-qmin;
        // printf("data->qmin=%i,i=%i, imin=%i, qmin=%i, j=%i\n",qmin,i,a->imin,a->qmin,j);
        
        if(j>-1 && j<N)
        {
            int i_re = 2*i;
            int i_im = i_re+1;
            int j_re = 2*j;
            int j_im = j_re+1;
            
            a_A[j_re] = a->tdi->A[i_re];
            a_A[j_im] = a->tdi->A[i_im];
            a_E[j_re] = a->tdi->E[i_re];
            a_E[j_im] = a->tdi->E[i_im];
        }//check that index is in range
    }//loop over waveform bins
    
    //Align waveforms into arrays for summing
    for(int i=0; i<b->BW; i++)
    {
        int j = i+b->qmin-qmin;
        
        if(j>-1 && j<N)
        {
            int i_re = 2*i;
            int i_im = i_re+1;
            int j_re = 2*j;
            int j_im = j_re+1;
            
            b_A[j_re] = b->tdi->A[i_re];
            b_A[j_im] = b->tdi->A[i_im];
            b_E[j_re] = b->tdi->E[i_re];
            b_E[j_im] = b->tdi->E[i_im];
        }//check that index is in range
    }//loop over waveform bins
    
    
    double aa = fourier_nwip(a_A,a_A,noise->invC[0][0],N) + fourier_nwip(a_E,a_E,noise->invC[1][1],N);
    double bb = fourier_nwip(b_A,b_A,noise->invC[0][0],N) + fourier_nwip(b_E,b_E,noise->invC[1][1],N);
    double ab = fourier_nwip(a_A,b_A,noise->invC[0][0],N) + fourier_nwip(a_E,b_E,noise->invC[1][1],N);
    
    match = ab/sqrt(aa*bb);
    
    free(a_A);
    free(a_E);
    free(b_A);
    free(b_E);
    
    return match;
}

double waveform_distance(struct Source *a, struct Source *b, struct Noise *noise)
{
  int N = a->tdi->N;
  int NFFT = 2*N;

  double *a_A = calloc(NFFT,sizeof(double));
  double *a_E = calloc(NFFT,sizeof(double));
  double *b_A = calloc(NFFT,sizeof(double));
  double *b_E = calloc(NFFT,sizeof(double));

  int qmin = a->qmin - a->imin;


  //Align waveforms into arrays for summing
  for(int i=0; i<a->BW; i++)
  {
    int j = i+a->qmin-qmin;
    
    if(j>-1 && j<N)
    {
      int i_re = 2*i;
      int i_im = i_re+1;
      int j_re = 2*j;
      int j_im = j_re+1;
      
      a_A[j_re] = a->tdi->A[i_re];
      a_A[j_im] = a->tdi->A[i_im];
      a_E[j_re] = a->tdi->E[i_re];
      a_E[j_im] = a->tdi->E[i_im];
    }//check that index is in range
  }//loop over waveform bins

  //Align waveforms into arrays for summing
  for(int i=0; i<b->BW; i++)
  {
    int j = i+b->qmin-qmin;
    
    if(j>-1 && j<N)
    {
      int i_re = 2*i;
      int i_im = i_re+1;
      int j_re = 2*j;
      int j_im = j_re+1;
      
      b_A[j_re] = b->tdi->A[i_re];
      b_A[j_im] = b->tdi->A[i_im];
      b_E[j_re] = b->tdi->E[i_re];
      b_E[j_im] = b->tdi->E[i_im];
    }//check that index is in range
  }//loop over waveform bins

  
    double aa = fourier_nwip(a_A,a_A,noise->invC[0][0],N) + fourier_nwip(a_E,a_E,noise->invC[1][1],N);
    double bb = fourier_nwip(b_A,b_A,noise->invC[0][0],N) + fourier_nwip(b_E,b_E,noise->invC[1][1],N);
    double ab = fourier_nwip(a_A,b_A,noise->invC[0][0],N) + fourier_nwip(a_E,b_E,noise->invC[1][1],N);

  double distance = (aa + bb - 2*ab)/4.;

  free(a_A);
  free(a_E);
  free(b_A);
  free(b_E);
  
  return distance;
}
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
        
    double epsilon    = 1.0e-7;
    //double invepsilon2= 1./(2.*epsilon);
    double invepsilon2= 1./(epsilon);
    double invstep;
    
    // Plus and minus parameters:
    double *params_p = calloc(UCB_MODEL_NP,sizeof(double));
    //double *params_m = calloc(NP,sizeof(double));
    
    // Plus and minus templates for each detector:
    struct Source *wave_p = malloc(sizeof(struct Source));
    //struct Source *wave_m = malloc(sizeof(struct Source));
    alloc_source(wave_p, data->N, data->Nchannel);
    //alloc_source(wave_m, data->N, data->Nchannel, NP);
    
    // TDI variables to hold derivatives of h
    struct TDI **dhdx = malloc(UCB_MODEL_NP*sizeof(struct TDI *));
    for(n=0; n<UCB_MODEL_NP; n++)
    {
        dhdx[n] = malloc(sizeof(struct TDI));
        alloc_tdi(dhdx[n], data->N, data->Nchannel);
    }
    
    /* assumes all the parameters are log or angle */
    int N2 = data->N*2;
    for(i=0; i<UCB_MODEL_NP; i++)
    {
        //step size for derivatives
        invstep = invepsilon2;
        
        // copy parameters
        for(j=0; j<UCB_MODEL_NP; j++)
        {
            wave_p->params[j] = source->params[j];
            //wave_m->params[j] = source->params[j];
        }
        
        // perturb parameters
        wave_p->params[i] += epsilon;
        //wave_m->params[i] -= epsilon;
        
	    // catch when cosine parameters get pushed out of bounds
        if(i==1 || i==4)
        {
            if(wave_p->params[i] > 1.0) wave_p->params[i] = 1.0;
            //if(wave_m->params[i] <-1.0) wave_m->params[i] =-1.0;
        }

        // complete info in source structure
        map_array_to_params(wave_p, wave_p->params, data->T);
        //map_array_to_params(wave_m, wave_m->params, data->T);
        
        // clean up TDI arrays, just in case
        for(j=0; j<N2; j++)
        {
            wave_p->tdi->X[j]=0.0;
            wave_p->tdi->Y[j]=0.0;
            wave_p->tdi->Z[j]=0.0;
            wave_p->tdi->A[j]=0.0;
            wave_p->tdi->E[j]=0.0;
        }
        
        // align perturbed waveforms in data array
        galactic_binary_alignment(orbit, data, wave_p);
        //galactic_binary_alignment(orbit, data, wave_m);
        
        // compute perturbed waveforms
        galactic_binary(orbit, data->format, data->T, data->t0, wave_p->params, UCB_MODEL_NP, wave_p->tdi->X, wave_p->tdi->Y, wave_p->tdi->Z, wave_p->tdi->A, wave_p->tdi->E, wave_p->BW, wave_p->tdi->Nchannel);
        
        // central differencing derivatives of waveforms w.r.t. parameters
        switch(source->tdi->Nchannel)
        {
            case 1:
                for(n=0; n<wave_p->BW*2; n++)
                {
                    dhdx[i]->X[n] = (wave_p->tdi->X[n] - source->tdi->X[n])*invstep;
                }
                break;
            case 2:
                for(n=0; n<wave_p->BW*2; n++)
                {
                    dhdx[i]->A[n] = (wave_p->tdi->A[n] - source->tdi->A[n])*invstep;
                    dhdx[i]->E[n] = (wave_p->tdi->E[n] - source->tdi->E[n])*invstep;
                }
                break;
            case 3:
                for(n=0; n<wave_p->BW*2; n++)
                {
                    /* Use AE channels for Fisher matrix cuz it's just a proposal*/
                    dhdx[i]->A[n] = (wave_p->tdi->A[n] - source->tdi->A[n])*invstep;
                    dhdx[i]->E[n] = (wave_p->tdi->E[n] - source->tdi->E[n])*invstep;
                     
                    
                    /*
                    dhdx[i]->X[n] = (wave_p->tdi->X[n] - source->tdi->X[n])*invstep;
                    dhdx[i]->Y[n] = (wave_p->tdi->Y[n] - source->tdi->Y[n])*invstep;
                    dhdx[i]->Z[n] = (wave_p->tdi->Z[n] - source->tdi->Z[n])*invstep;*/
                }
                break;

        }
    }
    
    // Calculate fisher matrix
    for(i=0; i<UCB_MODEL_NP; i++)
    {
        for(j=i; j<UCB_MODEL_NP; j++)
        {
            //source->fisher_matrix[i][j] = 10.0; //fisher gets a "DC" level to keep the inversion stable
            switch(source->tdi->Nchannel)
            {
                case 1:
                    source->fisher_matrix[i][j] = fourier_nwip(dhdx[i]->X, dhdx[j]->X, noise->invC[0][0], wave_p->BW);
                    break;
                case 2:
                    source->fisher_matrix[i][j] = fourier_nwip(dhdx[i]->A, dhdx[j]->A, noise->invC[0][0], wave_p->BW);
                    source->fisher_matrix[i][j] += fourier_nwip(dhdx[i]->E, dhdx[j]->E, noise->invC[1][1], wave_p->BW);
                    break;
                case 3:
                    source->fisher_matrix[i][j] = fourier_nwip(dhdx[i]->A, dhdx[j]->A, noise->invC[0][0], wave_p->BW);
                    source->fisher_matrix[i][j] += fourier_nwip(dhdx[i]->E, dhdx[j]->E, noise->invC[1][1], wave_p->BW);
                    break;

                    /*
                    source->fisher_matrix[i][j] = fourier_nwip(dhdx[i]->X, dhdx[j]->X, noise->invC[0][0], wave_p->BW);
                    source->fisher_matrix[i][j] += fourier_nwip(dhdx[i]->Y, dhdx[j]->Y, noise->invC[1][1], wave_p->BW);
                    source->fisher_matrix[i][j] += fourier_nwip(dhdx[i]->Z, dhdx[j]->Z, noise->invC[2][2], wave_p->BW);

                    source->fisher_matrix[i][j] += fourier_nwip(dhdx[i]->X, dhdx[j]->Y, noise->invC[0][1], wave_p->BW)*2.;
                    source->fisher_matrix[i][j] += fourier_nwip(dhdx[i]->X, dhdx[j]->Z, noise->invC[0][2], wave_p->BW)*2.;
                    source->fisher_matrix[i][j] += fourier_nwip(dhdx[i]->Y, dhdx[j]->Z, noise->invC[1][2], wave_p->BW)*2.;
                    break;*/
            }
            if(source->fisher_matrix[i][j]!=source->fisher_matrix[i][j])
            {
                fprintf(stderr,"WARNING: nan matrix element (line %d of file %s)\n",__LINE__,__FILE__);
                fprintf(stderr, "fisher_matrix[%i][%i], Snf=[%g,%g]\n",i,j,noise->C[0][0][data->N/2],noise->C[1][1][data->N/2]);
                for(int k=0; k<UCB_MODEL_NP; k++)
                {
                    fprintf(stderr,"source->params[%i]=%g\n",k,source->params[k]);
                }
                source->fisher_matrix[i][j] = 10.0;
            }
            source->fisher_matrix[j][i] = source->fisher_matrix[i][j];
        }
    }
    
    // Calculate eigenvalues and eigenvectors of fisher matrix
    matrix_eigenstuff(source->fisher_matrix, source->fisher_evectr, source->fisher_evalue, UCB_MODEL_NP);
    
    free(params_p);
    //free(params_m);
    free_source(wave_p);
    //free_source(wave_m);
    
    for(n=0; n<UCB_MODEL_NP; n++) free_tdi(dhdx[n]);
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
    double bw = 2*T*((4.+PI2*f*(AU/CLIGHT)*sintheta)/YEAR + fabs(fdot)*T);
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

void galactic_binary(struct Orbit *orbit, char *format, double T, double t0, double *params, int NParams, double *X, double *Y, double *Z, double *A, double *E, int BW, int NI)
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
    if(NParams>7)
        dfdt   = params[7]/(T*T);
    if(NParams>8)
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
    
    /*   Tensor construction for building slowly evolving LISA response   */
    //Gravitational Wave source basis vectors
    u[1] =  costh*cosph;  u[2] =  costh*sinph;  u[3] = -sinth;
    v[1] =  sinph;        v[2] = -cosph;        v[3] =  0.;
    k[1] = -sinth*cosph;  k[2] = -sinth*sinph;  k[3] = -costh;
    
    //GW polarization basis tensors
    /*
     * LDC convention:
     * https://gitlab.in2p3.fr/LISA/LDC/-/blob/develop/ldc/waveform/fastGB/GB.cc
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
            if(NParams>7) f[i] += dfdt*xi[i];
            
            //Second order in frequency
            if(NParams>8) f[i] += 0.5*d2fdt2*xi[i]*xi[i];
            
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
            if(NParams>7) arg2 += M_PI*dfdt*xi[i]*xi[i];
            
            //Second order frequency evolution
            if(NParams>8) arg2 += (M_PI/3.0)*d2fdt2*xi[i]*xi[i]*xi[i];
            
            //Evolution of amplitude
            double aevol = 1.0;
            
            //First order amplitude evolution
            if(NParams>7) aevol += 0.66666666666666666666*dfdt/f0*xi[i];
            
            //Second order amplitude evolution
            //if(NParams>8) aevol += const.*d2fdt2*xi[i]*xi[i]/f0;
            
            for(j=1; j<=3; j++)
            {
                if(i!=j)
                {
                    //Argument of transfer function
                    /*
                     * Set to match Radler LDC convention
                     *
                     https://gitlab.in2p3.fr/LISA/LDC/-/blob/develop/ldc/waveform/fastGB/GB.cc
                     */
                    double arg1 = 0.5*fonfs[i]*(1.0 + kdotr[i][j]);
                    
                    //Transfer function
                    double sinc = 0.25*sin(arg1)/arg1;
                    
                    //Real and imaginary pieces of time series (no complex exponential)
                    double tran1r = aevol*(dplus[i][j]*DPr + dcross[i][j]*DCr);
                    double tran1i = aevol*(dplus[i][j]*DPi + dcross[i][j]*DCi);
                    
                    /*
                     * Set to match Sangria LDC convention
                     * which defines the GW as e(-i Phi)
                    */
                    //Real and imaginry components of complex exponential
                    double tran2r = cos(arg1 - arg2);
                    double tran2i = sin(arg1 - arg2);

                    //Real & Imaginary part of the slowly evolving signal
                    TR[i][j] = sinc*( tran1r*tran2r + tran1i*tran2i);
                    TI[i][j] = sinc*(-tran1r*tran2i + tran1i*tran2r);
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
        LISA_tdi(orbit->L, orbit->fstar, T, d, f0, q, X-1, Y-1, Z-1, A-1, E-1, BW, NI);
    else if(strcmp("frequency",format) == 0)
        LISA_tdi_FF(orbit->L, orbit->fstar, T, d, f0, q, X-1, Y-1, Z-1, A-1, E-1, BW, NI);
    else if(strcmp("sangria",format) == 0)
        LISA_tdi_Sangria(orbit->L, orbit->fstar, T, d, f0, q, X-1, Y-1, Z-1, A-1, E-1, BW, NI);
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













static void RAantenna(struct Orbit *orbit, double *params, double Tobs, int NF, double *TF, double *FF, double *kdotx, double *FpAR, double *FpAI, double *FcAR, double *FcAI, double *FpER, double *FpEI, double *FcER, double *FcEI)
{
    
    /*   Indicies   */
    int i, j, k, n;
    
    /*   Gravitational Wave basis vectors   */
    double u[4],v[4],kv[4];
    
    /*   Polarization basis tensors   */
    double eplus[4][4], ecross[4][4];
    
    /*   Spacecraft position and separation vector   */
    double *x, *y, *z;
    double r12[4]={0},r13[4]={0},r23[4]={0};
    double r10[4]={0},r20[4]={0},r30[4]={0};
    
    double q1, q2, q3, q4;
    
    /*   Convenient quantities   */
    double dplus[4][4],dcross[4][4];

    /*   GW Source data   */
    double phi, psi;
    double costh, sinth, cosph, sinph, cosps, sinps;
    
    /*   Time and distance variables   */
    double t, xa, ya, za;
    
    /*   Miscellaneous  */
    double xi;
    
    double fpx, fcx, fpy, fcy, fpz, fcz;
    
    double TR[4][4], TI[4][4];
    double kdr[4][4];
    double kdg[4];
    
    double fr;
    
    double f0, fdot0, fddot0;
    
    /*   Allocating Arrays   */
    x = calloc(4,sizeof(double));
    y = calloc(4,sizeof(double));
    z = calloc(4,sizeof(double));

    
    /*
     params[0] // f*Tobs
     params[1]   // costh
     params[2] // phi
     params[3]  // log Amp
     params[4]  // cosi
     params[5]  // psi
     params[6]  // phi0
     params[7] // fdot
     params[8]   // fddot
     */
    
    
    f0 = params[0]/Tobs;
    fdot0 = params[7]/(Tobs*Tobs);
    fddot0 = params[8]/(Tobs*Tobs*Tobs);
    
    phi = params[2];   // EclipticLongitude
    
    // conventions are flipped relative to BH code
    psi = -params[5];   // polarization
    
    //Calculate cos and sin of sky position, inclination, polarization
    costh = params[1];   sinth = sqrt(1.0-costh*costh);
    cosph = cos(phi);     sinph = sin(phi);
    cosps = cos(2.*psi);  sinps = sin(2.*psi);
    
    
    /*   Tensor basis  */
    // Note that while this basis differs from the BH code, it amounts to a chain of sign for
    // eplus and ecross, which we correct for with an overall minus sign in the response
    u[1] =  costh*cosph;  u[2] =  costh*sinph;  u[3] = -sinth;
    v[1] =  sinph;        v[2] = -cosph;        v[3] =  0.;
    kv[1] = -sinth*cosph;  kv[2] = -sinth*sinph; kv[3] = -costh;
    
    
    for(i=1;i<=3;i++)
    {
        for(j=1;j<=3;j++)
        {
            eplus[i][j]  = u[i]*u[j] - v[i]*v[j];
            ecross[i][j] = u[i]*v[j] + v[i]*u[j];
        }
    }
    
    /*   Main Loop   */
    for(n=0; n< NF; n++)
    {
        // Barycenter time
        t = TF[n];
        
        //Calculate position of each spacecraft at time t
        (*orbit->orbit_function)(orbit, t, x, y, z);

        // guiding center
        xa = (x[1]+x[2]+x[3])/3.0;
        ya = (y[1]+y[2]+y[3])/3.0;
        za = (z[1]+z[2]+z[3])/3.0;
        
        kdotx[n] = (xa*kv[1]+ya*kv[2]+za*kv[3])/CLIGHT;
        
        // detector time and frequency
        xi  = t - kdotx[n];
        
        // frequency at the guiding center
        FF[n] = f0+fdot0*xi+0.5*fddot0*xi*xi;
        
        fr = FF[n]/(2.0*orbit->fstar);
        
        //Unit separation vector from spacecraft i to j
        r12[1] = (x[2] - x[1])/orbit->L;   r13[1] = (x[3] - x[1])/orbit->L;   r23[1] = (x[3] - x[2])/orbit->L;
        r12[2] = (y[2] - y[1])/orbit->L;   r13[2] = (y[3] - y[1])/orbit->L;   r23[2] = (y[3] - y[2])/orbit->L;
        r12[3] = (z[2] - z[1])/orbit->L;   r13[3] = (z[3] - z[1])/orbit->L;   r23[3] = (z[3] - z[2])/orbit->L;
        
        // These are not unit vectors. Just pulling out the L scaling
        r10[1] = (xa-x[1])/orbit->L;   r10[2] = (ya-y[1])/orbit->L;  r10[3] = (za-z[1])/orbit->L;
        r20[1] = (xa-x[2])/orbit->L;   r20[2] = (ya-y[2])/orbit->L;  r20[3] = (za-z[2])/orbit->L;
        r30[1] = (xa-x[3])/orbit->L;   r30[2] = (ya-y[3])/orbit->L;  r30[3] = (za-z[3])/orbit->L;
        
        kdr[1][2] = 0.0;
        for(k=1; k<=3; k++) kdr[1][2] += kv[k]*r12[k];
        kdr[1][3] = 0.0;
        for(k=1; k<=3; k++) kdr[1][3] += kv[k]*r13[k];
        kdr[2][3] = 0.0;
        for(k=1; k<=3; k++) kdr[2][3] += kv[k]*r23[k];
        
        kdr[2][1] = -kdr[1][2];  kdr[3][1] = -kdr[1][3];  kdr[3][2] = -kdr[2][3];
        
        kdg[1] = 0.0;
        for(k=1; k<=3; k++) kdg[1] += kv[k]*r10[k];
        kdg[2] = 0.0;
        for(k=1; k<=3; k++) kdg[2] += kv[k]*r20[k];
        kdg[3] = 0.0;
        for(k=1; k<=3; k++) kdg[3] += kv[k]*r30[k];
        
        for(i=1; i<=3; i++)
        {
            for(j=1; j<=3; j++)
            {
                if(i != j)
                {
                    q1 = fr*(1.0-kdr[i][j]);
                    q2 = fr*(1.0+kdr[i][j]);
                    q3 = -fr*(3.0+kdr[i][j]-2.0*kdg[i]);
                    q4 = -fr*(1.0+kdr[i][j]-2.0*kdg[i]);
                    q1 = (sin(q1)/q1);
                    q2 = (sin(q2)/q2);
                    TR[i][j] = 0.5*(q1*cos(q3)+q2*cos(q4));   // goes to 1 when f/fstar small
                    TI[i][j] = 0.5*(q1*sin(q3)+q2*sin(q4));   // goes to 0 when f/fstar small
                }
            }
        }
        
        dplus[1][2] = dplus[1][3] = dplus[2][1] = dplus[2][3] = dplus[3][1] = dplus[3][2] = 0.;
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
        
        dplus[2][1] = dplus[1][2];  dcross[2][1] = dcross[1][2];
        dplus[3][2] = dplus[2][3];  dcross[3][2] = dcross[2][3];
        dplus[3][1] = dplus[1][3];  dcross[3][1] = dcross[1][3];
        
        fpx = 0.5*( (dplus[1][2]*cosps+dcross[1][2]*sinps)*TR[1][2] - (dplus[1][3]*cosps+dcross[1][3]*sinps)*TR[1][3] );
        fcx = 0.5*( (-dplus[1][2]*sinps+dcross[1][2]*cosps)*TR[1][2] - (-dplus[1][3]*sinps + dcross[1][3]*cosps)*TR[1][3] );
        
        fpy = 0.5*( (dplus[2][3]*cosps+dcross[2][3]*sinps)*TR[2][3] - (dplus[2][1]*cosps+dcross[2][1]*sinps)*TR[2][1] );
        fcy = 0.5*( (-dplus[2][3]*sinps+dcross[2][3]*cosps)*TR[2][3] - (-dplus[2][1]*sinps + dcross[2][1]*cosps)*TR[2][1] );
        
        fpz = 0.5*( (dplus[3][1]*cosps+dcross[3][1]*sinps)*TR[3][1] - (dplus[3][2]*cosps+dcross[3][2]*sinps)*TR[3][2] );
        fcz = 0.5*( (-dplus[3][1]*sinps+dcross[3][1]*cosps)*TR[3][1] - (-dplus[3][2]*sinps + dcross[3][2]*cosps)*TR[3][2] );
        
        FpAR[n] = (2.0*fpx-fpy-fpz)/3.0;
        FcAR[n] = (2.0*fcx-fcy-fcz)/3.0;
        
        FpER[n] = (fpz-fpy)/SQ3;
        FcER[n] = (fcz-fcy)/SQ3;
        
        fpx = 0.5*( (dplus[1][2]*cosps+dcross[1][2]*sinps)*TI[1][2] - (dplus[1][3]*cosps+dcross[1][3]*sinps)*TI[1][3] );
        fcx = 0.5*( (-dplus[1][2]*sinps+dcross[1][2]*cosps)*TI[1][2] - (-dplus[1][3]*sinps + dcross[1][3]*cosps)*TI[1][3] );
        
        fpy = 0.5*( (dplus[2][3]*cosps+dcross[2][3]*sinps)*TI[2][3] - (dplus[2][1]*cosps+dcross[2][1]*sinps)*TI[2][1] );
        fcy = 0.5*( (-dplus[2][3]*sinps+dcross[2][3]*cosps)*TI[2][3] - (-dplus[2][1]*sinps + dcross[2][1]*cosps)*TI[2][1] );
        
        fpz = 0.5*( (dplus[3][1]*cosps+dcross[3][1]*sinps)*TI[3][1] - (dplus[3][2]*cosps+dcross[3][2]*sinps)*TI[3][2] );
        fcz = 0.5*( (-dplus[3][1]*sinps+dcross[3][1]*cosps)*TI[3][1] - (-dplus[3][2]*sinps + dcross[3][2]*cosps)*TI[3][2] );
        
        FpAI[n] = (2.0*fpx-fpy-fpz)/3.0;
        FcAI[n] = (2.0*fcx-fcy-fcz)/3.0;
        
        FpEI[n] = (fpz-fpy)/SQ3;
        FcEI[n] = (fcz-fcy)/SQ3;
        
        
    }
        
    
    free(x);
    free(y);
    free(z);
    
    return;
}
static void Extrinsic(struct Orbit *orbit, double *params, double Tobs, int NF, double *TF, double *AAmp, double *EAmp, double *APhase, double *EPhase)
{
    
    /*   Indicies   */
    int n;
    
    /*   Time and distance variables   */
    double *kdotx;
    
    /*   Miscellaneous  */
    double xi;
    
    double pa, pe;
    double paold, peold;
    double ja, je;
    
    double *FF;
    double *FpAR, *FcAR, *FpER, *FcER;
    double *FpAI, *FcAI, *FpEI, *FcEI;
    
    double t, Amp, Ampx, Phase, RR, II;
    
    double fonfs;
    
    double cosi;
    double Aplus, Across;
    
    double phi0, f0, fdot0, fddot0;
    
    
    kdotx = (double*)malloc(sizeof(double)* (NF));
    FF = (double*)malloc(sizeof(double)* (NF));
    FpAR = (double*)malloc(sizeof(double)* (NF));
    FcAR = (double*)malloc(sizeof(double)* (NF));
    FpER = (double*)malloc(sizeof(double)* (NF));
    FcER = (double*)malloc(sizeof(double)* (NF));
    FpAI = (double*)malloc(sizeof(double)* (NF));
    FcAI = (double*)malloc(sizeof(double)* (NF));
    FpEI = (double*)malloc(sizeof(double)* (NF));
    FcEI = (double*)malloc(sizeof(double)* (NF));
    
    RAantenna(orbit, params, Tobs, NF, TF, FF, kdotx, FpAR, FpAI, FcAR, FcAI, FpER, FpEI, FcER, FcEI);
    
    for(n=0; n< NF; n++)
    {
        AAmp[n] = 0.0;
        EAmp[n] = 0.0;
        APhase[n] = 0.0;
        EPhase[n] = 0.0;
    }
    
    /*
     params[0] // f*Tobs
     params[1]   // costh
     params[2] // phi
     params[3]  // log Amp
     params[4]  // cosi
     params[5]  // psi
     params[6]  // phi0
     params[7] // fdot
     params[8]   // fddot
     */
    
    // conventions are flipped relative to BH code
    cosi = -params[4];  // cos of inclination
    
    phi0 = params[6];
    
    f0 = params[0]/Tobs;
    fdot0 = params[7]/(Tobs*Tobs);
    fddot0 = params[8]/(Tobs*Tobs*Tobs);
    
    Amp = exp(params[3]);
    
    Aplus = 0.5*(1.+cosi*cosi);
    Across = -cosi;
    
    n = 0;
    RR = FpAR[n]*Aplus - FcAI[n]*Across;
    II = FcAR[n]*Across + FpAI[n]*Aplus;
    paold = atan2(II,RR);
    if(paold < 0.0) paold += PI2;
    RR = FpER[n]*Aplus - FcEI[n]*Across;
    II = FcER[n]*Across + FpEI[n]*Aplus;
    peold = atan2(II,RR);
    if(peold < 0.0) peold += PI2;
    
    ja = 0.0;
    je = 0.0;
    
    //out = fopen("map_fast.dat","w");
    for(n=0; n< NF; n++)
    {
        // Barycenter time
        t = TF[n];
        
        // detector time (guiding center)
        xi  = t - kdotx[n];
        
        fonfs = FF[n]/orbit->fstar;
        
        // including leading order amplitude evolution and TDI + fractional frequency modifiers
        Ampx = Amp*(1.0+0.66666666666666666666*fdot0/f0*xi)*(8.0*fonfs*sin(fonfs));
        
        //  - FF[n]/fstar  Derivation says this is needed in the phase. Doesn't seem to be.
        
        // This is the GW phase at the guiding center. The phase shifts at each spacecraft are
        // accounted for in the transfer function. We remove the carrier signal from the phase
        //so as not to overwhelm the spline interpolation. It is restored in the response
        Phase = PI2*(f0*xi - f0*t +0.5*fdot0*xi*xi + 0.1666666666666667*fddot0*xi*xi*xi) + phi0;
        
        RR = FpAR[n]*Aplus - FcAI[n]*Across;
        II = FcAR[n]*Across + FpAI[n]*Aplus;
        
        
        pa = atan2(II,RR);
        if(pa < 0.0) pa += PI2;
        
        if(pa-paold > 6.0) ja -= PI2;
        if(paold-pa > 6.0) ja += PI2;
        paold = pa;
        
        AAmp[n] = Ampx*sqrt(RR*RR+II*II);
        APhase[n] = Phase+pa+ja;
        
        RR = FpER[n]*Aplus - FcEI[n]*Across;
        II = FcER[n]*Across + FpEI[n]*Aplus;
        
        pe = atan2(II,RR);
        if(pe < 0.0) pe += PI2;
        
        if(pe-peold > 6.0) je -= PI2;
        if(peold-pe > 6.0) je += PI2;
        peold = pe;
        
        EAmp[n] = Ampx*sqrt(RR*RR+II*II);
        EPhase[n] = Phase+pe+je;
                
    }
    
    
    /*   Free Arrays   */
    
    free(kdotx);
    free(FF);
    
    free(FcAR);
    free(FpAR);
    free(FcER);
    free(FpER);
    free(FcAI);
    free(FpAI);
    free(FcEI);
    free(FpEI);
    
    return;
}

static void ResponseWavelet(struct Orbit *orbit, struct Wavelets *wdm, double Tobs, double *params, double *TF, double *AAmp, double *EAmp, double *APhase, double *EPhase, double *Af, double *Ef)
{
    double t;
    int i;
    int Nts;
    double px;
    double f0;
    double *TS;
    double *AAS, *EAS, *APS, *EPS;
    
    f0 = params[0]/Tobs;
    
    // pad out past ends to avoid edge effects with the spline
    Nts = (int)(Tobs/DAY)+11;
    
    TS = malloc(Nts*sizeof(double));
    AAS = malloc(Nts*sizeof(double));
    EAS = malloc(Nts*sizeof(double));
    APS = malloc(Nts*sizeof(double));
    EPS = malloc(Nts*sizeof(double));
    
    for (i=0; i< Nts; i++) TS[i] = (double)(i-5)*DAY;
    
    // The Amplitude is returned with the antenna patterns.
    Extrinsic(orbit, params, Tobs, Nts, TS, AAS, EAS, APS, EPS);
        
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    
    gsl_spline *Aspline = gsl_spline_alloc (gsl_interp_cspline, Nts);
    gsl_spline_init(Aspline, TS, AAS, Nts);
    
    gsl_spline *Espline = gsl_spline_alloc (gsl_interp_cspline, Nts);
    gsl_spline_init(Espline, TS, EAS, Nts);
    
    gsl_spline *PAspline = gsl_spline_alloc (gsl_interp_cspline, Nts);
    gsl_spline_init(PAspline, TS, APS, Nts);
    
    gsl_spline *PEspline = gsl_spline_alloc (gsl_interp_cspline, Nts);
    gsl_spline_init(PEspline, TS, EPS, Nts);
    
    free(TS);
    free(AAS);
    free(EAS);
    free(APS);
    free(EPS);
    
    for (i=0; i< wdm->NT; i++)
    {
        t = TF[i];
        
        // restore the Barycenter carrier phase
        px = PI2*f0*t;
        
        APhase[i] = px + gsl_spline_eval(PAspline, t, acc);
        EPhase[i] = px + gsl_spline_eval(PEspline, t, acc);
        
        AAmp[i] = gsl_spline_eval(Aspline, t, acc);
        EAmp[i] = gsl_spline_eval(Espline, t, acc);
        
        // instaneous frequency of the A channel signal
        Af[i] = f0 + gsl_spline_eval_deriv(PAspline, t, acc)/PI2;
        
        // instaneous frequency of the E channel signal
        Ef[i] = f0 + gsl_spline_eval_deriv(PEspline, t, acc)/PI2;
        
    }
    
    gsl_interp_accel_free(acc);
    gsl_spline_free(Aspline);
    gsl_spline_free(Espline);
    gsl_spline_free(PAspline);
    gsl_spline_free(PEspline);
    
    
    
}

void galactic_binary_wavelet(struct Orbit *orbit, struct Wavelets *wdm, double Tobs, double t0, double *params, int NM, double HBW, int *Alist, int *Elist, double *waveA, double *waveE)
{
    
    int i;
    double DT;
    
    double *AAmp, *EAmp, *APhase, *EPhase;
    double *Af, *Ef;
    double *TF;
    
    int Nt = wdm->NT;
    
    DT = Tobs/(double)(Nt);
    
    // indicates this pixel not used
    for(i=0; i< NM; i++)
    {
        Alist[i] = -1;
        Elist[i] = -1;
        waveA[i] = 0.0;
        waveE[i] = 0.0;
    }
    
    AAmp = (double*)malloc(sizeof(double)* (Nt));
    EAmp = (double*)malloc(sizeof(double)* (Nt));
    APhase = (double*)malloc(sizeof(double)* (Nt));
    EPhase = (double*)malloc(sizeof(double)* (Nt));
    Af = (double*)malloc(sizeof(double)* (Nt));
    Ef = (double*)malloc(sizeof(double)* (Nt));
    TF = (double*)malloc(sizeof(double)* (Nt));
    
    for (i=0; i< Nt; i++)
    {
        TF[i] = t0+((double)(i))*DT;  // center of the pixel
    }
    
    ResponseWavelet(orbit, wdm, Tobs, params, TF, AAmp, EAmp, APhase, EPhase, Af, Ef);
    
    wavemake(wdm, HBW, APhase, Af, AAmp, Alist, waveA);
    wavemake(wdm, HBW, EPhase, Ef, EAmp, Elist, waveE);
    
    free(AAmp);
    free(EAmp);
    free(APhase);
    free(EPhase);
    free(Af);
    free(Ef);
    free(TF);
    
    
}
