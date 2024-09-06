/*
 *  Copyright (C) 2024 Neil J. Cornish & Tyson B. Littenberg (MSFC-ST12)
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

#include "glass_utils.h"

static double phitilde(double om, double insDOM, double A, double B)
{
    double x, y, z;
    
    z = 0.0;
    
    if(fabs(om) >= A && fabs(om) < A+B)
    {
        x = (fabs(om)-A)/B;
        y = gsl_sf_beta_inc(WAVELET_FILTER_CONSTANT, WAVELET_FILTER_CONSTANT, x);
        z = insDOM*cos(y*M_PI/2.0);
    }
    
    if(fabs(om) < A) z = insDOM;
    
    return(z);
    
}

static void wavelet(struct Wavelets *wdm, int m, double *wave)
{
    
    int N = wdm->N;
    double A = wdm->A;
    double B = wdm->B;
    double insDOM = wdm->inv_root_dOmega;
    double dom = wdm->domega;
    double DOM = wdm->dOmega;
    
    double omega;
    double x, y, z;
    
    double *DE = (double*)malloc(sizeof(double)*(2*N));
    
    // zero and postive frequencies
    for(int i=0; i<=N/2; i++)
    {
        omega = (double)(i)*dom;
        
        y = phitilde(omega+(double)(m)*DOM, insDOM, A, B);
        z = phitilde(omega-(double)(m)*DOM, insDOM, A, B);
        
        x = y+z;
        
        REAL(DE,i) = M_SQRT1_2*x;
        IMAG(DE,i) = 0.0;
        
    }
    
    // negative frequencies
    for(int i=1; i< N/2; i++)
    {
        omega = -(double)(i)*dom;
        
        y = phitilde(omega+(double)(m)*DOM, insDOM, A, B);
        z = phitilde(omega-(double)(m)*DOM, insDOM, A, B);
        
        x = y+z;
        
        REAL(DE,N-i) = M_SQRT1_2*x;
        IMAG(DE,N-i) = 0.0;
        
    }
    
    gsl_fft_complex_radix2_backward(DE, 1, N);
        
    for(int i=0; i < N/2; i++)
    {
        wave[i] = REAL(DE,N/2+i)/wdm->norm;
        wave[i+N/2] = REAL(DE,i)/wdm->norm;
    }
    
    free(DE);
    
}

static void wavelet_window(struct Wavelets *wdm)
{
    double *DX = (double*)malloc(sizeof(double)*(2*wdm->N));
    
    //zero frequency
    REAL(DX,0) =  wdm->inv_root_dOmega;
    IMAG(DX,0) =  0.0;
    
    for(int i=1; i<= wdm->N/2; i++)
    {
        int j = wdm->N-i;
        double omega = (double)(i)*wdm->domega;
        
        // postive frequencies
        REAL(DX,i) = phitilde(omega, wdm->inv_root_dOmega, wdm->A, wdm->B);
        IMAG(DX,i) =  0.0;
        
        // negative frequencies
        REAL(DX,j) =  phitilde(-omega, wdm->inv_root_dOmega, wdm->A, wdm->B);
        IMAG(DX,j) =  0.0;
    }
    
    gsl_fft_complex_radix2_backward(DX, 1, wdm->N);
    
    
    wdm->window = (double*)malloc(sizeof(double)* (wdm->N));
    for(int i=0; i < wdm->N/2; i++)
    {
        wdm->window[i] = REAL(DX,wdm->N/2+i);
        wdm->window[wdm->N/2+i] = REAL(DX,i);
    }
    
    wdm->norm = 0.0;
    for(int i=0; i < wdm->N; i++) wdm->norm += wdm->window[i]*wdm->window[i]*wdm->cadence;
    wdm->norm = sqrt(wdm->norm);

    free(DX);
}

static void wavelet_lookup_table(struct Wavelets *wdm)
{
    double t,f,phase;
    double f0;
    char filename[128];
    FILE *outstream;
    
    // it turns out that all the wavelet layers are the same modulo a
    // shift in the reference frequency. Just have to do a single layer
    // we pick one far from the boundaries to avoid edge effects
    
    double *wave = (double*)malloc(sizeof(double)*(wdm->N));
    int ref_layer = wdm->NF/2;
    wavelet(wdm, ref_layer, wave);
    
    // The odd wavelets coefficienst can be obtained from the even.
    // odd cosine = -even sine, odd sine = even cosine
    // each wavelet covers a frequency band of width DW
    // execept for the first and last wasvelets
    // there is some overlap. The wavelet pixels are of width
    // DOM/PI, except for the first and last which have width
    // half that
    
    f0 = (double)(ref_layer)*wdm->df;
    

    for(int j=0; j<wdm->fdot_steps; j++)  // loop over f-dot slices
    {
        sprintf(filename, "WDMcoeffs%d.bin", j);
        outstream = fopen(filename,"wb");

        int NT = wdm->table[j]->size/2;
        
        
        for(int n=0; n<NT; n++)  // loop of frequency slices
        {
            f = f0 + ((double)(n-NT/2)+0.5)*wdm->deltaf;
            
            double real_coeff = 0.0;
            double imag_coeff = 0.0;
            
            for(int i=0; i<wdm->N; i++)
            {
                t = ((double)(i-wdm->N/2))*wdm->cadence;
                phase = PI2*f*t + M_PI*wdm->fdot[j]*t*t;
                real_coeff += wave[i]*cos(phase)*wdm->cadence;
                imag_coeff += wave[i]*sin(phase)*wdm->cadence;
            }
            gsl_vector_set(wdm->table[j],2*n,real_coeff);
            gsl_vector_set(wdm->table[j],2*n+1,imag_coeff);
        }
        
        gsl_vector_fwrite(outstream, wdm->table[j]);
        fclose(outstream);
    }
    
    free(wave);
}
void initialize_wavelet(struct Wavelets *wdm, double T)
{
    fprintf(stdout,"\n======= Initialize Wavelet Basis =======\n");
    //N = total number of data samples
    //NF = number of frequency layers
    //dt = sampling cadence
    
    wdm->NT = (int)ceil(T/WAVELET_DURATION);
    wdm->NF = WAVELET_DURATION/LISA_CADENCE;
    wdm->cadence = LISA_CADENCE;//7.5;
    wdm->df = WAVELET_BANDWIDTH;
    wdm->dt = WAVELET_DURATION;

    wdm->frequency_steps = 400;
    wdm->fdot_steps = 4;
    wdm->d_fdot = 0.1;
    wdm->oversample = 16.0;

    wdm->N = wdm->oversample * 2 * wdm->NF;
    wdm->T = wdm->N*wdm->cadence;

    wdm->Omega = M_PI/wdm->cadence;
    wdm->dOmega = wdm->Omega/(double)wdm->NF;
    wdm->domega = PI2/wdm->T;
    wdm->inv_root_dOmega = 1.0/sqrt(wdm->dOmega);
    wdm->B = wdm->Omega/(double)(2*wdm->NF);
    wdm->A = (wdm->dOmega-wdm->B)/2.0;
    wdm->BW = (wdm->A+wdm->B)/M_PI;
    wdm->deltaf = wdm->BW/(double)(wdm->frequency_steps);

    wdm->fdot = malloc(wdm->fdot_steps*sizeof(double));

    wdm->table   = malloc(wdm->fdot_steps*sizeof(gsl_vector *));
    wdm->n_table = malloc(wdm->fdot_steps*sizeof(int));

    double fdot_step = wdm->df/wdm->T*wdm->d_fdot; // sets the f-dot increment
            
    for(int n=0; n<wdm->fdot_steps; n++)
    {
        wdm->fdot[n] = (double)(n) * fdot_step;
        
        size_t N = (int)((wdm->BW+wdm->fdot[n]*wdm->T)/wdm->deltaf);
        if(N%2 != 0) N++; // makes sure it is an even number
        wdm->n_table[n] = N;
        wdm->table[n] =  gsl_vector_alloc(2*N);
    }
    
    //stores window function and normalization
    wavelet_window(wdm);
    
    //stores lookup table of wavelet basis functions
    wavelet_lookup_table(wdm);

    //set defaults for min and maximum pixels
    wavelet_pixel_to_index(wdm,0,1,&wdm->kmin);         //first pixel of second layer
    wavelet_pixel_to_index(wdm,0,wdm->NF-1,&wdm->kmax); //first pixel of last layer

    fprintf(stdout,"Number of time pixels:        %i\n", wdm->NT);
    fprintf(stdout,"Duration of time pixels:      %g [hr]\n", wdm->dt/3600);
    fprintf(stdout,"Number of frequency layers:   %i\n", wdm->NF);
    fprintf(stdout,"Bandwidth of frequency layer: %g [uHz]\n", wdm->df*1e6);
    fprintf(stdout,"\n========================================\n");
}

void wavelet_index_to_pixel(struct Wavelets *wdm, int *i, int *j, int k)
{
    int NT = wdm->NT;
    
    //which time
    *i = k%NT; 
    
    //which frequency
    *j = (k - (*i))/NT; //scary integer math
}

void wavelet_pixel_to_index(struct Wavelets *wdm, int i, int j, int *k)
{
    int NT = wdm->NT;
    
    *k = i + j*NT;
}

void wavelet_transform(struct Wavelets *wdm, double *data)
{
    //array index for tf pixel
    int k;
    
    //total data size
    int ND = wdm->NT*wdm->NF;
    
    //windowed data packets
    double *wdata = double_vector(wdm->N);

    //wavelet wavepacket transform of the signal
    double **wave = double_matrix(wdm->NT,wdm->NF);
    
    //normalization factor
    double fac = M_SQRT2*wdm->cadence/wdm->norm;
    
    //do the wavelet transform by convolving data w/ window and iFFT
    for(int i=0; i<wdm->NT; i++)
    {
        
        for(int j=0; j<wdm->N; j++)
        {
            int n = i*wdm->NF - wdm->N/2 + j;
            if(n < 0)   n += ND;  // periodically wrap the data
            if(n >= ND) n -= ND;  // periodically wrap the data
            wdata[j] = data[n] * wdm->window[j];  // apply the window
        }
        
        gsl_fft_real_radix2_transform(wdata, 1, wdm->N);
        
        //unpack Fourier transform
        wave[i][0] = M_SQRT2*fac*wdata[0];
        for(int j=1; j<wdm->NF; j++)
        {
            if((i+j)%2 ==0)
                wave[i][j] = fac*wdata[j*wdm->oversample];
            else
                wave[i][j] = -fac*wdata[wdm->N-j*wdm->oversample];
        }
    }
    
    //replace data vector with wavelet transform mapped from pixel to index
    for(int i=0; i<wdm->NT; i++)
    {
        for(int j=0; j<wdm->NF; j++)
        {
            //get index number k for tf pixel {i,j}
            wavelet_pixel_to_index(wdm,i,j,&k);
            
            //replace data array
            data[k] = wave[i][j];
        }
    }
    
    free_double_vector(wdata);
    free_double_matrix(wave,wdm->NT);
}

void wavelet_transform_from_table(struct Wavelets *wdm, double *phase, double *freq, double *freqd, double *amp, int *jmin, int *jmax, double *wave, int Nmax)
{
    
    int n, k, jj, kk;
    double dx, dy;
    double f, fdot;
    double fmid,fsam;
    double cos_phase, sin_phase, y, z, yy, zz;

    double df = wdm->deltaf;

    // maximum frequency and frequency derivative
    double f_max    = wdm->df*(wdm->NF-1);
    double fdot_max = wdm->fdot[wdm->fdot_steps-1];
    double d_fdot   = wdm->fdot[1]; // f-dot increment
    
    int wave_index = 0;

    for(int i=0; i<wdm->NT; i++)
    {
        f     = freq[i];
        fdot  = freqd[i];
        if(fdot < 0.0) fdot = 0.0;
        
        //skip this step if f or fdot violate bounds
        if(f>=f_max || fdot>=fdot_max)  continue;
        
        cos_phase = amp[i]*cos(phase[i]);
        sin_phase = amp[i]*sin(phase[i]);
        
        n = (int)floor(fdot/d_fdot);  // lower f-dot layer
        dy = fdot/d_fdot - n;         // where in the layer
                                    
        for(int j=jmin[i]; j<=jmax[i]; j++)
        {
        
            // central frequency
            fmid = j*wdm->df;
                
            kk = (int)floor( ( f - (fmid + 0.5*df) )/df );
            fsam = fmid + (kk + 0.5)*df;
            dx = (f - fsam)/df; // used for linear interpolation
                
            // interpolate over frequency
            y = 0.0;
            z = 0.0;
            yy = 0.0;
            zz = 0.0;

            jj = kk + wdm->n_table[n]/2;
            if(jj>=0 && jj< wdm->n_table[n]-1)
            {
                y = (1.0-dx)*gsl_vector_get(wdm->table[n],2*jj)   + dx*gsl_vector_get(wdm->table[n],2*(jj+1));
                z = (1.0-dx)*gsl_vector_get(wdm->table[n],2*jj+1) + dx*gsl_vector_get(wdm->table[n],2*(jj+1)+1);
            }

            jj = kk + wdm->n_table[n+1]/2;
            if(jj >=0 && jj < wdm->n_table[n]-1)
            {
                yy = (1.0-dx)*gsl_vector_get(wdm->table[n+1],2*jj)   + dx*gsl_vector_get(wdm->table[n+1],2*(jj+1));
                zz = (1.0-dx)*gsl_vector_get(wdm->table[n+1],2*jj+1) + dx*gsl_vector_get(wdm->table[n+1],2*(jj+1)+1);
            }
                
            // interpolate over fdot
            y = (1.0-dy)*y + dy*yy;
            z = (1.0-dy)*z + dy*zz;

            // make sure pixel is in range
            wavelet_pixel_to_index(wdm,i,j,&k);
            if(k>=wdm->kmin && k<wdm->kmax)
            {  
                if(wave_index<Nmax)
                {
                    if((i+j)%2 == 0) wave[wave_index] =  (cos_phase*y - sin_phase*z);
                    else             wave[wave_index] = -(cos_phase*z + sin_phase*y);
                    wave_index++;
                }
                else
                {
                    //fprintf(stderr,"Warning, wavelet_transform_from_table tried accessing array out of bounds\n");
                    //fflush(stderr);
                }
            }
            
        } //loop over frequency layers
        
    } //loop over time steps
}

void active_wavelet_list(struct Wavelets *wdm, double *freqX, double *freqY, double *freqZ, double *fdotX, double *fdotY, double *fdotZ, int *wavelet_list, int *Nwavelet, int *jmin, int *jmax)
{
    
    int n;
    int k;
    int N;
    double fdot;
    double fmx, fdmx, dfd, HBW;
    double fx, fy, fz;
    double fmax, fmin;
    
    double df = wdm->deltaf;
    double DF = wdm->df;
    double *fd = wdm->fdot;

    // maximum frequency and frequency derivative
    fmx  = (double)(wdm->NF-1)*DF;
    fdmx = fd[wdm->fdot_steps-1];
    dfd  = fd[1]; // f-dot increment
    
    N = 0;
    for(int i=0; i<wdm->NT; i++)
    {

        // find smallest fdot
        fdot = fdotX[i];
        if(fdotY[i] < fdot) fdot = fdotY[i];
        if(fdotZ[i] < fdot) fdot = fdotZ[i];
        if(fdot < 0.0) fdot = 0.0;
        
        fx = freqX[i];
        fy = freqY[i];
        fz = freqZ[i];
        
        // find the largest and smallest frequencies
        fmin = fx;
        fmax = fx;
        if(fy < fmin) fmin = fy;
        if(fy > fmax) fmax = fy;
        if(fz < fmin) fmin = fz;
        if(fz > fmax) fmax = fz;
        
        if(fmax < fmx && fdot < fdmx)
        {
            // lower f-dot layer
            n = (int)(floor(fdot/dfd));

            // half bandwidth of layer        
            HBW = 0.5*(double)(wdm->n_table[n]-1)*df;
            
            // lowest frequency layer
            jmin[i] = (int)ceil((fmin-HBW)/DF);
            
            // highest frequency layer
            jmax[i] = (int)floor((fmax+HBW)/DF);            
            
            for(int j=jmin[i]; j<=jmax[i]; j++)
            {
                
                wavelet_pixel_to_index(wdm,i,j,&k);
                
                //check that the pixel is in range
                if(k>=wdm->kmin && k<wdm->kmax)
                {
                    wavelet_list[N]=k-wdm->kmin;
                    N++;
                }
            }  // end loop over frequency layers 
        }  
    }
    *Nwavelet = N;
}




