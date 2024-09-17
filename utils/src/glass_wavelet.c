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
    
    double f0 = ref_layer*wdm->df;
    
    #pragma omp parallel for
    for(int j=0; j<wdm->fdot_steps; j++)  // loop over f-dot slices
    {
        //char filename[128];
        //sprintf(filename, "WDMcoeffs%d.bin", j);
        //FILE *outstream = fopen(filename,"wb");

        int NT = wdm->table[j]->size/2;
        
        
        for(int n=0; n<NT; n++)  // loop of frequency slices
        {
            double f = f0 + ((double)(n-NT/2)+0.5)*wdm->deltaf;
            
            double real_coeff = 0.0;
            double imag_coeff = 0.0;
            
            for(int i=0; i<wdm->N; i++)
            {
                double t = ((double)(i-wdm->N/2))*wdm->cadence;
                double phase = PI2*f*t + M_PI*wdm->fdot[j]*t*t;
                real_coeff += wave[i]*cos(phase)*wdm->cadence;
                imag_coeff += wave[i]*sin(phase)*wdm->cadence;
            }
            gsl_vector_set(wdm->table[j],2*n,real_coeff);
            gsl_vector_set(wdm->table[j],2*n+1,imag_coeff);
        }
        
        //gsl_vector_fwrite(outstream, wdm->table[j]);
        //fclose(outstream);
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
    wdm->fdot_steps = 50;
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
        wdm->fdot[n] = -fdot_step*wdm->fdot_steps/2 + n*fdot_step;
        size_t N = (int)((wdm->BW+fabs(wdm->fdot[n])*wdm->T)/wdm->deltaf);
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

void wavelet_transform_inverse(struct Wavelets *wdm, double *data)
{
    int k;
    int N = wdm->NT*wdm->NF;
    double *phit  = double_vector(wdm->NT/2+1);
    double *row   = double_vector(wdm->NT*2);
    double *dataf = double_vector(N); 
    double *datat = double_vector(N);
    double scale  = sqrt(M_PI)/wdm->cadence;
    double sign;
    
    for(int i=0; i<=wdm->NT/2; i++)
    {
        phit[i] = phitilde(i*wdm->domega, wdm->inv_root_dOmega, wdm->A, wdm->B);
    }

    for(int j=1; j<wdm->NF-1; j++)
    {
        if(j%2==0) sign =  1.0;
        else       sign = -1.0;

        for(int i=0; i<wdm->NT; i++)
        {
            wavelet_pixel_to_index(wdm,i,j,&k);
            k-=wdm->kmin;
            
            if((i+j)%2==0)
            {
                REAL(row,i) = data[k];
            }
            else
            {
                if(j%2==0) IMAG(row,i) = -data[k];
                else       IMAG(row,i) =  data[k];
            }
        }

        gsl_fft_complex_radix2_forward (row, 1, wdm->NT);
        int jj = j*(wdm->NT/2);

        // negative frequencies
        for(int i=wdm->NT/2-1; i>0; i--)
        {
            int kk = jj-i;
            dataf[kk]   += sign*scale*phit[i]*REAL(row,wdm->NT-i);
            dataf[N-kk] += sign*scale*phit[i]*IMAG(row,wdm->NT-i);
        }
        
        // positive frequencies
        for(int i=0; i<wdm->NT/2; i++)
        {
            int kk = jj+1;
            dataf[kk]   += sign*scale*phit[i]*REAL(row,i);
            dataf[N-kk] += sign*scale*phit[i]*IMAG(row,i);
        }

    }

    //print Fourier domain data
    FILE *out = fopen("dataf.dat","w");
    for(int n=1; n<N/2; n++) fprintf(out,"%e %e %e\n", n/wdm->T, dataf[n], dataf[N-n]);
    fclose(out);

    //print Time domain data
    memcpy(datat,dataf,N*sizeof(double));
    gsl_fft_halfcomplex_radix2_inverse(datat, 1, N);
    
    out = fopen("datat.dat","w");
    for(int n=0; n<N; n++)
    {
        fprintf(out,"%e %e\n", n*wdm->cadence, datat[n]);
    }
    fclose(out);

    free_double_vector(datat);
    free_double_vector(dataf);
    free_double_vector(row);
    free_double_vector(phit);
}

void wavelet_transform_from_table(struct Wavelets *wdm, double *phase, double *freq, double *freqd, double *amp, int *jmin, int *jmax, double *wave, int *list, int *rlist, int Nmax)
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
    double fdot_min = wdm->fdot[0];
    double d_fdot   = wdm->fdot[1]-wdm->fdot[0]; // f-dot increment
    
    for(int i=0; i<wdm->NT; i++)
    {
        f     = freq[i];
        fdot  = freqd[i];
        
        //skip this step if f or fdot violate bounds
        if(f>=f_max || fdot>=fdot_max || fdot<=fdot_min)  continue;
        
        cos_phase = amp[i]*cos(phase[i]);
        sin_phase = amp[i]*sin(phase[i]);
        
        n = (int)floor((fdot-fdot_min)/d_fdot);  // lower f-dot layer
        dy = (fdot-fdot_min)/d_fdot - n;         // where in the layer
                                    
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
                int n = rlist[k - wdm->kmin];
                if(n<Nmax)
                {
                    if((i+j)%2 == 0) wave[n] =  (cos_phase*y - sin_phase*z);
                    else             wave[n] = -(cos_phase*z + sin_phase*y);
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

void active_wavelet_list(struct Wavelets *wdm, double *freqX, double *freqY, double *freqZ, double *fdotX, double *fdotY, double *fdotZ, int *wavelet_list, int *reverse_list, int *Nwavelet, int *jmin, int *jmax)
{
    
    int n;
    int k;
    int N;
    double fmx, fdmx, fdmn, dfd, HBW;
    double fmax, fmin;
    double fdotmax, fdotmin;
    int Xflag, Yflag, Zflag;
    
    double df = wdm->deltaf;
    double DF = wdm->df;
    double *fd = wdm->fdot;

    // maximum frequency and frequency derivative
    fmx  = (double)(wdm->NF-1)*DF;
    fdmx = fd[wdm->fdot_steps-1];
    fdmn = fd[0];
    dfd  = fd[1]-fd[0]; // f-dot increment
    
    N = 0;
    for(int i=0; i<wdm->NT; i++)
    {
        // check to see if any of the channels are ok
        Xflag = Yflag = Zflag = 0;
        if(freqX[i] < fmx) Xflag = 1;
        if(freqY[i] < fmx) Yflag = 1;
        if(freqZ[i] < fmx) Zflag = 1;

        // shut off any channel that does not have valid fdots
        if(fdotX[i] < fdmn || fdotX[i] > fdmx) Xflag = 0;
        if(fdotY[i] < fdmn || fdotY[i] > fdmx) Yflag = 0;
        if(fdotZ[i] < fdmn || fdotZ[i] > fdmx) Zflag = 0;

        // skip if no channels have valid values
        if(!Xflag && !Yflag && !Zflag) continue;

        /*  find the largest and smallest frequencies and frequency derivatives
        but only using the valid channels */
        fmin = 1;
        fmax = 0;
        fdotmin =  1;
        fdotmax = -1;

        if(Xflag)
        {
            if(freqX[i]>fmax) fmax=freqX[i];
            if(freqX[i]<fmin) fmin=freqX[i];
            if(fdotX[i]>fdotmax) fdotmax=fdotX[i];
            if(fdotX[i]<fdotmin) fdotmin=fdotX[i];
        }
        
        if(Yflag)
        {
            if(freqY[i]>fmax) fmax=freqY[i];
            if(freqY[i]<fmin) fmin=freqY[i];
            if(fdotY[i]>fdotmax) fdotmax=fdotY[i];
            if(fdotY[i]<fdotmin) fdotmin=fdotY[i];
        }
        
        if(Zflag)
        {
            if(freqZ[i]>fmax) fmax=freqZ[i];
            if(freqZ[i]<fmin) fmin=freqZ[i];
            if(fdotZ[i]>fdotmax) fdotmax=fdotZ[i];
            if(fdotZ[i]<fdotmin) fdotmin=fdotZ[i];
        }
       
        //skip if max/min fdot go out of bounds 
        if(fdotmax >= fdmx || fdotmin <= fdmn) continue;

        // lowest f-dot layer
        n = (int)(floor(fdotmin-fdmn/dfd));
        int NL = wdm->n_table[n];

        // highest f-dot layer
        n = (int)(floor(fdotmax-fdmn/dfd));
        int NH = wdm->n_table[n];

        // find which has the largest number of samples
        if(NL > NH) NH = NL;

        // half bandwidth of layer        
        HBW = 0.5*(NH-1)*df;
        
        // lowest frequency layer
        jmin[i] = (int)ceil((fmin-HBW)/DF);
        
        // highest frequency layer
        jmax[i] = (int)floor((fmax+HBW)/DF);   

        // skip any out-of-bounds layers
        if(jmin[i] < 0) jmin[i] = 0;
        if(jmax[i] > wdm->NF-1) jmax[i] = wdm->NF-1;
        
        for(int j=jmin[i]; j<=jmax[i]; j++)
        {
            wavelet_pixel_to_index(wdm,i,j,&k);
            
            //check that the pixel is in range
            if(k>=wdm->kmin && k<wdm->kmax)
            {
                wavelet_list[N]=k-wdm->kmin;
                reverse_list[k-wdm->kmin]=N;
                N++;
            }
        }  
    }
    *Nwavelet = N;
}




