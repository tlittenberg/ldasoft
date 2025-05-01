//
//  waveform_benchmark.c
//  
//
//  Created by Tyson Littenberg on 3/18/25.
//

//gcc waveform_benchmark.c -lm -lgsl -lhdf5 -lgslcblas -lglass_utils -lglass_ucb -lopenblas -lomp -I/Users/tyson/opt/include -L/Users/tyson/opt/lib -o waveform_benchmark

#include <stdio.h>
#include <time.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_fft_complex.h>

#include <glass_utils.h>
#include <glass_ucb.h>
#include <glass_noise.h>

static void print_wavelet_pixels(struct Wavelets *wdm, struct TDI *tdi, FILE *fptr)
{
    for(int j=0; j<wdm->NF; j++)
    {
        for(int i=0; i<wdm->NT; i++)
        {
            int k;
            wavelet_pixel_to_index(wdm, i, j, &k);
            fprintf(fptr,"%.12g %.12g ",i*WAVELET_DURATION,j*WAVELET_BANDWIDTH + WAVELET_BANDWIDTH/2);
            fprintf(fptr,"%.12g ",tdi->X[k]);
            fprintf(fptr,"%.12g ",tdi->Y[k]);
            fprintf(fptr,"%.12g ",tdi->Z[k]);
            fprintf(fptr,"\n");
        }
        fprintf(fptr,"\n");
    }

}

int main(void)
{
    clock_t start, end;

    double t0=0.0;
    double Tobs = 31457280.0;
    int N = (int)(Tobs/LISA_CADENCE);

    /*
    Define and set up Orbit structure which contains spacecraft ephemerides
    */
    struct Orbit *orbit = malloc(sizeof(struct Orbit));
    initialize_analytic_orbit(orbit);

    /*
    Define and set up Wavelets structure which contains metadata for wavelet basis
    */
    struct Wavelets *wdm = malloc(sizeof(struct Wavelets));
    initialize_wavelet(wdm, Tobs);

    /*
     Set injection parameters for test signal
    */
    double *params = double_vector(8);
    
    double m1    = 0.6*TSUN;
    double m2    = 0.7*TSUN;
    double Mc    = pow(m1*m2,3.0/5.0)/pow(m1+m2,1.0/5.0);
    double f     = 5.0e-3;
    double fdot  = 96.0*pow(M_PI,8.0/3.0)/5.0*pow(Mc,5.0/3.0)*pow(f,11.0/3.0);
    double DL    = 1.0e3*PC/CLIGHT;
    double Amp   = 4.0*pow(Mc,5.0/3.0)*pow(M_PI*f,2.0/3.0)/DL;
    
    params[0] = f;    // f
    params[1] = 0.2;  // costh
    params[2] = 1.0;  // phi
    params[3] = Amp;  // Amp
    params[4] = -0.3; // cosi
    params[5] = 0.8;  // psi
    params[6] = 1.2;  // phi0
    params[7] = fdot; // fdot

    
    /*
    Generate wavelet domain waveform with fast LISA response
    */

    // get spacecraft ephemerides on the spline time grid
    initialize_interpolated_analytic_orbits(orbit, Tobs, t0);
    
    int Nwavelet;                          // number of non-zero wavelet pixels that hold the signal
    int *wavelet_list = int_vector(N);     // list of non-zero wavelet pixels that hold the signal
    struct TDI *tdi = malloc(sizeof(struct TDI));
    alloc_tdi(tdi,N,3);

    params[0] = params[0]*Tobs;
    params[3] = log(params[3]);
    params[7] = params[7]*(Tobs*Tobs);


    ucb_waveform_wavelet_het(orbit, wdm, Tobs, 0.0, params, wavelet_list, &Nwavelet, tdi->X, tdi->Y, tdi->Z);
    FILE *out = fopen("wavelet_het.dat","w");
    print_wavelet_pixels(wdm, tdi, out);
    fclose(out);

    
    
    start = clock();
    int Nwaveforms = 1000;
    for(int mc=0; mc<Nwaveforms; mc++)
    {
        ucb_waveform_wavelet_het(orbit, wdm, Tobs, 0.0, params, wavelet_list, &Nwavelet, tdi->X, tdi->Y, tdi->Z);
    }
    end = clock();
    double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("het wavelet calculation took %f seconds\n", cpu_time_used/(double)Nwaveforms);

    free_tdi(tdi);
    tdi = malloc(sizeof(struct TDI));
    alloc_tdi(tdi,N,3);

    ucb_waveform_wavelet(orbit, wdm, Tobs, 0.0, params, wavelet_list, &Nwavelet, tdi->X, tdi->Y, tdi->Z);
    out = fopen("wavelet_tab.dat","w");
    print_wavelet_pixels(wdm, tdi, out);
    fclose(out);

    start = clock();
    for(int mc=0; mc<Nwaveforms; mc++)
    {
        ucb_waveform_wavelet(orbit, wdm, Tobs, 0.0, params, wavelet_list, &Nwavelet, tdi->X, tdi->Y, tdi->Z);
    }
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("tab wavelet calculation took %f seconds\n", cpu_time_used/(double)Nwaveforms);

    
    return 0;
}
