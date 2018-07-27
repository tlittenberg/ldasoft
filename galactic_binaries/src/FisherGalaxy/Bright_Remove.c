/*********************************************************/
/*                                                       */
/*        Bright_Remove.c, Version 2.3, 4/28/2011        */                                                            
/*      Written by Neil Cornish & Tyson Littenberg       */                                          
/*                                                       */
/* gcc -O2 -o Bright_Remove Bright_Remove.c arrays.c -lm */
/*                                                       */
/*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "arrays.h"
#include "Constants.h"
#include "Detector.h"

struct lisa_orbit{                     
	int N;
	
	double L;
	double fstar;
	
	double *t;              
	double **x;             
	double **y;                     
	double **z;                  
	double **dx;
	double **dy;
	double **dz;
};
struct lisa_orbit orbit;

//void spacecraft(double t, double *x, double *y, double *z);
void spacecraft(struct lisa_orbit *orbit, double tint, double *xint, double *yint, double *zint);
void initialize_orbit(char OrbitFile[], struct lisa_orbit *orbit);
void LISA_spline(double *x, double *y, int n, double yp1, double ypn, double *y2);
void LISA_splint(double *xa, double *ya, double *y2a, int n, double x, double *y);

double quickselect(double *arr, int n, int k);
void medianX(long imin, long imax, double fstar, double L, double *XP, double *Xnoise, double *Xconf);
void medianAE(long imin, long imax, double fstar, double L, double *AEP, double *AEnoise, double *AEconf);
void convolve(long N, double *a, long M, double *b, double *cn);
void XYZ(double L, double fstar, double ***d, double f0, long q, long M, double *XLS, double *ALS, double *ELS);
void FAST_LISA(struct lisa_orbit *orbit, double *params, long N, long M, double *XLS, double *ALS, double *ELS);
void dfour1(double data[], unsigned long nn, int isign);
void instrument_noise(double f, double fstar, double L, double *SAE, double *SXYZ);
void instrument_noise2(double f, double fstar, double L, double *SAE, double *SXYZ);
double ran2(long *idum);
double gasdev2(long *idum);
void KILL(char*);

int main(int argc,char **argv)
{
	
	double f, fdot, theta, phi, A, iota, psi, phase;
	char Gfile[50];
	double *params;
	double *XfLS, *AALS, *EELS;
	double *XLS, *AA, *EE;
	double fonfs, Sn, Sm, Acut;
	long M, N, q;
	long i, j, k, cnt, cc1, mult, imax, imin;
	long rseed;
	double SAE, SXYZ, sqT;
	double XR, XI, AR, AI, ER, EI;
	double *Xnoise, *Xconf;
	double *AEnoise, *AEconf;
	double *XP, *AEP;
	double SNX, SNAE, SNRAE, SNRX;
	double SNRthres;
	
	FILE* Infile;
	FILE* Outfile;
	FILE* Xbright;
	FILE* Abright;
	
	if(argc !=5) KILL("Bright_Remove XAE.dat Noise.dat Bright.dat Orbit.dat\n");
	
	Xbright = fopen("BrightX.dat","w");
	Abright = fopen("BrightAE.dat","w");
	
	//Data structure for interpolating orbits from file	
	struct lisa_orbit *LISAorbit;
	LISAorbit = &orbit;
	
	//Set up orbit structure (allocate memory, read file, cubic spline)
	sprintf(Gfile,"%s",argv[4]);	
	initialize_orbit(Gfile, LISAorbit);
	double L     = LISAorbit->L;
	double fstar = LISAorbit->fstar;
	
	params = dvector(0,9);
	
	SNRthres = 7.0;
	
	if((T/year) <= 8.0) mult = 8;
	if((T/year) <= 4.0) mult = 4;
	if((T/year) <= 2.0) mult = 2;
	if((T/year) <= 1.0) mult = 1;
	
	XfLS = dvector(0,NFFT-1);  AALS = dvector(0,NFFT-1);  EELS = dvector(0,NFFT-1);
	
	for(i=0; i<NFFT; i++)
    {
		XfLS[i] = 0.0;
		AALS[i] = 0.0;
		EELS[i] = 0.0;
    }
	
	imax = (long)ceil(1.0e-2*T);
	imin = (long)floor(1.0e-4*T);
	sqT = sqrt(T);
	
	Infile = fopen(argv[1],"r");
	for(i=1; i< imax; i++)
	{
		fscanf(Infile,"%lf%lf%lf%lf%lf%lf%lf\n", &f, &XfLS[2*i], &XfLS[2*i+1],
			   &AALS[2*i], &AALS[2*i+1], &EELS[2*i], &EELS[2*i+1]);
	}
	fclose(Infile);
	
	XP = dvector(imin,imax);  AEP = dvector(imin,imax);
	Xnoise = dvector(imin,imax);  Xconf = dvector(imin,imax);
	AEnoise = dvector(imin,imax);  AEconf = dvector(imin,imax);
	
	Outfile = fopen("Power_0.dat","w");
	for(i=imin; i< imax; i++)
	{
		f = (double)(i)/T;
		if(noiseFlag==1)instrument_noise(f, fstar, L ,&SAE, &SXYZ);
		if(noiseFlag==2)instrument_noise2(f, fstar, L ,&SAE, &SXYZ);
		XP[i] = (2.0*(XfLS[2*i]*XfLS[2*i] + XfLS[2*i+1]*XfLS[2*i+1]));
		AEP[i] = (2.0*(AALS[2*i]*AALS[2*i]+AALS[2*i+1]*AALS[2*i+1])); 
		fprintf(Outfile,"%e %e %e %e %e\n", f, XP[i], SXYZ, AEP[i], SAE);
	}
	fclose(Outfile);
	
	Outfile = fopen(argv[2],"r");
	for(i=imin; i<= imax; i++)
	{
		fscanf(Outfile,"%lf%lf%lf%lf%lf\n", &f, &Xnoise[i], &Xconf[i], &AEnoise[i], &AEconf[i]);
	}
	fclose(Outfile);
	
	Infile = fopen(argv[3],"r");
	
	printf("Starting Removal\n");
	
	cnt = 0;
	cc1 = 0;
	
    while ( !feof(Infile) )
	{
		fscanf(Infile, "%lf%lf%lf%lf%lf%lf%lf%lf\n", &f, &fdot, &theta, &phi, &A, &iota, &psi, &phase);
		
		params[0] = f;
		params[1] = 0.5*pi-theta;
		params[2] = phi; 
		params[3] = A; 
		params[4] = iota;
		params[5] = psi;
		params[6] = phase;
		params[7] = fdot;
		params[8] = 11.0/3.0*fdot*fdot/f;
		
		if(f > 1.0e-2) cc1++;
		
		if(f < 1.0e-2)
		{
			
			N = 32*mult;
			if(f > 0.001) N = 64*mult;
			if(f > 0.01) N = 256*mult;
			if(f > 0.03) N = 512*mult;
			if(f > 0.1) N = 1024*mult;
			
			
			fonfs = f/fstar;
			
			q = (long)(f*T);
			
			if(noiseFlag==1)instrument_noise(f, fstar, L, &SAE, &SXYZ);
			if(noiseFlag==2)instrument_noise2(f, fstar, L, &SAE, &SXYZ);
			
			/*  calculate michelson noise  */
						
			/*inst2*/
			if(noiseFlag==1) Sm = SXYZ/(4.0*sin(f/fstar)*sin(f/fstar));
			if(noiseFlag==2) Sm = SXYZ/(4.0);//*sin(f/fstar)*sin(f/fstar));
			
			Acut = A*sqrt(T/Sm);
			
			M = (long)(pow(2.0,(rint(log(Acut)/log(2.0))+1.0)));
			
			if(M < N) M = N;
			if(N < M) N = M;
			if(M > 8192) M = 8192;
			
			N = M;
			//N = M = 2*8192;
			
			XLS = dvector(1,2*M);  
			AA = dvector(1,2*M);   EE = dvector(1,2*M);
			
			FAST_LISA(LISAorbit, params, N, M, XLS, AA, EE);
			
			/*inst2*/ 
			if(noiseFlag==1)
			{
				SNX = (SXYZ+Xconf[q]);//*sin(f/fstar)*sin(f/fstar);
				SNAE = (SAE+AEconf[q]);//*sin(f/fstar)*sin(f/fstar);
			}
			if(noiseFlag==2)
			{
				SNX = (SXYZ+Xconf[q])*sin(f/fstar)*sin(f/fstar);
				SNAE = (SAE+AEconf[q])*sin(f/fstar)*sin(f/fstar);
			}
			SNRAE = 0.0;
			SNRX = 0.0;
			for(i=1; i<=M; i++)
			{
				SNRX += 4.0*(XLS[2*i-1]*XLS[2*i-1]+XLS[2*i]*XLS[2*i]);
				SNRAE += 4.0*(AA[2*i-1]*AA[2*i-1]+AA[2*i]*AA[2*i]+EE[2*i-1]*EE[2*i-1]+EE[2*i]*EE[2*i]);
			}
			SNRAE *= T/SNAE;
			SNRX *= T/SNX;
			SNRAE = sqrt(SNRAE);
			SNRX = sqrt(SNRX);
			
			if(SNRX > SNRthres)
			{
				fprintf(Xbright, "%.16f %.10e %f %f %e %f %f %f\n", f, fdot, theta, phi, A, iota, psi, phase);
				for(i=1; i<=M; i++)
				{
					k = (q + i -1 - M/2);
					if(k>0)
					{
					if(noiseFlag==1)
					{
					XfLS[2*k] -= sqT*XLS[2*i-1];
					XfLS[2*k+1] -= sqT*XLS[2*i];
					}
					if(noiseFlag==2)
					{
						XfLS[2*k] -= sqT*XLS[2*i-1]/sin(f/fstar);
						XfLS[2*k+1] -= sqT*XLS[2*i]/sin(f/fstar);
					}
					}
				}
			}
			
			if(SNRAE > SNRthres)
			{
				fprintf(Abright, "%.16f %.10e %f %f %e %f %f %f\n", f, fdot, theta, phi, A, iota, psi, phase);
				for(i=1; i<=M; i++)
				{
					k = (q + i -1 - M/2);
					if(k>0)
					{
					if(noiseFlag==1)
					{
					AALS[2*k] -= sqT*AA[2*i-1];
					AALS[2*k+1] -= sqT*AA[2*i];
					EELS[2*k] -= sqT*EE[2*i-1];
					EELS[2*k+1] -= sqT*EE[2*i];
					}
					if(noiseFlag==2)
					{
						AALS[2*k] -= sqT*AA[2*i-1]/sin(f/fstar);
						AALS[2*k+1] -= sqT*AA[2*i]/sin(f/fstar);
						EELS[2*k] -= sqT*EE[2*i-1]/sin(f/fstar);
						EELS[2*k+1] -= sqT*EE[2*i]/sin(f/fstar);
					}
					}

				}
			}
			
			
			free_dvector(XLS,1,2*M);  free_dvector(AA,1,2*M);  free_dvector(EE,1,2*M);
			
		}
		else
		{
			// all the really high f sources will be detectable
			fprintf(Abright, "%.16f %.10e %f %f %e %f %f %f\n", f, fdot, theta, phi, A, iota, psi, phase);
			fprintf(Xbright, "%.16f %.10e %f %f %e %f %f %f\n", f, fdot, theta, phi, A, iota, psi, phase);
		}
		
		
		
	}
	
	printf("Removal Finished\n");
	
	printf("Number above 1e-2 Hz = %ld\n", cc1);
	
	
	Outfile = fopen("Galaxy_XAE_R1.dat","w");
	for(i=1; i< imax; i++)
	{
		f = (double)(i)/T;
		fprintf(Outfile,"%e %e %e %e %e %e %e\n", f, XfLS[2*i], XfLS[2*i+1],
				AALS[2*i], AALS[2*i+1], EELS[2*i], EELS[2*i+1]);
	}
	fclose(Outfile);
	
	
	for(i=imin; i< imax; i++)
	{
		XP[i] = (2.0*(XfLS[2*i]*XfLS[2*i] + XfLS[2*i+1]*XfLS[2*i+1]));
		AEP[i] = (2.0*(AALS[2*i]*AALS[2*i]+AALS[2*i+1]*AALS[2*i+1])); 
	}
	
	Outfile = fopen("Power_1.dat","w");
	for(i=imin; i< imax; i++)
	{
		f = (double)(i)/T;
		if(noiseFlag==1)instrument_noise(f, fstar, L, &SAE, &SXYZ);
		if(noiseFlag==2)instrument_noise2(f, fstar, L, &SAE, &SXYZ);
		XP[i] = (2.0*(XfLS[2*i]*XfLS[2*i] + XfLS[2*i+1]*XfLS[2*i+1]));
		AEP[i] = (2.0*(AALS[2*i]*AALS[2*i]+AALS[2*i+1]*AALS[2*i+1])); 
		fprintf(Outfile,"%e %e %e %e %e\n", f, XP[i], SXYZ, AEP[i], SAE);
	}
	fclose(Outfile);
	
	medianX(imin, imax, fstar, L, XP, Xnoise, Xconf);
	medianAE(imin, imax, fstar, L, AEP, AEnoise, AEconf);
	
	Outfile = fopen("Confusion_XAE_1.dat","w");
	for(i=imin; i<= imax; i++)
	{
		f = (double)(i)/T;
		fprintf(Outfile,"%e %e %e %e %e\n", f, Xnoise[i], Xconf[i], AEnoise[i], AEconf[i]);
	}
	fclose(Outfile);
	
	Outfile = fopen("Confusion_XAE_DS.dat","w");
	for(i=imin; i<= imax; i++)
	{
		if(i%100==0)
		{
			f = (double)(i)/T;
			fprintf(Outfile,"%e %e %e %e %e\n", f, Xnoise[i], Xconf[i], AEnoise[i], AEconf[i]);
		}
	}
	fclose(Outfile);
	
    return 0;
	
}



void FAST_LISA(struct lisa_orbit *orbit, double *params, long N, long M, double *XLS, double *ALS, double *ELS)
{
	
	/*   Indicies   */
	int i,j,n,m;
	/*   Carrier frequency bin  */
	long q;
	/*   Gravitational Wave basis vectors   */
	double *u,*v,*k;
	/*   Polarization basis tensors   */
	double **eplus, **ecross;
	/*   Spacecraft position and separation vector   */
	double *x, *y, *z;
	double *r12, *r13, *r21, *r23, *r31, *r32;
	/*   Dot products   */
	double *kdotx, **kdotr;
	/*   Convenient quantities   */
	double **dplus, **dcross;
	/*   GW Source data   */
	double Mc, theta, phi, psi, D, iota, Amp, Aplus, Across, f0, fdot, fddot, phio;
	double costh, sinth, cosph, sinph, cosi, cosps, sinps;
	/*   Time and distance variables   */
	double *xi, t;
	/*   Gravitational wave frequency & ratio of f and transfer frequency f*  */
	double L, fstar, *f, *fonfs;
	/*   LISA response to slow terms (Real & Imaginary pieces)   */
	//Static quantities (Re and Im)
	double DPr, DPi, DCr, DCi;
	//Time varrying quantities (Re & Im) broken up into convenient segments
	double **TR, **TI, tran1r, tran1i, tran2r, tran2i;
	//Miscellaneous constants used to speed up calculations
	double arg1, arg2, sinc, df;
	/*   Fourier coefficients before FFT and after convolution  */
	//Time series of slowly evolving terms at each vertex   
	double *data12, *data13, *data21, *data23, *data31, *data32;
	//Fourier coefficients of slowly evolving terms (numerical)
	double *a12, *a13, *a21, *a23, *a31, *a32;
	//Fourier coefficients of rapidly evolving terms (analytical)
	double *b;             
	//Fourier coefficients of entire response (convolution)
	double *c12, *c13, *c21, *c23, *c31, *c32;
	//Package cij's into proper form for TDI subroutines
	double ***d;
	/*   Miscellaneous  */
	double xm, fstep, power, aevol;
	
	
	/*   Allocating Arrays   */
	
	u = dvector(1,3); v = dvector(1,3); k = dvector(1,3);
	
	kdotx = dvector(1,3); kdotr = dmatrix(1,3,1,3);
	
	xi = dvector(1,3);
	
	f = dvector(1,3);
	
	fonfs = dvector(1,3);
	
	eplus  = dmatrix(1,3,1,3); ecross = dmatrix(1,3,1,3); 
	
	dplus  = dmatrix(1,3,1,3); dcross = dmatrix(1,3,1,3); 
	
	x = dvector(1,3); y = dvector(1,3); z = dvector(1,3);
	
	r12 = dvector(1,3); r21 = dvector(1,3); r31 = dvector(1,3);  
	r13 = dvector(1,3); r23 = dvector(1,3); r32 = dvector(1,3);
	
	TR = dmatrix(1,3,1,3); TI = dmatrix(1,3,1,3);
	
	data12 = dvector(1,2*N); data21 = dvector(1,2*N); data31 = dvector(1,2*N);
	data13 = dvector(1,2*N); data23 = dvector(1,2*N); data32 = dvector(1,2*N); 
	
	a12 = dvector(1,2*N+2); a21 = dvector(1,2*N+2); a31 = dvector(1,2*N+2);
	a13 = dvector(1,2*N+2); a23 = dvector(1,2*N+2); a32 = dvector(1,2*N+2);
	
	d = d3tensor(1,3,1,3,1,2*M);
	
	/*   Gravitational Wave source parameters   */
	
	f0 = params[0];
	costh =cos(params[1]);
	phi = params[2]; 
	Amp = params[3]; 
	cosi = cos(params[4]);
	psi = params[5];
	phio = params[6];
	fdot = params[7];
	fddot = params[8];
	
	//Calculate cos and sin of sky position, inclination, polarization
	sinth = sqrt(1.0-costh*costh);
	cosph = cos(phi);     sinph = sin(phi);
	cosps = cos(2.*psi);  sinps = sin(2.*psi);
	
	//Calculate GW polarization amplitudes
	Aplus = Amp*(1.+cosi*cosi);
	Across = -2.0*Amp*cosi;
	
	//Calculate carrier frequency bin
	q = (long)(f0*T);
	
	df = 2.0*pi*(((double)q)/T);
    
	//Calculate constant pieces of transfer functions
	DPr = Aplus*cosps;
	DPi = -Across*sinps;
	DCr = -Aplus*sinps;
	DCi = -Across*cosps;
	
	
	/*   Tensor construction for buildingslowly evolving LISA response   */
	//Gravitational Wave source basis vectors
	u[1] =  costh*cosph;  u[2] =  costh*sinph;  u[3] = -sinth;
	v[1] =  sinph;        v[2] = -cosph;        v[3] =  0.;
	k[1] = -sinth*cosph;  k[2] = -sinth*sinph;  k[3] = -costh;
	
	//GW polarization basis tensors
	for(i=1;i<=3;i++)
    {
		for(j=1;j<=3;j++)
	    {
			eplus[i][j]  = u[i]*u[j] - v[i]*v[j];
			ecross[i][j] = u[i]*v[j] + v[i]*u[j];
			
      	}
    }
	
	
	/*****************************   Main Loop   **********************************/
	L = orbit->L;
	fstar = orbit->fstar;
	for(n=1; n<=N; n++)
    {
		//First time sample must be at t=0 for phasing
		t = T*(double)(n-1)/(double)N;
		
		//Calculate position of each spacecraft at time t
		spacecraft(orbit, t, x, y, z);
		
		for(i=1; i<=3; i++) 
        {
			kdotx[i] = (x[i]*k[1]+y[i]*k[2]+z[i]*k[3])/clight;
			//Wave arrival time at spacecraft i
			xi[i]    = t - kdotx[i];
			//First order approximation to frequency at spacecraft i
			f[i]     = f0 + fdot*xi[i]+0.5*fddot*xi[i]*xi[i];
			//Ratio of true frequencto transfer frequency
			fonfs[i] = f[i]/orbit->fstar;
        }
		
		
		//Unit separation vector from spacecrafts i to j 
		r12[1] = (x[2] - x[1])/L;   r13[1] = (x[3] - x[1])/L;   r23[1] = (x[3] - x[2])/L;
		r12[2] = (y[2] - y[1])/L;   r13[2] = (y[3] - y[1])/L;   r23[2] = (y[3] - y[2])/L;
		r12[3] = (z[2] - z[1])/L;   r13[3] = (z[3] - z[1])/L;   r23[3] = (z[3] - z[2])/L;
		//Make use of symmetry
		for(i=1; i<=3; i++)
		{
			r21[i] = -r12[i];
			r31[i] = -r13[i];
			r32[i] = -r23[i];
		}
	    
		//Zero arrays to be summed
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
		//Makng use of symmetry
		dplus[2][1] = dplus[1][2];  dcross[2][1] = dcross[1][2];
		dplus[3][2] = dplus[2][3];  dcross[3][2] = dcross[2][3];
		dplus[3][1] = dplus[1][3];  dcross[3][1] = dcross[1][3];
		
		//Zero arrays to be summed
		kdotr[1][2] = kdotr[1][3] = kdotr[2][1] = kdotr[2][3] = kdotr[3][1] = kdotr[3][2] = 0.;
		for(i=1; i<=3; i++)
	    {
			kdotr[1][2] += k[i]*r12[i];   kdotr[1][3] += k[i]*r13[i];   kdotr[2][3] += k[i]*r23[i];	  
        }
		//Making use of antisymmetry  
		kdotr[2][1] = -kdotr[1][2];
		kdotr[3][1] = -kdotr[1][3];
		kdotr[3][2] = -kdotr[2][3];
		
		//Calculating Transfer function
		for(i=1; i<=3; i++)
        {
			for(j=1; j<=3; j++)
			{
				if(i!=j)
				{
					//Argument of transfer function
					arg1 = 0.5*fonfs[i]*(1 - kdotr[i][j]);
					//Argument of complex exponentials
					arg2 = 2.0*pi*f0*xi[i]+pi*fdot*xi[i]*xi[i]+pi*fddot*xi[i]*xi[i]*xi[i]/3.0 +phio - df*t;
					//Transfer function
					sinc = 0.25*sin(arg1)/arg1;
					//Evolution of amplitude
					aevol = 1.0 + 0.66666666666666666666*fdot/f0*xi[i];  
					///Real and imaginary pieces of time series (no complex exponential)
					tran1r = aevol*(dplus[i][j]*DPr + dcross[i][j]*DCr);
					tran1i = aevol*(dplus[i][j]*DPi + dcross[i][j]*DCi);
					//Real and imaginry components of complex exponential
					tran2r = cos(arg1 + arg2);
					tran2i = sin(arg1 + arg2);    
					//Real & Imaginary part of the slowly evolving signal
					TR[i][j] = sinc*(tran1r*tran2r - tran1i*tran2i);
					TI[i][j] = sinc*(tran1r*tran2i + tran1i*tran2r);
				}
			}
		}
		
		//Fill  time series data arrays with slowly evolving signal.  
		//dataij corresponds to fractional arm length difference yij
		data12[2*n-1] = TR[1][2];   data21[2*n-1] = TR[2][1];   data31[2*n-1] = TR[3][1];   
		data12[2*n]   = TI[1][2];   data21[2*n]   = TI[2][1];   data31[2*n]   = TI[3][1];
		data13[2*n-1] = TR[1][3];   data23[2*n-1] = TR[2][3];   data32[2*n-1] = TR[3][2];
		data13[2*n]   = TI[1][3];   data23[2*n]   = TI[2][3];   data32[2*n]   = TI[3][2];
    }
	
	/*   Numerical Fourier transform of slowly evolving signal   */
	dfour1(data12, N, -1);  dfour1(data21, N, -1);  dfour1(data31, N, -1);
	dfour1(data13, N, -1);  dfour1(data23, N, -1);  dfour1(data32, N, -1);  
	//Unpack arrays from dfour1.c and normalize
	for(i=1; i<=N; i++)
    {
		a12[i] = data12[N+i]/(double)N;  a21[i] = data21[N+i]/(double)N;  a31[i] = data31[N+i]/(double)N;
		a12[i+N] = data12[i]/(double)N;  a21[i+N] = data21[i]/(double)N;  a31[i+N] = data31[i]/(double)N;
		a13[i] = data13[N+i]/(double)N;  a23[i] = data23[N+i]/(double)N;  a32[i] = data32[N+i]/(double)N;
		a13[i+N] = data13[i]/(double)N;  a23[i+N] = data23[i]/(double)N;  a32[i+N] = data32[i]/(double)N;
    }
	
	a12[2*N+1] = data12[N+1]/(double)N;  a21[2*N+1] = data21[N+1]/(double)N;  a31[2*N+1] = data31[N+1]/(double)N;
	a12[2*N+2] = data12[N+2]/(double)N;  a21[2*N+2] = data21[N+2]/(double)N;  a31[2*N+2] = data31[N+2]/(double)N;
	a13[2*N+1] = data13[N+1]/(double)N;  a23[2*N+1] = data23[N+1]/(double)N;  a32[2*N+1] = data32[N+1]/(double)N;
	a13[2*N+2] = data13[N+2]/(double)N;  a23[2*N+2] = data23[N+2]/(double)N;  a32[2*N+2] = data32[N+2]/(double)N;
	
	/*   Renormalize so that the resulting time series is real   */
	for(i=1; i<=2*M; i++)
    {
		d[1][2][i] = 0.5*a12[i];  d[2][1][i] = 0.5*a21[i];  d[3][1][i] = 0.5*a31[i];
		d[1][3][i] = 0.5*a13[i];  d[2][3][i] = 0.5*a23[i];  d[3][2][i] = 0.5*a32[i];
    }
	
	/*   Call subroutines for synthesizing different TDI data channels  */
	/*   X Y Z-Channel   */
	XYZ(L, fstar, d, f0, q, M, XLS, ALS, ELS);
	
	
	/*   Deallocate Arrays   */
	
	free_dvector(u,1,3); free_dvector(v,1,3); free_dvector(k,1,3);
	
	free_dvector(kdotx,1,3); free_dmatrix(kdotr,1,3,1,3);
	
	free_dvector(xi,1,3);
	
	free_dvector(f,1,3);
	
	free_dvector(fonfs,1,3);
	
	free_dmatrix(eplus,1,3,1,3); free_dmatrix(ecross,1,3,1,3); 
	
	free_dmatrix(dplus,1,3,1,3); free_dmatrix(dcross,1,3,1,3); 
	
	free_dvector(x,1,3); free_dvector(y,1,3); free_dvector(z,1,3);
	
	free_dvector(r12,1,3); free_dvector(r21,1,3); free_dvector(r31,1,3);  
	free_dvector(r13,1,3); free_dvector(r23,1,3); free_dvector(r32,1,3);
	
	free_dmatrix(TR,1,3,1,3); free_dmatrix(TI,1,3,1,3);
	
	free_dvector(data12,1,2*N); free_dvector(data21,1,2*N); free_dvector(data31,1,2*N);
	free_dvector(data13,1,2*N); free_dvector(data23,1,2*N); free_dvector(data32,1,2*N); 
	
	free_dvector(a12,1,2*N+2); free_dvector(a21,1,2*N+2); free_dvector(a31,1,2*N+2);
	free_dvector(a13,1,2*N+2); free_dvector(a23,1,2*N+2); free_dvector(a32,1,2*N+2);
	
	free_d3tensor(d,1,3,1,3,1,2*M);
	
	
	return;
}

/*************************************************************************************/
/*                                                                                   */
/*                                    Subroutines                                    */
/*                                                                                   */
/*************************************************************************************/ 

void spacecraft(struct lisa_orbit *orbit, double tint, double *xint, double *yint, double *zint)
{
	int i;
	for(i=0; i<3; i++) 
	{
		LISA_splint(orbit->t, orbit->x[i], orbit->dx[i], orbit->N, tint, &xint[i+1]);
		LISA_splint(orbit->t, orbit->y[i], orbit->dy[i], orbit->N, tint, &yint[i+1]);
		LISA_splint(orbit->t, orbit->z[i], orbit->dz[i], orbit->N, tint, &zint[i+1]);
	}	
}

void initialize_orbit(char OrbitFile[], struct lisa_orbit *orbit)
{
	int n,i;
	double junk;
	
	FILE *infile = fopen(OrbitFile,"r");
	
	//how big is the file
	n=0;
	while(!feof(infile))
	{
		/*columns of orbit file:
		 t sc1x sc1y sc1z sc2x sc2y sc2z sc3x sc3y sc3z 
		 */
		n++;
		fscanf(infile,"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",&junk,&junk,&junk,&junk,&junk,&junk,&junk,&junk,&junk,&junk);
	}
	n--;
	rewind(infile);
	
	//allocate memory for local workspace
	int N = n;
	double *t  = dvector(0,N-1);
	double **x  = dmatrix(0,2,0,N-1);
	double **y  = dmatrix(0,2,0,N-1);
	double **z  = dmatrix(0,2,0,N-1);
	double **dx = dmatrix(0,2,0,N-1);
	double **dy = dmatrix(0,2,0,N-1);
	double **dz = dmatrix(0,2,0,N-1);
	
	//allocate memory for orbit structure
	orbit->N = n;
	
	orbit->t  = dvector(0,orbit->N-1);
	orbit->x  = dmatrix(0,2,0,orbit->N-1);
	orbit->y  = dmatrix(0,2,0,orbit->N-1);
	orbit->z  = dmatrix(0,2,0,orbit->N-1);
	orbit->dx = dmatrix(0,2,0,orbit->N-1);
	orbit->dy = dmatrix(0,2,0,orbit->N-1);
	orbit->dz = dmatrix(0,2,0,orbit->N-1);
	
	//read in orbits
	for(n=0; n<orbit->N; n++)
    {
		//First time sample must be at t=0 for phasing
		fscanf(infile,"%lg",&t[n]);
		for(i=0; i<3; i++) fscanf(infile,"%lg %lg %lg",&x[i][n],&y[i][n],&z[i][n]);
		
		orbit->t[n] = t[n];
		
		//Repackage orbit positions into arrays for interpolation
		for(i=0; i<3; i++) 
		{
			orbit->x[i][n] = x[i][n]; 
			orbit->y[i][n] = y[i][n]; 
			orbit->z[i][n] = z[i][n]; 			
		}
	}
	fclose(infile);
	
	//calculate derivatives for cubic spline
	for(i=0; i<3; i++) 
	{
		LISA_spline(t, orbit->x[i], N, 1.e30, 1.e30, orbit->dx[i]);
		LISA_spline(t, orbit->y[i], N, 1.e30, 1.e30, orbit->dy[i]);
		LISA_spline(t, orbit->z[i], N, 1.e30, 1.e30, orbit->dz[i]);
	}
	
	//calculate average arm length
	printf("estimating average armlengths -- assumes evenly sampled orbits\n");
	double L12=0.0;
	double L23=0.0;
	double L31=0.0;
	double x1,x2,x3,y1,y2,y3,z1,z2,z3;
	for(n=0; n<orbit->N; n++)
	{
		x1 = orbit->x[0][n];
		x2 = orbit->x[1][n];
		x3 = orbit->x[2][n];
		
		y1 = orbit->y[0][n];
		y2 = orbit->y[1][n];
		y3 = orbit->y[2][n];
		
		z1 = orbit->z[0][n];
		z2 = orbit->z[1][n];
		z3 = orbit->z[2][n];
		
		
		//L12
		L12 += sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) ); 
		
		//L23
		L23 += sqrt( (x3-x2)*(x3-x2) + (y3-y2)*(y3-y2) + (z3-z2)*(z3-z2) ); 
		
		//L31
		L31 += sqrt( (x1-x3)*(x1-x3) + (y1-y3)*(y1-y3) + (z1-z3)*(z1-z3) ); 
	}
	L12 /= (double)orbit->N;
	L31 /= (double)orbit->N;
	L23 /= (double)orbit->N;
	
	printf("Average arm lengths for the constellation:\n");
	printf("  L12 = %g\n",L12);
	printf("  L31 = %g\n",L31);
	printf("  L23 = %g\n",L23);
	
	//are the armlenghts consistent?
	double L = (L12+L31)/2.;//+L23)/3.;
	printf("Fractional deviation from average armlength for each side:\n");
	printf("  L12 = %g\n",fabs(L12-L)/L);
	printf("  L31 = %g\n",fabs(L31-L)/L);
	printf("  L23 = %g\n",fabs(L23-L)/L);
	
	//store armlenght & transfer frequency in orbit structure.
	orbit->L     = L;
	orbit->fstar = clight/(2.0*pi*L); 
	
	//free local memory
	free_dvector(t,0,N-1);
	free_dmatrix(x ,0,2,0,N-1);
	free_dmatrix(y ,0,2,0,N-1);
	free_dmatrix(z ,0,2,0,N-1);
	free_dmatrix(dx,0,2,0,N-1);
	free_dmatrix(dy,0,2,0,N-1);
	free_dmatrix(dz,0,2,0,N-1);
}

void LISA_spline(double *x, double *y, int n, double yp1, double ypn, double *y2)
// Unlike NR version, assumes zero-offset arrays.  CHECK THAT THIS IS CORRECT.
{
	int i, k;
	double p, qn, sig, un, *u;
	u = dvector(0, n-2);
	// Boundary conditions: Check which is best.
	if (yp1 > 0.99e30)
		y2[0] = u[0] = 0.0;
	else {
		y2[0] = -0.5;
		u[0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}
	
	
	for(i = 1; i < n-1; i++) {
		sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p = sig*y2[i-1] + 2.0;
		y2[i] = (sig-1.0)/p;
		u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	// More boundary conditions.
	if (ypn > 0.99e30)
		qn = un = 0.0;
	else {
		qn = 0.5;
		un = (3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	for (k = n-2; k >= 0; k--)
		y2[k] = y2[k]*y2[k+1]+u[k];
	free_dvector(u, 0, n-2);
}

void LISA_splint(double *xa, double *ya, double *y2a, int n, double x, double *y)
{
	// Unlike NR version, assumes zero-offset arrays.  CHECK THAT THIS IS CORRECT.
	int klo,khi,k;
	double h,b,a;
	
	klo=0;
	khi=n-1;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}


void XYZ(double L, double fstar, double ***d, double f0, long q, long M, double *XLS, double *ALS, double *ELS)
{ 
	int i;  
	double fonfs, sqT;
	double ReX, ImX, ReY, ImY, ReZ, ImZ;
	double c3, s3, c2, s2, c1, s1;
	double f,power;
	double *X, *Y, *Z;
	double *YLS, *ZLS;
	double phiLS, phiSL, cLS, sLS, cSL, sSL;
	
	X = dvector(1,2*M);  Y = dvector(1,2*M);  Z = dvector(1,2*M);
	YLS = dvector(1,2*M);  ZLS = dvector(1,2*M);
	
	phiLS = 2.0*pi*f0*(dt/2.0-L/clight);
	/* phiSL = pi/2.0-2.0*pi*f0*(L/clight); */
	cLS = cos(phiLS);
	sLS = sin(phiLS);
	/*cSL = cos(phiSL);
     sSL = sin(phiSL); */
	
	sqT = sqrt(T);
	for(i=1; i<=M; i++)
    {
		
		f = ((double)(q + i-1 - M/2))/T;
		fonfs = f/fstar;
		c3 = cos(3.*fonfs);  c2 = cos(2.*fonfs);  c1 = cos(1.*fonfs);
		s3 = sin(3.*fonfs);  s2 = sin(2.*fonfs);  s1 = sin(1.*fonfs);
		
		
		X[2*i-1] = (d[1][2][2*i-1]-d[1][3][2*i-1])*c3 + (d[1][2][2*i]-d[1][3][2*i])*s3 +
		(d[2][1][2*i-1]-d[3][1][2*i-1])*c2 + (d[2][1][2*i]-d[3][1][2*i])*s2 +
		(d[1][3][2*i-1]-d[1][2][2*i-1])*c1 + (d[1][3][2*i]-d[1][2][2*i])*s1 +
		(d[3][1][2*i-1]-d[2][1][2*i-1]);
		X[2*i]   = (d[1][2][2*i]-d[1][3][2*i])*c3 - (d[1][2][2*i-1]-d[1][3][2*i-1])*s3 +
		(d[2][1][2*i]-d[3][1][2*i])*c2 - (d[2][1][2*i-1]-d[3][1][2*i-1])*s2 +
		(d[1][3][2*i]-d[1][2][2*i])*c1 - (d[1][3][2*i-1]-d[1][2][2*i-1])*s1 +
		(d[3][1][2*i]-d[2][1][2*i]);
		
		Y[2*i-1] = (d[2][3][2*i-1]-d[2][1][2*i-1])*c3 + (d[2][3][2*i]-d[2][1][2*i])*s3 +
		(d[3][2][2*i-1]-d[1][2][2*i-1])*c2 + (d[3][2][2*i]-d[1][2][2*i])*s2+
		(d[2][1][2*i-1]-d[2][3][2*i-1])*c1 + (d[2][1][2*i]-d[2][3][2*i])*s1+
		(d[1][2][2*i-1]-d[3][2][2*i-1]);
		Y[2*i]   = (d[2][3][2*i]-d[2][1][2*i])*c3 - (d[2][3][2*i-1]-d[2][1][2*i-1])*s3+
		(d[3][2][2*i]-d[1][2][2*i])*c2 - (d[3][2][2*i-1]-d[1][2][2*i-1])*s2+
		(d[2][1][2*i]-d[2][3][2*i])*c1 - (d[2][1][2*i-1]-d[2][3][2*i-1])*s1+
		(d[1][2][2*i]-d[3][2][2*i]);
		
		Z[2*i-1] = (d[3][1][2*i-1]-d[3][2][2*i-1])*c3 + (d[3][1][2*i]-d[3][2][2*i])*s3+
		(d[1][3][2*i-1]-d[2][3][2*i-1])*c2 + (d[1][3][2*i]-d[2][3][2*i])*s2+
		(d[3][2][2*i-1]-d[3][1][2*i-1])*c1 + (d[3][2][2*i]-d[3][1][2*i])*s1+
		(d[2][3][2*i-1]-d[1][3][2*i-1]);
		
		Z[2*i]   = (d[3][1][2*i]-d[3][2][2*i])*c3 - (d[3][1][2*i-1]-d[3][2][2*i-1])*s3+
		(d[1][3][2*i]-d[2][3][2*i])*c2 - (d[1][3][2*i-1]-d[2][3][2*i-1])*s2+
		(d[3][2][2*i]-d[3][1][2*i])*c1 - (d[3][2][2*i-1]-d[3][1][2*i-1])*s1+
		(d[2][3][2*i]-d[1][3][2*i]);
        
		
		/* XSL[2*i-1] = 2.0*fonfs*(X[2*i-1]*cSL-X[2*i]*sSL);
		 XSL[2*i] = -2.0*fonfs*(X[2*i-1]*sSL+X[2*i]*cSL);
		 YSL[2*i-1] = 2.0*fonfs*(Y[2*i-1]*cSL-Y[2*i]*sSL);
		 YSL[2*i] = -2.0*fonfs*(Y[2*i-1]*sSL+Y[2*i]*cSL);
		 ZSL[2*i-1] = 2.0*fonfs*(Z[2*i-1]*cSL-Z[2*i]*sSL);
		 ZSL[2*i] = -2.0*fonfs*(Z[2*i-1]*sSL+Z[2*i]*cSL); */
		
		XLS[2*i-1] = (X[2*i-1]*cLS-X[2*i]*sLS);
		XLS[2*i] = -(X[2*i-1]*sLS+X[2*i]*cLS);
		YLS[2*i-1] = (Y[2*i-1]*cLS-Y[2*i]*sLS);
		YLS[2*i] = -(Y[2*i-1]*sLS+Y[2*i]*cLS);
		ZLS[2*i-1] = (Z[2*i-1]*cLS-Z[2*i]*sLS);
		ZLS[2*i] = -(Z[2*i-1]*sLS+Z[2*i]*cLS);
		
		
    }
	
	
	for(i=1; i<=2*M; i++)
	{
		// A channel
		ALS[i] = (2.0*XLS[i] - YLS[i] - ZLS[i])/3.0;
		// E channel
		ELS[i] = (ZLS[i]-YLS[i])/sq3;
	}
	
	free_dvector(X,1,2*M);  free_dvector(Y,1,2*M);  free_dvector(Z,1,2*M);
	free_dvector(YLS,1,2*M);  free_dvector(ZLS,1,2*M);
	
}

void instrument_noise(double f, double fstar, double L, double *SAE, double *SXYZ)
{
	//Power spectral density of the detector noise and transfer frequency
	double Sn, red, confusion_noise;
	double f1, f2;
	double A1, A2, slope;
	FILE *outfile;
	
	red = 0.0;//(2.0*pi*1.0e-4)*(2.0*pi*1.0e-4);
	
	
	// Calculate the power spectral density of the detector noise at the given frequency
	
	//No reddening
	if(redden==0)
	{
		*SAE = 16.0/3.0*pow(sin(f/fstar),2.0)*( ( (2.0+cos(f/fstar))*Sps + 2.0*(3.0+2.0*cos(f/fstar)+cos(2.0*f/fstar))*Sacc*(1.0/pow(2.0*pi*f,4)+ red/pow(2.0*pi*f,6))) / pow(2.0*L,2.0));
		*SXYZ = 4.0*pow(sin(f/fstar),2.0)*( ( 4.0*Sps + 8.0*(1.0+pow(cos(f/fstar),2.0))*Sacc*(1.0/pow(2.0*pi*f,4)+ red/pow(2.0*pi*f,6))) / pow(2.0*L,2.0));
	}
	
	//NGO reddening
	if(redden==1)
	{
		*SAE = 16.0/3.0*pow(sin(f/fstar),2.0)*( ( (2.0+cos(f/fstar))*Sps + 2.0*(3.0+2.0*cos(f/fstar)+cos(2.0*f/fstar))*Sacc*(1.+(0.0001/f))*(1.0/pow(2.0*pi*f,4)+ red/pow(2.0*pi*f,6))) / pow(2.0*L,2.0));
		*SXYZ = 4.0*pow(sin(f/fstar),2.0)*( ( 4.0*Sps + 8.0*(1.0+pow(cos(f/fstar),2.0))*Sacc*(1.+(0.0001/f))*(1.0/pow(2.0*pi*f,4)+ red/pow(2.0*pi*f,6))) / pow(2.0*L,2.0));
	}
	
	//LAGRANGE reddening
	if(redden==2)
	{
		*SAE = 16.0/3.0*pow(sin(f/fstar),2.0)*( ( (2.0+cos(f/fstar))*Sps + 2.0*(3.0+2.0*cos(f/fstar)+cos(2.0*f/fstar))*Sacc*pow(f,-3./2.)*(1.0/pow(2.0*pi*f,4)+ red/pow(2.0*pi*f,6))) / pow(2.0*L,2.0));
		*SXYZ = 4.0*pow(sin(f/fstar),2.0)*( ( 4.0*Sps + 8.0*(1.0+pow(cos(f/fstar),2.0))*Sacc*pow(f,-3./2.)*(1.0/pow(2.0*pi*f,4)+ red/pow(2.0*pi*f,6))) / pow(2.0*L,2.0));
	}
	
	if(*SAE < 1e-45) *SAE = 1e-45;
	if(*SXYZ< 1e-45) *SXYZ= 1e-45;
	
}



#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void dfour1(double data[], unsigned long nn, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	double tempr,tempi;
	
	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}
#undef SWAP

void KILL(char* Message)
{
	printf("\a\n");
	printf("%s",Message);
	printf("Terminating the program.\n\n\n");
	exit(1);
	
	
	return;
}

double gasdev2(long *idum)
{
	double ran2(long *idum);
	static int iset=0;
	static double gset;
	double fac, rsq, v1, v2;
	
	if(*idum < 0) iset = 0;
	if(iset == 0){
		do{
			v1 = 2.0 * ran2(idum)-1.0;
			v2 = 2.0 * ran2(idum)-1.0;
			rsq = v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset = v1 * fac;
		iset = 1;
		return(v2*fac);
	} else {
		iset = 0;
		return (gset);
	}
}

double ran2(long *idum)
{
    int j;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    double temp;
	
    if (*idum <= 0) {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        idum2=(*idum);
        for (j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if (*idum < 0) *idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}

void medianX(long imin, long imax, double fstar, double L, double *XP, double *Xnoise, double *Xconf)
{
	double f, logf;
	double XC, SAE, SXYZ;
	double chi, scale;
	long i, j, k;
	long segs;
	long rseed;
	int Npoly;
	double av;
	double *XX;
	double *fdata, *mdata, *pcx, *pcy, *inst;
	double *fsams;
	double **fpows;
	double chix, chiy, fit, alpha, beta, mul, conf;
	double lfmin, lfmax, x, dlf, lf, ln4;
	FILE *Xfile;
	
	XX = dvector(0,100);
	
	rseed = -546214;
	
	segs = (int)((double)(imax-imin)/101.0);
	
	lfmin = log((double)(imin-101)/T);
	lfmax = log((double)(imin+101*(segs))/T);
	
	
	Npoly = 30;
	
	dlf = (lfmax-lfmin)/(double)(Npoly);
	ln4 = log(1.0e-4);
	
	fdata = dvector(0,segs-1);
	mdata = dvector(0,segs-1);
	inst = dvector(0,segs-1);
	pcx = dvector(0,Npoly);
	pcy = dvector(0,Npoly);
	
	for(i=0; i < segs; i++)
	{
		for(j=0; j<=100; j++) XX[j] = XP[imin+101*i+j];
		f = (double)(imin+101*i-50)/T;
		if(noiseFlag==1)instrument_noise(f, fstar, L, &SAE, &SXYZ);
		if(noiseFlag==2)instrument_noise2(f, fstar, L, &SAE, &SXYZ);
		inst[i] = log(SXYZ*1.0e40);
		chi=quickselect(XX, 101, 51);
		//printf("%e %e\n", f, chi/0.72);
		fdata[i] = log(f);
		mdata[i] = log(chi/0.72*1.0e40);
	}
	
	// initialize fit
	for(i=1; i < Npoly; i++)
	{
		f = exp(lfmin+(double)(i)*dlf);
		j = (long)((f*T-(double)(imin-50))/101.0);
		//printf("%ld %ld\n", i, j);
		pcx[i] = mdata[j];
		//printf("%e %e\n", f, exp(pcx[i])*1.0e-40);
	}
	pcx[0] = pcx[1];
	pcx[Npoly] = pcx[Npoly-1];
	//printf("%ld\n", segs);
	
	
	chix = 0.0;
	for(i=0; i < segs; i++)
	{
		lf = log((double)(imin+101*i-50)/T);
        j = (long)floor((lf-lfmin)/dlf);
		fit = pcx[j]+((pcx[j+1]-pcx[j])/dlf)*(lf-(lfmin+(double)(j)*dlf));
        f = exp(-1.0*((lfmin+((double)(j)+0.5)*dlf)-ln4));
        chix += 500.0*(mdata[i] - fit)*(mdata[i] - fit)*f;
        //printf("%ld %e\n", i, (lf-(lfmin+(double)(j)*dlf))/dlf);
        //printf("%e %e %e\n", exp(lf), exp(mdata[i])*1.0e-40, exp(fit)*1.0e-40);
	}
	
	for(k=0; k < 10000; k++)
	{
		beta = pow(10.0,-0.0001*(10000.0-(double)(k))/10000.0);
		
		for(j=0; j <= Npoly; j++) pcy[j] = pcx[j];
		j = (long)((double)(Npoly+1)*ran2(&rseed));
		mul = 1.0;
		alpha = ran2(&rseed);
		if(alpha > 0.3) mul = 10.0;
		if(alpha > 0.7) mul = 0.01;
		pcy[j] = pcx[j]+mul*0.01*gasdev2(&rseed)/sqrt(beta);
		
        chiy = 0.0;
        for(i=0; i < segs; i++)
		{
			lf = log((double)(imin+101*i-50)/T);
			j = (long)floor((lf-lfmin)/dlf);
			fit = pcy[j]+((pcy[j+1]-pcy[j])/dlf)*(lf-(lfmin+(double)(j)*dlf));
			f = exp(-1.0*((lfmin+((double)(j)+0.5)*dlf)-ln4));
			chiy += 500.0*(mdata[i] - fit)*(mdata[i] - fit)*f;
        }
		
        alpha = log(ran2(&rseed));
		
        if(beta*(chix-chiy) > alpha)
		{
			chix = chiy;
            for(j=0; j <= Npoly; j++) pcx[j] = pcy[j];
		}
		//if(k%100==0) printf("%ld %.10e %.10e\n", k, chix, chiy);
	}
	
	//printf("%ld %.10e %.10e\n", k, chix, chiy);
	
	Xfile = fopen("Xfit.dat","w");
	for(i=0; i < segs; i++)
	{
		lf = log((double)(imin+101*i-50)/T);
        j = (long)floor((lf-lfmin)/dlf);
		fit = pcx[j]+((pcx[j+1]-pcx[j])/dlf)*(lf-(lfmin+(double)(j)*dlf));
        f = exp(lf);
		if(noiseFlag==1)instrument_noise(f, fstar, L, &SAE, &SXYZ);
		if(noiseFlag==2)instrument_noise2(f, fstar, L, &SAE, &SXYZ);
        fprintf(Xfile, "%e %e %e %e\n", exp(lf), exp(fit)*1.0e-40, exp(mdata[i])*1.0e-40, SXYZ);
	}
	fclose(Xfile);
	
    for(i=imin; i <= imax; i++)
	{
		f = (double)(i)/T;
		lf = log(f);
		j = (long)floor((lf-lfmin)/dlf);
		fit = pcx[j]+((pcx[j+1]-pcx[j])/dlf)*(lf-(lfmin+(double)(j)*dlf));
		if(noiseFlag==1)instrument_noise(f, fstar, L, &SAE, &SXYZ);
		if(noiseFlag==2)instrument_noise2(f, fstar, L, &SAE, &SXYZ);
		alpha = exp(fit)*1.0e-40;
		conf = alpha -SXYZ;
		if(conf < SXYZ/30.0) conf = 1.0e-46;
		Xnoise[i] = alpha;
		Xconf[i] = conf;
	}
	
    free_dvector(XX, 0,100);
    free_dvector(fdata, 0,segs-1);
    free_dvector(mdata, 0,segs-1);
    free_dvector(pcx, 0,Npoly);
    free_dvector(pcy, 0,Npoly);
	
	return;
}



void medianAE(long imin, long imax, double fstar, double L, double *AEP, double *AEnoise, double *AEconf)
{
	double f, logf;
	double XC, SAE, SXYZ;
	double chi, scale;
	long i, j, k;
	long segs;
	long rseed;
	int Npoly;
	double av;
	double *XX;
	double *fdata, *mdata, *pcx, *pcy, *inst;
	double *fsams;
	double **fpows;
	double chix, chiy, fit, alpha, beta, mul, conf;
	double lfmin, lfmax, x, dlf, lf, ln4;
	FILE *Xfile;
	
	XX = dvector(0,100);
	
	rseed = -546214;
	
	segs = (int)((double)(imax-imin)/101.0);
	
	lfmin = log((double)(imin-101)/T);
	lfmax = log((double)(imin+101*(segs))/T);
	
	
	Npoly = 30;
	
	dlf = (lfmax-lfmin)/(double)(Npoly);
	ln4 = log(1.0e-4);
	
	fdata = dvector(0,segs-1);
	mdata = dvector(0,segs-1);
	inst = dvector(0,segs-1);
	pcx = dvector(0,Npoly);
	pcy = dvector(0,Npoly);
	
	for(i=0; i < segs; i++)
	{
		for(j=0; j<=100; j++) XX[j] = AEP[imin+101*i+j];
		f = (double)(imin+101*i-50)/T;
		if(noiseFlag==1)instrument_noise(f, fstar, L, &SAE, &SXYZ);
		if(noiseFlag==2)instrument_noise2(f, fstar, L, &SAE, &SXYZ);
		inst[i] = log(SAE*1.0e40);
		chi=quickselect(XX, 101, 51);
		//printf("%e %e\n", f, chi/0.72);
		fdata[i] = log(f);
		mdata[i] = log(chi/0.72*1.0e40);
	}
	
	// initialize fit
	for(i=1; i < Npoly; i++)
	{
		f = exp(lfmin+(double)(i)*dlf);
		j = (long)((f*T-(double)(imin-50))/101.0);
		//printf("%ld %ld\n", i, j);
		pcx[i] = mdata[j];
		//printf("%e %e\n", f, exp(pcx[i])*1.0e-40);
	}
	pcx[0] = pcx[1];
	pcx[Npoly] = pcx[Npoly-1];
	//printf("%ld\n", segs);
	
	
	chix = 0.0;
	for(i=0; i < segs; i++)
	{
		lf = log((double)(imin+101*i-50)/T);
        j = (long)floor((lf-lfmin)/dlf);
		fit = pcx[j]+((pcx[j+1]-pcx[j])/dlf)*(lf-(lfmin+(double)(j)*dlf));
        f = exp(-1.0*((lfmin+((double)(j)+0.5)*dlf)-ln4));
        chix += 500.0*(mdata[i] - fit)*(mdata[i] - fit)*f;
        //printf("%ld %e\n", i, (lf-(lfmin+(double)(j)*dlf))/dlf);
        //printf("%e %e %e\n", exp(lf), exp(mdata[i])*1.0e-40, exp(fit)*1.0e-40);
	}
	
	for(k=0; k < 10000; k++)
	{
		beta = pow(10.0,-0.0001*(10000.0-(double)(k))/10000.0);
		
		for(j=0; j <= Npoly; j++) pcy[j] = pcx[j];
		j = (long)((double)(Npoly+1)*ran2(&rseed));
		mul = 1.0;
		alpha = ran2(&rseed);
		if(alpha > 0.3) mul = 10.0;
		if(alpha > 0.7) mul = 0.01;
		pcy[j] = pcx[j]+mul*0.01*gasdev2(&rseed)/sqrt(beta);
		
        chiy = 0.0;
        for(i=0; i < segs; i++)
		{
			lf = log((double)(imin+101*i-50)/T);
			j = (long)floor((lf-lfmin)/dlf);
			fit = pcy[j]+((pcy[j+1]-pcy[j])/dlf)*(lf-(lfmin+(double)(j)*dlf));
			f = exp(-1.0*((lfmin+((double)(j)+0.5)*dlf)-ln4));
			chiy += 500.0*(mdata[i] - fit)*(mdata[i] - fit)*f;
        }
		
        alpha = log(ran2(&rseed));
		
        if(beta*(chix-chiy) > alpha)
		{
			chix = chiy;
            for(j=0; j <= Npoly; j++) pcx[j] = pcy[j];
		}
		//if(k%100==0) printf("%ld %.10e %.10e\n", k, chix, chiy);
	}
	
	//printf("%ld %.10e %.10e\n", k, chix, chiy);
	
	Xfile = fopen("Afit.dat","w");
	for(i=0; i < segs; i++)
	{
		lf = log((double)(imin+101*i-50)/T);
        j = (long)floor((lf-lfmin)/dlf);
		fit = pcx[j]+((pcx[j+1]-pcx[j])/dlf)*(lf-(lfmin+(double)(j)*dlf));
        f = exp(lf);
		if(noiseFlag==1)instrument_noise(f, fstar, L, &SAE, &SXYZ);
		if(noiseFlag==2)instrument_noise2(f, fstar, L, &SAE, &SXYZ);
        fprintf(Xfile, "%e %e %e %e\n", exp(lf), exp(fit)*1.0e-40, exp(mdata[i])*1.0e-40, SAE);
	}
	fclose(Xfile);
	
    for(i=imin; i <= imax; i++)
	{
		f = (double)(i)/T;
		lf = log(f);
		j = (long)floor((lf-lfmin)/dlf);
		fit = pcx[j]+((pcx[j+1]-pcx[j])/dlf)*(lf-(lfmin+(double)(j)*dlf));
		if(noiseFlag==1)instrument_noise(f, fstar, L, &SAE, &SXYZ);
		if(noiseFlag==2)instrument_noise2(f, fstar, L, &SAE, &SXYZ);
		alpha = exp(fit)*1.0e-40;
		conf = alpha -SAE;
		if(conf < SAE/30.0) conf = 1.0e-46;
		AEnoise[i] = alpha;
		AEconf[i] = conf;
	}
	
    free_dvector(XX, 0,100);
    free_dvector(fdata, 0,segs-1);
    free_dvector(mdata, 0,segs-1);
    free_dvector(pcx, 0,Npoly);
    free_dvector(pcy, 0,Npoly);
	
	return;
}







#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

// for n odd, median given by k = (n+1)/2. for n even, median is average of k=n/2, k=n/2+1 elements

double quickselect(double *arr, int n, int k)
{
	unsigned long i,ir,j,l,mid;
	double a,temp;
	
	l=0;
	ir=n-1;
	for(;;) {
		if (ir <= l+1) { 
			if (ir == l+1 && arr[ir] < arr[l]) {
				SWAP(arr[l],arr[ir]);
			}
			return arr[k];
		}
		else {
			mid=(l+ir) >> 1; 
			SWAP(arr[mid],arr[l+1]);
			if (arr[l] > arr[ir]) {
				SWAP(arr[l],arr[ir]);
			}
			if (arr[l+1] > arr[ir]) {
				SWAP(arr[l+1],arr[ir]);
			}
			if (arr[l] > arr[l+1]) {
				SWAP(arr[l],arr[l+1]);
			}
			i=l+1; 
			j=ir;
			a=arr[l+1]; 
			for (;;) { 
				do i++; while (arr[i] < a); 
				do j--; while (arr[j] > a); 
				if (j < i) break; 
				SWAP(arr[i],arr[j]);
			} 
			arr[l+1]=arr[j]; 
			arr[j]=a;
			if (j >= k) ir=j-1; 
			if (j <= k) l=i;
		}
	}
}

void instrument_noise2(double f, double fstar, double L, double *SAE, double *SXYZ)
{
	//Power spectral density of the detector noise and transfer frequency
	double Sn, red, confusion_noise;
	double f1, f2;
	double A1, A2, slope;
	FILE *outfile;
	
	red = 0.0;//(2.0*pi*1.0e-4)*(2.0*pi*1.0e-4);
	
	
	// Calculate the power spectral density of the detector noise at the given frequency
	
	//No reddening
	if(redden==0)
	{
		*SAE = 16.0/3.0/*pow(sin(f/fstar),2.0)*/*( ( (2.0+cos(f/fstar))*Sps + 2.0*(3.0+2.0*cos(f/fstar)+cos(2.0*f/fstar))*Sacc*(1.0/pow(2.0*pi*f,4)+ red/pow(2.0*pi*f,6))) / pow(2.0*L,2.0));
		*SXYZ = 4.0/*pow(sin(f/fstar),2.0)*/*( ( 4.0*Sps + 8.0*(1.0+pow(cos(f/fstar),2.0))*Sacc*(1.0/pow(2.0*pi*f,4)+ red/pow(2.0*pi*f,6))) / pow(2.0*L,2.0));
	}
	
	//NGO reddening
	if(redden==1)
	{
		*SAE = 16.0/3.0/*pow(sin(f/fstar),2.0)*/*( ( (2.0+cos(f/fstar))*Sps + 2.0*(3.0+2.0*cos(f/fstar)+cos(2.0*f/fstar))*Sacc*(1.+(fred/f))*(1.0/pow(2.0*pi*f,4)+ red/pow(2.0*pi*f,6))) / pow(2.0*L,2.0));
		*SXYZ = 4.0/*pow(sin(f/fstar),2.0)*/*( ( 4.0*Sps + 8.0*(1.0+pow(cos(f/fstar),2.0))*Sacc*(1.+(fred/f))*(1.0/pow(2.0*pi*f,4)+ red/pow(2.0*pi*f,6))) / pow(2.0*L,2.0));
	}
	
	//LAGRANGE reddening
	if(redden==2)
	{
		*SAE = 16.0/3.0/*pow(sin(f/fstar),2.0)*/*( ( (2.0+cos(f/fstar))*Sps + 2.0*(3.0+2.0*cos(f/fstar)+cos(2.0*f/fstar))*Sacc*pow(f,-3./2.)*(1.0/pow(2.0*pi*f,4)+ red/pow(2.0*pi*f,6))) / pow(2.0*L,2.0));
		*SXYZ = 4.0/*pow(sin(f/fstar),2.0)*/*( ( 4.0*Sps + 8.0*(1.0+pow(cos(f/fstar),2.0))*Sacc*pow(f,-3./2.)*(1.0/pow(2.0*pi*f,4)+ red/pow(2.0*pi*f,6))) / pow(2.0*L,2.0));
	}
	
	if(*SAE < 1e-45) *SAE = 1e-45;
	if(*SXYZ< 1e-45) *SXYZ= 1e-45;
	
}
