// :tlittenb:20101218 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Declarations.h"
#include "Constants.h"
#include "Constants2.h"
#include "/Volumes/Apps_and_Docs/tlittenb/galaxy/codes/numrec.h"
#include "/Volumes/Apps_and_Docs/tlittenb/galaxy/codes/memory.h"


/* ********************************************************************************** */
/*                                                                                    */
/*                                Signal Generator Tools                              */
/*                                                                                    */
/* ********************************************************************************** */
int bandwidth(double params[])
{
	int     q = (int)floor(params[0]);
	double  A = exp(params[4]);
	double  f = params[0]/T;
	double sf = sin(f/fstar); //sin(f/f*)
	double sn = AEnoise(f);
	
	//Doppler spreading
	int DS = 64;
	if(q > 225000)	DS = 128;
	if(q > 629145)  DS = 512;
	if(q > 1887436) DS = 1024;
	if(q > 6291456) DS = 2048;
	
	//Sinc spreading
	double SNm  = sn/(4.*sf*sf);   //Michelson noise
	double SNRm = A*sqT/sqrt(SNm); //Michelson SNR (w/ no spread)
	
	int SS = (int)(pow(2.0,(rint(log(SNRm)/ln2)+1.0)));
	if(SS > N) SS = N;
	
	return (DS > SS) ? DS : SS; //return largest spread as bandwidth
}

double AEnoise(double f)
{
	//Power spectral density of the detector noise and transfer frequency
	double Sn, red, confusion_noise;
	double f1, f2;
	double A1, A2, slope;
	
	red = (2.0*pi*1.0e-4)*(2.0*pi*1.0e-4);
	
	confusion_noise = 0.0;
	
	// Calculate the power spectral density of the galactic background at the given frequency
	if((f > 1.0e-4) && (f <= 2.5e-4))
	{
		f1 = 1.0e-4;
		f2 =  2.5e-4;
		A1 = 1.65e-39;
		A2 = 2.0e-39;
		slope = (log10(A2)-log10(A1))/(log10(f2)-log10(f1));
		confusion_noise = A1*pow((f/f1),slope);
	}
	if((f > 2.5e-4) && (f <= 4.5e-4))
	{
		f1 = 2.5e-4;
		f2 =  4.5e-4;
		A1 = 2.0e-39;
		A2 = 1.5e-39;
		slope = (log10(A2)-log10(A1))/(log10(f2)-log10(f1));
		confusion_noise = A1*pow((f/f1),slope);
	}
	if((f > 4.5e-4) && (f <= 1.0e-3))
	{
		f1 = 4.5e-4;
		f2 =  1.0e-3;
		A1 = 1.5e-39;
		A2 = 7.0e-40;
		slope = (log10(A2)-log10(A1))/(log10(f2)-log10(f1));
		confusion_noise = A1*pow((f/f1),slope);
	}
	if((f > 1.0e-3) && (f <= 1.5e-3))
	{
		f1 = 1.0e-3;
		f2 =  1.5e-3;
		A1 = 7.0e-40;
		A2 = 3.5e-40;
		slope = (log10(A2)-log10(A1))/(log10(f2)-log10(f1));
		confusion_noise = A1*pow((f/f1),slope);
	}
	if((f > 1.5e-3) && (f <= 2.0e-3))
	{
		f1 = 1.5e-3;
		f2 =  2.0e-3;
		A1 = 3.5e-40;
		A2 = 1.5e-40;
		slope = (log10(A2)-log10(A1))/(log10(f2)-log10(f1));
		confusion_noise = A1*pow((f/f1),slope);
	}
	if((f > 2.0e-3) && (f <= 2.2e-3))
	{
		f1 = 2.0e-3;
		f2 =  2.2e-3;
		A1 = 1.5e-40;
		A2 = 9.2e-41;
		slope = (log10(A2)-log10(A1))/(log10(f2)-log10(f1));
		confusion_noise = A1*pow((f/f1),slope);
	}
	if((f > 2.2e-3) && (f <= 2.3e-3))
	{
		f1 = 2.2e-3;
		f2 =  2.3e-3;
		A1 = 9.2e-41;
		A2 = 4.1e-41;
		slope = (log10(A2)-log10(A1))/(log10(f2)-log10(f1));
		confusion_noise = A1*pow((f/f1),slope);
	}
	if((f > 2.3e-3) && (f <= 2.5e-3))
	{
		f1 = 2.3e-3;
		f2 =  2.5e-3;
		A1 = 4.1e-41;
		A2 = 9.1e-42;
		slope = (log10(A2)-log10(A1))/(log10(f2)-log10(f1));
		confusion_noise = A1*pow((f/f1),slope);
	}
	if((f > 2.5e-3))
	{
		f1 = 2.5e-3;
		f2 =  3.0e-3;
		A1 = 9.1e-42;
		A2 = 5.5e-43;
		slope = (log10(A2)-log10(A1))/(log10(f2)-log10(f1));
		confusion_noise = A1*pow((f/f1),slope);
	}
	
	// Calculate the power spectral density of the detector noise at the given frequency
	Sn = 16.0/3.0*pow(sin(f/fstar),2.0)*( ( (2.0+cos(f/fstar))*Sps + 2.0*(3.0+2.0*cos(f/fstar)+cos(2.0*f/fstar))*Sacc*(1.0/pow(2.0*pi*f,4)+ red/pow(2.0*pi*f,6))) / pow(2.0*L,2.0)) + confusion_noise;

	return Sn;
}

double overlap(double **hi, double **hj, double **Sn)
{
	int i;
	
	double hihi = 0.0;
	double hjhj = 0.0;
	double hihj = 0.0;
	
	for(i=0; i<NI; i++)
	{
		hihi += fourier_nwip(hi[i], hi[i], Sn[i], N);
		hjhj += fourier_nwip(hj[i], hj[i], Sn[i], N);
		hihj += fourier_nwip(hi[i], hj[i], Sn[i], N);
	}
	
	return hihj/sqrt(hihi*hjhj);
}

void waveforms(struct lisa_orbit *orbit, double *params, int b0, double **h)
{
	int i;
	int j;
	int k;
	int l;
	int m;
	int n;
	int p;
	int bi;
	int db;
	int bj;
	
	double *A;
	double *E;
	double P[NP];
	
	//copy signal parameters
	for(i=0; i<NP; i++) P[i] = params[i];
	
	/*****************************************/
	/* MLDC waveforms don't account for fdot */
	/*****************************************/	
	//Dynamically decide bandwidth
	int BW = bandwidth(P);
	int BW2  = 2*BW;
	int BWon2= BW/2;
	
	A = dvector(1,BW2);
	E = dvector(1,BW2);
	
	//zero waveform arrays
	for(i=0; i<N2;  i++) h[0][i] = h[1][i] = 0.0;
	for(i=1; i<=BW2; i++) A[i] = E[i] = 0.0;
	
	/*   Determine where in the datastream the about-to-be-generated template belongs   */ 
	bi = (long)P[0] - BWon2;
	
	if(b0 > bi) db = b0 - bi;
	else db = 0;
	
	/*   Wave-form generator produces templates A & E   */ 
	galactic_binary(orbit, P, A, E, BW); 
	
	bj = db;  
	for(i=1; i<=N; i++)
	{
		k = i-1;
		j = i-1 + b0;		
		
		l = 2*k;
		m = l+1;
		
		if(j>=bi && j<bi + BW) 
		{
			bj++;
			p = 2*bj;
			n = p-1;
			
			h[0][l]	= A[n];
			h[0][m] = A[p];
			h[1][l]	= E[n];
			h[1][m] = E[p];
		}     
	}
	free_dvector(A,1,BW2);
	free_dvector(E,1,BW2);	
}	

void galactic_binary(struct lisa_orbit *orbit, double params[], double *A, double *E, int BW)
{
	
	/*   Indicies   */
	int i,j,n;
	/*   Carrier frequency bin  */
	long q;
	/*   Bandwidth      */
	int BW2   = BW*2;
	int BWon2 = BW/2;	
	/*   Gravitational Wave basis vectors   */
	double u[4],v[4],k[4];
	/*   Polarization basis tensors   */
	double eplus[4][4], ecross[4][4];
	/*   Spacecraft position and separation vector   */
	double *x, *y, *z;
	double r12[4],r13[4],r21[4],r23[4],r31[4],r32[4];	
	/*   Dot products   */
	double kdotx[4],kdotr[4][4];
	/*   Convenient quantities   */
	double dplus[4][4],dcross[4][4];
	/*   GW source parameters   */
	double phi, psi, Amp, Aplus, Across, f0, fdot, phio, fd2;
	double costh, sinth, cosph, sinph, cosi, cosps, sinps;
	/*   Time and distance variables   */
	double t, xi[4];
	/*   Gravitational wave frequency & ratio of f and transfer frequency f*  */
	double f[4], fonfs[4];
	/*   LISA response to slow terms (Real & Imaginary pieces)   */
	//Static quantities (Re and Im)
	double DPr, DPi, DCr, DCi;
	//Time varrying quantities (Re & Im) broken up into convenient segments
	double tran1r, tran1i, tran2r, tran2i;
	double TR[4][4], TI[4][4];
	//Miscellaneous constants used to speed up calculations
	double arg1, arg2, sinc, df;
	/*   Fourier coefficients before FFT and after convolution  */
	//Time series of slowly evolving terms at each vertex   
	double *data12, *data13, *data21, *data23, *data31, *data32;
	//Fourier coefficients of slowly evolving terms (numerical)
	double a12[BW2+3], a13[BW2+3], a21[BW2+3], a23[BW2+3], a31[BW2+3], a32[BW2+3];
	//Package cij's into proper form for TDI subroutines
	double ***d;
	/*   Miscellaneous  */
	double aevol;
	
	/*   Allocating Arrays   */
	x = dvector(1,3); y = dvector(1,3); z = dvector(1,3);
	
	data12 = dvector(1,BW2); data21 = dvector(1,BW2); data31 = dvector(1,BW2);
	data13 = dvector(1,BW2); data23 = dvector(1,BW2); data32 = dvector(1,BW2); 
	
	d = d3tensor(1,3,1,3,1,BW2);
	
	/*   Gravitational Wave source parameters   */	
	//Calculate carrier frequency bin
	q = (long)(params[0]);
	
	//Intrinsic
	f0		= params[0]/T;
	fdot	= params[1]/Tsq;
	costh	= params[2];
	phi		= params[3];
	//Extrinisic
	Amp		= exp(params[4]);
	cosi	= params[5];
	psi		= params[6]*2.0;
	phio	= params[7];
	//second time derivative of f
	//	if(NP>8) fd2 = params[8]/T/Tsq;
	
	//Calculate cos and sin of sky position, inclination, polarization
	sinth	= sqrt(1.0 - costh*costh); //sin(theta) >= 0 (theta -> 0,pi)
	cosph	= cos(phi);     
	sinph	= sin(phi);
	cosps	= cos(psi);  
	sinps	= sin(psi);
	
	//Calculate GW polarization amplitudes
	Aplus  =  Amp*(1.+cosi*cosi);
	Across = -Amp*(2.0*cosi);
	
	df = pi2*(f0 - ((double)q)/T);											
	
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
	for(i=1;i<=3;i++)
	{
		for(j=1;j<=3;j++)
		{
			eplus[i][j]  = u[i]*u[j] - v[i]*v[j];
			ecross[i][j] = u[i]*v[j] + v[i]*u[j];
		}
	}
	
	
	/*****************************   Main Loop   **********************************/
	for(n=1; n<=BW; n++)
	{
		//First time sample must be at t=0 for phasing
		t = T*(double)(n-1)/(double)BW;
		
		//Calculate position of each spacecraft at time t
		spacecraft(orbit, t, x, y, z);
		
		for(i=1; i<=3; i++) 
		{
			kdotx[i] = (x[i]*k[1]+y[i]*k[2]+z[i]*k[3])/c;
			
			//Wave arrival time at spacecraft i
			xi[i] = t - kdotx[i];
			
			//First order approximation to frequency at spacecraft i
			f[i] = f0 + fdot*xi[i];
			
			//Second order in frequency
			//if(NP>8) f[i] += 0.5*fd2*xi[i]*xi[i];
			
			//Ratio of true frequency to transfer frequency
			fonfs[i] = f[i]/fstar;
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
			for(j=1; j<=3; j++)
			{
				if(i!=j)
				{
					//Argument of transfer function
					arg1 = 0.5*fonfs[i]*(1.0 - kdotr[i][j]);
					
					//Argument of complex exponentials
					arg2 = df*xi[i] + pi*fdot*xi[i]*xi[i] + phio - pi2*kdotx[i]*f0;
					
					//Second order frequency evolution
					//if(NP>8) arg2 += pi/3.0*fd2*xi[i]*xi[i]*x[i];
					
					//Transfer function
					sinc = 0.25*sin(arg1)/arg1;
					
					//Evolution of amplitude
					aevol = 1.0 + 0.66666666666666666666*fdot/f0*xi[i]; 
					
					//Second order amplitude evolution
					//if(NP>8) aevol += const.*fd2*xi[i]*xi[i]/f0;
					
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
	dfour1(data12, BW, -1);  dfour1(data21, BW, -1);  dfour1(data31, BW, -1);
	dfour1(data13, BW, -1);  dfour1(data23, BW, -1);  dfour1(data32, BW, -1);  
	
	//Unpack arrays from dfour1.c and normalize
	for(i=1; i<=BW; i++)
	{
		j = i + BW;
		a12[i] = data12[j]/(double)BW2;  a21[i] = data21[j]/(double)BW2;  a31[i] = data31[j]/(double)BW2;
		a12[j] = data12[i]/(double)BW2;  a21[j] = data21[i]/(double)BW2;  a31[j] = data31[i]/(double)BW2;
		a13[i] = data13[j]/(double)BW2;  a23[i] = data23[j]/(double)BW2;  a32[i] = data32[j]/(double)BW2;
		a13[j] = data13[i]/(double)BW2;  a23[j] = data23[i]/(double)BW2;  a32[j] = data32[i]/(double)BW2;
	}
	
	/*   Renormalize so that the resulting time series is real   */
	for(i=1; i<=BW2; i++)
	{
		d[1][2][i] = a12[i];  d[2][1][i] = a21[i];  d[3][1][i] = a31[i];
		d[1][3][i] = a13[i];  d[2][3][i] = a23[i];  d[3][2][i] = a32[i];
	}
	
	/*   Call subroutines for synthesizing different TDI data channels  */
	tdi(d, f0, q, A, E, BW);
	
	/*   Deallocate Arrays   */
	free_dvector(x,1,3); free_dvector(y,1,3); free_dvector(z,1,3);
	free_dvector(data12,1,BW2); free_dvector(data21,1,BW2); free_dvector(data31,1,BW2);
	free_dvector(data13,1,BW2); free_dvector(data23,1,BW2); free_dvector(data32,1,BW2); 
	free_d3tensor(d,1,3,1,3,1,BW2);
	
	return;
}


void spacecraft(struct lisa_orbit *orbit, double tint, double *xint, double *yint, double *zint)
{
	int i;
	for(i=0; i<3; i++) 
	{
		LISA_splint(orbit->t, orbit->x[i], orbit->dx[i], orbit->Norb, tint, &xint[i+1]);
		LISA_splint(orbit->t, orbit->y[i], orbit->dy[i], orbit->Norb, tint, &yint[i+1]);
		LISA_splint(orbit->t, orbit->z[i], orbit->dz[i], orbit->Norb, tint, &zint[i+1]);
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
	int Norb = n;
	double *t  = dvector(0,Norb-1);
	double **x  = dmatrix(0,2,0,Norb-1);
	double **y  = dmatrix(0,2,0,Norb-1);
	double **z  = dmatrix(0,2,0,Norb-1);
	double **dx = dmatrix(0,2,0,Norb-1);
	double **dy = dmatrix(0,2,0,Norb-1);
	double **dz = dmatrix(0,2,0,Norb-1);
	
	//allocate memory for orbit structure
	orbit->Norb = n;
	
	orbit->t  = dvector(0,orbit->Norb-1);
	orbit->x  = dmatrix(0,2,0,orbit->Norb-1);
	orbit->y  = dmatrix(0,2,0,orbit->Norb-1);
	orbit->z  = dmatrix(0,2,0,orbit->Norb-1);
	orbit->dx = dmatrix(0,2,0,orbit->Norb-1);
	orbit->dy = dmatrix(0,2,0,orbit->Norb-1);
	orbit->dz = dmatrix(0,2,0,orbit->Norb-1);
	
	//read in orbits
	for(n=0; n<orbit->Norb; n++)
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
		LISA_spline(t, orbit->x[i], Norb, 1.e30, 1.e30, orbit->dx[i]);
		LISA_spline(t, orbit->y[i], Norb, 1.e30, 1.e30, orbit->dy[i]);
		LISA_spline(t, orbit->z[i], Norb, 1.e30, 1.e30, orbit->dz[i]);
	}
	
	//free local memory
	free_dvector(t,0,Norb-1);
	free_dmatrix(x ,0,2,0,Norb-1);
	free_dmatrix(y ,0,2,0,Norb-1);
	free_dmatrix(z ,0,2,0,Norb-1);
	free_dmatrix(dx,0,2,0,Norb-1);
	free_dmatrix(dy,0,2,0,Norb-1);
	free_dmatrix(dz,0,2,0,Norb-1);
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



void tdi(double ***d, double f0, long q, double *A, double *E, int BW)
{ 
	int i,j,k; 
	int BW2   = BW*2;
	int BWon2 = BW/2;
	double fonfs;
	double c3, s3, c2, s2, c1, s1;
	double f;
	double X[BW2+1],Y[BW2+1],Z[BW2+1];
	double phiLS, cLS, sLS;
	
	//phiLS = pi2*f0*(0.5-LonC);//1 s sampling rate											
	phiLS = pi2*f0*(7.5-LonC);//15 s sampling rate
	cLS = cos(phiLS);
	sLS = sin(phiLS);
	
	for(i=1; i<=BW; i++)
	{
		k = 2*i;
		j = k-1;
		
		f = ((double)(q + i-1 - BWon2))/T;
		fonfs = f/fstar;
		c3 = cos(3.*fonfs);  c2 = cos(2.*fonfs);  c1 = cos(fonfs); 
		s3 = sin(3.*fonfs);  s2 = sin(2.*fonfs);  s1 = sin(fonfs);
		
		X[j] =	(d[1][2][j]-d[1][3][j])*c3 + (d[1][2][k]-d[1][3][k])*s3 +
		(d[2][1][j]-d[3][1][j])*c2 + (d[2][1][k]-d[3][1][k])*s2 +
		(d[1][3][j]-d[1][2][j])*c1 + (d[1][3][k]-d[1][2][k])*s1 +
		(d[3][1][j]-d[2][1][j]);
		
		X[k] =	(d[1][2][k]-d[1][3][k])*c3 - (d[1][2][j]-d[1][3][j])*s3 +
		(d[2][1][k]-d[3][1][k])*c2 - (d[2][1][j]-d[3][1][j])*s2 +
		(d[1][3][k]-d[1][2][k])*c1 - (d[1][3][j]-d[1][2][j])*s1 +
		(d[3][1][k]-d[2][1][k]);
		
		Y[j] =	(d[2][3][j]-d[2][1][j])*c3 + (d[2][3][k]-d[2][1][k])*s3 +
		(d[3][2][j]-d[1][2][j])*c2 + (d[3][2][k]-d[1][2][k])*s2+
		(d[2][1][j]-d[2][3][j])*c1 + (d[2][1][k]-d[2][3][k])*s1+
		(d[1][2][j]-d[3][2][j]);
		
		Y[k] =	(d[2][3][k]-d[2][1][k])*c3 - (d[2][3][j]-d[2][1][j])*s3+
		(d[3][2][k]-d[1][2][k])*c2 - (d[3][2][j]-d[1][2][j])*s2+
		(d[2][1][k]-d[2][3][k])*c1 - (d[2][1][j]-d[2][3][j])*s1+
		(d[1][2][k]-d[3][2][k]);
		
		Z[j] =	(d[3][1][j]-d[3][2][j])*c3 + (d[3][1][k]-d[3][2][k])*s3+
		(d[1][3][j]-d[2][3][j])*c2 + (d[1][3][k]-d[2][3][k])*s2+
		(d[3][2][j]-d[3][1][j])*c1 + (d[3][2][k]-d[3][1][k])*s1+
		(d[2][3][j]-d[1][3][j]);
		
		Z[k] =	(d[3][1][k]-d[3][2][k])*c3 - (d[3][1][j]-d[3][2][j])*s3+
		(d[1][3][k]-d[2][3][k])*c2 - (d[1][3][j]-d[2][3][j])*s2+
		(d[3][2][k]-d[3][1][k])*c1 - (d[3][2][j]-d[3][1][j])*s1+
		(d[2][3][k]-d[1][3][k]);
		
		A[j] =  sqT*((2.0*X[j]-Y[j]-Z[j])*cLS-(2.0*X[k]-Y[k]-Z[k])*sLS)/3.0;
		A[k] = -sqT*((2.0*X[j]-Y[j]-Z[j])*sLS+(2.0*X[k]-Y[k]-Z[k])*cLS)/3.0;
		
		E[j] =  sqT*((Z[j]-Y[j])*cLS-(Z[k]-Y[k])*sLS)/sq3;
		E[k] = -sqT*((Z[j]-Y[j])*sLS+(Z[k]-Y[k])*cLS)/sq3;
	}
	
}

/* ********************************************************************************** */
/*                                                                                    */
/*                                     Math Tools                                     */
/*                                                                                    */
/* ********************************************************************************** */

void phase_blind_time_shift(double *corr, double *corrf, double *data1, double *data2, double *psd)
{
	int i, l, k;
	
	for(i=0; i<N; i++)
	{
		l=2*i;
		k=l+1;
		
		corr[l]		= ( data1[l]*data2[l] + data1[k]*data2[k]) / psd[i];
		corr[k]		= ( data1[k]*data2[l] - data1[l]*data2[k]) / psd[i];
		corrf[l]	= ( data1[l]*data2[k] - data1[k]*data2[l]) / psd[i];
		corrf[k]	= ( data1[k]*data2[k] + data1[l]*data2[l]) / psd[i];
	}
	
	drealft(corr-1, N2, -1);
	drealft(corrf-1, N2, -1);
}

void max_array_element(double *max, int *index, double *array)
{
	int i;
	
	*max	= array[0];
	*index	= 0;
	
	for(i = 1; i < N2; i++)
	{
		if(array[i] > *max)
		{
			*max = array[i];
			*index = i;
		}
	}
}

void jacobi(double **a, int n, double *d, double **v, int *nrot)
{
	int j,iq,ip,i;
	double tresh,theta,tau,t,sm,s,h,g,cc,*b,*z;
	b=dvector(0,n-1);
	z=dvector(0,n-1);
	for (ip=0;ip<n;ip++) {
		for (iq=0;iq<n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=0;ip<n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	*nrot=0;
	for (i=0;i<100;i++) {
		sm=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0) {
			free_dvector(z,0,n-1);
			free_dvector(b,0,n-1);
			return;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
					&& (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((double)(fabs(h)+g) == (double)fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					cc=1.0/sqrt(1+t*t);
					s=t*cc;
					tau=s/(1.0+cc);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=0;j<ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=0;j<n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}
		for (ip=0;ip<n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	printf("  Too many iterations in routine jacobi\n");
	*nrot=-1;
	free_dvector(z,0,n-1);
	free_dvector(b,0,n-1);
	return;
}

void indexx(unsigned long n, double *arr, unsigned long *indx)
{
	unsigned long i,indxt,ir=n,j,k,l=1;
	int jstack=0,*istack,M=7;
	double a;
	
	istack=ivector(1,NSTACK);
	for (j=1;j<=n;j++) indx[j]=j;
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				indxt=indx[j];
				a=arr[indxt];
				for (i=j-1;i>=l;i--) {
					if (arr[indx[i]] <= a) break;
					indx[i+1]=indx[i];
				}
				indx[i+1]=indxt;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(indx[k],indx[l+1]);
			if (arr[indx[l]] > arr[indx[ir]]) {
				SWAP(indx[l],indx[ir])
			}
			if (arr[indx[l+1]] > arr[indx[ir]]) {
				SWAP(indx[l+1],indx[ir])
			}
			if (arr[indx[l]] > arr[indx[l+1]]) {
				SWAP(indx[l],indx[l+1])
			}
			i=l+1;
			j=ir;
			indxt=indx[l+1];
			a=arr[indxt];
			for (;;) {
				do i++; while (arr[indx[i]] < a);
				do j--; while (arr[indx[j]] > a);
				if (j < i) break;
				SWAP(indx[i],indx[j])
			}
			indx[l+1]=indx[j];
			indx[j]=indxt;
			jstack += 2;
			//if (jstack > NSTACK) nrerror("NSTACK too small in indexx.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free_ivector(istack,1,NSTACK);
}

void gaussj(double **a, int n, double **b, int m)
{
	int *indxc,*indxr,*ipiv;
	int i,icol=0,irow=0,j,k,l,ll;
	double big,dum,pivinv;
	
	indxc=ivector(1,n);
	indxr=ivector(1,n);
	ipiv=ivector(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					}
				}
		++(ipiv[icol]);
		if (irow != icol) { 
			for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
				for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
					}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) printf("gaussj: Singular Matrix\n"); 
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *= pivinv;
		for (l=1;l<=m;l++) b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n;l>=1;l--) {
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
}

double sign(double x)
{
	if(x > 0.0) 
	{
		return 1.0;
	}
	else 
	{
		if (x < 0) 
		{
			return -1.0;
		}
		else
		{
			return 0.0;
		}
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

double gauss(double x, double mean, double sigma)
{
	return -0.5*(x - mean)*(x - mean)/(sigma*sigma);
}

double fourier_nwip(double *a, double *b, double *Sn, int n)
{
	int i, j, k;
	double arg, product;
	double ReA, ReB, ImA, ImB;
	
	arg = 0.0;
	for(i=0; i<n; i++) 
	{
		j = i * 2;
		k = j + 1;
		ReA = a[j]; ImA = a[k];
		ReB = b[j]; ImB = b[k];
		product = ReA*ReB + ImA*ImB;      
		arg += product/Sn[i];
	}
	
	return(4.0*arg);
}

/* ********************************************************************************** */
/*																					  */
/*                                   Fourier Tools                                    */
/*																					  */
/* ********************************************************************************** */

void dfour1(double data[], unsigned long nn, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	double tempr,tempi, swap;
	
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

void drealft(double data[], unsigned long n, int isign)
{
	void dfour1(double data[], unsigned long nn, int isign);
	unsigned long i,i1,i2,i3,i4,np3;
	double c1=0.5,c2,h1r,h1i,h2r,h2i;
	double wr,wi,wpr,wpi,wtemp,theta;
	
	theta=3.141592653589793/(double) (n>>1);
	if (isign == 1) {
		c2 = -0.5;
		dfour1(data,n>>1,1);
	} else {
		c2=0.5;
		theta = -theta;
	}
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	np3=n+3;
	for (i=2;i<=(n>>2);i++) {
		i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
		h1r=c1*(data[i1]+data[i3]);
		h1i=c1*(data[i2]-data[i4]);
		h2r = -c2*(data[i2]+data[i4]);
		h2i=c2*(data[i1]-data[i3]);
		data[i1]=h1r+wr*h2r-wi*h2i;
		data[i2]=h1i+wr*h2i+wi*h2r;
		data[i3]=h1r-wr*h2r+wi*h2i;
		data[i4] = -h1i+wr*h2i+wi*h2r;
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1) {
		data[1] = (h1r=data[1])+data[2];
		data[2] = h1r-data[2];
	} else {
		data[1]=c1*((h1r=data[1])+data[2]);
		data[2]=c1*(h1r-data[2]);
		dfour1(data,n>>1,-1);
	}
}

double power_spectrum(double *data, int n)
{
	int i,j;
	double Re, Im;
	
	i = 2*n;
	j = i+1;
	Re = data[i];
	Im = data[j];
	
	return (Re*Re + Im*Im);
}

/* ********************************************************************************** */
/*                                                                                    */
/*                           Memory (de)allocation routines                           */
/*                                                                                    */
/* ********************************************************************************** */
int file_exist(char *filename)
{
	FILE *file;
    if(file = fopen(filename, "r"))
    {
        fclose(file);
        return 1;
    }
    return 0;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;
	
	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) fprintf(stderr,"allocation failure in ivector()");
	return v-nl+NR_END;
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v=0;
	
	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) fprintf(stderr,"allocation failure in dvector()");
	return v-nl+NR_END;
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v;
	
	v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) fprintf(stderr,"allocation failure in lvector()");
	return v-nl+NR_END;
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;
	
	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) fprintf(stderr, "allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;
	
	
	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m) fprintf(stderr, "allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
	/* return pointer to array of pointers to rows */
	return m;
}

unsigned long **lmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	unsigned long **m;
	
	/* allocate pointers to rows */
	m=(unsigned long **) malloc((size_t)((nrow+NR_END)*sizeof(long*)));
	if (!m) fprintf(stderr, "allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;
	
	
	/* allocate rows and set pointers to them */
	m[nrl]=(unsigned long *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(long)));
	if (!m) fprintf(stderr, "allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
	/* return pointer to array of pointers to rows */
	return m;
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_lmatrix(unsigned long **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;
	
	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) fprintf(stderr, "allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;
	
	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) fprintf(stderr,"allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
	/* return pointer to array of pointers to rows */
	return m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	double ***t;
	
	/* allocate pointers to pointers to rows */
	t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
	if (!t) fprintf(stderr,"allocation failure 1 in d3tensor()");
	t += NR_END;
	t -= nrl;
	
	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
	if (!t[nrl]) fprintf(stderr,"allocation failure 2 in d3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;
	
	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
	if (!t[nrl][ncl]) fprintf(stderr,"allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;
	
	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}
	
	/* return pointer to array of pointers to rows */
	return t;
}


void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* free a float d3tensor allocated by d3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

//int ***i3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
///* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
//{
//	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
//	int ***t;
//	
//	/* allocate pointers to pointers to rows */
//	t=(int ***) malloc((size_t)((nrow+NR_END)*sizeof(int**)));
//	if (!t) fprintf(stderr,"allocation failure 1 in d3tensor()");
//	t += NR_END;
//	t -= nrl;
//	
//	/* allocate pointers to rows and set pointers to them */
//	t[nrl]=(int **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int*)));
//	if (!t[nrl]) fprintf(stderr,"allocation failure 2 in d3tensor()");
//	t[nrl] += NR_END;
//	t[nrl] -= ncl;
//	
//	/* allocate rows and set pointers to them */
//	t[nrl][ncl]=(int *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(int)));
//	if (!t[nrl][ncl]) fprintf(stderr,"allocation failure 3 in f3tensor()");
//	t[nrl][ncl] += NR_END;
//	t[nrl][ncl] -= ndl;
//	
//	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
//	for(i=nrl+1;i<=nrh;i++) {
//		t[i]=t[i-1]+ncol;
//		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
//		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
//	}
//	
//	/* return pointer to array of pointers to rows */
//	return t;
//}
//
//void free_i3tensor(int ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
///* free an int d3tensor allocated by i3tensor() */
//{
//	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
//	free((FREE_ARG) (t[nrl]+ncl-NR_END));
//	free((FREE_ARG) (t+nrl-NR_END));
//}
//
//
//




