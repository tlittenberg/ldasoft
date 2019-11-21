/* gcc -O2 -o Catalogue Catalogue.c -lm */
/* 04/28/2011 -- convert GC to EC -- TBL */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define pi 3.141592653589793
#define pc 3.0856775807e16
#define clight 299792458.
#define Msun 1.9889e30
#define G 6.67259e-11
#define YEAR 31457280.0


double ran2(long *idum);

int main(int argc,char **argv)
{
  int i;
  long seed = -987654876;
  double junk;
  
  //Nelemns parameters
  double m1,m2,Porb,Porb_dot,m2_dot,l,b,DL,ecc;
  
  //Conversion parameters
  double Mchirp;
  
  //LISA parameters
  double f,fdot,theta,phi,Amp,iota,psi,phi0;
  
  double deg2rad = pi/180.0;
  double gc = 8500.0;
  double kpc2pc =  1000.0;
  double x,y,z;
  double xgc,ygc,zgc;
  double xec,yec,zec,rec;
  
  FILE *Infile;
  FILE *Outfile;
  FILE *Outfile2;
  
  int n,NS,index;
  
  /************************************************/
  /*              aCOSMIC Systems                 */
  /************************************************/
  
  printf("reading full_galaxy.dat...\n");
  Outfile = fopen("full_galaxy_GW.dat","w");
  Outfile2 = fopen("full_galaxy_Gal.dat","w");
  Infile = fopen("full_galaxy.dat","r");
  
  NS = 30000000;
  char header[128];
  fgets(header,128,Infile);
  for(n=0; n<NS; n++)
  {
    if(n%(NS/100)==0)printf("Converting Binaries: %i/%i\n",n,NS);
    
    /* COSMIC files */
    // index mass1 mass2 porb ecc xGx yGx zGx dist inc OMEGA omega
    fscanf(Infile, "%i%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg\n", &index, &m1, &m2, &Porb, &ecc, &xgc, &ygc, &zgc, &rec, &iota,&psi,&phi0);
    
    //convert orbital parameters to GW parameters
    f    =  2.0/(Porb*YEAR);
    fdot = 0.0;//-2.0*Porb_dot/Porb/Porb;
    
    if(f>.0004)
    {
      
      //convert galactic cartesian coordinates to ecliptic coordinates
      x = xgc+gc;
      y = ygc;
      z = zgc;
      
      fprintf(Outfile2,"%lg %lg %lg %lg\n",Porb,x,y,z);
      
      xec = -0.05487556043*x + 0.4941094278*y - 0.86766614920*z;
      yec = -0.99382137890*x - 0.1109907351*y - 0.00035159077*z;
      zec = -0.09647662818*x + 0.8622858751*y + 0.49714719180*z;
      
      DL = sqrt(xec*xec + yec*yec + zec*zec);
      
      theta = pi/2.0 - acos(zec/DL);
      phi   = atan2(yec,xec);
      while(phi<0.0) phi += 2.0*pi;
      
      
      //convert mass+distance to amplitude
      m1 *= Msun;
      m2 *= Msun;
      Amp    = 2.0 * G*G*m1*m2 * pow(pi*pi*f*f/(G*(m1+m2)),1./3.) / (DL*pc) /(clight*clight*clight*clight);
      
      fprintf(Outfile, "%.16g %.10g %.6g %.6g %.6g %.6g %.6g %.6g\n",f,fdot,theta,phi,Amp,iota,psi,phi0);
    }
  }
  fclose(Infile);
  fclose(Outfile);
  
  
  return 0;
  
}

/************************************************/
/*              Detached Systems                */
/************************************************/
//
//	printf("reading dwdSourceCatalogue.txt.1.1...\n");
//	Outfile = fopen("dwdSourceCatalogue.txt.1.1","w");
//	Infile = fopen("dwd_GWR_NewLISA.dat.1.1","r");
//	NS = 6521099;
//	for(n=0; n<NS; n++)
//	{
//		if(n%1000000==0)printf("Converting DWD: %i/%i\n",n,NS);
//
//		/*Nelemens files*/
//		fscanf(Infile, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%i\n", &m1, &m2, &Porb, &Porb_dot, &m2_dot, &l, &b, &DL, &junk, &junk, &junk, &i);
//
//		//convert galactic coordinates to radians
//		l = l*deg2rad;
//		b = (90.0 - b)*deg2rad;
//
//		//convert distance to pc
//		DL *= kpc2pc;
//
//		//change to ecliptic coordinates
//		x = DL*sin(b)*cos(l);
//		y = DL*sin(b)*sin(l);
//		z = DL*cos(b);
//
//		xgc = x-gc;
//		ygc = y;
//		zgc = z;
//
//		xec = -0.05487556043*x + 0.4941094278*y - 0.86766614920*z;
//		yec = -0.99382137890*x - 0.1109907351*y - 0.00035159077*z;
//		zec = -0.09647662818*x + 0.8622858751*y + 0.49714719180*z;
//
//		rec = sqrt(xec*xec + yec*yec + zec*zec);
//
//		theta = pi/2.0 - acos(zec/rec);
//		phi   = atan2(yec,xec);
//		while(phi<0.0) phi += 2.0*pi;
//
//		//convert orbital parameters to GW parameters
//		f    =  2.0/Porb;
//		fdot = -2.0*Porb_dot/Porb/Porb;
//
//		//convert mass+distance to amplitude
//		m1 *= Msun;
//		m2 *= Msun;
//		Amp    = 2.0 * G*G*m1*m2 * pow(pi*pi*f*f/(G*(m1+m2)),1./3.) / (DL*pc) /(clight*clight*clight*clight);
//
//		//randomly assign orientation parameters
//		iota = acos(-1.0 + 2.0*ran2(&seed));
//		psi  = ran2(&seed)*2.0*pi;
//		phi0 = ran2(&seed)*2.0*pi;
//
//		fprintf(Outfile, "%.16g %.10g %.6g %.6g %.6g %.6g %.6g %.6g\n",f,fdot,theta,phi,Amp,iota,psi,phi0);
//	}
//	fclose(Infile);
//	fclose(Outfile);
//
//	printf("reading dwdSourceCatalogue.txt.1.2...\n");
//	Outfile = fopen("dwdSourceCatalogue.txt.1.2","w");
//	Infile = fopen("dwd_GWR_NewLISA.dat.1.2","r");
//	NS = 6521098;
//	for(n=0; n<NS; n++)
//	{
//		if(n%1000000==0)printf("Converting DWD: %i/%i\n",n,NS);
//
//		/*Nelemens files*/
//		fscanf(Infile, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%i\n", &m1, &m2, &Porb, &Porb_dot, &m2_dot, &l, &b, &DL, &junk, &junk, &junk, &i);
//
//		//convert galactic coordinates to radians
//		l = l*deg2rad;
//		b = (90.0 - b)*deg2rad;
//
//		//convert distance to pc
//		DL *= kpc2pc;
//
//		//change to ecliptic coordinates
//		x = DL*sin(b)*cos(l);
//		y = DL*sin(b)*sin(l);
//		z = DL*cos(b);
//
//		xgc = x-gc;
//		ygc = y;
//		zgc = z;
//
//		xec = -0.05487556043*x + 0.4941094278*y - 0.86766614920*z;
//		yec = -0.99382137890*x - 0.1109907351*y - 0.00035159077*z;
//		zec = -0.09647662818*x + 0.8622858751*y + 0.49714719180*z;
//
//		rec = sqrt(xec*xec + yec*yec + zec*zec);
//
//		theta = pi/2.0 - acos(zec/rec);
//		phi   = atan2(yec,xec);
//		while(phi<0.0) phi += 2.0*pi;
//
//		//convert orbital parameters to GW parameters
//		f    =  2.0/Porb;
//		fdot = -2.0*Porb_dot/Porb/Porb;
//
//		//convert mass+distance to amplitude
//		m1 *= Msun;
//		m2 *= Msun;
//		Amp    = 2.0 * G*G*m1*m2 * pow(pi*pi*f*f/(G*(m1+m2)),1./3.) / (DL*pc) /(clight*clight*clight*clight);
//
//		//randomly assign orientation parameters
//		iota = acos(-1.0 + 2.0*ran2(&seed));
//		psi  = ran2(&seed)*2.0*pi;
//		phi0 = ran2(&seed)*2.0*pi;
//
//		fprintf(Outfile, "%.16g %.10g %.6g %.6g %.6g %.6g %.6g %.6g\n",f,fdot,theta,phi,Amp,iota,psi,phi0);
//	}
//	fclose(Infile);
//	fclose(Outfile);
//
//	printf("reading dwdSourceCatalogue.txt.2.1...\n");
//	Outfile = fopen("dwdSourceCatalogue.txt.2.1","w");
//	Infile = fopen("dwd_GWR_NewLISA.dat.2.1","r");
//	NS = 6521098;
//	for(n=0; n<NS; n++)
//	{
//		if(n%1000000==0)printf("Converting DWD: %i/%i\n",n,NS);
//
//		/*Nelemens files*/
//		fscanf(Infile, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%i\n", &m1, &m2, &Porb, &Porb_dot, &m2_dot, &l, &b, &DL, &junk, &junk, &junk, &i);
//
//		//convert galactic coordinates to radians
//		l = l*deg2rad;
//		b = (90.0 - b)*deg2rad;
//
//		//convert distance to pc
//		DL *= kpc2pc;
//
//		//change to ecliptic coordinates
//		x = DL*sin(b)*cos(l);
//		y = DL*sin(b)*sin(l);
//		z = DL*cos(b);
//
//		xgc = x-gc;
//		ygc = y;
//		zgc = z;
//
//		xec = -0.05487556043*x + 0.4941094278*y - 0.86766614920*z;
//		yec = -0.99382137890*x - 0.1109907351*y - 0.00035159077*z;
//		zec = -0.09647662818*x + 0.8622858751*y + 0.49714719180*z;
//
//		rec = sqrt(xec*xec + yec*yec + zec*zec);
//
//		theta = pi/2.0 - acos(zec/rec);
//		phi   = atan2(yec,xec);
//		while(phi<0.0) phi += 2.0*pi;
//
//		//convert orbital parameters to GW parameters
//		f    =  2.0/Porb;
//		fdot = -2.0*Porb_dot/Porb/Porb;
//
//		//convert mass+distance to amplitude
//		m1 *= Msun;
//		m2 *= Msun;
//		Amp    = 2.0 * G*G*m1*m2 * pow(pi*pi*f*f/(G*(m1+m2)),1./3.) / (DL*pc) /(clight*clight*clight*clight);
//
//		//randomly assign orientation parameters
//		iota = acos(-1.0 + 2.0*ran2(&seed));
//		psi  = ran2(&seed)*2.0*pi;
//		phi0 = ran2(&seed)*2.0*pi;
//
//		fprintf(Outfile, "%.16g %.10g %.6g %.6g %.6g %.6g %.6g %.6g\n",f,fdot,theta,phi,Amp,iota,psi,phi0);
//	}
//	fclose(Infile);
//	fclose(Outfile);
//
//	printf("reading dwdSourceCatalogue.txt.2.2...\n");
//	Outfile = fopen("dwdSourceCatalogue.txt.2.2","w");
//	Infile = fopen("dwd_GWR_NewLISA.dat.2.2","r");
//	NS = 6521098;
//	for(n=0; n<NS; n++)
//	{
//		if(n%1000000==0)printf("Converting DWD: %i/%i\n",n,NS);
//
//		/*Nelemens files*/
//		fscanf(Infile, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%i\n", &m1, &m2, &Porb, &Porb_dot, &m2_dot, &l, &b, &DL, &junk, &junk, &junk, &i);
//
//		//convert galactic coordinates to radians
//		l = l*deg2rad;
//		b = (90.0 - b)*deg2rad;
//
//		//convert distance to pc
//		DL *= kpc2pc;
//
//		//change to ecliptic coordinates
//		x = DL*sin(b)*cos(l);
//		y = DL*sin(b)*sin(l);
//		z = DL*cos(b);
//
//		xgc = x-gc;
//		ygc = y;
//		zgc = z;
//
//		xec = -0.05487556043*x + 0.4941094278*y - 0.86766614920*z;
//		yec = -0.99382137890*x - 0.1109907351*y - 0.00035159077*z;
//		zec = -0.09647662818*x + 0.8622858751*y + 0.49714719180*z;
//
//		rec = sqrt(xec*xec + yec*yec + zec*zec);
//
//		theta = pi/2.0 - acos(zec/rec);
//		phi   = atan2(yec,xec);
//		while(phi<0.0) phi += 2.0*pi;
//
//		//convert orbital parameters to GW parameters
//		f    =  2.0/Porb;
//		fdot = -2.0*Porb_dot/Porb/Porb;
//
//		//convert mass+distance to amplitude
//		m1 *= Msun;
//		m2 *= Msun;
//		Amp    = 2.0 * G*G*m1*m2 * pow(pi*pi*f*f/(G*(m1+m2)),1./3.) / (DL*pc) /(clight*clight*clight*clight);
//
//		//randomly assign orientation parameters
//		iota = acos(-1.0 + 2.0*ran2(&seed));
//		psi  = ran2(&seed)*2.0*pi;
//		phi0 = ran2(&seed)*2.0*pi;
//
//		fprintf(Outfile, "%.16g %.10g %.6g %.6g %.6g %.6g %.6g %.6g\n",f,fdot,theta,phi,Amp,iota,psi,phi0);
//	}
//	fclose(Infile);
//	fclose(Outfile);
//
//	/************************************************/
//	/*            Interacting Systems               */
//	/************************************************/
//
//	printf("reading AMCVn_GWR_MLDC_bulgefix_opt.dat.1.1...\n");
//	Infile = fopen("AMCVn_GWR_MLDC_bulgefix_opt.dat.1.1","r");
//	Outfile = fopen("AMCVnSourceCatalogue.txt.1.1","w");
//	NS = 8555876;
//	for(n=0; n<NS; n++)
//	{
//		/*Nelemens files*/
//		fscanf(Infile, "%lg %lg %lg %lg %lg %lg %lg %lg", &m1, &m2, &Porb, &Porb_dot, &m2_dot, &l, &b, &DL);
//
//		if(n%1000000==0)printf("Converting AM Cvn: %i/%i (m1=%g, m2=%g, DL=%g)\n",n,NS,m1,m2,DL);
//
//		//convert galactic coordinates to radians
//		l = l*deg2rad;
//		b = (90.0 - b)*deg2rad;
//
//		//convert distance to pc
//		DL *= kpc2pc;
//
//		//change to ecliptic coordinates
//		x = DL*sin(b)*cos(l);
//		y = DL*sin(b)*sin(l);
//		z = DL*cos(b);
//
//		xgc = x-gc;
//		ygc = y;
//		zgc = z;
//
//		xec = -0.05487556043*x + 0.4941094278*y - 0.86766614920*z;
//		yec = -0.99382137890*x - 0.1109907351*y - 0.00035159077*z;
//		zec = -0.09647662818*x + 0.8622858751*y + 0.49714719180*z;
//
//		rec = sqrt(xec*xec + yec*yec + zec*zec);
//
//		theta = pi/2.0 - acos(zec/rec);
//		phi   = atan2(yec,xec);
//		while(phi<0.0) phi += 2.0*pi;
//
//		//convert orbital parameters to GW parameters
//		f    =  2.0/Porb;
//		fdot = -2.0*Porb_dot/Porb/Porb;
//
//		//convert mass+distance to amplitude
//		m1 *= Msun;
//		m2 *= Msun;
//		Amp    = 2.0 * G*G*m1*m2 * pow(pi*pi*f*f/(G*(m1+m2)),1./3.) / (DL*pc) /(clight*clight*clight*clight);
//
//		//randomly assign orientation parameters
//		iota = acos(-1.0 + 2.0*ran2(&seed));
//		psi  = ran2(&seed)*2.0*pi;
//		phi0 = ran2(&seed)*2.0*pi;
//
//		fprintf(Outfile, "%.16g %.10g %.6g %.6g %.6g %.6g %.6g %.6g\n",f,fdot,theta,phi,Amp,iota,psi,phi0);
//	}
//	fclose(Infile);
//	fclose(Outfile);
//
//	printf("reading AMCVn_GWR_MLDC_bulgefix_opt.dat.1.2...\n");
//	Infile = fopen("AMCVn_GWR_MLDC_bulgefix_opt.dat.1.2","r");
//	Outfile = fopen("AMCVnSourceCatalogue.txt.1.2","w");
//	NS = 8555875;
//	for(n=0; n<NS; n++)
//	{
//		/*Nelemens files*/
//		fscanf(Infile, "%lg %lg %lg %lg %lg %lg %lg %lg", &m1, &m2, &Porb, &Porb_dot, &m2_dot, &l, &b, &DL);
//
//		if(n%1000000==0)printf("Converting AM Cvn: %i/%i (m1=%g, m2=%g, DL=%g)\n",n,NS,m1,m2,DL);
//
//		//convert galactic coordinates to radians
//		l = l*deg2rad;
//		b = (90.0 - b)*deg2rad;
//
//		//convert distance to pc
//		DL *= kpc2pc;
//
//		//change to ecliptic coordinates
//		x = DL*sin(b)*cos(l);
//		y = DL*sin(b)*sin(l);
//		z = DL*cos(b);
//
//		xgc = x-gc;
//		ygc = y;
//		zgc = z;
//
//		xec = -0.05487556043*x + 0.4941094278*y - 0.86766614920*z;
//		yec = -0.99382137890*x - 0.1109907351*y - 0.00035159077*z;
//		zec = -0.09647662818*x + 0.8622858751*y + 0.49714719180*z;
//
//		rec = sqrt(xec*xec + yec*yec + zec*zec);
//
//		theta = pi/2.0 - acos(zec/rec);
//		phi   = atan2(yec,xec);
//		while(phi<0.0) phi += 2.0*pi;
//
//		//convert orbital parameters to GW parameters
//		f    =  2.0/Porb;
//		fdot = -2.0*Porb_dot/Porb/Porb;
//
//		//convert mass+distance to amplitude
//		m1 *= Msun;
//		m2 *= Msun;
//		Amp    = 2.0 * G*G*m1*m2 * pow(pi*pi*f*f/(G*(m1+m2)),1./3.) / (DL*pc) /(clight*clight*clight*clight);
//
//		//randomly assign orientation parameters
//		iota = acos(-1.0 + 2.0*ran2(&seed));
//		psi  = ran2(&seed)*2.0*pi;
//		phi0 = ran2(&seed)*2.0*pi;
//
//		fprintf(Outfile, "%.16g %.10g %.6g %.6g %.6g %.6g %.6g %.6g\n",f,fdot,theta,phi,Amp,iota,psi,phi0);
//	}
//	fclose(Infile);
//	fclose(Outfile);
//
//	printf("reading AMCVn_GWR_MLDC_bulgefix_opt.dat.2.1...\n");
//	Infile = fopen("AMCVn_GWR_MLDC_bulgefix_opt.dat.2.1","r");
//	Outfile = fopen("AMCVnSourceCatalogue.txt.2.1","w");
//	NS = 8555873;
//	for(n=0; n<NS; n++)
//	{
//		/*Nelemens files*/
//		fscanf(Infile, "%lg %lg %lg %lg %lg %lg %lg %lg", &m1, &m2, &Porb, &Porb_dot, &m2_dot, &l, &b, &DL);
//
//		if(n%1000000==0)printf("Converting AM Cvn: %i/%i (m1=%g, m2=%g, DL=%g)\n",n,NS,m1,m2,DL);
//
//		//convert galactic coordinates to radians
//		l = l*deg2rad;
//		b = (90.0 - b)*deg2rad;
//
//		//convert distance to pc
//		DL *= kpc2pc;
//
//		//change to ecliptic coordinates
//		x = DL*sin(b)*cos(l);
//		y = DL*sin(b)*sin(l);
//		z = DL*cos(b);
//
//		xgc = x-gc;
//		ygc = y;
//		zgc = z;
//
//		xec = -0.05487556043*x + 0.4941094278*y - 0.86766614920*z;
//		yec = -0.99382137890*x - 0.1109907351*y - 0.00035159077*z;
//		zec = -0.09647662818*x + 0.8622858751*y + 0.49714719180*z;
//
//		rec = sqrt(xec*xec + yec*yec + zec*zec);
//
//		theta = pi/2.0 - acos(zec/rec);
//		phi   = atan2(yec,xec);
//		while(phi<0.0) phi += 2.0*pi;
//
//		//convert orbital parameters to GW parameters
//		f    =  2.0/Porb;
//		fdot = -2.0*Porb_dot/Porb/Porb;
//
//		//convert mass+distance to amplitude
//		m1 *= Msun;
//		m2 *= Msun;
//		Amp    = 2.0 * G*G*m1*m2 * pow(pi*pi*f*f/(G*(m1+m2)),1./3.) / (DL*pc) /(clight*clight*clight*clight);
//
//		//randomly assign orientation parameters
//		iota = acos(-1.0 + 2.0*ran2(&seed));
//		psi  = ran2(&seed)*2.0*pi;
//		phi0 = ran2(&seed)*2.0*pi;
//
//		fprintf(Outfile, "%.16g %.10g %.6g %.6g %.6g %.6g %.6g %.6g\n",f,fdot,theta,phi,Amp,iota,psi,phi0);
//	}
//	fclose(Infile);
//	fclose(Outfile);
//
//	printf("reading AMCVn_GWR_MLDC_bulgefix_opt.dat.2.2...\n");
//	Infile = fopen("AMCVn_GWR_MLDC_bulgefix_opt.dat.2.2","r");
//	Outfile = fopen("AMCVnSourceCatalogue.txt.2.2","w");
//	NS = 8555875;
//	for(n=0; n<NS; n++)
//	{
//		/*Nelemens files*/
//		fscanf(Infile, "%lg %lg %lg %lg %lg %lg %lg %lg", &m1, &m2, &Porb, &Porb_dot, &m2_dot, &l, &b, &DL);
//
//		if(n%1000000==0)printf("Converting AM Cvn: %i/%i (m1=%g, m2=%g, DL=%g)\n",n,NS,m1,m2,DL);
//
//		//convert galactic coordinates to radians
//		l = l*deg2rad;
//		b = (90.0 - b)*deg2rad;
//
//		//convert distance to pc
//		DL *= kpc2pc;
//
//		//change to ecliptic coordinates
//		x = DL*sin(b)*cos(l);
//		y = DL*sin(b)*sin(l);
//		z = DL*cos(b);
//
//		xgc = x-gc;
//		ygc = y;
//		zgc = z;
//
//		xec = -0.05487556043*x + 0.4941094278*y - 0.86766614920*z;
//		yec = -0.99382137890*x - 0.1109907351*y - 0.00035159077*z;
//		zec = -0.09647662818*x + 0.8622858751*y + 0.49714719180*z;
//
//		rec = sqrt(xec*xec + yec*yec + zec*zec);
//
//		theta = pi/2.0 - acos(zec/rec);
//		phi   = atan2(yec,xec);
//		while(phi<0.0) phi += 2.0*pi;
//
//		//convert orbital parameters to GW parameters
//		f    =  2.0/Porb;
//		fdot = -2.0*Porb_dot/Porb/Porb;
//
//		//convert mass+distance to amplitude
//		m1 *= Msun;
//		m2 *= Msun;
//		Amp    = 2.0 * G*G*m1*m2 * pow(pi*pi*f*f/(G*(m1+m2)),1./3.) / (DL*pc) /(clight*clight*clight*clight);
//
//		//randomly assign orientation parameters
//		iota = acos(-1.0 + 2.0*ran2(&seed));
//		psi  = ran2(&seed)*2.0*pi;
//		phi0 = ran2(&seed)*2.0*pi;
//
//		fprintf(Outfile, "%.16g %.10g %.6g %.6g %.6g %.6g %.6g %.6g\n",f,fdot,theta,phi,Amp,iota,psi,phi0);
//	}
//	fclose(Infile);
//	fclose(Outfile);
//
//	printf("combining SourceCatalogue.txt...\n");
//	system("rm SourceCatalogue.txt");
//	system("touch SourceCatalogue.txt");
//	system("cat dwd_GWR_NewLISA.dat.1.1 dwd_GWR_NewLISA.dat.1.2 dwd_GWR_NewLISA.dat.2.1 dwd_GWR_NewLISA.dat.2.2 AMCVnSourceCatalogue.txt.1.1 AMCVnSourceCatalogue.txt.1.2 AMCVnSourceCatalogue.txt.2.1 AMCVnSourceCatalogue.txt.2.2 >> SourceCatalogue.txt");

//	/************************************************/
//	/*            Verification Systems              */
//	/************************************************/
//
//	printf("reading Verification_binaries.dat...\n");
//	Outfile = fopen("VerificationBinaries.txt","w");
//	Infile = fopen("Verification_binaries.dat","r");
//	NS = 35;
//	for(n=0; n<NS; n++)
//	{
//		if(n%10==0)printf("Converting Verification Binaries: %i/%i\n",n,NS);
//		/*Nelemens files*/
//		fscanf(Infile, "%lf%lf%lf%lf%lf%i%lf%lf%lf%i\n", &m1, &m2, &Porb, &Porb_dot, &m2_dot, &i, &l, &b, &DL,&i);
//
//		//convert galactic coordinates to radians
//		l = l*deg2rad;
//		b = (90.0 - b)*deg2rad;
//
//		//convert distance to pc
//		DL *= kpc2pc;
//
//		//change to ecliptic coordinates
//		x = DL*sin(b)*cos(l);
//		y = DL*sin(b)*sin(l);
//		z = DL*cos(b);
//
//		xgc = x-gc;
//		ygc = y;
//		zgc = z;
//
//		xec = -0.05487556043*x + 0.4941094278*y - 0.86766614920*z;
//		yec = -0.99382137890*x - 0.1109907351*y - 0.00035159077*z;
//		zec = -0.09647662818*x + 0.8622858751*y + 0.49714719180*z;
//
//		rec = sqrt(xec*xec + yec*yec + zec*zec);
//
//		theta = pi/2.0 - acos(zec/rec);
//		phi   = atan2(yec,xec);
//		while(phi<0.0) phi += 2.0*pi;
//
//		//convert orbital parameters to GW parameters
//		f    =  2.0/Porb;
//		fdot = -2.0*Porb_dot/Porb/Porb;
//
//		//convert mass+distance to amplitude
//		m1 *= Msun;
//		m2 *= Msun;
//		Amp    = 2.0 * G*G*m1*m2 * pow(pi*pi*f*f/(G*(m1+m2)),1./3.) / (DL*pc) /(clight*clight*clight*clight);
//
//		//randomly assign orientation parameters
//		iota = acos(-1.0 + 2.0*ran2(&seed));
//		psi  = ran2(&seed)*2.0*pi;
//		phi0 = ran2(&seed)*2.0*pi;
//
//		fprintf(Outfile, "%.16g %.10g %.6g %.6g %.6g %.6g %.6g %.6g\n",f,fdot,theta,phi,Amp,iota,psi,phi0);
//	}
//	fclose(Infile);
//	fclose(Outfile);

//  /************************************************/
//  /*            Neutron Star Systems              */
//  /************************************************/
//
//  printf("reading NeutronStarBinariesEM.txt...\n");
//  Outfile = fopen("NeutronStarBinariesGW.txt","w");
//  Infile = fopen("NeutronStarBinariesEM.txt","r");
//  NS = 3;
//  for(n=0; n<NS; n++)
//  {
//    if(n%10==0)printf("Converting Neutron Star Binaries: %i/%i\n",n,NS);
//    /*Nelemens files*/
//    fscanf(Infile, "%lf%lf%lf%lf%lf%i%lf%lf%lf%i\n", &m1, &m2, &Porb, &Porb_dot, &m2_dot, &i, &l, &b, &DL,&i);
//
//    //convert galactic coordinates to radians
//    l = l*deg2rad;
//    b = (90.0 - b)*deg2rad;
//
//    //convert distance to pc
//    DL *= kpc2pc;
//
//    //change to ecliptic coordinates
//    x = DL*sin(b)*cos(l);
//    y = DL*sin(b)*sin(l);
//    z = DL*cos(b);
//
//    xgc = x-gc;
//    ygc = y;
//    zgc = z;
//
//    xec = -0.05487556043*x + 0.4941094278*y - 0.86766614920*z;
//    yec = -0.99382137890*x - 0.1109907351*y - 0.00035159077*z;
//    zec = -0.09647662818*x + 0.8622858751*y + 0.49714719180*z;
//
//    rec = sqrt(xec*xec + yec*yec + zec*zec);
//
//    theta = pi/2.0 - acos(zec/rec);
//    phi   = atan2(yec,xec);
//    while(phi<0.0) phi += 2.0*pi;
//
//    //convert orbital parameters to GW parameters
//    f    =  2.0/Porb;
//    fdot = -2.0*Porb_dot/Porb/Porb;
//
//    //convert mass+distance to amplitude
//    m1 *= Msun;
//    m2 *= Msun;
//    Amp    = 2.0 * G*G*m1*m2 * pow(pi*pi*f*f/(G*(m1+m2)),1./3.) / (DL*pc) /(clight*clight*clight*clight);
//
//    //randomly assign orientation parameters
//    iota = acos(-1.0 + 2.0*ran2(&seed));
//    psi  = ran2(&seed)*2.0*pi;
//    phi0 = ran2(&seed)*2.0*pi;
//
//    fprintf(Outfile, "%.16g %.10g %.6g %.6g %.6g %.6g %.6g %.6g\n",f,fdot,theta,phi,Amp,iota,psi,phi0);
//  }
//  fclose(Infile);
//  fclose(Outfile);
//
//
//    return 0;
//
//}

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define eps 1.2e-7
#define RNMX (1.0-eps)

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

