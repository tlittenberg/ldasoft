#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <math.h>

int main(int argc,char **argv)
{

	char ifilename[100];
	char WDWDfilename[100];
	char AMCVfilename[100];
	
	sprintf(ifilename,"%s",argv[1]);
	sprintf(WDWDfilename,"DetachedBinaries.txt");
	sprintf(AMCVfilename,"InteractingBinaries.txt");
	
	FILE *ifile = fopen(ifilename,"r");
	FILE *WDWDfile = fopen(WDWDfilename,"w");
	FILE *AMCVfile = fopen(AMCVfilename,"w");
	
	int i;
	double params[8];
	while(!feof(ifile))
	{
		for(i=0;i<8;i++)fscanf(ifile,"%lg",&params[i]);

		if(params[1]>=0) 
		{
			for(i=0;i<8;i++)fprintf(WDWDfile,"%lg ",params[i]);
			fprintf(WDWDfile,"\n");
		}
		if(params[1]<0) 
		{
			for(i=0;i<8;i++)fprintf(AMCVfile,"%lg ",params[i]);
			fprintf(AMCVfile,"\n");
		}
		
	}
	
	fclose(ifile);
	fclose(WDWDfile);
	fclose(AMCVfile);
	
	return 0;
	
	
	
	
	
	
	
//	char ifilename[100];
//	char ofilename[100];
//	
//	sprintf(ifilename,"%s",argv[1]);
//	sprintf(ofilename,"Reform%s",argv[1]);
//	
//	FILE *ifile = fopen(ifilename,"r");
//	FILE *ofile = fopen(ofilename,"w");
//	
//	int i;
//	double junk;
//	double sigmas[9];
//	while(!feof(ifile))
//	{
//		for(i=0;i<9;i++)fscanf(ifile,"%lg",&junk);
//		for(i=0;i<9;i++)fscanf(ifile,"%lg",&sigmas[i]);
//		fscanf(ifile,"%lg %lg",&junk,&junk);
//		for(i=0;i<9;i++)fprintf(ofile,"%lg ",log10(sigmas[i]));
//		fprintf(ofile,"\n");
//	}
//	
//	fclose(ifile);
//	fclose(ofile);
//	
//	return 0;
}
