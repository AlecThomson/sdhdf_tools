//  Copyright (C) 2019, 2020 George Hobbs

/*
 *    This file is part of sdhdfProc. 
 * 
 *    sdhdfProc is free software: you can redistribute it and/or modify 
 *    it under the terms of the GNU General Public License as published by 
 *    the Free Software Foundation, either version 3 of the License, or 
 *    (at your option) any later version. 
 *    sdhdfProc is distributed in the hope that it will be useful, 
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of 
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *    GNU General Public License for more details. 
 *    You should have received a copy of the GNU General Public License 
 *    along with sdhdfProc.  If not, see <http://www.gnu.org/licenses/>. 
 */

//
// Software that calculates statistical properties of spectra
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sdhdfProc.h"
#include "hdf5.h"
#include <cpgplot.h>
#include "TKfit.h"

#define VNUM "v0.1"
#define MAX_CHANS 262144
#define MAX_DUMPS 1024
#define MAX_TRANSMITTERS 1024

void loadBaselineRegions(char *fname, float *baselineX1,float *baselineX2,int *nBaseline);

int main(int argc,char *argv[])
{
  int        i,j,k,l,ii,kk,ll;
  int nAv;
  char       fname[MAX_FILES][64];
  sdhdf_fileStruct *inFile;
  int        idump,iband,ibeam;
  int av=0;
  int sump=0;
  int nx=1;
  int ny=1;
  int polPlot=1;
  double *highResX,*highResY1,*highResY2;
  float *zoomX,*zoomY1,*zoomY2;
  float *zoomY12;
  double *bX,*bY1,*bY2;
  int *bI;
  int nB;
  long writepos;
  int allocateMemory=0;
  float miny,maxy;
  float f0,f1;
  int inFiles;
  char text[1024];
  int sumPol_nPlot;
  char grDev[128];
  int nZoom;
  float xp1,xp2,yp1,yp2;
  float val1,val2;
  float baselineX1[1024],baselineX2[1024];
  int nBaseline;
  double params1[16],params2[16],v[16];
  int nfit;
  float fx[2],fy[2];
  float fc;
  double sx,sx2,sx_2,sx2_2,sdev1,sdev2,mean1,mean2;
  double freq;
  int np;
  int noBaseline=0,nc;
  float freq0 = 1660;
  float freq1 = 1665;
  
  //  loadBaselineRegions("cleanRegions.dat",baselineX1,baselineX2,&nBaseline);
  
  highResX = (double *)malloc(sizeof(double)*MAX_CHANS);
  highResY1 = (double *)calloc(sizeof(double),MAX_CHANS);
  highResY2 = (double *)calloc(sizeof(double),MAX_CHANS);

  // Defaults
  idump = iband = ibeam = 0;
  strcpy(grDev,"/xs");
  
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }

  //iband = 5;
  //  iband=0;
  iband = 0;
  
  sdhdf_initialiseFile(inFile);

  inFiles=0;
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-sb")==0)
	sscanf(argv[++i],"%d",&iband);
      else if (strcmp(argv[i],"-freqRange")==0)
	{
	  sscanf(argv[++i],"%f",&freq0);
	  sscanf(argv[++i],"%f",&freq1);
	}
      else
	strcpy(fname[inFiles++],argv[i]);
    }
  
  for (i=0;i<inFiles;i++)
    {
      sdhdf_initialiseFile(inFile);
      sdhdf_openFile(fname[i],inFile,1);

      sdhdf_loadMetaData(inFile);

      printf("Loading data for band %d\n",iband);
      sdhdf_loadBandData(inFile,ibeam,iband,1);
      
      // Process each band
      //      for (j=0;j<inFile->beam[ibeam].nBand;j++)
      j=iband;
      {
	printf("Processing band: %d\n",j);
	
	// Process each dump
	for (nAv = 1; nAv < inFile->beam[ibeam].bandHeader[j].ndump/2; nAv ++)
	  {
	    //	nAv = 1;
	    for (l=0;l<inFile->beam[ibeam].bandHeader[j].ndump;l+=nAv) 
	      {
		mean1=mean2=0.0;
		for (ll=0;ll<nAv;ll++)
		  {
		    np=0;
		    for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k++)
		      {
			freq = inFile->beam[ibeam].bandData[j].astro_data.freq[k];
			//		    if (freq > 1416 && freq < 1417)
			//		    if (freq > 1416.7 && freq < 1416.8)
			if (freq > freq0 && freq < freq1)
			  {
			    if (ll+l < inFile->beam[ibeam].bandHeader[j].ndump)
			      {
				val1 = inFile->beam[ibeam].bandData[j].astro_data.pol1[(ll+l)*inFile->beam[ibeam].bandHeader[j].nchan+k];
				val2 = inFile->beam[ibeam].bandData[j].astro_data.pol2[(ll+l)*inFile->beam[ibeam].bandHeader[j].nchan+k];
				if (ll==0)
				  {
				    highResX[np] = inFile->beam[ibeam].bandData[j].astro_data.freq[k];
				    highResY1[np] = val1;
				    highResY2[np] = val2;
				  }
				else
				  {
				    highResY1[np] += val1;
				    highResY2[np] += val2;
				  }
				mean1+=val1;
				mean2+=val2;
				np++;
			      }
			  }
		      }
		  }
		if (np > 0)
		  {
		    for (kk=0;kk<np;kk++)
		      {
			highResY1[kk]/=(double)nAv;
			highResY2[kk]/=(double)nAv;
		      }
		    mean1/=(double)(np*nAv);
		    mean2/=(double)(np*nAv);
		    nfit = 1;
		    
		TKremovePoly_d(highResX,highResY1,np,nfit);
		TKremovePoly_d(highResX,highResY2,np,nfit);
		
		// Get mean and variance in the baseline region		
		
		sx=sx2=0;
		sx_2=sx2_2=0;
		nc = 0;
		for (kk=0;kk<np;kk++)
		  {
		    sx    += highResY1[kk];
		    sx2   += pow(highResY1[kk],2);
		    sx_2  += highResY2[kk];
		    sx2_2 += pow(highResY2[kk],2);
		    nc++;
		  }	       		
		sdev1 = sqrt(1./(double)nc*sx2 - pow(1.0/(double)nc * sx,2));
		sdev2 = sqrt(1./(double)nc*sx2_2 - pow(1.0/(double)nc * sx_2,2));
		
		// Mean = before polynomial subtraction
		// Sdev = after polynomial subtraction
		printf("[stats] %d %d %d %g %g %g %g %.5f\n",i,l,ll,mean1,mean2,sdev1,sdev2,inFile->beam[0].bandData[j].astro_obsHeader[ll+l].mjd);	  
		  }
	      }
	  }
      }
      
           
      sdhdf_closeFile(inFile);

    }

  
  free(highResX); free(highResY1); free(highResY2);
}


void loadBaselineRegions(char *fname, float *baselineX1,float *baselineX2,int *nBaseline)
{
  FILE *fin;

  if (!(fin = fopen(fname,"r")))
    {
      printf("Unable to open >%s<\n",fname);
      exit(1);
    }
  *nBaseline = 0;
  while (!feof(fin))
    {
      if (fscanf(fin,"%f %f",&baselineX1[*nBaseline],&baselineX2[*nBaseline])==2)
	(*nBaseline)++;
    }
  fclose(fin);
}
