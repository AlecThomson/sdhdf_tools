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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sdhdfProc.h"
#include "hdf5.h"
#include <cpgplot.h>
#include "TKfit.h"


#define MAX_DUMPS 409600


typedef struct dumpStruct {
  char fname[128];
  int  dump;
  float az;
  float el;
  float scale;
} dumpStruct;


void makePlot(float *signalVal,int nchan,int ndump,float f0,float chbw,dumpStruct *dumpParams);

int main(int argc,char *argv[])
{
  int i,j,k,ii,jj,kk,b;
  sdhdf_fileStruct *inFile;
  int nchan,idump,nband;
  int totSize,chanPos;
  int ndump;
  //  spectralDumpStruct spectrum;
  int npol=4;
  int ibeam=0;
  int iband=0;
  int ichan;
  int setnband=-1;
  float tdump;
  int   np;
  int sdump=0;
  int sdump0=0;
  char aest0[MAX_STRLEN]="00:00:00";

  dumpStruct *dumpParams;
  double freq;
  float *timeVal;
  int min=1;
  float *elVal,*azVal;
  float *signalVal;

  FILE *fin;
  FILE *fout;
  FILE *fout2;
  char outName[1024];
  char header[1024];  
  char bandFile[1024] = "NULL";
  int  nFiles=0;
  char fname[MAX_FILES][MAX_STRLEN];

  int nsb=0;
  float *sbF0,*sbF1;
  char **bandHeader;
  int log=0;
  char labelFile[MAX_STRLEN] = "NULL";
  char xlabel[MAX_STRLEN];
  int nplot;
  float miny = -1,maxy = -1;
  float gap = 1800; // Gap length in seconds
  int i0=0;
  int sm=1;
  int nc;
  int div1=0;
  int sub0=0;

  float f0,f1,chbw;
  
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-log")==0)
	log=1;
      else
	strcpy(fname[nFiles++],argv[i]);
    }

  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }


  sdump=sdump0=0;
  for (i=0;i<nFiles;i++)
    {
      sdhdf_initialiseFile(inFile);
      sdhdf_openFile(fname[i],inFile,1);
      sdhdf_loadMetaData(inFile);
      
      // Find out which UWL bands we will need to load and process
      printf("Processing SB %d (%s)\n",j,inFile->fname);
      sdhdf_loadBandData(inFile,ibeam,iband,1);
      
      ndump = inFile->beam[ibeam].bandHeader[iband].ndump;
      nchan = inFile->beam[ibeam].bandHeader[iband].nchan;
      if (i==0)
	{
	  signalVal = (float *)malloc(sizeof(float)*MAX_DUMPS*nchan);
	  dumpParams = (dumpStruct *)malloc(sizeof(dumpStruct)*MAX_DUMPS);
	  f0 = inFile->beam[ibeam].bandData[iband].astro_data.freq[0];
	  f1 = inFile->beam[ibeam].bandData[iband].astro_data.freq[1];
	  chbw = f1-f0;
	}

      for (jj=0;jj<ndump;jj++)
	{
	  strcpy(dumpParams[jj+sdump].fname,inFile->fname);
	  dumpParams[jj+sdump].dump = jj;
	  dumpParams[jj+sdump].scale = 1;
	  //	  printf("Setting fname %s %d >%s<\n",inFile->fname,jj+sdump,dumpParams[jj+sdump].fname);
	  for (ii=0;ii<nchan;ii++)
	    {
	      signalVal[(jj+sdump)*nchan+ii] = (inFile->beam[ibeam].bandData[iband].astro_data.pol1[ii+jj*nchan]+inFile->beam[ibeam].bandData[iband].astro_data.pol2[ii+jj*nchan]);
	      signalVal[(jj+sdump)*nchan+ii] /= (double)inFile->beam[ibeam].bandHeader[iband].dtime;
	      //	      printf("dtim = %g\n",(double)inFile->beam[ibeam].bandHeader[iband].dtime);
	      signalVal[(jj+sdump)*nchan+ii] = log10(signalVal[(jj+sdump)*nchan+ii]);
	    }
	}
      sdhdf_closeFile(inFile);
      sdump+=ndump;
    }

  makePlot(signalVal,nchan,sdump,f0,chbw,dumpParams);
  

  free(signalVal);
  free(dumpParams);
}

void makePlot(float *signalVal,int nchan,int ndump,float f0,float chbw,dumpStruct *dumpParams)
{
  int i,j;
  float tr[6];
  float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
  float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
  float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
  float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};

  float specMiny,specMaxy;
  float *specX,*specY,*specY2,*specY3,*specY4;
  float val;

  float mx,my;
  char key;
  float *heatMap;
  
  int specSelect=-1;
  float minz,maxz;
  
  specX = (float *)malloc(sizeof(float)*nchan);
  specY = (float *)malloc(sizeof(float)*nchan);
  specY2 = (float *)malloc(sizeof(float)*nchan);
  specY3 = (float *)malloc(sizeof(float)*nchan);
  specY4 = (float *)malloc(sizeof(float)*nchan);
  heatMap = (float *)malloc(sizeof(float)*nchan*ndump);
  
  printf("Making plot with nchan = %d, ndump = %d\n",nchan,ndump);

  tr[0] = f0; tr[1] = chbw; tr[2] = 0;
  tr[3] = 0; tr[4] = 0; tr[5] = 1;
  cpgbeg(0,"/xs",1,1);
  cpgsch(1.4);
  cpgslw(2);
  cpgscf(2);

  do
    {
      
      for (i=0;i<nchan;i++)
	{
	  specX[i] = f0+chbw*i;
	  specY[i] = 0.0;
	  
	  
	  for (j=0;j<ndump;j++)
	    {
	      val = dumpParams[j].scale*pow(10,signalVal[j*nchan+i]);
	      heatMap[j*nchan+i] = log10(val);
	      if (i==0 && j==0)
		{
		  minz = maxz = heatMap[j*nchan+i];
		}
	      if (minz > heatMap[j*nchan+i]) minz = heatMap[j*nchan+i];
	      if (maxz < heatMap[j*nchan+i]) maxz = heatMap[j*nchan+i];
	      if (j==specSelect)
		specY4[i] = log10(val);
	      if (j==0)
		{
		  specY2[i] = specY3[i] = log10(val);
		}
	      else
		{
		  if (specY2[i] > log10(val)) specY2[i] = log10(val);
		  if (specY3[i] < log10(val)) specY3[i] = log10(val);
		}
	      specY[i] += val;
	    }
	  specY[i] = log10(specY[i]/ndump);
	  if (i==0)
	    {
	      specMiny = specY2[i];
	      specMaxy = specY3[i];
	    }
	  if (specMiny > specY2[i]) specMiny = specY2[i];
	  if (specMaxy < specY3[i]) specMaxy = specY3[i];      
	}
      cpgeras();
      cpgsvp(0.1,0.95,0.15,0.4);
      cpgswin(f0,f0+nchan*chbw,specMiny,specMaxy);
      cpgbox("ABCTSN",0,0,"ABCTSN",0,0);
      cpglab("Frequency (MHz)","","");
      cpgline(nchan,specX,specY);
      cpgsci(2); cpgline(nchan,specX,specY2);
      cpgsci(4); cpgline(nchan,specX,specY3);
      if (specSelect >= 0) {cpgsci(3); cpgline(nchan,specX,specY4);}
      cpgsci(1);
      
      cpgsvp(0.1,0.95,0.4,0.95);
      cpgswin(f0,f0+nchan*chbw,0,ndump);
      
      cpglab("","Spectral dump","");
      cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
      cpgimag(heatMap,nchan,ndump,1,nchan,1,ndump,minz,maxz,tr);
      cpgbox("ABCTS",0,0,"ABCTSN",0,0);

      cpgcurs(&mx,&my,&key);
      if (key=='A')
	{
	  printf("Mouse clicked at (%g,%g)\n",mx,my);
	  specSelect = (int)(my+0.5);
	  printf("Selected spectrum #%d\n",specSelect);
	  printf("Selecting spectra from %s, spectral dump number %d\n",dumpParams[specSelect].fname,dumpParams[specSelect].dump);
	}
      else if (key=='s') // Scale
	{
	  for (i=0;i<ndump;i++)
	    dumpParams[i].scale = 1.0/pow(10,signalVal[i*nchan+115]); // 800]); // FIX ***
	}
    } while (key != 'q');
      
  cpgend();

  free(specX);
  free(specY);
  free(specY2);
  free(specY3);
  free(specY4);
  free(heatMap);
}

