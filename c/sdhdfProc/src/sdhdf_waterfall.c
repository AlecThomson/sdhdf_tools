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


void makePlot(float *signalVal,int nchan,int ndump,float f0,float chbw,dumpStruct *dumpParams,float specMinx,float specMaxx,float cleanF0,float cleanF1,float *wt);

int main(int argc,char *argv[])
{
  int i,j,k,ii,jj,kk,b;
  sdhdf_fileStruct *inFile;
  int nchan,idump,nband;
  int totSize,chanPos;
  int ndump;
  int setNdump=MAX_DUMPS;
  float minx=-1,maxx=-1;
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
  float *wt;
  
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
  int poln=0;
  
  float f0,f1,chbw;
  float cleanF0=-1,cleanF1=-1;
  
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-log")==0)
	log=1;
      else if (strcmp(argv[i],"-pol")==0)
	sscanf(argv[++i],"%d",&poln);
      else if (strcmp(argv[i],"-ndump")==0)
	sscanf(argv[++i],"%d",&setNdump);
      else if (strcmp(argv[i],"-band")==0)
	sscanf(argv[++i],"%d",&iband);
      else if (strcasecmp(argv[i],"-frange")==0 || strcasecmp(argv[i],"-freqrange")==0)
	{
	  sscanf(argv[++i],"%f",&minx);
	  sscanf(argv[++i],"%f",&maxx);
	}
      else if (strcasecmp(argv[i],"-crange")==0 || strcasecmp(argv[i],"-cleanrange")==0)
	{
	  sscanf(argv[++i],"%f",&cleanF0);
	  sscanf(argv[++i],"%f",&cleanF1);
	}
      else
	strcpy(fname[nFiles++],argv[i]);
    }

  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      printf("Try and reduce the total number of spectral dumps using -ndump\n");
      printf("or reduce the frequency resolution of your data\n");
      exit(1);
    }


  sdump=sdump0=0;
  for (i=0;i<nFiles;i++)
    {
      sdhdf_initialiseFile(inFile);
      sdhdf_openFile(fname[i],inFile,1);
      sdhdf_loadMetaData(inFile);
      
      // Find out which UWL bands we will need to load and process
      printf("Processing SB %d (%s)\n",iband,inFile->fname);
      sdhdf_loadBandData(inFile,ibeam,iband,1);
      
      ndump = inFile->beam[ibeam].bandHeader[iband].ndump;
      nchan = inFile->beam[ibeam].bandHeader[iband].nchan;
      if (i==0)
	{
	  if (!(signalVal = (float *)malloc(sizeof(float)*setNdump*nchan)))
	    {
	      printf("ERROR: unable to allocate enough memory for the data: %g\n",(float)(setNdump*nchan));
	      exit(1);	      
	    }
	  wt = (float *)malloc(sizeof(float)*nchan*setNdump);
	  dumpParams = (dumpStruct *)malloc(sizeof(dumpStruct)*setNdump);
	  f0 = inFile->beam[ibeam].bandData[iband].astro_data.freq[0];
	  f1 = inFile->beam[ibeam].bandData[iband].astro_data.freq[1];
	  chbw = f1-f0;
	}
    
      for (jj=0;jj<ndump;jj++)
	{
	  printf("Processing dump %d\n",jj);
	  strcpy(dumpParams[jj+sdump].fname,inFile->fname);
	  dumpParams[jj+sdump].dump = jj;
	  dumpParams[jj+sdump].scale = 1;
	  for (ii=0;ii<nchan;ii++)
	    {
	      if (poln == 1)		
		signalVal[(jj+sdump)*nchan+ii] = (inFile->beam[ibeam].bandData[iband].astro_data.pol1[ii+jj*nchan]);
	      else
		signalVal[(jj+sdump)*nchan+ii] = (inFile->beam[ibeam].bandData[iband].astro_data.pol1[ii+jj*nchan]+inFile->beam[ibeam].bandData[iband].astro_data.pol2[ii+jj*nchan]);

	      wt[(jj+sdump)*nchan+ii] = inFile->beam[ibeam].bandData[iband].astro_data.dataWeights[ii+jj*nchan];
	      //	      if (jj==1)
	      //		printf("flag1: %d %d %g\n",ii,jj,wt[(jj*sdump)*nchan+ii]);
	      
	      signalVal[(jj+sdump)*nchan+ii] /= (double)inFile->beam[ibeam].bandHeader[iband].dtime;
	      //	      printf("dtim = %g\n",(double)inFile->beam[ibeam].bandHeader[iband].dtime);
	      signalVal[(jj+sdump)*nchan+ii] = log10(signalVal[(jj+sdump)*nchan+ii]);
	    }
	}
      sdhdf_closeFile(inFile);
      sdump+=ndump;
    }

  /*
  for (ii=0;ii<nchan;ii++)
    {
      for (jj=0;jj<sdump;jj++)
	{
	  if (jj==1)
	    printf("flag2: %d %d %g\n",ii,jj,wt[(jj*sdump)*nchan+ii]);

	}
    }
  */  
  makePlot(signalVal,nchan,sdump,f0,chbw,dumpParams,minx,maxx,cleanF0,cleanF1,wt);
  
  free(wt);
  free(signalVal);
  free(dumpParams);
}

void makePlot(float *signalVal,int nchan,int ndump,float f0,float chbw,dumpStruct *dumpParams,float specMinx,float specMaxx,float cleanF0,float cleanF1,float *wt)
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
  float minz=0,maxz=0;

  float heatMiny = 0;
  float heatMaxy = ndump;
  

  int t=0,t_z=0;
  int log=-1;
  int colourScale=1;
  int imax,jmax;
  
  if (specMinx == -1) specMinx = f0;
  if (specMaxx == -1) specMaxx = f0+nchan*chbw;
  
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
	}
      if (cleanF0 > 0 &&  cleanF1 > 0)
	{	  
	  float scl;
	  int nv=0;
	  printf("Cleaning the data\n");
	  
	  for (i=0;i<ndump;i++)
	    {
	      nv=0;
	      scl=0;
	      for (j=0;j<nchan;j++)
		{
		  if (specX[j] > cleanF0 && specX[j] <= cleanF1 && wt[i*nchan+j] != 0)
		    {
		      scl += pow(10,signalVal[i*nchan+j]);
		      nv++;
		    }	
	  if (nv > 0)
		    dumpParams[i].scale = 1.0/(scl/nv);
		  else
		    dumpParams[i].scale = 1.0;
		}
	    }
	}

      t=0;
      t_z=0;
      for (i=0;i<nchan;i++)
	{
	  specY[i] = 0.0;
	  
	  
	  for (j=0;j<ndump;j++)
	    {
	      if (wt[j*nchan+i] == 0)
		val = 0;
	      else		
		val = dumpParams[j].scale*pow(10,signalVal[j*nchan+i]);
	      //	      if (j==1)
		//		printf("Have %d %d %g\n",i,j,wt[j*nchan+i]);
	      //	      if (val > 0 && j > 0)
	      //		printf("%d %d Val = %g\n",i,j,val);
	      if (log==1 && val > 0)
		heatMap[j*nchan+i] = log10(val);
	      else
		heatMap[j*nchan+i] = val;
	      
	      if (specX[i] > specMinx && specX[i] < specMaxx  && wt[j*nchan+i] > 0)
		{
		  if (t_z==0)
		    {
		      imax = i; jmax = j;
		      minz = maxz = heatMap[j*nchan+i];
		      printf("In here setting %g %g\n",minz,maxz);
		      t_z=1;
		    }
		  else
		    {
		      if (minz > heatMap[j*nchan+i]) minz = heatMap[j*nchan+i];
		      if (maxz < heatMap[j*nchan+i]) {maxz = heatMap[j*nchan+i]; imax = i; jmax = j;}
		    }
		    }
	    
	      if (j==specSelect)
		{
		  if (log==1)
		    specY4[i] = log10(val);
		  else
		    specY4[i] = (val);
		}
	      if (j==0)
		{
		  if (log==1)
		    specY2[i] = specY3[i] = log10(val);
		  else
		    specY2[i] = specY3[i] = val;
		}
	      else
		{
		  if (log==1)
		    {
		      if (specY2[i] > log10(val)) specY2[i] = log10(val);
		      if (specY3[i] < log10(val)) specY3[i] = log10(val);
		    }
		  else
		    {
		      if (specY2[i] > log10(val)) specY2[i] = (val);
		      if (specY3[i] < log10(val)) specY3[i] = (val);
		    }
		}
	      specY[i] += val;
	    }
	  if (log==1)
	    specY[i] = log10(specY[i]/ndump);
	  else
	    specY[i] = (specY[i]/ndump);
	
	  if (specX[i] > specMinx && specX[i] < specMaxx)
	    {
	      if (t==0)
		{
		  specMiny = specY2[i];
		  specMaxy = specY3[i];
		  t=1;
		}
	      
	      if (specMiny > specY2[i]) specMiny = specY2[i];
	      if (specMaxy < specY3[i]) specMaxy = specY3[i];      
	    }
	}
      //      maxz = 2000;
      printf("Have colour range from %g to %g\n",minz,maxz);
      printf("Brightest pixel at (%d,%d) corresponding to frequency %.6f\n",imax,jmax,specX[imax]);
      printf("X range from %g to %g\n",specMinx,specMaxx);
      printf("Y range from %g to %g\n",specMiny,specMaxy);
      cpgeras();
      cpgsvp(0.1,0.95,0.15,0.4);
      cpgswin(specMinx,specMaxx,specMiny,specMaxy);
      cpgbox("ABCTSN",0,0,"ABCTSN",0,0);
      cpglab("Frequency (MHz)","","");
      cpgline(nchan,specX,specY);
      cpgsci(2); cpgline(nchan,specX,specY2);
      cpgsci(4); cpgline(nchan,specX,specY3);
      if (specSelect >= 0) {cpgsci(3); cpgline(nchan,specX,specY4);}
      cpgsci(1);
      
      cpgsvp(0.1,0.95,0.4,0.95);
      cpgswin(specMinx,specMaxx,heatMiny,heatMaxy+0.5);
      
      cpglab("","Spectral dump","");
      cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
      if (colourScale==1)
	cpgimag(heatMap,nchan,ndump,1,nchan,1,ndump,minz,maxz,tr);
      else if (colourScale==2)
	cpgimag(heatMap,nchan,ndump,1,nchan,1,ndump,maxz,minz,tr);
      else if (colourScale==3)
	cpggray(heatMap,nchan,ndump,1,nchan,1,ndump,maxz,minz,tr);
      else if (colourScale==4)
	cpggray(heatMap,nchan,ndump,1,nchan,1,ndump,minz,maxz,tr);
	cpgbox("ABCTS",0,0,"ABCTSN",0,0);

      cpgcurs(&mx,&my,&key);
      if (key=='A')
	{
	  printf("Mouse clicked at (%g,%g)\n",mx,my);
	  specSelect = (int)(my-0.5);
	  if (specSelect < 0) specSelect = 0;
	  if (specSelect >= ndump) specSelect = ndump-1;
	  
	  printf("Selected spectrum #%d\n",specSelect);
	  printf("Selecting spectra from %s, spectral dump number %d\n",dumpParams[specSelect].fname,dumpParams[specSelect].dump);
	}
      else if (key=='c')
	{
	  colourScale++;
	  if (colourScale==5) colourScale=1;
	}
      else if (key=='l')
	log*=-1;
      else if (key=='y')
	{
	  float mx2,my2;
	  
	  cpgband(3,0,mx,my,&mx2,&my2,&key);
	  if (my < my2) {heatMiny = my; heatMaxy = my2;}
	  else {heatMiny = my2; heatMaxy = my;}
	}
      else if (key=='x')
	{
	  float mx2,my2;
	  
	  cpgband(4,0,mx,my,&mx2,&my2,&key);
	  if (mx < mx2) {specMinx = mx; specMaxx = mx2;}
	  else {specMinx = mx2; specMaxx = mx;}
	}
      else if (key=='f')
	{
	  printf("Enter frequency range: ");
	  scanf("%f %f",&specMinx,&specMaxx);
	}
      else if (key=='u')
	{
	  heatMiny = 0;
	  heatMaxy = ndump;
	  specMinx = f0;
	  specMaxx = f0+nchan*chbw;	  
	}
      else if (key=='a')
	specSelect=-1;
      else if (key=='s') // Scale
	{
	  printf("Enter a clean frequency range: f0 f1 ");
	  scanf("%f %f",&cleanF0,&cleanF1);
	  
	  
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

