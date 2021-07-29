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
// Software to plot a spectrum
//
// Usage:
// sdhdf_plotSpectrum -f <filename>
//
// Compilation
// gcc -lm -o sdhdf_plotMultiSpec sdhdf_plotMultiSpec.c sdhdfProc.c -I../hdf5-1.10.4/src/ ../hdf5-1.10.4/src/.libs/libhdf5.a -ldl -lz -lcpgplot -L/pulsar/psr/software/stable/stretch/lib/ -I//pulsar/psr/software/stable/stretch/include -lcalceph -Isofa/20190722/c/src/ -Lsofa/20190722/c/src -lsofa_c

//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sdhdfProc.h"
#include "hdf5.h"
#include <cpgplot.h>

#define VNUM "v0.1"

void drawIncludeWeights(int nchan,float *freq,float *pol,float *wt);
void drawMolecularLine(float freq,char *label,float minX,float maxX,float minY,float maxY);
void plotSpectrum(sdhdf_fileStruct *inFile,int ibeam, int iband,int idump,char *grDev,char *fname,float f0,int setf0,float f1,int setf1,int av,int sump,int nx,int ny,int polPlot,float chSize,float locky1,float locky2,int join,double fref,float bl_f0,float bl_f1,int setBaseline,int setLog);


void help()
{
  printf("sdhdf_plotMultiSpec: %s\n",VNUM);
  printf("sdhfProc version:   %s\n",SOFTWARE_VER);
  printf("author:             George Hobbs\n");
  printf("\n");
  printf("Software to plot multiple spectra\n");
  printf("\n\nCommand line arguments:\n\n");


  printf("-4pol               Plot 4 polarisations (default is just 1 or 2 polarisation)\n");
  printf("-av <av>            Frequency average this number of channels\n");
  printf("-ch <chSize>        Character size\n");
  printf("-h                  This help\n");
  printf("-f0 <freq0>         Start frequency (should be within specified band)\n");
  printf("-f1 <freq1>         End frequency (should be within specified band)\n");
  printf("-f <filename>       SDHDF file corresponding to observations\n");
  printf("-fref <freq>        Set a reference frequency in MHz\n");
  printf("-g <display>        Set display type\n");
  printf("-join               Join plots together\n");
  printf("-locky <y1> <y2>    Fix the y-axis to be between y1 and y2\n");
  printf("-log                Set logarithmic plotting\n");
  printf("-nx <nx>            Set nx panels in the x-direction\n");  
  printf("-ny <ny>            Set ny panels in the y-direction\n");
  printf("-p                  Sum polarisations\n");
  printf("-sb <bandNumber>    Plot data within specified band number\n");
  printf("-sd <dumpNumber>    Plot data within specified spectral dump\n");
  
    
  
  
  
  printf("\nExample:\n\n");
  printf("sdhdf_plotMultiSpec -sb 0 -f0 810 -f1 812 uwl_210722_214441.hdf.T.flag\n");
  printf("---------------------\n");
}


int main(int argc,char *argv[])
{
  int        i,j,ii;
  char       fname[MAX_STRLEN];
  char       grDev[MAX_STRLEN];
  sdhdf_fileStruct *inFile;
  int        idump,iband,ibeam;
  float chSize=1.4;
  float f0=-1,f1=-1;
  int setf0=0,setf1=0;
  int av=0;
  int sump=0;
  int nx=1;
  int ny=1;
  int polPlot=1;
  float locky1=-1,locky2=-1;
  int join=0;
  double fref=-1;
  float bl_f0=0,bl_f1=0;
  int setBaseline=0;
  int setLog=-1;
  
  //  help();
  
  // Defaults
  ibeam = 0;
  idump = iband = 0;
  strcpy(grDev,"/xs");
  
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }

  sdhdf_initialiseFile(inFile);
  
  for (i=1;i<argc;i++)
    {      
      if (strcmp(argv[i],"-h")==0)
	{help(); exit(1);}
      else if (strcmp(argv[i],"-log")==0)
	setLog=1;      
      else if (strcmp(argv[i],"-g")==0)
	strcpy(grDev,argv[++i]);
      else if (strcmp(argv[i],"-fref")==0)
	sscanf(argv[++i],"%lf",&fref);
      else if (strcmp(argv[i],"-nx")==0) 
	sscanf(argv[++i],"%d",&nx);
      else if (strcmp(argv[i],"-ny")==0)
	sscanf(argv[++i],"%d",&ny);
      else if (strcmp(argv[i],"-p")==0)
	sump=1;
      else if (strcmp(argv[i],"-join")==0)
	join=1;
      else if (strcmp(argv[i],"-av")==0)
	sscanf(argv[++i],"%d",&av);
      else if (strcmp(argv[i],"-bl")==0)
	{
	  sscanf(argv[++i],"%f",&bl_f0);
	  sscanf(argv[++i],"%f",&bl_f1);
	  setBaseline=1;
	}
      else if (strcmp(argv[i],"-f0")==0)
	{sscanf(argv[++i],"%f",&f0); setf0=1;}
      else if (strcmp(argv[i],"-f1")==0)
	{sscanf(argv[++i],"%f",&f1); setf1=1;}
      else if (strcmp(argv[i],"-locky")==0)
	{
	  sscanf(argv[++i],"%f",&locky1);
	  sscanf(argv[++i],"%f",&locky2);
	}
      else if (strcmp(argv[i],"-sd")==0)
	sscanf(argv[++i],"%d",&idump);
      else if (strcmp(argv[i],"-4pol")==0)
	polPlot=2;
      else if (strcmp(argv[i],"-ch")==0) // Character height
	sscanf(argv[++i],"%f",&chSize);
      else if (strcmp(argv[i],"-sb")==0)
	sscanf(argv[++i],"%d",&iband);
	//	i++;
      else
	{
	  strcpy(fname,argv[i]);

	  sdhdf_initialiseFile(inFile);
	  if (sdhdf_openFile(fname,inFile,1)==-1)
	    {
	      printf("Unable to open input file >%s<\n",fname);
	      free(inFile);
	      exit(1);
	    }
	  sdhdf_loadMetaData(inFile);
	  printf("%-22.22s %-22.22s %-5.5s %-6.6s %-20.20s %-10.10s %-10.10s %-7.7s %3d\n",fname,inFile->primary[0].utc0,
		 inFile->primary[0].hdr_defn_version,inFile->primary[0].pid,inFile->beamHeader[ibeam].source,inFile->primary[0].telescope,
		 inFile->primary[0].observer, inFile->primary[0].rcvr,inFile->beam[ibeam].nBand);

	  // FIX ME
	  //	  for (ii=0;ii<argc;ii++)
	  //	    {
	  //	      if (strcmp(argv[ii],"-sb")==0)
	  //		iband = sdhdf_getBandID(inFile,argv[++ii]);
	  //	    }


	  plotSpectrum(inFile,ibeam, iband,idump,grDev,fname,f0,setf0,f1,setf1,av,sump,nx,ny,polPlot,chSize,locky1,locky2,join,fref,bl_f0,bl_f1,setBaseline,setLog);	  	  

	  sdhdf_closeFile(inFile);

	}
    }
  cpgend();  
  free(inFile);
}


void plotSpectrum(sdhdf_fileStruct *inFile,int ibeam, int iband,int idump,char *grDev,char *fname,float f0,int setf0,float f1,int setf1,int av,int sump,int nx,int ny,int polPlot,float chSize,float locky1,float locky2,int join,double fref,float bl_f0,float bl_f1,int setBaseline,int setLog)
{
  static int entry=0;
  static int entryX=0;
  static int entryY=0;
  //  spectralDumpStruct spectrum;
  char key;
  float mx,my;
  int i,j,k;
  float *aa,*bb,*ab,*abs;
  float minx,maxx,miny,maxy,minz,maxz;
  float ominx,omaxx,ominy,omaxy;
  int t=0;
  char title[1024];
  char xlabel[1024];
  char ylabel[1024];
  int sumAll=-1;
  int plot=1;
  int molecularLines=1;
  int nchan,maxNchan;
  float *pol1,*pol2,*pol3,*pol4;
  float *wt;
  int *useP1,*useP2;
  float *freq;
  float xp1,xp2,yp1,yp2;
  float wx1 = 0.1, wx2 = 0.9;
  float wy1 = 0.15, wy2 = 0.9;
  float mean1,mean2,mean3,mean4;
  int nc=0;
  int npol;


  mean1=mean2=mean3=mean4=0.0;
  npol = inFile->beam[ibeam].bandHeader[iband].npol;
  
  if (entry==0)
    {
      if (join==0)
	cpgbeg(0,grDev,nx,ny);
      else
	cpgbeg(0,grDev,1,1);
      cpgscf(2);
      cpgsch(chSize);
      cpgslw(2);
    }

  maxNchan=0;
  for (i=0;i<inFile->beam[ibeam].nBand;i++)
    {
      if (inFile->beam[ibeam].bandHeader[iband].nchan > maxNchan)
	maxNchan = inFile->beam[ibeam].bandHeader[iband].nchan;
    }

  nchan = inFile->beam[ibeam].bandHeader[iband].nchan;
  freq = (float *)malloc(sizeof(float)*nchan);
  useP1 = (int *)malloc(sizeof(int)*nchan);
  useP2 = (int *)malloc(sizeof(int)*nchan);
  wt = (float *)malloc(sizeof(float)*nchan);
  pol1 = (float *)malloc(sizeof(float)*nchan);
  pol2 = (float *)malloc(sizeof(float)*nchan);
  pol3 = (float *)malloc(sizeof(float)*nchan);
  pol4 = (float *)malloc(sizeof(float)*nchan);
  
  sdhdf_loadBandData(inFile,ibeam,iband,1);
  
  if (av < 1)
    {
      for (i=0;i<nchan;i++)
	{
	  useP1[i] = 1;
	  useP2[i] = 1;
	  if (fref < 0)
	    freq[i] = inFile->beam[ibeam].bandData[iband].astro_data.freq[i];
	  else
	    freq[i] = (1.0-inFile->beam[ibeam].bandData[iband].astro_data.freq[i]/(fref))*SPEED_LIGHT/1000.; // km/s
	  pol1[i] = inFile->beam[ibeam].bandData[iband].astro_data.pol1[i+idump*nchan];
	  wt[i]   = inFile->beam[ibeam].bandData[iband].astro_data.dataWeights[i+idump*nchan];
	  if (wt[i] == 0) {useP1[i] = 0; useP2[i] = 0;}
	  if (npol > 1)
	    pol2[i] = inFile->beam[ibeam].bandData[iband].astro_data.pol2[i+idump*nchan];
	  if (polPlot==2)
	    {
	      pol3[i] = inFile->beam[ibeam].bandData[iband].astro_data.pol3[i+idump*nchan];
	      pol4[i] = inFile->beam[ibeam].bandData[iband].astro_data.pol4[i+idump*nchan];
	    }
	  if (sump==1 && npol > 1) pol1[i] = pol1[i]+pol2[i];
	  
	  if (setLog == 1)
	    {
	      if (pol1[i] > 0)
		pol1[i]=log10(pol1[i]);
	      else
		{pol1[i] = -100; useP1[i] = 0;}
	      
	      if (npol > 1)
		{
		  if (pol2[i] > 0)
		    pol2[i]=log10(pol2[i]);
		  else
		    {pol2[i] = -100; useP2[i] = 0;}
		}
	    }
	  if (freq[i] > bl_f0 && freq[i] <= bl_f1)
	    {
	      mean1+=pol1[i];
	      mean2+=pol2[i];
	      mean3+=pol3[i];
	      mean4+=pol4[i];
	      nc++;
	    }
	}
    }
  else
    {
      float avF,av1,av2,av3,av4;
      int n=0;
      for (i=0;i<nchan;i+=av)
	{
	  avF=av1=av2=av3=av4=0;
	  for (j=i;j<i+av;j++)
	    {
	      avF += inFile->beam[ibeam].bandData[iband].astro_data.freq[j];
	      av1 += inFile->beam[ibeam].bandData[iband].astro_data.pol1[j+idump*nchan];
	      if (npol > 1)
		av2 += inFile->beam[ibeam].bandData[iband].astro_data.pol2[j+idump*nchan];
	      if (polPlot==2)
		{
		  av3 += inFile->beam[ibeam].bandData[iband].astro_data.pol3[j+idump*nchan];
		  av4 += inFile->beam[ibeam].bandData[iband].astro_data.pol4[j+idump*nchan];
		}
	    }
	  freq[n] = avF/av;
	  pol1[n] = av1/av;
	  if (npol > 1) pol2[n] = av2/av;
	  if (polPlot==2)
	    {
	      pol3[n] = av3/av;
	      pol4[n] = av4/av;
	    }
	  if (sump==1 && npol > 1)
	    pol1[n] = pol1[n]+pol2[n];
	  
	  if (setLog == 1)
	    {
	      pol1[n]=log10(pol1[n]);
	      if (npol > 1) pol2[n]=log10(pol2[n]);
	    }
	  if (freq[i] > bl_f0 && freq[i] <= bl_f1)
	    {
	      
	      mean1+=pol1[i];
	      mean2+=pol2[i];
	      mean3+=pol3[i];
	      mean4+=pol4[i];
	      nc++;
	    }
	  n++;
	}
      nchan = n;
    }

  if (setBaseline==1)
    {
      for (i=0;i<nchan;i++)
	{
	  pol1[i] -= mean1/(float)nc;
	  pol2[i] -= mean2/(float)nc;
	  pol3[i] -= mean3/(float)nc;
	  pol4[i] -= mean4/(float)nc;
	}
    }

  
  if (setf0 == 1)
    minx = f0;
  else
    minx = freq[0];
  
  if (setf1 == 1)
    maxx = f1;
  else
    maxx = freq[nchan-1];
  if (minx > maxx)
    {
      float tm = minx;
      minx = maxx;
      maxx = tm;
    }

  for (i=0;i<nchan;i++)
    {
      if (freq[i] >= minx && freq[i] <= maxx)
	{
	  if (useP1[i] == 1 && t==0)
	    {
	      miny = maxy = pol1[i];
	      t=1;
	    }
	  else
	    {
	      if (useP1[i] == 1 && pol1[i] > maxy) maxy = pol1[i];
	      if (useP1[i] == 1 && pol1[i] < miny) miny = pol1[i];
	      if (npol > 1)
		{
		  if (useP2[i] == 1 && pol2[i] > maxy) maxy = pol2[i];
		  if (useP2[i] == 1 && pol2[i] < miny) miny = pol2[i];
		}
	      if (polPlot==2)
		{
		  if (pol3[i] > maxy && useP1[i]==1) maxy = pol3[i];
		  if (pol4[i] > maxy && useP1[i]==1) maxy = pol4[i];
		  if (pol3[i] < miny && useP1[i]==1) miny = pol3[i];
		  if (pol4[i] < miny && useP1[i]==1) miny = pol4[i];		  
		}
	    }
	}
    }

  if (locky1 != locky2)
    {
      miny = locky1;
      maxy = locky2;
    }

  if (join==0)
    {
      if (setLog==-1)
	cpgenv(minx,maxx,miny,maxy,0,1);
      else
	cpgenv(minx,maxx,miny,maxy,0,20);
      sprintf(title,"%s: Spectral dump %d",inFile->beam[ibeam].bandHeader[iband].label,idump);

      // FIX ME
      /*
      sdhdf_loadFrequencyAttributes(inFile,iband);
      sdhdf_loadDataAttributes(inFile,iband);
      if (fref < 0)
	sprintf(xlabel,"%s frequency (%s)",inFile->frequency_attr.frame,inFile->frequency_attr.unit);
      else
	sprintf(xlabel,"%s velocity (km/s)",inFile->frequency_attr.frame);
      */
      if (fref < 0)
	sprintf(xlabel,"Frequency (MHz)");
      else
	sprintf(xlabel,"Velocity (km/s)");

      cpglab(xlabel,"Signal strength (arbitrary)",fname);
      
    }
  else
    {      
      float dh,dw;

      // NX must be 1 for join to work
      if (nx != 1)
	{
	  printf("WARNING: Plot not correct with -join if nx != 1\n");
	}
      if (entry==0)
	{
	  cpgsvp(wx1,wx2,wy1,wy2);
	  cpgswin(0,1,0,1);
	  sprintf(title,"%s: Spectral dump %d",inFile->beam[ibeam].bandHeader[iband].label,idump);

	  // FIX ME
	  //	  sdhdf_loadFrequencyAttributes(inFile,iband);
	  //	  sdhdf_loadDataAttributes(inFile,iband);

	  /*
	  if (fref < 0)
	    sprintf(xlabel,"%s frequency (%s)",inFile->frequency_attr.frame,inFile->frequency_attr.unit);
	  else
	    sprintf(xlabel,"%s velocity (km/s)",inFile->frequency_attr.frame);
	  */
	  if (fref < 0)
	    sprintf(xlabel,"frequency (MHz)");
	  else
	    sprintf(xlabel,"velocity (km/s)");

	  cpglab(xlabel,"Signal strength (arbitrary)","");
	}
      
      dh = (wy2-wy1)/ny;
      dw = (wx2-wx1)/nx;
      xp1 = wx1+dw*(entryX);
      xp2 = xp1 + dw;;
      yp1 = wy2-(dh)*(entryY+1);
      yp2 = wy2-(dh)*entryY;

      cpgsvp(xp1,xp2,yp1,yp2);
      cpgswin(minx,maxx,miny,maxy);
      if (entry >= nx*ny-nx)
	cpgbox("BCNTS",0,0,"BCTSN",0,0);
      else
	cpgbox("BCTS",0,0,"BCTSN",0,0);
      cpgtext(minx+(maxx-minx)*0.05,maxy-(maxy-miny)*0.15,inFile->fname);
    }

  if (molecularLines==1)
    {
      drawMolecularLine(1420.405752,"HI",minx,maxx,miny,maxy);
      drawMolecularLine(1612.2310,"OH",minx,maxx,miny,maxy);
      drawMolecularLine(1665.4018,"OH",minx,maxx,miny,maxy);
      drawMolecularLine(1667.3590,"OH",minx,maxx,miny,maxy);
      drawMolecularLine(1720.5300,"OH",minx,maxx,miny,maxy);
    }

  if (polPlot==1)
    {
      // Do not plot flagged channels
      cpgsci(1); drawIncludeWeights(nchan,freq,pol1,wt); cpgsci(1);
      if (sump==0 && npol > 1)
	{
	  cpgsci(2); drawIncludeWeights(nchan,freq,pol2,wt); cpgsci(1);
	}
    }
  else if (polPlot==2)
    {
      cpgsci(2); drawIncludeWeights(nchan,freq,pol2,wt); cpgsci(1);
      cpgsci(4); drawIncludeWeights(nchan,freq,pol3,wt); cpgsci(1);
      cpgsci(5); drawIncludeWeights(nchan,freq,pol4,wt); cpgsci(1);
      cpgsci(1); drawIncludeWeights(nchan,freq,pol1,wt); cpgsci(1);
    }

  free(wt);
  free(useP1);
  free(useP2);
  free(pol1);
  free(pol2);
  free(pol3);
  free(pol4);
  free(freq);
  sdhdf_releaseBandData(inFile,ibeam,iband,1);
  //  sdhdf_freeSpectrumMemory(&spectrum);
  /*
  for (i=0;i<inFile->bandHeader[iband].nchan;i++)
    {
      printf("%d %d %d %.5f %g %g %g %g %g %g %g %g\n",iband,idump,i,inFile->bandHeader[iband].topoFreq[i],
	     spectrum.pol1[i].val,spectrum.pol2[i].val,spectrum.pol3[i].val,spectrum.pol4[i].val,
	     spectrum.pol1[i].weight,spectrum.pol2[i].weight,spectrum.pol3[i].weight,spectrum.pol4[i].weight);
    }
  */  

  entry++;
  entryX++;
  if (entryX == nx)
    {
      entryX = 0;
      entryY ++;
    }
}

void drawIncludeWeights(int nchan,float *freq,float *pol,float *wt)
{
  int i,j,i0;
  int drawIt=-1;

  for (i=0;i<nchan;i++)
    {
      if (wt[i] != 0 && drawIt==-1)
	{
	  drawIt=1;
	  i0=i;
	}
      if (wt[i] == 0 && drawIt==1)
	{
	  cpgline(i-1-i0,freq+i0,pol+i0);
	  drawIt=-1;
	  
	}	      
    }
  if (drawIt==1) cpgline(i-1-i0,freq+i0,pol+i0); 
  
  //  cpgline(nchan,freq,pol);
}

void drawMolecularLine(float freq,char *label,float minX,float maxX,float minY,float maxY)
{
  float fx[2],fy[2];

  fx[0] = fx[1] = freq;
  fy[0] = minY;
  fy[1] = maxY;
  cpgsci(3); cpgsls(4); cpgline(2,fx,fy); cpgsls(1); cpgsci(1);

}

