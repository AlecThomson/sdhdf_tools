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
// gcc -lm -o sdhdf_plotSpectrum sdhdf_plotSpectrum.c sdhdfProc.c -I../hdf5-1.10.4/src/ ../hdf5-1.10.4/src/.libs/libhdf5.a -ldl -lz -lcpgplot -L/pulsar/psr/software/stable/stretch/lib/ -I//pulsar/psr/software/stable/stretch/include -lcalceph -Isofa/20190722/c/src/ -Lsofa/20190722/c/src -lsofa_c

//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sdhdfProc.h"
#include "hdf5.h"
#include <cpgplot.h>

void plotSpectrum(sdhdf_fileStruct *inFile,int ibeam,int iband,int idump,double fref,char *yUnit,char *freqFrame,char *freqUnit);
void drawMolecularLine(float freq,char *label,float minX,float maxX,float minY,float maxY);


#define VNUM "v0.1"

void help()
{
  printf("sdhdf_plotSpectrum: %s\n",VNUM);
  printf("sdhfProc version:   %s\n",SOFTWARE_VER);
  printf("author:             George Hobbs\n");
  printf("\n");
  printf("Software to plot a spectrum in an interactive matter\n");
  printf("\n\nCommand line arguments:\n\n");
  printf("-h                  This help\n");
  printf("-f <filename>       SDHDF file corresponding to observations\n");
  
  printf("\nExample:\n\n");
  printf("sdhdf_plotSpectrum -f diff.hdf -sb 0\n");
  printf("---------------------\n");
}


int main(int argc,char *argv[])
{
  int i,j;
  int nAttr;
  char fname[MAX_STRLEN]="unset";
  sdhdf_fileStruct *inFile;
  double fref=-1;
  char yUnit[MAX_STRLEN] = "not set";
  char freqFrame[MAX_STRLEN] = "[unknown]";
  char freqUnit[MAX_STRLEN] = "unknown";
  char att_yUnit[MAX_STRLEN] = "arbitrary";
  int idump,iband;
  char dataName[MAX_STRLEN];
  int ibeam=0;
  
  help();
  
  // Defaults
  idump = iband = 0;
  
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }

  for (i=1;i<argc;i++)
    {      
      if (strcmp(argv[i],"-f")==0)
	strcpy(fname,argv[++i]);	
      else if (strcmp(argv[i],"-yunit")==0)
	strcpy(yUnit,argv[++i]);	
      else if (strcmp(argv[i],"-fref")==0)
	sscanf(argv[++i],"%lf",&fref);
      else if (strcmp(argv[i],"-h")==0)
	{help(); exit(1);}
      else if (strcmp(argv[i],"-sd")==0)
	sscanf(argv[++i],"%d",&idump);
      else if (strcmp(argv[i],"-sb")==0) // FIX ME - ENABLE BAND DESCRIPTORS
	sscanf(argv[++i],"%d",&iband);
    }

  sdhdf_initialiseFile(inFile);
    
  if (strcmp(fname,"unset")==0)
    {
      printf("Please define an input file name using the -f option\n");
      free(inFile);
      exit(1);
    }
  
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

  // FIX THIS
  /*
    for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-sb")==0)
	iband = sdhdf_getBandID(inFile,argv[++i]);
    }
  */
  
  // FIX THESE
  //  sdhdf_loadFrequencyAttributes(inFile,iband);
  //  sdhdf_loadDataAttributes(inFile,iband);

  //  strcpy(freqFrame,inFile->frequency_attr.frame);
  //  strcpy(freqUnit,inFile->frequency_attr.unit);
  //  strcpy(att_yUnit,inFile->data_attr.unit);
  
  if (strcmp(yUnit,"not set")==0)
    strcpy(yUnit,att_yUnit);
  
  plotSpectrum(inFile,ibeam,iband,idump,fref,yUnit,freqFrame,freqUnit);

  
  sdhdf_closeFile(inFile);
  free(inFile);
}

void plotSpectrum(sdhdf_fileStruct *inFile,int ibeam,int iband,int idump,double fref,char *yUnit,char *freqFrame,char *freqUnit)
{
  char key;
  float mx,my;
  int i,j;
  float *aa,*bb,*ab,*abs;
  float minx,maxx,miny,maxy,minz,maxz;
  float ominx,omaxx,ominy,omaxy;
  int t=0;
  int setLog=-1;
  char title[1024];
  char ylabel[1024];
  char xlabel[1024];
  char legend[1024];
  int plot=1;
  int molecularLines=-1;
  int nchan,maxNchan;
  float *pol1,*pol2,*pol3,*pol4;
  float *freq;
  int xplot=1;
  int npol=4;
  int reload=1;
  int flagIt=2;
  
  cpgbeg(0,"/xs",1,1);
  cpgask(0);
  cpgscf(2);
  cpgsch(1.4);
  cpgslw(2);

  maxNchan=0;
  for (i=0;i<inFile->beam[ibeam].nBand;i++)
    {
      if (inFile->beam[ibeam].bandHeader[i].nchan > maxNchan)
	maxNchan = inFile->beam[ibeam].bandHeader[i].nchan;
    }
  pol1 = (float *)malloc(sizeof(float)*maxNchan);
  pol2 = (float *)malloc(sizeof(float)*maxNchan);
  pol3 = (float *)malloc(sizeof(float)*maxNchan);
  pol4 = (float *)malloc(sizeof(float)*maxNchan);
  freq = (float *)malloc(sizeof(float)*maxNchan);

  do
    {
      if (plot==1)
	{
	  nchan = inFile->beam[ibeam].bandHeader[iband].nchan;
	  npol  = inFile->beam[ibeam].bandHeader[iband].npol;

	  if (reload==1) // Should reload if band or beam changes
	    {
	      printf("Loading band data %d\n",iband);
	      // If already loaded then should release data **
	      // SHOULD ONLY LOAD IF NOT LOADED YET
	      if (inFile->beam[ibeam].bandData[iband].astro_data.pol1AllocatedMemory == 0)
		sdhdf_loadBandData(inFile,ibeam,iband,1);
	      reload=0;
	    }
	  for (i=0;i<nchan;i++)
	    {
	      if (xplot==1)
		{
		  if (fref < 0)
		      freq[i] = inFile->beam[ibeam].bandData[iband].astro_data.freq[i];
		  else
		    {
		      freq[i] = (1.0-inFile->beam[ibeam].bandData[iband].astro_data.freq[i]/(fref))*SPEED_LIGHT/1000.; // km/s
		    }
		}
	      else
		freq[i] = i;

	      pol1[i] = inFile->beam[ibeam].bandData[iband].astro_data.pol1[i+idump*nchan];
	      if (npol > 1) pol2[i] = inFile->beam[ibeam].bandData[iband].astro_data.pol2[i+idump*nchan];
	      if (setLog == 1)
		{
		  pol1[i]=log10(pol1[i]);
		  if (npol > 1) pol2[i]=log10(pol2[i]);
		}
	    }

	
	  if (t==0 || t==1)
	    {
	      int setMiny=0;
	      for (i=0;i<nchan;i++)
		{
		  if (inFile->beam[ibeam].bandData[iband].astro_data.flag[i] == 0 || flagIt==0)
		    {
		      if (setMiny==0)
			{
			  miny = maxy = pol1[i];
			  setMiny=1;
			}
		      else
			{
			  if (pol1[i] > maxy) maxy = pol1[i];
			  if (pol1[i] < miny) miny = pol1[i];
			  
			  if (npol > 1)
			    {
			      if (pol2[i] > maxy) maxy = pol2[i];
			      if (pol2[i] < miny) miny = pol2[i];
			    }
			}
		    }
		}
	      if (t==0)
		{
		  if (xplot==1)
		    {
		      minx = freq[0]; maxx = freq[nchan-1];
		    }
		  else
		    {
		      minx = 0;
		      maxx = nchan;
		    }
		      ominx = minx; omaxx = maxx;
		}
	      ominy = miny; omaxy = maxy;
	      t=2;
	      
	    }
	}

      
      if (setLog==-1)
	cpgenv(minx,maxx,miny,maxy,0,0);
      else
	cpgenv(minx,maxx,miny,maxy,0,20);
      sprintf(ylabel,"Signal strength (%s)",yUnit);
      if (inFile->beam[ibeam].bandHeader[iband].ndump==1)       
	sprintf(title,"%s, %s",inFile->fname,inFile->beam[ibeam].bandHeader[iband].label);
      else
	sprintf(title,"%s, %s, spectral dump %d",inFile->fname,inFile->beam[ibeam].bandHeader[iband].label,idump);
      if (fref < 0)
	{
	  sprintf(xlabel,"%s frequency (%s)",freqFrame,freqUnit);
	}
      else
	{
	  sprintf(xlabel,"Velocity (km s\\u-1\\d) [f\\dref\\u = %.2f MHz]",fref);
	}
      cpglab(xlabel,ylabel,title);
	  
      cpgsch(0.9);
      sprintf(legend,"Source: %s",inFile->beamHeader[ibeam].source);
      cpgtext(minx+0.05*(maxx-minx),maxy-0.05*(maxy-miny),legend);
      sprintf(legend,"Tobs: %.2f sec",inFile->beam[ibeam].bandHeader[iband].dtime);
      cpgtext(minx+0.05*(maxx-minx),maxy-0.1*(maxy-miny),legend);
      sprintf(legend,"RA:  %s",inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].raStr);
      cpgtext(minx+0.05*(maxx-minx),maxy-0.15*(maxy-miny),legend);
      sprintf(legend,"Dec: %s",inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].decStr);
      cpgtext(minx+0.05*(maxx-minx),maxy-0.2*(maxy-miny),legend);
      cpgsch(1.4);
      
      if (molecularLines==1)
	{
	  drawMolecularLine(1420.405752,"HI",minx,maxx,miny,maxy);
	  drawMolecularLine(1612.2310,"OH",minx,maxx,miny,maxy);
	  drawMolecularLine(1665.4018,"OH",minx,maxx,miny,maxy);
	  drawMolecularLine(1667.3590,"OH",minx,maxx,miny,maxy);
	  drawMolecularLine(1720.5300,"OH",minx,maxx,miny,maxy);
	}
      if (flagIt==2)
	{
	  int i0=0;
	  int drawIt=-1;
	  for (i=0;i<nchan;i++)
	    {
	      if (inFile->beam[ibeam].bandData[iband].astro_data.flag[i] == 0 && drawIt==-1)
		{
		  drawIt=1;
		  i0=i;
		}
	      if (inFile->beam[ibeam].bandData[iband].astro_data.flag[i] != 0 && drawIt==1)
		{
		  cpgsci(1); cpgline(i-1-i0,freq+i0,pol1+i0); cpgsci(1);
		  if (npol > 1)
		    {
		      cpgsci(2); cpgline(i-1-i0,freq+i0,pol2+i0); cpgsci(1);
		    }
		  drawIt=-1;

		}	      
	    }
	  if (drawIt==1)
	    {
	      cpgsci(1); cpgline(i-1-i0,freq+i0,pol1+i0); cpgsci(1);
	      if (npol > 1)
		{
		  cpgsci(2); cpgline(i-1-i0,freq+i0,pol2+i0); cpgsci(1);
		}
	    }
	}
      else
	{
	  cpgsci(1); cpgline(nchan,freq,pol1); cpgsci(1);
	  if (npol > 1)
	    {
	      cpgsci(2); cpgline(nchan,freq,pol2); cpgsci(1);
	    }
	}
      cpgcurs(&mx,&my,&key);
      if (key=='z')
	{
	  float mx2,my2;
	  cpgband(2,0,mx,my,&mx2,&my2,&key);
	  if (mx != mx2 && my != my2)
	    {
	      if (mx < mx2)
		{minx = mx; maxx = mx2;}
	      else
		{minx = mx2; maxx = mx;}
	      
	      if (my < my2)
		{miny = my; maxy = my2;}
	      else
		{miny = my2; maxy = my;}
	    }
	  else
	    printf("Please press 'z' and then move somewhere and click left mouse button\n");
	}
      else if (key=='f')
	{
	  flagIt++;
	  if (flagIt==3) flagIt=0;
	  t=0;
	}
      else if (key=='A')
	printf("Mouse cursor = (%g,%g)\n",mx,my);
      else if (key=='1') {plot=1; t=0;}
      else if (key=='2') {plot=2; t=1;}
      else if (key=='l')
	{setLog*=-1; t=0;}
      else if (key=='>')
	{	  
	  iband++;
	  if (iband >= inFile->beam[ibeam].nBand)
	    iband = inFile->beam[ibeam].nBand-1;
	  reload=1;
	  t=0;
	  
	}
      else if (key=='<')
	{
	  iband--;
	  if (iband < 0) iband = 0;
	  reload=1;
	  t=0;
	}
      else if (key=='+')
	{
	  if (idump < inFile->beam[ibeam].bandHeader[iband].ndump-1)
	    idump++;
	}
      else if (key=='x')
	{
	  xplot*=-1;
	  t=0;
	}
      else if (key=='m')
	molecularLines*=-1;
      else if (key=='-')
	{
	  if (idump > 0)
	    idump--;
	}
      else if (key=='o')
	{
	  char fname[1024];
	  FILE *fout;
	  printf("Enter output filename ");
	  scanf("%s",fname);
	  fout = fopen(fname,"w");
	  for (i=0;i<nchan;i++)
	    {
	      if (freq[i] > minx && freq[i] <= maxx)
		{
		  if (npol > 1)
		    fprintf(fout,"%.6f %g %g\n",freq[i],pol1[i],pol2[i]);
		  else
		    fprintf(fout,"%.6f %g\n",freq[i],pol1[i]);
		}
	    }
	  fclose(fout);
	}
    else if (key=='u')
      {
	minx = ominx;
	maxx = omaxx;
	miny = ominy;
	maxy = omaxy;
      }
  } while (key!='q');
  
  cpgend();


  free(pol1);
  free(pol2);
  free(pol3);
  free(pol4);

  
  /*
  for (i=0;i<inFile->bandHeader[iband].nchan;i++)
    {
      printf("%d %d %d %.5f %g %g %g %g %g %g %g %g\n",iband,idump,i,inFile->bandHeader[iband].topoFreq[i],
	     spectrum.pol1[i].val,spectrum.pol2[i].val,spectrum.pol3[i].val,spectrum.pol4[i].val,
	     spectrum.pol1[i].weight,spectrum.pol2[i].weight,spectrum.pol3[i].weight,spectrum.pol4[i].weight);
    }
  */  

}

void drawMolecularLine(float freq,char *label,float minX,float maxX,float minY,float maxY)
{
  float fx[2],fy[2];

  fx[0] = fx[1] = freq;
  fy[0] = minY;
  fy[1] = maxY;
  cpgsci(3); cpgsls(4); cpgline(2,fx,fy); cpgsls(1); cpgsci(1);

}

