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
// gcc -lm -o sdhdf_plotWide sdhdf_plotWide.c sdhdfProc.c -I../hdf5-1.10.4/src/ ../hdf5-1.10.4/src/.libs/libhdf5.a -ldl -lz -lcpgplot -L/pulsar/psr/software/stable/stretch/lib/ -I//pulsar/psr/software/stable/stretch/include -lcalceph -Isofa/20190722/c/src/ -Lsofa/20190722/c/src -lsofa_c

//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sdhdfProc.h"
#include "hdf5.h"
#include <cpgplot.h>

#define VNUM "v0.1"

void drawMolecularLine(float freq,char *label,float minX,float maxX,float minY,float maxY);
void plotSpectrum(sdhdf_fileStruct *inFile,int iband,int idump,char *grDev,char *fname,float f0,float f1,int av,int sump,int nx,int ny,int polPlot);
void showTransmitter(float freq,float bw,char *label,float miny,float maxy);

int main(int argc,char *argv[])
{
  int        i,j,k,l,ii;
  char       fname[MAX_FILES][64];
  sdhdf_fileStruct *inFile;
  int        idump,iband,ibeam;
  float f0=-1,f1=-1;
  int av=0;
  int sump=0;
  int nx=1;
  int ny=1;
  int polPlot=1;
  long nVals=0;
  float *px,*py1,*py2,*py1_max,*py2_max;
  int *psum;
  long writepos;
  int allocateMemory=0;
  float miny,maxy,val1,val2,minx,maxx;
  float ominy,omaxy,ominx,omaxx;
  float setMinX=-1;
  float setMaxX=-1;
  float plotTransmitters=-1;
  int plotMaxHold=-1;
  int molecularLines=-1;
  int inFiles;
  int sb=-1;
  float mx,my;
  char key;
  char grDev[128]="/xs";
  
  printf("Starting\n");
  // Defaults
  idump = iband = ibeam = 0;
  strcpy(grDev,"/xs");
  
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }

  sdhdf_initialiseFile(inFile);

  inFiles=0;
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-minx")==0)
	sscanf(argv[++i],"%f",&setMinX);
      else if (strcmp(argv[i],"-maxx")==0)
	sscanf(argv[++i],"%f",&setMaxX);
      else if (strcmp(argv[i],"-transmitters")==0)
	plotTransmitters=1;
      else if (strcmp(argv[i],"-maxhold")==0)
	plotMaxHold=1;
      else if (strcmp(argv[i],"-g")==0)
	strcpy(grDev,argv[++i]);
      else if (strcmp(argv[i],"-sb")==0)
	sscanf(argv[++i],"%d",&sb);
      else
	strcpy(fname[inFiles++],argv[i]);
    }
  
  printf("Have inFiles=  %d\n",inFiles);
  for (i=0;i<inFiles;i++)
    {
      printf("Loading %s\n",fname[i]);
      
      sdhdf_initialiseFile(inFile);
      sdhdf_openFile(fname[i],inFile,1);

      sdhdf_loadMetaData(inFile);
      
      if (allocateMemory ==0)
	{
	  nVals=0;
	  for (j=0;j<inFile->beam[ibeam].nBand;j++)
	    {
	      sdhdf_loadBandData(inFile,ibeam,j,1);
	      nVals+=(inFile->beam[ibeam].bandHeader[j].nchan);
	    }
	  px = (float *)malloc(sizeof(float)*nVals);
	  py1 = (float *)calloc(sizeof(float),nVals);
	  py2 = (float *)calloc(sizeof(float),nVals);
	  py1_max = (float *)calloc(sizeof(float),nVals);
	  py2_max = (float *)calloc(sizeof(float),nVals);
	  psum = (int *)calloc(sizeof(int),nVals);
	  allocateMemory=1;
	  writepos=0;
	  for (j=0;j<inFile->beam[ibeam].nBand;j++)
	    {
	      for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k++)
		px[writepos+k] = inFile->beam[ibeam].bandData[j].astro_data.freq[k];
	      writepos+=inFile->beam[ibeam].bandHeader[j].nchan;
	    }
	}
      printf("Total number of values to plot =  %d\n",(int)nVals);
      writepos=0;
      for (j=0;j<inFile->beam[ibeam].nBand;j++)
	{
	  for (k=0;k<inFile->beam[ibeam].bandHeader[j].nchan;k++)
	    {
	      for (l=0;l<inFile->beam[ibeam].bandHeader[j].ndump;l++)
		{
		  val1 = inFile->beam[ibeam].bandData[j].astro_data.pol1[l*inFile->beam[ibeam].bandHeader[j].nchan+k];
		  val2 = inFile->beam[ibeam].bandData[j].astro_data.pol2[l*inFile->beam[ibeam].bandHeader[j].nchan+k];
		  py1[writepos+k] += val1;
		  py2[writepos+k] += val2;
		  if (val1 > py1_max[writepos+k]) py1_max[writepos+k] = val1;
		  if (val2 > py2_max[writepos+k]) py2_max[writepos+k] = val2;
		  psum[writepos+k] ++;
		}
	    }
	  writepos+=inFile->beam[ibeam].bandHeader[j].nchan;
	}
      

      
      printf("Closing\n");
      sdhdf_closeFile(inFile);
      printf("Done close\n");
    }

  miny = 1e99; maxy = -1e99;
  
  for (i=0;i<nVals;i++)
    {
      py1[i]/=(float)psum[i];
      py2[i]/=(float)psum[i];

      py1[i] = log10(py1[i]);
      py2[i] = log10(py2[i]);

      py1_max[i] = log10(py1_max[i]);
      py2_max[i] = log10(py2_max[i]);
      
      if (py1[i] > maxy) maxy = py1[i];
      if (py2[i] > maxy) maxy = py2[i];
      if (py1_max[i] > maxy) maxy = py1_max[i];
      if (py2_max[i] > maxy) maxy = py2_max[i];
      if (py1[i] < miny) miny = py1[i];
      if (py2[i] < miny) miny = py2[i];
    }

  if (sb > -1)
    {
      setMinX = 704+sb*128;
      setMaxX = 704+(sb+1)*128;
    }
  
  
  if (setMinX >= 0)
    minx = setMinX;
  else
    minx = 700;

  if (setMaxX >= 0)
    maxx = setMaxX;
  else
    maxx = 4040;

  ominy = miny;
  omaxy = maxy;
  ominx = minx;
  omaxx = maxx;
  
  cpgbeg(0,grDev,1,1);  
  cpgsch(1.4);
  cpgslw(2);
  cpgscf(2);
  cpgask(0);

  do {
    cpgenv(minx,maxx,miny,maxy,0,20);
    cpglab("Observing frequency (MHz)","Signal strength","");
    cpgsci(2);
    cpgline(nVals,px,py1);
    if (plotMaxHold==1)
      {
	cpgsls(2); cpgline(nVals,px,py1_max); cpgsls(1);
      }
    cpgsci(5);
    cpgline(nVals,px,py2);
    if (plotMaxHold==1)
      {
	cpgsls(2); cpgline(nVals,px,py2_max); cpgsls(1);
      }
    cpgsci(1);
    
    if (molecularLines==1)
      {
	drawMolecularLine(1420.405752,"HI",minx,maxx,miny,maxy);
	drawMolecularLine(1612.2310,"OH",minx,maxx,miny,maxy);
	drawMolecularLine(1665.4018,"OH",minx,maxx,miny,maxy);
	drawMolecularLine(1667.3590,"OH",minx,maxx,miny,maxy);
	drawMolecularLine(1720.5300,"OH",minx,maxx,miny,maxy);
      }

    if (plotTransmitters==1)
      {
	cpgsch(0.8);
	showTransmitter(763,10,"Optus",miny,maxy);
	showTransmitter(778,20,"Telstra",miny,maxy);
	showTransmitter(847.8,0.4,"Parkes radio",miny,maxy);
	showTransmitter(849.5,0.2,"Police",miny,maxy);
	showTransmitter(872.5,5,"Vodafone",miny,maxy);
	showTransmitter(882.5,14.9,"Telstra",miny,maxy);
	showTransmitter(905,3,"Amateur radio/walkie talkie",miny,maxy);
	showTransmitter(930.5,0.2,"Police",miny,maxy);
	showTransmitter(947.6,8.4,"Optus",miny,maxy);
	showTransmitter(956.2,5,"Vodafone",miny,maxy);
	showTransmitter(1018,0.5,"PKS airport",miny,maxy);
	showTransmitter(1575.42,32.736,"Beidou",miny,maxy);
	showTransmitter(1176.45,20.46,"Beidou",miny,maxy);
	showTransmitter(1268.52,20.46,"Beidou",miny,maxy);
	showTransmitter(1542,34,"Inmarsat",miny,maxy);
	showTransmitter(1643.25,33.5,"Inmarsat",miny,maxy);
	showTransmitter(1575,32,"Galileo",miny,maxy);
	showTransmitter(1280,40,"Galileo",miny,maxy);
	showTransmitter(1189,50,"Galileo",miny,maxy);
	showTransmitter(1575.42,15.456,"GPS",miny,maxy);
	showTransmitter(1227.6,11,"GPS",miny,maxy);
	showTransmitter(1176.45,12.5,"GPS",miny,maxy);
	showTransmitter(1626.25,0.5,"Iridium",miny,maxy);
	showTransmitter(1542.5,35,"Thurya",miny,maxy);
	showTransmitter(1597.21875,16.3125,"Glonass",miny,maxy);
	showTransmitter(1245.78125,5.6875,"Glonass",miny,maxy);
	showTransmitter(1815,20,"Telstra",miny,maxy);
	showTransmitter(1835,20,"Telstra",miny,maxy);
	showTransmitter(1857.5,25,"Optus",miny,maxy);
	showTransmitter(2115,10,"Vodafone",miny,maxy);
	showTransmitter(2167.5,4.8,"Vodafone",miny,maxy);
	showTransmitter(2122.5,5,"Telstra",miny,maxy);
	showTransmitter(2127.5,5,"Telstra",miny,maxy);
	showTransmitter(2142.5,5,"Optus",miny,maxy);
	showTransmitter(2147.5,5,"Optus",miny,maxy);
	showTransmitter(2152.5,5,"Optus",miny,maxy);
	showTransmitter(2160.0,9.9,"Telstra",miny,maxy);
	showTransmitter(2312,19.2,"NBN",miny,maxy);
	showTransmitter(2331.2,19.2,"NBN",miny,maxy);
	showTransmitter(2350.4,19.2,"NBN",miny,maxy);
	showTransmitter(2369.6,19.2,"NBN",miny,maxy);
	showTransmitter(2386.7,15,"NBN",miny,maxy);

	showTransmitter(2650.0,40,"Telstra",miny,maxy);
	showTransmitter(2680.0,20,"Optus",miny,maxy);

	showTransmitter(3455.0,19.9,"NBN",miny,maxy);
	showTransmitter(3560.0,19.9,"NBN",miny,maxy);
	
	showTransmitter(921.5,13,"WiFi 802.11",miny,maxy);
	showTransmitter(2442.0,82,"WiFi",miny,maxy);

	showTransmitter(1090,50./1000.,"ADS-B",miny,maxy);

	cpgsch(1.4);
      }

    if (strcmp(grDev,"/xs")==0)
      {
	cpgcurs(&mx,&my,&key);
	
	if (key=='z')
	  {
	    float mx2,my2;
	    cpgband(2,0,mx,my,&mx2,&my2,&key);
	    printf("Zooming in on %g %g %g %g\n",mx,my,mx2,my2);
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
	else if (key=='m')
	  plotMaxHold*=-1;
	else if (key=='M')
	  molecularLines*=-1;
	else if (key=='t')
	  plotTransmitters*=-1;
	else if (key=='u')
	  {
	    miny = ominy;
	    maxy = omaxy;
	    minx = ominx;
	    maxx = omaxx;
	  }
      }
    else
      key='q';

  } while (key!='q');

  
  cpgend();  
  free(px);
  free(py1);
  free(py2);
  free(py1_max);
  free(py2_max);      
  free(psum);
  free(inFile);
}




void showTransmitter(float freq,float bw,char *label,float miny,float maxy)
{
  float fx[2],fy[2];
  cpgptxt(freq,maxy-(maxy-miny)*0.05,0,0.5,label);
  cpgsls(4);
  fx[0] = fx[1] = freq-bw/2;
  fy[0] = miny;  fy[1] = maxy;
  cpgline(2,fx,fy);
  fx[0] = fx[1] = freq+bw/2;
  fy[0] = miny;  fy[1] = maxy;
  cpgline(2,fx,fy);
  cpgsls(1);
}

void drawMolecularLine(float freq,char *label,float minX,float maxX,float minY,float maxY)
{
  float fx[2],fy[2];

  fx[0] = fx[1] = freq;
  fy[0] = minY;
  fy[1] = maxY;
  cpgsci(3); cpgsls(4); cpgline(2,fx,fy); cpgsls(1); cpgsci(1);

}

