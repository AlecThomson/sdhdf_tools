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
// Software to process and view the calibrator
//
// Usage:
// sdhdf_cal
//
// Compilation
// gcc -lm -o sdhdf_cal sdhdf_cal.c sdhdfProc.c -I../hdf5-1.10.4/src/ ../hdf5-1.10.4/src/.libs/libhdf5.a -ldl -lz -lcpgplot -L/pulsar/psr/software/stable/stretch/lib/ -I//pulsar/psr/software/stable/stretch/include -lcalceph -Isofa/20190722/c/src/ -Lsofa/20190722/c/src -lsofa_c 
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sdhdfProc.h"
#include "hdf5.h"
#include <cpgplot.h>


void doPlot(sdhdf_fileStruct *inFile,int beam,int totChan);

int main(int argc,char *argv[])
{
  int i,j;
  char fname[MAX_STRLEN];
  sdhdf_fileStruct *inFile;
  int iband=0,nchan,idump;
  int ndump=0;
  int beam=0;
  int totChan=0;
  
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }

  sdhdf_initialiseFile(inFile);
    
  for (i=1;i<argc;i++)
    {      
      if (strcmp(argv[i],"-f")==0)
	strcpy(fname,argv[++i]);	
    }
  
  sdhdf_openFile(fname,inFile,1);
  printf("File opened\n");
  sdhdf_loadMetaData(inFile);



  ndump = inFile->beam[beam].calBandHeader[0].ndump;
  printf("ndumps = %d\n",ndump);


  for (i=0;i<inFile->beam[beam].nBand;i++)
    {
      nchan = inFile->beam[beam].calBandHeader[i].nchan;
      totChan+=nchan;
      printf("Cal nchan = %d\n",nchan);
      sdhdf_loadBandData(inFile,beam,i,2);
      sdhdf_loadBandData(inFile,beam,i,3);
    }
  doPlot(inFile,beam,totChan);

  sdhdf_closeFile(inFile);
  free(inFile);
}

void doPlot(sdhdf_fileStruct *inFile,int beam,int totChan)
{
  float fx[totChan];
  float fy1[totChan];
  float fy2[totChan];
  float fy3[totChan];
  float fy4[totChan];
  int np=0;
  int i,j,k;
  float mx,my;
  char key;
  int idump=0;
  int nchan;
  int plot=1;
  int nLine=2;
  float onP1,offP1,onP2,offP2;
  float onP3,offP3,onP4,offP4;
  float miny,maxy;
  float uminx,umaxx;
  float uminy,umaxy;
  float set_miny=-1,set_maxy=-1;
  float set_minx=-1,set_maxx=-1;
  int log=0;
  
  cpgbeg(0,"/xs",1,1);
  cpgask(0);
  cpgslw(2);
  cpgscf(2);
  cpgsch(1.4);
  do {
    np=0;
    for (i=0;i<inFile->beam[beam].nBand;i++)
      {
	nchan = inFile->beam[beam].calBandHeader[i].nchan;
	for (j=0;j<nchan;j++)
	  {
	    fx[np] = inFile->beam[beam].bandData[i].cal_on_data.freq[j];
	    onP1  = inFile->beam[beam].bandData[i].cal_on_data.pol1[j+idump*nchan];
	    offP1 = inFile->beam[beam].bandData[i].cal_off_data.pol1[j+idump*nchan];
 	    onP2  = inFile->beam[beam].bandData[i].cal_on_data.pol2[j+idump*nchan];
	    offP2 = inFile->beam[beam].bandData[i].cal_off_data.pol2[j+idump*nchan];
 	    onP3  = inFile->beam[beam].bandData[i].cal_on_data.pol3[j+idump*nchan];
	    offP3 = inFile->beam[beam].bandData[i].cal_off_data.pol3[j+idump*nchan];
 	    onP4  = inFile->beam[beam].bandData[i].cal_on_data.pol4[j+idump*nchan];
	    offP4 = inFile->beam[beam].bandData[i].cal_off_data.pol4[j+idump*nchan];

	    if (plot==1)
	      {
		fy1[np] = onP1;
		fy2[np] = onP2;
		
		fy1[np] = log10(fy1[np]);
		fy2[np] = log10(fy2[np]);
		nLine=2;
		log=1;
	      }
	    else if (plot==2)
	      {
		fy1[np] = onP1;
		fy2[np] = offP1;
		fy1[np] = log10(fy1[np]);
		fy2[np] = log10(fy2[np]);
		nLine=2;
		log=1;
	      }
	    else if (plot==3)
	      {
		fy1[np] = offP1/(onP1-offP1);
		fy2[np] = offP2/(onP2-offP2);
		printf("Have %g %g\n",fy1[np],fy2[np]);
		nLine=2;
		log=0;
	      }
	    else if (plot==4)
	      {
		float ical = (onP1-offP1) + (onP2-offP2);
		float qcal = (onP1-offP1) - (onP2-offP2);

		fy1[np] = 2*qcal/ical;
		printf("Have %g %g\n",fy1[np],fy2[np]);
		nLine=1;
		log=0;
	      }
	    else if (plot==5)
	      {
		float ucal = 2*(onP3-offP3);
		float vcal = 2*(onP4-offP4);
		fy1[np] = atan2(vcal,ucal)*180/M_PI;
		printf("Have %g %g\n",fy1[np],fy2[np]);
		nLine=1;
		log=0;
	      }
	    np++;
	  }
      }
    for (i=0;i<np;i++)
      {
	if (i==0)
	  {
	    miny = maxy = fy1[i];
	  }
	if (miny > fy1[i]) miny = fy1[i];
	if (maxy < fy1[i]) maxy = fy1[i];
	if (nLine > 1)
	  {
	    if (miny > fy2[i]) miny = fy2[i];
	    if (maxy < fy2[i]) maxy = fy2[i];
	  }
	  }
    if (set_miny==-1 && set_maxy==-1)
      {
	uminy = miny;
	umaxy = maxy;
	uminx = 704;
	umaxx = 4032;
      }
    else
      {
	uminy = set_miny;
	umaxy = set_maxy;
	uminx = set_minx;
	umaxx = set_maxx;
      }
    if (log==1)
      cpgenv(uminx,umaxx,uminy,umaxy,0,20);
    else
      cpgenv(uminx,umaxx,uminy,umaxy,0,1);
    if (plot==3)
      cpglab("Frequency (MHz)","OFF/(ON-OFF)","");
    else
      cpglab("Frequency (MHz)","","");

    cpgsci(1); cpgline(np,fx,fy1);
    if (nLine > 1)
      {
	cpgsci(2); cpgline(np,fx,fy2);
      }
    
    cpgsci(1);


    cpgcurs(&mx,&my,&key);
    if (key=='1') {plot=1; set_miny = set_maxy = -1;}
    else if (key=='2') {plot=2;set_miny = set_maxy = -1;}
    else if (key=='3') {plot=3;set_miny = set_maxy = -1;}
    else if (key=='4') {plot=4;set_miny = set_maxy = -1;}
    else if (key=='5') {plot=5;set_miny = set_maxy = -1;}
    else if (key=='u')
      {
	set_miny = set_maxy = -1;
      }
    else if (key=='z')
      {
	float mx2,my2;
	cpgband(2,0,mx,my,&mx2,&my2,&key);
	if (mx != mx2 && my != my2)
	  {
	    if (mx < mx2)
	      {set_minx = mx; set_maxx = mx2;}
	    else
	      {set_minx = mx2; set_maxx = mx;}
	    
	    if (my < my2)
	      {set_miny = my; set_maxy = my2;}
	    else
	      {set_miny = my2; set_maxy = my;}
	  }
	else
	  printf("Please press 'z' and then move somewhere and click left mouse button\n");
	
      }
  } while (key!='q');

  
}
