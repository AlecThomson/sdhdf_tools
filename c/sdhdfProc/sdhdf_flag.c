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
// Software to flag data
//
// Usage:
// sdhdf_cal
//
// Compilation
// gcc -lm -o sdhdf_flag sdhdf_flag.c sdhdfProc.c -I../hdf5-1.10.4/src/ ../hdf5-1.10.4/src/.libs/libhdf5.a -ldl -lz -lcpgplot -L/pulsar/psr/software/stable/stretch/lib/ -I//pulsar/psr/software/stable/stretch/include -lcalceph -Isofa/20190722/c/src/ -Lsofa/20190722/c/src -lsofa_c 
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sdhdfProc.h"
#include "hdf5.h"
#include <cpgplot.h>
#include "TKfit.h"

void doPlot(sdhdf_fileStruct *inFile,int ibeam);
void saveFile(sdhdf_fileStruct *inFile,int ibeam);

int main(int argc,char *argv[])
{
  int i,j;
  char fname[MAX_STRLEN];
  sdhdf_fileStruct *inFile;
  int iband=0,ibeam=0,nchan,idump,nband;
  int ndump=0,totSize,chanPos;
  //  spectralDumpStruct spectrum;
  int npol=4;
  int setnband=-1;
  
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }

  sdhdf_initialiseFile(inFile);
    
  for (i=1;i<argc;i++)
    {      
      if (strcmp(argv[i],"-f")==0)	       strcpy(fname,argv[++i]);	
      else if (strcmp(argv[i],"-setnband")==0) sscanf(argv[++i],"%d",&setnband);
      else if (strcmp(argv[i],"-sb")==0)       sscanf(argv[++i],"%d",&iband);

    }
  
  sdhdf_openFile(fname,inFile,1);
  printf("File opened\n");
  sdhdf_loadMetaData(inFile);
  
  if (setnband > 0)
    inFile->beam[ibeam].nBand = setnband;
  

  // FIX ME
  //  sdhdf_loadFlagData(inFile);
  
  totSize=0;
  nband = inFile->beam[ibeam].nBand;
   
  npol=4; // FIX HARDCODING OF NPOL
  for (i=0;i<nband;i++)
    {
      totSize+=(inFile->beam[ibeam].bandHeader[i].nchan*inFile->beam[ibeam].bandHeader[i].ndump*npol);
    }
  

  for (i=0;i<nband;i++)
    {
      nchan = inFile->beam[ibeam].bandHeader[i].nchan;
      ndump = inFile->beam[ibeam].bandHeader[i].ndump;
      sdhdf_loadBandData(inFile,ibeam,i,1);
    }

  doPlot(inFile,ibeam);

  sdhdf_closeFile(inFile);

  free(inFile);
}

void doPlot(sdhdf_fileStruct *inFile,int ibeam)
{
  int plot=1,i,j,k;
  float *px,*py1,*py2,*pf1,*pf2;
  float pointX[2],pointY[2];
  float maxPosX;
  float mx,my,mx2,my2;
  float minx,maxx,miny,maxy;
  char key;
  char title[MAX_STRLEN];
  int iband=0;
  int nchan,ndump;
  long max_nchan;
  int npol=4;
  int recalc=1;
  int nFlagRegion=0;
  int n;
  int c0;
  int regionType=1;
  
  cpgbeg(0,"/xs",1,1);
  cpgslw(2);
  cpgscf(2);
  cpgsch(1.4);
  cpgask(0);
  printf("Doing plot\n");


  max_nchan=0;
  for (i=0;i<inFile->beam[ibeam].nBand;i++)
    {
      if (max_nchan < inFile->beam[ibeam].bandHeader[i].nchan)
	max_nchan = inFile->beam[ibeam].bandHeader[i].nchan;
    }
  printf("Maximum channels = %d\n",max_nchan);
  
  px  = (float *)malloc(sizeof(float)*max_nchan);
  py1 = (float *)malloc(sizeof(float)*max_nchan);
  py2 = (float *)malloc(sizeof(float)*max_nchan);
  pf1 = (float *)malloc(sizeof(float)*max_nchan);
  pf2 = (float *)malloc(sizeof(float)*max_nchan);
  
  do {    
    nchan = inFile->beam[ibeam].bandHeader[iband].nchan;
    ndump = inFile->beam[ibeam].bandHeader[iband].ndump;

    if (recalc > 0)
      {
	regionType=1;
	nFlagRegion=0;
	c0 = 0;
	for (i=0;i<iband;i++)
	  c0+=(inFile->beam[ibeam].bandHeader[iband].nchan*inFile->beam[ibeam].bandHeader[iband].ndump*npol);

	if (recalc!=2) {minx = 1e99; maxx = -1e99;}
	miny = 1e99; maxy = -1e99;

	for (i=0;i<nchan;i++)
	  {
	    px[i] = inFile->beam[ibeam].bandData[iband].astro_data.freq[i];
	    py1[i] = 0.0;
	    py2[i] = 0.0;
	    for (j=0;j<ndump;j++)
	      {
		py1[i] += inFile->beam[ibeam].bandData[iband].astro_data.pol1[i+j*nchan];
		py2[i] += inFile->beam[ibeam].bandData[iband].astro_data.pol2[i+j*nchan];
	      }

	    if (recalc!=2)
	      {
		if (px[i] > maxx) maxx = px[i];
		if (px[i] < minx) minx = px[i];
	      }


	    if (inFile->beam[ibeam].bandData[iband].astro_data.flag[i] == 0)
	      {
 		if (py1[i] > maxy) {maxy = py1[i]; maxPosX = px[i];}
		if (py1[i] < miny) miny = py1[i];
		if (py2[i] > maxy) {maxy = py2[i]; maxPosX = px[i];}
		if (py2[i] < miny) miny = py2[i];
	      }
	      
	    
	    if (regionType==1 && inFile->beam[ibeam].bandData[iband].astro_data.flag[i] == 1) // No region set
	      {
		regionType=2;
		pf1[nFlagRegion] = px[i];
	      }	     
	    else if (regionType==2 && inFile->beam[ibeam].bandData[iband].astro_data.flag[i] == 0) 
	      {
		regionType=1;
		pf2[nFlagRegion++] = px[i];
	      }
	      
	  }
	if (regionType==2)
	    pf2[nFlagRegion++] = px[i-1];

	n = nchan;
	recalc=0;
      }
    
    cpgenv(minx,maxx,miny,maxy,0,1);
    sprintf(title,"Band: %s",inFile->beam[ibeam].bandHeader[iband].label);
    cpglab("Frequency (MHz)","Signal strength (arbitrary)",title);
    cpgline(n,px,py1);
    cpgsci(2); cpgline(n,px,py2); cpgsci(1);

    cpgsfs(3);
    cpgsci(5);
    for (i=0;i<nFlagRegion;i++)
      {
	cpgrect(pf1[i],pf2[i],miny,maxy);
      }
    cpgsci(1);
    cpgsfs(1);

    pointX[0] = maxPosX;
    pointY[0] = maxy;
    cpgsci(3);
    cpgpt(1,pointX,pointY,31);
    cpgsci(1);
    
    cpgcurs(&mx,&my,&key);
    if (key=='u')
      recalc=1;
    else if (key=='b') // Zap 5% of the band edges -- assuming Parkes UWL data
      {
	int ii,jj,zz;
	for (ii=0;ii<inFile->beam[ibeam].nBand;ii++)
	  {
	    for (jj=0;jj<inFile->beam[ibeam].bandHeader[ii].nchan;jj++)
	      {
		for (zz=0;zz<=26;zz++)
		  {
		    if (inFile->beam[ibeam].bandData[ii].astro_data.freq[jj] >= 704-6.4+zz*128.0 &&
			inFile->beam[ibeam].bandData[ii].astro_data.freq[jj] < 704+6.4+zz*128.0)
		      {			
			inFile->beam[ibeam].bandData[ii].astro_data.flag[jj] = 1;
		      }
		  }
	      }	
	  }
	recalc=1;
      }
    else if (key=='s')
      {
	saveFile(inFile,ibeam);
      }
    else if (key=='+')
      {
	iband++;
	if (iband >= inFile->beam[ibeam].nBand)
	  iband = inFile->beam[ibeam].nBand-1;
	recalc=1;
      }
    else if (key=='-')
      {
	iband--;
	recalc=1;
	if (iband < 0) iband=0;
      }
    else if (key=='z')
      {
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
    else if (key=='Z')
      {
	float lowX,highX;
	cpgband(4,0,mx,my,&mx2,&my2,&key);
	if (mx != mx2)
	  {
	    if (mx < mx2)
	      { lowX = mx; highX = mx2;}
	    else
	      { lowX = mx2; highX = mx;}
	    for (i=0;i<inFile->beam[ibeam].nBand;i++)
	      {
		for (j=0;j<inFile->beam[ibeam].bandHeader[i].nchan;j++)
		  { 
		    
		    if (inFile->beam[ibeam].bandData[i].astro_data.freq[j] >= lowX && inFile->beam[ibeam].bandData[i].astro_data.freq[j] <= highX)
		      inFile->beam[ibeam].bandData[i].astro_data.flag[j] = 1;
		  }
	      }
	    recalc=2;
	  }
	else
	  printf("Please press 'z' and then move somewhere and click left mouse button\n");
	
      }
  } while (key != 'q');

  free(px);
  free(py1);
  free(py2);
  free(pf1);
  free(pf2);
}

void saveFile(sdhdf_fileStruct *inFile,int ibeam)
{
  char oname[MAX_STRLEN];
  char flagName[MAX_STRLEN];
  sdhdf_fileStruct *outFile;
  int i,j;
  hsize_t dims[1];
  hid_t dset_id,dataspace_id;
  herr_t status;
  int *outFlags;

  // Should check if already .flag extension
  //
  sprintf(oname,"%s.flag",inFile->fname);
  printf("Saving to %s, please wait ... \n",oname);

  
  if (!(outFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >outFile<\n");
      exit(1);
    }

  sdhdf_initialiseFile(outFile);
  sdhdf_openFile(oname,outFile,3);

  sdhdf_copyRemainder(inFile,outFile,0);
  // Now add the flag table
  for (i=0;i<inFile->beam[ibeam].nBand;i++)
    {
      sdhdf_writeFlags(outFile,ibeam,i,inFile->beam[ibeam].bandData[i].astro_data.flag,inFile->beam[ibeam].bandHeader[i].nchan,inFile->beam[ibeam].bandHeader[i].label);
    }


  sdhdf_closeFile(outFile);

  free(outFile);
  printf("Save completed\n");
}
