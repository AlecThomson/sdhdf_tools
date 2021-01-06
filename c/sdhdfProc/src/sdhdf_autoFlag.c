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
// Software to flag data automatically
//
// Usage:
// sdhdf_autoFlag
//
// Compilation
// gcc -lm -o sdhdf_autoFlag sdhdf_autoFlag.c sdhdfProc.c -I../hdf5-1.10.4/src/ ../hdf5-1.10.4/src/.libs/libhdf5.a -ldl -lz -lcpgplot -L/pulsar/psr/software/stable/stretch/lib/ -I//pulsar/psr/software/stable/stretch/include -lcalceph -Isofa/20190722/c/src/ -Lsofa/20190722/c/src -lsofa_c 
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sdhdfProc.h"
#include "hdf5.h"
#include <cpgplot.h>
#include "TKfit.h"

void saveFile(sdhdf_fileStruct *inFile);

int main(int argc,char *argv[])
{
  int i,j,ii,b;
  char fname[MAX_STRLEN];
  char fromFileName[MAX_STRLEN];
  sdhdf_fileStruct *inFile,*fromFile;
  int iband=0,nchan,idump,nband;
  int ndump=0,totSize,chanPos;
  //  spectralDumpStruct spectrum;
  int npol=4;
  int setnband=-1;
  
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }
  if (!(fromFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >fromFile<\n");
      exit(1);
    }

    
  for (i=1;i<argc;i++)
    {      
      if (strcmp(argv[i],"-from")==0)	strcpy(fromFileName,argv[++i]);	
    }

  sdhdf_initialiseFile(fromFile);
  sdhdf_openFile(fromFileName,fromFile,1);
  sdhdf_loadMetaData(fromFile);
  printf("Loading the bands\n");
  for (b=0;b<fromFile->nBeam;b++)
    {
      // This is a waste as don't need to load in all the data: FIX ME
      for (i=0;i<fromFile->beam[b].nBand;i++)
	sdhdf_loadBandData(fromFile,b,i,1);
    }
  //

  for (ii=1;ii<argc;ii++)
    {
      if (strcmp(argv[ii],"-from")==0)
	ii++;
      else
	{
	  strcpy(fname,argv[ii]);
	  sdhdf_initialiseFile(inFile);
	  sdhdf_openFile(fname,inFile,1);
	  sdhdf_loadMetaData(inFile);
	  for (b=0;b<inFile->nBeam;b++)
	    {
	      nband = inFile->beam[b].nBand;
	      for (i=0;i<nband;i++)
		{
		  // Again a waste as don't need to load in the data: FIX ME
		  sdhdf_loadBandData(inFile,b,i,1);
		  for (j=0;j<inFile->beam[b].bandHeader[i].nchan;j++)
		    inFile->beam[b].bandData[i].astro_data.dataWeights[j] =
		      fromFile->beam[b].bandData[i].astro_data.dataWeights[j];
		}
	    }
	  saveFile(inFile);
	  sdhdf_closeFile(inFile);
	}
    }
  sdhdf_closeFile(fromFile);
  free(inFile);
  free(fromFile);
  
}


// Should make a generic saveFile function and also update sdhdf_flag.c

void saveFile(sdhdf_fileStruct *inFile)
{
  char oname[MAX_STRLEN];
  char flagName[MAX_STRLEN];
  sdhdf_fileStruct *outFile;
  int i,j,b;
  hsize_t dims[1];
  hid_t dset_id,dataspace_id;
  herr_t status;
  int *outFlags;

  // Should check if already .flag extension
  //
  sprintf(oname,"%s.autoflag",inFile->fname);

  
  if (!(outFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >outFile<\n");
      exit(1);
    }

  sdhdf_initialiseFile(outFile);
  sdhdf_openFile(oname,outFile,3);
  sdhdf_copyRemainder(inFile,outFile,0);
  
  // Now add the flag table
  for (b=0;b<inFile->nBeam;b++)
    {
      for (i=0;i<inFile->beam[b].nBand;i++)
	{
	  sdhdf_writeFlags(outFile,b,i,inFile->beam[b].bandData[i].astro_data.flag,inFile->beam[b].bandHeader[i].nchan,inFile->beam[b].bandHeader[i].label);
	}
    }

  sdhdf_closeFile(outFile);

  free(outFile);
}

