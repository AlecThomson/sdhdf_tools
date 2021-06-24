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

#define MAX_RFI 128

void saveFile(sdhdf_fileStruct *inFile);

int main(int argc,char *argv[])
{
  int i,j,ii,b,k,s;
  char fname[MAX_STRLEN];
  char fromFileName[MAX_STRLEN];
  sdhdf_fileStruct *inFile,*fromFile;
  int iband=0,nchan,idump,nband;
  int ndump=0,totSize,chanPos;
  //  spectralDumpStruct spectrum;
  int npol=4;
  int setnband=-1;
  int flagType=-1;
  sdhdf_rfi rfi[MAX_RFI];
  int nRFI;
  
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

  // Two types of flagging.  Type 1 is copying a flag table
  // from one file to another.  Type 2 is running an algorithm to choose what to flag  
  flagType=2;
  for (i=1;i<argc;i++)
    {      
      if (strcmp(argv[i],"-from")==0)	{strcpy(fromFileName,argv[++i]); flagType=1;}
    }

  if (flagType==1)
    {
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
	  if (flagType==2)
	    {
	      printf("Observatory site = %s\n",inFile->primary[0].telescope);
	      if (strcmp(inFile->primary[0].telescope,"Parkes")!=0)
		printf("WARNING: ONLY IMPLEMENTED PARKES OBSERVATORY\n"); // FIX ME	      
	      sdhdf_loadPersistentRFI(rfi,&nRFI,MAX_RFI,"parkes");
	      printf("Have loaded %d RFI signals\n",nRFI);
	    }
	  for (b=0;b<inFile->nBeam;b++)
	    {
	      nband = inFile->beam[b].nBand;
	      for (i=0;i<nband;i++)
		{
		  // A waste as don't need to load in the data: FIX ME
		  // NOT APPLYING TO ALL SPECTRAL DUMPS **** FIX ME
		  sdhdf_loadBandData(inFile,b,i,1);
		  if (flagType==1)
		    {
 		      for (j=0;j<inFile->beam[b].bandHeader[i].nchan;j++)
			inFile->beam[b].bandData[i].astro_data.dataWeights[j] =
			  fromFile->beam[b].bandData[i].astro_data.dataWeights[j];
		    }
		  else if (flagType==2)
		    {
		      double freq;
		      int flagIt;

		      for (s=0;s<inFile->beam[b].bandHeader[i].ndump;s++)
			{
			  printf("Processing spectral dump %d\n",s);
			  for (j=0;j<inFile->beam[b].bandHeader[i].nchan;j++)
			    {
			      freq = inFile->beam[b].bandData[i].astro_data.freq[j];
			      for (k=0;k<nRFI;k++) // CAN DO THIS MUCH QUICKER			   
				{
				  //			      printf("Checking %.5f %g %g\n",freq,rfi[k].f0,rfi[k].f1);
				  if (freq > rfi[k].f0 && freq < rfi[k].f1)
				    {
				      inFile->beam[b].bandData[i].astro_data.dataWeights[j+s*inFile->beam[b].bandHeader[i].nchan] = 0;
				      break;
				    }
				}
			    }
			}
		    }
		}
	    }
	  saveFile(inFile);
	  sdhdf_closeFile(inFile);
	}
    }
  free(inFile);

  if (flagType==1)
    {
      sdhdf_closeFile(fromFile);     
      free(fromFile);
    }
  
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
	  sdhdf_writeDataWeights(outFile,b,i,inFile->beam[b].bandData[i].astro_data.dataWeights,inFile->beam[b].bandHeader[i].nchan,inFile->beam[b].bandHeader[i].ndump,inFile->beam[b].bandHeader[i].label);
	}
    }

  sdhdf_closeFile(outFile);

  free(outFile);
}

