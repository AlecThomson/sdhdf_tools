//  Copyright (C) 2019, 2020, 2021 George Hobbs

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sdhdfProc.h"

#define VERSION "v0.5"
#define MAX_POL_CAL_CHAN 4096    // FIX ME -- SHOULD SET DYNAMICALLY

void help()
{
  printf("sdhdf_calibrate %s (SDHDFProc %s)\n",VERSION,SOFTWARE_VER);
  printf("Authors: G. Hobbs\n");
  printf("Purpose: to calibrate the data sets\n");
  printf("\n");
  printf("Command line arguments:\n\n");
  exit(1);
}


int main(int argc,char *argv[])
{
  int i,j,k,ii,beam,band;
  int nFiles=0;
  char fname[MAX_FILES][MAX_STRLEN];
  char runtimeDir[MAX_STRLEN];
  float pol1,pol2,pol3,pol4;
  sdhdf_fileStruct *inFile;
  sdhdf_calibration *polCal;
  int nPolCalChan=0;
  FILE *fout;

  int npol,nchan,ndump,nchanCal,ndumpCal;



  
  // Load in the file names
  for (i=1;i<argc;i++)
    {
      strcpy(fname[nFiles++],argv[i]);
    }
  
  inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct));
  polCal = (sdhdf_calibration *)malloc(sizeof(sdhdf_calibration)*MAX_POL_CAL_CHAN);



  if (getenv("SDHDF_RUNTIME")==0)
    {
      printf("=======================================================================\n");
      printf("Error: sdhdf_calibrate requires that the SDHDF_RUNTIME directory is set\n");
      printf("=======================================================================\n");
      exit(1);
    }
  strcpy(runtimeDir,getenv("SDHDF_RUNTIME"));
  
  for (i=0;i<nFiles;i++)
    {
      sdhdf_initialiseFile(inFile);
      if (sdhdf_openFile(fname[i],inFile,1)==-1)
	printf("Warning: unable to open file >%s<. Skipping\n",fname[i]);
      else
	{
	  sdhdf_loadMetaData(inFile);
	  sdhdf_loadPCM(polCal,&nPolCalChan,"parkes","UWL","uwl_181105_105441_b4.pcm"); // REMOVE HARDCODE
	  for (j=0;j<nPolCalChan;j++)
	    printf("%d %g %g %g %g %g %g %g %g %g %g %g\n",
		   j,polCal[j].freq,polCal[j].noiseSource_QoverI,polCal[j].noiseSource_UoverI,polCal[j].noiseSource_VoverI,
		   polCal[j].constant_gain,polCal[j].constant_diff_gain,polCal[j].constant_diff_phase,
		   polCal[j].constant_b1,polCal[j].constant_b2,polCal[j].constant_r1,polCal[j].constant_r2);
	  // Load the relevant PCM file
	  // FIX ME: Should check if the user wants polarisation calibration	
	  
	  for (beam=0;beam<inFile->nBeam;beam++)
	    {
	      for (band=0;band<inFile->beam[beam].nBand;band++)
		{
		  printf("Loading beam %d and band %d\n",beam,band);


		  npol = inFile->beam[beam].bandHeader[band].npol;
		  nchan = inFile->beam[beam].bandHeader[band].nchan;
		  ndump = inFile->beam[beam].bandHeader[band].ndump;
		  sdhdf_loadBandData(inFile,beam,band,1);
		  sdhdf_loadBandData(inFile,beam,band,2);
		  sdhdf_loadBandData(inFile,beam,band,3);

		  printf("Have loaded data\n");

		  nchanCal = inFile->beam[beam].calBandHeader[band].nchan;
		  ndumpCal = inFile->beam[beam].calBandHeader[band].ndump;


		  sdhdf_releaseBandData(inFile,beam,band,1);
		  sdhdf_releaseBandData(inFile,beam,band,2);
		  sdhdf_releaseBandData(inFile,beam,band,3);
		  
		}
	    }
	  sdhdf_closeFile(inFile);
	}
    }

  free(inFile);
  free(polCal);
}



