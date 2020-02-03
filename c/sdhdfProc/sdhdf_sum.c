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
// Software to sum spectra together.
//
// Usage:
// sdhdf_sum -o <outfilename> ... list of input filenames ...
//
// Compilation
// gcc -lm -o sdhdf_sum sdhdf_sum.c sdhdfProc.c -I../hdf5-1.10.4/src/ ../hdf5-1.10.4/src/.libs/libhdf5.a -ldl -lz 
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sdhdfProc.h"
#include "hdf5.h"

#define VNUM "v0.1"

void help()
{
  printf("sdhdf_sum:        %s\n",VNUM);
  printf("sdhfProc version: %s\n",SOFTWARE_VER);
  printf("author:           George Hobbs\n");
  printf("\n");
  printf("Software to sum spectra together\n");
  printf("\n\nCommand line arguments:\n\n");
  printf("-h                This help\n");
  
  printf("\nExample:\n\n");
  printf("TO DO\n");
  printf("---------------------\n");
}

int main(int argc,char *argv[])
{
  int i,j,k,p,b;
  int timeThrough=1;
  sdhdf_fileStruct *inFile;
  sdhdf_fileStruct *file0;
  sdhdf_fileStruct *outFile;
  //  spectralDumpStruct spectrum;
  float *out_data;
  float in_pol1;
  float in_pol2;
  float in_pol3;
  float in_pol4;
  double tav=0;
  sdhdf_bandHeaderStruct *inBandParams;  
  long nchan,npol,ndump;
  
  char outFileName[MAX_STRLEN];
  int nFiles=0;
  char **fname;
  FILE *fout;

  //  help();
  
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >onFile<\n");
      exit(1);
    }
  if (!(file0 = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >file0<\n");
      exit(1);
    }
  if (!(outFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >outFile<\n");
      exit(1);
    }

  sdhdf_initialiseFile(inFile);
  sdhdf_initialiseFile(file0);
  sdhdf_initialiseFile(outFile);

  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-o")==0)
	strcpy(outFileName,argv[++i]);
      else if (strcmp(argv[i],"-h")==0) // Should actually free memory
	exit(1);
      else
	nFiles++;
    }

  fname = (char **)malloc(sizeof(char *)*nFiles);
  for (i=0;i<nFiles;i++)
    fname[i] = (char *)malloc(sizeof(char)*MAX_STRLEN);

  nFiles=0;
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-o")==0)
	i++;
      else if (strcmp(argv[i],"-h")==0)
	{} // Do nothing
      else
	strcpy(fname[nFiles++],argv[i]);
    }
  printf("Have %d input files\n",nFiles);
  // Get information from the first file

  sdhdf_openFile(fname[0],file0,1);
  sdhdf_loadMetaData(file0);
  printf("Number of beams = %d\n",file0->nBeam);
  // Open the output file
  sdhdf_openFile(outFileName,outFile,3);
  printf("File opened for output >%s<\n",outFileName);
  
  for (b=0;b<file0->nBeam;b++)
    {
      inBandParams = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*file0->beam[b].nBand);      
      sdhdf_copyBandHeaderStruct(file0->beam[b].bandHeader,inBandParams,file0->beam[b].nBand);
                
      for (i=0;i<file0->beam[b].nBand;i++)
	{
	  tav=0;
	  sdhdf_loadBandData(file0,b,i,1);
	  
	  printf("Processing band: %s\n",file0->beam[b].bandHeader[i].label);
	  nchan = file0->beam[b].bandHeader[i].nchan;
	  npol  = file0->beam[b].bandHeader[i].npol;
	  
	  // GH: FIX ME
	  //      sdhdf_loadFrequencyAttributes(file0,i);
	  //      sdhdf_loadDataAttributes(file0,i);
	  //      strcpy(outFile->frequency_attr.frame,file0->frequency_attr.frame);
	  //      strcpy(outFile->frequency_attr.unit,file0->frequency_attr.unit);
	  //      strcpy(outFile->data_attr.unit,file0->data_attr.unit);
	  	  
	  ndump = 1;
	  out_data = (float *)calloc(sizeof(float),nchan*npol*ndump);      
	  
	  // OPEN AND CLOSE ALL THE FILES LOTS OF TIMES -- CAN IMPROVE THE CODING HERE
	  for (j=0;j<nFiles;j++)
	    {
	      printf(" ... %s\n",fname[j]);
	      sdhdf_initialiseFile(inFile);	  
	      sdhdf_openFile(fname[j],inFile,1);
	      sdhdf_loadMetaData(inFile);
	      
	      nchan = inFile->beam[b].bandHeader[i].nchan;
	      npol  = inFile->beam[b].bandHeader[i].npol;
	      ndump = 1;
	      
	      sdhdf_loadBandData(inFile,b,i,1);
	      tav += inFile->beam[b].bandHeader[i].dtime;
	      
	      for (k=0;k<nchan;k++)
		{
		  // 0 here should be the spectral dump number
		  if (npol == 1)
		    {
		      in_pol1 = inFile->beam[b].bandData[i].astro_data.pol1[k];
		      out_data[k+0*npol*nchan]         += (in_pol1)/nFiles;
		    }
		  else if (npol==2)
		    {
		      in_pol1 = inFile->beam[b].bandData[i].astro_data.pol1[k];
		      in_pol2 = inFile->beam[b].bandData[i].astro_data.pol2[k];
		      out_data[k+0*npol*nchan]         += (in_pol1)/nFiles;
		      out_data[k+0*npol*nchan+nchan]   += (in_pol2)/nFiles;
		      
		    }
		  else
		    {
		      in_pol1 = inFile->beam[b].bandData[i].astro_data.pol1[k];
		      in_pol2 = inFile->beam[b].bandData[i].astro_data.pol2[k];
		      in_pol3 = inFile->beam[b].bandData[i].astro_data.pol3[k];
		      in_pol4 = inFile->beam[b].bandData[i].astro_data.pol4[k];

		      out_data[k+0*npol*nchan]         += (in_pol1)/nFiles;
		      out_data[k+0*npol*nchan+nchan]   += (in_pol2)/nFiles;
		      out_data[k+0*npol*nchan+2*nchan] += (in_pol3)/nFiles;
		      out_data[k+0*npol*nchan+3*nchan] += (in_pol4)/nFiles;		      
		    }
		}
	      sdhdf_releaseBandData(inFile,b,i,1);
	      sdhdf_closeFile(inFile);
	    }

	  inBandParams[i].ndump=1;
	  inBandParams[i].dtime = tav;
	  sdhdf_writeSpectrumData(outFile,file0->beam[b].bandHeader[i].label,b,i,out_data,file0->beam[b].bandData[i].astro_data.freq,nchan,npol,ndump,0);
	  // GH: FIX ME
	  //      sdhdf_writeFrequencyAttributes(outFile,file0->bandHeader[i].label);
	  //      sdhdf_writeDataAttributes(outFile,file0->bandHeader[i].label);

	  free(out_data);
	}
      sdhdf_writeBandHeader(outFile,inBandParams,b,file0->beam[b].nBand,1);
      free(inBandParams);
    }
  sdhdf_writeHistory(outFile,file0->history,file0->nHistory);  
  sdhdf_copyRemainder(file0,outFile,0);
  
  sdhdf_closeFile(outFile);
  sdhdf_closeFile(file0);

  free(inFile);
  free(outFile);
  free(file0);

  for (i=0;i<nFiles;i++)
    free(fname[i]);
  free(fname);

  printf("Complete: have written file %s\n",outFileName);
}

