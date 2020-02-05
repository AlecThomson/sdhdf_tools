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
// Software to output text files relating to the data sets in the file
//
// Usage:
// sdhdf_quickdump <filename.hdf>  <filename.hdf> ...
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sdhdfProc.h"

#define VERSION "v0.5"

void help()
{
  printf("sdhdf_quickdump %s (SDHDFProc %s)\n",VERSION,SOFTWARE_VER);
  printf("Authors: G. Hobbs\n");
  printf("Purpose: to present data within a file\n");
  printf("\n");
  printf("Command line arguments:\n\n");
  exit(1);

}

int main(int argc,char *argv[])
{
  int i,j,k,beam,band;
  int display=0;
  int nFiles=0;
  char fname[MAX_FILES][MAX_STRLEN];
  float freq;
  int   flag;
  float pol1,pol2,pol3,pol4;
  sdhdf_fileStruct *inFile;
  int npol=0;
  int outFile=0;
  float *tsys;
  char extFileName[MAX_STRLEN];
  char outFileName[MAX_STRLEN];
  FILE *fout;
  float freq0,freq1;
  int   setFreqRange=0;
  int sd0,sd1;
  int   setDumpRange=0;
  int nchan;
  int dataType=1;
  
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-h")==0)
	help();
      else if (strcmp(argv[i],"-tsys")==0)
	dataType=2;
      else if (strcmp(argv[i],"-freqRange")==0)
	{
	  sscanf(argv[++i],"%f",&freq0);
	  sscanf(argv[++i],"%f",&freq1);
	  setFreqRange=1;
	}
      else if (strcmp(argv[i],"-dumpRange")==0)
	{
	  sscanf(argv[++i],"%d",&sd0);
	  sscanf(argv[++i],"%d",&sd1);
	  setDumpRange=1;
	}
      else if (strcmp(argv[i],"-e")==0)
	{
	  strcpy(extFileName,argv[++i]);
	  outFile=1;
	}
      else
	strcpy(fname[nFiles++],argv[i]);
    }

  inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct));
  
  for (i=0;i<nFiles;i++)
    {
      sdhdf_initialiseFile(inFile);
      if (sdhdf_openFile(fname[i],inFile,1)==-1)
	printf("Warning: unable to open file >%s<. Skipping\n",fname[i]);
      else
	{
	  sdhdf_loadMetaData(inFile);
	  if (outFile==1)
	    {
	      sprintf(outFileName,"%s.%s",inFile->fname,extFileName);
	      if (!(fout = fopen(outFileName,"w")))
		{
		  printf("ERROR: Unable to open file >%s< for writing\n",outFileName);
		  exit(1);
		}
	    }
	  
	  for (beam=0;beam<inFile->nBeam;beam++)
	    {
	      for (band=0;band<inFile->beam[beam].nBand;band++)
		{
		  npol = inFile->beam[beam].bandHeader[band].npol;
		  nchan = inFile->beam[beam].bandHeader[band].nchan;
		  sdhdf_loadBandData(inFile,beam,band,1);
		  if (dataType==2)
		    {
		      sdhdf_loadBandData(inFile,beam,band,2);
		      sdhdf_loadBandData(inFile,beam,band,3);
		    }
		  if (dataType==1) // Astronomy data
		    {
		      if (setDumpRange==0)
			{
			  sd0 = 0;
			  sd1 = inFile->beam[beam].bandHeader[band].ndump;
			}
		      
		      for (k=sd0;k<sd1;k++)
			{
			  //		      printf("Here with %d\n",k);
			  for (j=0;j<nchan;j++)
			    {
			      display=1;
			      //			  printf("Loading freq\n");
			      freq = inFile->beam[beam].bandData[band].astro_data.freq[j];
			      //			  printf("Loading pol1 %g\n",freq);
			      pol1 = inFile->beam[beam].bandData[band].astro_data.pol1[j+k*nchan];
			      //			  printf("Done load\n");
			      if (setFreqRange==1 && (freq < freq0 || freq > freq1))
				display=0;
			      
			      if (npol > 1)
				{
				  pol2 = inFile->beam[beam].bandData[band].astro_data.pol2[j+k*nchan];
				  if (npol > 2)
				    {
				      pol3 = inFile->beam[beam].bandData[band].astro_data.pol3[j+k*nchan];
				      pol4 = inFile->beam[beam].bandData[band].astro_data.pol4[j+k*nchan];
				    }
				  else
				    pol3 = pol4 = 0;
				}
			      else
				pol2 = pol3 = pol4 = 0;
			      
			      flag = inFile->beam[beam].bandData[band].astro_data.flag[j];
			      
			      if (display==1)
				{
				  if (outFile==1)
				    fprintf(fout,"%s %d %d %d %d %.6f %g %g %g %g %d\n",inFile->fname,beam,band,k,j,freq,pol1,pol2,pol3,pol4,flag);
				  else
				    printf("%s %d %d %d %d %.6f %g %g %g %g %d\n",inFile->fname,beam,band,k,j,freq,pol1,pol2,pol3,pol4,flag);
				}
			    }
			}
		    }
		  else if (dataType==2)
		    {
		      tsys = (float *)malloc(sizeof(float)*inFile->beam[beam].calBandHeader[band].ndump*inFile->beam[beam].calBandHeader[band].nchan*2); // * 2 = AA and BB
		      sdhdf_loadCalProc(inFile,beam,band,"cal_proc_tsys",tsys);
		      for (k=0;k<inFile->beam[beam].calBandHeader[band].ndump;k++)
			{
			  //		      printf("Here with %d\n",k);
			  for (j=0;j<inFile->beam[beam].calBandHeader[band].nchan;j++)
			    {			      
			      if (outFile==1)
				fprintf(fout,"%s %d %d %d %d %.6f %g %g\n",inFile->fname,beam,band,k,j,inFile->beam[beam].bandData[band].cal_on_data.freq[j],
					tsys[k*2*inFile->beam[beam].calBandHeader[band].nchan + j],
					tsys[k*2*inFile->beam[beam].calBandHeader[band].nchan+inFile->beam[beam].calBandHeader[band].nchan + j]);
			      else
				printf("%s %d %d %d %d %.6f %g %g\n",inFile->fname,beam,band,k,j,inFile->beam[beam].bandData[band].cal_on_data.freq[j],
					tsys[k*2*inFile->beam[beam].calBandHeader[band].nchan + j],
					tsys[k*2*inFile->beam[beam].calBandHeader[band].nchan+inFile->beam[beam].calBandHeader[band].nchan + j]);

			    }
			  if (outFile==1)
			    fprintf(fout,"\n");
			}
		      free(tsys);
		    }
		  sdhdf_releaseBandData(inFile,beam,band,1);
		  
		}
	    }
	  if (outFile==1)
	    fclose(fout);
	  sdhdf_closeFile(inFile);
	}
    }

  free(inFile);
}
