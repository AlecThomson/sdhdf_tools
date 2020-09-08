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
// Software to produce a new SDHDF file based on parts of an existing one
//
// Usage:
// sdhdf_extract -f <filename.hdf> -o <outputFile.hdf>
//
// Compilation
// gcc -lm -o sdhdf_extract sdhdf_extract.c sdhdfProc.c -I../hdf5-1.10.4/src/ ../hdf5-1.10.4/src/.libs/libhdf5.a -ldl -lz -L/pulsar/psr/software/stable/stretch/lib/ -I//pulsar/psr/software/stable/stretch/include -lcalceph -Isofa/20190722/c/src/ -Lsofa/20190722/c/src -lsofa_c
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sdhdfProc.h"
#include "hdf5.h"

#define MAX_BANDS 26     // FIX THIS

//
// Note that the spectral dump number may be different for different subbands
// Should use the band extract first
//

int main(int argc,char *argv[])
{
  int ii,i,j,k,kk,jj,l,nchan,totNchan,b;
  char fname[MAX_FILES][64];
  int nFiles=0;
  char ext[MAX_STRLEN]="extract";
  char oname[MAX_STRLEN];
  sdhdf_fileStruct *inFile,*outFile;
  herr_t status;
  int selectDump[MAX_BANDS];
  float *outVals,*freqVals,*inData;
  int  nSelectDumps=0;
  int  copySD=0;
  float fSelect0[MAX_BANDS];
  float fSelect1[MAX_BANDS];
  int zoomBand=0;
  sdhdf_obsParamsStruct *outObsParams;
  sdhdf_bandHeaderStruct *outBandParams,*outCalBandParams,*inBandParams,*inCalBandParams;
  int nBand=0;
  char zoomLabel[MAX_STRLEN];
  char groupName[MAX_STRLEN];

  long npol,ndump,out_ndump,out_ndump_num;
  float *out_data,*out_freq;
  
  strcpy(oname,"sdhdf_extract_output.hdf");
  
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }

  if (!(outFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >outFile<\n");
      exit(1);
    }
  
  
  for (i=1;i<argc;i++)
    {      
      if (strcmp(argv[i],"-e")==0)
	strcpy(ext,argv[++i]);
      else if (strcmp(argv[i],"-d")==0)
	sscanf(argv[++i],"%d",&selectDump[nSelectDumps++]);
      else
	{
	  strcpy(fname[nFiles],argv[i]);
	  nFiles++;
	}
    }

  for (ii=0;ii<nFiles;ii++)
    {
      sdhdf_initialiseFile(inFile);
      sdhdf_initialiseFile(outFile);
      
      sprintf(oname,"%s.%s",fname[ii],ext);
            
      sdhdf_openFile(fname[ii],inFile,1);
      sdhdf_openFile(oname,outFile,3);
      
      if (inFile->fileID!=-1) // Did we successfully open the file?
	{
	  sdhdf_loadMetaData(inFile);
	  printf("%-22.22s %-22.22s %-5.5s %-6.6s %-20.20s %-10.10s %-10.10s %-7.7s %3d\n",fname,inFile->primary[0].utc0,
		 inFile->primary[0].hdr_defn_version,inFile->primary[0].pid,inFile->beamHeader[0].source,inFile->primary[0].telescope,
		 inFile->primary[0].observer, inFile->primary[0].rcvr,inFile->beam[0].nBand);

	  sdhdf_allocateBeamMemory(outFile,inFile->nBeam);
	  printf("Allocated memory\n");
	  for (b=0;b<inFile->nBeam;b++)
	    {
	      printf("Processing beam %d\n",b);
	      inBandParams = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*inFile->beam[b].nBand);      
	      sdhdf_copyBandHeaderStruct(inFile->beam[b].bandHeader,inBandParams,inFile->beam[b].nBand);
	      
	      for (i=0;i<inFile->beam[b].nBand;i++)
		{
		  nchan = inFile->beam[b].bandHeader[i].nchan;
		  npol  = inFile->beam[b].bandHeader[i].npol;
		  ndump = inFile->beam[b].bandHeader[i].ndump;
		  printf("Processing subband %d (number of spectral dumps = %d)\n",i,ndump);
		  out_ndump = nSelectDumps;
		  outObsParams = (sdhdf_obsParamsStruct *)malloc(sizeof(sdhdf_obsParamsStruct)*out_ndump);      
		  for (k=0;k<out_ndump;k++)
		    sdhdf_copySingleObsParams(inFile,b,i,k,&outObsParams[k]);


		  out_data = (float *)malloc(sizeof(float)*nchan*npol*out_ndump);
		  out_freq = (float *)malloc(sizeof(float)*nchan);

		  sdhdf_loadBandData(inFile,b,i,1);
		  for (j=0;j<nchan;j++)
		    out_freq[j] = inFile->beam[b].bandData[i].astro_data.freq[j];
		  out_ndump_num=0;
		  for (j=0;j<ndump;j++)
		    {
		      copySD=0;
		      for (jj=0;jj<nSelectDumps;jj++)
			{
			  if (selectDump[jj]==j)
			    {
			      printf("Copying subband %d\n",j);
			      for (kk=0;kk<nchan;kk++)
				{
				  if (npol==1)
				    out_data[kk+out_ndump_num*nchan]         = inFile->beam[b].bandData[i].astro_data.pol1[kk+j*nchan];
				  else if (npol==2)
				    {
				      out_data[kk+out_ndump_num*nchan*2]         = inFile->beam[b].bandData[i].astro_data.pol1[kk+j*nchan];
				      out_data[kk+nchan+out_ndump_num*nchan*2]   = inFile->beam[b].bandData[i].astro_data.pol2[kk+j*nchan];
				    }
				  else if (npol==4)
				    {
				      out_data[kk        +out_ndump_num*nchan*4] = inFile->beam[b].bandData[i].astro_data.pol1[kk+j*nchan];
				      out_data[kk+nchan  +out_ndump_num*nchan*4] = inFile->beam[b].bandData[i].astro_data.pol2[kk+j*nchan];
				      out_data[kk+2*nchan+out_ndump_num*nchan*4] = inFile->beam[b].bandData[i].astro_data.pol3[kk+j*nchan];
				      out_data[kk+3*nchan+out_ndump_num*nchan*4] = inFile->beam[b].bandData[i].astro_data.pol4[kk+j*nchan];
				    }
				}
			      out_ndump_num++;
			    }
			}
		    }
		  //		  sdhdf_writeBandHeader(outFile,outBandParams,b,nSelectBands,1);
		  //		  sdhdf_writeBandHeader(outFile,outCalBandParams,b,nSelectBands,2);
		  printf("Writing spectral data\n");
		  sdhdf_writeSpectrumData(outFile,inFile->beam[b].bandHeader[i].label,b,i,out_data,out_freq,nchan,npol,out_ndump,0);
		  printf("Writing obs params\n");
		  sdhdf_writeObsParams(outFile,inFile->beam[b].bandHeader[i].label,b,i,outObsParams,out_ndump);
		  printf("Releasing data\n");
		  sdhdf_releaseBandData(inFile,b,ii,1);		  

		  free(outObsParams);
		  inBandParams[i].ndump = out_ndump;
		}
	      printf("Writing band header\n");
	      sdhdf_writeBandHeader(outFile,inBandParams,b,inFile->beam[b].nBand,1);
	      free(inBandParams);
	    }
	  printf("Freeing\n");
	    

	      

	  printf("Writing other tables\n");
	  // Copy other primary tables
	  sdhdf_writeHistory(outFile,inFile->history,inFile->nHistory);
	  printf("Copy remainder\n");
	  sdhdf_copyRemainder(inFile,outFile,0);
	  printf("Closing files\n");
	  sdhdf_closeFile(inFile);
	  sdhdf_closeFile(outFile);
	}
    }

  free(inFile);
  free(outFile);
}


