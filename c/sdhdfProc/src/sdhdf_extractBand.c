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

int main(int argc,char *argv[])
{
  int ii,i,j,k,l,nchan,totNchan,b;
  char fname[MAX_FILES][64];
  int nFiles=0;
  char ext[MAX_STRLEN]="extract";
  char oname[MAX_STRLEN];
  sdhdf_fileStruct *inFile,*outFile;
  herr_t status;
  char selectBand[MAX_BANDS][MAX_STRLEN];
  float *outVals,*freqVals,*inData;
  int  nSelectBands=0;
  int  copyBand=0;
  float fSelect0[MAX_BANDS];
  float fSelect1[MAX_BANDS];
  int zoomBand=0;
  sdhdf_bandHeaderStruct *outBandParams,*outCalBandParams,*inBandParams,*inCalBandParams;
  int nBand=0;
  char zoomLabel[MAX_STRLEN];
  char groupName[MAX_STRLEN];
  int cal=0;
  
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
      else if (strcmp(argv[i],"-zoom")==0)
	{
	  sscanf(argv[++i],"%f",&fSelect0[zoomBand]);
	  sscanf(argv[++i],"%f",&fSelect1[zoomBand]);
	  zoomBand++;
	}
      else if (strcmp(argv[i],"-b")==0)
	strcpy(selectBand[nSelectBands++],argv[++i]);
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
      
      if (zoomBand>0)
	nSelectBands=zoomBand;
      
      sdhdf_openFile(fname[ii],inFile,1);
      sdhdf_openFile(oname,outFile,3);
      
      if (inFile->fileID!=-1) // Did we successfully open the file?
	{
	  sdhdf_loadMetaData(inFile);
	  if (strcmp(inFile->primary[0].cal_mode,"OFF")==0)
	    cal=0;
	  else
	    cal=1;
	  
	  printf("%-22.22s %-22.22s %-5.5s %-6.6s %-20.20s %-10.10s %-10.10s %-7.7s %3d\n",fname,inFile->primary[0].utc0,
		 inFile->primary[0].hdr_defn_version,inFile->primary[0].pid,inFile->beamHeader[0].source,inFile->primary[0].telescope,
		 inFile->primary[0].observer, inFile->primary[0].rcvr,inFile->beam[0].nBand);
	  
	  for (b=0;b<inFile->nBeam;b++)
	    {
	      inBandParams  = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*inFile->beam[b].nBand);
	      outBandParams = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*nSelectBands);
	      sdhdf_copyBandHeaderStruct(inFile->beam[b].bandHeader,inBandParams,inFile->beam[b].nBand);
	      if (cal==1)
		{
		  inCalBandParams  = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*inFile->beam[b].nBand);         
		  outCalBandParams = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*nSelectBands);
		  sdhdf_copyBandHeaderStruct(inFile->beam[b].calBandHeader,inCalBandParams,inFile->beam[b].nBand);
		}
	      nBand=0;

	      if (zoomBand==0)
		{
		  for (i=0;i<inFile->beam[b].nBand;i++)
		    {
		      copyBand=0;
		      for (j=0;j<nSelectBands;j++)
			{
			  if (strcmp(inFile->beam[b].bandHeader[i].label,selectBand[j])==0)
			    {
			      copyBand=1;
			      break;
			    }
			}
		      if (copyBand==1)
			{
			  sprintf(groupName,"beam_%d",b);

			  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
			    {
			      hid_t group_id;
			      group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
			      status = H5Gclose(group_id);
			    }
			  sprintf(groupName,"beam_%d/%s",b,inFile->beam[b].bandHeader[i].label);
			  sdhdf_copyEntireGroup(groupName,inFile,outFile);

			  strcpy(outBandParams[nBand].label,inBandParams[i].label);
			  outBandParams[nBand].fc = inBandParams[i].fc;
			  outBandParams[nBand].f0 = inBandParams[i].f0;
			  outBandParams[nBand].f1 = inBandParams[i].f1;
			  outBandParams[nBand].nchan = inBandParams[i].nchan;
			  outBandParams[nBand].dtime = inBandParams[i].dtime;
			  outBandParams[nBand].npol = inBandParams[i].npol;
			  strcpy(outBandParams[nBand].pol_type, inBandParams[i].pol_type);
			  outBandParams[nBand].ndump = inBandParams[i].ndump;

			  if (cal==1)
			    {
			      strcpy(outCalBandParams[nBand].label,inCalBandParams[i].label);
			      outCalBandParams[nBand].fc = inCalBandParams[i].fc;
			      outCalBandParams[nBand].f0 = inCalBandParams[i].f0;
			      outCalBandParams[nBand].f1 = inCalBandParams[i].f1;
			      outCalBandParams[nBand].nchan = inCalBandParams[i].nchan;
			      outCalBandParams[nBand].dtime = inCalBandParams[i].dtime;
			      outCalBandParams[nBand].npol = inCalBandParams[i].npol;
			      strcpy(outCalBandParams[nBand].pol_type, inCalBandParams[i].pol_type);
			      outCalBandParams[nBand].ndump = inCalBandParams[i].ndump;
			    }
			  nBand++;
			}
		    }
		}
	      else
		{
		  // FIX ME -- THIS ISN'T CORRECTLY DEALING WITH OBS_PARAMS NOR CALS
		  for (j=0;j<zoomBand;j++)
		    {
		      sprintf(zoomLabel,"band_zoom%03d",j);
		      strcpy(outBandParams[j].label,zoomLabel);
		      outBandParams[j].fc = (fSelect0[j] + fSelect1[j])/2.;
		      outBandParams[j].f0 = fSelect0[j];
		      outBandParams[j].f1 = fSelect1[j];
		      strcpy(outBandParams[j].pol_type,inBandParams[0].pol_type);
		      outBandParams[j].npol = inBandParams[0].npol;
		      outBandParams[j].ndump = 1; // FIX THIS

		      if (cal==1)
			{
			  strcpy(outCalBandParams[j].label,zoomLabel);
			  outCalBandParams[j].fc = (fSelect0[j] + fSelect1[j])/2.;
			  outCalBandParams[j].f0 = fSelect0[j];
			  outCalBandParams[j].f1 = fSelect1[j];
			  strcpy(outCalBandParams[j].pol_type,inCalBandParams[0].pol_type);
			  outCalBandParams[j].npol = inCalBandParams[0].npol;
			  outCalBandParams[j].ndump = 1; // FIX THIS
			}
		      
		      // Count nchan (note that we may go across band edges)
		      totNchan=0;
		      
		      for (k=0;k<inFile->beam[b].nBand;k++)
			{
			  sdhdf_loadBandData(inFile,b,k,1);
			  for (l=0;l<inFile->beam[b].bandHeader[k].nchan;l++)
			    {
			      if (inFile->beam[b].bandData[k].astro_data.freq[l] >= fSelect0[j] &&
				  inFile->beam[b].bandData[k].astro_data.freq[l] <= fSelect1[j])
				totNchan++;
			    }
			}
		      
		      outBandParams[j].nchan = totNchan;
		      outBandParams[j].dtime = 1;// FIX	   	      
		      outVals  = (float *)malloc(sizeof(float)*outBandParams[j].nchan*4*1); // Fix 4 = npol, 1 = ndump
		      freqVals = (float *)malloc(sizeof(float)*outBandParams[j].nchan); 
		      
		      nchan=0;
		      for (k=0;k<inFile->beam[b].nBand;k++)
			{
			  // inData = (float *)malloc(sizeof(float)*inFile->beam[b].bandHeader[k].nchan*inFile->beam[b].bandHeader[k].ndump*4); // 4 = npol FIX
			  //  sdhdf_loadEntireDump(inFile,k,inData);
			  for (l=0;l<inFile->beam[b].bandHeader[k].nchan;l++)
			    {
			      if (inFile->beam[b].bandData[k].astro_data.freq[l] >= fSelect0[j] && inFile->beam[b].bandData[k].astro_data.freq[l] <= fSelect1[j])
				{
				  outVals[nchan] = inFile->beam[b].bandData[k].astro_data.pol1[l]; // FIX FOR NDUMP
				  // FIX FOR NPOL
				  outVals[nchan+totNchan] = inFile->beam[b].bandData[k].astro_data.pol2[l]; // FIX FOR NDUMP
				  outVals[nchan+2*totNchan] = inFile->beam[b].bandData[k].astro_data.pol3[l]; // FIX FOR NDUMP
				  outVals[nchan+3*totNchan] = inFile->beam[b].bandData[k].astro_data.pol4[l]; // FIX FOR NDUMP
				  freqVals[nchan] = inFile->beam[b].bandData[k].astro_data.freq[l];
				  nchan++;
				}
			    }
			  //			  free(inData);
			}
		      printf("Output nchan = %d (%d)\n",nchan,totNchan);
		      //		      sdhdf_writeSpectrumData(outFile,inFile,b,j,outVals,freqVals,nchan,4,1,0); // FIX 4,1,0
		      sdhdf_writeSpectrumData(outFile,outBandParams[j].label,b,j,outVals,freqVals,nchan,4,1,0); // FIX 4,1,0
		      //		      sdhdf_writeNewBand(outFile,j,outVals,freqVals,&(outBandParams[j]));
		    
		      
		      free(outVals);
		      free(freqVals);
		    }
	      
		}
	      sdhdf_writeBandHeader(outFile,outBandParams,b,nSelectBands,1);
	      if (cal==1)
		sdhdf_writeBandHeader(outFile,outCalBandParams,b,nSelectBands,2);

	      
	      //	      sprintf(groupName,"beam_%d/metadata",b);
	      //	      sdhdf_copyEntireGroup(groupName,inFile,outFile);	      	      

	      inFile->beamHeader[b].nBand = nSelectBands;

	    }
       
	  // Copy other primary tables

	  sdhdf_copyEntireGroup("metadata",inFile,outFile);	      	      
	  sdhdf_writeBeamHeader(outFile,inFile->beamHeader,inFile->nBeam); 

	  // Don't want to "copyRemainder" as have only selected specific bands
	  sdhdf_copyEntireGroup("config",inFile,outFile);
	  	  
	  free(outBandParams);
	  free(inBandParams);
	  if (cal==1)
	    {
	      free(outCalBandParams);
	      free(inCalBandParams);
	    }
	  sdhdf_closeFile(inFile);
	  sdhdf_closeFile(outFile);
	}
    }

  free(inFile);
  free(outFile);
}



