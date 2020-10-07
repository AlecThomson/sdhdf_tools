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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "sdhdfProc.h"



//
// Initialise the sdhdf_fileStruct structure
// This should be run as soon as the memory for the file structure has been allocated
//
void sdhdf_initialiseFile(sdhdf_fileStruct *inFile)
{
  inFile->fileOpen = 0;
  strcpy(inFile->fname,"unset");
  inFile->nPrimary = -1;
  inFile->primaryAllocatedMemory = 0;
  inFile->nSoftware = -1;
  inFile->softwareAllocatedMemory = 0;
  inFile->nHistory = -1;
  inFile->historyAllocatedMemory = 0;
  inFile->nBeam = -1;
  inFile->beamAllocatedMemory = 0;  
}

// Open an SDHDF file
// openType =
//            1 Open for read only
//            2 Open for read/write
//            3 Open for read/write as truncate
// Returns -1 if unsuccessfully opened
//
int sdhdf_openFile(char *fname,sdhdf_fileStruct *inFile,int openType)
{
  if (openType==1)
    inFile->fileID  = H5Fopen(fname,H5F_ACC_RDONLY,H5P_DEFAULT);
  else if (openType==2)
    inFile->fileID  = H5Fopen(fname,H5F_ACC_RDWR,H5P_DEFAULT);
  else if (openType==3)
    inFile->fileID  = H5Fcreate(fname,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  if (inFile->fileID != -1)
    {
      strcpy(inFile->fname,fname);
      inFile->fileOpen = 1;
    }
  return inFile->fileID;
}

// Routine to close a file that is already open
// it also frees memory used by the sdhdf_fileStruct
//
void sdhdf_closeFile(sdhdf_fileStruct *inFile)
{
  herr_t status;
  int i,j;

  if (inFile->fileOpen==1)
    {
      if (inFile->beamAllocatedMemory > 0)
	{
	  for (i=0;i<inFile->nBeam;i++)
	    {
	      if (inFile->beam[i].bandAllocatedMemory > 0)
		{
		  for (j=0;j<inFile->beam[i].nBand;j++)
		    {
		      if (inFile->beam[i].bandData[j].astro_data.freqAllocatedMemory > 0)
			{
			  free(inFile->beam[i].bandData[j].astro_data.freq);
			  inFile->beam[i].bandData[j].astro_data.freqAllocatedMemory = 0;
			}
		      if (inFile->beam[i].bandData[j].astro_data.flagAllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].astro_data.flag); inFile->beam[i].bandData[j].astro_data.flagAllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].astro_data.pol1AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].astro_data.pol1); inFile->beam[i].bandData[j].astro_data.pol1AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].astro_data.pol2AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].astro_data.pol2); inFile->beam[i].bandData[j].astro_data.pol2AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].astro_data.pol3AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].astro_data.pol3); inFile->beam[i].bandData[j].astro_data.pol3AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].astro_data.pol4AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].astro_data.pol4); inFile->beam[i].bandData[j].astro_data.pol4AllocatedMemory = 0; }

		      if (inFile->beam[i].bandData[j].cal_on_data.freqAllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_on_data.freq); inFile->beam[i].bandData[j].cal_on_data.freqAllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_on_data.flagAllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_on_data.flag); inFile->beam[i].bandData[j].cal_on_data.flagAllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_on_data.pol1AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_on_data.pol1); inFile->beam[i].bandData[j].cal_on_data.pol1AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_on_data.pol2AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_on_data.pol2); inFile->beam[i].bandData[j].cal_on_data.pol2AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_on_data.pol3AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_on_data.pol3); inFile->beam[i].bandData[j].cal_on_data.pol3AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_on_data.pol4AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_on_data.pol4); inFile->beam[i].bandData[j].cal_on_data.pol4AllocatedMemory = 0; }

		      if (inFile->beam[i].bandData[j].cal_off_data.freqAllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_off_data.freq); inFile->beam[i].bandData[j].cal_off_data.freqAllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_off_data.flagAllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_off_data.flag); inFile->beam[i].bandData[j].cal_off_data.flagAllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_off_data.pol1AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_off_data.pol1); inFile->beam[i].bandData[j].cal_off_data.pol1AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_off_data.pol2AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_off_data.pol2); inFile->beam[i].bandData[j].cal_off_data.pol2AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_off_data.pol3AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_off_data.pol3); inFile->beam[i].bandData[j].cal_off_data.pol3AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_off_data.pol4AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_off_data.pol4); inFile->beam[i].bandData[j].cal_off_data.pol4AllocatedMemory = 0; }

		      if (inFile->beam[i].bandData[j].cal_proc_tsys.freqAllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_tsys.freq); inFile->beam[i].bandData[j].cal_proc_tsys.freqAllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_tsys.flagAllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_tsys.flag); inFile->beam[i].bandData[j].cal_proc_tsys.flagAllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_tsys.pol1AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_tsys.pol1); inFile->beam[i].bandData[j].cal_proc_tsys.pol1AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_tsys.pol2AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_tsys.pol2); inFile->beam[i].bandData[j].cal_proc_tsys.pol2AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_tsys.pol3AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_tsys.pol3); inFile->beam[i].bandData[j].cal_proc_tsys.pol3AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_tsys.pol4AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_tsys.pol4); inFile->beam[i].bandData[j].cal_proc_tsys.pol4AllocatedMemory = 0; }

		      if (inFile->beam[i].bandData[j].cal_proc_diff_gain.freqAllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_diff_gain.freq); inFile->beam[i].bandData[j].cal_proc_diff_gain.freqAllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_diff_gain.flagAllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_diff_gain.flag); inFile->beam[i].bandData[j].cal_proc_diff_gain.flagAllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_diff_gain.pol1AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_diff_gain.pol1); inFile->beam[i].bandData[j].cal_proc_diff_gain.pol1AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_diff_gain.pol2AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_diff_gain.pol2); inFile->beam[i].bandData[j].cal_proc_diff_gain.pol2AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_diff_gain.pol3AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_diff_gain.pol3); inFile->beam[i].bandData[j].cal_proc_diff_gain.pol3AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_diff_gain.pol4AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_diff_gain.pol4); inFile->beam[i].bandData[j].cal_proc_diff_gain.pol4AllocatedMemory = 0; }

		      if (inFile->beam[i].bandData[j].cal_proc_diff_phase.freqAllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_diff_phase.freq); inFile->beam[i].bandData[j].cal_proc_diff_phase.freqAllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_diff_phase.flagAllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_diff_phase.flag); inFile->beam[i].bandData[j].cal_proc_diff_phase.flagAllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_diff_phase.pol1AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_diff_phase.pol1); inFile->beam[i].bandData[j].cal_proc_diff_phase.pol1AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_diff_phase.pol2AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_diff_phase.pol2); inFile->beam[i].bandData[j].cal_proc_diff_phase.pol2AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_diff_phase.pol3AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_diff_phase.pol3); inFile->beam[i].bandData[j].cal_proc_diff_phase.pol3AllocatedMemory = 0; }
		      if (inFile->beam[i].bandData[j].cal_proc_diff_phase.pol4AllocatedMemory > 0)
			{ free(inFile->beam[i].bandData[j].cal_proc_diff_phase.pol4); inFile->beam[i].bandData[j].cal_proc_diff_phase.pol4AllocatedMemory = 0; }
				      

		      if (inFile->beam[i].bandData[j].astro_obsHeaderAllocatedMemory == 1)
			free(inFile->beam[i].bandData[j].astro_obsHeader);
		      if (inFile->beam[i].bandData[j].cal_obsHeaderAllocatedMemory == 1)
			free(inFile->beam[i].bandData[j].cal_obsHeader);


		    }
		  free(inFile->beam[i].bandData);
		  free(inFile->beam[i].bandHeader);

		  if (inFile->beam[i].calBandAllocatedMemory==1)
		    {
		      free(inFile->beam[i].calBandHeader);
		      inFile->beam[i].calBandAllocatedMemory = 0;
		    }
		  inFile->beam[i].bandAllocatedMemory = 0;
		}
	    }

	  free(inFile->beam);
	  inFile->beamAllocatedMemory = 0;
	}
      if (inFile->historyAllocatedMemory > 0)
	{free(inFile->history); inFile->historyAllocatedMemory = 0;}

      if (inFile->softwareAllocatedMemory > 0)
	{free(inFile->software); inFile->softwareAllocatedMemory = 0;}

      if (inFile->primaryAllocatedMemory > 0)
	{free(inFile->primary); inFile->primaryAllocatedMemory = 0;}

      status = H5Fclose(inFile->fileID);
      inFile->fileOpen = 0;
    }
}

// Type = 1: astro_data
// Type = 2: cal_on
// Type = 3: cal_off
void sdhdf_loadBandData(sdhdf_fileStruct *inFile,int beam,int band,int type)
{
  int i,j,n;
  int nchan;
  int ndump;
  int npol;
  char dataName[MAX_STRLEN];
  hid_t dataset_id;
  herr_t status;
  float *allData;

  // astro_data
  if (type==1)
    {
      nchan = inFile->beam[beam].bandHeader[band].nchan;
      ndump = inFile->beam[beam].bandHeader[band].ndump;
      npol  = inFile->beam[beam].bandHeader[band].npol;

      sdhdf_allocateBandData(&(inFile->beam[beam].bandData[band].astro_data),nchan,ndump,npol);
      
      sprintf(dataName,"beam_%d/%s/astronomy_data/frequency",beam,inFile->beam[beam].bandHeader[band].label);      
      dataset_id   = H5Dopen2(inFile->fileID,dataName,H5P_DEFAULT);
      status = H5Dread(dataset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,inFile->beam[beam].bandData[band].astro_data.freq);  

      sprintf(dataName,"beam_%d/%s/astronomy_data/data",beam,inFile->beam[beam].bandHeader[band].label);      
      dataset_id   = H5Dopen2(inFile->fileID,dataName,H5P_DEFAULT);
      allData = (float *)malloc(sizeof(float)*nchan*npol*ndump);
      status = H5Dread(dataset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,allData);  
      sdhdf_extractPols(&(inFile->beam[beam].bandData[band].astro_data),allData,nchan,ndump,npol);
      free(allData);
      status = H5Dclose(dataset_id);
      // Flags
      
      sprintf(dataName,"beam_%d/%s/astronomy_data/flag",beam,inFile->beam[beam].bandHeader[band].label);      
      if (sdhdf_checkGroupExists(inFile,dataName)==0)
	{
	  dataset_id   = H5Dopen2(inFile->fileID,dataName,H5P_DEFAULT);
	  status = H5Dread(dataset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,inFile->beam[beam].bandData[band].astro_data.flag);        
	  status = H5Dclose(dataset_id);
	}
    }
  // calibrator data
  else if (type==2 || type==3)
    {
      nchan = inFile->beam[beam].calBandHeader[band].nchan;
      ndump = inFile->beam[beam].calBandHeader[band].ndump;
      npol  = inFile->beam[beam].calBandHeader[band].npol;

      if (type==2)
	sdhdf_allocateBandData(&(inFile->beam[beam].bandData[band].cal_on_data),nchan,ndump,npol);
      else if (type==3)
	sdhdf_allocateBandData(&(inFile->beam[beam].bandData[band].cal_off_data),nchan,ndump,npol);
      
      sprintf(dataName,"beam_%d/%s/calibrator_data/cal_frequency",beam,inFile->beam[beam].bandHeader[band].label);      
      dataset_id   = H5Dopen2(inFile->fileID,dataName,H5P_DEFAULT);
      if (type==2)
	status = H5Dread(dataset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,inFile->beam[beam].bandData[band].cal_on_data.freq);  
      else if (type==3)
	status = H5Dread(dataset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,inFile->beam[beam].bandData[band].cal_off_data.freq);  

      if (type==2)
	sprintf(dataName,"beam_%d/%s/calibrator_data/cal_data_on",beam,inFile->beam[beam].bandHeader[band].label);      
      else if (type==3)
	sprintf(dataName,"beam_%d/%s/calibrator_data/cal_data_off",beam,inFile->beam[beam].bandHeader[band].label);
      
      dataset_id   = H5Dopen2(inFile->fileID,dataName,H5P_DEFAULT);
      allData = (float *)malloc(sizeof(float)*nchan*npol*ndump);
      status = H5Dread(dataset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,allData);  
      if (type==2)
	{
	  printf("Extrating with %d %d %d\n",nchan,ndump,npol);
	  sdhdf_extractPols(&(inFile->beam[beam].bandData[band].cal_on_data),allData,nchan,ndump,npol);
	}
	  else if (type==3)
	sdhdf_extractPols(&(inFile->beam[beam].bandData[band].cal_off_data),allData,nchan,ndump,npol);
      free(allData);

      status = H5Dclose(dataset_id);
    }
}



void sdhdf_extractPols(sdhdf_spectralDumpsStruct *spec,float *in,int nchan,int ndump,int npol)
{
  int i,j;

  // SHOULD USE MEMCPY
  for (i=0;i<ndump;i++)
    {
      for (j=0;j<nchan;j++)
	{
	  spec->pol1[j+nchan*i] = in[j+nchan*i*npol];
	  if (npol > 1)
	    {
	      spec->pol2[j+nchan*i] = in[j+nchan*i*npol+nchan];
	      if (npol > 2)
		{
		  spec->pol3[j+nchan*i] = in[j+nchan*i*npol+2*nchan];
		  spec->pol4[j+nchan*i] = in[j+nchan*i*npol+3*nchan];
		}
	    }
	}
    }

}

void sdhdf_allocateBandData(sdhdf_spectralDumpsStruct *spec,int nchan,int ndump,int npol)
{
  if (spec->freqAllocatedMemory == 0)
    {
      spec->freq = (float *)malloc(sizeof(float)*nchan);
      spec->freqAllocatedMemory = 1;
    }
  if (spec->flagAllocatedMemory == 0)
    {
      spec->flag = (int *)calloc(sizeof(int),nchan); // Initialise to 0
      spec->flagAllocatedMemory = 1;
    }
  if (spec->pol1AllocatedMemory == 0)
    {
      spec->pol1 = (float *)malloc(sizeof(float)*nchan*npol*ndump); 
      spec->pol1AllocatedMemory = 1;
    }
  if (npol > 1 && spec->pol2AllocatedMemory == 0)
    {
      spec->pol2 = (float *)malloc(sizeof(float)*nchan*npol*ndump); 
      spec->pol2AllocatedMemory = 1;
    }
  if (npol > 3 && spec->pol3AllocatedMemory == 0 && spec->pol4AllocatedMemory == 0)
    {
      spec->pol3 = (float *)malloc(sizeof(float)*nchan*ndump);
      spec->pol4 = (float *)malloc(sizeof(float)*nchan*ndump); 
      spec->pol3AllocatedMemory = 1;
      spec->pol4AllocatedMemory = 1;
    }
}
  

// Type = 1: astro_data
// FIX ME __ OTHER TYPES
void sdhdf_releaseBandData(sdhdf_fileStruct *inFile,int beam,int band,int type)
{
  if (type==1)
    {
      if (inFile->beam[beam].bandData[band].astro_data.freqAllocatedMemory == 1)
	{free(inFile->beam[beam].bandData[band].astro_data.freq); inFile->beam[beam].bandData[band].astro_data.freqAllocatedMemory = 0;}
      if (inFile->beam[beam].bandData[band].astro_data.flagAllocatedMemory == 1)
	{free(inFile->beam[beam].bandData[band].astro_data.flag); inFile->beam[beam].bandData[band].astro_data.flagAllocatedMemory = 0;}
      if (inFile->beam[beam].bandData[band].astro_data.pol1AllocatedMemory == 1)
	{free(inFile->beam[beam].bandData[band].astro_data.pol1); inFile->beam[beam].bandData[band].astro_data.pol1AllocatedMemory = 0;}
      if (inFile->beam[beam].bandData[band].astro_data.pol2AllocatedMemory == 1)
	{free(inFile->beam[beam].bandData[band].astro_data.pol2); inFile->beam[beam].bandData[band].astro_data.pol2AllocatedMemory = 0;}
      if (inFile->beam[beam].bandData[band].astro_data.pol3AllocatedMemory == 1)
	{free(inFile->beam[beam].bandData[band].astro_data.pol3); inFile->beam[beam].bandData[band].astro_data.pol3AllocatedMemory = 0;}
      if (inFile->beam[beam].bandData[band].astro_data.pol4AllocatedMemory == 1)
	{free(inFile->beam[beam].bandData[band].astro_data.pol4); inFile->beam[beam].bandData[band].astro_data.pol4AllocatedMemory = 0;}	      
    }
}

void sdhdf_add1arg(char *args,char *add)
{
  if (strlen(args) < MAX_STRLEN - strlen(add) - 1)
    {
      strcat(args,add);
      strcat(args," ");
    }
}


void sdhdf_add2arg(char *args,char *add1,char *add2)
{
  if (strlen(args) < MAX_STRLEN - strlen(add1) - 1)
    {
      strcat(args,add1);
      strcat(args," ");
    }
  if (strlen(args) < MAX_STRLEN - strlen(add2) - 1)
    {
      strcat(args,add2);
      strcat(args," ");
    }
}


// ** SHOULD COPY FLAG TABLE AS WELL **

// type = 0: update DATA and frequency, but copy the cal datasets
// type = 1: update data and frequency only
// type = 2: update cal_on and frequency only
// type = 3: update cal_off only
// type = 4: update flag, data and frequency only

void sdhdf_writeSpectrumData(sdhdf_fileStruct *outFile,char *blabel, int ibeam,int iband,  float *out,float *freq,long nchan,long npol,long nsub,int type)
{
  int i;
  hid_t dset_id,datatype_id,group_id,ocpl_id;
  herr_t status;
  hid_t dataspace_id,stid;
  char groupName[MAX_STRLEN];
  char dsetName[MAX_STRLEN],label[MAX_STRLEN];
  hsize_t dims[5];

  if (type==0 || type==1 || type==4)
    {
      sprintf(groupName,"beam_%d",ibeam);
      if (sdhdf_checkGroupExists(outFile,groupName) == 1)
	{
	  group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	  status = H5Gclose(group_id);
	}
      sprintf(groupName,"beam_%d/%s",ibeam,blabel);
      if (sdhdf_checkGroupExists(outFile,groupName) == 1)
	{
	  group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	  status = H5Gclose(group_id);
	}
      sprintf(groupName,"beam_%d/%s/astronomy_data",ibeam,blabel);
      if (sdhdf_checkGroupExists(outFile,groupName) == 1)
	{
	  group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	  status = H5Gclose(group_id);
	}
      // FIX ME -- SHOULD MAKE CAL DATA SETS
    }

  dims[0] = nsub;
  dims[1] = 1;
  dims[2] = npol;
  dims[3] = nchan;
  dims[4] = 1;
  dataspace_id = H5Screate_simple(5,dims,NULL);

  /* FIX ME -- SHOULD COPY CAL INFO
  if (type==0)    
    {
      ocpl_id = H5Pcreate(H5P_OBJECT_COPY);
      sprintf(label,"%s/cal_data_on",groupName);
      status = H5Ocopy(inFile->fileID, label, outFile->fileID, label, ocpl_id, H5P_DEFAULT);
      sprintf(label,"%s/cal_data_off",groupName);
      status = H5Ocopy(inFile->fileID, label, outFile->fileID, label, ocpl_id, H5P_DEFAULT);
      sprintf(label,"%s/cal_frequency",groupName);
      status = H5Ocopy(inFile->fileID, label, outFile->fileID, label, ocpl_id, H5P_DEFAULT);

      H5Pclose(ocpl_id);
    }
  */

  // FIX ME
  /*
  if (type==4)
    {
      ocpl_id = H5Pcreate(H5P_OBJECT_COPY);
      sprintf(label,"%s/flag",groupName);
      status = H5Ocopy(inFile->fileID, label, outFile->fileID, label, ocpl_id, H5P_DEFAULT);
      H5Pclose(ocpl_id);
    }
  */

  if (type==0 || type==1 || type==4)
    {
      sprintf(groupName,"beam_%d/%s/astronomy_data",ibeam,blabel);
      sprintf(dsetName,"%s/data",groupName);
      dset_id = H5Dcreate2(outFile->fileID,dsetName,H5T_NATIVE_FLOAT,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);    
    }
  else if (type==2)
    {
      sprintf(groupName,"beam_%d/%s/calibrator_data",ibeam,blabel);
      sprintf(dsetName,"%s/cal_data_on",groupName);
      dset_id = H5Dcreate2(outFile->fileID,dsetName,H5T_NATIVE_FLOAT,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);    
    }
  else if (type==3)
    {
      sprintf(groupName,"beam_%d/%s/calibrator_data",ibeam,blabel);
      sprintf(dsetName,"%s/cal_data_off",groupName);
      dset_id = H5Dcreate2(outFile->fileID,dsetName,H5T_NATIVE_FLOAT,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);    
    }
  
  status  = H5Dwrite(dset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,out);  
  status = H5Dclose(dset_id);
  status = H5Sclose(dataspace_id);

  dims[0] = nchan;
  dataspace_id = H5Screate_simple(1,dims,NULL);

  if (type==0 || type==1 || type==4)
    {
      sprintf(groupName,"beam_%d/%s/astronomy_data",ibeam,blabel);
      sprintf(dsetName,"%s/frequency",groupName);
      dset_id = H5Dcreate2(outFile->fileID,dsetName,H5T_NATIVE_FLOAT,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);     
      status  = H5Dwrite(dset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,freq);  
      status  = H5Dclose(dset_id);
    }
  else if (type==2)
    {
      sprintf(groupName,"beam_%d/%s/calibrator_data",ibeam,blabel);
      sprintf(dsetName,"%s/cal_frequency",groupName);
      dset_id = H5Dcreate2(outFile->fileID,dsetName,H5T_NATIVE_FLOAT,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);     
      status  = H5Dwrite(dset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,freq);  
      status  = H5Dclose(dset_id);      
    }
  status = H5Sclose(dataspace_id);
  if (type == 0 || type == 1 || type == 4)
    status = H5Gclose(group_id);


}

// type = 0: copy all
void sdhdf_copyRemainder(sdhdf_fileStruct *inFile,sdhdf_fileStruct *outFile,int type)
{
  char groupName[1024];
  hid_t group_id;
  herr_t status;
  int b,i;
  
  // File information
  sprintf(groupName,"/metadata");
  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
    sdhdf_copyEntireGroup(groupName,inFile,outFile);
  else
    {
      sprintf(groupName,"/metadata/beam_params");
      if (sdhdf_checkGroupExists(outFile,groupName) == 1) sdhdf_copyEntireGroup(groupName,inFile,outFile);
      sprintf(groupName,"/metadata/history");
      if (sdhdf_checkGroupExists(outFile,groupName) == 1) sdhdf_copyEntireGroup(groupName,inFile,outFile);
      sprintf(groupName,"/metadata/primary_header");
      if (sdhdf_checkGroupExists(outFile,groupName) == 1) sdhdf_copyEntireGroup(groupName,inFile,outFile);
      sprintf(groupName,"/metadata/software_versions");
      if (sdhdf_checkGroupExists(outFile,groupName) == 1) sdhdf_copyEntireGroup(groupName,inFile,outFile);
    }
  
  // Config
  sprintf(groupName,"/config");
  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
    sdhdf_copyEntireGroup(groupName,inFile,outFile);


  // Beams
  for (b=0;b<inFile->nBeam;b++)
    {
      sprintf(groupName,"/beam_%d",b);
      if (sdhdf_checkGroupExists(outFile,groupName) == 1)
	sdhdf_copyEntireGroup(groupName,inFile,outFile);
      else
	{
	  sprintf(groupName,"/beam_%d/metadata",b);
	  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
	    sdhdf_copyEntireGroup(groupName,inFile,outFile);
	  else
	    {
	      sprintf(groupName,"/beam_%d/metadata/band_params",b);
	      if (sdhdf_checkGroupExists(outFile,groupName) == 1)
		sdhdf_copyEntireGroup(groupName,inFile,outFile);
	      sprintf(groupName,"/beam_%d/metadata/cal_band_params",b);
	      if (sdhdf_checkGroupExists(inFile,groupName) == 0 && sdhdf_checkGroupExists(outFile,groupName) == 1)
		sdhdf_copyEntireGroup(groupName,inFile,outFile);
	    }
	  
	  for (i=0;i<inFile->beam[b].nBand;i++)
	    {
	      sprintf(groupName,"/beam_%d/%s",b,inFile->beam[b].bandHeader[i].label);
	      if (sdhdf_checkGroupExists(outFile,groupName) == 1)
		sdhdf_copyEntireGroup(groupName,inFile,outFile);
	      else
		{
		  // Metadata
		  sprintf(groupName,"/beam_%d/%s/metadata",b,inFile->beam[b].bandHeader[i].label);
		  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
		    sdhdf_copyEntireGroup(groupName,inFile,outFile);
		  else
		    {
		      sprintf(groupName,"/beam_%d/%s/metadata/cal_obs_params",b,inFile->beam[b].bandHeader[i].label);
		      if (sdhdf_checkGroupExists(inFile,groupName) == 0 && sdhdf_checkGroupExists(outFile,groupName) == 1)
			sdhdf_copyEntireGroup(groupName,inFile,outFile);
		      sprintf(groupName,"/beam_%d/%s/metadata/obs_params",b,inFile->beam[b].bandHeader[i].label);
		      if (sdhdf_checkGroupExists(inFile,groupName) == 0 && sdhdf_checkGroupExists(outFile,groupName) == 1)
			sdhdf_copyEntireGroup(groupName,inFile,outFile);
		      
		    }
		  // calibrator_data
		  sprintf(groupName,"/beam_%d/%s/calibrator_data",b,inFile->beam[b].bandHeader[i].label);
		  if (sdhdf_checkGroupExists(inFile,groupName) == 0 && sdhdf_checkGroupExists(outFile,groupName) == 1)
		    sdhdf_copyEntireGroup(groupName,inFile,outFile);
		  
		  // astronomy_data
		  sprintf(groupName,"/beam_%d/%s/astronomy_data",b,inFile->beam[b].bandHeader[i].label);
		  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
		    sdhdf_copyEntireGroup(groupName,inFile,outFile);

		  // Flagging information
		  sprintf(groupName,"beam_%d/%s/astronomy_data/flag",b,inFile->beam[b].bandHeader[i].label);
		  if (sdhdf_checkGroupExists(inFile,groupName) == 0 && sdhdf_checkGroupExists(outFile,groupName) == 1)
		    sdhdf_copyEntireGroup(groupName,inFile,outFile);
		  
		}
	    }
	}
    }

}

void sdhdf_loadCalProc(sdhdf_fileStruct *inFile,int ibeam,int iband,char *cal_label,float *vals)
{
  int i;
  hid_t dataset_id;
  herr_t status;
  char dataName[1024];

  sprintf(dataName,"/beam_%d/%s/calibrator_proc/%s",ibeam,inFile->beam[ibeam].bandHeader[iband].label,cal_label);
  dataset_id   = H5Dopen2(inFile->fileID,dataName,H5P_DEFAULT);
  status = H5Dread(dataset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,vals);  
  status = H5Dclose(dataset_id);
}

void sdhdf_writeCalProc(sdhdf_fileStruct *outFile,int ibeam,int iband,char *band_label,char *cal_label,float *vals,int nchan,int npol,int ndumps)
{
  int i;
  hid_t dset_id,datatype_id,group_id;
  herr_t status;
  hid_t dataspace_id,stid;
  char groupName[MAX_STRLEN];
  char dsetName[MAX_STRLEN],label[MAX_STRLEN];
  hsize_t dims[5];

  // Create group if needed
  sprintf(groupName,"/beam_%d/%s/calibrator_proc",ibeam,band_label);

  // Check if group exists
  status = H5Eset_auto1(NULL,NULL);
  status = H5Gget_objinfo(outFile->fileID,groupName,0,NULL);
  if (status == 0)  
    {
      //printf("cal_proc already exists\n");
    }
  else
    {
      status=0;
      group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    }
  dims[0] = ndumps;
  dims[1] = 1;
  dims[2] = npol;
  dims[3] = nchan;
  dims[4] = 1;
  dataspace_id = H5Screate_simple(5,dims,NULL);
  
  sprintf(dsetName,"beam_%d/%s/calibrator_proc/%s",ibeam,band_label,cal_label);
  dset_id = H5Dcreate2(outFile->fileID,dsetName,H5T_NATIVE_FLOAT,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);    
  status  = H5Dwrite(dset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,vals);  
  status = H5Dclose(dset_id);
  status = H5Sclose(dataspace_id);
  status = H5Gclose(group_id);      

}


void sdhdf_copyEntireGroup(char *bandLabel,sdhdf_fileStruct *inFile,sdhdf_fileStruct *outFile)
{
  hid_t loc_id,ocpl_id,lcpl_id;
  herr_t status;

  loc_id = inFile->fileID;

  ocpl_id = H5Pcreate(H5P_OBJECT_COPY);
  lcpl_id = H5Pcreate(H5P_LINK_CREATE);


  status = H5Ocopy(loc_id, bandLabel, outFile->fileID, bandLabel, ocpl_id, lcpl_id);

  H5Pclose(ocpl_id);
  H5Pclose(lcpl_id);
}

void sdhdf_writeFlags(sdhdf_fileStruct *outFile,int ibeam,int iband,int *flag,int nchan,char *bandLabel)
{
  char dSetName[MAX_STRLEN];
  char groupName[MAX_STRLEN];
  hid_t dataspace_id,stid,dset_id,datatype_id,group_id;
  hsize_t dims[1];
  herr_t status;
  
  dims[0] = nchan;
  dataspace_id = H5Screate_simple(1,dims,NULL);

  sprintf(dSetName,"beam_%d/%s/astronomy_data/flag",ibeam,bandLabel);
  if (sdhdf_checkGroupExists(outFile,dSetName) == 1)
    dset_id = H5Dcreate2(outFile->fileID,dSetName,H5T_NATIVE_INT,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);    
  else
    dset_id = H5Dopen2(outFile->fileID,dSetName,H5P_DEFAULT);
  
  status  = H5Dwrite(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,flag);  
  status = H5Dclose(dset_id);
  status = H5Sclose(dataspace_id);
}
