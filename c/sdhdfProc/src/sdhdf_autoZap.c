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

void saveFile(fileStruct *inFile);

int main(int argc,char *argv[])
{
  int i,j,ii;
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
  sdhdf_loadPrimaryMetadata(fromFile);  
  sdhdf_loadBandMetadata(fromFile,1);
  sdhdf_loadSpectralDumpMetadata(fromFile);
  sdhdf_loadFlagData(fromFile);

  //
  for (ii=1;ii<argc;ii++)
    {
      if (strcmp(argv[ii],"-from")!=0)
	{
	  strcpy(fname,argv[ii]);
	  printf("Processing %s\n",fname);
	  sdhdf_initialiseFile(inFile);
	  sdhdf_openFile(fname,inFile,1);
	  sdhdf_loadPrimaryMetadata(inFile);  
	  sdhdf_loadBandMetadata(inFile,1);	  
	  sdhdf_loadSpectralDumpMetadata(inFile);
	  sdhdf_loadFlagData(inFile);

	  
	  nband = inFile->nBands;
	  for (i=0;i<nband;i++)
	    {
	      for (j=0;j<inFile->bandHeader[i].nchan;j++)
		inFile->bandData[i].flag[j] = fromFile->bandData[i].flag[j];
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

void saveFile(fileStruct *inFile)
{
  char oname[MAX_STRLEN];
  char flagName[MAX_STRLEN];
  fileStruct *outFile;
  int i,j;
  hsize_t dims[1];
  hid_t dset_id,dataspace_id;
  herr_t status;
  int *outFlags;

  // Should check if already .flag extension
  //
  sprintf(oname,"%s.autoflag",inFile->fname);

  
  if (!(outFile = (fileStruct *)malloc(sizeof(fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >outFile<\n");
      exit(1);
    }

  sdhdf_initialiseFile(outFile);
  sdhdf_openFile(oname,outFile,3);

  sdhdf_copyEntireGroup("band_header",inFile,outFile);
  sdhdf_copyEntireGroup("obs_config",inFile,outFile);	      	      
  sdhdf_copyEntireGroup("cal_obs_params",inFile,outFile);
  sdhdf_copyEntireGroup("history",inFile,outFile); // SHOULD UPDATE THIS GROUP
  sdhdf_copyEntireGroup("obs_params",inFile,outFile);
  sdhdf_copyEntireGroup("primary_header",inFile,outFile);
  sdhdf_copyEntireGroup("software_versions",inFile,outFile);
  for (i=0;i<inFile->nBands;i++)
    {
      sdhdf_copyEntireGroup(inFile->bandHeader[i].label,inFile,outFile);
      // Add in the flags here if not already set - if so then update

      outFlags = (int *)malloc(sizeof(int)*inFile->bandHeader[i].nchan);
      for (j=0;j<inFile->bandHeader[i].nchan;j++)
	outFlags[j] = (int)inFile->bandData[i].flag[j];
      
      dims[0] = inFile->bandHeader[i].nchan;
      dataspace_id = H5Screate_simple(1,dims,NULL);

      sprintf(flagName,"%s/flag",inFile->bandHeader[i].label);
      if (H5Lexists(outFile->fileID,flagName,H5P_DEFAULT) <= 0)
	dset_id = H5Dcreate2(outFile->fileID,flagName,H5T_NATIVE_INT,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      else
	dset_id = H5Dopen2(outFile->fileID,flagName,H5P_DEFAULT);
      status  = H5Dwrite(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,outFlags);
      status  = H5Dclose(dset_id);
      status  = H5Sclose(dataspace_id);  
      free(outFlags);
    }
  sdhdf_closeFile(outFile);

  free(outFile);
}

