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
// Software to join files together
//
// Usage:
//
// Compilation
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

void help()
{
  printf("sdhdf_join: routine to join sdhdf files together\n");
  printf("\n\n");
  printf("Command line arguments:\n\n");
  printf("-buildFreq         Build a new file with specific frequency bands (defined by -dsb)\n");
  printf("-dsb <bandNum>     Specify band number to join into output file (note this can be used multiple times)\n");
  printf("-h                 This help\n");
  printf("-joinFreq          Join in the frequency direction (by default we join in the time direction)\n");
  printf("-o <filename>      Output filename\n");
  printf("\n\n");
  printf("Input filenames are listed on the command line:\n");
  printf("sdhdf_join -o join.hdf uwl*.hdf\n");
  
}

int main(int argc,char *argv[])
{
  int joinFreq=0;
  int ii,i,j,k,kk,jj,l,totNchan,b;
  int buildFreq=0;
  long nchan_v;
  int npol_v;
  int nchan[MAX_BANDS],npol[MAX_BANDS];
  char fname[MAX_FILES][64];
  int nFiles=0;
  char ext[MAX_STRLEN]="extract";
  char oname[MAX_STRLEN];
  sdhdf_fileStruct *inFile,*outFile;
  herr_t status;
  int selectDump[MAX_BANDS];
  float *outVals,*freqVals,*inData;
  int outNdump[MAX_BANDS];
  int outNdumpCal[MAX_BANDS];
  int  copySD=0;
  float fSelect0[MAX_BANDS];
  float fSelect1[MAX_BANDS];
  int zoomBand=0;
  sdhdf_obsParamsStruct **outObsParams;
  sdhdf_obsParamsStruct *singleObsParams;
  sdhdf_obsParamsStruct **outObsParamsCal;
  sdhdf_bandHeaderStruct *outBandParams,*outCalBandParams;
  int nBand=0;
  int nBeam=0;
  char zoomLabel[MAX_STRLEN];
  char groupName[MAX_STRLEN];

  int defineSubband[MAX_BANDS];
  char defineSubbandFile[MAX_BANDS][MAX_STRLEN];
  int nDefineBands=0;
  
  sdhdf_attributes_struct dataAttributes[MAX_ATTRIBUTES];
  sdhdf_attributes_struct freqAttributes[MAX_ATTRIBUTES];
  int nDataAttributes=0;
  int nFreqAttributes=0;


  sdhdf_beamHeaderStruct *beamHeader;

  int haveCal;
  long ndump,ndumpCal,out_ndump,out_ndump_num;
  float **out_data,**out_freq;
  float *single_out_data,*single_out_freq;
  float **dataWts;
  int offsetDump[MAX_BANDS];
  int offsetPos[MAX_BANDS];
  int bandOffset;
  
  // Assume that all the input files have identical structure
  printf("Warning: assuming that all the input files have identical structure\n");
  
  strcpy(oname,"sdhdf_join_output.hdf");
  
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
      if (strcmp(argv[i],"-o")==0)
	strcpy(oname,argv[++i]);
      else if (strcmp(argv[i],"-h")==0)
	{help(); exit(1);}
      else if (strcmp(argv[i],"-buildFreq")==0)
	buildFreq=1;
      else if (strcasecmp(argv[i],"-joinFreq")==0)
	joinFreq=1;
      else if (strcmp(argv[i],"-dsb")==0)
	{
	  sscanf(argv[++i],"%d",&defineSubband[nDefineBands]);
	  strcpy(defineSubbandFile[nDefineBands],argv[++i]);
	  nDefineBands++;
	}      
      else
	{
	  strcpy(fname[nFiles],argv[i]);
	  nFiles++;
	}
    }


  if (buildFreq==1)
    {
      // Join file by extracting specific sub-bands from other files
      nBand = nDefineBands;

      
      //      for (i=0;i<nBand;i++)
      //	{
	  //      printf("band %d, output spectral dumps %d, nchan = %d, npol = %d\n",i,outNdump[i],nchan[i],npol[i]);
	  //	  outObsParams[i] = (sdhdf_obsParamsStruct *)malloc(sizeof(sdhdf_obsParamsStruct)*outNdump[i]);      
	  //	  out_data[i] = (float *)malloc(sizeof(float)*nchan[i]*npol[i]*outNdump[i]);
	  //	  out_freq[i] = (float *)malloc(sizeof(float)*nchan[i]);
	  //	}
      nBeam = 1; // HARDCODE -- FIX ME
      sdhdf_initialiseFile(outFile);
      sdhdf_openFile(oname,outFile,3);
      sdhdf_allocateBeamMemory(outFile,nBeam);
    
      for (b=0;b<nBeam;b++)
	{
	
	  outBandParams = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*nDefineBands);
	  for (ii=0;ii<nDefineBands;ii++)
	    {
	      printf("Loading file %s\n",defineSubbandFile[ii]);
	      sdhdf_initialiseFile(inFile);                
	      sdhdf_openFile(defineSubbandFile[ii],inFile,1); 
	      sdhdf_loadMetaData(inFile);

	      printf("%-30.30s %-22.22s %-5.5s %-6.6s %-20.20s %-10.10s %-10.10s %-7.7s %3d\n",fname[ii],inFile->primary[0].utc0,
		     inFile->primary[0].hdr_defn_version,inFile->primary[0].pid,inFile->beamHeader[0].source,
		     inFile->primary[0].telescope, inFile->primary[0].observer, inFile->primary[0].rcvr,inFile->beam[0].nBand);

	      
	      ndump = inFile->beam[b].bandHeader[defineSubband[ii]].ndump;
	      ndumpCal = inFile->beam[b].calBandHeader[defineSubband[ii]].ndump;
	      nchan_v = inFile->beam[b].bandHeader[defineSubband[ii]].nchan;
	      npol_v = inFile->beam[b].bandHeader[defineSubband[ii]].npol;

	      singleObsParams = (sdhdf_obsParamsStruct *)malloc(sizeof(sdhdf_obsParamsStruct)*ndump);
	      //	      singleObsParamsCal = (sdhdf_obsParamsStruct *)malloc(sizeof(sdhdf_obsParamsStruct)*ndumpCal);
	      single_out_data = (float *)malloc(sizeof(float)*nchan_v*npol_v*ndump);
	      single_out_freq = (float *)malloc(sizeof(float)*nchan_v);

	      
	      sdhdf_copyBandHeaderStruct(&(inFile->beam[b].bandHeader[defineSubband[ii]]),&outBandParams[ii],1);  
	      outBandParams[ii].ndump = ndump;
	      for (k=0;k<ndump;k++)
		{
		  sdhdf_copySingleObsParams(inFile,b,defineSubband[ii],k,&singleObsParams[k]);
		}
	      // Directly read the data into an array
	      sdhdf_loadBandData2Array(inFile,b,defineSubband[ii],1,single_out_data);
	      sdhdf_loadFrequency2Array(inFile,b,defineSubband[ii],single_out_freq);

	      // SHOULD COPY ATTRIBUTES -- FIX ME
	      if (ii!=nDefineBands-1)
		sdhdf_closeFile(inFile);

	      sdhdf_writeObsParams(outFile,outBandParams[ii].label,inFile->beamHeader[b].label,ii,singleObsParams,ndump,1);
	      sdhdf_writeSpectrumData(outFile,inFile->beamHeader[b].label,outBandParams[ii].label,b,ii,single_out_data,single_out_freq,nchan_v,npol_v,ndump,0,dataAttributes,nDataAttributes,freqAttributes,nFreqAttributes); // FIX THE ATTRIBUTES STUFF
	      free(singleObsParams);
	      free(single_out_data);
	      free(single_out_freq);

	    }
	  sdhdf_writeBandHeader(outFile,outBandParams,inFile->beamHeader[b].label,nDefineBands,1);
	
	  // SHOULD DO SOMETHING WITH THE CAL
	
	  
	      
	  free(outBandParams);

	  // Close last file
	  sdhdf_copyRemainder(inFile,outFile,0);
	  sdhdf_closeFile(inFile);
	
	}
     
    }
  else
    {
      // Identify how many spectral line dumps we'll have
      printf("Number of files: %d\n",nFiles);
      for (i=0;i<MAX_BANDS;i++)
	{
	  outNdump[i]=0;  
	  outNdumpCal[i]=0;
	  offsetDump[i]=0;
	  offsetPos[i]=0;
	}

      // Go through the files once to work out the structure of the output file
      nBand=0;
      bandOffset=0;
      for (i=0;i<nFiles;i++)
	{
	  printf("Processing file: %d %s\n",i,fname[i]);
	  sdhdf_initialiseFile(inFile);                
	  sdhdf_openFile(fname[i],inFile,1); // Open first file      
	  sdhdf_loadMetaData(inFile);
	  if (strcmp(inFile->primary[0].cal_mode,"OFF")==0)
	    haveCal=0;
	  else
	    haveCal=1;
	    
    
	  printf("Loaded metadata\n");
	  if (i==0)	
	    {
	      if (joinFreq==0)
		nBand = inFile->beam[0].nBand;
	      nBeam = inFile->nBeam;
	    }
	  if (joinFreq==1)
	    nBand += inFile->beam[0].nBand;
	  for (b=0;b<inFile->beam[0].nBand;b++)	
	    {
	      if (joinFreq==1)
		{
		  outNdump[b+bandOffset] = inFile->beam[0].bandHeader[b].ndump;
		  if (haveCal==1) outNdumpCal[b+bandOffset] = inFile->beam[0].calBandHeader[b].ndump;		  
		}
	      else
		{
		  outNdump[b] += inFile->beam[0].bandHeader[b].ndump;
		  if (haveCal==1) outNdumpCal[b] += inFile->beam[0].calBandHeader[b].ndump;
		}
	      if (i==0 && joinFreq==0)
		{
		  npol[b]  = inFile->beam[0].bandHeader[b].npol;
		  nchan[b] = inFile->beam[0].bandHeader[b].nchan;
		}
	      else if (joinFreq==1)
		{
		  npol[b+bandOffset]  = inFile->beam[0].bandHeader[b].npol;
		  nchan[b+bandOffset] = inFile->beam[0].bandHeader[b].nchan;
		}
	    }
	  bandOffset+=inFile->beam[0].nBand;
	  printf("Closing file\n");
	  sdhdf_closeFile(inFile);
	}
      //       for (b=0;b<inFile->beam[0].nBand;b++)
      //	 printf("outNdum for band %d is %d\n",b,outNdump[b]);
      outBandParams = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*nBand);
      if (haveCal==1)
	{
	  outCalBandParams = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*nBand);
	  outObsParamsCal = (sdhdf_obsParamsStruct **)malloc(sizeof(sdhdf_obsParamsStruct *)*nBand);
	}
      outObsParams = (sdhdf_obsParamsStruct **)malloc(sizeof(sdhdf_obsParamsStruct *)*nBand);
      out_data = (float **)malloc(sizeof(float *)*nBand);
      out_freq = (float **)malloc(sizeof(float *)*nBand);
      dataWts  = (float **)malloc(sizeof(float *)*nBand);
      beamHeader = (sdhdf_beamHeaderStruct *)malloc(sizeof(sdhdf_beamHeaderStruct)*nBeam);


      
      for (i=0;i<nBand;i++)
	{
	  outObsParams[i] = (sdhdf_obsParamsStruct *)malloc(sizeof(sdhdf_obsParamsStruct)*outNdump[i]);
	  if (haveCal==1)
	    outObsParamsCal[i] = (sdhdf_obsParamsStruct *)malloc(sizeof(sdhdf_obsParamsStruct)*outNdumpCal[i]);      
	  out_data[i] = (float *)malloc(sizeof(float)*nchan[i]*npol[i]*outNdump[i]);
	  out_freq[i] = (float *)malloc(sizeof(float)*nchan[i]);
	  dataWts[i]  = (float *)malloc(sizeof(float)*nchan[i]*outNdump[i]);
	}
      
      sdhdf_initialiseFile(outFile);
      sdhdf_openFile(oname,outFile,3);
      sdhdf_allocateBeamMemory(outFile,nBeam);
      //
      // Now produce the output file
      //
      for (b=0;b<nBeam;b++)
	{	 	  
	  bandOffset=0;
	  for (ii=0;ii<nFiles;ii++)
	    {
	      sdhdf_initialiseFile(inFile);                
	      sdhdf_openFile(fname[ii],inFile,1); 
	      sdhdf_loadMetaData(inFile);
	      printf("%-30.30s %-22.22s %-5.5s %-6.6s %-20.20s %-10.10s %-10.10s %-7.7s %3d\n",fname[ii],inFile->primary[0].utc0,
		     inFile->primary[0].hdr_defn_version,inFile->primary[0].pid,inFile->beamHeader[0].source,
		     inFile->primary[0].telescope, inFile->primary[0].observer, inFile->primary[0].rcvr,inFile->beam[0].nBand);
	      if (ii==0)
		{
		  strcpy(beamHeader[b].label,inFile->beamHeader[b].label);
		  beamHeader[b].nBand = nBand;
		  strcpy(beamHeader[b].source,inFile->beamHeader[b].source);
		}

	      
	      if (nBeam != inFile->nBeam)
		printf("WARNING: nBeam = %d != nBeam for file %s\n",nBeam,fname[ii]);
	      
	      if (ii==0 && joinFreq==0)
		{	       
		  sdhdf_copyBandHeaderStruct(inFile->beam[b].bandHeader,outBandParams,nBand);
		  if (haveCal==1)
		    sdhdf_copyBandHeaderStruct(inFile->beam[b].calBandHeader,outCalBandParams,nBand);
		  for (i=0;i<nBand;i++)
		    {
		      outBandParams[i].ndump = outNdump[i];
		      if (haveCal==1) outCalBandParams[i].ndump = outNdumpCal[i];
		    }
		      // FIX ME -- MUST UPDATE CAL HERE
		}
	      else if (joinFreq==1)
		{
		  sdhdf_copyBandHeaderStruct(inFile->beam[b].bandHeader,outBandParams+bandOffset,inFile->beam[b].nBand);
		  if (haveCal==1)
		    sdhdf_copyBandHeaderStruct(inFile->beam[b].calBandHeader,outCalBandParams+bandOffset,inFile->beam[b].nBand);
		  for (i=0;i<inFile->beam[b].nBand;i++)
		    {
		      outBandParams[i+bandOffset].ndump = outNdump[i];
		      if (haveCal==1)
			outCalBandParams[i+bandOffset].ndump = outNdumpCal[i];
		    }
		}
	      
	      for (i=0;i<inFile->beam[b].nBand;i++)
		{
		  ndump        = inFile->beam[b].bandHeader[i].ndump;
		  if (haveCal==1)
		    ndumpCal     = inFile->beam[b].calBandHeader[i].ndump;
		  
		  for (k=0;k<ndump;k++)
		    {

		      if (joinFreq==0)
			{
			  sdhdf_copySingleObsParams(inFile,b,i,k,&outObsParams[i][k+offsetDump[i]]);
			}
		      else
			{
			  sdhdf_copySingleObsParams(inFile,b,i,k,&outObsParams[i+bandOffset][k+offsetDump[i]]);
			}
			  // offsetPos is wrong here **
		      //		      dataWts[i][offsetPos[i]+k] = inFile->beam[b].bandData[i].astro_data.dataWeights[k*nchan[i]+];
		    }
		  if (haveCal==1)
		    {
		      for (k=0;k<ndumpCal;k++)
			{
			  if (joinFreq==0)
			    {
			      sdhdf_copySingleObsParamsCal(inFile,b,i,k,&outObsParamsCal[i][k+offsetDump[i]]);			  
			    }
			  else
			    {
			      sdhdf_copySingleObsParamsCal(inFile,b,i,k,&outObsParamsCal[i+bandOffset][k+offsetDump[i]]);
			    }
			}
		    }
		  // Directly read the data into an array
		  if (joinFreq==0)
		    sdhdf_loadBandData2Array(inFile,b,i,1,&out_data[i][offsetPos[i]]);
		  //		  else
		  //		    sdhdf_loadBandData2Array(inFile,b,i,1,&out_data[i+bandOffset][0]);		 // FIX ME -- 0 is odd here
		  //		  for (j=0;j<20;j++)
		  //		    printf("Loaded %g\n",out_data[i][offsetPos[i]+j]);
		  if (ii==0 && joinFreq==0)
		    {
		      // Directly read the frequency data into an array
		      sdhdf_loadFrequency2Array(inFile,b,i,out_freq[i]);
		    }
		  else if (joinFreq==1)
		    {
		      sdhdf_loadFrequency2Array(inFile,b,i,out_freq[i+bandOffset]);
		    }
		  if (joinFreq==0)
		    {
		      offsetPos[i] +=ndump*nchan[i]*npol[i];
		      offsetDump[i]+=ndump;
		    }
		  else
		    offsetPos[i+bandOffset] +=ndump*nchan[i+bandOffset]*npol[i+bandOffset];
		  if (ii==0 && joinFreq==0) // Note that this is getting the attributes from the first file -- really should check if they are the same in all the files
		    {
		      
		      // Need to do this separately for each band *** FIX ME **** CURRENTLY IT WILL ONLY TAKE THE LAST ONE
		      sdhdf_copyAttributes(inFile->beam[b].bandData[i].astro_obsHeaderAttr,inFile->beam[b].bandData[i].nAstro_obsHeaderAttributes,dataAttributes,&nDataAttributes);
		      sdhdf_copyAttributes(inFile->beam[b].bandData[i].astro_obsHeaderAttr_freq,inFile->beam[b].bandData[i].nAstro_obsHeaderAttributes_freq,freqAttributes,&nFreqAttributes);	  
		      //		  printf("Have copied %d attributes\n",nFreqAttributes);
		    }
		  else if (joinFreq==1)
		    {
		      sdhdf_copyAttributes(inFile->beam[b].bandData[i].astro_obsHeaderAttr,inFile->beam[b].bandData[i].nAstro_obsHeaderAttributes,dataAttributes,&nDataAttributes);
		      sdhdf_copyAttributes(inFile->beam[b].bandData[i].astro_obsHeaderAttr_freq,inFile->beam[b].bandData[i].nAstro_obsHeaderAttributes_freq,freqAttributes,&nFreqAttributes);

		    }
		}
	      bandOffset+=inFile->beam[b].nBand;
	      sdhdf_closeFile(inFile);
	    }
	  sdhdf_writeBeamHeader(outFile,beamHeader,nBeam);		
	  sdhdf_writeBandHeader(outFile,outBandParams,b,nBand,1);
	  if (haveCal==1)
	    sdhdf_writeBandHeader(outFile,outCalBandParams,b,nBand,2);

	  for (i=0;i<nBand;i++)
	    {
	      sdhdf_writeObsParams(outFile,outBandParams[i].label,inFile->beamHeader[b].label,i,outObsParams[i],outNdump[i],1);
	      if (haveCal==1)
		sdhdf_writeObsParams(outFile,outBandParams[i].label,inFile->beamHeader[b].label,i,outObsParamsCal[i],outNdumpCal[i],2);

	      // NOTE THAT THE ATTRIBUTES ARE INCORRECT HERE -- SEE COMMENT ABOVE ABOUT ONLY TAKING THE LAST ONE ***** FIX ME
	      //	      printf("Writing data for band %d (%s) ndump = %d, npol = %d, nchan = %d\n",i,outBandParams[i].label,outNdump[i],npol[i],nchan[i]);
	      if (joinFreq==0)
		sdhdf_writeSpectrumData(outFile,inFile->beamHeader[b].label,outBandParams[i].label,b,i,out_data[i],out_freq[i],nchan[i],npol[i],outNdump[i],0,dataAttributes,nDataAttributes,freqAttributes,nFreqAttributes);
	    }
	  // *********** DON'T FORGET THE CAL *********
	  
	}
      
      // Re-load the first file and copy the remainder
      if (joinFreq==1)
	{
	  for (i=0;i<nFiles;i++)
	    {
	      sdhdf_initialiseFile(inFile);                
	      
	      sdhdf_openFile(fname[i],inFile,1); 
	      sdhdf_loadMetaData(inFile);
	      sdhdf_copyRemainder(inFile,outFile,0);
	      sdhdf_closeFile(inFile);
	    }
	    }
      else
	{
	  sdhdf_initialiseFile(inFile);                
	  
	  sdhdf_openFile(fname[0],inFile,1); 
	  sdhdf_loadMetaData(inFile);
	  sdhdf_copyRemainder(inFile,outFile,0);
	  sdhdf_closeFile(inFile);

	}
      sdhdf_closeFile(outFile);
      free(outBandParams);
      if (haveCal==1)
	free(outCalBandParams);
      for (i=0;i<nBand;i++)
	{
	  free(outObsParams[i]);
	  if (haveCal==1)
	    free(outObsParamsCal[i]);
	  free(out_data[i]);
	  free(out_freq[i]);
	  
	}
      free(outObsParams);
      if (haveCal==1)
	free(outObsParamsCal);
      free(out_data);
      free(out_freq);
      free(beamHeader);
    }

  free(inFile);
  free(outFile);
}



