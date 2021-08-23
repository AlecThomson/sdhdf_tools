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

#define MAX_RFI 512

void saveFile(sdhdf_fileStruct *inFile,char *extension);

void help()
{
  printf("sdhdf_autoFlag:   code to carry out automatic flagging of a data file\n");
  printf("Note 1: sdhdf_flag is used for manual flagging\n");
  printf("Note 2: here we use the term 'flagging', but it is the weights table in the SDHDF file that is set to zero\n");
  printf("\n\n");
  printf("-e <extension>            Set the output file extension\n");
  printf("-from <filename.hdf>      Input file that has flagging information that should be copied to other files\n");
  printf("-edge <value>             Percent of the sub-band boundaries that should be flagged\n");
  printf("-persistent               Remove persistent RFI specific to the observatory\n");
  printf("-transient                Remove transient RFI specific to the observatory\n");
  printf("\n");
  printf("Example: sdhdf_autoFlag -from manualZap.hdf uwl*.hdf\n");
	 

  exit(1);
}


int main(int argc,char *argv[])
{
  int i,j,ii,b,k,s;
  char fname[MAX_STRLEN];
  char fromFileName[MAX_STRLEN];
  char extension[MAX_STRLEN]="autoFlag";
  sdhdf_fileStruct *inFile,*fromFile;
  int iband=0,nchan,idump,nband;
  int ndump=0,totSize,chanPos;
  //  spectralDumpStruct spectrum;
  int npol=4;
  int setnband=-1;
  int flagType=-1;
  sdhdf_rfi rfi[MAX_RFI];
  sdhdf_transient_rfi transient_rfi[MAX_RFI];
  int nTransientRFI=0;
  int nRFI;
  float bandEdge=0;
  int fromFlag=0;
  int bandEdgeFlag=0;
  int persistentFlag=0;
  int transientFlag=0;
  
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
  flagType=0;
  for (i=1;i<argc;i++)
    {      
      if (strcmp(argv[i],"-from")==0)	{strcpy(fromFileName,argv[++i]); flagType=1; fromFlag=1;}
      else if (strcmp(argv[i],"-e")==0) strcpy(extension,argv[++i]);
      else if (strcmp(argv[i],"-edge")==0) {sscanf(argv[++i],"%f",&bandEdge); bandEdgeFlag=1;}
      else if (strcmp(argv[i],"-persistent")==0) {persistentFlag=1; flagType=2;}
      else if (strcmp(argv[i],"-transient")==0) {transientFlag=1; flagType=2;}
      else if (strcmp(argv[i],"-h")==0) help();
    }

  if (fromFlag==0 && bandEdgeFlag==0 && persistentFlag==0 && transientFlag==0)
    {
      printf("ERROR: You need to use either the -from, -edge, -transient or -persistent command line options\n");
      exit(1);
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
      if (strcmp(argv[ii],"-from")==0 || strcmp(argv[ii],"-edge")==0 || strcmp(argv[ii],"-e")==0)
	ii++;
      else if (strcmp(argv[ii],"-persistent")==0 || strcmp(argv[ii],"-transient")==0)
	{
	  // Do nothing
	}
      else // Processing file
	{
	  printf("ii = %d\n",ii);
	  strcpy(fname,argv[ii]);
	  printf("Starting %s\n",fname);  	  
	  sdhdf_initialiseFile(inFile);
	  sdhdf_openFile(fname,inFile,1);
	  sdhdf_loadMetaData(inFile);
	  if (flagType==2)
	    {
	      printf("Observatory site = %s\n",inFile->primary[0].telescope);
	      if (strcmp(inFile->primary[0].telescope,"Parkes")!=0)
		printf("WARNING: ONLY IMPLEMENTED PARKES OBSERVATORY\n"); // FIX ME	      
	      if (transientFlag==1)
		  sdhdf_loadTransientRFI(transient_rfi,&nTransientRFI,MAX_RFI,"parkes");
	      if (persistentFlag==1)
		sdhdf_loadPersistentRFI(rfi,&nRFI,MAX_RFI,"parkes");

	      if (bandEdge > 0)
		{
		  double freq0,freq1,chbw;
		  printf("In band edge\n");
		  nband = inFile->beam[0].nBand; // FIX ME - USING BEAM 0
		  printf("nband = %d\n",nband);
		  for (i=0;i<nband;i++)
		    {
		      // FIX ME -- THIS IS CRAZY TO LOAD IN THE DATA EACH TIME JUST TO SET THE FLAGS **** FIX ME *****
		      
		      sdhdf_loadBandData(inFile,0,i,1); // FIX ME - USING BEAM 0
		  
		      printf("Setting RFI %d\n",nRFI);
		      freq0 = inFile->beam[0].bandData[i].astro_data.freq[0];
		      freq1 = inFile->beam[0].bandData[i].astro_data.freq[inFile->beam[0].bandHeader[i].nchan-1];
		      chbw = inFile->beam[0].bandData[i].astro_data.freq[1]-inFile->beam[0].bandData[i].astro_data.freq[0];
		      printf("Got %g %g\n",freq0,freq1);
		      rfi[nRFI].type=1;  rfi[nRFI].f0 = freq0-chbw/2.; rfi[nRFI].f1 = freq0+(freq1-freq0)*bandEdge/100.+chbw/2.; nRFI++;
		      rfi[nRFI].type=1;  rfi[nRFI].f0 = freq1 -chbw/2.-(freq1-freq0)*bandEdge/100.; rfi[nRFI].f1 = freq1+chbw/2.; nRFI++;
		      printf("Adding RFI at %g %g\n",freq0,freq0+(freq1-freq0)*bandEdge/100.);
		      sdhdf_releaseBandData(inFile,b,i,1);

		    }
		}
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

		      if (persistentFlag==1)
			{
			  // ii, b, i (band)
			  for (k=0;k<nRFI;k++)
			    {
			      // Is this RFI within the band?
			      if ((rfi[k].f0 > inFile->beam[b].bandHeader[i].f0 && rfi[k].f0 < inFile->beam[b].bandHeader[i].f1)
				  || (rfi[k].f1 > inFile->beam[b].bandHeader[i].f0 && rfi[k].f1 < inFile->beam[b].bandHeader[i].f1))
				{
				  // Now zap from each sub-band
				  printf("Removing persistent RFI: %s\n",rfi[k].description);
				  for (s=0;s<inFile->beam[b].bandHeader[i].ndump;s++)
				    {
				      for (j=0;j<inFile->beam[b].bandHeader[i].nchan;j++)
					{
					  freq = inFile->beam[b].bandData[i].astro_data.freq[j];
					  if (freq > rfi[k].f0 && freq < rfi[k].f1)
					    inFile->beam[b].bandData[i].astro_data.dataWeights[j+s*inFile->beam[b].bandHeader[i].nchan] = 0;
					}
				    }
				}
			    }
			}

		      if (transientFlag==1)
			{
			  double v1,v2,v3;
			  int n1,n2;
			  int flagIt=0;
			  double sumPol=0;
			  // ii, b, i (band)
			  for (k=0;k<nTransientRFI;k++)
			    {
			      // Is this RFI within the band?
			      if ((transient_rfi[k].f0 > inFile->beam[b].bandHeader[i].f0 && transient_rfi[k].f0 < inFile->beam[b].bandHeader[i].f1)
				  || (transient_rfi[k].f1 > inFile->beam[b].bandHeader[i].f0 && transient_rfi[k].f1 < inFile->beam[b].bandHeader[i].f1))
				{
				  // Now zap from each sub-band
				  printf("Checking transient RFI: %s comparing %g -> %g with %g -> %g\n",transient_rfi[k].description,
					 transient_rfi[k].f0,transient_rfi[k].f1,transient_rfi[k].f2,transient_rfi[k].f3);
				  for (s=0;s<inFile->beam[b].bandHeader[i].ndump;s++)
				    {
				      v1 = v2 = v3 = 0;
				      
				      n1=n2=0;
				      flagIt=0;
				      for (j=0;j<inFile->beam[b].bandHeader[i].nchan;j++)
					{
					  freq = inFile->beam[b].bandData[i].astro_data.freq[j];
					  if (freq > transient_rfi[k].f0 && freq < transient_rfi[k].f1)
					    {
					      sumPol = (inFile->beam[b].bandData[i].astro_data.pol1[j+s*inFile->beam[b].bandHeader[i].nchan] 
							+ inFile->beam[b].bandData[i].astro_data.pol2[j+s*inFile->beam[b].bandHeader[i].nchan]);
					      if (v1 < sumPol) v1 = sumPol;
					      n1 ++;
					    }
					  if (freq > transient_rfi[k].f2 && freq < transient_rfi[k].f3)
					    {
					      sumPol = (inFile->beam[b].bandData[i].astro_data.pol1[j+s*inFile->beam[b].bandHeader[i].nchan] 
							+ inFile->beam[b].bandData[i].astro_data.pol2[j+s*inFile->beam[b].bandHeader[i].nchan]);
					      v2 +=  sumPol;
					      n2 ++;
					    }
					}
				      v2/=(double)n2; // SHOULD CHECK n2 > 0 --- FIX ME


				      if (v1 > transient_rfi[k].threshold*v2)
					{
					  flagIt=1;
					  printf("Flagging dump %d, comparing %g and %g (%d %d)\n",s,v1,v2,n1,n2);
					}
				      else
					printf("**Not** flagging dump %d, comparing %g and %g (%d %d)\n",s,v1,v2,n1,n2);
				      if (flagIt==1) // Do the flagging if necessary
					{

					  for (j=0;j<inFile->beam[b].bandHeader[i].nchan;j++)
					    {					      
					      freq = inFile->beam[b].bandData[i].astro_data.freq[j];
					      if (freq > transient_rfi[k].f0 && freq < transient_rfi[k].f1)
						inFile->beam[b].bandData[i].astro_data.dataWeights[j+s*inFile->beam[b].bandHeader[i].nchan] = 0;
					    }					  
					}
				    }
				}
			    }
			}
/*
		      for (s=0;s<inFile->beam[b].bandHeader[i].ndump;s++)
			{
			  printf("Processing spectral dump %d\n",s);
			  for (j=0;j<inFile->beam[b].bandHeader[i].nchan;j++)
			    {
			      freq = inFile->beam[b].bandData[i].astro_data.freq[j];
			      if (persistentFlag==1)
				{
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
			      if (transientFlag==1)
				{
				  for (k=0;k<nTransientRFI;k++) // CAN DO THIS MUCH QUICKER			   
				    {
				      // WHAT ABOUT RFI ACROSS SUB-BAND BOUNDARIES?
				      //			      printf("Checking %.5f %g %g\n",freq,rfi[k].f0,rfi[k].f1);
				      if (freq > transient_rfi[k].f0 && freq < transient_rfi[k].f1)
					{
					  inFile->beam[b].bandData[i].astro_data.dataWeights[j+s*inFile->beam[b].bandHeader[i].nchan] = 0;
					  break;
					}
				    }

				}
			    }
			}
		      */
		    }
		  sdhdf_releaseBandData(inFile,b,i,1);
		}
	    }
	  saveFile(inFile,extension);
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

void saveFile(sdhdf_fileStruct *inFile,char *extension)
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
  sprintf(oname,"%s.%s",inFile->fname,extension);

  
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

