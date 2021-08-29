//  Copyright (C) 2021 George Hobbs

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

/* To include

- Background sky from Haslam or equivalent
- HI from GASS or equivalent
- Able to include bright continuum sources (Haslam update or NVSS/SUMSS?)
/*
# DON'T FORGET:
#
# BANDPASS - SLOPES
# RIPPLES
# GAIN VARIATIONS
# TONES AND COMBS
# PARALLACTIFY??
*/

- Include bandpass shape/slope across band
- Correctly model gain in different pols and allow gain variations
- Include ripples
- Include RFI
      - impulsive (aircraft/satellites)
- Have ability to convert to counts (or Kelvin)
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sdhdfProc.h"
#include "T2toolkit.h"
#include "hdf5.h"


typedef struct dumpParameterStruct {
  float tdump; // Time of spectral dump
  float raj0;   // Pointing position in degrees at centre of spectral dump (RA)
  float decj0;   // Pointing position in degrees at centre of spectral dump (DEC)

  double timeFromStart;
  double mjd;
} dumpParameterStruct;

typedef struct bandParameterStruct {
  char  label[MAX_STRLEN];    // Band label
  float f0;
  float f1;
  int   nchan;
} bandParameterStruct;

typedef struct beamParameterStruct {
  char  label[MAX_STRLEN];   // Beam label
  float tsys_aa;             // System temperature for this beam for AA polarisation (K)
  float tsys_bb;             // System temperature for this beam for BB polarisation (K)
  float delta_ra;            // Offset of beam from nominal pointing direction in RA (deg)
  float delta_dec;           // Offset of beam from nominal pointing direction in DEC (deg)
  char  src[MAX_STRLEN];     // Source name
} beamParameterStruct;

typedef struct parameterStruct {
  int nbeams;
  int nbands;
  int ndumps;
  int npols;
  char fname[MAX_STRLEN];
  double mjd0;
  char pid[MAX_STRLEN];
  char observer[MAX_STRLEN];
  char telescope[MAX_STRLEN];
  char rcvr[MAX_STRLEN];
  beamParameterStruct *beam;
  bandParameterStruct *band;
  dumpParameterStruct *dump;

  int persistentRFI;
  int galacticHI;
} parameterStruct;

int main(int argc,char *argv[])
{
  int i,j,k,ii,jj;
  long iseed = TKsetSeed();
  char fname[MAX_STRLEN]="unset";
  char paramFile[MAX_STRLEN]="unset";
  char labelStr[MAX_STRLEN];
  sdhdf_fileStruct          *outFile;
  sdhdf_beamHeaderStruct    *beamHeader;
  sdhdf_primaryHeaderStruct *primaryHeader;
  sdhdf_bandHeaderStruct    *bandHeader;
  sdhdf_obsParamsStruct     *obsParams;
  double dt;
  int nbeam=0;
  int ibeam=0;
  int ndump=0;
  int idump,iband;
  int nband;
  int nchan,nsub,npol;
  int nDataAttributes=0;
  int nFreqAttributes=0;
  sdhdf_attributes_struct *freqAttributes;
  sdhdf_attributes_struct *dataAttributes;
  float *freq,*data;
  float ssys_aa = 22;
  float ssys_bb = 23;
  float chbw;
  float req_chbw = 0.5e3; // in Hz
  float band_bw  = 192; // in MHz
  float rcvr_f0  = 700.0; // MHz
  float rcvr_f1  = 1852.0; // MHz. How will the band be split ??
  float signal;
  FILE *fin;
  char line[MAX_STRLEN];
  char word1[MAX_STRLEN];
  char word2[MAX_STRLEN];
  parameterStruct *params;
  
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-p")==0)
	strcpy(paramFile,argv[++i]);
    }

  if (strcmp(paramFile,"unset")==0)
    {
      printf("ERROR: must set a parameter filename using the -p option\n");
      exit(1);
    }

  if (!(fin = fopen(paramFile,"r")))
    {
      printf("ERROR: unable to open parameter file >%s<\n",paramFile);
      exit(1);
    }
  params = (parameterStruct *)malloc(sizeof(parameterStruct));
  // initialise
  params->persistentRFI=0;
  params->galacticHI=0;
  
  while (!feof(fin))
    {
      if (fgets(line,MAX_STRLEN,fin)!=NULL)
	{
	  if (line[0]!='#' && strlen(line) > 1) // Not a comment line
	    {
	      strcpy(word2,"NULL");
	      if (sscanf(line,"%s %s",word1,word2)==2)
		{
		  if (strcmp(word1,"nbeams:")==0)
		    {
		      sscanf(word2,"%d",&(params->nbeams));
		      // Allocate memory for this number of beams
		      params->beam = (beamParameterStruct *)malloc(sizeof(beamParameterStruct)*params->nbeams);
		      ibeam=0;
		    }
		  else if (strcmp(word1,"persistent_rfi:")==0 && strcmp(word2,"on")==0)
		    params->persistentRFI=1;
		  else if (strcmp(word1,"galactic_hi:")==0 && strcmp(word2,"gaussian")==0)
		    params->galacticHI=1;
		  else if (strcmp(word1,"beam:")==0)
		    {
		      if (ibeam == params->nbeams)
			printf("WARNING: TOO MANY BEAMS SPECIFIED. IGNORING >%s<\n",line);
		      else
			{
			  sscanf(line,"%s %s %s %f %f %f %f",word1,params->beam[ibeam].label,params->beam[ibeam].src,
				 &(params->beam[ibeam].tsys_aa),
				 &(params->beam[ibeam].tsys_bb),&(params->beam[ibeam].delta_ra),&(params->beam[ibeam].delta_dec));
			  ibeam++;
			}
		    }
		  else if (strcmp(word1,"nbands:")==0)
		    {
		      sscanf(word2,"%d",&(params->nbands));
		      params->band = (bandParameterStruct *)malloc(sizeof(bandParameterStruct)*params->nbands);
		      iband=0;
		    }
		  else if (strcmp(word1,"band:")==0)
		    {
		      if (iband == params->nbands)
			{
			  printf("WARNING: TOO MANY BANDS SPECIFIED. IGNORING >%s<\n",line);
			  printf("word1 = %s, word2 = %s\n",word1,word2);
			  printf("iband = %d, nbands = %d\n",iband,params->nbands);
			}
			  else
			{
			  sscanf(line,"%s %s %f %f %d",word1,params->band[iband].label,
				 &(params->band[iband].f0), &(params->band[iband].f1),&(params->band[iband].nchan));
			  iband++;
			}
		    }
		  else if (strcmp(word1,"ndumps:")==0)
		    {
		      sscanf(word2,"%d",&(params->ndumps));
		      params->dump = (dumpParameterStruct *)malloc(sizeof(dumpParameterStruct)*params->ndumps);
		      idump=0;
		    }
		  else if (strcmp(word1,"dump:")==0)
		    {
		      if (idump == params->ndumps)
			  printf("WARNING: TOO MANY DUMPS SPECIFIED. IGNORING >%s<\n",line);
		      else
			{
			  sscanf(line,"%s %f %f %f",word1,
				 &(params->dump[idump].tdump), &(params->dump[idump].raj0),&(params->dump[idump].decj0));
			  idump++;
			}
		    }
		  else if (strcmp(word1,"npols:")==0)
		    sscanf(word2,"%d",&(params->npols));
		  else if (strcmp(word1,"start_mjd:")==0)
		    sscanf(word2,"%lf",&(params->mjd0));
		  else if (strcmp(word1,"output:")==0)
		    strcpy(params->fname,word2);
		  else if (strcmp(word1,"pid:")==0)
		    strcpy(params->pid,word2);
		  else if (strcmp(word1,"observer:")==0)
		    strcpy(params->observer,line+10);
		  else if (strcmp(word1,"receiver:")==0)
		    strcpy(params->rcvr,word2);
		  else if (strcmp(word1,"telescope:")==0)
		    strcpy(params->telescope,word2);
		}
	    }
	}
    }
  fclose(fin);

  // Calculate times
  dt=0;
  for (i=0;i<params->ndumps;i++)
    {
      dt+=params->dump[i].tdump/2.;
      params->dump[i].timeFromStart = dt;
      params->dump[i].mjd = params->mjd0 + dt/86400.0;
      dt+=params->dump[i].tdump/2.;
    }

  
  
  if (!(outFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }
  sdhdf_initialiseFile(outFile);
  if (sdhdf_openFile(params->fname,outFile,3)==-1)
    {
      printf("Unable to open output file >%s<\n",fname);
      free(outFile);
      exit(1);
    }

  nbeam = params->nbeams;
  nband = params->nbands; 
  printf("Simulating %d bands\n",nband);
  
  beamHeader = (sdhdf_beamHeaderStruct *)malloc(sizeof(sdhdf_beamHeaderStruct)*nbeam);
  primaryHeader = (sdhdf_primaryHeaderStruct *)malloc(sizeof(sdhdf_primaryHeaderStruct));

  for (i=0;i<nbeam;i++)
    {
      strcpy(beamHeader[i].label,params->beam[i].label);
      beamHeader[i].nBand = nband;
      strcpy(beamHeader[i].source,params->beam[i].src); 
    }
  sdhdf_writeBeamHeader(outFile,beamHeader,nbeam);

  // Write the primary header
  strcpy(primaryHeader[0].date,"unknown");
  strcpy(primaryHeader[0].hdr_defn,"5");
  strcpy(primaryHeader[0].hdr_defn_version,"1.9");
  strcpy(primaryHeader[0].file_format,"unknown");
  strcpy(primaryHeader[0].file_format_version,"unknown");
  primaryHeader[0].sched_block_id = 1;
  strcpy(primaryHeader[0].cal_mode,"OFF");
  strcpy(primaryHeader[0].instrument,"simulate");
  strcpy(primaryHeader[0].observer,params->observer);
  strcpy(primaryHeader[0].pid,params->pid);
  strcpy(primaryHeader[0].rcvr,params->rcvr);
  strcpy(primaryHeader[0].telescope,params->telescope);
  strcpy(primaryHeader[0].utc0,"unknown");
  primaryHeader[0].nbeam = nbeam;
  sdhdf_writePrimaryHeader(outFile,primaryHeader);


  // Simulate data for each beam
  for (i=0;i<nbeam;i++)
    {
      bandHeader = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*nband); // FIX ME -- may have different numbers of bands 
      for (j=0;j<nband;j++)
	{
	  ndump = params->ndumps;
	  npol  = params->npols;
	  
	  bandHeader[j].f0 = params->band[j].f0;
	  bandHeader[j].f1 = params->band[j].f1;
	  bandHeader[j].fc = (bandHeader[j].f0+bandHeader[j].f1)/2.0;

	  nchan = params->band[j].nchan;

	  obsParams = (sdhdf_obsParamsStruct *)malloc(sizeof(sdhdf_obsParamsStruct)*ndump);
	  freq = (float *)malloc(sizeof(float)*nchan);
	  data = (float *)malloc(sizeof(float)*nchan*ndump*npol);
	  strcpy(bandHeader[j].label,params->band[j].label);

	  bandHeader[j].nchan = nchan;
	  bandHeader[j].npol = npol;
	  if (npol==4)
	    strcpy(bandHeader[j].pol_type,"AABBCRCI");
	  else if (npol==2)
	    strcpy(bandHeader[j].pol_type,"AA,BB");
	  else if (npol==1)
	    strcpy(bandHeader[j].pol_type,"AA+BB");
	  bandHeader[j].dtime = params->dump[0].tdump; 
	  bandHeader[j].ndump = ndump;	 
	  chbw = (bandHeader[j].f1 - bandHeader[j].f0)/(float)nchan;
	  // Setup the observation parameters for each band
	  for (k=0;k<nchan;k++)
	    freq[k] = bandHeader[j].f0+k*(bandHeader[j].f1-bandHeader[j].f0)/(float)nchan + chbw/2.0;
	  
	  for (k=0;k<ndump;k++)
	    {
	      obsParams[k].timeElapsed = params->dump[k].timeFromStart;
	      strcpy(obsParams[k].timedb,"unknown");
	      obsParams[k].mjd = params->dump[k].mjd;
	      strcpy(obsParams[k].utc,"unknown");
	      strcpy(obsParams[k].ut_date,"unknown");
	      strcpy(obsParams[k].aest,"aest_unknown");
	      strcpy(obsParams[k].raStr,"unknown");
	      strcpy(obsParams[k].decStr,"unknown");
	      obsParams[k].raOffset = 0;
	      obsParams[k].decOffset = 0;
	      obsParams[k].gl = 0;
	      obsParams[k].gb = 0;
	      obsParams[k].az = 0;
	      obsParams[k].el = 0;
	      obsParams[k].az_drive_rate = 0;
	      obsParams[k].ze_drive_rate = 0;
	      obsParams[k].hourAngle = 0;
	      obsParams[k].paraAngle = 0;
	      obsParams[k].windDir = 0;
	      obsParams[k].windSpd = 0;		      
	      for (ii=0;ii<nchan;ii++)
		{
		  if (params->galacticHI==1)
		    signal = 40*exp(-pow(freq[ii]-1421.0,2)/2./0.2/0.2);
		  if (params->persistentRFI==1)
		    {
		      if (freq[ii] > 754 && freq[ii] < 768)
			signal += 1e6;
		      else if (freq[ii] > 768 && freq[ii] < 788)
			signal += 8e5;
		      else if (freq[ii] > 869.95 && freq[ii] < 875.05)
			signal += 4e5;
		      else if (freq[ii] > 875.05 && freq[ii] < 889.95)
			signal += 3e5;
		      else if (freq[ii] > 943.4 && freq[ii] < 951.8)
			signal += 3e5;
		      else if (freq[ii] > 953.7 && freq[ii] < 960)
			signal += 3e5;
		      else if (freq[ii] > 1017 && freq[ii] < 1019)
			signal += 3e5;
		      else if (freq[ii] > 1023 && freq[ii] < 1025)
			signal += 9e5;
		      else if (freq[ii] > 1029 && freq[ii] < 1031)
			signal += 8e5;
		      else if (freq[ii] > 1805 && freq[ii] < 1865)
			signal += 9e5;
		    }
		  data[k*nchan*npol + ii]           = TKgaussDev(&iseed) + ssys_aa + signal;
		  data[k*nchan*npol + nchan + ii]   = TKgaussDev(&iseed) + ssys_bb + signal;
		  data[k*nchan*npol + 2*nchan + ii] = TKgaussDev(&iseed);
		  data[k*nchan*npol + 3*nchan + ii] = TKgaussDev(&iseed);
		}
	    }
	  sdhdf_writeObsParams(outFile,bandHeader[j].label,i,j,obsParams,ndump,1);			       
	  free(obsParams);

	  // SHOULD SET UP ATTRIBUTES
	  sdhdf_writeSpectrumData(outFile,bandHeader[j].label,i,j,data,freq,nchan,npol,ndump,1,dataAttributes,nDataAttributes,freqAttributes,nFreqAttributes);
	  free(freq);
	  free(data);

	}
      sdhdf_writeBandHeader(outFile,bandHeader,i,nband,1);
      free(bandHeader);

    }

  
  sdhdf_closeFile(outFile);
  free(outFile);
  free(primaryHeader);
  free(beamHeader);
  free(params->beam);
  free(params->band);
  free(params->dump);
  free(params);
}

