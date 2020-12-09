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
#include <string.h>
#include <math.h>
#include "fitsio.h"
#include "sdhdfProc.h"
#include "hdf5.h"

#define VNUM "v0.1"

int main(int argc,char *argv[])
{
  int i,j,k,b;

  int fitsType=1; // 1 = SDFITS, 2 = PSRFITS
  
  sdhdf_fileStruct *outFile;
  sdhdf_softwareVersionsStruct *softwareVersions;
  sdhdf_historyStruct *history;
  sdhdf_primaryHeaderStruct *primaryHeader;
  sdhdf_beamHeaderStruct *beamHeader;
  sdhdf_bandHeaderStruct *bandHeader;
  sdhdf_obsParamsStruct  *obsParams;
  sdhdf_obsParamsStruct  *calObsParams;
  sdhdf_bandHeaderStruct *calBandHeader;
  
  
  // Input FITS file
  int status=0;
  fitsfile *fptr;
  int maxdim = 5;
  int naxis;
  long naxes[maxdim];
  int type;
  long nrows;
  long maxsize;
  int    nchan;
  int    nchan_cal;
  int    npol;
  int    npol_cal;
  int    ndump;
  int    ndump_cal;
  long   long_ndump;
  double freq0;
  double chan_bw;
  float  *floatVals;
  float  *floatCalValsOn;
  float  *floatCalValsOff;
  float  *infloatVals;
  float  *freqVals;
  float  *freqCalVals;
  float *datOffs;
  float *datScl;
  double *doubleVals;
  short int *shortVals;
  short int n_sval=0;
  int    binVal = 0; 
  int    nVals;
  float  n_fval=0.;
  double n_dval=0.;
  long   n_lval=0;
  int    colnum;
  int    colnum_datoffs;
  int    colnum_datscl;

  int    initflag=0;
  double mean=0;

  char fname[1024];
  char outname[1024]="convert.hdf";

  int cal=0;
  char calFile[1024];
  int  calOn1,calOn2;
  int  calOff1,calOff2;

  int pulseOff1,pulseOff2;
  double *baseline1,*baseline2,*baseline3,*baseline4;
  pulseOff1=pulseOff2=-1;
  
  primaryHeader    = (sdhdf_primaryHeaderStruct *)malloc(sizeof(sdhdf_primaryHeaderStruct));
  beamHeader       = (sdhdf_beamHeaderStruct *)malloc(sizeof(sdhdf_beamHeaderStruct));
  bandHeader       = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct));
  softwareVersions = (sdhdf_softwareVersionsStruct *)malloc(sizeof(sdhdf_softwareVersionsStruct));
  history          = (sdhdf_historyStruct *)malloc(sizeof(sdhdf_historyStruct));
  
  sdhdf_setMetadataDefaults(primaryHeader,beamHeader,bandHeader,softwareVersions,history,1,1);

  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	strcpy(fname,argv[++i]);
      else if (strcmp(argv[i],"-cal")==0)
	{
	  cal=1;
	  strcpy(calFile,argv[++i]);
	}
      else if (strcmp(argv[i],"-calOn")==0)
	{
	  sscanf(argv[++i],"%d",&calOn1);
	  sscanf(argv[++i],"%d",&calOn2);
	}
      else if (strcmp(argv[i],"-calOff")==0)
	{
	  sscanf(argv[++i],"%d",&calOff1);
	  sscanf(argv[++i],"%d",&calOff2);
	}
      else if (strcmp(argv[i],"-pulseOff")==0)
	{
	  sscanf(argv[++i],"%d",&pulseOff1);
	  sscanf(argv[++i],"%d",&pulseOff2);
	}
      else if (strcmp(argv[i],"-bin")==0)
	sscanf(argv[++i],"%d",&binVal);
      else if (strcmp(argv[i],"-type")==0)
	sscanf(argv[++i],"%d",&fitsType);
      else if (strcmp(argv[i],"-o")==0)
	strcpy(outname,argv[++i]);
    }
  
  

  if (fitsType==1) // SDFITS
    {    
      printf("Opening file >%s<\n",fname);
      fits_open_file(&fptr,fname,READONLY,&status);
      fits_report_error(stderr,status);
      fits_movnam_hdu(fptr,BINARY_TBL,"SINGLE DISH",1,&status);
      
      fits_get_num_rows(fptr,&long_ndump,&status);
      ndump = (int)long_ndump;
      printf("ndump = %d\n",ndump);
      doubleVals = (double *)malloc(sizeof(double)*ndump);
      
      fits_get_colnum(fptr,CASEINSEN,"DATA",&colnum,&status);
      fits_read_tdim(fptr,colnum,maxdim,&naxis,naxes,&status);
      
      fits_get_colnum(fptr,CASEINSEN,"NCHAN",&colnum,&status);
      npol = naxes[0];
      nchan = naxes[1];
      printf("nchan = %d, npol = %d, status = %d\n",nchan,npol,status);
      
      infloatVals = (float *)malloc(sizeof(float)*nchan*npol);
      floatVals = (float *)malloc(sizeof(float)*nchan*npol*ndump);
      freqVals = (float *)malloc(sizeof(float)*nchan);
      
      fits_get_colnum(fptr,CASEINSEN,"FREQ",&colnum,&status);
      fits_read_col(fptr,TDOUBLE,colnum,1,1,ndump,&n_dval,doubleVals,&initflag,&status);
      // Note that the frequency should change based on the LSR correction for different sub-ints -- but it doesn't seem to do this
      freq0 = doubleVals[0];
      printf("freq0 = %g status= %d\n",freq0,status);
      
      fits_get_colnum(fptr,CASEINSEN,"CHAN_BW",&colnum,&status);
      fits_read_col(fptr,TDOUBLE,colnum,1,1,ndump,&n_dval,doubleVals,&initflag,&status);
      // Note that the frequency should change based on the LSR correction for different sub-ints -- but it doesn't seem to do this
      chan_bw = doubleVals[0];

      printf("Freq0 = %g chan_bw = %g\n",freq0,chan_bw);
      
      fits_get_colnum(fptr,CASEINSEN,"DATA",&colnum,&status);

      for (i=0;i<nchan;i++)
	freqVals[i] = freq0+i*chan_bw;

      for (i=0;i<ndump;i++)
	{
	  printf("Loading spectrum dump %d/%d\n",i,ndump-1);
	  fits_read_col(fptr,TFLOAT,colnum,i+1,1,nchan*npol,&n_fval,infloatVals,&initflag,&status);
	  // Rearrange ordering
	  for (j=0;j<npol;j++)
	    {
	      for (k=0;k<nchan;k++)
		{
		  floatVals[i*nchan*npol + j*nchan + k] = infloatVals[k*npol+j];
		}
	    }
	}
    }
  else if (fitsType==2)
    {
      int nbin;
      int nbin_cal;
      int nOn,nOff;
      if (cal==1)
	{
	  calBandHeader  = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct));
	  calObsParams   = (sdhdf_obsParamsStruct *)malloc(sizeof(sdhdf_obsParamsStruct)*ndump);

	  strcpy(primaryHeader->cal_mode,"ON");
	  
	  
	  // First process cal
	  printf("Opening file >%s<\n",calFile);
	  fits_open_file(&fptr,calFile,READONLY,&status);
	  fits_report_error(stderr,status);      

	  printf("Reading PSRFITS file\n");
	  fits_movnam_hdu(fptr,BINARY_TBL,"SUBINT",1,&status);
	  fits_get_num_rows(fptr,&long_ndump,&status);
	  ndump_cal = (int)long_ndump;
	  
	  doubleVals = (double *)malloc(sizeof(double)*ndump_cal);
	  
	  fits_get_colnum(fptr,CASEINSEN,"DATA",&colnum,&status);
	  fits_read_tdim(fptr,colnum,maxdim,&naxis,naxes,&status);
	  fits_report_error(stderr,status);
	  
	  nbin_cal  = naxes[0];
	  nchan_cal = naxes[1];
	  npol_cal  = naxes[2];
	  printf("nbin = %d, nchan = %d, npol = %d, status = %d\n",nbin_cal,nchan_cal,npol_cal,status);
	  
	  floatCalValsOn  = (float *)malloc(sizeof(float)*nchan_cal*npol_cal*ndump_cal);
	  floatCalValsOff = (float *)malloc(sizeof(float)*nchan_cal*npol_cal*ndump_cal);
	  freqCalVals    = (float *)malloc(sizeof(float)*nchan_cal);
	  shortVals   = (short int *)malloc(sizeof(short int)*nchan_cal*npol_cal*nbin_cal);
	  
	  fits_get_colnum(fptr,CASEINSEN,"DAT_OFFS",&colnum_datoffs,&status);
	  fits_get_colnum(fptr,CASEINSEN,"DAT_SCL",&colnum_datscl,&status);
	  datOffs = (float *)malloc(sizeof(float)*nchan_cal*npol_cal);
	  datScl = (float *)malloc(sizeof(float)*nchan_cal*npol_cal);	

	  
	  fits_get_colnum(fptr,CASEINSEN,"DAT_FREQ",&colnum,&status);
	  fits_read_col(fptr,TFLOAT,colnum,1,1,nchan_cal,&n_fval,freqCalVals,&initflag,&status);


	  strcpy(calBandHeader->label,"band0");
	  calBandHeader->fc = (freqCalVals[0]+freqCalVals[nchan_cal-1])/2.0; 
	  calBandHeader->f0 = freqCalVals[0];
	  calBandHeader->f1 = freqCalVals[nchan_cal-1]; 
	  calBandHeader->nchan = nchan_cal;
	  calBandHeader->npol = npol_cal;
	  strcpy(calBandHeader->pol_type,"AABBCRCI");
	  calBandHeader->dtime = 100; // FIX
	  calBandHeader->ndump = ndump_cal;

	  
	  fits_get_colnum(fptr,CASEINSEN,"DATA",&colnum,&status);
	  for (i=0;i<ndump_cal;i++)
	    {
	      fits_read_col(fptr,TFLOAT,colnum_datoffs,i+1,1,nchan_cal*npol_cal,&n_fval,datOffs,&initflag,&status);
	      fits_read_col(fptr,TFLOAT,colnum_datscl,i+1,1,nchan_cal*npol_cal,&n_fval,datScl,&initflag,&status);


	      calObsParams[i].timeElapsed = i; // FIX
	      strcpy(calObsParams[i].timedb,"UNKNOWN"); // FIX
	      calObsParams[i].mjd = 56000; // FIX
	      strcpy(calObsParams[i].utc,"UNKNOWN"); // FIX
	      strcpy(calObsParams[i].ut_date,"UNKNOWN"); // FIX
	      strcpy(calObsParams[i].aest,"UNKNOWN"); // FIX  --- AND NOT AEST -- HAVE LOCAL TIME
	      strcpy(calObsParams[i].raStr,"UNKNOWN");
	      strcpy(calObsParams[i].decStr,"UNKNOWN");
	      calObsParams[i].raDeg = 0;
	      calObsParams[i].decDeg = 0;
	      calObsParams[i].raOffset = 0;
	      calObsParams[i].decOffset = 0;
	      calObsParams[i].gl = 0;
	      calObsParams[i].gb = 0;
	      calObsParams[i].az = 0;
	      calObsParams[i].ze = 0;
	      calObsParams[i].el = 0;
	      calObsParams[i].az_drive_rate = 0;
	      calObsParams[i].ze_drive_rate = 0;
	      calObsParams[i].hourAngle = 0;
	      calObsParams[i].paraAngle = 0;
	      calObsParams[i].windDir = 0;
	      calObsParams[i].windSpd = 0;      


	      
	      printf("Loading subintegration %d/%d\n",i,ndump_cal-1);
	      fits_read_col(fptr,TSHORT,colnum,i+1,1,nchan_cal*npol_cal*nbin_cal,&n_sval,shortVals,&initflag,&status);
	      // Select bin and scale by DAT_OFFS and DAT_SCL
	      for (j=0;j<npol_cal;j++)
		{
		  for (k=0;k<nchan_cal;k++)
		    {
		      floatCalValsOn[i*nchan_cal*npol_cal + j*nchan_cal + k] = 0;
		      floatCalValsOff[i*nchan_cal*npol_cal + j*nchan_cal + k] = 0;
		      nOn = nOff = 0;
		      for (b=0;b<nbin_cal;b++)
			{
			  if (b >= calOn1 && b <= calOn2)
			    {			    
			      floatCalValsOn[i*nchan_cal*npol_cal + j*nchan_cal + k] += (shortVals[k*nbin_cal+j*nchan_cal*nbin_cal + b])*datScl[j*nchan_cal+k]+datOffs[j*nchan_cal+k]; // WHAT ABOUT zeroOff?? -- CHECK
			      nOn++;
			    }
			  if (b>= calOff1 && b <= calOff2)
			    {
			      floatCalValsOff[i*nchan_cal*npol_cal + j*nchan_cal + k] += (shortVals[k*nbin_cal+j*nchan_cal*nbin_cal + b])*datScl[j*nchan_cal+k]+datOffs[j*nchan_cal+k]; // WHAT ABOUT zeroOff?? -- CHECK
			      nOff++;
			    }
			}
		      //		      printf("%d %d nOn = %d nOff = %d %d %d %d %d %g %g\n",j,k,nOn,nOff,calOn1,calOn2,calOff1,calOff2, floatCalValsOn[i*nchan_cal*npol_cal + j*nchan_cal + k],floatCalValsOff[i*nchan_cal*npol_cal + j*nchan_cal + k]);
		      floatCalValsOn[i*nchan_cal*npol_cal + j*nchan_cal + k] /= (double)nOn;
		      floatCalValsOff[i*nchan_cal*npol_cal + j*nchan_cal + k] /= (double)nOff;
		    }
		}
	    }
	
	  printf("Closing cal file\n");
	  fits_close_file(fptr,&status);
	  free(shortVals); free(datOffs); free(datScl);
	}
      // Now process the astronomy file
      
      printf("Opening file >%s<\n",fname);
      fits_open_file(&fptr,fname,READONLY,&status);
      fits_report_error(stderr,status);
      
      printf("Reading PSRFITS file\n");
      fits_movnam_hdu(fptr,BINARY_TBL,"SUBINT",1,&status);
      fits_get_num_rows(fptr,&long_ndump,&status);
      ndump = (int)long_ndump;

      doubleVals = (double *)malloc(sizeof(double)*ndump);
     
      fits_get_colnum(fptr,CASEINSEN,"DATA",&colnum,&status);
      fits_read_tdim(fptr,colnum,maxdim,&naxis,naxes,&status);
      fits_report_error(stderr,status);

      nbin  = naxes[0];
      nchan = naxes[1];
      npol  = naxes[2];
      printf("nbin = %d, nchan = %d, npol = %d, status = %d\n",nbin,nchan,npol,status);
      
      infloatVals = (float *)malloc(sizeof(float)*nchan*npol); // DON'T NEED THIS ONE .. BUT FREE AT THE END .. FIX ME
      floatVals   = (float *)malloc(sizeof(float)*nchan*npol*ndump*nbin);
      freqVals    = (float *)malloc(sizeof(float)*nchan);
      shortVals   = (short int *)malloc(sizeof(short int)*nchan*npol*nbin);
      baseline1   = (double *)malloc(sizeof(double)*nchan);
      baseline2   = (double *)malloc(sizeof(double)*nchan);
      baseline3   = (double *)malloc(sizeof(double)*nchan);
      baseline4   = (double *)malloc(sizeof(double)*nchan);
      
      
      fits_get_colnum(fptr,CASEINSEN,"DAT_OFFS",&colnum_datoffs,&status);
      fits_get_colnum(fptr,CASEINSEN,"DAT_SCL",&colnum_datscl,&status);
      datOffs = (float *)malloc(sizeof(float)*nchan*npol);
      datScl = (float *)malloc(sizeof(float)*nchan*npol);	
      
      fits_get_colnum(fptr,CASEINSEN,"DAT_FREQ",&colnum,&status);
      fits_read_col(fptr,TFLOAT,colnum,1,1,nchan,&n_fval,freqVals,&initflag,&status);

      fits_get_colnum(fptr,CASEINSEN,"DATA",&colnum,&status);

      // Do we want a baseline
      if (pulseOff1 >= 0 && pulseOff2 >= 0)
	{
	  int kk;
	  for (k=0;k<nchan;k++)
	    baseline1[k]=baseline2[k]=baseline3[k]=baseline4[k]=0;

	  for (i=0;i<ndump;i++)
	    {
	      fits_read_col(fptr,TFLOAT,colnum_datoffs,i+1,1,nchan*npol,&n_fval,datOffs,&initflag,&status);
	      fits_read_col(fptr,TFLOAT,colnum_datscl,i+1,1,nchan*npol,&n_fval,datScl,&initflag,&status);
	      
	      printf("Loading subintegration %d/%d\n",i,ndump-1);
	      fits_read_col(fptr,TSHORT,colnum,i+1,1,nchan*npol*nbin,&n_sval,shortVals,&initflag,&status);
	      // Select bin and scale by DAT_OFFS and DAT_SCL
	      for (k=0;k<nchan;k++)
		{
		  for (kk=pulseOff1;kk<pulseOff2;kk++)
		    {
		      j=0; baseline1[k] += (shortVals[k*nbin+j*nchan*nbin + kk])*datScl[j*nchan+k]+datOffs[j*nchan+k]; // WHAT ABOUT zeroOff?? -- CHECK
		      j=1; baseline2[k] += (shortVals[k*nbin+j*nchan*nbin + kk])*datScl[j*nchan+k]+datOffs[j*nchan+k]; // WHAT ABOUT zeroOff?? -- CHECK
		      j=2; baseline3[k] += (shortVals[k*nbin+j*nchan*nbin + kk])*datScl[j*nchan+k]+datOffs[j*nchan+k]; // WHAT ABOUT zeroOff?? -- CHECK
		      j=3; baseline4[k] += (shortVals[k*nbin+j*nchan*nbin + kk])*datScl[j*nchan+k]+datOffs[j*nchan+k]; // WHAT ABOUT zeroOff?? -- CHECK
		    }
		}
	    }

	  for (k=0;k<nchan;k++)
	    {
	      baseline1[k] /= (double)(ndump*(pulseOff2-pulseOff1));
	      baseline2[k] /= (double)(ndump*(pulseOff2-pulseOff1));
	      baseline3[k] /= (double)(ndump*(pulseOff2-pulseOff1));
	      baseline4[k] /= (double)(ndump*(pulseOff2-pulseOff1));
	      printf("baseline %d %g %g %g %g\n",k,baseline1[k],baseline2[k],baseline3[k],baseline4[k]);
	    }
	}
      else
	{
	  for (k=0;k<nchan;k++)
	    baseline1[k]=baseline2[k]=baseline3[k]=baseline4[k]=0;
	}

      for (i=0;i<ndump;i++)
	{
	  fits_read_col(fptr,TFLOAT,colnum_datoffs,i+1,1,nchan*npol,&n_fval,datOffs,&initflag,&status);
	  fits_read_col(fptr,TFLOAT,colnum_datscl,i+1,1,nchan*npol,&n_fval,datScl,&initflag,&status);

	  printf("Loading subintegration %d/%d\n",i,ndump-1);
	  fits_read_col(fptr,TSHORT,colnum,i+1,1,nchan*npol*nbin,&n_sval,shortVals,&initflag,&status);
	  // Select bin and scale by DAT_OFFS and DAT_SCL
	  for (j=0;j<npol;j++)
	    {
	      for (k=0;k<nchan;k++)
		{
		  floatVals[i*nchan*npol + j*nchan + k] = (shortVals[k*nbin+j*nchan*nbin + binVal])*datScl[j*nchan+k]+datOffs[j*nchan+k]; // WHAT ABOUT zeroOff?? -- CHECK
		  if (j==0)      floatVals[i*nchan*npol + j*nchan + k] -= baseline1[k];
		  else if (j==1) floatVals[i*nchan*npol + j*nchan + k] -= baseline2[k];
		  else if (j==2) floatVals[i*nchan*npol + j*nchan + k] -= baseline3[k];
		  else if (j==3) floatVals[i*nchan*npol + j*nchan + k] -= baseline4[k];
		}
	    }
	}
      free(shortVals);
    }
  printf("Closing file\n");
  fits_close_file(fptr,&status);
  obsParams = (sdhdf_obsParamsStruct *)malloc(sizeof(sdhdf_obsParamsStruct)*ndump);
  for (i=0;i<ndump;i++)
    {
      obsParams[i].timeElapsed = i; // FIX
      strcpy(obsParams[i].timedb,"UNKNOWN"); // FIX
      obsParams[i].mjd = 56000; // FIX
      strcpy(obsParams[i].utc,"UNKNOWN"); // FIX
      strcpy(obsParams[i].ut_date,"UNKNOWN"); // FIX
      strcpy(obsParams[i].aest,"UNKNOWN"); // FIX  --- AND NOT AEST -- HAVE LOCAL TIME
      strcpy(obsParams[i].raStr,"UNKNOWN");
      strcpy(obsParams[i].decStr,"UNKNOWN");
      obsParams[i].raDeg = 0;
      obsParams[i].decDeg = 0;
      obsParams[i].raOffset = 0;
      obsParams[i].decOffset = 0;
      obsParams[i].gl = 0;
      obsParams[i].gb = 0;
      obsParams[i].az = 0;
      obsParams[i].ze = 0;
      obsParams[i].el = 0;
      obsParams[i].az_drive_rate = 0;
      obsParams[i].ze_drive_rate = 0;
      obsParams[i].hourAngle = 0;
      obsParams[i].paraAngle = 0;
      obsParams[i].windDir = 0;
      obsParams[i].windSpd = 0;      
    }

  // Open the output HDF5 file
  if (!(outFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >outFile<\n");
      exit(1);
    }
  sdhdf_initialiseFile(outFile);

  // Set up the primary header information
      
  // Set up the beam information

  // Set up the band information
  strcpy(bandHeader->label,"band0");
  bandHeader->fc = (freqVals[0]+freqVals[nchan-1])/2.0; 
  bandHeader->f0 = freqVals[0];
  bandHeader->f1 = freqVals[nchan-1]; 
  bandHeader->nchan = nchan;
  bandHeader->npol = npol;
  strcpy(bandHeader->pol_type,"AABBCRCI");
  bandHeader->dtime = 100; // FIX
  bandHeader->ndump = ndump;

  printf("Writing output file\n");
  sdhdf_openFile(outname,outFile,3);
  sdhdf_writePrimaryHeader(outFile,primaryHeader);
  sdhdf_writeBeamHeader(outFile,beamHeader,1);
  sdhdf_writeBandHeader(outFile,bandHeader,0,1,1);
  sdhdf_writeObsParams(outFile,bandHeader[0].label,0,0,obsParams,ndump,1);
  sdhdf_writeSpectrumData(outFile,bandHeader->label,0,0,floatVals,freqVals,nchan,npol,ndump,1);
  if (cal==1)
    {
      sdhdf_writeBandHeader(outFile,calBandHeader,0,1,2);
      sdhdf_writeObsParams(outFile,bandHeader[0].label,0,0,obsParams,ndump_cal,2);
      printf("In here\n");
      // GEORGE HERE:
      // ... need to setup and then write the cal metadata
      // 

      
      sdhdf_writeSpectrumData(outFile,bandHeader->label,0,0,floatCalValsOn,freqCalVals,nchan_cal,npol_cal,ndump_cal,2);
      sdhdf_writeSpectrumData(outFile,bandHeader->label,0,0,floatCalValsOff,freqCalVals,nchan_cal,npol_cal,ndump_cal,3);
      printf("Completed writing cal\n");      
    }
  sdhdf_writeSoftwareVersions(outFile,softwareVersions);
  sdhdf_writeHistory(outFile,history,1);
  
  sdhdf_closeFile(outFile);
  free(outFile);

  if (cal==1)
    {
      free(floatCalValsOn);
      free(floatCalValsOff);
      free(freqCalVals);
      free(calBandHeader);
      free(calObsParams);	    
    }
if (type==2)
  {
    free(baseline1); free(baseline2); free(baseline3); free(baseline4);
  }
  free(primaryHeader);
  free(beamHeader);
  free(bandHeader);
  free(obsParams);
  free(softwareVersions);
  free(history);
  free(floatVals);
  free(infloatVals);
  free(freqVals);
  free(doubleVals);
} 