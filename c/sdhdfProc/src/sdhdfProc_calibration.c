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
#include "fitsio.h"

//
// Routine to load a PCM file (FITS format)
// to model cross-coupling within the receiver and noise source system
//
void sdhdf_loadPCM(sdhdf_calibration *polCal,int *nPolCalChan,char *observatory, char *rcvr,char *pcmFile)
{
  int status=0;
  fitsfile *fptr;
  char fname[1024];
  char runtimeDir[1024];
  int nchan_cal_poln;
  int nchan_pcm;
  
  int colnum_data,colnum_freq;
  float *data;
  double *dataFreq;
  
  float n_fval=0;
  double n_dval=0;
  
  int initflag=0;
  int i,j;

  
  if (getenv("SDHDF_RUNTIME")==0)
    {
      printf("=======================================================================\n");
      printf("Error: sdhdfProc_calibration requires that the SDHDF_RUNTIME directory is set\n");
      printf("=======================================================================\n");
      exit(1);
    }
  strcpy(runtimeDir,getenv("SDHDF_RUNTIME"));
  sprintf(fname,"%s/observatory/%s/calibration/%s/%s",runtimeDir,observatory,rcvr,pcmFile);
  printf("Opening: %s\n",fname);
  
  fits_open_file(&fptr,fname,READONLY,&status);
  fits_report_error(stderr,status);

  // Obtain the Stokes parameters of the noise source
  // NOTE THAT THIS DOESN'T HAVE A REQUENCY AXIS - That's probably a mistake in the FITS file definition
  //
  fits_movnam_hdu(fptr,BINARY_TBL,"CAL_POLN",1,&status);
  fits_read_key(fptr,TINT,"NCHAN",&nchan_cal_poln,NULL,&status);
  printf("Number of channels = %d\n",nchan_cal_poln);
  fits_get_colnum(fptr,CASEINSEN,"DATA",&colnum_data,&status);

  // Each entry has Q/I, U/I and V/I (i.e., 3 parameters)
  data = (float *)malloc(sizeof(float)*nchan_cal_poln*3);
  fits_read_col(fptr,TFLOAT,colnum_data,1,1,nchan_cal_poln*3,&n_fval,data,&initflag,&status);
  for (i=0;i<nchan_cal_poln;i++)
    {
      polCal[i].noiseSource_QoverI = data[i*3];
      polCal[i].noiseSource_UoverI = data[i*3+1];
      polCal[i].noiseSource_VoverI = data[i*3+2];      
    }
  *nPolCalChan = nchan_cal_poln;

  free(data);

  
  // Now load the PCM parameteres
  fits_movnam_hdu(fptr,BINARY_TBL,"FEEDPAR",1,&status);
  fits_report_error(stderr,status);
  
  fits_read_key(fptr,TINT,"NCHAN",&nchan_pcm,NULL,&status);
  if (nchan_pcm != nchan_cal_poln)
    {
      printf("ERROR: In sdhdfProc_calibration: nchan_pcm != nchan_cal_poln\n");
      printf("nchan_pcm = %d\n",nchan_pcm);
      printf("nchan_cal_poln = %d\n",nchan_cal_poln);
      exit(1);      
    }
  //
  // Have 7 parameters to define the feed and noise source system
  //
  data = (float *)malloc(sizeof(float)*nchan_pcm*7);
  fits_get_colnum(fptr,CASEINSEN,"DATA",&colnum_data,&status);

  printf("data column = %d\n",colnum_data);
  fits_report_error(stderr,status);
  // NOTE SOME OF THESE SEEM TO BE NAN - BE CAREFUL READING IN -- FIX ME
  fits_read_col(fptr,TFLOAT,colnum_data,1,1,nchan_pcm*7,&n_fval,data,&initflag,&status);
  fits_report_error(stderr,status);
  for (i=0;i<nchan_cal_poln;i++)
    {
      polCal[i].constant_gain = data[i*7];
      polCal[i].constant_diff_gain = data[i*7+1];
      polCal[i].constant_diff_phase = data[i*7+2];
      polCal[i].constant_b1 = data[i*7+3];
      polCal[i].constant_b2 = data[i*7+4];
      polCal[i].constant_r1 = data[i*7+5];
      polCal[i].constant_r2 = data[i*7+6];	    
    }
  free(data);

  
  // Read the frequency axis
  dataFreq = (double *)malloc(sizeof(double)*nchan_pcm);
  fits_get_colnum(fptr,CASEINSEN,"DAT_FREQ",&colnum_freq,&status);
  fits_read_col(fptr,TDOUBLE,colnum_freq,1,1,nchan_pcm,&n_dval,dataFreq,&initflag,&status);
  fits_report_error(stderr,status);
  for (i=0;i<nchan_cal_poln;i++)
    polCal[i].freq = dataFreq[i];
  
  free(dataFreq);


  

  
  
  
  fits_close_file(fptr,&status);
  
}

int sdhdf_loadTcal(sdhdf_tcal_struct *tcalData,char *fname)
{
  FILE *fin;
  int n=0;
  if (!(fin = fopen(fname,"r")))
    {
      printf("Unable to open filename >%s<\n",fname);
      exit(1);
    }
  while (!feof(fin))
    {
      if (fscanf(fin,"%lf %lf %lf",&(tcalData[n].freq),&(tcalData[n].tcalA),&(tcalData[n].tcalB))==3)
	n++;
    }
  fclose(fin);
  return n;
}

void sdhdf_get_tcal(sdhdf_tcal_struct *tcalData,int n,double f0,double *tcalA,double *tcalB)
{
  int i;
  double m,c;
  for (i=0;i<n-1;i++)
    {
      // Should do an interpolation here
      if (f0 >= tcalData[i].freq && f0 < tcalData[i+1].freq)
	{
	  m = (tcalData[i].tcalA-tcalData[i+1].tcalA)/(tcalData[i].freq-tcalData[i+1].freq);
	  c = tcalData[i].tcalA-m*tcalData[i].freq;
	  *tcalA = m*f0+c;

	  m = (tcalData[i].tcalB-tcalData[i+1].tcalB)/(tcalData[i].freq-tcalData[i+1].freq);
	  c = tcalData[i].tcalB-m*tcalData[i].freq;
	  *tcalB = m*f0+c;

	  break;
	}
    }
}

//
// This agrees with compute_stokes in simpol.C in PSRCHIVE
//
void sdhdf_convertStokes(float p1,float p2,float p3,float p4,float *stokesI,float *stokesQ,float *stokesU,float *stokesV)
{
  *stokesI = p1 + p2;
  *stokesQ = p1 - p2;
  *stokesU = 2*p3;
  *stokesV = 2*p4;
}
