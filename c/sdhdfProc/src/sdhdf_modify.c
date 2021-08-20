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
// Software to produce a new SDHDF file based processing an existing one
//
// Usage:
// sdhdf_modify <filename.hdf> -o <outputFile.hdf>
//
// Compilation
// gcc -lm -o sdhdf_modify sdhdf_modify.c sdhdfProc.c -I../hdf5-1.10.4/src/ ../hdf5-1.10.4/src/.libs/libhdf5.a -ldl -lz -L/pulsar/psr/software/stable/stretch/lib/ -I//pulsar/psr/software/stable/stretch/include -lcalceph -Isofa/20190722/c/src/ -Lsofa/20190722/c/src -lsofa_c
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sdhdfProc.h"
#include "hdf5.h"
#include "TKfit.h"
#include <libgen.h>

#define MAX_IGNORE_SDUMP 4096

double linearInterpolate(double *freq,double *scaleAA,double in_freq,int kk,int nScale);
void   getScal(double sysGain_freq,double *scalFreq,double *scalAA,double *scalBB,int nScal,double *scalAA_val,double *scalBB_val);
void   TKspline_interpolate(int n,double *x,double *y,double yd[][4],double *interpX,
			  double *interpY,int nInterp);
void   TKcmonot (int n, double x[], double y[], double yd[][4]);


// SHOULD ALLOCATE THIS PROPERLY **** FIX ME
double yd1[4096][4];
double yd2[4096][4];

int main(int argc,char *argv[])
{
  int i,ii,j,k,p,l,kp,b;
  int nchan,npol;
  int interp=1; // 1 = linear, 2 = cubic spline
  char fname[MAX_STRLEN];
  char bname[MAX_STRLEN];
  char oname[MAX_STRLEN];
  sdhdf_fileStruct *inFile,*outFile;
  herr_t status;
  sdhdf_bandHeaderStruct *inBandParams;
  sdhdf_obsParamsStruct  *outObsParams;

  sdhdf_attributes_struct dataAttributes[MAX_ATTRIBUTES];
  sdhdf_attributes_struct freqAttributes[MAX_ATTRIBUTES];
  int nDataAttributes=0;
  int nFreqAttributes=0;
  int polAv=1;

  int nBeam=0;
  int nBand=0;
  int tScrunch=0;
  int fScrunch=0;
  int pol=0;
  int stokes=0;
  int fAv=0;
  float fchbw=-1;
  int nchanRequested=-1;
  int bary=0,lsr=0;
  int regrid=0;
  int nh;
  float *in_data,*out_data,*out_freq;
  double *in_freq;
  float *out_Tdata,*out_Fdata;
  float *dataWts;
  
  long out_nchan,out_npol,out_ndump,nsd;
  double av1,av2,av3,av4,avFreq;
  double tav=0;
  char extension[MAX_STRLEN];
  sdhdf_eopStruct *eop;
  int nEOP=0;
  char args[MAX_STRLEN]="";
  char scaleFile[1024]="NULL";
  int ignoreSD[MAX_IGNORE_SDUMP];
  int nIgnoreSD=0;
  double freq[MAX_CHAN_SCAL],scaleAA[MAX_CHAN_SCAL],scaleBB[MAX_CHAN_SCAL];

  int nScale=0;
  double *scalFreq,*scalAA,*scalBB;
  float *cal_p1,*cal_p2,*cal_freq;

  double sysGain_freq[4096],sysGain_p1[4096],sysGain_p2[4096],scalAA_val,scalBB_val;
  double lastGood_p1=1e8,lastGood_p2=1e8;
  char   scalFname[1024]="NULL";
  int    nScal=-1;
  int    c_nchan;
  int    nc,nc2;
  int    nRF1,nRF2,nRF3;
  int    calBin = 32;  // *** MUST SET PROPERLY **

  double scA_rf1_val,scB_rf1_val;
  double scA_rf2_val,scB_rf2_val;
  double scA_rf3_val,scB_rf3_val;
  double scA_rf1_x[4096],scA_rf1_y[4096];
  double scB_rf1_x[4096],scB_rf1_y[4096];
  double scA_rf2_x[4096],scA_rf2_y[4096];
  double scB_rf2_x[4096],scB_rf2_y[4096];
  double scA_rf3_x[4096],scA_rf3_y[4096];
  double scB_rf3_x[4096],scB_rf3_y[4096];
  double scA_rf1[3],scB_rf1[3]; // Parameters modelling the bands			  
  double scA_rf2[3],scB_rf2[3];
  double scA_rf3[3],scB_rf3[3];

  double diffGain2[4096],diffPhase[4096],g_e[4096];
  double scale1[4096],scale2[4096],scale3[4096],scale4[4096];
  double useScale1,useScale2,useScale3,useScale4;
  double t_phase,t_len;
  double calAA,calBB,calReAB,calImAB;
  double usePara=0;
  int    setPara=0;
  int    tcal=0;
  
  double divideCal=1;
  double divideAstro=1;

  double cc_a=0,cc_b=0,cc_c=0,cc_d=0; // Cross coupling parameters
  double cc_phi1,cc_phi2;
  double cc_eps1,cc_eps2;
  int muellerI = 0;
  int origNchan = -1;
  double chbw;
  float tdumpRequest = -1;

  char ephemName[1024] = "DE436.1950.2050";

  
  cc_eps1=cc_eps2=0;
  
  scA_rf1[0]=scA_rf1[1]=scA_rf1[2] = 0.0;
  scB_rf1[0]=scB_rf1[1]=scB_rf1[2] = 0.0;
  scA_rf2[0]=scA_rf2[1]=scA_rf2[2] = 0.0;
  scB_rf2[0]=scB_rf2[1]=scB_rf2[2] = 0.0;
  scA_rf3[0]=scA_rf3[1]=scA_rf3[2] = 0.0;
  scB_rf3[0]=scB_rf3[1]=scB_rf3[2] = 0.0;

  
  int calFitN=0;
  
  strcpy(oname,"sdhdf_modify_output.hdf");
  strcpy(extension,"modify");
  
  eop = (sdhdf_eopStruct *)malloc(sizeof(sdhdf_eopStruct)*MAX_EOP_LINES);
  
  scalFreq = (double *)malloc(sizeof(double)*MAX_CHAN_SCAL);
  scalAA = (double *)malloc(sizeof(double)*MAX_CHAN_SCAL);
  scalBB = (double *)malloc(sizeof(double)*MAX_CHAN_SCAL);
  
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
	{sdhdf_add2arg(args,argv[i],argv[i+1]); strcpy(extension,argv[++i]);}
      else if (strcmp(argv[i],"-ephem")==0)
	strcpy(ephemName,argv[++i]);
      else if (strcmp(argv[i],"-origNchan")==0)
	sscanf(argv[++i],"%d",&origNchan);
      else if (strcmp(argv[i],"-divideCal")==0)
	sscanf(argv[++i],"%lf",&divideCal);
      else if (strcmp(argv[i],"-divideAstro")==0)
	sscanf(argv[++i],"%lf",&divideAstro);
      else if (strcmp(argv[i],"-calFitN")==0) // Currently 1, 2 or 3 (1 = default) for fitting to the ON-OFF cal values
	  sscanf(argv[++i],"%d",&calFitN);
      else if (strcmp(argv[i],"-calBin")==0)
	sscanf(argv[++i],"%d",&calBin);
      else if (strcmp(argv[i],"-muellerI")==0)
	muellerI=1;
      else if (strcasecmp(argv[i],"-cc_phi1")==0)
	sscanf(argv[++i],"%lf",&cc_phi1);
      else if (strcasecmp(argv[i],"-cc_phi2")==0)
	sscanf(argv[++i],"%lf",&cc_phi2);
      else if (strcasecmp(argv[i],"-cc_eps1")==0)
	sscanf(argv[++i],"%lf",&cc_eps1);
      else if (strcasecmp(argv[i],"-cc_eps2")==0)
	sscanf(argv[++i],"%lf",&cc_eps2);
      else if (strcmp(argv[i],"-scal")==0)
	{strcpy(scalFname,argv[++i]); tcal=0;}
      else if (strcmp(argv[i],"-tcal")==0)
	{strcpy(scalFname,argv[++i]); tcal=1;}
      else if (strcmp(argv[i],"-psum")==0)
	polAv=0;
      else if (strcmp(argv[i],"-para")==0)
	{
	  setPara=1;
	  sscanf(argv[++i],"%lf",&usePara);
	}
      else if (strcmp(argv[i],"-interp")==0) // 1 = linear, 2 = cubic spline
	sscanf(argv[++i],"%d",&interp);
      else if (strcmp(argv[i],"-regrid")==0)
	regrid=1;
      else if (strcmp(argv[i],"-T")==0)
	{sdhdf_add1arg(args,argv[i]); tScrunch=1;}
      else if (strcmp(argv[i],"-F")==0)
	{sdhdf_add1arg(args,argv[i]); fScrunch=1;}
      else if (strcmp(argv[i],"-tdump")==0)
	{sdhdf_add1arg(args,argv[i]); sscanf(argv[++i],"%f",&tdumpRequest);}
      else if (strcmp(argv[i],"-scaleFile")==0)
	{sdhdf_add2arg(args,argv[i],argv[i+1]); strcpy(scaleFile,argv[++i]);}
      else if (strcmp(argv[i],"-id")==0)
	{sdhdf_add2arg(args,argv[i],argv[i+1]); sscanf(argv[++i],"%d",&ignoreSD[nIgnoreSD++]);}
      else if (strcmp(argv[i],"-p1")==0)
	{sdhdf_add1arg(args,argv[i]); pol=1;}
      else if (strcmp(argv[i],"-p2")==0)
	{sdhdf_add1arg(args,argv[i]); pol=2;}
      else if (strcmp(argv[i],"-p3")==0)
	{sdhdf_add1arg(args,argv[i]); pol=3;}
      else if (strcmp(argv[i],"-S")==0)
	{sdhdf_add1arg(args,argv[i]); stokes=1;}
      else if (strcmp(argv[i],"-bary")==0)
	{
	  bary=1;
	  nEOP = sdhdf_loadEOP(eop);
	  sdhdf_add1arg(args,argv[i]); 
	}
      else if (strcmp(argv[i],"-lsr")==0)
	{
	  lsr=1;
	  nEOP = sdhdf_loadEOP(eop);
	  sdhdf_add1arg(args,argv[i]); 
	}
      else if (strcasecmp(argv[i],"-fav")==0)
	{sdhdf_add2arg(args,argv[i],argv[i+1]); sscanf(argv[++i],"%d",&fAv);}
      else if (strcasecmp(argv[i],"-fchbw")==0)
	{sdhdf_add2arg(args,argv[i],argv[i+1]); sscanf(argv[++i],"%f",&fchbw);}
      else if (strcasecmp(argv[i],"-nch")==0)
	{sdhdf_add2arg(args,argv[i],argv[i+1]); sscanf(argv[++i],"%d",&nchanRequested);}
      else
	{
	  strcpy(fname,argv[i]);
	  printf("Processing %s\n",fname);
	  //	  strcpy(bname,basename(fname));
	  //	  printf("Bname = %s\n",bname);
	  if (strcmp(scaleFile,"NULL")!=0)
	    {
	      FILE *fScale;
	      nScale=0;
	      fScale = fopen(scaleFile,"r");
	      while (!feof(fScale))
		{
		  if (fscanf(fScale,"%lf %lf %lf",&freq[nScale],&scaleAA[nScale],&scaleBB[nScale])==3)
		    {
		      nScale++;
		      if (nScale == MAX_CHAN_SCAL+1)
			{
			  printf("ERROR: Must increase MAX_CHAN_SCAL\n");
			  fclose(fScale);
			  exit(1); // Should exit more gracefully
			}
		    }
		}
	      fclose(fScale);
	    }
	
	  if (strcmp(scalFname,"NULL")!=0)
	    {
	      FILE *fin;

	      if (strcmp(scalFname,"auto")==0)
		{
		  char runtimeDir[1024];
		  if (getenv("SDHDF_RUNTIME")==0)
		    {
		      printf("=======================================================================\n");
		      printf("Error: sdhdf_modify requires that the SDHDF_RUNTIME directory is set\n");
		      printf("=======================================================================\n");
		      exit(1);
		    }
		  strcpy(runtimeDir,getenv("SDHDF_RUNTIME"));
		  sprintf(scalFname,"%s/observatory/parkes/calibration/scalAverage.dat",runtimeDir); // FIX ME -- TOO MUCH HARDCODED HERE
		  printf("Using SCAL = %s\n",scalFname);
		}
	      
	      nScal=0;
	      if ((fin = fopen(scalFname,"r"))==NULL)
		{
		  printf("Unable to open file >%s<\n",scalFname);
		  exit(1); // Should close files before exiting ::: FIX ME
		}
	      printf("Loading %s\n",scalFname);
	      while (!feof(fin))
		{
		  if (fscanf(fin,"%lf %lf %lf",&scalFreq[nScal],&scalAA[nScal],&scalBB[nScal])==3)
		    {
		      if (tcal==1)
			{
			  scalAA[nScal]/=2.0; // THIS IS BECAUSE OF A DIFFERENCE IN HOW TCAL AND SCAL ARE STORED
			  scalBB[nScal]/=2.0;
			}
		      nScal++;
		      if (nScal == MAX_CHAN_SCAL+1)
			{
			  printf("ERROR: Must increase MAX_CHAN_SCAL (currently set to %d)\n",MAX_CHAN_SCAL);
			  fclose(fin);
			  exit(1); // Should exit more gracefully
			}

		    }
		}
	      fclose(fin);
	      printf("Complete loading %s\n",scalFname);
	    }
	  

	  sdhdf_initialiseFile(inFile);
	  sdhdf_initialiseFile(outFile);

	  sdhdf_openFile(fname,inFile,1);
	  //	  inFile->fileID  = H5Fopen(fname,H5F_ACC_RDONLY,H5P_DEFAULT);
	  sprintf(oname,"%s.%s",fname,extension);
	  printf("Opening file >%s<\n",oname);
	  sdhdf_openFile(oname,outFile,3);
	  if (inFile->fileID!=-1) // Did we successfully open the file?
	    {
	      sdhdf_loadMetaData(inFile);
	      
	      printf("%-22.22s %-22.22s %-5.5s %-6.6s %-20.20s %-10.10s %-10.10s %-7.7s %3d\n",fname,inFile->primary[0].utc0,
		     inFile->primary[0].hdr_defn_version,inFile->primary[0].pid,inFile->beamHeader[0].source,inFile->primary[0].telescope,
		     inFile->primary[0].observer, inFile->primary[0].rcvr,inFile->beam[0].nBand);
	      
	      nBeam = inFile->nBeam;

	      sdhdf_allocateBeamMemory(outFile,nBeam);

	      for (b=0; b<nBeam; b++)
		{
		  printf("Processing beam %d/%d\n",b,nBeam-1);
		  inBandParams = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*inFile->beam[b].nBand);      
		  sdhdf_copyBandHeaderStruct(inFile->beam[b].bandHeader,inBandParams,inFile->beam[b].nBand);

		  // Do we need to extract the calibrator information?
		  if (nScal > 0)
		    {
		      int ci,cj,ck;
		      double onP1,offP1,onP2,offP2;
		      double onOffP1,onOffP2;
		      double freq;
		      int totNchan_cal=0;

		      for (ci=0;ci<inFile->beam[b].nBand;ci++)
			totNchan_cal+=inFile->beam[b].calBandHeader[ci].nchan;
		      
		      
		      printf("Calibrating\n");
		      // Could calibrate dump by dump, or average the entire data file
		      // FOR NOW **** just averaging all the calibration information in the file
		      cal_p1   = (float *)malloc(sizeof(float)*totNchan_cal);
		      cal_p2   = (float *)malloc(sizeof(float)*totNchan_cal);
		      cal_freq = (float *)malloc(sizeof(float)*totNchan_cal);
		      nc=0;

		      nRF1=nRF2=nRF3=0;
		      
		      
		      for (ci=0;ci<inFile->beam[b].nBand;ci++)
			{
			  sdhdf_loadBandData(inFile,b,ci,2);
			  sdhdf_loadBandData(inFile,b,ci,3);
		       
			  c_nchan = inFile->beam[b].calBandHeader[ci].nchan;
			  // printf("Number of channels = %d, number of dumps = %d\n",c_nchan,inFile->beam[b].calBandHeader[ci].ndump);
			  for (cj=0;cj<c_nchan;cj++)
			    {
			      onP1=offP1=onP2=offP2=onOffP1=onOffP2=0.0;
		
			      for (ck=0;ck<inFile->beam[b].calBandHeader[ci].ndump;ck++)
				{
				  onP1  += inFile->beam[b].bandData[ci].cal_on_data.pol1[cj+ck*c_nchan];
				  offP1 += inFile->beam[b].bandData[ci].cal_off_data.pol1[cj+ck*c_nchan];
				  onP2  += inFile->beam[b].bandData[ci].cal_on_data.pol2[cj+ck*c_nchan];
				  offP2 += inFile->beam[b].bandData[ci].cal_off_data.pol2[cj+ck*c_nchan];

				  onOffP1 += (inFile->beam[b].bandData[ci].cal_on_data.pol1[cj+ck*c_nchan] - inFile->beam[b].bandData[ci].cal_off_data.pol1[cj+ck*c_nchan]);		  
				  onOffP2 += (inFile->beam[b].bandData[ci].cal_on_data.pol2[cj+ck*c_nchan] - inFile->beam[b].bandData[ci].cal_off_data.pol2[cj+ck*c_nchan]);
				}
			      //			      printf("cal: %.5f %g %g %g %g\n",inFile->beam[b].bandData[ci].cal_on_data.freq[cj],onP1,offP1,onP2,offP2);

			      // Note: summing, not averaging above ***
			      //
			      freq = cal_freq[nc] = inFile->beam[b].bandData[ci].cal_on_data.freq[cj];

			      //			      cal_p1[nc] = (onP1-offP1);  
			      //			      cal_p2[nc] = (onP2-offP2);
			      cal_p1[nc] = onOffP1;
			      cal_p2[nc] = onOffP2;
			      //			      printf("freq = %g\n",freq);

			      if ((freq > 814 && freq < 818) ||
				  (freq > 994 && freq < 998) ||
				  (freq > 1068 && freq < 1073) ||
				  (freq > 1194 && freq < 1198) ||
				  (freq > 1295 && freq < 1305))
				{
				  if (cal_p1[nc] > 0 && cal_p2[nc] > 0)
				    {
				      scA_rf1_x[nRF1] = scB_rf1_x[nRF1] = log10(freq);
				      scA_rf1_y[nRF1] = log10(cal_p1[nc]);
				      scB_rf1_y[nRF1] = log10(cal_p2[nc]);
				      nRF1++;
				    }
				  //				  scA_rf1+=cal_p1[nc]; scB_rf1+=cal_p2[nc]; nRF1++;
				}
			      if ((freq > 1309 && freq < 1400) ||
				  (freq > 1444 && freq < 1454) ||
				  (freq > 1510 && freq < 1525) ||
				  (freq > 1635 && freq < 1660) ||
				  (freq > 1770 && freq < 1785) ||
				  (freq > 1830 && freq < 1840) ||
				  (freq > 1930 && freq < 1960) ||
				  (freq > 2020 && freq < 2040) ||
				  (freq > 2060 && freq < 2070) ||
				  (freq > 2175 && freq < 2195))
				{
				  if (cal_p1[nc] > 0 && cal_p2[nc] > 0)
				    {
				      scA_rf2_x[nRF2] = scB_rf2_x[nRF2] = log10(freq);
				      scA_rf2_y[nRF2] = log10(cal_p1[nc]);
				      scB_rf2_y[nRF2] = log10(cal_p2[nc]);
				      nRF2++;
				    }
				  
				  //				  scA_rf2+=cal_p1[nc]; scB_rf2+=cal_p2[nc]; nRF2++;
				}
			      if ((freq > 2520 && freq < 2550) ||
				  (freq > 2654 && freq < 2660) ||
				  (freq > 2710 && freq < 2730) ||
				  (freq > 2780 && freq < 2860) ||
				  (freq > 2900 && freq < 2980) ||
				  (freq > 3030 && freq < 3060) ||
				  (freq > 3080 && freq < 3120) ||
				  (freq > 3180 && freq < 3240) ||
				  (freq > 3290 && freq < 3320) ||
				  (freq > 3340 && freq < 3370) ||
				  (freq > 3415 && freq < 3435) ||
				  (freq > 3485 && freq < 3495) ||
				  (freq > 3590 && freq < 3605) ||
				  (freq > 3680 && freq < 3690) ||
				  (freq > 3826 && freq < 3829) ||
				  (freq > 3860 && freq < 3880) ||
				  (freq > 4010 && freq < 4015))
				{
				  if (cal_p1[nc] > 0 && cal_p2[nc] > 0)
				    {
				      scA_rf3_x[nRF3] = scB_rf3_x[nRF3] = log10(freq);
				      scA_rf3_y[nRF3] = log10(cal_p1[nc]);
				      scB_rf3_y[nRF3] = log10(cal_p2[nc]);
				      nRF3++;
				    }
				  //  scA_rf3+=cal_p1[nc]; scB_rf3+=cal_p2[nc]; nRF3++;
				}
			      
			      nc++;
			    }
			}

		      // Fit straight lines to the cal data
		      if (calFitN>0)
			{

			  if (nRF1 == 0)
			    printf("WARNING: Attempting to fit to RF1, but no data\n");
			  else
			    {
			      TKleastSquares_svd_noErr(scA_rf1_x,scA_rf1_y,nRF1,scA_rf1,calFitN,TKfitPoly);
			      TKleastSquares_svd_noErr(scB_rf1_x,scB_rf1_y,nRF1,scB_rf1,calFitN,TKfitPoly);
			    }

			  if (nRF2 == 0)
			    printf("WARNING: Attempting to fit to RF2, but no data\n");
			  else
			    {
			      TKleastSquares_svd_noErr(scA_rf2_x,scA_rf2_y,nRF2,scA_rf2,calFitN,TKfitPoly);
			      TKleastSquares_svd_noErr(scB_rf2_x,scB_rf2_y,nRF2,scB_rf2,calFitN,TKfitPoly);
			    }
			  if (nRF3 == 0)
			    printf("WARNING: Attempting to fit to RF3, but no data\n");
			  else
			    {
			      TKleastSquares_svd_noErr(scA_rf3_x,scA_rf3_y,nRF3,scA_rf3,calFitN,TKfitPoly);
			      TKleastSquares_svd_noErr(scB_rf3_x,scB_rf3_y,nRF3,scB_rf3,calFitN,TKfitPoly);
			    }
			}

		      {
			FILE *fout;
			int i0;
			fout = fopen("calOutput_rf1.dat","w");
			for (i0=0;i0<nRF1;i0++)
			  fprintf(fout,"calOutput: %g %g %g %g %g %g %g %g\n",scA_rf1_x[i0],scA_rf1_y[i0],scB_rf1_x[i0],scB_rf1_y[i0],scA_rf1[0],scB_rf1[0],scA_rf1[1],scB_rf1[1]);
			fclose(fout);

			fout = fopen("calOutput_rf2.dat","w");
			for (i0=0;i0<nRF2;i0++)
			  fprintf(fout,"calOutput: %.6f %.6f %.6f %.6f %.4f %.4f %.4f %.4f %.4f %.4f\n",scA_rf2_x[i0],scA_rf2_y[i0],scB_rf2_x[i0],scB_rf2_y[i0],scA_rf2[0],scB_rf2[0],scA_rf2[1],scB_rf2[1],scA_rf2[2],scB_rf2[2]);
			fclose(fout);

			
			fout = fopen("calOutput_rf3.dat","w");
			for (i0=0;i0<nRF3;i0++)
			  fprintf(fout,"calOutput: %g %g %g %g %g %g %g %g\n",scA_rf3_x[i0],scA_rf3_y[i0],scB_rf3_x[i0],scB_rf3_y[i0],scA_rf3[0],scB_rf3[0],scA_rf3[1],scB_rf3[1]);
			fclose(fout);
		      }
		      
		    
		      printf("Loaded cal %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %d %d %d\n",scA_rf1[0],scB_rf1[0],scA_rf1[1],scB_rf1[1],scA_rf1[2],scB_rf1[2],
			     scA_rf2[0],scB_rf2[0],scA_rf2[1],scB_rf2[1],scA_rf2[2],scB_rf2[2],scA_rf3[0],scB_rf3[0],scA_rf3[1],scA_rf3[2],scB_rf3[2],
			     scB_rf3[1],nRF1,nRF2,nRF3);
		      //			  exit(1);

		    }

	       
		  // Copy the bands
		  nc2=0;
		  for (ii=0;ii<inFile->beam[b].nBand;ii++)
		    {
		      printf(" ... Processing subband %d/%d\n",ii,inFile->beam[b].nBand-1);
		      // If "timeScrunch" is set then need to
		      // Replace the first spectral dump in the file with the scrunched information (note integration times
		      //  in case one spectral dump is shorter than the other ones??)
		      // Remove the other spectral dumps
		      // Update all relevant meta-data
		      //   ... integration time
		      //   ... number of spectral dumps
		      //   ... time of spectral dump
		      
		      // I don't want to copy the data column as I cannot resize it.
		      nchan = inFile->beam[b].bandHeader[ii].nchan;
		      npol  = inFile->beam[b].bandHeader[ii].npol;

		      // FIX ME
		      //		      sdhdf_loadFrequencyAttributes(inFile,ii);
		      //		      sdhdf_loadDataAttributes(inFile,ii);
		      
		      //		      strcpy(outFile->frequency_attr.frame,inFile->frequency_attr.frame);
		      //		      strcpy(outFile->frequency_attr.unit,inFile->frequency_attr.unit);
		      //		      strcpy(outFile->data_attr.unit,inFile->data_attr.unit);		     		      


		      if (nchanRequested > 0) // Define fav
			{
			  fAv = (int)(nchan/nchanRequested);
			  printf("Setting number of channels to average to %d (original number of channels = %d)\n",fAv,nchan);
			}
		      
		      if (fScrunch==1)
			fAv = nchan;
		   
		      if (fchbw > 0)
			{
			  double chbw;
			  sdhdf_loadBandData(inFile,b,ii,1);
			  chbw = fabs(inFile->beam[b].bandData[ii].astro_data.freq[1]-inFile->beam[b].bandData[ii].astro_data.freq[0]); // Assuming time ordered, evenly spaced channels -- FIX ME
			  fAv = (int)(fabs(fchbw)/chbw);
			  printf("Frequency averaging by %d channels\n",fAv);
			}
		      
		      if (fAv < 2)
			out_nchan = nchan;
		      else
			out_nchan = (int)(nchan/fAv);
		      
		      if (pol == 1)
			out_npol=1;
		      else
			out_npol  = npol;


		      if (tdumpRequest > 0)
			{
			  float tdump;
			  tdump = inFile->beam[b].bandHeader[ii].dtime;
			  printf("tdump = %g request = %g\n",tdump,tdumpRequest);
			}


		      
		      if (tScrunch==1)
			out_ndump = 1;
		      else
			out_ndump = inFile->beam[b].bandHeader[ii].ndump;
		    

		      if (nScal > 0)
			{
			  int cj,ck;
			  double onP1,offP1,onP2,offP2,onP3,offP3,onP4,offP4;
			  double freq;
			  
			  printf("Calibrating\n");
			  // Could calibrate dump by dump, or average the entire data file
			  // FOR NOW **** just averaging all the calibration information in the file
			  sdhdf_loadBandData(inFile,b,ii,2);
			  sdhdf_loadBandData(inFile,b,ii,3);
			  
			  c_nchan = inFile->beam[b].calBandHeader[ii].nchan;
			  printf("Number of channels = %d\n",c_nchan);
			  for (cj=0;cj<c_nchan;cj++)
			    {
			      onP1=offP1=onP2=offP2=onP3=offP3=onP4=offP4=0.0;
			      //			      printf("ndump = %d\n",inFile->beam[b].calBandHeader[ii].ndump);
			      for (ck=0;ck<inFile->beam[b].calBandHeader[ii].ndump;ck++)
				{
				  onP1  += inFile->beam[b].bandData[ii].cal_on_data.pol1[cj+ck*c_nchan];
				  offP1 += inFile->beam[b].bandData[ii].cal_off_data.pol1[cj+ck*c_nchan];
				  onP2  += inFile->beam[b].bandData[ii].cal_on_data.pol2[cj+ck*c_nchan];
				  offP2 += inFile->beam[b].bandData[ii].cal_off_data.pol2[cj+ck*c_nchan];
 				  onP3  += inFile->beam[b].bandData[ii].cal_on_data.pol3[cj+ck*c_nchan];
				  offP3 += inFile->beam[b].bandData[ii].cal_off_data.pol3[cj+ck*c_nchan];
				  onP4  += inFile->beam[b].bandData[ii].cal_on_data.pol4[cj+ck*c_nchan];
				  offP4 += inFile->beam[b].bandData[ii].cal_off_data.pol4[cj+ck*c_nchan];
				}
			      
			      onP1 /= (inFile->beam[b].calBandHeader[ii].ndump*divideCal);
			      offP1 /= (inFile->beam[b].calBandHeader[ii].ndump*divideCal);
			      onP2 /= (inFile->beam[b].calBandHeader[ii].ndump*divideCal);
			      offP2 /= (inFile->beam[b].calBandHeader[ii].ndump*divideCal);
			      onP3 /= (inFile->beam[b].calBandHeader[ii].ndump*divideCal);
			      offP3 /= (inFile->beam[b].calBandHeader[ii].ndump*divideCal);
			      onP4 /= (inFile->beam[b].calBandHeader[ii].ndump*divideCal);
			      offP4 /= (inFile->beam[b].calBandHeader[ii].ndump*divideCal);
			      
			      // Note: summing, not averaging above ***
			      //
			      calAA   = onP1-offP1;
			      calBB   = onP2-offP2;
			      calReAB = onP3-offP3;
			      calImAB = onP4-offP4;

			      sysGain_freq[cj] = inFile->beam[b].bandData[ii].cal_on_data.freq[cj];
			      diffGain2[cj]    = (calAA-calBB)/(calAA+calBB);
			      diffPhase[cj]    = atan2(calImAB,calReAB);
			      scale1[cj] = calAA;
			      scale2[cj] = calBB;
			      scale3[cj] = calReAB;
			      scale4[cj] = calImAB;
			      
			      getScal(sysGain_freq[cj],scalFreq,scalAA,scalBB,nScal,&scalAA_val,&scalBB_val);
			      //  g_e[cj] = (calAA+calBB)/(scalAA_val+scalBB_val);

			      // Note should use -tcal isntead of -scal if using TCAL measurements
			      g_e[cj] = (1.0)/(scalAA_val+scalBB_val); // Check THIS VERY CAREFULLY -- NOW IN cals/Jy
			      //			      g_e[cj] = (8.0)/(scalAA_val+scalBB_val); // Check THIS VERY CAREFULLY -- NOW IN cals/Jy
			      //			      printf("WARNING SCALING BY A FACTOR OF 8 -- CHECK THIS *****\n");
			      
			      //			      printf("g_e = %g\n",g_e[cj]);
			      if (ii==5)
				printf("Processing freq = %.6f %g %g %g\n",sysGain_freq[cj],diffGain2[cj],diffPhase[cj]*180/M_PI,g_e[cj]);

			      if (calFitN>0)
				{
				  scA_rf1_val = pow(10,scA_rf1[0] + scA_rf1[1]*log10(sysGain_freq[cj]) + scA_rf1[2]*pow(log10(sysGain_freq[cj]),2));
				  scB_rf1_val = pow(10,scB_rf1[0] + scB_rf1[1]*log10(sysGain_freq[cj]) + scB_rf1[2]*pow(log10(sysGain_freq[cj]),2));
				  scA_rf2_val = pow(10,scA_rf2[0] + scA_rf2[1]*log10(sysGain_freq[cj]) + scA_rf2[2]*pow(log10(sysGain_freq[cj]),2));
				  scB_rf2_val = pow(10,scB_rf2[0] + scB_rf2[1]*log10(sysGain_freq[cj]) + scB_rf2[2]*pow(log10(sysGain_freq[cj]),2));
				  scA_rf3_val = pow(10,scA_rf3[0] + scA_rf3[1]*log10(sysGain_freq[cj]) + scA_rf3[2]*pow(log10(sysGain_freq[cj]),2));
				  scB_rf3_val = pow(10,scB_rf3[0] + scB_rf3[1]*log10(sysGain_freq[cj]) + scB_rf3[2]*pow(log10(sysGain_freq[cj]),2));

				  //				  printf("[%d] Comparison %.6g %.6g %.6g %.6g %.6g %.6g\n",calFitN,sysGain_freq[cj],onP1-offP1,cal_p1[nc2],scA_rf1_val,scA_rf2_val,scA_rf3_val);
				  							      
				  if (sysGain_freq[cj] < 1344) // HARDCODE TO PARKES
				    {
				      sysGain_p1[cj] = scA_rf1_val/scalAA_val;
				      sysGain_p2[cj] = scB_rf1_val/scalBB_val;
				    }
				  else if (sysGain_freq[cj] < 2368) // HARDCODE TO PARKES
				    {
				      sysGain_p1[cj] = scA_rf2_val/scalAA_val;
				      sysGain_p2[cj] = scB_rf2_val/scalBB_val;
				    }
				  else if (sysGain_freq[cj] <= 4032) // HARDCODE TO PARKES
				    {
				      sysGain_p1[cj] = scA_rf3_val/scalAA_val;
				      sysGain_p2[cj] = scB_rf3_val/scalBB_val;
				    }				  
				}
			      else
				{
				  sysGain_p1[cj] = cal_p1[nc2]/scalAA_val;
				  sysGain_p2[cj] = cal_p2[nc2]/scalBB_val;				  
				}
			  
			      nc2++;

			      /*
			      printf("Gain: %.6f %g %g %g %g %g %g %g %g %g %g [ %g %g ]\n",
				     sysGain_freq[cj],sysGain_p1[cj],sysGain_p2[cj],scalAA_val,scalBB_val,
				     scA_rf1_val,scB_rf1_val, scA_rf2_val,scB_rf2_val,scA_rf3_val,scB_rf3_val,
				     scA_rf1[0],scA_rf1[1]);
			      */
			    }
			}
		    		    
		      outObsParams = (sdhdf_obsParamsStruct *)malloc(sizeof(sdhdf_obsParamsStruct)*out_ndump);      
		      if (tScrunch==1)
			sdhdf_copySingleObsParams(inFile,b,ii,0,&outObsParams[0]);
		      else
			{
			  int kk;
			  for (kk=0;kk<out_ndump;kk++)
			    sdhdf_copySingleObsParams(inFile,b,ii,kk,&outObsParams[kk]);
			}
		      nsd = inFile->beam[b].bandHeader[ii].ndump;
		      if (!(in_data = (float *)malloc(sizeof(float)*nchan*npol*inFile->beam[b].bandHeader[ii].ndump)))
			{
			  printf("ERROR: Unable to allocate memory for in_data: %d\n",nchan*npol*nsd);
			  exit(1);
			}
		      in_freq   = (double *)malloc(sizeof(double)*nchan);
		      out_freq  = (float *)malloc(sizeof(float)*out_nchan);
		      out_data  = (float *)calloc(sizeof(float),out_nchan*out_npol*out_ndump);
		      out_Tdata = (float *)calloc(sizeof(float),nchan*npol*out_ndump);
		      out_Fdata = (float *)calloc(sizeof(float),out_nchan*npol*out_ndump);
		      dataWts   = (float *)calloc(sizeof(float),out_nchan*out_ndump);
		      tav=0;
		   
		      sdhdf_loadBandData(inFile,b,ii,1);
		      for (j=0;j<nchan;j++)
			{
			  in_freq[j] = inFile->beam[b].bandData[ii].astro_data.freq[j];
			  if (fAv < 2)
			    out_freq[j] = in_freq[j];
			}
		      if (strcmp(scaleFile,"NULL")!=0)
			{
			  float sclAA,sclBB;
			  int kk,k0;
			  //			  printf("Applying scaling\n");
			  // Apply the scaling
			  for (j=0;j<inFile->beam[b].bandHeader[ii].ndump;j++)
			    {
			      k0=0;
			      for (k=0;k<nchan;k++)
				{
				  sclAA = sclBB = 0.0;
				  for (kk=k0;kk<nScale;kk++)
				    {
				      if (in_freq[k] <= freq[kk])
					{
					  sclAA = 1./linearInterpolate(freq,scaleAA,in_freq[k],kk,nScale);
					  sclBB = 1./linearInterpolate(freq,scaleBB,in_freq[k],kk,nScale);
					  if (k==0) k0=kk;
					  break;
					}				     
				    }
				  //				  printf("WARNING: ASSUMING CAL HAS 128 CHANNELS/SUBBAND\n");
				  //				  sclAA *= (nchan/128);
				  //				  sclBB *= (nchan/128);
				  inFile->beam[b].bandData[ii].astro_data.pol1[k+j*nchan]*=sclAA;
				  inFile->beam[b].bandData[ii].astro_data.pol2[k+j*nchan]*=sclBB;
				}
			    }
			  //			  printf("Fixed the scaling\n");
			}
		      // Time averaging
		      if (tScrunch==1)
			{
			  int ignore=0;
			  int jj;
			  int nsum=0;
			  double swt=0,wt;
			  
			  for (k=0;k<nchan;k++)
			    {
			      swt=0;
			      tav=0;
			      for (j=0;j<inFile->beam[b].bandHeader[ii].ndump;j++)
				{
				  ignore=0;
				  
				  for (jj=0;jj<nIgnoreSD;jj++)
				    {
				      if (ignoreSD[jj] == j)
					{ignore=1; break;}
				    }
				  if (ignore==0)
				    {
				      tav += inFile->beam[b].bandHeader[ii].dtime;
				      nsum++;
				      //				  printf("Updating time averaged to %g\n",tav);
				      
				      wt  = inFile->beam[b].bandData[ii].astro_data.dataWeights[j*nchan+k];
				      
				      if (wt!=0)
					{
					  if (npol==1)
					    out_Tdata[k] += wt*inFile->beam[b].bandData[ii].astro_data.pol1[k+j*nchan]; ///(float)nsd;
					  else if (npol==2)
					    {
					      out_Tdata[k]         += wt*inFile->beam[b].bandData[ii].astro_data.pol1[k+j*nchan]; // /(float)nsd;
					      out_Tdata[k+nchan]   += wt*inFile->beam[b].bandData[ii].astro_data.pol2[k+j*nchan]; // /(float)nsd;
					    }
					  else
					    {
					      out_Tdata[k]         += wt*inFile->beam[b].bandData[ii].astro_data.pol1[k+j*nchan]; // /(float)nsd;
					      out_Tdata[k+nchan]   += wt*inFile->beam[b].bandData[ii].astro_data.pol2[k+j*nchan]; // /(float)nsd;
					      out_Tdata[k+2*nchan] += wt*inFile->beam[b].bandData[ii].astro_data.pol3[k+j*nchan]; // /(float)nsd;
					      out_Tdata[k+3*nchan] += wt*inFile->beam[b].bandData[ii].astro_data.pol4[k+j*nchan]; // /(float)nsd;  
					    }
					  
					  
					  //				      printf("Updating dataWts %d %d %g\n",j,k,inFile->beam[b].bandData[ii].astro_data.dataWeights[j*nchan+k]);
					  // Update the data weights
					  // GEORGE: SHOULD CHECK IF THE ORIGINAL FILE DIDN'T HAVE WEIGHTS ***
					  // Note that we don't have j*out_nchan here as we're time scrunching
					  dataWts[k] += wt;
					  swt += wt;
					}
					  //				      printf("Done\n");
				    }
				}
			      if (npol==1)
				{if (swt > 0) {out_Tdata[k] /= swt;} else {out_Tdata[k] = 0;}}
			      else if (npol==2)
				{
				  if (swt > 0)
				    {
				      out_Tdata[k]         /= swt;
				      out_Tdata[k+nchan]   /= swt;
				    }
				  else
				    {
				      out_Tdata[k]         = 0;
				      out_Tdata[k+nchan]   = 0;				      
				    }
				}
			      else
				{
				  if (swt > 0)
				    {
				      out_Tdata[k]         /= swt;
				      out_Tdata[k+nchan]   /= swt;
				      out_Tdata[k+2*nchan] /= swt;
				      out_Tdata[k+3*nchan] /= swt;
				    }
				  else
				    {
				      out_Tdata[k]         = 0;
				      out_Tdata[k+nchan]   = 0;
				      out_Tdata[k+2*nchan] = 0;
				      out_Tdata[k+3*nchan] = 0;
				    }
				}
			      
			      //			      printf("data_wts = %g\n",dataWts[k]);
			    }


			
			  
			  // Scale the output
			  // GEORGE **** THIS SHOULD BE EARLIER AS YOU MAY NOT CHOOSE TO TIME SCRUNCH *****
			  if (nScal > 0 && stokes==0)
			    {
			      double gainVal1,gainVal2;
			      double *interpY1,*interpY2;
			      
			      if (interp==2) // Cubic spline
				{
				  interpY1 = (double *)malloc(sizeof(double)*nchan);
				  interpY2 = (double *)malloc(sizeof(double)*nchan);
				  printf("Setting up cubic spline\n");
				  TKcmonot(c_nchan,sysGain_freq,sysGain_p1,yd1);
				  TKcmonot(c_nchan,sysGain_freq,sysGain_p2,yd2);
				  TKspline_interpolate(c_nchan,sysGain_freq,sysGain_p1,yd1,in_freq,interpY1,nchan);
				  TKspline_interpolate(c_nchan,sysGain_freq,sysGain_p2,yd2,in_freq,interpY2,nchan);
				  printf("Complete setting up cubic spline\n");
				}
			      for (k=0;k<nchan;k++)
				{
				  //
				  // Note that this is using the getScal interpolation routine, but is getting sysGain instead
				  //
				  if (interp==1)
				    getScal(in_freq[k],sysGain_freq,sysGain_p1,sysGain_p2,c_nchan,&gainVal1,&gainVal2);
				  else if (interp==2)
				    {
				      gainVal1 = interpY1[k];
				      gainVal2 = interpY2[k];
				    }
				  //				  printf("Scale factors: %.6f %g %g\n",in_freq[k],gainVal1,gainVal2);
				  // Only scaling 2 polarisations *** <<<
				  
				  if (gainVal1 <= 0) gainVal1=1e9; // SHOULD SET MORE SENSIBLY ****
				  if (gainVal2 <= 0) gainVal2=1e9; // SHOULD SET MORE SENSIBLY ****
				  
				  out_Tdata[k]       *= nchan/c_nchan/gainVal1/calBin; ///nsum; // Note scaling by number of spectral dumps as each one now should be in Jy
				  out_Tdata[k+nchan] *= nchan/c_nchan/gainVal2/calBin; ///nsum;				  
				}
			      if (interp==2)
				{
				  free(interpY1);
				  free(interpY2);
				}
			    }
			}
		      else
			{
			  for (j=0;j<inFile->beam[b].bandHeader[ii].ndump;j++)
			    {
			      for (k=0;k<nchan;k++)
				{
				  if (npol==1)
				    out_Tdata[k+j*nchan*npol] = inFile->beam[b].bandData[ii].astro_data.pol1[k+j*nchan]; // /(float)nsd;
				  else if (npol==2)
				    {
				      out_Tdata[k+j*nchan*npol]         = inFile->beam[b].bandData[ii].astro_data.pol1[k+j*nchan]; ///(float)nsd;
				      out_Tdata[k+j*nchan*npol+nchan]   = inFile->beam[b].bandData[ii].astro_data.pol2[k+j*nchan]; ///(float)nsd;
				    }
				  else if (npol==4)
				    {
				      out_Tdata[k+j*nchan*npol]         = inFile->beam[b].bandData[ii].astro_data.pol1[k+j*nchan]; ///(float)nsd;
				      out_Tdata[k+j*nchan*npol+nchan]   = inFile->beam[b].bandData[ii].astro_data.pol2[k+j*nchan]; ///(float)nsd;
				      out_Tdata[k+j*nchan*npol+2*nchan] = inFile->beam[b].bandData[ii].astro_data.pol3[k+j*nchan]; ///(float)nsd;
				      out_Tdata[k+j*nchan*npol+3*nchan] = inFile->beam[b].bandData[ii].astro_data.pol4[k+j*nchan]; ///(float)nsd; 
				    }
				}	
			    }
			}
		    
		      // Frequency averaging (note: weighted mean)
		      if (fAv >= 2)
			{
			  double wt;
			  double swt;
			  printf("In averaging\n");
			  for (j=0;j<out_ndump;j++)
			    {
			      kp=0;
			      for (k=0;k<nchan;k+=fAv)
				{				  
				  av1 = av2 = av3 = av4 = avFreq = 0.0;
				  swt=0;
				  dataWts[j*out_nchan+kp] = 0; 
				  for (l=k;l<k+fAv;l++)
				    {
				      // This needs fixing -- because dataWeights not set for the output ndump if time scrunched
				      wt = inFile->beam[b].bandData[ii].astro_data.dataWeights[j*nchan+l];
				      //				      printf("wt = %g\n",wt);
				      avFreq += in_freq[l]; // Should be a weighted sum  *** FIX ME
				      if (npol==1)
					av1 += wt*out_Tdata[l+j*npol*nchan];
				      else if (npol==2)
					{
					  av1 += wt*out_Tdata[l+j*npol*nchan];
					  av2 += wt*out_Tdata[l+j*npol*nchan+nchan];
					}
				      else
					{
					  //					  printf("Got here with %g %g %g\n",wt,out_Tdata[l+j*npol*nchan],av1);				  
					  av1 += wt*out_Tdata[l+j*npol*nchan];
					  av2 += wt*out_Tdata[l+j*npol*nchan+1*nchan];
					  av3 += wt*out_Tdata[l+j*npol*nchan+2*nchan];
					  av4 += wt*out_Tdata[l+j*npol*nchan+3*nchan];
					}
				      //  printf("Writing to %d %d nchan = %d l = %d k = %d\n",j,kp,nchan,l,k);
				      //  printf("Input weight = %g\n",inFile->beam[b].bandData[ii].astro_data.dataWeights[j*nchan+l]);
				      dataWts[j*out_nchan+kp] += wt; //inFile->beam[b].bandData[ii].astro_data.dataWeights[j*nchan+l];
				      swt += wt;
				      //				      printf("Done write\n");
				    }
				
				  //				  printf("Updating the sum\n");
				  if (j==0) out_freq[kp] = avFreq/fAv;
				  if (npol==1)
				    {
				      if (swt > 0)
					out_Fdata[kp+out_nchan*j*npol]               = av1/swt; 
				      else
					out_Fdata[kp+out_nchan*j*npol]               = 0; 
				    }
				  else if (npol==2)
				    {
				      if (swt > 0)
					{
					  out_Fdata[kp+out_nchan*j*npol]             = av1/swt; 
					  out_Fdata[kp+out_nchan*j*npol+out_nchan]   = av2/swt; 
					}
				      else
					{
					  out_Fdata[kp+out_nchan*j*npol]             = 0; 
					  out_Fdata[kp+out_nchan*j*npol+out_nchan]   = 0; 
					}
				    }
				  else
				    {
				      if (swt > 0)
					{
					  out_Fdata[kp+out_nchan*j*npol]             = av1/swt; 
					  out_Fdata[kp+out_nchan*j*npol+out_nchan]   = av2/swt; 
					  out_Fdata[kp+out_nchan*j*npol+2*out_nchan] = av3/swt; 
					  out_Fdata[kp+out_nchan*j*npol+3*out_nchan] = av4/swt; 
					}
				      else
					{
					  out_Fdata[kp+out_nchan*j*npol]             = 0; 
					  out_Fdata[kp+out_nchan*j*npol+out_nchan]   = 0; 
					  out_Fdata[kp+out_nchan*j*npol+2*out_nchan] = 0; 
					  out_Fdata[kp+out_nchan*j*npol+3*out_nchan] = 0; 
					}
				    }
				  kp++;
				  
				}
			      printf("Processed dump %d\n",j);
			    }
			  printf("Finished frequency averaging\n");
			}
		      else
			{
			  // Should use memcpy
			  for (j=0;j<out_ndump;j++)
			    {
			      for (k=0;k<out_nchan;k++)
				{
				  // We've already set the weights above ???  FIX ME *** CHECK *****
				  if (tScrunch!=1)
				    dataWts[j*out_nchan+k] = inFile->beam[b].bandData[ii].astro_data.dataWeights[j*nchan+j];
				  if (npol==1)
				    out_Fdata[k+out_nchan*j*4]             = out_Tdata[k+out_nchan*j*4];  // *4 NEED FIXING
				  else if (npol==2)
				    {
				      out_Fdata[k+out_nchan*j*4]             = out_Tdata[k+out_nchan*j*4];  // *4 NEED FIXING
				      out_Fdata[k+out_nchan*j*4+out_nchan]   = out_Tdata[k+out_nchan*j*4+out_nchan];
				    }
				  else if (npol==4)
				    {
				      out_Fdata[k+out_nchan*j*4]             = out_Tdata[k+out_nchan*j*4];
				      out_Fdata[k+out_nchan*j*4+out_nchan]   = out_Tdata[k+out_nchan*j*4+out_nchan];
				      out_Fdata[k+out_nchan*j*4+2*out_nchan] = out_Tdata[k+out_nchan*j*4+2*out_nchan];
				      out_Fdata[k+out_nchan*j*4+3*out_nchan] = out_Tdata[k+out_nchan*j*4+3*out_nchan];
				    }
				}
			    }
			}
		      // This was moved lower because of the need to read the weights
		      printf("Releasing band data\n");
		      sdhdf_releaseBandData(inFile,b,ii,1); 		      
		      printf("Creating final data\n");
		      // Now create final data
		      // Do polarisation summing as required
		      if (pol==0)
			{
			  for (j=0;j<out_ndump;j++)
			    {
			      for (k=0;k<out_nchan;k++)
				{
				  if (npol==1)
				    out_data[k+out_nchan*j*out_npol]               = out_Fdata[k+out_nchan*j*npol];
				  else if (npol==2)
				    {
				      out_data[k+out_nchan*j*out_npol]             = out_Fdata[k+out_nchan*j*npol];
				      out_data[k+out_nchan*j*out_npol+out_nchan]   = out_Fdata[k+out_nchan*j*npol+out_nchan];
				    }
				  else if (npol==4) // Remember npol = input npol 
				    {
				      out_data[k+out_nchan*j*out_npol]             = out_Fdata[k+out_nchan*j*npol];
				      out_data[k+out_nchan*j*out_npol+out_nchan]   = out_Fdata[k+out_nchan*j*npol+out_nchan];
				      out_data[k+out_nchan*j*out_npol+2*out_nchan] = out_Fdata[k+out_nchan*j*npol+2*out_nchan];
				      out_data[k+out_nchan*j*out_npol+3*out_nchan] = out_Fdata[k+out_nchan*j*npol+3*out_nchan];
				    }
				}
			    }
			}
		      else if (pol==1) 
			{
			  for (j=0;j<out_ndump;j++)
			    {
			      for (k=0;k<out_nchan;k++)
				{
				  if (polAv==1)
				    out_data[k+out_nchan*j] = (out_Fdata[k+out_nchan*j*4]+out_Fdata[k+out_nchan*j*4+out_nchan])/2.;
				  else
				    out_data[k+out_nchan*j] = (out_Fdata[k+out_nchan*j*4]+out_Fdata[k+out_nchan*j*4+out_nchan]); 
				}
			      // out_data[k+out_nchan*j*4] = (out_Fdata[k+out_nchan*j*4]+out_Fdata[k+out_nchan*j+4+out_nchan])/2.;
			    }
			}
		      else if (pol==2)
			{
			  for (j=0;j<out_ndump;j++)
			    {
			      for (k=0;k<out_nchan;k++)
				{
				  out_data[k+out_nchan*j*4]             = out_Fdata[k+out_nchan*j*4];
				  out_data[k+out_nchan*j*4+out_nchan]   = out_Fdata[k+out_nchan*j*4+out_nchan];
				}
			    }
			}
		      else if (pol==3)
			{
			  for (j=0;j<out_ndump;j++)
			    {
			      for (k=0;k<out_nchan;k++)
				out_data[k+out_nchan*j*4] = (out_Fdata[k+out_nchan*j*4]+out_Fdata[k+out_nchan*j+4+out_nchan]);
			    }
			}
		      // Conversion to Stokes
		   
		      if (stokes==1)
			{
			  float p1,p2,p3,p4;
			  double t_p3,t_p4;
			  float mI,mQ,mU,mV;
			  float aI,aQ,aU,aV;
			  
			  double phi,dgain2;
			  float marray[4][4],mueller[4][4],mueller_inv[4][4];
			  float mastro[4][4];
			  float mcc[4][4];
			  float mpa[4][4];
			  float mfeed[4][4];
			  float mhand[4][4];
			  float tmatrix[4][4];
			  double paraAng;
			  //			  double feedGamma=-90*M_PI/180;
			  double feedGamma=0;      //-90*M_PI/180;
			  double ge;
			  double temp;

			  int sI,sJ;
			  float stokesM[4],stokesA[4];
			  
			  printf("Stokes In here with %d\n",out_ndump);
			  printf("WARNING: cal bin set to: %d\n",calBin);

			  // Check for Nchannel changes because of an extraction


			  if (origNchan == -1)
			    {
			      chbw = in_freq[1]-in_freq[0];
			      if (fabs(128.0 / chbw - nchan) > 1e-3) // hardcoding band width and 1e-3 values here ****
				{
				  printf("WARNING: The number of channels seems to have changed since the original file produced\n");
				  printf("Setting nchan to %d\n",(int)(128/chbw+0.5));
				  origNchan = (int)(128/chbw+0.5);				    
				}
				
			    }
			  
			  for (j=0;j<out_ndump;j++)
			    {
			      
			      for (k=0;k<out_nchan;k++)
				{
				  // FD_HAND = -1;
				  // Q1: WHAT SHOULD THIS BE
				  /*
				  
				  p2 = out_data[k+out_nchan*j*4];   
				  p1 = out_data[k+out_nchan*j*4+out_nchan];				  
				  p3 = out_data[k+out_nchan*j*4+2*out_nchan];
				  p4 = out_data[k+out_nchan*j*4+3*out_nchan];
				  */

				
				  p1 = out_data[k+out_nchan*j*4];   
				  p2 = out_data[k+out_nchan*j*4+out_nchan];
				  p3 = out_data[k+out_nchan*j*4+2*out_nchan];
				  p4 = out_data[k+out_nchan*j*4+3*out_nchan];

				  p1 /= divideAstro;
				  p2 /= divideAstro;
				  p3 /= divideAstro;
				  p4 /= divideAstro;
				  
				  //				  sdhdf_convertStokes(p1,p2,p3,p4,&out_data[k],&out_data[k+out_nchan],&out_data[k+2*out_nchan],&out_data[k+3*out_nchan]);

				  // Pre-stabilise the system				  
				  getScal(in_freq[k],sysGain_freq,scale1,scale2,c_nchan,&useScale1,&useScale2);
				  getScal(in_freq[k],sysGain_freq,scale3,scale4,c_nchan,&useScale3,&useScale4);
				  p1/=useScale1;
				  p2/=useScale2;

				  t_phase = atan2(useScale4,useScale3);
				  t_len = sqrt(pow(useScale3,2) + pow(useScale4,2));
				  t_p3 = p3;
				  t_p4 = p4;
				  p3 = (t_p3*cos(t_phase) + t_p4*sin(t_phase))/t_len;
				  p4 = (t_p4*cos(t_phase) - t_p3*sin(t_phase))/t_len;
				  //				  printf("Scaling: %g %g %g %g %g %g %g %g %g\n",sysGain_freq[cj],scale1[cj],scale2[cj],scale3[cj],scale4[cj],t_phase,t_len,calReAB,calImAB);

				  

				  //				  printf("Converting Stokes %g %g %g %g %g %g %g %g\n",p1,p2,p3,p4,useScale1,useScale2,useScale3,useScale4);
				  sdhdf_convertStokes(p1,p2,p3,p4,&mI,&mQ,&mU,&mV);

				  
				  stokesM[0] = mI; ///mI /8/25.;  // REMOVE THE /8 *******
				  stokesM[1] = mQ; // /8/25.;
				  stokesM[2] = mU; // /8/25.;
				  stokesM[3] = mV; // /8/25.; 


				  
				  if (in_freq[k] > 1400 && in_freq[k] < 1401)
				    {
				      printf("Stokes in: (mI,mQ,mU,mV) = (%g,%g,%g,%g)\n",stokesM[0],stokesM[1],stokesM[2],stokesM[3]);
				    }
				  
				  //				  getScal(in_freq[k],sysGain_freq,diffGain2,diffPhase,c_nchan,&dgain2,&phi);

				  // Just need g_e
				  getScal(in_freq[k],sysGain_freq,g_e,diffPhase,c_nchan,&ge,&phi); // Note interpolating two together and so need to pass phi again - FIX THIS
				  // ge = 1; // FIX THIS
				  // printf("diffGain: %.6f %g %g %g\n",in_freq[k],dgain2,phi,ge);

				  //
				  // ge/=(double)calBin;

				  // Should set the outFile, or search the correct place in the inFile???? **** FIX
				  if (j==0 && k==0)
				    printf("WARNING: using Parallactic angle: %g\n",inFile->beam[b].bandData[ii].astro_obsHeader[0].paraAngle);
				  paraAng = inFile->beam[b].bandData[ii].astro_obsHeader[0].paraAngle*M_PI/180.;

				  if (setPara==1)
				    {
				      printf("Using command line parallactic angle %g degrees\n",usePara);
				      paraAng = usePara*M_PI/180.0;
				    }


				  // Now apply the receiver Mueller matrix
				  
				  //				  sdhdf_setGain2Phase(marray,dgain2,phi);

				  //				  if (in_freq[k] > 3200 && in_freq[k] < 3201)
				  //				    {
				  //				      printf("amplifiers: %g (%g %g)\n",in_freq[k],dgain2*2,phi);
				  //				      displayMatrix_4x4(marray);
				  //				    }
				  sdhdf_setIdentity_4x4(marray);
				  // Parallactic angle
				  sdhdf_setParallacticAngle(mpa,paraAng);
				  sdhdf_setParallacticAngle(mastro,M_PI/4.);

				  //				  if (in_freq[k] > 1400 && in_freq[k] < 1401)
				  //				    {
				  //				      printf("msky: %g (%g)\n",in_freq[k],paraAng);
				  //				      displayMatrix_4x4(mpa);
				  //				    }
				
				  sdhdf_setIdentity_4x4(mhand);
				  sdhdf_setIdentity_4x4(mcc);
				  sdhdf_setIdentity_4x4(mfeed);
				  //				  mhand[1][1] = -1;  // FD_HAND = -1
				  //				  mhand[3][3] = -1;
				  
				  // Feed
				  //				  sdhdf_setFeed(mfeed,+45*M_PI/180.0);

				  // Cross coupling
				  //				  sdhdf_setIdentity_4x4(mcc);

				  // CHECK THESE

				  /*
				  cc_a = cc_eps1*cos(cc_phi1) + cc_eps2*cos(cc_phi2);
				  cc_b = cc_eps1*sin(cc_phi1) + cc_eps2*sin(cc_phi2);
				  cc_c = cc_eps1*cos(cc_phi1) - cc_eps2*cos(cc_phi2);
				  cc_d = cc_eps1*sin(cc_phi1) - cc_eps2*sin(cc_phi2);

				  mcc[2][0] = cc_a;
				  mcc[3][0] = cc_b;
				  mcc[2][1] = cc_c;
				  mcc[3][1] = cc_d;

				  mcc[0][2] = cc_a;
				  mcc[1][2] = -cc_c;
				  mcc[0][3] = cc_b;
				  mcc[1][3] = -cc_d;

				  
				  if (in_freq[k] > 1400 && in_freq[k] < 1401)
				    {
				      printf("mfeed: %g \n",in_freq[k]);
				      displayMatrix_4x4(mfeed);
				    }
				  */

				  
				  sdhdf_copy_mat4(marray,mueller);
				  sdhdf_mult4x4_replace(mueller,mcc);
				  sdhdf_mult4x4_replace(mueller,mfeed);
				  sdhdf_mult4x4_replace(mueller,mhand);
				  sdhdf_mult4x4_replace(mueller,mpa);
				  sdhdf_mult4x4_replace(mueller,mastro);
				  

				  if (muellerI == 1) // Do not do any calibration
				    sdhdf_setIdentity_4x4(mueller);
				  
				  sdhdf_inv4x4(mueller,mueller_inv);
				  
				  if (in_freq[k] > 1400 && in_freq[k] < 1401)
				    {
				      printf("Mueller:\n");
				      displayMatrix_4x4(mueller);
				      printf("Mueller^-1:\n");
				      displayMatrix_4x4(mueller_inv);
				    }

				  sdhdf_multMat_vec_replace(mueller_inv,stokesM);
				  if (in_freq[k] > 1400 && in_freq[k] < 1401)
				    {
				      printf("\n\nStokes after inversion\n");
				      sdhdf_display_vec4(stokesM);
				    }
				  
				  aI = stokesM[0];
				  aQ = stokesM[1];
				  aU = stokesM[2];
				  aV = stokesM[3];
				  
				  if (in_freq[k] > 1400 && in_freq[k] < 1401)
				    {
				      printf("(aI,aQ,aU,aV) = (%g,%g,%g,%g)\n",aI,aQ,aU,aV);
				    }


				  //				  *= nchan/c_nchan/gainVal2/calBin; ///nsum;				  

				  /*
				  out_data[k]             = aI/ge;
				  out_data[k+out_nchan]   = aQ/ge;
				  out_data[k+2*out_nchan] = aU/ge;
				  out_data[k+3*out_nchan] = aV/ge;
				  */
				  //				  printf("nchan, c_nchan = %d %d\n",nchan,c_nchan);
				  // 262144, 128
				  //				  out_Tdata[k+nchan] *= nchan/c_nchan/gainVal2/calBin; ///nsum;				  


				  //  
				  // GEORGE: SHOULD GET THE SCALING VIA CHANNEL BANDWIDTHS IN THE ASTRO AND CAL CHANNEL
				  //
				  if (origNchan > -1)
				    {
				      //	      printf("Using the original Nchan %d instead of %d\n",origNchan,nchan);
				      out_data[k]             = aI/ge/calBin*origNchan/c_nchan;
				      out_data[k+out_nchan]   = aQ/ge/calBin*origNchan/c_nchan;
				      out_data[k+2*out_nchan] = aU/ge/calBin*origNchan/c_nchan;
				      out_data[k+3*out_nchan] = aV/ge/calBin*origNchan/c_nchan;
				    }
				  else
				    {				      
				      out_data[k]             = aI/ge/calBin*nchan/c_nchan;
				      out_data[k+out_nchan]   = aQ/ge/calBin*nchan/c_nchan;
				      out_data[k+2*out_nchan] = aU/ge/calBin*nchan/c_nchan;
				      out_data[k+3*out_nchan] = aV/ge/calBin*nchan/c_nchan;
				    }
				}
			    }
			}
		    
		      // update frequencies if required
		      if (bary==1 || lsr==1)
			{
			  double vOverC;
			  double iX1,iY1,iX2,iY2,m,c;
			  double iX,iY;
			  for (j=0;j<out_ndump;j++)
			    {
			      vOverC=sdhdf_calcVoverC(inFile,b,ii,j,eop,nEOP,lsr,ephemName);
			      if (regrid==0)
				{
				  printf("Changing frequency axis => not regridding. V/c = %g. Out_nchan = %d\n",vOverC,out_nchan);				  
				  for (k=0;k<out_nchan;k++)
				    out_freq[k]*=(1.0-vOverC);
				}
			      else
				{
				  float * temp_data = (float *)calloc(sizeof(float),out_nchan*out_npol);
				  double freqNew,freqOld;
				  int kk;
				  int deltaI;
				  double fracDeltaI;
				  double df;
				
				  printf("Re-gridding: out_npol = %d out_nchan = %d out_ndump = %d\n",out_npol,out_nchan,out_ndump);
				  memcpy(temp_data,out_data+j*out_nchan*out_npol,sizeof(float)*out_nchan*out_npol);
				  for (kk=0;kk<out_npol;kk++)
				    {
				      for (k=0;k<out_nchan;k++)
					out_data[j*out_nchan*out_npol + kk*out_nchan + k] = 0.0;
				    }
			       
				  // Should do a proper gridding or a memcpy here
				  for (k=0;k<out_nchan;k++)
				    {
				      freqNew = out_freq[k]*(1.0-vOverC);
				      freqOld = out_freq[k];
				      // Should first check if we're still in topocentric frequencies -- DO THIS ***
				      df = ((double)(freqOld-freqNew)/(double)(out_freq[1]-out_freq[0])); // Check if frequency channelisation changes - e.g., at subband boundaries ** FIX ME
				      
				      if (df > 0)
					{
					  deltaI = (int)(df+0.5);
					  fracDeltaI = deltaI - df;
					}
				      else
					{
					  deltaI = -(int)(fabs(df)+0.5);
					  fracDeltaI = deltaI - df;  // CHECK MINUS SIGN
					}
				    
				      //				      printf("Here with %.6f %.6f %g %d %g %d %d\n",freqOld,freqNew,df,deltaI,fracDeltaI,k,out_nchan);
				      //				      deltaI  = 0;
				      if (k + deltaI >= 0 && k + deltaI < out_nchan)
					{
					  for (kk=0;kk<out_npol;kk++)
					    {
					      // ** Do an interpolation here **
					      iX1 = 0; iX2 = 1;
					      iY1 = temp_data[kk*out_nchan + k + deltaI];
					      iY2 = temp_data[kk*out_nchan + k + deltaI + 1]; // SHOULD CHECK IF AT EDGE
					      m   = (iY2-iY1)/(iX2-iX1);
					      c   = iY1;
					      iX  = fracDeltaI;
					      iY  = m*iX+c;

					      // out_data[j*out_nchan*out_npol + kk*out_nchan + k] = temp_data[kk*out_nchan + k + deltaI];
					      out_data[j*out_nchan*out_npol + kk*out_nchan + k] = iY;

					      // GEORGE HACK
					      //					      out_data[j*out_nchan*out_npol + kk*out_nchan + k] = iY1;
					    }
					}
				    }
				  free(temp_data);
				}
			    }
			}
		      
		      
		      sdhdf_copyAttributes(inFile->beam[b].bandData[ii].astro_obsHeaderAttr,inFile->beam[b].bandData[ii].nAstro_obsHeaderAttributes,dataAttributes,&nDataAttributes);
		      sdhdf_copyAttributes(inFile->beam[b].bandData[ii].astro_obsHeaderAttr_freq,inFile->beam[b].bandData[ii].nAstro_obsHeaderAttributes_freq,freqAttributes,&nFreqAttributes);

		    
		      if (bary==1 || lsr == 1)
			{
			  int it;
			  for (it=0;it<nFreqAttributes;it++)
			    {
			      if (strcmp(freqAttributes[it].key,"FRAME")==0)
				{
				  if (bary==1)
				    strcpy(freqAttributes[it].value,"barycentric");
				  else if (lsr==1)
				    strcpy(freqAttributes[it].value,"LSR");				    
				}
			    }
			}
		   
		      sdhdf_writeSpectrumData(outFile,inFile->beam[b].bandHeader[ii].label,b,ii,out_data,out_freq,out_nchan,out_npol,out_ndump,0,dataAttributes,nDataAttributes,freqAttributes,nFreqAttributes);
		      // Write out the obs_params file
		      sdhdf_writeObsParams(outFile,inFile->beam[b].bandHeader[ii].label,b,ii,outObsParams,out_ndump,1);
		      free(outObsParams);
		      
		      printf("Writing out the data weights for band %d\n",ii);
		      sdhdf_writeDataWeights(outFile,b,ii,dataWts,out_nchan,out_ndump,inFile->beam[b].bandHeader[ii].label);
		      printf("Complete writing data weights\n");
		   
		    

		      // FIX ME
		      //		      sdhdf_writeFrequencyAttributes(outFile,inFile->bandHeader[ii].label);
		      //		      sdhdf_writeDataAttributes(outFile,inFile->bandHeader[ii].label);
		      
		      
		      inBandParams[ii].nchan  = out_nchan;
		      inBandParams[ii].npol   = out_npol;
		      inBandParams[ii].f0     = out_freq[0];
		      inBandParams[ii].f1     = out_freq[out_nchan-1];
		      inBandParams[ii].ndump  = out_ndump;
		      if (tScrunch==1)
			inBandParams[ii].dtime = tav;
		      
		      free(in_data);
		      free(out_data);
		      free(out_Tdata);
		      free(out_Fdata);
		      free(in_freq);
		      free(out_freq);
		      free(dataWts);
		      // NOW NEED TO UPDATE META-DATA
		      // Subintegration change
		      //  band_header
		      //  obs_params
		      // Polarisation change
		      //  
		      
		    }
		  sdhdf_writeBandHeader(outFile,inBandParams,b,inFile->beam[b].nBand,1);
		  free(inBandParams);

		  if (nScal > 0)
		    {
		      free(cal_p1);
		      free(cal_p2);
		      free(cal_freq);
		    }

		}
	      // Copy other primary tables
	      
	      // FIX ME -- NEED TO INCREASE SIZE OF HISTORY ARRAY
	      //	      sdhdf_addHistory(inFile->history,inFile->nHistory,"sdhdf_modify","sdhdfProc software to modify a file",args);	     
	      //	      inFile->nHistory++;
	      sdhdf_writeHistory(outFile,inFile->history,inFile->nHistory);

	   

	      
	      sdhdf_copyRemainder(inFile,outFile,0);
	      
	      /*
		FIX ME
		sdhdf_copyEntireGroup("obs_config",inFile,outFile);	      	      
		sdhdf_copyEntireGroup("cal_band_header",inFile,outFile);
		sdhdf_copyEntireGroup("cal_obs_params",inFile,outFile);
	      */
	      // SHOULD BE UPDATING THE HEADER INFORMATION
	      
	      
	      // FIX ME
	      /*
		sdhdf_writeSpectralDumpHeader(outFile,inFile->spectralDumpHeader,out_ndump);
		sdhdf_copyEntireGroup("primary_header",inFile,outFile);
		sdhdf_copyEntireGroup("software_versions",inFile,outFile);
	      */
	      sdhdf_closeFile(inFile);
	      sdhdf_closeFile(outFile);
	      printf("Finished file %s %s\n",fname,oname);
	    }
	}
    }
  free(scalFreq); free(scalAA); free(scalBB);
  free(inFile);
  free(outFile);
  free(eop);
  
}


double linearInterpolate(double *freq, double *scaleAA,double in_freq,int kk,int nScale)
{
  double y1,y2;
  double x1,x2;
  double m,c;
  if (kk==0)
    {
      y1 = scaleAA[kk];
      y2 = scaleAA[kk+1];
      x1 = freq[kk];
      x2 = freq[kk+1];
    }
  else 
    {
      y1 = scaleAA[kk-1];
      y2 = scaleAA[kk];
      x1 = freq[kk-1];
      x2 = freq[kk];
    }
  m = (y2-y1)/(x2-x1);
  c = y1-m*x1;
  return m*in_freq+c;
}


void getScal(double sysGain_freq,double *scalFreq,double *scalAA,double *scalBB,int nScal,double *scalAA_val,double *scalBB_val)
{
  int i;
  double m,c,y1,y2,x1,x2;
  int set=0;
  *scalAA_val = *scalBB_val = 0.0;
  
  for (i=0;i<nScal-1;i++)
    {
      if (sysGain_freq <= scalFreq[i])
	{
	  if (i!=0)
	    {
	      x1 = scalFreq[i-1]; x2 = scalFreq[i]; y1 = scalAA[i-1]; y2 = scalAA[i];
	      m = (y2-y1)/(x2-x1); c = y1-m*x1;
	      *scalAA_val = m*sysGain_freq+c; // Should do an interpolation
	      
	      x1 = scalFreq[i-1]; x2 = scalFreq[i]; y1 = scalBB[i-1]; y2 = scalBB[i];
	      m = (y2-y1)/(x2-x1); c = y1-m*x1;
	      *scalBB_val = m*sysGain_freq+c;
	    }
	  else
	    {
	      x1 = scalFreq[i+1]; x2 = scalFreq[i]; y1 = scalAA[i+1]; y2 = scalAA[i];
	      m = (y2-y1)/(x2-x1); c = y1-m*x1;
	      *scalAA_val = m*sysGain_freq+c; // Should do an interpolation
	      
	      x1 = scalFreq[i+1]; x2 = scalFreq[i]; y1 = scalBB[i+1]; y2 = scalBB[i];
	      m = (y2-y1)/(x2-x1); c = y1-m*x1;
	      *scalBB_val = m*sysGain_freq+c;
	    }
	  set=1;
	      
	  break;
	}
    }
  if (set==0)
    {
      x1 = scalFreq[nScal-2]; x2 = scalFreq[nScal-1]; y1 = scalAA[nScal-2]; y2 = scalAA[nScal-1];
      m = (y2-y1)/(x2-x1); c = y1-m*x1;
      *scalAA_val = m*sysGain_freq+c; // Should do an interpolation
      
      x1 = scalFreq[nScal-2]; x2 = scalFreq[nScal-1]; y1 = scalBB[nScal-2]; y2 = scalBB[nScal-1];
      m = (y2-y1)/(x2-x1); c = y1-m*x1;
      *scalBB_val = m*sysGain_freq+c;
    }
}

/* ************************************************************** *
 * Spline routines                                                *
 * ************************************************************** */
#define ABS(x) ((x) < 0 ? -(x) : (x))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))


void TKcmonot (int n, double x[], double y[], double yd[][4])
/*****************************************************************************
compute cubic interpolation coefficients via the Fritsch-Carlson method,
which preserves monotonicity
******************************************************************************
Input:
n		number of samples
x  		array[n] of monotonically increasing or decreasing abscissae
y		array[n] of ordinates

Output:
yd		array[n][4] of cubic interpolation coefficients (see notes)
******************************************************************************
Notes:
The computed cubic spline coefficients are as follows:
yd[i][0] = y(x[i])    (the value of y at x = x[i])
yd[i][1] = y'(x[i])   (the 1st derivative of y at x = x[i])
yd[i][2] = y''(x[i])  (the 2nd derivative of y at x = x[i])
yd[i][3] = y'''(x[i]) (the 3rd derivative of y at x = x[i])

To evaluate y(x) for x between x[i] and x[i+1] and h = x-x[i],
use the computed coefficients as follows:
y(x) = yd[i][0]+h*(yd[i][1]+h*(yd[i][2]/2.0+h*yd[i][3]/6.0))

The Fritsch-Carlson method yields continuous 1st derivatives, but 2nd
and 3rd derivatives are discontinuous.  The method will yield a
monotonic interpolant for monotonic data.  1st derivatives are set
to zero wherever first divided differences change sign.

For more information, see Fritsch, F. N., and Carlson, R. E., 1980,
Monotone piecewise cubic interpolation:  SIAM J. Numer. Anal., v. 17,
n. 2, p. 238-246.

Also, see the book by Kahaner, D., Moler, C., and Nash, S., 1989,
Numerical Methods and Software, Prentice Hall.  This function was
derived from SUBROUTINE PCHEZ contained on the diskette that comes
with the book.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 09/30/89
Modified:  Dave Hale, Colorado School of Mines, 02/28/91
	changed to work for n=1.
Modified:  Dave Hale, Colorado School of Mines, 08/04/91
	fixed bug in computation of left end derivative
*****************************************************************************/
{
  int i;
  double h1,h2,del1,del2,dmin,dmax,hsum,hsum3,w1,w2,drat1,drat2,divdf3;
  
  /* copy ordinates into output array */
  for (i=0; i<n; i++)
    yd[i][0] = y[i];
  
  /* if n=1, then use constant interpolation */
  if (n==1) {
    yd[0][1] = 0.0;
    yd[0][2] = 0.0;
    yd[0][3] = 0.0;
    return;
    
    /* else, if n=2, then use linear interpolation */
  } else if (n==2) {
    yd[0][1] = yd[1][1] = (y[1]-y[0])/(x[1]-x[0]);
    yd[0][2] = yd[1][2] = 0.0;
    yd[0][3] = yd[1][3] = 0.0;
    return;
  }
  
  /* set left end derivative via shape-preserving 3-point formula */
  h1 = x[1]-x[0];
  h2 = x[2]-x[1];
  hsum = h1+h2;
  del1 = (y[1]-y[0])/h1;
  del2 = (y[2]-y[1])/h2;
  w1 = (h1+hsum)/hsum;
  w2 = -h1/hsum;
  yd[0][1] = w1*del1+w2*del2;
  if (yd[0][1]*del1<=0.0)
    yd[0][1] = 0.0;
  else if (del1*del2<0.0) {
    dmax = 3.0*del1;
    if (ABS(yd[0][1])>ABS(dmax)) yd[0][1] = dmax;
  }
  
  /* loop over interior points */
  for (i=1; i<n-1; i++) {
    
    /* compute intervals and slopes */
    h1 = x[i]-x[i-1];
    h2 = x[i+1]-x[i];
    hsum = h1+h2;
    del1 = (y[i]-y[i-1])/h1;
    del2 = (y[i+1]-y[i])/h2;
    
    /* if not strictly monotonic, zero derivative */
    if (del1*del2<=0.0) {
      yd[i][1] = 0.0;
      
      /*
       * else, if strictly monotonic, use Butland's formula:
       *      3*(h1+h2)*del1*del2
       * -------------------------------
       * ((2*h1+h2)*del1+(h1+2*h2)*del2)
       * computed as follows to avoid roundoff error
       */
    } else {
      hsum3 = hsum+hsum+hsum;
      w1 = (hsum+h1)/hsum3;
      w2 = (hsum+h2)/hsum3;
      dmin = MIN(ABS(del1),ABS(del2));
      dmax = MAX(ABS(del1),ABS(del2));
      drat1 = del1/dmax;
      drat2 = del2/dmax;
      yd[i][1] = dmin/(w1*drat1+w2*drat2);
    }
  }
  
  /* set right end derivative via shape-preserving 3-point formula */
  w1 = -h2/hsum;
  w2 = (h2+hsum)/hsum;
  yd[n-1][1] = w1*del1+w2*del2;
  if (yd[n-1][1]*del2<=0.0)
    yd[n-1][1] = 0.0;
  else if (del1*del2<0.0) {
    dmax = 3.0*del2;
    if (ABS(yd[n-1][1])>ABS(dmax)) yd[n-1][1] = dmax;
  }
  
  /* compute 2nd and 3rd derivatives of cubic polynomials */
  for (i=0; i<n-1; i++) {
    h2 = x[i+1]-x[i];
    del2 = (y[i+1]-y[i])/h2;
    divdf3 = yd[i][1]+yd[i+1][1]-2.0*del2;
    yd[i][2] = 2.0*(del2-yd[i][1]-divdf3)/h2;
    yd[i][3] = (divdf3/h2)*(6.0/h2);
  }
  yd[n-1][2] = yd[n-2][2]+(x[n-1]-x[n-2])*yd[n-2][3];
  yd[n-1][3] = yd[n-2][3];
}

// Interpolate the spline fit on to a given set of x-values
void TKspline_interpolate(int n,double *x,double *y,double yd[][4],double *interpX,
			  double *interpY,int nInterp)
{
  double h;
  int jpos;
  int i,j;

  //To evaluate y(x) for x between x[i] and x[i+1] and h = x-x[i],
  //use the computed coefficients as follows:
  //y(x) = yd[i][0]+h*(yd[i][1]+h*(yd[i][2]/2.0+h*yd[i][3]/6.0))
  
  // Assume the data are sorted so x[i] < x[i+1]
  for (i=0;i<nInterp;i++)
    {
      jpos=-1;
      for (j=0;j<n-1;j++)
	{
	  if (interpX[i]>=x[j] && interpX[i]<x[j+1])
	    {
	      jpos=j;
	      break;
	    }
	}
      if (jpos!=-1)
	{
	  h = interpX[i]-x[jpos];
	  interpY[i] = yd[jpos][0]+h*(yd[jpos][1]+h*(yd[jpos][2]/2.0+h*yd[jpos][3]/6.0));
	}
      else
	interpY[i] = 0.0;
    }

}
