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
void getScal(double sysGain_freq,double *scalFreq,double *scalAA,double *scalBB,int nScal,double *scalAA_val,double *scalBB_val);

void TKspline_interpolate(int n,double *x,double *y,double yd[][4],double *interpX,
			  double *interpY,int nInterp);
void TKcmonot (int n, double x[], double y[], double yd[][4]);


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
  int nBeam=0;
  int nBand=0;
  int tScrunch=0;
  int fScrunch=0;
  int pol=0;
  int stokes=0;
  int fAv=0;
  int bary=0,lsr=0;
  int regrid=0;
  int nh;
  float *in_data,*out_data,*out_freq;
  double *in_freq;
  float *out_Tdata,*out_Fdata;
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
  double freq[128*26],scaleAA[128*26],scaleBB[128*26];

  int nScale=0;
  double *scalFreq,*scalAA,*scalBB;
  float *cal_p1,*cal_p2,*cal_freq;

  double sysGain_freq[4096],sysGain_p1[4096],sysGain_p2[4096],scalAA_val,scalBB_val;
  double lastGood_p1=1e8,lastGood_p2=1e8;
  char scalFname[1024]="NULL";
  int  nScal=-1;
  int c_nchan;
  int nc;
  int nRF1,nRF2,nRF3;
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
  
  strcpy(oname,"sdhdf_modify_output.hdf");

  eop = (sdhdf_eopStruct *)malloc(sizeof(sdhdf_eopStruct)*MAX_EOP_LINES);
  
  scalFreq = (double *)malloc(sizeof(double)*3330);
  scalAA = (double *)malloc(sizeof(double)*3330);
  scalBB = (double *)malloc(sizeof(double)*3330);


  
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
      else if (strcmp(argv[i],"-scal")==0)
	strcpy(scalFname,argv[++i]);
      else if (strcmp(argv[i],"-interp")==0) // 1 = linear, 2 = cubic spline
	sscanf(argv[++i],"%d",&interp);
      else if (strcmp(argv[i],"-regrid")==0)
	regrid=1;
      else if (strcmp(argv[i],"-T")==0)
	{sdhdf_add1arg(args,argv[i]); tScrunch=1;}
      else if (strcmp(argv[i],"-F")==0)
	{sdhdf_add1arg(args,argv[i]); fScrunch=1;}
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
      else
	{
	  strcpy(fname,argv[i]);
	  printf("Here with %s\n",fname);
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
		    nScale++;
		}
	      fclose(fScale);
	    }

	  if (strcmp(scalFname,"NULL")!=0)
	    {
	      FILE *fin;
	      
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
		    nScal++;
		}
	      fclose(fin);
	      printf("Complete loading %s\n",scalFname);
	    }
	  

	  sdhdf_initialiseFile(inFile);
	  sdhdf_initialiseFile(outFile);

	  printf("Opening in file\n");
	  sdhdf_openFile(fname,inFile,1);
	  //	  inFile->fileID  = H5Fopen(fname,H5F_ACC_RDONLY,H5P_DEFAULT);
	  printf("Complete opening file with %d\n",inFile->fileID);
	  sprintf(oname,"%s.%s",fname,extension);
	  printf("Opening output file >%s<\n",oname);
	  sdhdf_openFile(oname,outFile,3);
	  printf("Complete open\n");
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
		  
		  inBandParams = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*inFile->beam[b].nBand);      
		  sdhdf_copyBandHeaderStruct(inFile->beam[b].bandHeader,inBandParams,inFile->beam[b].nBand);

		  // Do we need to extract the calibrator information?
		  if (nScal > 0)
		    {
		      int ci,cj,ck;
		      double onP1,offP1,onP2,offP2;
		      double freq;
		      
		      printf("Calibrating\n");
		      // Could calibrate dump by dump, or average the entire data file
		      // FOR NOW **** just averaging all the calibration information in the file
		      cal_p1 = (float *)malloc(sizeof(float)*128*26); // Should set properly
		      cal_p2 = (float *)malloc(sizeof(float)*128*26);
		      cal_freq = (float *)malloc(sizeof(float)*128*26);
		      nc=0;

		      nRF1=nRF2=nRF3=0;
		      
		      
		      for (ci=0;ci<inFile->beam[b].nBand;ci++)
			{
			  sdhdf_loadBandData(inFile,b,ci,2);
			  sdhdf_loadBandData(inFile,b,ci,3);
			  
			  
			  c_nchan = inFile->beam[b].calBandHeader[ci].nchan;
			  printf("Number of channels = %d\n",c_nchan);
			  for (cj=0;cj<c_nchan;cj++)
			    {
			      onP1=offP1=onP2=offP2=0.0;
			      //			      printf("ndump = %d\n",inFile->beam[b].calBandHeader[ii].ndump);
			      for (ck=0;ck<inFile->beam[b].calBandHeader[ci].ndump;ck++)
				{
				  onP1  += inFile->beam[b].bandData[ci].cal_on_data.pol1[cj+ck*c_nchan];
				  offP1 += inFile->beam[b].bandData[ci].cal_off_data.pol1[cj+ck*c_nchan];
				  onP2  += inFile->beam[b].bandData[ci].cal_on_data.pol2[cj+ck*c_nchan];
				  offP2 += inFile->beam[b].bandData[ci].cal_off_data.pol2[cj+ck*c_nchan];
				}
			      // Note: summing, not averaging above ***
			      //
			      freq = cal_freq[nc] = inFile->beam[b].bandData[ci].cal_on_data.freq[cj];
			      cal_p1[nc] = (onP1-offP1);  
			      cal_p2[nc] = (onP2-offP2);
			      //			      printf("freq = %g\n",freq);

			      if ((freq > 814 && freq < 818) ||
				  (freq > 994 && freq < 998) ||
				  (freq > 1068 && freq < 1073) ||
				  (freq > 1194 && freq < 1198) ||
				  (freq > 1295 && freq < 1305))
				{
				  scA_rf1_x[nRF1] = scB_rf1_x[nRF1] = log10(freq);
				  scA_rf1_y[nRF1] = log10(cal_p1[nc]);
				  scB_rf1_y[nRF1] = log10(cal_p2[nc]);
				  nRF1++;
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
				  scA_rf2_x[nRF2] = scB_rf2_x[nRF2] = log10(freq);
				  scA_rf2_y[nRF2] = log10(cal_p1[nc]);
				  scB_rf2_y[nRF2] = log10(cal_p2[nc]);
				  nRF2++;
				  
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
				  scA_rf3_x[nRF3] = scB_rf3_x[nRF3] = log10(freq);
				  scA_rf3_y[nRF3] = log10(cal_p1[nc]);
				  scB_rf3_y[nRF3] = log10(cal_p2[nc]);
				  nRF3++;
				  //  scA_rf3+=cal_p1[nc]; scB_rf3+=cal_p2[nc]; nRF3++;
				}
			      
			      
			      nc++;
			    }
			}
		      // Fit straight lines to the cal data
		      TKleastSquares_svd_noErr(scA_rf1_x,scA_rf1_y,nRF1,scA_rf1,3,TKfitPoly);
		      TKleastSquares_svd_noErr(scB_rf1_x,scB_rf1_y,nRF1,scB_rf1,3,TKfitPoly);
		      
		      TKleastSquares_svd_noErr(scA_rf2_x,scA_rf2_y,nRF2,scA_rf2,3,TKfitPoly);
		      TKleastSquares_svd_noErr(scB_rf2_x,scB_rf2_y,nRF2,scB_rf2,3,TKfitPoly);

		      TKleastSquares_svd_noErr(scA_rf3_x,scA_rf3_y,nRF3,scA_rf3,3,TKfitPoly);
		      TKleastSquares_svd_noErr(scB_rf3_x,scB_rf3_y,nRF3,scB_rf3,3,TKfitPoly);

		      {
			FILE *fout;
			fout = fopen("calOutput.dat","w");
			for (i=0;i<nRF3;i++)
			  fprintf(fout,"calOutput: %g %g %g %g %g %g %g %g\n",scA_rf3_x[i],scA_rf3_y[i],scB_rf3_x[i],scB_rf3_y[i],scA_rf3[0],scB_rf3[0],scA_rf3[1],scB_rf3[1]);
			fclose(fout);
		      }

		      
		      printf("Loaded cal %g %g %g %g %g %g %g %g %g %g %g %g %d %d %d\n",scA_rf1[0],scB_rf1[0],scA_rf1[1],scB_rf1[1],
			     scA_rf2[0],scB_rf2[0],scA_rf2[1],scB_rf2[1],scA_rf3[0],scB_rf3[0],scA_rf3[1],scB_rf3[1],nRF1,nRF2,nRF3);

		    }
	       		
		  
		  // Copy the bands
		  for (ii=0;ii<inFile->beam[b].nBand;ii++)
		    {
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
		      
		      if (fScrunch==1)
			fAv = nchan;
		      
		      if (fAv < 2)
			out_nchan = nchan;
		      else
			out_nchan = (int)(nchan/fAv);
		      
		      if (pol == 1)
			out_npol=1;
		      else
			out_npol  = npol;
		      
		      if (tScrunch==1)
			out_ndump = 1;
		      else
			out_ndump = inFile->beam[b].bandHeader[ii].ndump;


		      if (nScal > 0)
			{
			  int cj,ck;
			  double onP1,offP1,onP2,offP2;
			  double freq;
			  
			  printf("Calibrating\n");
			  // Could calibrate dump by dump, or average the entire data file
			  // FOR NOW **** just averaging all the calibration information in the file
			  sdhdf_loadBandData(inFile,b,ii,2);
			  //			  printf("Loaded cal 1\n");
			  sdhdf_loadBandData(inFile,b,ii,3);
			  //			  printf("Loaded cal 2\n");
			  
			  c_nchan = inFile->beam[b].calBandHeader[ii].nchan;
			  printf("Number of channels = %d\n",c_nchan);
			  for (cj=0;cj<c_nchan;cj++)
			    {
			      //			      onP1=offP1=onP2=offP2=0.0;
			      //			      printf("ndump = %d\n",inFile->beam[b].calBandHeader[ii].ndump);
			      //			      for (ck=0;ck<inFile->beam[b].calBandHeader[ii].ndump;ck++)
			      //				{
			      //				  onP1  += inFile->beam[b].bandData[ii].cal_on_data.pol1[cj+ck*c_nchan];
			      //				  offP1 += inFile->beam[b].bandData[ii].cal_off_data.pol1[cj+ck*c_nchan];
			      //				  onP2  += inFile->beam[b].bandData[ii].cal_on_data.pol2[cj+ck*c_nchan];
			      //				  offP2 += inFile->beam[b].bandData[ii].cal_off_data.pol2[cj+ck*c_nchan];
			      //				}
			      // Note: summing, not averaging above ***
			      //
			      sysGain_freq[cj] = inFile->beam[b].bandData[ii].cal_on_data.freq[cj];
			      getScal(sysGain_freq[cj],scalFreq,scalAA,scalBB,nScal,&scalAA_val,&scalBB_val);
			      //			      sysGain_p1[cj] = (onP1-offP1);  // NOTE /scalAA_val for true gain;
			      //			      sysGain_p2[cj] = (onP2-offP2);

						  
			  
			      // GEORGE TESTING
			      //sysGain_p1[cj] = (onP1-offP1)/scalAA_val;
			      //sysGain_p2[cj] = (onP2-offP2)/scalBB_val;
			      //			      scA_rf1 /= nRF1;
			      //			      scB_rf1 /= nRF1;

			      scA_rf1_val = pow(10,scA_rf1[0] + scA_rf1[1]*log10(sysGain_freq[cj]) + scA_rf1[2]*pow(log10(sysGain_freq[cj]),2));
			      scB_rf1_val = pow(10,scB_rf1[0] + scB_rf1[1]*log10(sysGain_freq[cj]) + scB_rf1[2]*pow(log10(sysGain_freq[cj]),2));
			      scA_rf2_val = pow(10,scA_rf2[0] + scA_rf2[1]*log10(sysGain_freq[cj]) + scA_rf2[2]*pow(log10(sysGain_freq[cj]),2));
			      scB_rf2_val = pow(10,scB_rf2[0] + scB_rf2[1]*log10(sysGain_freq[cj]) + scB_rf2[2]*pow(log10(sysGain_freq[cj]),2));
			      scA_rf3_val = pow(10,scA_rf3[0] + scA_rf3[1]*log10(sysGain_freq[cj]) + scA_rf3[2]*pow(log10(sysGain_freq[cj]),2));
			      scB_rf3_val = pow(10,scB_rf3[0] + scB_rf3[1]*log10(sysGain_freq[cj]) + scB_rf3[2]*pow(log10(sysGain_freq[cj]),2));
			      
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
			      
			      /*			      
			      if (sysGain_p1[cj] < 1e6) // THIS IS A MADE-UP NUMBER -- GEORGE FIX
				{sysGain_p1[cj] = lastGood_p1; printf("WARNING 1\n");}
			      if (sysGain_p2[cj] < 1e6) // THIS IS A MADE-UP NUMBER -- GEORGE FIX
				{sysGain_p2[cj] = lastGood_p2; printf("WARNING 2\n");}
			      if (sysGain_p1[cj] > 1e12) // THIS IS A MADE-UP NUMBER -- GEORGE FIX
				{sysGain_p1[cj] = lastGood_p1; printf("WARNING 3\n");}
			      if (sysGain_p2[cj] > 1e12) // THIS IS A MADE-UP NUMBER -- GEORGE FIX
				{sysGain_p2[cj] = lastGood_p2; printf("WARNING 4\n");}
			      lastGood_p1 = sysGain_p1[cj];
			      lastGood_p2 = sysGain_p2[cj];
			      */
			      printf("Gain: %.6f %g %g %g %g %g %g %g %g %g %g %g %g\n",sysGain_freq[cj],sysGain_p1[cj],sysGain_p2[cj],scalAA_val,scalBB_val,
				     scA_rf1_val,scB_rf1_val, scA_rf2_val,scB_rf2_val,scA_rf3_val,scB_rf3_val,scA_rf1[0],scA_rf1[1]);
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
		      
		      in_data = (float *)malloc(sizeof(float)*nchan*npol*inFile->beam[b].bandHeader[ii].ndump);
		      in_freq = (double *)malloc(sizeof(double)*nchan);
		      out_freq = (float *)malloc(sizeof(float)*out_nchan);
		      out_data = (float *)calloc(sizeof(float),out_nchan*out_npol*out_ndump);
		      out_Tdata = (float *)calloc(sizeof(float),nchan*npol*out_ndump);
		      out_Fdata = (float *)calloc(sizeof(float),out_nchan*npol*out_ndump);
		      
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
				  // FIX ME -- ****
				  sclAA *= (nchan/128);
				  sclBB *= (nchan/128);
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
				  for (k=0;k<nchan;k++)
				    {
				      // NOT ACCOUNTING FOR WEIGHTINGS			  
				      if (npol==1)
					out_Tdata[k] += inFile->beam[b].bandData[ii].astro_data.pol1[k+j*nchan]; ///(float)nsd;
				      else if (npol==2)
					{
					  out_Tdata[k]         += inFile->beam[b].bandData[ii].astro_data.pol1[k+j*nchan]; // /(float)nsd;
					  out_Tdata[k+nchan]   += inFile->beam[b].bandData[ii].astro_data.pol2[k+j*nchan]; // /(float)nsd;
					}
				      else
					{
					  out_Tdata[k]         += inFile->beam[b].bandData[ii].astro_data.pol1[k+j*nchan]; // /(float)nsd;
					  out_Tdata[k+nchan]   += inFile->beam[b].bandData[ii].astro_data.pol2[k+j*nchan]; // /(float)nsd;
					  out_Tdata[k+2*nchan] += inFile->beam[b].bandData[ii].astro_data.pol3[k+j*nchan]; // /(float)nsd;
					  out_Tdata[k+3*nchan] += inFile->beam[b].bandData[ii].astro_data.pol4[k+j*nchan]; // /(float)nsd;  
					}
				      
				    }
				}
			    }
			  // Scale the output
			  // GEORGE **** THIS SHOULD BE EARLIER AS YOU MAY NOT CHOOSE TO TIME SCRUNCH *****
			  if (nScal > 0)
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
				  printf("Scale factors: %.6f %g %g\n",in_freq[k],gainVal1,gainVal2);
				  // Only scaling 2 polarisations *** <<<
				  
				  if (gainVal1 <= 0) gainVal1=1e9; // SHOULD SET MORE SENSIBLY ****
				  if (gainVal2 <= 0) gainVal2=1e9; // SHOULD SET MORE SENSIBLY ****
				  
				  out_Tdata[k]         *= nchan/c_nchan/gainVal1/nsum; // Note scaling by number of spectral dumps as each one now should be in Jy
				  out_Tdata[k+nchan]   *= nchan/c_nchan/gainVal2/nsum;				  
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
				    out_Tdata[k+j*nchan*npol] = inFile->beam[b].bandData[ii].astro_data.pol1[k+j*nchan]/(float)nsd;
				  else if (npol==2)
				    {
				      out_Tdata[k+j*nchan*npol]         = inFile->beam[b].bandData[ii].astro_data.pol1[k+j*nchan]/(float)nsd;
				      out_Tdata[k+j*nchan*npol+nchan]   = inFile->beam[b].bandData[ii].astro_data.pol2[k+j*nchan]/(float)nsd;
				    }
				  else if (npol==4)
				    {
				      out_Tdata[k+j*nchan*npol]         = inFile->beam[b].bandData[ii].astro_data.pol1[k+j*nchan]/(float)nsd;
				      out_Tdata[k+j*nchan*npol+nchan]   = inFile->beam[b].bandData[ii].astro_data.pol2[k+j*nchan]/(float)nsd;
				      out_Tdata[k+j*nchan*npol+2*nchan] = inFile->beam[b].bandData[ii].astro_data.pol3[k+j*nchan]/(float)nsd;
				      out_Tdata[k+j*nchan*npol+3*nchan] = inFile->beam[b].bandData[ii].astro_data.pol4[k+j*nchan]/(float)nsd; 
				    }
				}	
			    }
			}
		      sdhdf_releaseBandData(inFile,b,ii,1);
		      // Frequency averaging
		      if (fAv >= 2)
			{
			  for (j=0;j<out_ndump;j++)
			    {
			      kp=0;
			      for (k=0;k<nchan;k+=fAv)
				{
				  av1 = av2 = av3 = av4 = avFreq = 0.0;
				  for (l=k;l<k+fAv;l++)
				    {
				      avFreq += in_freq[l];
				      if (npol==1)
					av1 += out_Tdata[l+j*npol*nchan];
				      else if (npol==2)
					{
					  av1 += out_Tdata[l+j*npol*nchan];
					  av2 += out_Tdata[l+j*npol*nchan+nchan];
					}
				      else
					{
					  av1 += out_Tdata[l+j*npol*nchan];
					  av2 += out_Tdata[l+j*npol*nchan+nchan];
					  av3 += out_Tdata[l+j*npol*nchan+2*nchan];
					  av4 += out_Tdata[l+j*npol*nchan+3*nchan];
					}
				    }
				  
				  if (j==0) out_freq[kp] = avFreq/fAv;
				  
				  if (npol==1)
				    out_Fdata[kp+out_nchan*j*npol]               = av1/fAv;
				  else if (npol==2)
				    {
				      out_Fdata[kp+out_nchan*j*npol]             = av1/fAv;
				      out_Fdata[kp+out_nchan*j*npol+out_nchan]   = av2/fAv;
				    }
				  else
				    {
				      out_Fdata[kp+out_nchan*j*npol]             = av1/fAv;
				      out_Fdata[kp+out_nchan*j*npol+out_nchan]   = av2/fAv;
				      out_Fdata[kp+out_nchan*j*npol+2*out_nchan] = av3/fAv;
				      out_Fdata[kp+out_nchan*j*npol+3*out_nchan] = av4/fAv;
				    }
				  kp++;
				  
				}
			    }
			}
		      else
			{
			  // Should use memcpy
			  printf("Not frequency scrunching\n");
			  for (j=0;j<out_ndump;j++)
			    {
			      for (k=0;k<out_nchan;k++)
				{
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
				out_data[k+out_nchan*j] = (out_Fdata[k+out_nchan*j*4]+out_Fdata[k+out_nchan*j*4+out_nchan]); 
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
			  for (j=0;j<out_ndump;j++)
			    {
			      float p1,p2,p3,p4;
			      for (k=0;k<out_nchan;k++)
				{
				  p1 = out_data[k];   // SHOULD BE A J IN HERE SOMEWHERE *** FIX
				  p2 = out_data[k+out_nchan];
				  p3 = out_data[k+2*out_nchan];
				  p4 = out_data[k+3*out_nchan];
				  sdhdf_convertStokes(p1,p2,p3,p4,&out_data[k],&out_data[k+out_nchan],&out_data[k+2*out_nchan],&out_data[k+3*out_nchan]);
				}
			    }
			}
		      
		      // update frequencies if required
		      if (bary==1 || lsr==1)
			{
			  double vOverC;
			  for (j=0;j<out_ndump;j++)
			    {
			      vOverC=sdhdf_calcVoverC(inFile,b,ii,j,eop,nEOP,lsr);
			      if (regrid==0)
				{
				  printf("Changing frequency axis\n");
				  for (k=0;k<out_nchan;k++)
				    out_freq[k]*=(1.0-vOverC);
				}
			      else
				{
				  float * temp_data = (float *)calloc(sizeof(float),out_nchan*out_npol);
				  double freqNew,freqOld;
				  int kk;
				  int deltaI;
				  
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
				      deltaI  = (int)((freqOld-freqNew)/(out_freq[1]-out_freq[0])+0.5); // Check if frequency channelisation changes - e.g., at subband boundaries ** FIX ME
				      //				      printf("Here with %.6f %.6f %d\n",freqOld,freqNew,deltaI);
				      //				      deltaI  = 0;
				      if (k + deltaI >= 0 && k + deltaI < out_nchan)
					{
					  for (kk=0;kk<out_npol;kk++)
					    out_data[j*out_nchan*out_npol + kk*out_nchan + k] = temp_data[kk*out_nchan + k + deltaI];
					}
				    }
				  free(temp_data);
				}
			    }
			}
		      
		      
		      // FIX ME
		      //			  if (bary==1)
		      //			    strcpy(outFile->frequency_attr.frame,"barycentric");
		      //			  else if (lsr == 1)
		      //			    strcpy(outFile->frequency_attr.frame,"LSR");
		      sdhdf_writeSpectrumData(outFile,inFile->beam[b].bandHeader[ii].label,b,ii,out_data,out_freq,out_nchan,out_npol,out_ndump,0);
		      // Write out the obs_params file
		      sdhdf_writeObsParams(outFile,inFile->beam[b].bandHeader[ii].label,b,ii,outObsParams,out_ndump);
		      free(outObsParams);
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
		      
		      // NOW NEED TO UPDATE META-DATA
		      // Subintegration change
		      //  band_header
		      //  obs_params
		      // Polarisation change
		      //  
		      
		    }
		  sdhdf_writeBandHeader(outFile,inBandParams,b,inFile->beam[b].nBand,1);
		  free(inBandParams);
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
  if (nScal > 0)
    {
      free(cal_p1);
      free(cal_p2);
      free(cal_freq);
    }
  
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
