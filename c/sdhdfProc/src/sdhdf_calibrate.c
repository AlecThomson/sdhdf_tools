//  Copyright (C) 2019, 2020, 2021, 2022 George Hobbs

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


// To fix:
//
// Should check if the user actually wants polarisation and flux calibration, or just 1
// Need a way to read the calibration information from an earlier file - or inject an earlier cal into a different file
// Need a way to specify which fluxcal/pcm file to use (and how the user can know what options exist)
// Currrently hardcode alpha = -45 degrees in sdhdfProc_calibration.c -- FIX THIS!
// Currently always time-averaging the noise source signal
// FIX ME: IMPORTANT: ISSUE IN THAT SHOULDN'T INTERPOLATE USING BAD ENTIRES (AS THEY ARE ZERO) (important for noise source stokes)
// Must work with older-style, non-normalised files.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sdhdfProc.h"

#define VERSION "v0.5"
#define MAX_POL_CAL_CHAN 4096    // FIX ME -- SHOULD SET DYNAMICALLY

void processFile(char *fname,char *oname, int stabiliseType, int out_npol, char *fluxCalFile, int verbose);

void help()
{
  printf("sdhdf_calibrate %s (SDHDFProc %s)\n",VERSION,SOFTWARE_VER);
  printf("Authors: G. Hobbs\n");
  printf("Purpose: to calibrate the data sets\n");
  printf("\n");
  printf("Command line arguments:\n\n");
  printf("-e <ext>    Output file extension\n");
  printf("-h          This help\n");
  printf("\n\n");
  printf("Example: sdhdf_calibrate -e cal uwl_*.hdf.T.f1024\n");
  exit(1);
}


int main(int argc,char *argv[])
{
  int i;
  int nFiles=0;
  char fname[MAX_FILES][MAX_STRLEN];
  char oname[MAX_STRLEN];
  char extension[MAX_STRLEN];
  int verbose=0;
  int stabiliseType=1;
  int out_npol = 4;
  char fluxCalFile[1024] = "uwl_220705_132107.fluxcal";
  
  strcpy(extension,"calibrate");
  
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-v")==0)
	verbose=1;    
      else if (strcmp(argv[i],"-V")==0)
	verbose=2;    
      else if (strcmp(argv[i],"-e")==0)
	strcpy(extension,argv[++i]);
      else if (strcmp(argv[i],"-fluxcal")==0)
	strcpy(fluxCalFile,argv[++i]);
      else if (strcmp(argv[i],"-s")==0)
	sscanf(argv[++i],"%d",&stabiliseType);
      else if (strcmp(argv[i],"-p")==0)
	sscanf(argv[++i],"%d",&out_npol);
      else
	strcpy(fname[nFiles++],argv[i]);      
    }
  printf("Processing: %d files\n",nFiles);
  for (i=0;i<nFiles;i++)
    {
      printf("Processing file: %s\n",fname[i]);
      sdhdf_formOutputFilename(fname[i],extension,oname);
      processFile(fname[i],oname,stabiliseType,out_npol,fluxCalFile,verbose);
    }

}

void processFile(char *fname,char *oname, int stabiliseType,int out_npol,char *fluxCalFile, int verbose)
{
  int ii,i,c,j,k,b;
  sdhdf_fileStruct *inFile,*outFile;
  int nbeam,npol,nchan,ndump,nchanFreq,load_nchanCal;
  int use_nchanCal;
  float *load_calFreq,*load_calAA,*load_calBB,*load_calRe,*load_calIm;
  float **interpCoeff_AA,**interpCoeff_BB,**interpCoeff_Re,**interpCoeff_Im;
  float **interpCoeff_fluxCalAA,**interpCoeff_fluxCalBB;
  float *out_data,*out_freq;
  float final_aa,final_bb,final_rab,final_iab;
  float stabilised_aa,stabilised_bb,stabilised_rab,stabilised_iab;
  float cal_aa,cal_bb,cal_rab,cal_iab;
  float measured_aa,measured_bb,measured_rab,measured_iab;
  float freq;
  sdhdf_attributes_struct dataAttributes[MAX_ATTRIBUTES];
  sdhdf_attributes_struct freqAttributes[MAX_ATTRIBUTES];
  int nDataAttributes=0;
  int nFreqAttributes=0;
  int nFluxCalChan=0;
  double stabilise_normFactor;
  float *fluxCalFreq,*fluxCalValAA,*fluxCalValBB;
  float fluxValAA,fluxValBB;
  sdhdf_fluxCalibration *fluxCal;

  double complex Jast[2][2];
  double complex rho[2][2];
  double paraAng;
  double complex R_feed_pa[2][2];
  double complex R_feed_pa_dag[2][2];
  double complex finalJ[2][2];
  double alpha = -45.0*M_PI/180.0; // FIX ME
  double fluxScale;

  FILE *debugOut1,*debugOut2;
  
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

  
  sdhdf_initialiseFile(inFile);
  sdhdf_initialiseFile(outFile);
  
  sdhdf_openFile(fname,inFile,1);
  sdhdf_openFile(oname,outFile,3);
  sdhdf_loadMetaData(inFile);

  nbeam = inFile->nBeam;
  sdhdf_allocateBeamMemory(outFile,nbeam);

  nFluxCalChan=0;
  fluxCal = (sdhdf_fluxCalibration *)malloc(sizeof(sdhdf_fluxCalibration)*MAX_POL_CAL_CHAN); //  Should change to MAX_FLUXCAL OR SIMILAR ** FIX ME
  sdhdf_loadFluxCal(fluxCal,&nFluxCalChan,"parkes","UWL",fluxCalFile); // REMOVE HARDCODE

  // FIX ME: SHOULD KNOW IF DATA FLAGGED OR NOT HERE IN THE FLUX CAL
  // SHOULD SET FLAGS AT END TO SAY WHAT HAS NOT BEEN CORRECTLY FLAGGED
  
  if (verbose==2)
    {
      FILE *fout;
      fout = fopen("sdhdf_calibrate.fluxcalInput","w");
      for (i=0;i<nFluxCalChan;i++)
	fprintf(fout,"%d %g %g %g\n",i,fluxCal[i].freq,fluxCal[i].scalAA,fluxCal[i].scalBB);
      fclose(fout);
    }
  interpCoeff_fluxCalAA = (float **)malloc(sizeof(float *)*nFluxCalChan);
  interpCoeff_fluxCalBB = (float **)malloc(sizeof(float *)*nFluxCalChan);
  fluxCalFreq = (float *)malloc(sizeof(float)*nFluxCalChan);
  fluxCalValAA = (float *)malloc(sizeof(float)*nFluxCalChan);
  fluxCalValBB = (float *)malloc(sizeof(float)*nFluxCalChan);

  for (i=0;i<nFluxCalChan;i++)
    {
      interpCoeff_fluxCalAA[i] = (float *)malloc(sizeof(float)*4);
      interpCoeff_fluxCalBB[i] = (float *)malloc(sizeof(float)*4);

      fluxCalFreq[i]  = fluxCal[i].freq;
      fluxCalValAA[i] = fluxCal[i].scalAA;
      fluxCalValBB[i] = fluxCal[i].scalBB;
    }
  TKcmonot(nFluxCalChan,fluxCalFreq,fluxCalValAA,interpCoeff_fluxCalAA);
  TKcmonot(nFluxCalChan,fluxCalFreq,fluxCalValBB,interpCoeff_fluxCalBB);
  //  for (freq=704;freq<4032;freq+=0.1)
  //    printf("fluxcalInterp: %g %g %g\n",freq,sdhdf_splineValue(freq,nFluxCalChan,fluxCalFreq,interpCoeff_fluxCalAA),sdhdf_splineValue(freq,nFluxCalChan,fluxCalFreq,interpCoeff_fluxCalBB));

  if (verbose==2)
    {
      debugOut1 = fopen("sdhdf_calibrate.interpolations","w");
    }
  printf("Complete interpolation\n");
  for (b=0;b<nbeam;b++)
    {
      for (j=0;j<inFile->beam[b].nBand;j++)
	{
	  printf("Processing band %d\n",j);
	  sdhdf_loadBandData(inFile,b,j,1);
	  nchan  = inFile->beam[b].bandHeader[j].nchan;
	  npol   = inFile->beam[b].bandHeader[j].npol;
	  ndump  = inFile->beam[b].bandHeader[j].ndump;
	  nchanFreq = inFile->beam[b].bandData[j].astro_data.nFreqDumps;
	  sdhdf_copyAttributes(inFile->beam[b].bandData[j].astro_obsHeaderAttr,inFile->beam[b].bandData[j].nAstro_obsHeaderAttributes,dataAttributes,&nDataAttributes);
	  sdhdf_copyAttributes(inFile->beam[b].bandData[j].astro_obsHeaderAttr_freq,inFile->beam[b].bandData[j].nAstro_obsHeaderAttributes_freq,freqAttributes,&nFreqAttributes);

	  // SHOULD CHECK IF THE CAL IS ON OR OFF -- FIX ME

	  // Load in the noise source data
	  sdhdf_loadBandData(inFile,b,j,2);
	  sdhdf_loadBandData(inFile,b,j,3);
	  load_nchanCal = inFile->beam[b].calBandHeader[j].nchan;
	  load_calFreq = (float *)malloc(sizeof(float)*load_nchanCal);
	  load_calAA = (float *)malloc(sizeof(float)*load_nchanCal);
	  load_calBB = (float *)malloc(sizeof(float)*load_nchanCal);
	  load_calRe = (float *)malloc(sizeof(float)*load_nchanCal);
	  load_calIm = (float *)malloc(sizeof(float)*load_nchanCal);

	  interpCoeff_AA = (float **)malloc(sizeof(float *)*load_nchanCal);
	  interpCoeff_BB = (float **)malloc(sizeof(float *)*load_nchanCal);
	  interpCoeff_Re = (float **)malloc(sizeof(float *)*load_nchanCal);
	  interpCoeff_Im = (float **)malloc(sizeof(float *)*load_nchanCal);
	  for (k=0;k<load_nchanCal;k++)
	    {
	      interpCoeff_AA[k] = (float *)malloc(sizeof(float)*4);
	      interpCoeff_BB[k] = (float *)malloc(sizeof(float)*4);
	      interpCoeff_Re[k] = (float *)malloc(sizeof(float)*4);
	      interpCoeff_Im[k] = (float *)malloc(sizeof(float)*4);
	    }
	  
	  // Note this averages the noise source in the entire observation
	  sdhdf_noiseSourceOnOff(inFile,b,j,load_calFreq,load_calAA,load_calBB,load_calRe,load_calIm);
	  for (k=0;k<load_nchanCal;k++)
	    printf("MEASURECAL %g %g %g %g %g\n",load_calFreq[k],load_calAA[k],load_calBB[k],load_calRe[k],load_calIm[k]);
	  // NOTE: THIS INTERPOLATION INCLUDES ALL THE RFI! FIX ME
	  // EXTRACT_BAND SHOULD EXTRACT THE CAL AS WELL
	  // NEED TO CHECK THIS INTERPOLATION MORE
	  //
	  TKcmonot(load_nchanCal,load_calFreq,load_calAA,interpCoeff_AA);
	  TKcmonot(load_nchanCal,load_calFreq,load_calBB,interpCoeff_BB);
	  TKcmonot(load_nchanCal,load_calFreq,load_calRe,interpCoeff_Re);
	  TKcmonot(load_nchanCal,load_calFreq,load_calIm,interpCoeff_Im);
	  
	  if (verbose==1 && j==0)
	    {
	      for (freq=1600;freq<1728;freq+=0.1)
		printf("Using noise source: %g %g %g %g %g\n",freq,sdhdf_splineValue(freq,load_nchanCal,load_calFreq,interpCoeff_AA),sdhdf_splineValue(freq,load_nchanCal,load_calFreq,interpCoeff_BB),sdhdf_splineValue(freq,load_nchanCal,load_calFreq,interpCoeff_Re),sdhdf_splineValue(freq,load_nchanCal,load_calFreq,interpCoeff_Im));
	      //	      for (k=0;k<use_nchanCal;k++)
		//
	    }

	  out_freq  = (float *)malloc(sizeof(float)*nchan*nchanFreq);
	  out_data  = (float *)calloc(sizeof(float),nchan*npol*ndump);
	  	  
	  for (k=0;k<inFile->beam[b].bandHeader[j].ndump;k++)
	    {
	      
	      for (ii=0;ii<nchan;ii++)
		{
		  freq = inFile->beam[b].bandData[j].astro_data.freq[k*nchan+ii];
		  final_aa = final_bb = final_rab = final_iab = 0;
		  out_freq[k*nchan+ii] = freq;  // NOTE: Assuming that the frequency is constant in all bands - not true if Doppler corrected -- FIX ME

		  measured_aa = inFile->beam[b].bandData[j].astro_data.pol1[ii+k*nchan];
		  measured_bb = inFile->beam[b].bandData[j].astro_data.pol2[ii+k*nchan];
		  measured_rab = inFile->beam[b].bandData[j].astro_data.pol3[ii+k*nchan];
		  measured_iab = inFile->beam[b].bandData[j].astro_data.pol4[ii+k*nchan];

		  cal_aa  = sdhdf_splineValue(freq,load_nchanCal,load_calFreq,interpCoeff_AA);
		  cal_bb  = sdhdf_splineValue(freq,load_nchanCal,load_calFreq,interpCoeff_BB);
		  cal_rab = sdhdf_splineValue(freq,load_nchanCal,load_calFreq,interpCoeff_Re);
		  cal_iab = sdhdf_splineValue(freq,load_nchanCal,load_calFreq,interpCoeff_Im);

		  printf("testCAL: %g %g\n",freq,(pow(cal_rab,2)+pow(cal_iab,2))/cal_aa/cal_bb);
		  
		  // Stabilise the data set using the noise source		  
		  if (stabiliseType==1 || stabiliseType == 2)
		    {
		      stabilised_aa = measured_aa/cal_aa;
		      stabilised_bb = measured_bb/cal_bb;
		    }
		  else if (stabiliseType==3) // Using gain from AB
		    {
		      stabilised_aa = measured_aa*cal_bb/(pow(cal_rab,2)+pow(cal_iab,2));
		      stabilised_bb = measured_bb*cal_aa/(pow(cal_rab,2)+pow(cal_iab,2));
		    }
		  // SHOULD CHECK |AB_cal|^2 = AA_cal*BB_cal *** DO THIS ****
		  if (stabiliseType==1 || stabiliseType==3) // Not assuming perfect correlation of the cal
		    stabilise_normFactor = (pow(cal_rab,2)+pow(cal_iab,2));
		  else if (stabiliseType==2)
		    stabilise_normFactor = sqrt(pow(cal_rab,2)+pow(cal_iab,2))*sqrt(cal_aa*cal_bb);
		  else
		    stabilise_normFactor = 1;
		  
		  stabilised_rab = (measured_rab*cal_rab + measured_iab*cal_iab)/stabilise_normFactor;
		  stabilised_iab = (measured_iab*cal_rab - measured_rab*cal_iab)/stabilise_normFactor;

		  fluxValAA = sdhdf_splineValue(freq,nFluxCalChan,fluxCalFreq,interpCoeff_fluxCalAA);
		  fluxValBB = sdhdf_splineValue(freq,nFluxCalChan,fluxCalFreq,interpCoeff_fluxCalBB);
		  fluxScale = fluxValAA + fluxValBB;
		  if (verbose==2)
		    {
		      fprintf(debugOut1,"%.6f %g %g\n",freq,fluxValAA,fluxValBB);
		    }
		  // Form Jones matrix
		  sdhdf_complex_matrix_2x2(rho,stabilised_aa,stabilised_rab-I*stabilised_iab,
					   stabilised_rab+I*stabilised_iab,stabilised_bb);
		  
		  paraAng = inFile->beam[b].bandData[j].astro_obsHeader[k].paraAngle * M_PI/180.0;
		  sdhdf_complex_matrix_2x2(R_feed_pa,
					   I*sin(alpha+paraAng),
					   -I*cos(alpha+paraAng),
					   -I*cos(alpha+paraAng),
					   -I*sin(alpha+paraAng));
		  sdhdf_complex_matrix_2x2_dagger(R_feed_pa,R_feed_pa_dag);
		  sdhdf_copy_complex_matrix_2x2(finalJ,R_feed_pa);
		  sdhdf_multiply_complex_matrix_2x2(finalJ,rho);
		  sdhdf_multiply_complex_matrix_2x2(finalJ,R_feed_pa_dag);

		  final_aa  = fluxScale*creal(finalJ[0][0]);
		  final_bb  = fluxScale*creal(finalJ[1][1]);
		  final_rab = fluxScale*creal(finalJ[1][0]);
		  final_iab = fluxScale*cimag(finalJ[0][1]);
		  
		  out_data[ii+nchan*k*npol]         = final_aa;
		  out_data[ii+nchan*k*npol+nchan]   = final_bb;
		  out_data[ii+nchan*k*npol+2*nchan] = final_rab;
		  out_data[ii+nchan*k*npol+3*nchan] = final_iab;	      
		}
	    }
       
	  // Do the calibration
	  sdhdf_writeSpectrumData(outFile,inFile->beamHeader[b].label,inFile->beam[b].bandHeader[j].label,b,j, out_data,out_freq,ndump,nchan,npol,ndump,0,dataAttributes,nDataAttributes,freqAttributes,nFreqAttributes);	  
	  sdhdf_releaseBandData(inFile,b,j,1); 		      
	
	  
	  free(out_freq);
	  free(out_data);
	  free(load_calFreq);
	  free(load_calAA);
	  free(load_calBB);
	  free(load_calRe);
	  free(load_calIm);
	  //
	  for (k=0;k<load_nchanCal;k++)
	    {
	      free(interpCoeff_AA[k]);
	      free(interpCoeff_BB[k]);
	      free(interpCoeff_Re[k]);
	      free(interpCoeff_Im[k]); 
	    }
	  free(interpCoeff_AA);
	  free(interpCoeff_BB);
	  free(interpCoeff_Re);
	  free(interpCoeff_Im); 	  
	}
    }
  if (verbose==2)
    fclose(debugOut1);

  sdhdf_writeHistory(outFile,inFile->history,inFile->nHistory);
  sdhdf_copyRemainder(inFile,outFile,0);
  
  sdhdf_closeFile(inFile);
  sdhdf_closeFile(outFile);


  for (i=0;i<nFluxCalChan;i++)
    {
      free(interpCoeff_fluxCalAA[i]);
      free(interpCoeff_fluxCalBB[i]);
    }

  free(interpCoeff_fluxCalAA);
  free(interpCoeff_fluxCalBB);

  
  free(inFile); free(outFile);
  free(fluxCal);  
  free(fluxCalFreq); free(fluxCalValAA); free(fluxCalValBB);
}

