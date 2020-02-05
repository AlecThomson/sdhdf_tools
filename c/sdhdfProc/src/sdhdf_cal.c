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
// Software to calibrate data files to form Stokes (I,Q,U,V)
//
// Usage:
// sdhdf_cal <filename.hdf>  <filename.hdf> ...
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sdhdfProc.h"

#define VERSION "v0.5"

void help()
{
  printf("sdhdf_cal %s (SDHDFProc %s)\n",VERSION,SOFTWARE_VER);
  printf("Authors: G. Hobbs\n");
  printf("Purpose: to calibrate data in a file to Stokes parameters (I,Q,U,V) \n");
  printf("\n");
  printf("Command line arguments:\n\n");
  exit(1);

}

int main(int argc,char *argv[])
{
  int i,j,k,l;
  int beam,band,dump;
  int npol,nchanAst,ndumpAst,nchanCal,ndumpCal;
  float tdumpAst,tdumpCal;
  float scale;
  int nFiles=0;
  char fname[MAX_FILES][MAX_STRLEN];
  float freq;
  float pol1,pol2,pol3,pol4;
  float astI_m,astQ_m,astU_m,astV_m;
  float astI,astQ,astU,astV;
  float con1,con2,con3,con4;  
  float coff1,coff2,coff3,coff4;
  float dc1,dc2,dc3,dc4;
  float dcI,dcQ,dcU,dcV;
  
  float i0 = 1;
  float *deltaG_cal,*deltaPhi_cal,*alphaB_cal;
  float deltaG_ast,deltaPhi_ast,alphaB_ast;
  float freq1,freq2;
  int   idump;
  int   ii0,ii1;
  int res;
  sdhdf_fileStruct *inFile;
  float chanbw_cal;
  float m,c;

  float m1[4][4],mInv[4][4],m2[4][4];

  // REMOVE LATER
  float *resultI,*resultQ,*resultU,*resultV;

  // Calibrator information
  sdhdf_tcal_struct *tcalData;
  int nTcal;
  double tcalA,tcalB;
  
  tcalData = (sdhdf_tcal_struct *)malloc(sizeof(sdhdf_tcal_struct)*3328); // 3328 = number of channels in Tcal files  

  // HARDCODED
  nTcal = sdhdf_loadTcal(tcalData,"/u/hob044/software/new_c/sdhdfProc/runtime/Tcal_result.12.6.6.dat");


  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-h")==0)
	help();
      else
	strcpy(fname[nFiles++],argv[i]);
    }

  inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct));

  for (i=0;i<nFiles;i++)
    {
      sdhdf_initialiseFile(inFile);
      if (sdhdf_openFile(fname[i],inFile,1)==-1)
	printf("Warning: unable to open file >%s<. Skipping\n",fname[i]);
      else
	{
	  sdhdf_loadMetaData(inFile);

	  for (beam=0;beam<inFile->nBeam;beam++)
	    {
	      //	      for (band=0;band<inFile->beam[beam].nBand;band++)
	      for (band=1;band<2;band++)
		{
		  npol     = inFile->beam[beam].bandHeader[band].npol;
		  nchanAst = inFile->beam[beam].bandHeader[band].nchan;
		  nchanCal = inFile->beam[beam].calBandHeader[band].nchan;
		  ndumpAst = inFile->beam[beam].bandHeader[band].ndump;
		  ndumpCal = inFile->beam[beam].calBandHeader[band].ndump;

		  tdumpAst = inFile->beam[beam].bandHeader[band].dtime;
		  tdumpCal = inFile->beam[beam].calBandHeader[band].dtime;

		  sdhdf_loadBandData(inFile,beam,band,1);
		  sdhdf_loadBandData(inFile,beam,band,2);
		  sdhdf_loadBandData(inFile,beam,band,3);

		  deltaG_cal = (float *)malloc(sizeof(float)*nchanCal*ndumpCal);
		  deltaPhi_cal = (float *)malloc(sizeof(float)*nchanCal*ndumpCal);
		  alphaB_cal = (float *)malloc(sizeof(float)*nchanCal*ndumpCal);
		  
		  resultI = (float *)calloc(sizeof(float),nchanAst);
		  resultQ = (float *)calloc(sizeof(float),nchanAst);
		  resultU = (float *)calloc(sizeof(float),nchanAst);
		  resultV = (float *)calloc(sizeof(float),nchanAst);
		  
		  
		  for (dump=0;dump<ndumpCal;dump++)
		    {
		      for (j=0;j<nchanCal;j++)
			{
 			  freq  = inFile->beam[beam].bandData[band].cal_on_data.freq[j];
			  con1  = inFile->beam[beam].bandData[band].cal_on_data.pol1[j+dump*nchanCal];
			  con2  = inFile->beam[beam].bandData[band].cal_on_data.pol2[j+dump*nchanCal];
			  con3  = inFile->beam[beam].bandData[band].cal_on_data.pol3[j+dump*nchanCal];
			  con4  = inFile->beam[beam].bandData[band].cal_on_data.pol4[j+dump*nchanCal];
			  coff1 = inFile->beam[beam].bandData[band].cal_off_data.pol1[j+dump*nchanCal];
			  coff2 = inFile->beam[beam].bandData[band].cal_off_data.pol2[j+dump*nchanCal];
			  coff3 = inFile->beam[beam].bandData[band].cal_off_data.pol3[j+dump*nchanCal];
			  coff4 = inFile->beam[beam].bandData[band].cal_off_data.pol4[j+dump*nchanCal];
			  dc1 = con1-coff1;
			  dc2 = con2-coff2;
			  dc3 = con3-coff3;
			  dc4 = con4-coff4;

			  // Should check these
			  dcI = dc1+dc2;
			  dcQ = dc1-dc2;
			  dcU = 2*dc3;
			  dcV = 2*dc4;
			  deltaG_cal[dump*nchanCal+j]   = 2*dcQ/dcI; 
			  deltaPhi_cal[dump*nchanCal+j] = atan2(-dcV,dcU);

			  // UPDATE I0 SHOULD BE SET FROM TCAL OR SCAL *****
			  sdhdf_get_tcal(tcalData,nTcal,freq,&tcalA,&tcalB);
			  //
			  // FIX ME -- WHAT SHOULD THIS BE
			  //
			  //			  printf("Here with %g %g\n",tcalA,tcalB);
			  i0 = (tcalA+tcalB)/2.; 
			  alphaB_cal[dump*nchanCal+j]   = dcI/i0;

			  //			  if (dump==1)
			  //			    printf("result1: %.6f %g %g %g %g %g %g %g %g %g %d %d\n",freq,con1,con2,con3,con4,coff1,coff2,coff3,coff4,deltaG_cal[dump*nchanCal+j],nchanCal,ndumpCal);
			    //			    printf("result1: %.6f %g %g %g %g %g %g %g\n",freq,dcI,dcQ,dcU,dcV,
			    //				   deltaG_cal[dump*nchanCal+j],deltaPhi_cal[dump*nchanCal+j],alphaB_cal[dump*nchanCal+j]);


			}
		    }

		  // Now determine the calibration parameters on the resolution of the astronomy spectrum
		  chanbw_cal = inFile->beam[beam].bandData[band].cal_on_data.freq[1]-inFile->beam[beam].bandData[band].cal_on_data.freq[0];
		  for (dump=0;dump<ndumpAst;dump++)
		    {
		      for (j=0;j<nchanAst;j++)
			{
			  freq = inFile->beam[beam].bandData[band].astro_data.freq[j];
			  pol1 = inFile->beam[beam].bandData[band].astro_data.pol1[j+dump*nchanAst];
			  pol2 = inFile->beam[beam].bandData[band].astro_data.pol2[j+dump*nchanAst];
			  pol3 = inFile->beam[beam].bandData[band].astro_data.pol3[j+dump*nchanAst];
			  pol4 = inFile->beam[beam].bandData[band].astro_data.pol4[j+dump*nchanAst];

			  // Rescaling to cal 
			  scale = ((nchanAst/nchanCal/2.)*tdumpCal/tdumpAst);
			  pol1 *= scale;
			  pol2 *= scale;
			  pol3 *= scale;
			  pol4 *= scale;
			  
			  astI_m = pol1+pol2;
			  astQ_m = pol1-pol2;
			  astU_m = 2*pol3;
			  astV_m = 2*pol4;

			  ii0 = (int)((freq-inFile->beam[beam].bandData[band].cal_on_data.freq[0])/chanbw_cal);
			  if (ii0 >= nchanCal-1)
			    {
			      ii0 = nchanCal-2;
			      ii1 = nchanCal-1;
			    }
			  else
			    ii1 = ii0+1;
			  freq1 = inFile->beam[beam].bandData[band].cal_on_data.freq[ii0];
			  freq2 = inFile->beam[beam].bandData[band].cal_on_data.freq[ii1];

			  //
			  // SHOULD DO A MULTI POINT INTERPOLATION ***
			  //
			  // Need to do a much better smoothing/interpolation than done here
			  // to ensure the functions are smooth
			  //

			  // This doesn't need to be within the channel loop
			  idump = (int)((dump*tdumpAst)/tdumpCal+0.5);
			  if (idump > ndumpCal-1) idump = ndumpCal-1;
			  //			  printf("idump = %d %d\n",idump,ndumpCal);
			  m = (deltaG_cal[idump*nchanCal+ii1]-deltaG_cal[idump*nchanCal+ii0])/(freq2-freq1);
			  c = (deltaG_cal[idump*nchanCal+ii1] - m*freq2);
			  deltaG_ast = m*freq+c;

			  m = (deltaPhi_cal[idump*nchanCal+ii1]-deltaPhi_cal[idump*nchanCal+ii0])/(freq2-freq1);
			  c = (deltaPhi_cal[idump*nchanCal+ii1] - m*freq2);
			  deltaPhi_ast = m*freq+c;
			  
			  m = (alphaB_cal[idump*nchanCal+ii1]-alphaB_cal[idump*nchanCal+ii0])/(freq2-freq1);
			  c = (alphaB_cal[idump*nchanCal+ii1] - m*freq2);
			  alphaB_ast = m*freq+c;
			  
			  m1[0][0] = 1;
			  m1[0][1] = deltaG_ast/2.;
			  m1[0][2] = 0;
			  m1[0][3] = 0;
			  m1[1][0] = deltaG_ast/2.;
			  m1[1][1] = 1;
			  m1[1][2] = 0;
			  m1[1][3] = 0;
			  m1[2][0] = 0;
			  m1[2][1] = 0;
			  m1[2][2] = cos(deltaPhi_ast);
			  m1[2][3] = -sin(deltaPhi_ast);
			  m1[3][0] = 0;
			  m1[3][1] = 0;
			  m1[3][2] = sin(deltaPhi_ast);
			  m1[3][3] = cos(deltaPhi_ast);
			
			  //			  printf("matrix in: %d %d %d (%d %d %d [%.5f %.5f %.5f %g %g]) %.5f %g %g\n",band,dump,j,idump,ii0,ii1,freq,freq1,freq2,deltaG_cal[idump*nchanCal+ii0],deltaG_cal[idump*nchanCal+ii1], deltaG_ast/2,cos(deltaPhi_ast),sin(deltaPhi_ast));

			  res = sdhdf_inv4x4(m1,mInv);

			  if (res==-1)
			    {
			      printf("ERROR: unable to invert matrix\n");
			      exit(1);
			    }
			  // FIX: SHOULD INTERPOLATE THE MUELLER MATRIX ON THE SCALE OF THE CAL - NOT OF THE ASTRONOMY SIGNAL
			  /*
			    for (k=0;k<4;k++)
			    {
			    for (l=0;l<4;l++)
			    printf("(%g,%g) ",m1[k][l],mInv[k][l]);
			    printf("\n");
			      
			    }
			    printf("\n"); 
			  exit(1);
			  */
			  
			  // Do the calibration
			  astI = astI_m*mInv[0][0]+astQ_m*mInv[0][1]+astU_m*mInv[0][2]+astV_m*mInv[0][3]; 
			  astQ = astI_m*mInv[1][0]+astQ_m*mInv[1][1]+astU_m*mInv[1][2]+astV_m*mInv[1][3]; 
			  astU = astI_m*mInv[2][0]+astQ_m*mInv[2][1]+astU_m*mInv[2][2]+astV_m*mInv[2][3]; 
			  astV = astI_m*mInv[3][0]+astQ_m*mInv[3][1]+astU_m*mInv[3][2]+astV_m*mInv[3][3]; 

			  astI/=alphaB_ast;
			  astQ/=alphaB_ast;
			  astU/=alphaB_ast;
			  astV/=alphaB_ast;

			  resultI[j] += astI;
			  resultQ[j] += astQ;
			  resultU[j] += astU;
			  resultV[j] += astV;
			  //			  printf("result2: %.6f %g %g %g %g %d %d %g %g\n",freq,pol1,pol2,pol3,pol4,nchanAst,nchanCal,tdumpAst,tdumpCal);
			  //			  printf("result2: %.6f %g %g %g %g delta: %g %g %g measured: %g %g %g %g minv: %g %g %g %g\n",freq,astI,astQ,astU,astV,deltaG_ast,deltaPhi_ast,alphaB_ast,astI_m,astQ_m,astU_m,astV_m,mInv[0][0],mInv[0][1],mInv[0][2],mInv[0][3]);
			}
		    }
		  for (j=0;j<nchanAst;j++)
		    {
		      printf("result3: %.6f %g %g %g %g\n",inFile->beam[beam].bandData[band].astro_data.freq[j],resultI[j]/ndumpAst,resultQ[j]/ndumpAst,resultU[j]/ndumpAst,resultV[j]/ndumpAst);
		    }
		
		  
		  free(deltaG_cal);
		  free(deltaPhi_cal);
		  free(alphaB_cal);
		  free(resultI);
		  free(resultQ);
		  free(resultU);
		  free(resultV);
		}
	    }


	  sdhdf_closeFile(inFile);
	}
    }

  free(inFile);
  free(tcalData);
}
