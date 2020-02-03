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

int main(int argc,char *argv[])
{
  int i,ii,j,k,p,l,kp,b;
  int nchan,npol;
  char fname[MAX_STRLEN];
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
  int nh;
  float *in_data,*out_data,*in_freq,*out_freq;
  float *out_Tdata,*out_Fdata;
  long out_nchan,out_npol,out_ndump,nsd;
  double av1,av2,av3,av4,avFreq;
  double tav=0;
  char extension[MAX_STRLEN];
  sdhdf_eopStruct *eop;
  int nEOP=0;
  char args[MAX_STRLEN]="";
   
  strcpy(oname,"sdhdf_modify_output.hdf");

  eop = (sdhdf_eopStruct *)malloc(sizeof(sdhdf_eopStruct)*MAX_EOP_LINES);
  
  
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
      else if (strcmp(argv[i],"-T")==0)
	{sdhdf_add1arg(args,argv[i]); tScrunch=1;}
      else if (strcmp(argv[i],"-F")==0)
	{sdhdf_add1arg(args,argv[i]); fScrunch=1;}
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

	  sdhdf_initialiseFile(inFile);
	  sdhdf_initialiseFile(outFile);

	  sdhdf_openFile(fname,inFile,1);
	  sprintf(oname,"%s.%s",fname,extension);
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
		  
		  inBandParams = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*inFile->beam[b].nBand);      
		  sdhdf_copyBandHeaderStruct(inFile->beam[b].bandHeader,inBandParams,inFile->beam[b].nBand);

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
		      in_freq = (float *)malloc(sizeof(float)*nchan);
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
	      // Time averaging
		      if (tScrunch==1)
			{
			  for (j=0;j<inFile->beam[b].bandHeader[ii].ndump;j++)
			    {
			      for (k=0;k<nchan;k++)
				{
				  // NOT ACCOUNTING FOR WEIGHTINGS			  
				  if (npol==1)
				    out_Tdata[k] += inFile->beam[b].bandData[ii].astro_data.pol1[k+j*nchan]/(float)nsd;
				  else if (npol==2)
				    {
				      out_Tdata[k]         += inFile->beam[b].bandData[ii].astro_data.pol1[k+j*nchan]/(float)nsd;
				      out_Tdata[k+nchan]   += inFile->beam[b].bandData[ii].astro_data.pol2[k+j*nchan]/(float)nsd;
				    }
				  else
				    {
				      out_Tdata[k]         += inFile->beam[b].bandData[ii].astro_data.pol1[k+j*nchan]/(float)nsd;
				      out_Tdata[k+nchan]   += inFile->beam[b].bandData[ii].astro_data.pol2[k+j*nchan]/(float)nsd;
				      out_Tdata[k+2*nchan] += inFile->beam[b].bandData[ii].astro_data.pol3[k+j*nchan]/(float)nsd;
				      out_Tdata[k+3*nchan] += inFile->beam[b].bandData[ii].astro_data.pol4[k+j*nchan]/(float)nsd;		      
				    }
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
				    out_Tdata[k+j*nchan*npol]         = inFile->beam[b].bandData[ii].astro_data.pol1[k+j*nchan]/(float)nsd;
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
			  for (j=0;j<out_ndump;j++)
			    {
			      for (k=0;k<out_nchan;k++)
				{
				  if (npol==1)
				    out_Fdata[k+out_nchan*j*4]             = out_Tdata[k+out_nchan*j*4];
				  else if (npol==2)
				    {
				      out_Fdata[k+out_nchan*j*4]             = out_Tdata[k+out_nchan*j*4];
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
				  else if (npol==4)
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
				out_data[k+out_nchan*j*4] = (out_Fdata[k+out_nchan*j*4]+out_Fdata[k+out_nchan*j+4+out_nchan])/2.;
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
				  p1 = out_data[k];
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
			      
			      for (k=0;k<out_nchan;k++)
				out_freq[k]*=(1.0-vOverC);
			    }
			  tav += inFile->beam[b].bandHeader[ii].dtime;

			  // FIX ME
			  //			  if (bary==1)
			  //			    strcpy(outFile->frequency_attr.frame,"barycentric");
			  //			  else if (lsr == 1)
			  //			    strcpy(outFile->frequency_attr.frame,"LSR");
			}
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

  free(inFile);
  free(outFile);
  free(eop);
}



