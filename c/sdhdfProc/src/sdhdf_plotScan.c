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

// gcc -lm -o sdhdf_plotScan sdhdf_plotScan.c sdhdfProc.c -I../hdf5-1.10.4/src/ ../hdf5-1.10.4/src/.libs/libhdf5.a -ldl -lz -lcpgplot  -L/pulsar/psr/software/stable/stretch/lib/ -I//pulsar/psr/software/stable/stretch/include -lcalceph -Isofa/20190722/c/src/ -Lsofa/20190722/c/src -lsofa_c


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "hdf5.h"
#include "sdhdfProc.h"
#include <cpgplot.h>
#include "TKnonlinearFit.h"

#define MAX_CHAN 262144 // SHOULD UPDATE THIS!!
#define MAX_SUBINT 512 // SHOULD UPDATE THIS
#define MAX_BANDS 26   // SHOULD UPDATE THIS

typedef struct modelFitStruct {
  double amp;
  double offset;
  double width;
  double baseline; 
} modelFitStruct;


double haversine(double centre_long,double centre_lat,double src_long,double src_lat);


int main(int argc,char *argv[])
{
  int i,j,sb;
  sdhdf_fileStruct *inFile;
  char fname[MAX_STRLEN];
  int nchan,npol,ndumps;
  int ibeam=0;
  int band0 = 0; // Used to identify number of subbands
  float gx[2],gy[2];
  int setMinMax=1;
  float minx,maxx,miny,maxy;
  float mx,my;
  char key;

  int fitted=0;
  
  float timeVal[MAX_BANDS][MAX_SUBINT];
  float fluxAA[MAX_BANDS][MAX_SUBINT],fluxBB[MAX_BANDS][MAX_SUBINT];
  float fluxI[MAX_BANDS][MAX_SUBINT];
  float fx[MAX_BANDS][MAX_SUBINT],fy[MAX_BANDS][MAX_SUBINT];
  float lineX[2],lineY[2];
  float fModelX[MAX_BANDS][1024];
  float fModelY[MAX_BANDS][1024];
  modelFitStruct modelFit[MAX_BANDS];
  float maxFlux=0,minFlux=1e99;
  int nFile=0;
  int useWeights=1;
  float wSum;
  int sbHighlight=-1;
  int okay;
  int sub;
  int recalc=1;
  int logy=-1;
  int label=-1;
  int normalise=-1;
  int offset=-1;
  int xaxis=0;
  float maxy_norm,miny_norm;
  float offsetVal = 1;
  float onAA,onBB;
  char offsetStr[1024];
  char title[1024];
  int select=-1;
  double sigma;
  
  double ra0 = 0;
  double dec0 = 0;
  
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }

  sdhdf_initialiseFile(inFile);


  for (i=0;i<argc;i++)
    {
     if (strcmp(argv[i],"-f")==0)
       strcpy(fname,argv[++i]);
    }

  if (sdhdf_openFile(fname,inFile,1)==-1)
    {
      printf("Unable to open file >%s<\n",fname);
      free(inFile);
      exit(1);
    }

  sdhdf_loadMetaData(inFile);

  ndumps = inFile->beam[ibeam].bandHeader[band0].ndump;

  // FIX ME: Should read these from a file
  if (strcmp(inFile->beamHeader[ibeam].source,"0407-658")==0)
    {
      printf("Setting default position to 0407-657\n");
      //      04 08 20.37884 -65 45 09.0806
      ra0 = (4.+8/60.0+20.37884/60./60.)*180./12.;
      dec0 = -65.0-45.0/60.-9.0806/60./60.;
    }
  else if (strcmp(inFile->beamHeader[ibeam].source,"1934_RASCAN")==0)
    {
      printf("Setting default position to 1934-638\n");
      //  19:39:25.01 -63:42:45.7
      ra0 = (19.+39/60.0+25.01/60./60.)*180./12.;
      dec0 = -63.0-42.0/60.-45.7/60./60.;
    }
  else if (strcmp(inFile->beamHeader[ibeam].source,"1934_DCSCAN")==0)
    {
      printf("Setting default position to 1934-638\n");
      //  19:39:25.01 -63:42:45.7
      ra0 = (19.+39/60.0+25.01/60./60.)*180./12.;
      dec0 = -63.0-42.0/60.-45.7/60./60.;
    }
  else
    {
      printf("Setting default position to start of observation. Do not know source >%s<\n",
	     inFile->beamHeader[ibeam].source);
      ra0 = inFile->beam[ibeam].bandData[band0].astro_obsHeader[0].raDeg;
      dec0 = inFile->beam[ibeam].bandData[band0].astro_obsHeader[0].decDeg;
    }
  // Define beam model
  for (i=0;i<inFile->beam[ibeam].nBand;i++)
    {
      modelFit[i].amp = 1;
      modelFit[i].offset = 0;
      modelFit[i].baseline = 0;
      // FWHM (Gaussian) = 2.355 sigma
      modelFit[i].width = 1.02*SPEED_LIGHT/(inFile->beam[ibeam].bandHeader[i].fc*1e6)*180./M_PI/64.0/2.355; // HARD CODE TO 64m
      printf("fc setting is %g\n",inFile->beam[ibeam].bandHeader[i].fc);
    }
  
  for (i=0;i<inFile->beam[ibeam].nBand;i++)
    {
      nchan = inFile->beam[ibeam].bandHeader[i].nchan;
      npol  = inFile->beam[ibeam].bandHeader[i].npol;

      sdhdf_loadBandData(inFile,ibeam,i,1);
      
      for (sub=0;sub<inFile->beam[ibeam].bandHeader[i].ndump;sub++)
 	{
	  fluxAA[i][sub] = fluxBB[i][sub] = 0;
	  wSum=0;
	  for (j=0;j<inFile->beam[ibeam].bandHeader[i].nchan;j++)
	    {
	      if (inFile->beam[ibeam].bandData[i].astro_data.flag[j] == 0)
		{
		  onAA = inFile->beam[ibeam].bandData[i].astro_data.pol1[j+sub*nchan]; 
		  onBB = inFile->beam[ibeam].bandData[i].astro_data.pol2[j+sub*nchan]; 
		  fluxAA[i][sub] += onAA;
		  fluxBB[i][sub] += onBB;
		  wSum+=1;
		}
	    }
	  fluxAA[i][sub]/=(double)wSum;
	  fluxBB[i][sub]/=(double)wSum;
	  fluxI[i][sub] = fluxAA[i][sub] + fluxBB[i][sub];
	  
	  if (maxFlux < fluxI[i][sub]) maxFlux = fluxI[i][sub];
	  if (minFlux > fluxI[i][sub]) minFlux = fluxI[i][sub];
	  //	  printf("Band: %d, subint %d aa = %g bb = %g\n",i,sub,fluxAA[i][sub],fluxBB[i][sub]);

	}
      sdhdf_releaseBandData(inFile,ibeam,i,1);
    }
  
// Now make a picture
  cpgbeg(0,"/xs",1,1);
  cpgsch(1.4);
  cpgscf(2);
  cpgslw(2);
  cpgask(0);

  
  do {
    if (recalc==1)
      {
	int t=0;
	for (i=0;i<inFile->beam[ibeam].nBand;i++)
	  {	    
	    for (sub=0;sub<ndumps;sub++)
	      {       	
		if (xaxis==0)
		  timeVal[i][sub] = sub;
		else if (xaxis==1)
		  timeVal[i][sub] = inFile->beam[ibeam].bandData[band0].astro_obsHeader[sub].timeElapsed;
		else if (xaxis==2)
		  timeVal[i][sub] = inFile->beam[ibeam].bandData[band0].astro_obsHeader[sub].raDeg;
		else if (xaxis==3)
		  timeVal[i][sub] = inFile->beam[ibeam].bandData[band0].astro_obsHeader[sub].decDeg;
		else if (xaxis==4)
		  timeVal[i][sub] = inFile->beam[ibeam].bandData[band0].astro_obsHeader[sub].az;
		else if (xaxis==5)
		  timeVal[i][sub] = inFile->beam[ibeam].bandData[band0].astro_obsHeader[sub].el;
		else if (xaxis==6)
		  timeVal[i][sub] = haversine(ra0,dec0,inFile->beam[ibeam].bandData[band0].astro_obsHeader[sub].raDeg,inFile->beam[ibeam].bandData[band0].astro_obsHeader[sub].decDeg);
		else if (xaxis==7)
		  {
		    timeVal[i][sub] = inFile->beam[ibeam].bandData[band0].astro_obsHeader[sub].raDeg-ra0;
		  }
		else if (xaxis==8)
		  {
		    timeVal[i][sub] = inFile->beam[ibeam].bandData[band0].astro_obsHeader[sub].decDeg-dec0;
		  }
		fx[i][sub] = timeVal[i][sub];
		if (logy==1)
		  fy[i][sub] = log10(fluxI[i][sub]);
		else
		  fy[i][sub] = fluxI[i][sub];		
		if (sub==0)
		  miny_norm = maxy_norm = fy[i][sub];
		else
		  {
		    if (miny_norm > fy[i][sub]) miny_norm = fy[i][sub];
		    if (maxy_norm < fy[i][sub]) maxy_norm = fy[i][sub];
		  }
	      }

	    if (xaxis == 7 || xaxis == 8)
	      {
	      
		for (j=0;j<1024;j++)
		  {		    
		    fModelX[i][j]  = minx + j*(maxx-minx)/1024.;
		    fModelY[i][j]  = modelFit[i].amp*exp(-pow((fModelX[i][j]-modelFit[i].offset),2)/2./pow(modelFit[i].width,2))+modelFit[i].baseline;
		  }
	      }
	    
	    if (normalise==1)
	      {
		for (sub=0;sub<ndumps;sub++)
		  fy[i][sub] = (fy[i][sub] - miny_norm)/(maxy_norm-miny_norm);
	      }
	    if (offset==1 && select==-1)
	      {
		for (sub=0;sub<ndumps;sub++)
		  fy[i][sub] += (offsetVal*i);
	      }
	    
	    // Find max/min ranges
	    if (select==-1 || i==select)
	      {
		for (sub=0;sub<ndumps;sub++)
		  {
		    if (t==0)
		      {
			minx = maxx = fx[i][sub];
			miny = maxy = fy[i][sub];
			t=1;
		      }
		    else
		      {
			if (minx > fx[i][sub]) minx = fx[i][sub];
			if (maxx < fx[i][sub]) maxx = fx[i][sub];
			if (miny > fy[i][sub]) miny = fy[i][sub];
			if (maxy < fy[i][sub]) maxy = fy[i][sub];
		      }
		  }
	      }

	  }
	recalc=0;
      }

    cpgenv(minx,maxx,miny,maxy,0,0);
    sprintf(title,"%s",inFile->beamHeader[ibeam].source);
    if (xaxis==0)
      cpglab("Spectral dump number","Signal strength (arbitrary)",title);
    else if (xaxis==1)
      cpglab("Time since observation start (seconds)","Signal strength (arbitrary)",title);
    else if (xaxis==2)
      cpglab("Right ascension (degrees)","Signal strength (arbitrary)",title);
    else if (xaxis==3)
      cpglab("Declination (degrees)","Signal strength (arbitrary)",title);
    else if (xaxis==4)
      cpglab("Azimuth (degrees)","Signal strength (arbitrary)",title);
    else if (xaxis==5)
      cpglab("Elevation (degrees)","Signal strength (arbitrary)",title);
    else if (xaxis==6)
      {
	sprintf(offsetStr,"Angular offset from (%.2f,%.2f) (degrees)",ra0,dec0);
	cpglab(offsetStr,"Signal strength (arbitrary)",title);
      }
    else if (xaxis==7)
      {
	sprintf(offsetStr,"Right ascension - [%.4f] (degrees)",ra0);
	cpglab(offsetStr,"Signal strength (arbitrary)",title);
      }
    else if (xaxis==8)
      {
	sprintf(offsetStr,"Declination - [%.4f] (degrees)",dec0);
	cpglab(offsetStr,"Signal strength (arbitrary)",title);
      }
    if (select==-1)
      {
	for (i=0;i<inFile->beam[ibeam].nBand;i++)
	  {
	    if (sbHighlight == i)
	      cpgsci(2);
	    cpgline(ndumps,fx[i],fy[i]);	    
	    cpgsci(1);
	    if (label==1)
	      {
		cpgsch(0.6);
		cpgtext(minx+(maxx-minx)*0.1,fy[i][ndumps/4],inFile->beam[ibeam].bandHeader[i].label);
		cpgsch(1.4);
	      }
	  }
      }
    else
      {
	cpgsch(0.8);
	cpgtext(minx+(maxx-minx)*0.1,maxy-(maxy-miny)*0.1,inFile->beam[ibeam].bandHeader[select].label);
	cpgsch(1.4);
	printf("In here %g %g\n",(minx+(maxx-minx)*0.1),maxy-(maxy-miny)*0.1);
	cpgline(ndumps,fx[select],fy[select]);
	cpgpt(ndumps,fx[select],fy[select],6);
	if (xaxis == 7 || xaxis == 8)
	  {
	    cpgsci(2);
	    cpgline(1024,fModelX[select],fModelY[select]);
	    cpgsci(1);
	    lineX[0] = 0; lineX[1]=0;
	    lineY[0] = miny; lineY[1] = maxy;
	    cpgsls(4); cpgline(2,lineX,lineY); cpgsls(1);
	  }
      }
    
    cpgcurs(&mx,&my,&key);
    if (key=='+')
      {
	sbHighlight++;
	if (sbHighlight > inFile->beam[ibeam].nBand-1) sbHighlight = 0;
	if (select!=-1) {select = sbHighlight; recalc=1;}
      }
    else if (key=='-')
      {
	sbHighlight--;
	if (sbHighlight == -1) sbHighlight = inFile->beam[ibeam].nBand-1;
	if (select!=-1) {select = sbHighlight; recalc=1;}
      }
    else if (key=='f')
      {
	int nfit=4;
	double pval[nfit];
	int npts = ndumps;
	double px[npts];
	double py[npts];
	lm_status_struct status;
	lm_control_struct control = lm_control_double;
	int b;
	if (select==-1) select=0;
	printf("Band | Fc (MHz) | Amp | Angular offset | Width | Baseline\n");

	if (fitted==0 && normalise==-1)
	  {
	    for (b=0;b<inFile->beam[ibeam].nBand;b++)
	      {
		modelFit[b].baseline = fx[b][0];         // Have a more reasonable baseline if not normalised
		modelFit[b].amp = maxy-miny;
	      }
	  }
	
	for (b=0;b<inFile->beam[ibeam].nBand;b++)
	  {
	    for (i=0;i<npts;i++)
	      {
		px[i] = fx[b][i];
		py[i] = fy[b][i];
	      }
	    
	    pval[0] = modelFit[b].amp;
	    pval[1] = modelFit[b].offset;
	    pval[2] = modelFit[b].width;
	    pval[3] = modelFit[b].baseline;
	    
	    lmcurve_fit(nfit,pval,npts,px,py,nonlinearFunc,&control,&status);  
	    printf("%s %.4f %g %g %g %g\n",inFile->beam[ibeam].bandHeader[b].label,inFile->beam[ibeam].bandHeader[b].fc,pval[0],pval[1],pval[2],pval[3]);

	    modelFit[b].amp = pval[0];
	    modelFit[b].offset = pval[1];
	    modelFit[b].width = pval[2];
	    modelFit[b].baseline = pval[3];
	  }
	recalc=1;
	fitted=1;

      }     
    else if (key=='x')
      {
	xaxis++;
	if (xaxis==9)
	  xaxis=0;
	recalc=1;
      }
    else if (key=='z')
      {
	float mx2,my2;
	cpgband(2,0,mx,my,&mx2,&my2,&key);
	minx = mx;
	maxx = mx2;
	miny = my;
	maxy = my2;

      }
    else if (key=='u')
      recalc=1;
    else if (key=='s')
      {
	if (select==-1)
	  select = sbHighlight;
	else
	  select = -1;
	recalc=1;
      }
    else if (key=='p')
      {
	printf("Current (ra0,dec0) = (%.2f,%.2f)\n",ra0,dec0);
	printf("Enter new values:  ");
	scanf("%lf %lf",&ra0,&dec0);
	recalc=1;
      }
    else if (key=='n')
      {
	normalise*=-1;
	recalc=1;
      }
    else if (key=='o')
      {
	offset*=-1;
	recalc=1;
      }
    else if (key=='L')
      label*=-1;
    else if  (key=='l')
      {
	logy*=-1;
	recalc=1;
      }
  } while (key!='q');
    
  cpgend();


  sdhdf_closeFile(inFile);
  
  free(inFile);
  
} 

double haversine(double centre_long,double centre_lat,double src_long,double src_lat)
{
  double dlon,dlat,a,c;
  double deg2rad = M_PI/180.0;

  centre_long*=deg2rad;
  centre_lat*=deg2rad;
  src_long*=deg2rad;
  src_lat*=deg2rad;
  
  /* Apply the Haversine formula */
  dlon = (src_long - centre_long);
  dlat = (src_lat  - centre_lat);
  a = pow(sin(dlat/2.0),2) + cos(centre_lat) *
    cos(src_lat)*pow(sin(dlon/2.0),2);
  if (a==1)
    c = M_PI;
  else
    c = 2.0 * atan2(sqrt(a),sqrt(1.0-a));
  return c/deg2rad;
}
