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
// Software to provide a quick-look at the metadata in a SDHDF file.
//
// Usage:
// sdhdf_identify <filename.hdf> <filename2.hdf> ...
//
// Compilation
// gcc -lm -o sdhdf_identify sdhdf_identify.c sdhdfProc.c -I../hdf5-1.10.4/src/ ../hdf5-1.10.4/src/.libs/libhdf5.a -ldl -lz -L/pulsar/psr/software/stable/stretch/lib/ -I//pulsar/psr/software/stable/stretch/include -lcalceph -Isofa/20190722/c/src/ -Lsofa/20190722/c/src -lsofa_c
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sdhdfProc.h"
#include "hdf5.h"

#define MAX_SRC 128

typedef struct infoStruct {
  char fname[MAX_STRLEN];
  char source[MAX_STRLEN];
  float ra;
  float dec;
  double mjd;
} infoStruct;

double haversine(double centre_long,double centre_lat,double src_long,double src_lat);

int main(int argc,char *argv[])
{
  int i,j,k,l;
  char fname[MAX_FILES][64];
  int nFiles=0;
  int ibeam=0;
  int iband=0;
  infoStruct *info;
  sdhdf_fileStruct *inFile;
  int useCoord=0;
  float ra0,dec0,raddist;
  double dist;
  int useSource=0;
  int infoDisp=0;
  char **src;
  char **onSrc;
  char **offSrc;
  int nSrc=0;
  int nOnSrc=0;
  int nOffSrc=0;
  int exactMatch=1;
  int openFile=0;
  int pair=0;
  char outFile[MAX_STRLEN];
  int writeOut=0;
  FILE *fout;
  
  src = (char **)malloc(sizeof(char *)*MAX_SRC);
  onSrc = (char **)malloc(sizeof(char *)*MAX_SRC);
  offSrc = (char **)malloc(sizeof(char *)*MAX_SRC);
  for (i=0;i<MAX_SRC;i++)
    {
      src[i] = (char *)malloc(sizeof(char)*MAX_STRLEN);
      onSrc[i] = (char *)malloc(sizeof(char)*MAX_STRLEN);
      offSrc[i] = (char *)malloc(sizeof(char)*MAX_STRLEN);
    }
  
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }

  for (i=1;i<argc;i++)
    {      
      if (strcmp(argv[i],"-coord")==0)
	{
	  useCoord=1;
	  sscanf(argv[++i],"%f",&ra0);
	  sscanf(argv[++i],"%f",&dec0);
	  sscanf(argv[++i],"%f",&raddist);
	}
      else if (strcmp(argv[i],"-o")==0)
	{
	  strcpy(outFile,argv[++i]);
	  writeOut=1;
	}
      else if (strcmp(argv[i],"-pair")==0)
	pair=1;
      else if (strcmp(argv[i],"-partial")==0)
	exactMatch=0;
      else if (strcmp(argv[i],"-info")==0)
	infoDisp=1;
      else if (strcmp(argv[i],"-src")==0)
	{
	  useSource=1;
	  strcpy(src[nSrc++],argv[++i]);
	}
      else if (strcmp(argv[i],"-on")==0)
	strcpy(onSrc[nOnSrc++],argv[++i]);
      else if (strcmp(argv[i],"-off")==0)
	strcpy(offSrc[nOffSrc++],argv[++i]);
      else
	{
	  strcpy(fname[nFiles++],argv[i]);
	}
    }


  info = (infoStruct *)malloc(sizeof(infoStruct)*nFiles);

  for (i=0;i<nFiles;i++)
    {
      sdhdf_initialiseFile(inFile);
      openFile = sdhdf_openFile(fname[i],inFile,1);
      if (openFile==-1)
	{
	  printf("ERROR: Unable to open file %s\n",fname[i]);
	}
      else
	{
	  sdhdf_loadMetaData(inFile);
	  strcpy(info[i].fname,fname[i]);
	  strcpy(info[i].source,inFile->beamHeader[ibeam].source);
	  info[i].ra = inFile->beam[ibeam].bandData[iband].astro_obsHeader[0].raDeg; // Note 0 here
	  info[i].dec = inFile->beam[ibeam].bandData[iband].astro_obsHeader[0].decDeg;
	  info[i].mjd = inFile->beam[ibeam].bandData[iband].astro_obsHeader[0].mjd;
	}
      sdhdf_closeFile(inFile);
    }
  printf("Loaded data for %d files\n",nFiles);


  if (writeOut==1)
    fout = fopen(outFile,"w");
    
  if (pair==1)
    {
      int foundOn=0;
      int foundOff=0;
      int bestOff;
      double timeDiff;
      for (i=0;i<nFiles;i++)
	{
	  // Find the on source
	  foundOn=0;
	  for (j=0;j<nOnSrc;j++)
	    {
	      if (strcmp(onSrc[j],info[i].source)==0)
		{foundOn=1; break;}
	    }
	  if (foundOn==1)
	    {
	      // Now find a matching off source
	      bestOff=-1;
	      for (j=0;j<nFiles;j++)
		{
		  foundOff=0;

		  for (k=0;k<nOffSrc;k++)
		    {
		      if (strcmp(offSrc[k],info[j].source)==0)
			{foundOff=1; break;}
		    }
		  if (foundOff==1)
		    {
		      if (bestOff == -1)
			{
			  bestOff = j;
			  timeDiff = fabs(info[j].mjd-info[i].mjd);
			}
		      else if (fabs(info[j].mjd-info[i].mjd) < timeDiff)
			{
			  bestOff = j;
			  timeDiff = fabs(info[j].mjd-info[i].mjd);
			}
		    }
		}
	      printf("[MATCH] %s %s # SRC_ON: %s SRC_OFF: %s TIMEDIFF: %g\n",info[i].fname,info[bestOff].fname,info[i].source,info[bestOff].source,timeDiff);
	      if (writeOut==1)
		fprintf(fout,"%s %s # SRC_ON: %s SRC_OFF: %s TIMEDIFF: %g\n",info[i].fname,info[bestOff].fname,info[i].source,info[bestOff].source,timeDiff);
	      //	      printf("Found %d as off source. Time difference = %g\n",bestOff,timeDiff);
	    }
	}
    }
  else
    {  
      for (i=0;i<nFiles;i++)
	{
	  if (useSource==1)
	    {
	      for (j=0;j<nSrc;j++)
		{
		  if (exactMatch==1)
		    {
		      if (strcmp(info[i].source,src[j])==0)
			{
			  if (infoDisp==1) printf("[SRC MATCH] %s %s\n",fname[i],info[i].source);
			  else printf("%s\n",fname[i]);
			}
		    }
		  else
		    {
		      if (strstr(info[i].source,src[j])!=NULL)
			{
			  if (infoDisp==1) printf("[SRC MATCH] %s %s\n",fname[i],info[i].source);
			  else printf("%s\n",fname[i]);
			}
		    }
		}
	      
	    }
	  if (useCoord==1)
	    {
	      printf("CURRENTLY NOT CHECKING COORDINATES IN SDHDF_IDENTIFY -- FIX ME\n");
	      /*  UPDATE TO USE INFO STUCT
		  for (j=0;j<inFile->beam[ibeam].bandHeader[iband].ndump;j++)
		  {
		  dist = haversine(ra0,dec0,inFile->beam[ibeam].bandData[iband].astro_obsHeader[j].raDeg,
		  inFile->beam[ibeam].bandData[iband].astro_obsHeader[j].decDeg);
		  if (dist < raddist)
		  {
		  if (infoDisp==1) printf("[DIST MATCH] %s %d %g\n",fname[i],j,dist);
		  else printf("%s\n",fname[i]);		    
		  }
		  }
	      */
	    }
	}
    }
  if (writeOut==1)
    fclose(fout);
    
  free(inFile);
  free(info);
  for (i=0;i<MAX_SRC;i++)
    {
      free(src[i]);
      free(onSrc[i]);
      free(offSrc[i]);
    }
  free(src);
  free(onSrc);
  free(offSrc);

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
