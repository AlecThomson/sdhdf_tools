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

double haversine(double centre_long,double centre_lat,double src_long,double src_lat);

int main(int argc,char *argv[])
{
  int i,j;
  char fname[MAX_FILES][64];
  int nFiles=0;
  int ibeam=0;
  int iband=0;
  sdhdf_fileStruct *inFile;
  int useCoord=0;
  float ra0,dec0,raddist;
  double dist;
  int useSource=0;
  int info=0;
  char src[MAX_STRLEN];
  int exactMatch=1;
  
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
      else if (strcmp(argv[i],"-partial")==0)
	exactMatch=0;
      else if (strcmp(argv[i],"-info")==0)
	info=1;
      else if (strcmp(argv[i],"-src")==0)
	{
	  useSource=1;
	  strcpy(src,argv[++i]);
	}
      else
	{
	  strcpy(fname[nFiles++],argv[i]);
	}
    }

  for (i=0;i<nFiles;i++)
    {
      sdhdf_initialiseFile(inFile);
      sdhdf_openFile(fname[i],inFile,1);
      sdhdf_loadMetaData(inFile);
      if (useSource==1)
	{
	  if (exactMatch==1)
	    {
	      if (strcmp(inFile->beamHeader[ibeam].source,src)==0)
		{
		  if (info==1) printf("[SRC MATCH] %s %s\n",fname[i],inFile->beamHeader[ibeam].source);
		  else printf("%s\n",fname[i]);
		}
	    }
	  else
	    {
	      if (strstr(inFile->beamHeader[ibeam].source,src)!=NULL)
		{
		  if (info==1) printf("[SRC MATCH] %s %s\n",fname[i],inFile->beamHeader[ibeam].source);
		  else printf("%s\n",fname[i]);
		}
	    }
	      
	}
      if (useCoord==1)
	{
	  for (j=0;j<inFile->beam[ibeam].bandHeader[iband].ndump;j++)
	    {
	      dist = haversine(ra0,dec0,inFile->beam[ibeam].bandData[iband].astro_obsHeader[j].raDeg,
			       inFile->beam[ibeam].bandData[iband].astro_obsHeader[j].decDeg);
	      if (dist < raddist)
		{
		  if (info==1) printf("[DIST MATCH] %s %d %g\n",fname[i],j,dist);
		  else printf("%s\n",fname[i]);		    
		}
	    }
	}
      sdhdf_closeFile(inFile);
    }
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
