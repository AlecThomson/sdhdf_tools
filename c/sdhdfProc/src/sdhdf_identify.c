//  Copyright (C) 2019, 2020, 2021 George Hobbs

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
#define MAX_BANDS 26

typedef struct infoStruct {
  char fname[MAX_STRLEN];
  char source[MAX_STRLEN];
  float ra;
  float dec;
  double mjd;
  int nband;
  int sb;
} infoStruct;

double haversine(double centre_long,double centre_lat,double src_long,double src_lat);
double dms_turn(char *line);
double hms_turn(char *line);
int turn_hms(double turn, char *hms);
int turn_dms(double turn, char *dms);
double turn_deg(double turn);


void help()
{
  printf("sdhdf_identify\n\n");
  printf("Command line options\n\n");
  printf("-allBand                        - look for matches amongst all sub-bands in the files (otherwise use band defined by -sb)\n");
  printf("-coord <ra0> <dec0> <raddist>   - selects file within a particular angular radius of those coordinates (RA = HH:MM:SS.SSS, DEC = DD:MM:SS.SSS, raddist = cone radius in degrees\n");
  printf("-h                              - this help\n");
  printf("-info                           - display more information output lines\n");
  printf("-o <outFile>                    - write the output to an output file\n");
  printf("-pair                           - produce paired 'on' and 'off' output - need -on and -off set\n");
  printf("-partial                        - does not require an exact match when checking source names\n");
  printf("-sb <subband>                   - select sub-band number \n");
  printf("-src <name>                     - selects based on source name\n");
  printf("-nband <num>                    - selects based on the number of sub-bands in the file\n");
  printf("-on <source name>               - defines 'on' source to be the specified name\n");
  printf("-off <source name>              - defines 'off' source to be the specified name\n");
  printf("\n\n");
  printf("Example: sdhdf_identify -pair -on 0407-658 -off 0407-658_R *.hdf\n");

}

int main(int argc,char *argv[])
{
  int i,j,k,l;
  char fname[MAX_FILES][64];
  int nFiles=0;
  int ibeam=0;
  int iband=0;
  int allBand=0;
  infoStruct *info;
  sdhdf_fileStruct *inFile;
  int useCoord=0;
  char ra0Str[1024],dec0Str[1024];
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
  int useNband=0,selNband;
  int writeOut=0;
  FILE *fout;
  int disp;
  int nEntry=0;
  
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
	  sscanf(argv[++i],"%s",ra0Str);
	  sscanf(argv[++i],"%s",dec0Str);
	  sscanf(argv[++i],"%f",&raddist);
	}
      else if (strcasecmp(argv[i],"-allband")==0)
	allBand=1;
      else if (strcmp(argv[i],"-sb")==0)
	sscanf(argv[++i],"%d",&iband);
      else if (strcmp(argv[i],"-h")==0)
	help();
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
      else if (strcmp(argv[i],"-nband")==0)
	{
	  useNband = 1;
	  sscanf(argv[++i],"%d",&selNband);
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


  info = (infoStruct *)malloc(sizeof(infoStruct)*nFiles*MAX_BANDS); // SHOULD FIX THE MAX_BANDS HERE

  nEntry=0;
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
	  if (allBand == 0)
	    {
	      strcpy(info[nEntry].fname,fname[i]);
	      strcpy(info[nEntry].source,inFile->beamHeader[ibeam].source);
	      info[nEntry].ra = inFile->beam[ibeam].bandData[iband].astro_obsHeader[0].raDeg; // Note 0 here
	      info[nEntry].dec = inFile->beam[ibeam].bandData[iband].astro_obsHeader[0].decDeg;
	      info[nEntry].mjd = inFile->beam[ibeam].bandData[iband].astro_obsHeader[0].mjd;
	      info[nEntry].nband = inFile->beam[ibeam].nBand;
	      info[nEntry].sb = iband;
	      nEntry++;
	    }
	  else
	    {
	      for (j=0;j<inFile->beam[ibeam].nBand;j++)
		{
		  strcpy(info[nEntry].fname,fname[i]);
		  strcpy(info[nEntry].source,inFile->beamHeader[ibeam].source);
		  info[nEntry].ra = inFile->beam[ibeam].bandData[j].astro_obsHeader[0].raDeg; // Note 0 here
		  info[nEntry].dec = inFile->beam[ibeam].bandData[j].astro_obsHeader[0].decDeg;
		  info[nEntry].mjd = inFile->beam[ibeam].bandData[j].astro_obsHeader[0].mjd;
		  info[nEntry].nband = inFile->beam[ibeam].nBand;
		  info[nEntry].sb = j;
		  nEntry++;
		}
	    }
	}
      sdhdf_closeFile(inFile);
    }
  //    printf("Loaded data for %d files\n",nFiles);


  if (writeOut==1)
    fout = fopen(outFile,"w");
    
  if (pair==1)
    {
      int foundOn=0;
      int foundOff=0;
      int bestOff;
      double timeDiff;
      for (i=0;i<nEntry;i++)
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
	      for (j=0;j<nEntry;j++)
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
		fprintf(fout,"%s %s %s.%s.onoff # SRC_ON: %s SRC_OFF: %s TIMEDIFF: %g\n",info[i].fname,info[bestOff].fname,info[i].fname,info[bestOff].fname,info[i].source,info[bestOff].source,timeDiff);
	      //	      printf("Found %d as off source. Time difference = %g\n",bestOff,timeDiff);
	    }
	}
    }
  else
    {
      for (i=0;i<nEntry;i++)
	{
	  // Check source names
	  if (useSource==1)
	    {
	      for (j=0;j<nSrc;j++)
		{
		  disp=0;
		  if (exactMatch==1)
		    {
		      if (strcmp(info[i].source,src[j])==0) disp=1;			
		    }
		  else
		    {
		      if (strstr(info[i].source,src[j])!=NULL) disp=1;
		    }
		  if (disp==1)
		    {
		      if (infoDisp==1) printf("[SRC MATCH] %s %s sub-band %d\n",info[i].fname,info[i].source,info[i].sb);
		      else printf("%s\n",info[i].fname);
		    }
		}
	      
	    }
	  // Check coordinates
	  if (useCoord==1)
	    {	      
	      ra0 = turn_deg(hms_turn(ra0Str));
	      dec0 = turn_deg(dms_turn(dec0Str));
	      dist = haversine(ra0,dec0,info[i].ra,info[i].dec);
	      
	      if (dist < raddist)
		{
		  if (infoDisp==1) printf("[DIST MATCH] %s %g %g %g sub-band %d\n",info[i].fname,info[i].ra,info[i].dec,dist,info[i].sb);
		  else printf("%s\n",info[i].fname);		    
		}
	    }

	  // Check nband
	  if (useNband == 1)
	    {
	      if (info[i].nband == selNband)
		{
		  if (infoDisp==1) printf("[NBAND MATCH] %s sub-band %d\n",fname[i],info[i].sb);
		  else
		    printf("%s\n",fname[i]);
		}
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



// Input in degrees
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


double turn_deg(double turn){
 
  /* Converts double turn to string "sddd.ddd" */
  return turn*360.0;
}


int turn_dms(double turn, char *dms){
  
  /* Converts double turn to string "sddd:mm:ss.sss" */
  
  int dd, mm, isec;
  double trn, sec;
  char sign;
  
  sign=' ';
  if (turn < 0.){
    sign = '-';
    trn = -turn;
  }
  else{
    sign = '+';
    trn = turn;
  }
  dd = trn*360.;
  mm = (trn*360.-dd)*60.;
  sec = ((trn*360.-dd)*60.-mm)*60.;
  isec = (sec*1000. +0.5)/1000;
    if(isec==60){
      sec=0.;
      mm=mm+1;
      if(mm==60){
        mm=0;
        dd=dd+1;
      }
    }
  sprintf(dms,"%c%02d:%02d:%010.7f",sign,dd,mm,sec);
 
}


int turn_hms(double turn, char *hms){
 
  /* Converts double turn to string " hh:mm:ss.ssss" */
  
  int hh, mm, isec;
  double sec;

  hh = turn*24.;
  mm = (turn*24.-hh)*60.;
  sec = ((turn*24.-hh)*60.-mm)*60.;
  isec = (sec*10000. +0.5)/10000;
    if(isec==60){
      sec=0.;
      mm=mm+1;
      if(mm==60){
        mm=0;
        hh=hh+1;
        if(hh==24){
          hh=0;
        }
      }
    }

  sprintf(hms," %02d:%02d:%010.7f",hh,mm,sec);
 
}


double hms_turn(char *line){

  /* Converts string " hh:mm:ss.ss" or " hh mm ss.ss" to double turn */
  
  int i;int turn_hms(double turn, char *hms);
  double hr, min, sec, turn=0;
  char hold[MAX_STRLEN];

  strcpy(hold,line);

  /* Get rid of ":" */
  for(i=0; *(line+i) != '\0'; i++)if(*(line+i) == ':')*(line+i) = ' ';

  i = sscanf(line,"%lf %lf %lf", &hr, &min, &sec);
  if(i > 0){
    turn = hr/24.;
    if(i > 1)turn += min/1440.;
    if(i > 2)turn += sec/86400.;
  }
  if(i == 0 || i > 3)turn = 1.0;


  strcpy(line,hold);

  return turn;
}

double dms_turn(char *line){

  /* Converts string "-dd:mm:ss.ss" or " -dd mm ss.ss" to double turn */
  
  int i;
  char *ic, ln[40];
  double deg, min, sec, sign, turn=0;

  /* Copy line to internal string */
  strcpy(ln,line);

  /* Get rid of ":" */
  for(i=0; *(ln+i) != '\0'; i++)if(*(ln+i) == ':')*(ln+i) = ' ';

  /* Get sign */
  if((ic = strchr(ln,'-')) == NULL)
     sign = 1.;
  else {
     *ic = ' ';
     sign = -1.;
  }

  /* Get value */
  i = sscanf(ln,"%lf %lf %lf", &deg, &min, &sec);
  if(i > 0){
    turn = deg/360.;
    if(i > 1)turn += min/21600.;
    if(i > 2)turn += sec/1296000.;
    if(turn >= 1.0)turn = turn - 1.0;
    turn *= sign;
  }
  if(i == 0 || i > 3)turn =1.0;

  return turn;
}
