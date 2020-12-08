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
