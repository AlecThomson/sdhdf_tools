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
#include <string.h>
#include <stdlib.h>
#include "sdhdfProc.h"

#define MAX_REST_FREQUENCIES 1024

void checkLine(float f0,float f1,char *str,float refFreq);

int main(int argc,char *argv[])
{
  float f0,f1;
  int nFreq=0;
  sdhdf_restfrequency_struct *restFrequencies;
  int i;
  
  restFrequencies = (sdhdf_restfrequency_struct *)malloc(sizeof(sdhdf_restfrequency_struct)*MAX_REST_FREQUENCIES);
  sdhdf_loadRestFrequencies(restFrequencies,&nFreq);


  if (argc == 3)
    {
      sscanf(argv[1],"%f",&f0);
      sscanf(argv[2],"%f",&f1);
    }
  else
    {
      f0 = 704;
      f1 = 4032;
    }
  for (i=0;i<nFreq;i++)
    checkLine(f0,f1,restFrequencies[i].label,restFrequencies[i].f0);
  
  free(restFrequencies);
}

// Should include number of decimal places to print
void checkLine(float f0,float f1,char *str,float refFreq)
{
  if (refFreq >= f0 && refFreq <= f1)
    printf(" %-15.4f %-30.30s\n",refFreq,str);
}
