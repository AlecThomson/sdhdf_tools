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

void checkLine(float f0,float f1,char *str,float refFreq);

int main(int argc,char *argv[])
{
  float f0,f1;

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
  checkLine(f0,f1,"CH",704.175);
  checkLine(f0,f1,"CH",722.303);
  checkLine(f0,f1,"CH",724.791);
  checkLine(f0,f1,"CH_3OH",834.285);
  checkLine(f0,f1,"CH_3CHO",1065.076);
  checkLine(f0,f1,"CH_2CHCN",1371.723);
  checkLine(f0,f1,"CH_2CHCN",1371.796);
  checkLine(f0,f1,"CH_2CHCN",1371.931);
  checkLine(f0,f1,"HI: Neutral hydrogen",1420.405752);
  checkLine(f0,f1,"NH_2CHO",1538.108);
  checkLine(f0,f1,"NH_2CHO",1538.676);
  checkLine(f0,f1,"NH_2CHO",1539.264);
  checkLine(f0,f1,"NH_2CHO",1539.527);
  checkLine(f0,f1,"NH_2CHO",1539.832);
  checkLine(f0,f1,"NH_2CHO",1540.998);
  checkLine(f0,f1,"NH_2^13CHO",1570.805);
  checkLine(f0,f1,"^18OH",1584.274);
  checkLine(f0,f1,"CH_3OCHO",1610.247);
  checkLine(f0,f1,"CH_3OCHO",1610.900);
  checkLine(f0,f1,"OH: Hydroxyl radical",1612.2310);
  checkLine(f0,f1,"^17OH",1624.518);
  checkLine(f0,f1,"^17OH",1626.161);
  checkLine(f0,f1,"^18OH",1637.564);
  checkLine(f0,f1,"HCOOH",1638.805);
  checkLine(f0,f1,"^18OH",1639.503);  
  checkLine(f0,f1,"OH: Hydroxyl radical",1665.4018);
  checkLine(f0,f1,"OH: Hydroxyl radical",1667.3590);
  checkLine(f0,f1,"^18OH",1692.795);  
  checkLine(f0,f1,"OH: Hydroxyl radical",1720.5300);
  checkLine(f0,f1,"CH_3OCHO",1849.634);
  checkLine(f0,f1,"CNCHO",2078.068);
  checkLine(f0,f1,"HC_5N",2661.604);
  checkLine(f0,f1,"HC_5N",2662.877);
  checkLine(f0,f1,"HC_5N",2664.786);
  checkLine(f0,f1,"H_2CS",3139.405);
  checkLine(f0,f1,"CH_3CHO",3195.162);
  checkLine(f0,f1,"CH",3263.794);
  checkLine(f0,f1,"CH",3335.481);
  checkLine(f0,f1,"CH",3349.193);
  
  

}

// Should include number of decimal places to print
void checkLine(float f0,float f1,char *str,float refFreq)
{
  if (refFreq >= f0 && refFreq <= f1)
    printf(" %-15.4f %-20.20s\n",refFreq,str);
}
