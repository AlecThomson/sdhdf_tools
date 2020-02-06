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

// Multiple 4x4 matrix
static inline void sdhdf_mult4x4(float src1[4][4], float src2[4][4], float dest[4][4])
{
  dest[0][0] = src1[0][0] * src2[0][0] + src1[0][1] * src2[1][0] + src1[0][2] * src2[2][0] + src1[0][3] * src2[3][0];
  dest[0][1] = src1[0][0] * src2[0][1] + src1[0][1] * src2[1][1] + src1[0][2] * src2[2][1] + src1[0][3] * src2[3][1];
  dest[0][2] = src1[0][0] * src2[0][2] + src1[0][1] * src2[1][2] + src1[0][2] * src2[2][2] + src1[0][3] * src2[3][2];
  dest[0][3] = src1[0][0] * src2[0][3] + src1[0][1] * src2[1][3] + src1[0][2] * src2[2][3] + src1[0][3] * src2[3][3];
  dest[1][0] = src1[1][0] * src2[0][0] + src1[1][1] * src2[1][0] + src1[1][2] * src2[2][0] + src1[1][3] * src2[3][0];
  dest[1][1] = src1[1][0] * src2[0][1] + src1[1][1] * src2[1][1] + src1[1][2] * src2[2][1] + src1[1][3] * src2[3][1];
  dest[1][2] = src1[1][0] * src2[0][2] + src1[1][1] * src2[1][2] + src1[1][2] * src2[2][2] + src1[1][3] * src2[3][2];
  dest[1][3] = src1[1][0] * src2[0][3] + src1[1][1] * src2[1][3] + src1[1][2] * src2[2][3] + src1[1][3] * src2[3][3];
  dest[2][0] = src1[2][0] * src2[0][0] + src1[2][1] * src2[1][0] + src1[2][2] * src2[2][0] + src1[2][3] * src2[3][0];
  dest[2][1] = src1[2][0] * src2[0][1] + src1[2][1] * src2[1][1] + src1[2][2] * src2[2][1] + src1[2][3] * src2[3][1];
  dest[2][2] = src1[2][0] * src2[0][2] + src1[2][1] * src2[1][2] + src1[2][2] * src2[2][2] + src1[2][3] * src2[3][2];
  dest[2][3] = src1[2][0] * src2[0][3] + src1[2][1] * src2[1][3] + src1[2][2] * src2[2][3] + src1[2][3] * src2[3][3];
  dest[3][0] = src1[3][0] * src2[0][0] + src1[3][1] * src2[1][0] + src1[3][2] * src2[2][0] + src1[3][3] * src2[3][0];
  dest[3][1] = src1[3][0] * src2[0][1] + src1[3][1] * src2[1][1] + src1[3][2] * src2[2][1] + src1[3][3] * src2[3][1];
  dest[3][2] = src1[3][0] * src2[0][2] + src1[3][1] * src2[1][2] + src1[3][2] * src2[2][2] + src1[3][3] * src2[3][2];
  dest[3][3] = src1[3][0] * src2[0][3] + src1[3][1] * src2[1][3] + src1[3][2] * src2[2][3] + src1[3][3] * src2[3][3];
};


// Invert 4x4 matrix
int sdhdf_inv4x4(float m[4][4],float invOut[4][4])
{
  float det;
  float inv[4][4];
  int i,j;

  inv[0][0] = m[1][1]  * m[2][2] * m[3][3] -
    m[1][1]  * m[2][3] * m[3][2] -
    m[2][1]  * m[1][2]  * m[3][3] +
    m[2][1]  * m[1][3]  * m[3][2] +
    m[3][1]  * m[1][2]  * m[2][3] -
    m[3][1]  * m[1][3]  * m[2][2];
  
  inv[1][0] = -m[1][0]  * m[2][2] * m[3][3] +
    m[1][0]  * m[2][3] * m[3][2] +
    m[2][0]  * m[1][2]  * m[3][3] -
    m[2][0]  * m[1][3]  * m[3][2] -
    m[3][0] * m[1][2]  * m[2][3] +
    m[3][0] * m[1][3]  * m[2][2];
  
  inv[2][0] = m[1][0]  * m[2][1] * m[3][3] -
    m[1][0]  * m[2][3] * m[3][1] -
    m[2][0]  * m[1][1] * m[3][3] +
    m[2][0]  * m[1][3] * m[3][1] +
    m[3][0] * m[1][1] * m[2][3] -
    m[3][0] * m[1][3] * m[2][1];
  
  inv[3][0] = -m[1][0]  * m[2][1] * m[3][2] +
    m[1][0]  * m[2][2] * m[3][1] +
    m[2][0]  * m[1][1] * m[3][2] -
    m[2][0]  * m[1][2] * m[3][1] -
    m[3][0] * m[1][1] * m[2][2] +
    m[3][0] * m[1][2] * m[2][1];
  
  inv[0][1] = -m[0][1]  * m[2][2] * m[3][3] +
    m[0][1]  * m[2][3] * m[3][2] +
    m[2][1]  * m[0][2] * m[3][3] -
    m[2][1]  * m[0][3] * m[3][2] -
    m[3][1] * m[0][2] * m[2][3] +
    m[3][1] * m[0][3] * m[2][2];
  
  inv[1][1] = m[0][0]  * m[2][2] * m[3][3] -
    m[0][0]  * m[2][3] * m[3][2] -
    m[2][0]  * m[0][2] * m[3][3] +
    m[2][0]  * m[0][3] * m[3][2] +
    m[3][0] * m[0][2] * m[2][3] -
    m[3][0] * m[0][3] * m[2][2];
  
  inv[2][1] = -m[0][0]  * m[2][1] * m[3][3] +
    m[0][0]  * m[2][3] * m[3][1] +
    m[2][0]  * m[0][1] * m[3][3] -
    m[2][0]  * m[0][3] * m[3][1] -
    m[3][0] * m[0][1] * m[2][3] +
    m[3][0] * m[0][3] * m[2][1];
  
  inv[3][1] = m[0][0]  * m[2][1] * m[3][2] -
    m[0][0]  * m[2][2] * m[3][1] -
    m[2][0]  * m[0][1] * m[3][2] +
    m[2][0]  * m[0][2] * m[3][1] +
    m[3][0] * m[0][1] * m[2][2] -
    m[3][0] * m[0][2] * m[2][1];
  
  inv[0][2] = m[0][1]  * m[1][2] * m[3][3] -
    m[0][1]  * m[1][3] * m[3][2] -
    m[1][1]  * m[0][2] * m[3][3] +
    m[1][1]  * m[0][3] * m[3][2] +
    m[3][1] * m[0][2] * m[1][3] -
    m[3][1] * m[0][3] * m[1][2];
  
  inv[1][2] = -m[0][0]  * m[1][2] * m[3][3] +
    m[0][0]  * m[1][3] * m[3][2] +
    m[1][0]  * m[0][2] * m[3][3] -
    m[1][0]  * m[0][3] * m[3][2] -
    m[3][0] * m[0][2] * m[1][3] +
    m[3][0] * m[0][3] * m[1][2];
  
  inv[2][2] = m[0][0]  * m[1][1] * m[3][3] -
    m[0][0]  * m[1][3] * m[3][1] -
    m[1][0]  * m[0][1] * m[3][3] +
    m[1][0]  * m[0][3] * m[3][1] +
    m[3][0] * m[0][1] * m[1][3] -
    m[3][0] * m[0][3] * m[1][1];
  
  inv[3][2] = -m[0][0]  * m[1][1] * m[3][2] +
    m[0][0]  * m[1][2] * m[3][1] +
    m[1][0]  * m[0][1] * m[3][2] -
    m[1][0]  * m[0][2] * m[3][1] -
    m[3][0] * m[0][1] * m[1][2] +
    m[3][0] * m[0][2] * m[1][1];
  
  inv[0][3] = -m[0][1] * m[1][2] * m[2][3] +
    m[0][1] * m[1][3] * m[2][2] +
    m[1][1] * m[0][2] * m[2][3] -
    m[1][1] * m[0][3] * m[2][2] -
    m[2][1] * m[0][2] * m[1][3] +
    m[2][1] * m[0][3] * m[1][2];
  
  inv[1][3] = m[0][0] * m[1][2] * m[2][3] -
    m[0][0] * m[1][3] * m[2][2] -
    m[1][0] * m[0][2] * m[2][3] +
    m[1][0] * m[0][3] * m[2][2] +
    m[2][0] * m[0][2] * m[1][3] -
    m[2][0] * m[0][3] * m[1][2];
  
  inv[2][3] = -m[0][0] * m[1][1] * m[2][3] +
    m[0][0] * m[1][3] * m[2][1] +
    m[1][0] * m[0][1] * m[2][3] -
    m[1][0] * m[0][3] * m[2][1] -
    m[2][0] * m[0][1] * m[1][3] +
    m[2][0] * m[0][3] * m[1][1];
  
  inv[3][3] = m[0][0] * m[1][1] * m[2][2] -
    m[0][0] * m[1][2] * m[2][1] -
    m[1][0] * m[0][1] * m[2][2] +
    m[1][0] * m[0][2] * m[2][1] +
    m[2][0] * m[0][1] * m[1][2] -
    m[2][0] * m[0][2] * m[1][1];
  
  det = m[0][0] * inv[0][0] + m[0][1] * inv[1][0] + m[0][2] * inv[2][0] + m[0][3] * inv[3][0];
  //  printf("det = %g\n",det);
  if (det == 0)
    return -1;

  det = 1.0 / det;
  
  for (i = 0; i < 4; i++)
    {
      for (j=0;j<4;j++)
	invOut[i][j] = inv[i][j] * det;
    }
  return 0;
}

/* Computes the dot product of two vectors */
double sdhdf_dotproduct(double *v1,double *v2)
{
  double dot;
  dot = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
  return dot;
}

// Function based on spot.py 
// FIX ME: MUST CHECK THIS CAREFULLY
void sdhdf_para(double dxd,double ddc,double q,double *axd,double *eld)
{
  // Determine position of beam in (RA,DEC) frame
  // given its RA (dxd) and DEC (ddc) offsets, plus
  // parallactic angle (q).
  q = q * M_PI/180.0;
  *axd = -dxd*cos(q) + ddc*sin(q);
  *eld =  ddc*cos(q) + dxd*sin(q);
  
}

