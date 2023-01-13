//  Copyright (C) 2021, 2022 George Hobbs

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
#include <string.h>
#include <math.h>
#include "sdhdfProc.h"
#include "hdf5.h"
#include <cpgplot.h>


// Data modelling structure
typedef struct dataModelStruct {
  char bandLabel[MAX_STRLEN];
  int original_nsub;
  int original_nchan;
  int original_npol;
  float *freq;
  float *originalData;
  float *model;
  float *wts;
} dataModelStruct;


//
// Spline interpolation
//
void TKspline_interpolate(int n,float *x,float **yd,float *interpX,
			  float *interpY,int nInterp);
void TKcmonot (int n, float *x, float *y, float **yd);
void writeModel(sdhdf_fileStruct *outFile,dataModelStruct *model,int nband);

// Smoothing spline algorithms
void SPLINE(int N, float *Z, float *ZF, float ZSP, float ZPV);
void **MATRIX(int nrows, int ncols, int first_row, int first_col,int element_size);
void errorAction(int N, double *Y, float *ZF);
void modelSpectrum(dataModelStruct *model,int nband,sdhdf_fileStruct *outFile);

void help()
{
  printf("sdhdf_model: routine to model spectra (and baselines)\n");
  printf("\n\n");
  printf("Command line arguments:\n\n");
  printf("-f <filename>     Input filename\n");
  printf("-sb <band>        Input band number\n");
  printf("\n\n");
  printf("The output is an interactive display. The following commands are available\n\n");
  printf(">                 Move to next band\n");
  printf("<                 Move to previous band\n");
  printf("d                 Toggle plotting the raw data\n");
  printf("h                 This help\n");
  printf("l                 Define spline length\n");
  printf("m                 Toggle plotting the model\n");
  printf("u                 Unzoom\n");
  printf("q                 Quit\n");
  printf("s                 Write a new file containing the model\n");
  printf("v                 Define spline variance\n");
  printf("z                 Zoom into a region\n");
    
}

int main(int argc,char *argv[])
{
  int i,j,k;
  char fname[MAX_STRLEN]="unset";
  sdhdf_fileStruct *inFile;
  dataModelStruct *model;
  
  int ibeam=0;
  int idump,iband;
  int dataSize;
  int nband;
  int nchan,nsub,npol;
  char oname[1024]="model.hdf";
  // Open the output file
  sdhdf_fileStruct *outFile;

  idump = iband = 0;

  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }
    if (!(outFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for the output file\n");
      exit(1);
    }
  sdhdf_initialiseFile(outFile);

  
  for (i=1;i<argc;i++)
    {      
      if (strcmp(argv[i],"-f")==0)   	     strcpy(fname,argv[++i]);	
      else if (strcmp(argv[i],"-sb")==0) sscanf(argv[++i],"%d",&iband);
      else if (strcmp(argv[i],"-h")==0) {help(); exit(1);}
    }
  sdhdf_openFile(oname,outFile,3);
  sdhdf_initialiseFile(inFile);

      
  if (strcmp(fname,"unset")==0)
    {
      printf("Please define an input file name using the -f option\n");
      free(inFile);
      exit(1);
    }
  if (sdhdf_openFile(fname,inFile,1)==-1)
    {
      printf("Unable to open input file >%s<\n",fname);
      free(inFile);
      exit(1);
    }

  sdhdf_loadMetaData(inFile);

  // Allocate memory
  nband = inFile->beam[ibeam].nBand;
  model = (dataModelStruct *)malloc(sizeof(dataModelStruct)*nband);
  for (i=0;i<nband;i++) // Fix ME -- hardcode beam 0
    {
      nchan = inFile->beam[ibeam].bandHeader[i].nchan;
      nsub  = inFile->beam[ibeam].bandHeader[i].ndump;
      npol  = inFile->beam[ibeam].bandHeader[i].npol;
      strcpy(model[i].bandLabel,inFile->beam[ibeam].bandHeader[i].label);
      model[i].original_nsub  = nsub;
      model[i].original_nchan = nchan;
      model[i].original_npol  = npol;
      dataSize = nchan*nsub*npol;
      model[i].freq = (float *)malloc(sizeof(float)*nchan);
      model[i].originalData = (float *)malloc(sizeof(float)*dataSize);
      model[i].model = (float *)malloc(sizeof(float)*dataSize);
      model[i].wts = (float *)malloc(sizeof(float)*nchan*nsub);

      sdhdf_loadBandData(inFile,ibeam,i,1);

      // SHOULD CHECK FOR WEIGHTINGS       
      for (k=0;k<nsub;k++)
	{
	  for (j=0;j<nchan;j++)
	    {
	      // FIX ME: using [0] for frequency dump
	      model[i].freq[j]         = inFile->beam[ibeam].bandData[i].astro_data.freq[j];
	      model[i].model[k*nchan*npol+j] = model[i].originalData[k*nchan*npol+j] = inFile->beam[ibeam].bandData[i].astro_data.pol1[j+k*nchan];
	      model[i].model[k*nchan*npol+nchan + j] = model[i].originalData[k*nchan*npol+nchan+j] = inFile->beam[ibeam].bandData[i].astro_data.pol2[j+k*nchan];
	      // SHOULD GET POL 2 and POL 3	      

	    }
	}            
      // COULD RELEASE DATA HERE
    }  
  sdhdf_copyRemainder(inFile,outFile,0);
  sdhdf_closeFile(inFile);
  free(inFile);

  
  // Do the modelling
  modelSpectrum(model,nband,outFile);


  for (i=0;i<nband;i++)
    {
      free(model[i].freq);
      free(model[i].originalData);
      free(model[i].model);
      free(model[i].wts);
    }

  free(model);
  sdhdf_closeFile(outFile);

}


void modelSpectrum(dataModelStruct *model,int nband,sdhdf_fileStruct *outFile) // sdhdf_fileStruct *inFile,int ibeam,int iband,int idump)
{
  float *fx;
  float *fy1,*fy2,*fy3,*fy4;
  float *fy1_model,*fy2_model;
  float minx,maxx,miny,maxy;
  int nplot,nchan;
  int i,j,k;
  float mx,my;
  char key;
  int   plotData=1;
  int   plotModel=1;

  int recalcMinMax=1;
  int iband=0;
  
  float splineLength = 100; // Length (rigidity) of spline to be used to model series   
  float splineVariance = 0.5;  // Portion of variance at wavelength ZSP contained in spline (0<ZPV<1)

  int maxNchan=0;
  // Produce the data for the plot

  nplot=0;
  for (i=0;i<nband;i++)
    {
      if (maxNchan < model[i].original_nchan) maxNchan = model[i].original_nchan;
    }

  
  
  fx  = (float *)malloc(sizeof(float)*maxNchan);
  fy1 = (float *)malloc(sizeof(float)*maxNchan);
  fy1_model = (float *)malloc(sizeof(float)*maxNchan);
  fy2 = (float *)malloc(sizeof(float)*maxNchan);
  fy2_model = (float *)malloc(sizeof(float)*maxNchan);
  fy3 = (float *)malloc(sizeof(float)*maxNchan);
  fy4 = (float *)malloc(sizeof(float)*maxNchan);
    
  // Make the plot
  cpgbeg(0,"/xs",1,1);
  cpgask(0);
  do {
    nplot=0;
    for (i=0;i<model[iband].original_nchan;i++)
      {
	// SHOULD CHECK FOR WEIGHTINGS       
	// CHANGE HERE

	fx[nplot]  = model[iband].freq[i];
	fy1[nplot] = model[iband].originalData[i];  // FIX ME -- GET DUMPS IN HERE
	fy2[nplot] = model[iband].originalData[i+model[iband].original_nchan];
	// SHOULD GET POL 2 and POL 3
	if (recalcMinMax==1)
	  {
	    if (nplot==0)
	      {
		minx=maxx = fx[nplot];
		if (fy1[nplot] > fy2[nplot])
		  {
		    miny = fy2[nplot]; maxy = fy1[nplot];
		  }
		else
		  {
		    miny = fy1[nplot]; maxy = fy2[nplot];
		  }
	      }
	    else
	      {
		if (minx > fx[nplot]) minx = fx[nplot];
		if (maxx < fx[nplot]) maxx = fx[nplot];
		if (miny > fy1[nplot]) miny = fy1[nplot];
		if (maxy < fy1[nplot]) maxy = fy1[nplot];
		if (miny > fy2[nplot]) miny = fy2[nplot];
		if (maxy < fy2[nplot]) maxy = fy2[nplot];
	      }
	  }
	nplot++;
      }
    recalcMinMax=0;
	        
    // Interpolate using smoothing spline
    // Note that this requires regularly sampled data -- or needs to be split up between zapped regions
    // FIX AND CHECK HERE
    
    for (i=0;i<nband;i++)
      {
	SPLINE(model[i].original_nchan,&(model[i].originalData[0]),&(model[i].model[0]),splineLength,splineVariance);
	SPLINE(model[i].original_nchan,&(model[i].originalData[model[i].original_nchan]),
	       &(model[i].model[model[i].original_nchan]),splineLength,splineVariance);
      }
    for (j=0;j<model[iband].original_nchan;j++)
      {
	fy1_model[j] = model[iband].model[j];
	fy2_model[j] = model[iband].model[j+model[iband].original_nchan];
	//	printf("model = %g %g  (%d %d %d)\n",fy1_model[j],fy2_model[j],i,model[i].original_nchan,maxNchan);
      } 
    
    cpgenv(minx,maxx,miny,maxy,0,1);
    cpglab("Frequency (MHz)","Signal strength","");
    if (plotData==1)
      {
	cpgsci(2); cpgline(nplot,fx,fy1);
	cpgsci(4); cpgline(nplot,fx,fy2);
      }
    if (plotModel==1)
      {      
	cpgslw(3); cpgsci(1); cpgline(model[iband].original_nchan,fx,fy1_model); cpgslw(1);
	cpgslw(3); cpgsci(1); cpgline(model[iband].original_nchan,fx,fy2_model); cpgslw(1);
      }
    cpgsci(1);
    cpgcurs(&mx,&my,&key);
    //    
    // Should have an option to split the plot into each polarisation
    //
    if (key=='d') plotData*=-1;
    else if (key=='s') writeModel(outFile,model,nband);
    else if (key=='m') plotModel*=-1;
    else if (key=='l') {
      printf("Current spline length (rigidity) = %f. Enter new value: ",splineLength);
      scanf("%f",&splineLength);
    }
    else if (key=='v') {
      printf("Current spline variance = %f. Enter new value (0 < v < 1): ",splineVariance);
      scanf("%f",&splineVariance);
    }
    else if (key=='u')
      recalcMinMax=1;
    else if (key=='h')
      help();
    else if (key=='>')
      {
	if (iband < nband-1)
	  {
	    iband++; recalcMinMax=1;
	  }
      }
    else if (key=='<')
      {
	if (iband > 0)
	  {
	    iband--; recalcMinMax=1;
	  }
      }
    else if (key=='z') {
      float mx2,my2;
      cpgband(2,0,mx,my,&mx2,&my2,&key);
      if (mx != mx2 && my != my2)
	{
	  if (mx < mx2)
	    {minx = mx; maxx = mx2;}
	  else
	    {minx = mx2; maxx = mx;}
	  
	  if (my < my2)
	    {miny = my; maxy = my2;}
	  else
	    {miny = my2; maxy = my;}
	}
      else
	printf("Please press 'z' and then move somewhere and click left mouse button\n");
    }
  } while (key != 'q');
  cpgend();
  
  free(fx); free(fy1); free(fy2); free(fy3); free(fy4);
  free(fy1_model); free(fy2_model);
}

/* ************************************************************** *
 * Spline routines                                                *
 * ************************************************************** */
#define ABS(x) ((x) < 0 ? -(x) : (x))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))


void TKcmonot (int n, float *x, float *y, float **yd)
/*****************************************************************************
compute cubic interpolation coefficients via the Fritsch-Carlson method,
which preserves monotonicity
******************************************************************************
Input:
n		number of samples
x  		array[n] of monotonically increasing or decreasing abscissae
y		array[n] of ordinates

Output:
yd		array[n][4] of cubic interpolation coefficients (see notes)
******************************************************************************
Notes:
The computed cubic spline coefficients are as follows:
yd[i][0] = y(x[i])    (the value of y at x = x[i])
yd[i][1] = y'(x[i])   (the 1st derivative of y at x = x[i])
yd[i][2] = y''(x[i])  (the 2nd derivative of y at x = x[i])
yd[i][3] = y'''(x[i]) (the 3rd derivative of y at x = x[i])

To evaluate y(x) for x between x[i] and x[i+1] and h = x-x[i],
use the computed coefficients as follows:
y(x) = yd[i][0]+h*(yd[i][1]+h*(yd[i][2]/2.0+h*yd[i][3]/6.0))

The Fritsch-Carlson method yields continuous 1st derivatives, but 2nd
and 3rd derivatives are discontinuous.  The method will yield a
monotonic interpolant for monotonic data.  1st derivatives are set
to zero wherever first divided differences change sign.

For more information, see Fritsch, F. N., and Carlson, R. E., 1980,
Monotone piecewise cubic interpolation:  SIAM J. Numer. Anal., v. 17,
n. 2, p. 238-246.

Also, see the book by Kahaner, D., Moler, C., and Nash, S., 1989,
Numerical Methods and Software, Prentice Hall.  This function was
derived from SUBROUTINE PCHEZ contained on the diskette that comes
with the book.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 09/30/89
Modified:  Dave Hale, Colorado School of Mines, 02/28/91
	changed to work for n=1.
Modified:  Dave Hale, Colorado School of Mines, 08/04/91
	fixed bug in computation of left end derivative
*****************************************************************************/
{
  int i;
  float h1,h2,del1,del2,dmin,dmax,hsum,hsum3,w1,w2,drat1,drat2,divdf3;
  
  /* copy ordinates into output array */
  for (i=0; i<n; i++)
    yd[i][0] = y[i];
  
  /* if n=1, then use constant interpolation */
  if (n==1) {
    yd[0][1] = 0.0;
    yd[0][2] = 0.0;
    yd[0][3] = 0.0;
    return;
    
    /* else, if n=2, then use linear interpolation */
  } else if (n==2) {
    yd[0][1] = yd[1][1] = (y[1]-y[0])/(x[1]-x[0]);
    yd[0][2] = yd[1][2] = 0.0;
    yd[0][3] = yd[1][3] = 0.0;
    return;
  }
  
  /* set left end derivative via shape-preserving 3-point formula */
  h1 = x[1]-x[0];
  h2 = x[2]-x[1];
  hsum = h1+h2;
  del1 = (y[1]-y[0])/h1;
  del2 = (y[2]-y[1])/h2;
  w1 = (h1+hsum)/hsum;
  w2 = -h1/hsum;
  yd[0][1] = w1*del1+w2*del2;
  if (yd[0][1]*del1<=0.0)
    yd[0][1] = 0.0;
  else if (del1*del2<0.0) {
    dmax = 3.0*del1;
    if (ABS(yd[0][1])>ABS(dmax)) yd[0][1] = dmax;
  }
  
  /* loop over interior points */
  for (i=1; i<n-1; i++) {
    
    /* compute intervals and slopes */
    h1 = x[i]-x[i-1];
    h2 = x[i+1]-x[i];
    hsum = h1+h2;
    del1 = (y[i]-y[i-1])/h1;
    del2 = (y[i+1]-y[i])/h2;
    
    /* if not strictly monotonic, zero derivative */
    if (del1*del2<=0.0) {
      yd[i][1] = 0.0;
      
      /*
       * else, if strictly monotonic, use Butland's formula:
       *      3*(h1+h2)*del1*del2
       * -------------------------------
       * ((2*h1+h2)*del1+(h1+2*h2)*del2)
       * computed as follows to avoid roundoff error
       */
    } else {
      hsum3 = hsum+hsum+hsum;
      w1 = (hsum+h1)/hsum3;
      w2 = (hsum+h2)/hsum3;
      dmin = MIN(ABS(del1),ABS(del2));
      dmax = MAX(ABS(del1),ABS(del2));
      drat1 = del1/dmax;
      drat2 = del2/dmax;
      yd[i][1] = dmin/(w1*drat1+w2*drat2);
    }
  }
  
  /* set right end derivative via shape-preserving 3-point formula */
  w1 = -h2/hsum;
  w2 = (h2+hsum)/hsum;
  yd[n-1][1] = w1*del1+w2*del2;
  if (yd[n-1][1]*del2<=0.0)
    yd[n-1][1] = 0.0;
  else if (del1*del2<0.0) {
    dmax = 3.0*del2;
    if (ABS(yd[n-1][1])>ABS(dmax)) yd[n-1][1] = dmax;
  }
  
  /* compute 2nd and 3rd derivatives of cubic polynomials */
  for (i=0; i<n-1; i++) {
    h2 = x[i+1]-x[i];
    del2 = (y[i+1]-y[i])/h2;
    divdf3 = yd[i][1]+yd[i+1][1]-2.0*del2;
    yd[i][2] = 2.0*(del2-yd[i][1]-divdf3)/h2;
    yd[i][3] = (divdf3/h2)*(6.0/h2);
  }
  yd[n-1][2] = yd[n-2][2]+(x[n-1]-x[n-2])*yd[n-2][3];
  yd[n-1][3] = yd[n-2][3];
}

// Interpolate the spline fit on to a given set of x-values
void TKspline_interpolate(int n,float *x,float **yd,float *interpX,
			  float *interpY,int nInterp)
{
  float h;
  int jpos;
  int i,j;

  //To evaluate y(x) for x between x[i] and x[i+1] and h = x-x[i],
  //use the computed coefficients as follows:
  //y(x) = yd[i][0]+h*(yd[i][1]+h*(yd[i][2]/2.0+h*yd[i][3]/6.0))
  
  // Assume the data are sorted so x[i] < x[i+1]
  for (i=0;i<nInterp;i++)
    {
      jpos=-1;
      for (j=0;j<n-1;j++)
	{
	  if (interpX[i]>=x[j] && interpX[i]<x[j+1])
	    {
	      jpos=j;
	      break;
	    }
	}
      if (jpos!=-1)
	{
	  h = interpX[i]-x[jpos];
	  interpY[i] = yd[jpos][0]+h*(yd[jpos][1]+h*(yd[jpos][2]/2.0+h*yd[jpos][3]/6.0));
	}
      else
	interpY[i] = 0.0;
    }

}


//
// Smoothing spline routines
// These are obtained from: https://www.ltrr.arizona.edu/pub/trees/doc/old_html/spline_8c-source.html
//
/*******************************************************************
  RTN SPLINE: Fits cubic smoothing spline to time series

  Derived from IMSL routines by Edward R Cook, Tree Ring Laboratory,
  Lamont-Doherty Earth Observatory, Palisades, New York, USA

  Four routines combined into one by
  Richard L Holmes, University of Arizona, Tucson, Arizona, USA
  Modified copyright (C) 10 AUG 1998

********************************************************************/


/* DYNAMICALLY ALLOCATE A 2-D ARRAY */
/* Assumption:  nrows*ncols*element_size, rounded up to a multiple   */
/* of sizeof(long double), must fit in a long type.  If not, then    */
/* the "i += ..." step could overflow.                               */

void **MATRIX(int nrows, int ncols, int first_row, int first_col,int element_size)
{
    void **p;
    int alignment;
    long i;
    
    if(nrows < 1 || ncols < 1) return(NULL);
    i = nrows*sizeof(void *);
    /* align the addr of the data to be a multiple of sizeof(long double) */
    alignment = i % sizeof(long double);
    if(alignment != 0) alignment = sizeof(long double) - alignment;
    i += nrows*ncols*element_size+alignment;
    if((p = (void **)malloc((size_t)i)) != NULL)
    {
        /* compute the address of matrix[first_row][0] */
        p[0] = (char *)(p+nrows)+alignment-first_col*element_size;
        for(i = 1; i < nrows; i++)
            /* compute the address of matrix[first_row+i][0] */
            p[i] = (char *)(p[i-1])+ncols*element_size;
        /* compute the address of matrix[0][0] */
        p -= first_row;
    }
    return(p);
}



/* This function is called from SPLINE when :  */
/* 1. Series is too short to compute spline    */
/* 2. Matrix A is not positive definite        */

void errorAction(int N, double *Y, float *ZF)
{
    int k;
    double ZN;

    ZN = 0.0;
    for(k = 1; k <= N; k++)
        ZN = ZN + Y[k];

    if (N > 0)
    {
        ZN = ZN/(float)N;
        for(k = 1; k <= N; k++)
            ZF[k-1] = ZN;
    }

    return;
}



/* Function  SPLINE: Fits cubic smoothing spline to time series               */
/* Arguments:                                                                 */
/*                                                                            */
/* N:   Number of values in time series                                       */
/* Z:   Time series array to be modeled with spline                           */
/* ZF:  Computed cubic spline function array fit to time series               */
/* ZSP: Length (rigidity) of spline to be used to model series                */
/* ZPV: Portion of variance at wavelength ZSP contained in spline (0<ZPV<1)   */
/*                                                                            */
/* Arguments Z, ZF, ZSP and ZPV are single precision;                         */
/* computation is done entirely in double-precision arithmetic                */

void SPLINE(int N, float *Z, float *ZF, float ZSP, float ZPV)
{
    int i, j, k, l, m;
    int NC, NC1, NCP1, IMNCP1, I1, I2, JM1, IW, KL, N1, K1;
    double RSP, VPV, PD, RN, D1, D2, SUM;
    double **A, *F, *Y, C1[5], C2[4];

    C1[0] =  0.0;
    C1[1] =  1.0;
    C1[2] = -4.0;
    C1[3] =  6.0;
    C1[4] = -2.0;
    
    C2[0] =  0.0;
    C2[1] =  0.0;
    C2[2] =  0.33333333333333;
    C2[3] =  1.33333333333333;
    
    /* Allocate arrays to store intermeediate results */
    A = (double **)MATRIX(N+1, 5, 0, 0, sizeof(double));
    F = (double *)malloc((N+1)*sizeof(double));
    Y = (double *)malloc((N+1)*sizeof(double));
    if (A == NULL || F == NULL || Y == NULL)
    {
        printf("\nSPLINE >> Unable to allocate memory\n");
        return;
    }
    
    /* Check whether series is too short to compute spline */
    if (N < 4)
    {
        errorAction(N, Y, ZF);
        return;
    }

    /* Copy time series into double precision array */
    for(j = 1; j <= N; j++)
        Y[j] = (double)Z[j-1];
        
    /* Compute Lagrange multiplier, which defines frequency response of spline */
    RSP = (double)ZSP;
    VPV = (double)ZPV;
    PD = ((1.0/(1.0-VPV)-1.0)*6.0* pow((cos(M_PI*2.0/RSP)-1.0),2.0))/(cos(M_PI*2.0/RSP)+2.0);
    for(i = 1; i <= N-2; i++)
        for(j = 1; j <= 3; j++)
        {
            A[i][j] = C1[j] + PD * C2[j];
            A[i][4] = Y[i] + C1[4] * Y[i+1] + Y[i+2];
        }

    A[1][1] = C2[1];
    A[1][2] = C2[1];
    A[2][1] = C2[1];
    NC = 2;

    /* Begin LUDAPB */
    RN = 1.0/((double)(N-2) * 16.0);
    D1 = 1.0;
    D2 = 0.0;
    NCP1 = NC + 1;

    /* Initialize zero elements */
    if (NC != 0)
    {
        for(i = 1; i <= NC; i++)
            for(j = i; j <= NC; j++)
            {
                k = NCP1 - j;
                A[i][k] = 0.0;
            }
    }

    /* i: row index of element being computed */
    /* j: column index of element being computed */
    /* l: row index of previously computed vector being used to compute inner product */
    /* m: column index */
    for(i = 1; i <= N-2; i++)
    {
        IMNCP1 = i - NCP1;
        I1 = (1 < 1 - IMNCP1? 1 - IMNCP1: 1);
        for(j = I1; j <= NCP1; j++)
        {
            l = IMNCP1 + j;
            I2 = NCP1 - j;
            SUM = A[i][j];
            JM1 = j - 1;

            if (JM1 > 0)
            {
                for(k = 1; k <= JM1; k++)
                {
                    m = I2 + k;
                    SUM = SUM - (A[i][k]*A[l][m]);
                }
            }

            /* Matrix not positive definite */
            if (j == NCP1)
            {
                if (A[i][j]+SUM*RN <= A[i][j])
                {
                    printf("\nSPLINE >> Matrix not positive definite\n");
                    errorAction(N, Y, ZF);
                    return;             
                }
                
                A[i][j] = 1.0/sqrt(SUM);
                
                /* Update determinant */
                D1 = D1*SUM;
                while (fabs(D1) > 1.0)
                {
                    D1 = D1*0.0625;
                    D2 = D2+4.0;
                }
                
                while (fabs(D1) <= 0.0625)
                {
                    D1 = D1*16.0;
                    D2 = D2-4.0;
                }
                continue;
            }
            A[i][j] = SUM*A[l][NCP1];                           
        }
    }
    /* End LUDAPB */

    /* Begin LUELPB */
    /* Solution LY = B */
    NC1 = NC + 1;
    IW = 0;
    l = 0;

    for(i = 1; i <= N-2; i++)
    {
        SUM = A[i][4];
        if (NC > 0) {
            if (IW != 0) {
                l = l + 1;
                if (l > NC)
                {
                    l = NC;
                }
                k = NC1 - l;
                KL = i - l;

                for(j = k; j <= NC; j++)
                {
                    SUM = SUM - (A[KL][4]*A[i][j]);
                    KL = KL + 1;
                }
            } else if (SUM != 0.0)
            {
                IW = 1;
            }
        }
        A[i][4] = SUM*A[i][NC1];
    }

    /* Solution UX = Y */
    A[N-2][4] = A[N-2][4]*A[N-2][NC1];
    N1 = N-2+1;

    for(i = 2; i <= N-2; i++)
    {
        k = N1 - i;
        SUM = A[k][4];
        if (NC > 0)
        {
            KL = k +1;
            K1 = (N-2 < k+NC? N-2: k+NC);
            l = 1;
            for(j = KL; j <= K1; j++)
            {
                SUM = SUM - A[j][4]*A[j][NC1-l];
                l = l + 1;
            }
        }
        A[k][4] = SUM*A[k][NC1];
    }
    /* End LUELPB */

    /* Calculate Spline Curve */
    for(i = 3; i <= N-2; i++)
        F[i] = A[i-2][4]+C1[4]*A[i-1][4]+A[i][4];

    F[1] = A[1][4];
    F[2] = C1[4]*A[1][4]+A[2][4];
    F[N-1] = A[N-2-1][4]+C1[4]*A[N-2][4];
    F[N] = A[N-2][4];

    for(i = 1; i <= N; i++)
        ZF[i-1] = Y[i] - F[i];

    return;
}

void writeModel(sdhdf_fileStruct *outFile,dataModelStruct *model,int nband) 
{
  int nbeam=1;
  int i,j,k,b,ii;
  sdhdf_bandHeaderStruct *inBandParams;
  int nsub;
  int npol;
  int nchan;
  
  // Now replace the relevant spectra
  b=0; // FIX ME
  for (i=0;i<nband;i++)
    {
      nchan = model[i].original_nchan;
      nsub  = model[i].original_nsub;
      npol  = model[i].original_npol;

      sdhdf_replaceSpectrumData(outFile,model[i].bandLabel,b,i,model[i].model,nsub,npol,nchan);
    }
}
