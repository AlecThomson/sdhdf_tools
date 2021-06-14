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
#include <string.h>
#include <math.h>
#include <hdf5.h>
#include "sdhdf_v1.9.h"

#define SOFTWARE_VER "v0.1"
#define MAX_STRLEN    512
#define MAX_FILES     8192         // Maximum number of files to be processed in batch processing
#define MAX_CHAN_SCAL 65536        // Maximum number of channels in the SCAL measurements

#define SPEED_LIGHT   299792458.0  // Speed of light
#define BOLTZMANN 1.38064852e-23   // Boltzmann constant
#define MAX_EOP_LINES 208480       // Maximum number of lines in the Earth Orientation parameter file
#define MAX_POLY_COEFF 16          // Maximum number of coefficients in baseline fits

// Handling of geodetic coordinates
#define GRS80_A 6378137.0          // semi-major axis (m) 
#define GRS80_F 1.0/298.257222101  // flattening 


// Calibrator temperature structure
typedef struct sdhdf_tcal_struct {
  double freq;
  double tcalA;
  double tcalB;
} sdhdf_tcal_struct;


// Earth orientation parameters
typedef struct sdhdf_eopStruct {
  double mjd;
  double x;
  double y;
  double dut1;
} sdhdf_eopStruct;


// Parameters relating to baselines
typedef struct baselineFitStruct {
  float  applicableRange_f0; // Minimum frequency (MHz) for where this baseline fitting is applicable
  float  applicableRange_f1; // Maximum frequency (MHz) for where this baseline fitting is applicable
  double fitX0;
  double fitCoeff_pol1[MAX_POLY_COEFF];
  double fitCoeff_pol2[MAX_POLY_COEFF];
  double fitCoeff_pol3[MAX_POLY_COEFF];
  double fitCoeff_pol4[MAX_POLY_COEFF];
  int    nCoeff;
  int    set;                // 0 = not set, 1 = set
} baselineFitStruct;


// Function definitions

// File input/output
void sdhdf_initialiseFile(sdhdf_fileStruct *inFile);
int  sdhdf_openFile(char *fname,sdhdf_fileStruct *inFile,int openType);
void sdhdf_closeFile(sdhdf_fileStruct *inFile);

// Data 
void sdhdf_releaseBandData(sdhdf_fileStruct *inFile,int beam,int band,int type);
void sdhdf_loadBandData(sdhdf_fileStruct *inFile,int beam,int band,int type);
void sdhdf_loadBandData2Array(sdhdf_fileStruct *inFile,int beam,int band,int type,float *arr);
void sdhdf_loadFrequency2Array(sdhdf_fileStruct *inFile,int beam,int band,float *arr);
void sdhdf_allocateBandData(sdhdf_spectralDumpsStruct *spec,int nchan,int ndump,int npol);
void sdhdf_extractPols(sdhdf_spectralDumpsStruct *spec,float *in,int nchan,int ndump,int npol);

// Command line arguments
void sdhdf_add1arg(char *args,char *add);
void sdhdf_add2arg(char *args,char *add1,char *add2);


// Loading metadata information
void sdhdf_loadMetaData(sdhdf_fileStruct *inFile);  // Include loading attributes
void sdhdf_loadPrimaryHeader(sdhdf_fileStruct *inFile);
void sdhdf_loadBeamHeader(sdhdf_fileStruct *inFile);
void sdhdf_loadBandHeader(sdhdf_fileStruct *inFile,int type);
void sdhdf_loadObsHeader(sdhdf_fileStruct *inFile,int type);
void sdhdf_loadHistory(sdhdf_fileStruct *inFile);
void sdhdf_loadSoftware(sdhdf_fileStruct *inFile);
void sdhdf_loadSpectrum(sdhdf_fileStruct *inFile,int ibeam, int iband, sdhdf_spectralDumpsStruct *spectrum);
void sdhdf_initialise_spectralDumps(sdhdf_spectralDumpsStruct *in);
void sdhdf_initialise_bandHeader(sdhdf_bandHeaderStruct *header);
void sdhdf_initialise_obsHeader(sdhdf_obsParamsStruct *obs);
void sdhdf_copyBandHeaderStruct(sdhdf_bandHeaderStruct *in,sdhdf_bandHeaderStruct *out,int n);
void sdhdf_allocateBeamMemory(sdhdf_fileStruct *inFile,int nbeam);
int sdhdf_checkGroupExists(sdhdf_fileStruct *inFile,char *groupName);
void sdhdf_copySingleObsParams(sdhdf_fileStruct *inFile,int ibeam,int iband,int idump,sdhdf_obsParamsStruct *obsParam);
void sdhdf_loadCalProc(sdhdf_fileStruct *inFile,int ibeam,int iband,char *cal_label,float *vals);

// Attributes
int sdhdf_getNattributes(sdhdf_fileStruct *inFile,char *dataName);
void sdhdf_readAttributeFromNum(sdhdf_fileStruct *inFile,char *dataName,int num,sdhdf_attributes_struct *attribute);
void sdhdf_copyAttributes(sdhdf_attributes_struct *in,int n_in,sdhdf_attributes_struct *out,int *n_out);
void sdhdf_writeAttribute(sdhdf_fileStruct *outFile,char *dataName,char *attrName,char *result);

// Loading data
void sdhdf_loadFrequencies(sdhdf_fileStruct *inFile,int ibeam,int iband,int type);
void sdhdf_loadData(sdhdf_fileStruct *inFile,int ibeam,int iband,float *in_data,int type);
int  sdhdf_loadFlagData(sdhdf_fileStruct *inFile);
int  sdhdf_loadTcal(sdhdf_tcal_struct *tcalData,char *fname);
int  sdhdf_loadTsys(sdhdf_fileStruct *inFile,int band,float *tsys);
int  sdhdf_loadPhase(sdhdf_fileStruct *inFile,int band,float *phase);
int  sdhdf_loadGain(sdhdf_fileStruct *inFile,int band,float *gain);

// SDHDF identification
int sdhdf_doWeHave(sdhdf_fileStruct *inFile,char *groupLabel);
int sdhdf_getBandID(sdhdf_fileStruct *inFile,char *input);

// Manipulation and processing
void sdhdf_convertStokes(float p1,float p2,float p3,float p4,float *stokesI,float *stokesQ,float *stokesU,float *stokesV);
void sdhdf_get_tcal(sdhdf_tcal_struct *tcalData,int n,double f0,double *tcalA,double *tcalB);
void sdhdf_addHistory(sdhdf_historyStruct *history,int n,char *procName,char *descr,char *args);

// Output information
void sdhdf_copyEntireGroup(char *groupLabel,sdhdf_fileStruct *inFile,sdhdf_fileStruct *outFile); 
void sdhdf_copyEntireGroupDifferentLabels(char *bandLabelIn,sdhdf_fileStruct *inFile,char *bandLabelOut,sdhdf_fileStruct *outFile);
void sdhdf_setMetadataDefaults(sdhdf_primaryHeaderStruct *primaryHeader,sdhdf_beamHeaderStruct *beamHeader,
			       sdhdf_bandHeaderStruct *bandHeader,sdhdf_softwareVersionsStruct *softwareVersions,sdhdf_historyStruct *history,
			       int nbeam,int nband);
void sdhdf_writeHistory(sdhdf_fileStruct *outFile,sdhdf_historyStruct *outParams,int n);
void sdhdf_writeSoftwareVersions(sdhdf_fileStruct *outFile,sdhdf_softwareVersionsStruct *outParams);
void sdhdf_writePrimaryHeader(sdhdf_fileStruct *outFile,sdhdf_primaryHeaderStruct *primaryHeader);
void sdhdf_writeBandHeader(sdhdf_fileStruct *outFile,sdhdf_bandHeaderStruct *outBandParams,int ibeam,int outBands,int type);
void sdhdf_writeBeamHeader(sdhdf_fileStruct *outFile,sdhdf_beamHeaderStruct *beamHeader,int nBeams);
void sdhdf_writeSpectrumData(sdhdf_fileStruct *outFile,char *blabel, int ibeam,int iband,  float *out,float *freq,long nchan,long npol,long nsub,int type,sdhdf_attributes_struct *dataAttributes,int nDataAttributes,sdhdf_attributes_struct *freqAttributes,int nFreqAttributes);
void sdhdf_copyRemainder(sdhdf_fileStruct *inFile,sdhdf_fileStruct *outFile,int type);
void sdhdf_writeObsParams(sdhdf_fileStruct *outFile,char *bandLabel,int ibeam,int iband,sdhdf_obsParamsStruct *obsParams,int ndump,int type);
void sdhdf_writeFlags(sdhdf_fileStruct *outFile,int ibeam,int iband,int *flag,int nchan,char *bandLabel);
void sdhdf_writeDataWeights(sdhdf_fileStruct *outFile,int ibeam,int iband,float *flag,int nchan,int ndump,char *bandLabel);
void sdhdf_writeCalProc(sdhdf_fileStruct *outFile,int ibeam,int iband,char *band_label,char *cal_label,float *vals,int nchan,int npol,int ndumps);

// HDF5 reading functions
void sdhdf_loadHeaderString(hid_t header_id,char *parameter,char *outStr);
void sdhdf_loadHeaderInt(hid_t header_id,char *parameter,int *outInt);
void sdhdf_loadColumnDouble(char *param,hid_t header_id,hid_t headerT,double *dVals);
void sdhdf_loadColumnString(char *param,hid_t header_id,hid_t headerT,char **vals,int nrows);
void sdhdf_loadColumnInt(char *param,hid_t header_id,hid_t headerT,long int *dVals);
void sdhdf_loadHeaderDouble(hid_t header_id,char *parameter,double *outDouble);

// Ephemerides
long   sdhdf_loadEOP(sdhdf_eopStruct *eop);
double sdhdf_calcVoverC(sdhdf_fileStruct *inFile,int ibeam,int iband,int idump,sdhdf_eopStruct *eop,int nEOP,int lsr);
void sdhdf_obsCoord_IAU2000B(double observatory_trs[3],
			     double zenith_trs[3],
			     long double tt_mjd, long double utc_mjd,
			     double observatory_crs[3],
			     double zenith_crs[3],
			     double observatory_velocity_crs[3],sdhdf_eopStruct *eop,int nEOP);
void sdhdf_interpolate_EOP(double mjd, double *xp, double *yp, double *dut1, double *dut1dot, sdhdf_eopStruct *eop,int nEOP);
void sdhdf_ITRF_to_GRS80(double x,double y,double z,double *long_grs80,double *lat_grs80,double *height_grs80);

// Mathematics
//int sdhdf_inv4x4(float m[4][4],float inv[4][4]);

void sdhdf_setIdentity_4x4(float mat[4][4]);
void sdhdf_display_vec4(float *vec);
void sdhdf_setFeed(float mf[4][4],float gamma);
void sdhdf_multMat_vec_replace(float mat[4][4],float *vec);
void sdhdf_copy_mat4(float in[4][4],float out[4][4]);
void sdhdf_mult4x4_replace(float src1[4][4], float src2[4][4]);
void sdhdf_copy_vec4(float *in,float *out);
void sdhdf_setGainPhase(float ma[4][4],float diffGain,float diffPhase);
void sdhdf_setGain2Phase(float ma[4][4],float diffGain,float diffPhase);
void sdhdf_setParallacticAngle(float msky[4][4],float pa);

int sdhdf_inv4x4(float m[4][4],float inv[4][4]);
void sdhdf_mult4x4(float src1[4][4], float src2[4][4], float dest[4][4]);
double sdhdf_dotproduct(double *v1,double *v2);
void sdhdf_para(double dxd,double ddc,double q,double *axd,double *eld);
void displayMatrix_4x4(float matrix[4][4]);

// Function definitions (OLD)

//void sdhdf_loadPrimaryMetadata(sdhdf_fileStruct *inFile);
//void sdhdf_loadBandMetadata(sdhdf_fileStruct *inFile,int type);
//void sdhdf_loadSpectralDumpMetadata(sdhdf_fileStruct *inFile);

//void sdhdf_initialiseSpectrum(spectralDumpStruct *spectrum);
//void sdhdf_load_calBandFrequencies(sdhdf_fileStruct *inFile,int iband);


/*
void sdhdf_writeBandHeader(sdhdf_fileStruct *outFile,sdhdf_bandHeaderStruct *outBandParams,int outBands,int type);
void sdhdf_writeSpectralDumpHeader(sdhdf_fileStruct *outFile,sdhdf_obsParamsStruct *outParams,int outDumps);
//void sdhdf_loadEntireDump(sdhdf_fileStruct *inFile,int band,float *in_data);
//void sdhdf_load_calEntireDump(sdhdf_fileStruct *inFile,int band,float *in_data,int onOff);
void sdhdf_writeEntireDump(sdhdf_fileStruct *outFile,sdhdf_fileStruct *inFile, int iband,  float *out,float *freq,long nchan,long npol,long nsub,int type);
void sdhdf_writeNewBand(sdhdf_fileStruct *outFile, int iband, float *vals,float *freq,sdhdf_bandHeaderStruct *bandHeader);
long sdhdf_getNdumpsFromData(sdhdf_fileStruct *inFile,int iband);
void sdhdf_freeSpectrumMemory(spectralDumpStruct *spectrum);
//void sdhdf_loadHistory(sdhdf_fileStruct *inFile,int extraMemory);
void sdhdf_writeHistory(sdhdf_fileStruct *outFile,sdhdf_historyStruct *outParams,int n);
//void sdhdf_loadSoftware(sdhdf_fileStruct *file);
void sdhdf_writeSoftware(sdhdf_fileStruct *outFile,sdhdf_softwareVersionsStruct *outParams,int n);
//void sdhdf_loadPrimary(sdhdf_fileStruct *file);
void sdhdf_writePrimary(sdhdf_fileStruct *outFile,sdhdf_primaryHeaderStruct *outParams,int n);



void sdhdf_writeDataSet(sdhdf_fileStruct *outFile,int band,float *freq,float *data,int nchan,int npol,int ndump,int type);
int sdhdf_getNattributes(sdhdf_fileStruct *inFile,char *dataName);
void sdhdf_readAttribute(sdhdf_fileStruct *inFile,char *dataName,char *attrName,char *out);
void sdhdf_loadFrequencyAttributes(sdhdf_fileStruct *inFile,int iband);
void sdhdf_loadDataAttributes(sdhdf_fileStruct *inFile,int iband);
void sdhdf_writeFrequencyAttributes(sdhdf_fileStruct *outFile,char *bandName);
void sdhdf_writeDataAttributes(sdhdf_fileStruct *outFile,char *bandName);
void sdhdf_writeCalProc(sdhdf_fileStruct *outFile,int iband,char *band_label,char *cal_label,float *vals,int nchan,int npol,int ndumps);
*/
