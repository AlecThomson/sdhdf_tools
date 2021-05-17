#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sdhdfProc.h"
#include "hdf5.h"
#include "fitsio.h"


// Compilation
// gcc -lm -o sdhdf_convertTo sdhdf_convertTo.c -lcfitsio -I../hdf5-1.10.4/src/ ../hdf5-1.10.4/src/.libs/libhdf5.a -ldl -lz -L/pulsar/psr/software/stable/stretch/lib/ -I//pulsar/psr/software/stable/stretch/include -lcalceph 


int main(int argc,char *argv[])
{
  int i,j,ii,k,kk;
  char fname[MAX_STRLEN];
  char oname[MAX_STRLEN];
  char ext[MAX_STRLEN] = "fits";

  // SDHDF file
  sdhdf_fileStruct *inFile;
  
  // Parameters for fits file
  fitsfile *fptr;
  int fitsStatus=0;
  float fval;
  float *fdata;
  float hr,min,sec;
  int   colnum;
  int   rowNum=0,ival;
  int   scan=0;
  long naxes[4];
  float intTime;  // Total integration time
  char **strArray;
  char timeStr[1024];

  int ibeam = 0;
  
  if (!(inFile = (sdhdf_fileStruct *)malloc(sizeof(sdhdf_fileStruct))))
    {
      printf("ERROR: unable to allocate sufficient memory for >inFile<\n");
      exit(1);
    }
      
  strArray = (char **)malloc(sizeof(char *));
  strArray[0] = (char *)malloc(sizeof(char)*1024);

  printf("WARNING: Various parameters not being set:\n"); //Should set FOCUSTAN correctly\n");
  printf("FOCUSROT, FOCUSTAN, SCANRATE, TSYS, CALFCTR, TCAL, TCALTIME, AZIMUTH, ELEVATIO, TAMBIENT, PRESSURE, HUMIDITY, WINDSPEE, WINDDIRE\n");

  
  
  for (i=1;i<argc;i++)
    {
      strcpy(fname,argv[i]);

      sdhdf_initialiseFile(inFile);
      sdhdf_openFile(fname,inFile,1);
      sdhdf_loadMetaData(inFile);
      
      // NOTE HARDCODE HERE
      sprintf(oname,"!%s.%s(/u/hob044/software/new_c/sdhdfProc/sdHeader.fits)",fname,ext);
      printf("opening  %s %s\n",fname,oname);

      fits_create_file(&fptr,oname,&fitsStatus);
      fits_report_error(stderr, fitsStatus);  /* print out any error messages */

      // Primary header information
      fits_movabs_hdu(fptr,1,NULL,&fitsStatus);
      fits_write_date(fptr,&fitsStatus);

      // Header information in SINGLE_DISH table
      fits_movnam_hdu(fptr, BINARY_TBL, (char *)"SINGLE DISH", 0, &fitsStatus);
      fits_update_key(fptr,TSTRING,(char *)"TELESCOP",(char *)"PARKES",NULL,&fitsStatus);
      //       fits_update_key(fptr,TSTRING,(char *)"TELESCOP",(char *)"ATPKSMB",NULL,&fitsStatus);
      fits_update_key(fptr,TSTRING,(char *)"INSTRUME",(char *)"Medusa",NULL,&fitsStatus);
      fits_update_key(fptr,TSTRING,(char *)"DATE-OBS",inFile->primary[0].utc0,NULL,&fitsStatus);
      fits_update_key(fptr,TSTRING,(char *)"PROJID",inFile->primary[0].pid,NULL,&fitsStatus);
      fits_update_key(fptr,TSTRING,(char *)"FRONTEND",(char *)"UWL",NULL,&fitsStatus);
      fits_update_key(fptr,TSTRING,(char *)"HDRVER",inFile->primary[0].hdr_defn_version,NULL,&fitsStatus);
      fits_update_key(fptr,TSTRING,(char *)"OBSERVER",inFile->primary[0].observer,NULL,&fitsStatus);
      fits_update_key(fptr,TSTRING,(char *)"EQUINOX",(char *)"2000",NULL,&fitsStatus);
      fits_update_key(fptr,TSTRING,(char *)"RADESYS",(char *)"FK5",NULL,&fitsStatus);
      fits_update_key(fptr,TSTRING,(char *)"SYSOBS",(char *)"TOPOCENT",NULL,&fitsStatus);
      
      //       fits_update_key(fp,TSTRING,(char *)"SPECSYS",(char *)"TOPO",NULL,&fitsStatus);
      printf("Setting SPECSYS to LSRK -- I think this should be TOPO, but that doesn't work in ASAP\n");
      fits_update_key(fptr,TSTRING,(char *)"SPECSYS",(char *)"LSRK",NULL,&fitsStatus);
      //
       /*
	 CTYPE2          'STOKES  '               DATA array axis 2: polarization code
       CRPIX2          1.0                      Polarization code reference pixel
	 CRVAL2          -5.0                     Polarization code at reference pixel (XX)
	 CDELT2          -1.0                     Polarization code axis increment
	 CTYPE3          'RA      '               DATA array axis 3 (degenerate): RA (mid-int)
	 CRPIX3          1.0                      RA reference pixel
	 TTYPE18         'CRVAL3  '               label for field
	 TFORM18         '1D      '               format of field
	   TUNIT18         'deg     '               units of field
	   CDELT3          -1.0                     RA axis increment
	   CTYPE4          'DEC     '               DATA array axis 4 (degenerate): Dec (mid-int)
	   CRPIX4          1.0                      Dec reference pixel
	   TTYPE19         'CRVAL4  '               label for field
	   TFORM19         '1D      '               format of field
	     TUNIT19         'deg     '               units of field
	     CDELT4          1.0                      Dec axis increment
       */    


       
       fits_update_key(fptr,TSTRING,(char *)"CTYPE2",(char *)"STOKES",NULL,&fitsStatus); // SHOULDn'T HARD CODE 17 HERE
       fval = 1.0; fits_update_key(fptr,TFLOAT,(char *)"CRPIX2",&fval,NULL,&fitsStatus);
       fval = -5.0; fits_update_key(fptr,TFLOAT,(char *)"CRVAL2",&fval,NULL,&fitsStatus); // See https://casa.nrao.edu/aips2_docs/notes/236/node14.html
       fval = -1.0; fits_update_key(fptr,TFLOAT,(char *)"CDELT2",&fval,NULL,&fitsStatus);
       printf("Incorrectly setting CDELT3 and CDELT4\n");
       fval = -1.0; fits_update_key(fptr,TFLOAT,(char *)"CDELT3",&fval,NULL,&fitsStatus);
       fval = 1.0; fits_update_key(fptr,TFLOAT,(char *)"CDELT4",&fval,NULL,&fitsStatus);

       // These are taken from an old SDFITS file
       //
       fval = -4.554232087E+06; fits_update_key(fptr,TFLOAT,(char *)"OBSGEO-X",&fval,NULL,&fitsStatus);
       fval = 2.816759046E+06; fits_update_key(fptr,TFLOAT,(char *)"OBSGEO-Y",&fval,NULL,&fitsStatus);
       fval = -3.454035950E+06; fits_update_key(fptr,TFLOAT,(char *)"OBSGEO-Z",&fval,NULL,&fitsStatus);

       //
       // ******************* Process each sub-band
       for (ii=0;ii<inFile->beam[ibeam].nBand;ii++)
	 {
	   printf("Processing sub-band %d (%s)\n",ii,inFile->beam[ibeam].bandHeader[ii].label);
	   naxes[0] = inFile->beam[ibeam].bandHeader[ii].nchan;
	   naxes[1] = inFile->beam[ibeam].bandHeader[ii].npol;
	   naxes[2] = 1;
	   naxes[3] = 1;
	   if (strlen(inFile->beam[ibeam].bandHeader[ii].label)>0)
	     {
	       printf("Checking subintegrations\n");
	       // ******************* Now process each subintegration
	       fdata = (float *)malloc(sizeof(float)*inFile->beam[ibeam].bandHeader[ii].nchan*inFile->beam[ibeam].bandHeader[ii].npol);
	       sdhdf_loadBandData(inFile,ibeam,ii,1);

	       
	       // Need to use the existing arrays of copy the data into fdata **** GEORGE TO DO *****
	       //	       sdhdf_loadEntireDump(inFile,ii,fdata);

	       intTime = inFile->beam[ibeam].bandHeader[ii].dtime * inFile->beam[ibeam].bandHeader[ii].ndump; /// GEORGE CHECK WHAT THIS SHOULD BE
	       
	       for (j=0;j<inFile->beam[ibeam].bandHeader[ii].ndump;j++)
		 {
		   for (k=0;k<inFile->beam[ibeam].bandHeader[ii].npol;k++)
		     {
		       for (kk=0;kk<inFile->beam[ibeam].bandHeader[ii].nchan;kk++)
			 {
			   if (k==0)
			     {
			       fdata[k*inFile->beam[ibeam].bandHeader[ii].nchan+kk] = inFile->beam[ibeam].bandData[ii].astro_data.pol1[kk+j*inFile->beam[ibeam].bandHeader[ii].nchan];
			       //			       printf("Got %g\n",inFile->beam[ibeam].bandData[ii].astro_data.pol1[kk+j*inFile->beam[ibeam].bandHeader[ii].nchan]);
				      
			     }
			   else if (k==1)
			     fdata[k*inFile->beam[ibeam].bandHeader[ii].nchan+kk] = inFile->beam[ibeam].bandData[ii].astro_data.pol2[kk+j*inFile->beam[ibeam].bandHeader[ii].nchan];
			   else if (k==2)
			     fdata[k*inFile->beam[ibeam].bandHeader[ii].nchan+kk] = inFile->beam[ibeam].bandData[ii].astro_data.pol3[kk+j*inFile->beam[ibeam].bandHeader[ii].nchan];
			   else if (k==4)
			     fdata[k*inFile->beam[ibeam].bandHeader[ii].nchan+kk] = inFile->beam[ibeam].bandData[ii].astro_data.pol4[kk+j*inFile->beam[ibeam].bandHeader[ii].nchan];
			 }
		     }

		   //		   printf("Checking fdata\n");
		   //		   for (k=0;k<inFile->beam[ibeam].bandHeader[ii].npol*inFile->beam[ibeam].bandHeader[ii].nchan;k++)
		   //		     printf("%g\n",fdata[k]);
		   
		   //		   printf("Subintegration %d/%d (sb %d/%d)\n",j,inFile->beam[ibeam].bandHeader[ii].ndump,ii,inFile->beam[ibeam].nBand);
		   fits_insert_rows(fptr,rowNum,1,&fitsStatus);
		   fits_get_colnum(fptr,CASEINSEN,(char *)"SCAN",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   ival = scan+1;  fits_write_col(fptr,TINT,colnum,rowNum+1,1,1,&ival,&fitsStatus);
		   fits_get_colnum(fptr,CASEINSEN,(char *)"CYCLE",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   ival = j;  fits_write_col(fptr,TINT,colnum,rowNum+1,1,1,&ival,&fitsStatus);
		   fits_get_colnum(fptr,CASEINSEN,(char *)"IF",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   printf("Setting IF to %d %d %d\n",colnum,rowNum+1,ii+1);
		   ival = ii+1;  fits_write_col(fptr,TINT,colnum,rowNum+1,1,1,&ival,&fitsStatus);
		   fits_get_colnum(fptr,CASEINSEN,(char *)"BEAM",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   ival = 1;  fits_write_col(fptr,TINT,colnum,rowNum+1,1,1,&ival,&fitsStatus);
		   fits_get_colnum(fptr,CASEINSEN,(char *)"DATE-OBS",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   strcpy(strArray[0],inFile->primary[ibeam].utc0);
		   fits_write_col(fptr,TSTRING,colnum,rowNum+1,1,1,strArray,&fitsStatus);
		
		   // TIME
		   sprintf(timeStr,inFile->primary[ibeam].utc0+11); // CHECK THAT THIS IS CORRECT *** GEORGE
		   sscanf(timeStr,"%f:%f:%f",&hr,&min,&sec);
		   fval = (hr*60*60+min*60.+sec);  
		   fits_get_colnum(fptr,CASEINSEN,(char *)"TIME",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);
		   fits_get_colnum(fptr,CASEINSEN,(char *)"EXPOSURE",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval = intTime; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);
		   fits_get_colnum(fptr,CASEINSEN,(char *)"OBJECT",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   strcpy(strArray[0],inFile->beamHeader[ibeam].source);
		   fits_write_col(fptr,TSTRING,colnum,rowNum+1,1,1,strArray,&fitsStatus);
		   
		   // OBJ-RA
		   fits_get_colnum(fptr,CASEINSEN,(char *)"OBJ-RA",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval = inFile->beam[ibeam].bandData[ii].astro_obsHeader[j].raDeg; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);		   
		
		   // OBJ-DEC
		   fits_get_colnum(fptr,CASEINSEN,(char *)"OBJ-DEC",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval = inFile->beam[ibeam].bandData[ii].astro_obsHeader[j].decDeg; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);

		   // RESTFRQ
		   fits_get_colnum(fptr,CASEINSEN,(char *)"RESTFRQ",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval = 704e6; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);
		   
		   
		   // OBSMODE
		   strcpy(strArray[0],"POINTING");
		   fits_get_colnum(fptr,CASEINSEN,(char *)"OBSMODE",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fits_write_col(fptr,TSTRING,colnum,rowNum+1,1,1,strArray,&fitsStatus);



		   // FREQRES
		   fits_get_colnum(fptr,CASEINSEN,(char *)"FREQRES",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval = 1e6*fabs(inFile->beam[ibeam].bandHeader[ii].f1-inFile->beam[ibeam].bandHeader[ii].f0)/inFile->beam[ibeam].bandHeader[ii].nchan; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);
		   
		   // BANDWID
		   fits_get_colnum(fptr,CASEINSEN,(char *)"BANDWID",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval = 1e6*fabs(inFile->beam[ibeam].bandHeader[ii].f1-inFile->beam[ibeam].bandHeader[ii].f0); fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);

		   // CRPIX1 - must check if this is NCHAN - "array location of the reference pixel along axis 1"
		   fits_get_colnum(fptr,CASEINSEN,(char *)"CRPIX1",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   //		   fval = (float)inFile->beam[ibeam].bandHeader[ii].nchan; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);
		   fval = 0; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);
		   
		   // CRVAL1
		   fits_get_colnum(fptr,CASEINSEN,(char *)"CRVAL1",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval = 1e6*inFile->beam[ibeam].bandHeader[ii].f0; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);

		   // CDELT1
		   fits_get_colnum(fptr,CASEINSEN,(char *)"CDELT1",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval = 1e6*fabs(inFile->beam[ibeam].bandHeader[ii].f1-inFile->beam[ibeam].bandHeader[ii].f0)/inFile->beam[ibeam].bandHeader[ii].nchan; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);

		   // CRVAL3
		   fits_get_colnum(fptr,CASEINSEN,(char *)"CRVAL3",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval =  inFile->beam[ibeam].bandData[ii].astro_obsHeader[j].raDeg; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);
		   
		   // CRVAL4
		   fits_get_colnum(fptr,CASEINSEN,(char *)"CRVAL4",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval =  inFile->beam[ibeam].bandData[ii].astro_obsHeader[j].decDeg; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);



		   // SCANRATE
		   // TSYS
		   // CALFCTR
		   // TCAL
		   // TCALTIME
		   // AZIMUTH
		   // ELEVATIO
		   // PARANGLE
		   
		   fits_get_colnum(fptr,CASEINSEN,(char *)"PARANGLE",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval =  inFile->beam[ibeam].bandData[ii].astro_obsHeader[j].paraAngle; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);
		   
		   // FOCUSTAN
		   fits_get_colnum(fptr,CASEINSEN,(char *)"FOCUSTAN",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval = 0; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);
		   
		   // FOCUSROT
		   fits_get_colnum(fptr,CASEINSEN,(char *)"FOCUSROT",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fval = 0; fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,1,&fval,&fitsStatus);
		   
		   // TAMBIENT
		   // PRESUURE
		   // HUMIDITY
		   // WINDSPEE
		   // WINDDIRE
		   
		   // Data
		   		   
		   fits_get_colnum(fptr,CASEINSEN,(char *)"DATA",&colnum,&fitsStatus); fits_report_error(stderr,fitsStatus);
		   fits_modify_vector_len (fptr, colnum, inFile->beam[ibeam].bandHeader[ii].nchan*inFile->beam[ibeam].bandHeader[ii].npol, &fitsStatus); fits_report_error(stderr,fitsStatus);
		   fits_delete_key(fptr, (char *)"TDIM23", &fitsStatus); // THIS SHOULD NOT BE HARDCODED
		   if (fitsStatus) { fits_report_error(stderr,fitsStatus); exit(1);}
		   fits_write_tdim(fptr,colnum,4,naxes,&fitsStatus);fits_report_error(stderr,fitsStatus);
		   if (fitsStatus) { fits_report_error(stderr,fitsStatus); exit(1);}
		   fits_write_col(fptr,TFLOAT,colnum,rowNum+1,1,inFile->beam[ibeam].bandHeader[ii].nchan*inFile->beam[ibeam].bandHeader[ii].npol,fdata,&fitsStatus);
		   rowNum++;
		 }
	       free(fdata);
	     }
 
	 }

       
       

      fits_report_error(stderr, fitsStatus);  /* print out any error messages */
      fits_close_file(fptr,&fitsStatus);
      printf("Complete\n");


      sdhdf_closeFile(inFile);
    }
  free(inFile);
  free(strArray[0]);
  free(strArray);

}
