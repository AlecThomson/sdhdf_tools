//  Copyright (C) 2019, 2020, 2021, 2022, 2023, 2024 George Hobbs

/*
 *    This file is part of INSPECTA.
 *
 *    INSPECTA is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *    INSPECTA is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *    You should have received a copy of the GNU General Public License
 *    along with INSPECTA.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "inspecta.h"
#include <locale.h>
#include <wchar.h>
#include <time.h>
#include <sys/utsname.h>

void sdhdf_loadConfig(sdhdf_fileStruct *inFile);

void sdhdf_storeArguments(char *args,int maxLen,int argc,char *argv[])
{
  int i;
  int len1,len2;
  for (i=1;i<argc;i++)
    {
      len1 = strlen(args);
      len2 = strlen(argv[i]);
      if (len1 + len2 > maxLen-5)
	{strcat(args,"..."); break;}
      else
	{
	  strcat(args,argv[i]);
	  strcat(args," ");
	}
    }
}


// Sort out the filename when adding on an extension
// If ending of filename = hdf then remove that, add the extension and then add .hdf at the end
// If ending of filename = hdf5 then remove that, add the extension and then add .hdf5 at the end
// Otherwise just append to the end
void sdhdf_formOutputFilename(char *inFile,char *extension,char *oname)
{
  if (strcmp(inFile+strlen(inFile)-4,".hdf")==0)
    {
      char tmp[1024];
      strcpy(tmp,inFile);
      tmp[strlen(tmp)-4]='\0';
      sprintf(oname,"%s.%s.hdf",tmp,extension);
    }
  else if (strcmp(inFile+strlen(inFile)-5,".hdf5")==0)
    {
      char tmp[1024];
      strcpy(tmp,inFile);
      tmp[strlen(tmp)-5]='\0';
      sprintf(oname,"%s.%s.hdf",tmp,extension);
    }
  else if (strcmp(inFile+strlen(inFile)-6,".sdhdf")==0)
    {
      char tmp[1024];
      strcpy(tmp,inFile);
      tmp[strlen(tmp)-6]='\0';
      sprintf(oname,"%s.%s.hdf",tmp,extension);
    }
  else
    sprintf(oname,"%s.%s",inFile,extension);
}


// Fix underscores for GIZA
void sdhdf_fixUnderscore(char *input,char *output)
{
  char *tok;
  char tstr[1024];
  int pos=0;
  strcpy(tstr,input);
  tok = strtok(tstr,"_");
  strcpy(output,"");
  while (tok!=NULL)
    {
      if (pos>0)
	strcat(output,"\\_");
      strcat(output,tok);
      tok = strtok(NULL,"_");
      pos=1;
    }
}

//
// Load all available metadata including attributes from the SDHDF file
//
void sdhdf_loadMetaData(sdhdf_fileStruct *inFile)
{
	char beam_path[MAX_STRLEN];
	char band_path[MAX_STRLEN];
	char meta_path[MAX_STRLEN];

	// **NOTE: sdhdf_loadGroup currently does not assign
	// memory properly and causes segfaults**

  // file
	//if (DEBUG == 1) printf(">LOADING FILE\n");
	//sdhdf_loadGroup(inFile, "/", &inFile->fileAttr);
	//if (DEBUG == 1) printf("<FINISHED FILE\n\n");

	// /metadata group
	//if (DEBUG == 1) printf(">LOADING /METADATA\n");
	// TODO: currently not getting attributres from group
  //sdhdf_loadGroup(inFile, METADATA_GRP, &inFile->metaAttr);
	//if (DEBUG == 1) printf("<FINISHED /METADATA\n\n");

	// abort trap errors from above. Memory allocation issue?


	if (DEBUG == 1) printf(">LOADING PRIMARY HEADER\n");
  sdhdf_loadPrimaryHeader(inFile);
	if (DEBUG == 1) printf("<FINISHED PRIMARY HEADER\n\n");

	//exit(1); // TEMP

	if (DEBUG == 1) printf(">LOADING BEAM HEADER\n");
	sdhdf_loadBeamHeader(inFile);
	if (DEBUG == 1) printf("<FINISHED BEAM HEADER\n\n");

	if (DEBUG == 1) printf(">LOADING HISTORY\n");
  sdhdf_loadHistory(inFile);
	if (DEBUG == 1) printf("<FINISHED HISTORY\n\n");

	if (DEBUG == 1) printf(">LOADING SOFTWARE\n");
  sdhdf_loadSoftware(inFile);
	if (DEBUG == 1) printf("<FINISHED SOFTWARE\n\n");

	if (DEBUG == 1) printf(">LOADING SCHEDULE\n");
  sdhdf_loadSchedule(inFile);
	if (DEBUG == 1) printf("<FINISHED SCHEDULE\n\n");

	// /configuration group
	if (DEBUG == 1) printf(">LOADING CONFIG\n");
	// TODO: currently not getting attributres from group
  //sdhdf_loadGroup(inFile, CONFIG_GRP, &inFile->configAttr);
  //sdhdf_loadConfig(inFile); // done
	if (DEBUG == 1) printf("<FINISHED CONFIG\n\n");

	// /beam group
	//if (DEBUG == 1) printf(">LOADING BEAM\n");
	// TODO: currently not getting attributres from group
	// default: beam_0
	//sprintf(beam_path, "%s_0", BEAM_GRP);
  //sdhdf_loadGroup(inFile, beam_path, &inFile->beamAttr);
	//if (DEBUG == 1) printf("<FINISHED BEAM\n\n");

	// /beam/metadata group
	//if (DEBUG == 1) printf(">LOADING BEAM/METADATA\n");
	// TODO: currently not getting attributres from group
	//sprintf(beam_path, "%s_0/%s", BEAM_GRP, METADATA_GRP);
  //sdhdf_loadGroup(inFile, beam_path, &inFile->metaAttr);
	//if (DEBUG == 1) printf("<FINISHED BEAM/METADATA\n\n");

	if (DEBUG == 1) printf(">LOADING BAND HEADER (ASTRO)\n");
	sdhdf_loadBandHeader(inFile,1);
	if (DEBUG == 1) printf("<FINISHED BAND HEADER (ASTRO)\n\n");

	// /beam/band group
	//if (DEBUG == 1) printf(">LOADING BEAM/BAND\n");
	// TODO: currently not getting attributres from group
	//sprintf(band_path, "%s_0/%s_SB0", BEAM_GRP, BAND_GRP);
  //sdhdf_loadGroup(inFile, band_path, &inFile->bandAttr); // todo
	//if (DEBUG == 1) printf("<FINISHED BEAM/BAND\n");

	// /beam/band/metadata group
	//if (DEBUG == 1) printf(">LOADING BEAM/BAND/METADATA\n");
	// TODO: currently not getting attributres from group
	//sprintf(meta_path, "%s_0/%s_SB0/%s", BEAM_GRP, BAND_GRP, METADATA_GRP);
  //sdhdf_loadGroup(inFile, meta_path, &inFile->metaAttr); // todo
	//if (DEBUG == 1) printf("<FINISHED BEAM/BAND/METADATA\n");

	// /beam/band/metadata/ datasets
  if (DEBUG == 1) printf(">LOADING OBS HEADER (ASTRO)\n");
	sdhdf_loadObsHeader(inFile,1);
	if (DEBUG == 1) printf("<FINISHED OBS HEADER (ASTRO)\n\n");

  if (inFile->beam[0].bandData[0].astro_obsHeader[0].dtime == -1) // Check if dump time not in the observation parameters
    {
      int i,j,b;
      for (b=0;b<inFile->nBeam;b++)
	{
	  for (i=0;i<inFile->beam[b].nBand;i++)
	    {
	      for (j=0;j<inFile->beam[b].bandHeader[i].ndump;j++)
		inFile->beam[b].bandData[i].astro_obsHeader[j].dtime = inFile->beam[b].bandHeader[i].dtime;
	    }
	}
    }

  if (strcmp(inFile->primary[0].cal_mode,"ON")==0)
    {
			if (DEBUG == 1) printf(">LOADING BAND HEADER (CAL)\n");
      sdhdf_loadBandHeader(inFile,2);
			if (DEBUG == 1) printf("<FINISHED BAND HEADER (CAL)\n\n");

			if (DEBUG == 1) printf(">LOADING OBS HEADER (CAL)\n");
      sdhdf_loadObsHeader(inFile,2);
			if (DEBUG == 1) printf("<FINISHED LOADING OBS HEADER (CAL)\n\n");
    }

}

void sdhdf_loadGroup(sdhdf_fileStruct *inFile, char *grp, sdhdf_attributes_struct *attrStruct)
{
	char groupName[128];

	sprintf(groupName,"%s", grp);

	if (sdhdf_checkGroupExists(inFile,groupName) == 1)
    {
      printf("WARNING: No group %s in SDHDF file\n", groupName);
    }
  else
    {
			// load attributes
			sdhdf_loadGroupAttributes(inFile, groupName, &attrStruct);
		}
		//printf("finished loadgroupattr\n");
}

/*void sdhdf_loadFile(sdhdf_fileStruct *inFile)
{
  hid_t header_id, stid, val_tid;
	herr_t status;
  int ndims;
  hsize_t dims[2];
	char groupName[128];

	sprintf(groupName,"%s","/");

	// load attributes
	sdhdf_loadAttributes(inFile, groupName, &inFile->fileAttr);

}*/

void sdhdf_loadConfig(sdhdf_fileStruct *inFile)
{
  hid_t header_id, backend_header_id, frontend_header_id, telescope_header_id, stid, val_tid;
  herr_t status;
  int ndims;
  hsize_t dims[2];
  sdhdf_backendConfigStruct *backend_configVals;
	sdhdf_frontendConfigStruct *frontend_configVals;
	sdhdf_telescopeConfigStruct *telescope_configVals;
  char groupName[128];
	char backend_cfg[128];
	char frontend_cfg[128];
	char telescope_cfg[128];

  sprintf(groupName,"%s",CONFIG_GRP);
	sprintf(backend_cfg,"%s/%s",CONFIG_GRP,BACKEND_CONFIG);
	sprintf(frontend_cfg,"%s/%s",CONFIG_GRP,FRONTEND_CONFIG);
	sprintf(telescope_cfg,"%s/%s",CONFIG_GRP,TELESCOPE_CONFIG);

  if (sdhdf_checkGroupExists(inFile,groupName) == 1)
    {
      printf("WARNING: No configuration group in SDHDF file\n");
    }
  else
    {

      backend_header_id = H5Dopen2(inFile->fileID,backend_cfg,H5P_DEFAULT);
			frontend_header_id = H5Dopen2(inFile->fileID,frontend_cfg,H5P_DEFAULT);
			telescope_header_id = H5Dopen2(inFile->fileID,telescope_cfg,H5P_DEFAULT);

      backend_configVals = (sdhdf_backendConfigStruct *)malloc(sizeof(sdhdf_backendConfigStruct)*dims[0]);
			frontend_configVals = (sdhdf_frontendConfigStruct *)malloc(sizeof(sdhdf_frontendConfigStruct)*dims[0]);
			telescope_configVals = (sdhdf_telescopeConfigStruct *)malloc(sizeof(sdhdf_telescopeConfigStruct)*dims[0]);

			// instrument_configuration
      val_tid = H5Tcreate(H5T_COMPOUND,sizeof(sdhdf_backendConfigStruct));
      stid    = H5Tcopy(H5T_C_S1);
      status  = H5Tset_size(stid,64); // Should set to value defined in sdhdf_v1.9.h
      H5Tinsert(val_tid,"BACKEND_PHASE",HOFFSET(sdhdf_backendConfigStruct,backend_phase),stid);
      H5Tinsert(val_tid,"CAL_FREQ",HOFFSET(sdhdf_backendConfigStruct,cal_freq),stid);
      H5Tinsert(val_tid,"CAL_EPOCH",HOFFSET(sdhdf_backendConfigStruct,cal_epoch),stid);
      H5Tinsert(val_tid,"CAL_DUTY_CYCLE",HOFFSET(sdhdf_backendConfigStruct,cal_duty_cycle),stid);
      H5Tinsert(val_tid,"CAL_PHASE",HOFFSET(sdhdf_backendConfigStruct,cal_phase),stid);
      status  = H5Dread(backend_header_id,val_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,backend_configVals);

      strcpy(inFile->cal_epoch,backend_configVals[0].cal_epoch);
      sscanf(backend_configVals[0].cal_freq,"%lf",&(inFile->cal_freq));
      sscanf(backend_configVals[0].cal_phase,"%lf",&(inFile->cal_phase));
      sscanf(backend_configVals[0].cal_duty_cycle,"%lf",&(inFile->cal_duty_cycle));

			status = H5Tclose(val_tid);
      status = H5Tclose(stid);
      status = H5Dclose(backend_header_id);

			// receiver_configuration
			val_tid = H5Tcreate(H5T_COMPOUND,sizeof(sdhdf_frontendConfigStruct));
      stid    = H5Tcopy(H5T_C_S1);
      status  = H5Tset_size(stid,64); // Should set to value defined in sdhdf_v1.9.h
      H5Tinsert(val_tid,"RECEIVER",HOFFSET(sdhdf_frontendConfigStruct,receiver),stid);
      H5Tinsert(val_tid,"LOW_FREQUENCY",HOFFSET(sdhdf_frontendConfigStruct,low_freq),stid);
      H5Tinsert(val_tid,"HIGH_FREQUENCY",HOFFSET(sdhdf_frontendConfigStruct,high_freq),stid);
      H5Tinsert(val_tid,"POLARISATION",HOFFSET(sdhdf_frontendConfigStruct,polarisation),stid);
      H5Tinsert(val_tid,"HAND",HOFFSET(sdhdf_frontendConfigStruct,rx_hand),stid);
      status  = H5Dread(frontend_header_id,val_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,frontend_configVals);

			status = H5Tclose(val_tid);
      status = H5Tclose(stid);
      status = H5Dclose(frontend_header_id);

			// telescope_configuration
			val_tid = H5Tcreate(H5T_COMPOUND,sizeof(sdhdf_telescopeConfigStruct));
      stid    = H5Tcopy(H5T_C_S1);
      status  = H5Tset_size(stid,64); // Should set to value defined in sdhdf_v1.9.h
      H5Tinsert(val_tid,"TELESCOPE",HOFFSET(sdhdf_telescopeConfigStruct,telescope),stid);
      H5Tinsert(val_tid,"ITRF_X_COORDINATE",HOFFSET(sdhdf_telescopeConfigStruct,itrf_x_coord),stid);
      H5Tinsert(val_tid,"ITRF_Y_COORDINATE",HOFFSET(sdhdf_telescopeConfigStruct,itrf_y_coord),stid);
      H5Tinsert(val_tid,"ITRF_Z_COORDINATE",HOFFSET(sdhdf_telescopeConfigStruct,itrf_z_coord),stid);
      H5Tinsert(val_tid,"ITRF_REFERENCE_EPOCH",HOFFSET(sdhdf_telescopeConfigStruct,itrf_ref),stid);
			H5Tinsert(val_tid,"NATIVE_COORDINATE_SYSTEM",HOFFSET(sdhdf_telescopeConfigStruct,native_coord),stid);
      status = H5Dread(telescope_header_id,val_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,telescope_configVals);

			status = H5Tclose(val_tid);
      status = H5Tclose(stid);
      status = H5Dclose(telescope_header_id);

			// load attributes
			sdhdf_loadAttributes(inFile, groupName, &inFile->configAttr);
			sdhdf_loadAttributes(inFile, backend_cfg, &inFile->backend_configAttr);
			sdhdf_loadAttributes(inFile, frontend_cfg, &inFile->frontend_configAttr);
			sdhdf_loadAttributes(inFile, telescope_cfg, &inFile->telescope_configAttr);

      free(backend_configVals);
			free(frontend_configVals);
			free(telescope_configVals);

    }
}

void sdhdf_loadPrimaryHeader(sdhdf_fileStruct *inFile)
{
  int i,j;
  hid_t header_id,headerT,space,stid,val_tid;
  herr_t status;
  int ndims;
  hsize_t dims[2];
	char dsetPthName[MAX_STRLEN];

	sprintf(dsetPthName,"%s/%s",METADATA_GRP,PRIMARY_HEADER);
  if (sdhdf_checkGroupExists(inFile,dsetPthName) == 1)
    {
      printf("WARNING: We do not have %s in the data file %s\n",dsetPthName,inFile->fname);
    }
  else
    {
      header_id  = H5Dopen2(inFile->fileID,dsetPthName,H5P_DEFAULT);
      headerT    = H5Dget_type(header_id);
      space      = H5Dget_space(header_id);
      ndims      = H5Sget_simple_extent_dims(space,dims,NULL);

      if (inFile->primaryAllocatedMemory == 0)
	{
	  inFile->primary = (sdhdf_primaryHeaderStruct *)malloc(sizeof(sdhdf_primaryHeaderStruct)*dims[0]);
	  inFile->nPrimary = dims[0];
	  inFile->primaryAllocatedMemory = 1;
	  if (inFile->nPrimary != 1)
	    {
	      printf("ERROR: UNEXPECTED SIZE IN PRIMARY HEADER\n");
	      exit(1);
	    }
	}

      val_tid = H5Tcreate(H5T_COMPOUND,sizeof(sdhdf_primaryHeaderStruct));
      stid    = H5Tcopy(H5T_C_S1);
      status  = H5Tset_size(stid,20); // Should set to value defined in sdhdf_v1.9.h

      H5Tinsert(val_tid,DATE,HOFFSET(sdhdf_primaryHeaderStruct,date),stid);
      H5Tinsert(val_tid,HDR_DEFN,HOFFSET(sdhdf_primaryHeaderStruct,hdr_defn),stid);
      H5Tinsert(val_tid,HDR_DEFN_VERSION,HOFFSET(sdhdf_primaryHeaderStruct,hdr_defn_version),stid);
      H5Tinsert(val_tid,FILE_FORMAT,HOFFSET(sdhdf_primaryHeaderStruct,file_format),stid);
      H5Tinsert(val_tid,FILE_FORMAT_VERSION,HOFFSET(sdhdf_primaryHeaderStruct,file_format_version),stid);
      H5Tinsert(val_tid,SCHED_BLOCK_ID,HOFFSET(sdhdf_primaryHeaderStruct,sched_block_id),H5T_NATIVE_INT);
      H5Tinsert(val_tid,CAL_MODE,HOFFSET(sdhdf_primaryHeaderStruct,cal_mode),stid);
      H5Tinsert(val_tid,INSTRUMENT,HOFFSET(sdhdf_primaryHeaderStruct,instrument),stid);
      H5Tinsert(val_tid,OBSERVER,HOFFSET(sdhdf_primaryHeaderStruct,observer),stid);
      H5Tinsert(val_tid,PID,HOFFSET(sdhdf_primaryHeaderStruct,pid),stid);
      H5Tinsert(val_tid,RECEIVER,HOFFSET(sdhdf_primaryHeaderStruct,rcvr),stid);

      status = H5Tset_size(stid,64); // Should set to value defined in sdhdf_v1.9.h
      H5Tinsert(val_tid,TELESCOPE,HOFFSET(sdhdf_primaryHeaderStruct,telescope),stid);
      H5Tinsert(val_tid,UTC_START,HOFFSET(sdhdf_primaryHeaderStruct,utc0),stid);

      H5Tinsert(val_tid,N_BEAMS,HOFFSET(sdhdf_primaryHeaderStruct,nbeam),H5T_NATIVE_INT);

      status  = H5Dread(header_id,val_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,inFile->primary);
      inFile->nBeam    = inFile->primary[0].nbeam;

      // Load attributes
			sdhdf_loadAttributes(inFile, dsetPthName, &inFile->primaryAttr);

      status = H5Tclose(val_tid);
      status = H5Tclose(stid);
      status = H5Dclose(header_id);
    }
}


void sdhdf_allocateBeamMemory(sdhdf_fileStruct *inFile,int nbeam)
{
  int i;
  if (inFile->beamAllocatedMemory == 0)
    {
      inFile->beam = (sdhdf_beamStruct *)malloc(sizeof(sdhdf_beamStruct)*nbeam);
      inFile->beamHeader = (sdhdf_beamHeaderStruct *)malloc(sizeof(sdhdf_beamHeaderStruct)*nbeam);
      inFile->nBeam = nbeam;
      inFile->beamAllocatedMemory = 1;
      inFile->nBeamHeaderAttributes=0;
      for (i=0;i<nbeam;i++)
	{
	  inFile->beam[i].nBand=0;
	  inFile->beam[i].bandAllocatedMemory=0;
	  inFile->beam[i].nCalBand=0;
	  inFile->beam[i].calBandAllocatedMemory=0;

	  inFile->beam[i].nBandHeaderAttributes=0;
	  inFile->beam[i].nBandDataAttributes=0;
	}
    }
}

void sdhdf_loadBeamHeader(sdhdf_fileStruct *inFile)
{
  int   i,j;
  hid_t header_id,headerT,space,stid,val_tid;
  herr_t status;
  int ndims;
  hsize_t dims[2];
  int nbeam;
  char label[MAX_STRLEN];

  nbeam = inFile->nBeam;

  if (inFile->beamAllocatedMemory == 0)
    {
      if (nbeam < 0)
	{
	  printf("ERROR: Unexpected number of beams >%d<.\n",nbeam);
	  exit(1);
	}
      inFile->beam       = (sdhdf_beamStruct *)malloc(sizeof(sdhdf_beamStruct)*nbeam);
      inFile->beamHeader = (sdhdf_beamHeaderStruct *)malloc(sizeof(sdhdf_beamHeaderStruct)*nbeam);
      inFile->nBeam      = nbeam;
      inFile->beamAllocatedMemory = 1;
    }

  sprintf(label,"%s/%s",METADATA_GRP,BEAM_PARAMS);
  header_id  = H5Dopen2(inFile->fileID,label,H5P_DEFAULT);
  headerT    = H5Dget_type(header_id);
  space      = H5Dget_space(header_id);
  ndims      = H5Sget_simple_extent_dims(space,dims,NULL);

  if (dims[0] != nbeam)
    {
      printf("ERROR: wrong number of beams in %s\n",label);
      printf("nbeam = %d\n",nbeam);
      printf("ndims = %d\n",ndims);
      exit(1);
    }

  val_tid = H5Tcreate(H5T_COMPOUND,sizeof(sdhdf_beamHeaderStruct));
  stid    = H5Tcopy(H5T_C_S1);
  status  = H5Tset_size(stid,MAX_STRLEN); // Should set to value defined in sdhdf_v1.9.h

  H5Tinsert(val_tid,LABEL,HOFFSET(sdhdf_beamHeaderStruct,label),stid);
  H5Tinsert(val_tid,N_BANDS,HOFFSET(sdhdf_beamHeaderStruct,nBand),H5T_NATIVE_INT);
  H5Tinsert(val_tid,SOURCE,HOFFSET(sdhdf_beamHeaderStruct,source),stid);

  status  = H5Dread(header_id,val_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,(inFile->beamHeader));
  for (i=0;i<nbeam;i++)
    {
      inFile->beam[i].bandAllocatedMemory = 0;
      inFile->beam[i].calBandAllocatedMemory = 0;
      inFile->beam[i].nBand = inFile->beamHeader[i].nBand;
    }

  status = H5Tclose(val_tid);
  status = H5Tclose(stid);
  status = H5Dclose(header_id);

  // Load attributes
	sdhdf_loadAttributes(inFile, label, &inFile->beamHeaderAttr);


}

void sdhdf_setMetadataDefaults(sdhdf_primaryHeaderStruct *primaryHeader,sdhdf_beamHeaderStruct *beamHeader,
			       sdhdf_bandHeaderStruct *bandHeader,sdhdf_softwareVersionsStruct *softwareVersions,sdhdf_historyStruct *history,
			       int nbeam,int nband)
{
  int i;

  // History
  strcpy(history->date,"UNKNOWN");
  strcpy(history->proc_name,"UNKNOWN");
  strcpy(history->proc_descr,"UNKNOWN");
  strcpy(history->proc_args,"UNKNOWN");
  strcpy(history->proc_host,"UNKNOWN");

  // Software versions
  strcpy(softwareVersions->proc_name,"UNKNOWN");
  strcpy(softwareVersions->software,"UNKNONW");
  strcpy(softwareVersions->software_descr,"UNKNOWN");
  strcpy(softwareVersions->software_version,"UNKNOWN");

	// Schedule
  /*strcpy(schedHeader->date,"UNKNOWN");
  strcpy(schedHeader->host,"UNKNOWN");
	strcpy(schedHeader->host_name,"UNKNOWN");
	strcpy(schedHeader->api,"UNKNOWN");
	strcpy(schedHeader->version,"UNKNOWN");
	strcpy(schedHeader->meta,"UNKNOWN");*/

  // Primary header
  strcpy(primaryHeader->date,"UNKNOWN");
  strcpy(primaryHeader->hdr_defn,"UNKNOWN");
  strcpy(primaryHeader->hdr_defn_version,"UNKNOWN");
  strcpy(primaryHeader->file_format,"UNKNOWN");
  strcpy(primaryHeader->file_format_version,"UNKNOWN");
  primaryHeader->sched_block_id=-1;
  strcpy(primaryHeader->cal_mode,"UNKNOWN");
  strcpy(primaryHeader->instrument,"UNKNOWN");
  strcpy(primaryHeader->observer,"UNKNOWN");
  strcpy(primaryHeader->pid,"UNKNOWN");
  strcpy(primaryHeader->rcvr,"UNKNOWN");
  strcpy(primaryHeader->telescope,"UNKNOWN");
  strcpy(primaryHeader->utc0,"UNKNOWN");
  primaryHeader->nbeam=nbeam;

  // Should setup multiple beams correctly
  strcpy(beamHeader->label,"beam_0");
  beamHeader->nBand=1;
  strcpy(beamHeader->source,"SOURCE NAME");

  strcpy(bandHeader->label,"band_0");
  bandHeader->fc = 0;
  bandHeader->f0 = 0;
  bandHeader->f1 = 0;
  bandHeader->nchan = 0;
  bandHeader->npol = 0;
  strcpy(bandHeader->pol_type,"UNKNOWN");
  bandHeader->dtime = 0;
  bandHeader->ndump = 0;

}

void sdhdf_writeSoftwareVersions(sdhdf_fileStruct *outFile,sdhdf_softwareVersionsStruct *outParams)
{
  hid_t dset_id,datatype_id,group_id;
  herr_t status;
  hid_t dataspace_id,stid;
  hsize_t dims[1];
  char name[1024];
  char groupName[1024];
	char dsetPthName[1024];

  dims[0] = 1;

  datatype_id = H5Tcreate (H5T_COMPOUND, sizeof(sdhdf_softwareVersionsStruct));
  stid = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(stid,64); // Should set to value defined in sdhdf_v1.9.h

  H5Tinsert(datatype_id,PROC,HOFFSET(sdhdf_softwareVersionsStruct,proc_name),stid);
  H5Tinsert(datatype_id,SOFTWARE,HOFFSET(sdhdf_softwareVersionsStruct,software),stid);
  H5Tinsert(datatype_id,SOFTWARE_DESCR,HOFFSET(sdhdf_softwareVersionsStruct,software_descr),stid);
  H5Tinsert(datatype_id,SOFTWARE_VERSION,HOFFSET(sdhdf_softwareVersionsStruct,software_version),stid);

  dataspace_id = H5Screate_simple(1,dims,NULL);

  sprintf(groupName,"%s",METADATA_GRP);
	sprintf(dsetPthName,"%s/%s",METADATA_GRP,SOFTWARE_VERSIONS);
  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
    {
      group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Gclose(group_id);
    }
  sprintf(name,dsetPthName);
  if (sdhdf_checkGroupExists(outFile,name) == 1)
      dset_id = H5Dcreate2(outFile->fileID,name,datatype_id,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  else
      dset_id  = H5Dopen2(outFile->fileID,name,H5P_DEFAULT);
  printf("Created %s\n",SOFTWARE_VERSIONS);
  status  = H5Dwrite(dset_id,datatype_id,H5S_ALL,H5S_ALL,H5P_DEFAULT,outParams);
  status  = H5Dclose(dset_id);

  status  = H5Sclose(dataspace_id);
  status  = H5Tclose(datatype_id);

  status  = H5Tclose(stid);
}

void sdhdf_writeBeamHeader(sdhdf_fileStruct *outFile,sdhdf_beamHeaderStruct *beamHeader,int nBeams)
{
  hid_t dset_id,datatype_id,group_id;
  herr_t status;
  hid_t dataspace_id,stid;
  hsize_t dims[1];
  char name[1024];
  char groupName[1024];
	char dsetPthName[1024];

  dims[0] = nBeams;

  datatype_id = H5Tcreate (H5T_COMPOUND, sizeof(sdhdf_beamHeaderStruct));
  stid = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(stid,MAX_STRLEN);
  printf("Setting size to %d with status %d\n",MAX_STRLEN,status);
  H5Tinsert(datatype_id,LABEL,HOFFSET(sdhdf_beamHeaderStruct,label),stid);
  H5Tinsert(datatype_id,N_BANDS,HOFFSET(sdhdf_beamHeaderStruct,nBand),H5T_NATIVE_INT);
  H5Tinsert(datatype_id,SOURCE,HOFFSET(sdhdf_beamHeaderStruct,source),stid);
  printf("In beam write with %s\n",beamHeader[0].label);
  dataspace_id = H5Screate_simple(1,dims,NULL);

  sprintf(groupName,"%s",METADATA_GRP);
	sprintf(dsetPthName,"%s/%s",METADATA_GRP,BEAM_PARAMS);
  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
    {
      group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Gclose(group_id);
    }
  sprintf(name,dsetPthName);
  if (sdhdf_checkGroupExists(outFile,name) == 1)
    dset_id = H5Dcreate2(outFile->fileID,name,datatype_id,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  else
    {
      dset_id  = H5Dopen2(outFile->fileID,name,H5P_DEFAULT);
    }

  status  = H5Dwrite(dset_id,datatype_id,H5S_ALL,H5S_ALL,H5P_DEFAULT,beamHeader);

  status  = H5Dclose(dset_id);
  status  = H5Sclose(dataspace_id);
  status  = H5Tclose(datatype_id);
  status  = H5Tclose(stid);
}

void sdhdf_writePrimaryHeader(sdhdf_fileStruct *outFile,sdhdf_primaryHeaderStruct *primaryHeader)
{
  hid_t dset_id,datatype_id,group_id;
  herr_t status;
  hid_t dataspace_id,stid;
  hsize_t dims[1];
  char name[1024];
  char groupName[1024];
	char dsetPthName[1024];

  dims[0] = 1;

  datatype_id = H5Tcreate (H5T_COMPOUND, sizeof(sdhdf_primaryHeaderStruct));
  stid = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(stid,20); // Should set to value defined in sdhdf_v1.9.h

  H5Tinsert(datatype_id,DATE,HOFFSET(sdhdf_primaryHeaderStruct,date),stid);
  H5Tinsert(datatype_id,HDR_DEFN,HOFFSET(sdhdf_primaryHeaderStruct,hdr_defn),stid);
  H5Tinsert(datatype_id,HDR_DEFN_VERSION,HOFFSET(sdhdf_primaryHeaderStruct,hdr_defn_version),stid);
  H5Tinsert(datatype_id,FILE_FORMAT,HOFFSET(sdhdf_primaryHeaderStruct,file_format),stid);
  H5Tinsert(datatype_id,FILE_FORMAT_VERSION,HOFFSET(sdhdf_primaryHeaderStruct,file_format_version),stid);
  H5Tinsert(datatype_id,SCHED_BLOCK_ID,HOFFSET(sdhdf_primaryHeaderStruct,sched_block_id),H5T_NATIVE_INT);
  H5Tinsert(datatype_id,CAL_MODE,HOFFSET(sdhdf_primaryHeaderStruct,cal_mode),stid);
  H5Tinsert(datatype_id,INSTRUMENT,HOFFSET(sdhdf_primaryHeaderStruct,instrument),stid);
  H5Tinsert(datatype_id,OBSERVER,HOFFSET(sdhdf_primaryHeaderStruct,observer),stid);
  H5Tinsert(datatype_id,PID,HOFFSET(sdhdf_primaryHeaderStruct,pid),stid);
  H5Tinsert(datatype_id,RECEIVER,HOFFSET(sdhdf_primaryHeaderStruct,rcvr),stid);

  status = H5Tset_size(stid,64); // Should set to value defined in sdhdf_v1.9.h
  H5Tinsert(datatype_id,TELESCOPE,HOFFSET(sdhdf_primaryHeaderStruct,telescope),stid);
  H5Tinsert(datatype_id,UTC_START,HOFFSET(sdhdf_primaryHeaderStruct,utc0),stid);
  H5Tinsert(datatype_id,N_BEAMS,HOFFSET(sdhdf_primaryHeaderStruct,nbeam),H5T_NATIVE_INT);

  dataspace_id = H5Screate_simple(1,dims,NULL);

  sprintf(groupName,"%s",METADATA_GRP);
	sprintf(dsetPthName,"%s/%s",METADATA_GRP,PRIMARY_HEADER);
  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
    {
      group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Gclose(group_id);
    }
  sprintf(name,dsetPthName);
  if (sdhdf_checkGroupExists(outFile,name) == 1)
      dset_id = H5Dcreate2(outFile->fileID,name,datatype_id,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  else
      dset_id  = H5Dopen2(outFile->fileID,name,H5P_DEFAULT);
  printf("Created %s\n",PRIMARY_HEADER);
  status  = H5Dwrite(dset_id,datatype_id,H5S_ALL,H5S_ALL,H5P_DEFAULT,primaryHeader);
  status  = H5Dclose(dset_id);

  status  = H5Sclose(dataspace_id);
  status  = H5Tclose(datatype_id);

  status  = H5Tclose(stid);


}

//
// Load the band header meta data in the beam_xx/beam_metadata group
// type = 1 for astronomy
// type = 2 for cal
//
void sdhdf_loadBandHeader(sdhdf_fileStruct *inFile,int type)
{
  int   i,j;
  hid_t header_id,headerT,space,stid,val_tid;
  herr_t status;
  int ndims;
  hsize_t dims[2];
  int nbeam,nband;
  char label[MAX_STRLEN];
  char beamLabel[MAX_STRLEN];

  nbeam = inFile->nBeam;

  if (inFile->beamAllocatedMemory == 0 || nbeam < 1)
    {
      printf("ERROR: Trying to load bands, but no beams allocated\n");
      exit(1);
    }

  for (i=0;i<nbeam;i++)
    {
      strcpy(beamLabel,inFile->beamHeader[i].label);
      nband = inFile->beam[i].nBand;
      if (nband < 1)
	{
	  printf("ERROR: No bands allocated\n");
	  exit(1);
 	}

      if (type==1)
	{
	  if (inFile->beam[i].bandAllocatedMemory == 0)
	    {
	      inFile->beam[i].bandData       = (sdhdf_bandStruct *)malloc(sizeof(sdhdf_bandStruct)*nband);
	      inFile->beam[i].bandHeader     = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*nband);
	      inFile->beam[i].bandAllocatedMemory = 1;
	    }
	}
      else if (type==2)
	{
	  if (inFile->beam[i].calBandAllocatedMemory == 0)
	    {
	      inFile->beam[i].calBandHeader     = (sdhdf_bandHeaderStruct *)malloc(sizeof(sdhdf_bandHeaderStruct)*nband);
	      inFile->beam[i].calBandAllocatedMemory = 1;
	    }
	}
      for (j=0;j<nband;j++)
	{
	  if (type==1)
	    {
	      sdhdf_initialise_bandHeader(&(inFile->beam[i].bandHeader[j]));
	      inFile->beam[i].bandData[j].astro_obsHeaderAllocatedMemory = 0;
	      inFile->beam[i].bandData[j].cal_obsHeaderAllocatedMemory   = 0;
	      inFile->beam[i].bandData[j].nAstro_obsHeader               = -1;
	      inFile->beam[i].bandData[j].nCal_obsHeader                 = -1;
	      sdhdf_initialise_spectralDumps(&(inFile->beam[i].bandData[j].astro_data));
	      sdhdf_initialise_spectralDumps(&(inFile->beam[i].bandData[j].cal_on_data));
	      sdhdf_initialise_spectralDumps(&(inFile->beam[i].bandData[j].cal_off_data));
	      sdhdf_initialise_spectralDumps(&(inFile->beam[i].bandData[j].cal_proc_tsys));
	      sdhdf_initialise_spectralDumps(&(inFile->beam[i].bandData[j].cal_proc_diff_gain));
	      sdhdf_initialise_spectralDumps(&(inFile->beam[i].bandData[j].cal_proc_diff_phase));
	    }
	  else if (type==2)
	    {
	      sdhdf_initialise_bandHeader(&(inFile->beam[i].calBandHeader[j]));
	    }
	}

      if (type==1)
	sprintf(label,"%s/%s/%s",beamLabel,METADATA_GRP,BAND_PARAMS);
      else if (type==2)
	sprintf(label,"%s/%s/%s",beamLabel,METADATA_GRP,CAL_BAND_PARAMS);

      header_id  = H5Dopen2(inFile->fileID,label,H5P_DEFAULT);
      headerT    = H5Dget_type(header_id);
      space      = H5Dget_space(header_id);
      ndims      = H5Sget_simple_extent_dims(space,dims,NULL);

      if (dims[0] != nband)
	{
	  printf("ERROR: incorrect number of bands.  In beam header have %d band. Dimension of data is %d\n",nband,(int)dims[0]);
	  printf("ndims = %d\n",ndims);
	  printf("Label = %s\n",label);
	  printf("File = %s\n",inFile->fname);
	  exit(1);
	}

      val_tid = H5Tcreate(H5T_COMPOUND,sizeof(sdhdf_bandHeaderStruct));
      stid = H5Tcopy(H5T_C_S1);
      status = H5Tset_size(stid,64); // Should set to value defined in sdhdf_v1.9.h

      H5Tinsert(val_tid,LABEL,HOFFSET(sdhdf_bandHeaderStruct,label),stid);
      H5Tinsert(val_tid,C_FREQ,HOFFSET(sdhdf_bandHeaderStruct,fc),H5T_NATIVE_DOUBLE);
      H5Tinsert(val_tid,LOW_FREQ,HOFFSET(sdhdf_bandHeaderStruct,f0),H5T_NATIVE_DOUBLE);
      H5Tinsert(val_tid,HIGH_FREQ,HOFFSET(sdhdf_bandHeaderStruct,f1),H5T_NATIVE_DOUBLE);
      H5Tinsert(val_tid,N_CHANS,HOFFSET(sdhdf_bandHeaderStruct,nchan),H5T_NATIVE_INT);
      H5Tinsert(val_tid,N_POLS,HOFFSET(sdhdf_bandHeaderStruct,npol),H5T_NATIVE_INT);
      status = H5Tset_size(stid,20); // Should set to value defined in sdhdf_v1.9.h
      H5Tinsert(val_tid,POL_TYPE,HOFFSET(sdhdf_bandHeaderStruct,pol_type),stid);
      H5Tinsert(val_tid,DUMP_TIME,HOFFSET(sdhdf_bandHeaderStruct,dtime),H5T_NATIVE_DOUBLE);
      H5Tinsert(val_tid,N_DUMPS,HOFFSET(sdhdf_bandHeaderStruct,ndump),H5T_NATIVE_INT);

      if (type==1)
			{
			  status  = H5Dread(header_id,val_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,inFile->beam[i].bandHeader);
				// Load attributes
			  sdhdf_loadAttributes(inFile, label, &inFile->beam[i].bandHeaderAttr);
			}
			else if (type==2)
			{
				status  = H5Dread(header_id,val_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,inFile->beam[i].calBandHeader);
				// Load attributes
			  sdhdf_loadAttributes(inFile, label, &inFile->beam[i].bandHeaderAttr);
			}

      status  = H5Tclose(val_tid);
      status  = H5Tclose(stid);
      status  = H5Dclose(header_id);
    }
}


void sdhdf_writeBandHeader(sdhdf_fileStruct *outFile,sdhdf_bandHeaderStruct *outBandParams,char *beamLabel,int outBands,int type)
{
  hid_t dset_id,datatype_id,group_id;
  herr_t status;
  hid_t dataspace_id,stid;
  hsize_t dims[1];
  char name[1024];
  char groupName[1024];

  dims[0] = outBands;

  datatype_id = H5Tcreate (H5T_COMPOUND, sizeof(sdhdf_bandHeaderStruct));
  stid = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(stid,12);

  H5Tinsert(datatype_id,LABEL,HOFFSET(sdhdf_bandHeaderStruct,label),stid);
  H5Tinsert(datatype_id,C_FREQ,HOFFSET(sdhdf_bandHeaderStruct,fc),H5T_NATIVE_DOUBLE);
  H5Tinsert(datatype_id,LOW_FREQ,HOFFSET(sdhdf_bandHeaderStruct,f0),H5T_NATIVE_DOUBLE);
  H5Tinsert(datatype_id,HIGH_FREQ,HOFFSET(sdhdf_bandHeaderStruct,f1),H5T_NATIVE_DOUBLE);
  H5Tinsert(datatype_id,N_CHANS,HOFFSET(sdhdf_bandHeaderStruct,nchan),H5T_NATIVE_INT);
  H5Tinsert(datatype_id,N_POLS,HOFFSET(sdhdf_bandHeaderStruct,npol),H5T_NATIVE_INT);
  status = H5Tset_size(stid,20); // Should set to value defined in sdhdf_v1.9.h
  H5Tinsert(datatype_id,POL_TYPE,HOFFSET(sdhdf_bandHeaderStruct,pol_type),stid);
  H5Tinsert(datatype_id,DUMP_TIME,HOFFSET(sdhdf_bandHeaderStruct,dtime),H5T_NATIVE_DOUBLE);
  H5Tinsert(datatype_id,N_DUMPS,HOFFSET(sdhdf_bandHeaderStruct,ndump),H5T_NATIVE_INT);
  dataspace_id = H5Screate_simple(1,dims,NULL);

  // Do we need to create the groups
  sprintf(groupName,"%s",beamLabel);
  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
    {
      group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Gclose(group_id);
    }

  sprintf(groupName,"%s/%s",beamLabel,METADATA_GRP);
  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
    {
      group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Gclose(group_id);
    }
  if (type==1)
    {
      sprintf(name,"%s/%s/%s",beamLabel,METADATA_GRP,BAND_PARAMS);
      if (sdhdf_checkGroupExists(outFile,name) == 1)
	{
	  dset_id = H5Dcreate2(outFile->fileID,name,datatype_id,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	}
      else
	dset_id  = H5Dopen2(outFile->fileID,name,H5P_DEFAULT);
    }
  else
    {
      sprintf(name,"%s/%s/%s",beamLabel,METADATA_GRP,CAL_BAND_PARAMS);
      if (sdhdf_checkGroupExists(outFile,name) == 1)
	dset_id = H5Dcreate2(outFile->fileID,name,datatype_id,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      else
	dset_id  = H5Dopen2(outFile->fileID,name,H5P_DEFAULT);
    }
  status  = H5Dwrite(dset_id,datatype_id,H5S_ALL,H5S_ALL,H5P_DEFAULT,outBandParams);
  status  = H5Dclose(dset_id);
  status  = H5Sclose(dataspace_id);
  status  = H5Tclose(datatype_id);
  status  = H5Tclose(stid);
}

// 0 = already exists
int sdhdf_checkGroupExists(sdhdf_fileStruct *inFile,char *groupName)
{
  herr_t status;
  // Check if group exists
  status = H5Eset_auto1(NULL,NULL);
  status = H5Gget_objinfo(inFile->fileID,groupName,0,NULL);
  if (status==0)
    return 0;
  else
    return 1;
}

//
// Load the band observation meta data in the beam_xx/band_yy/band_metadata group
// type = 1: astronomy meta-data
// type = 2: cal meta-data
//
void sdhdf_loadObsHeader(sdhdf_fileStruct *inFile,int type)
{
  int   i,j,k;
  hid_t header_id,headerT,space,stid,val_tid;
  herr_t status;
  int ndims;
  hsize_t dims[2];
  int nbeam,nband,ndump;
  char label[MAX_STRLEN];
  char beamLabel[MAX_STRLEN];

  nbeam = inFile->nBeam;

  if (inFile->beamAllocatedMemory == 0 || nbeam < 1)
    {
      printf("ERROR: Trying to load observation metadata, but no beams allocated\n");
      exit(1);
    }

  for (i=0;i<nbeam;i++)
    {
      strcpy(beamLabel,inFile->beamHeader[i].label);
      nband = inFile->beam[i].nBand;

      if (inFile->beam[i].bandAllocatedMemory == 0 || nband < 1)
	{
	  printf("ERROR: Trying to load observation metadata, but no bands allocated\n");
	  exit(1);
	}

      for (j=0;j<nband;j++)
	{
	  if (type==1)
	    ndump = inFile->beam[i].bandHeader[j].ndump;
	  else if (type==2)
	    ndump = inFile->beam[i].calBandHeader[j].ndump;


	  if (ndump < 1)
	    {
	      printf("ERROR: Trying to load observation metadata, but ndump = 0\n");
	      exit(1);
	    }

	  if (type==1)
	    {
	      if (inFile->beam[i].bandData[j].astro_obsHeaderAllocatedMemory == 0)
		{
		  inFile->beam[i].bandData[j].astro_obsHeader = (sdhdf_obsParamsStruct *)malloc(sizeof(sdhdf_obsParamsStruct)*ndump);
		  inFile->beam[i].bandData[j].astro_obsHeaderAllocatedMemory = 1;
		  inFile->beam[i].bandData[j].nAstro_obsHeader = ndump;
		}
	    }
	  else if (type==2)
	    {
	      if (inFile->beam[i].bandData[j].cal_obsHeaderAllocatedMemory == 0)
		{
		  inFile->beam[i].bandData[j].cal_obsHeader = (sdhdf_obsParamsStruct *)malloc(sizeof(sdhdf_obsParamsStruct)*ndump);
		  inFile->beam[i].bandData[j].cal_obsHeaderAllocatedMemory = 1;
		  inFile->beam[i].bandData[j].nCal_obsHeader = ndump;
		}
	    }

	  for (k=0;k<ndump;k++)
	    {
	      if (type==1)
		sdhdf_initialise_obsHeader(&(inFile->beam[i].bandData[j].astro_obsHeader[k]));
	      else
		sdhdf_initialise_obsHeader(&(inFile->beam[i].bandData[j].cal_obsHeader[k]));
	    }

	  if (type==1)
	    sprintf(label,"%s/%s/%s/%s",beamLabel,inFile->beam[i].bandHeader[j].label,METADATA_GRP,OBS_PARAMS);
	  else if (type==2)
	    sprintf(label,"%s/%s/%s/%s",beamLabel,inFile->beam[i].bandHeader[j].label,METADATA_GRP,CAL_OBS_PARAMS);

	  header_id  = H5Dopen2(inFile->fileID,label,H5P_DEFAULT);
	  headerT    = H5Dget_type(header_id);
	  space      = H5Dget_space(header_id);
	  ndims      = H5Sget_simple_extent_dims(space,dims,NULL);
	  if (ndims == -1)
	    {
	      printf("ERROR: number of dimensions is -1 in %s\n",label);
	      exit(1);
	    }
	  if (dims[0] != ndump)
	    {
	      printf("ERROR: [%s] Missing data detected. In %s number of dumps %d (expected) !== %d (actual) (band = %d) [%s]. ndims = %d, dims[0] = %d, type = %d\n",inFile->fname,inFile->beam[i].bandHeader[j].label,ndump,(int)dims[0],j,label,ndims,(int)dims[0],(int)type);
	      printf("ATTEMPTING TO CONTINUE ....\n");
	    }
	  else
	    {
	      val_tid = H5Tcreate(H5T_COMPOUND,sizeof(sdhdf_obsParamsStruct));
	      stid = H5Tcopy(H5T_C_S1);
	      status = H5Tset_size(stid,64); // Should set to value defined in sdhdf_v1.9.h

	      H5Tinsert(val_tid,"ELAPSED_TIME",HOFFSET(sdhdf_obsParamsStruct,timeElapsed),H5T_NATIVE_DOUBLE);
	      H5Tinsert(val_tid,"INTEGRATION_TIME",HOFFSET(sdhdf_obsParamsStruct,dtime),H5T_NATIVE_DOUBLE);
				// deprecated in v4
				//H5Tinsert(val_tid,"TIME_DB",HOFFSET(sdhdf_obsParamsStruct,timedb),stid);
	      H5Tinsert(val_tid,"MJD",HOFFSET(sdhdf_obsParamsStruct,mjd),H5T_NATIVE_DOUBLE);
	      H5Tinsert(val_tid,"UTC",HOFFSET(sdhdf_obsParamsStruct,utc),stid);
				// deprecated in v4
	      //H5Tinsert(val_tid,"UT_DATE",HOFFSET(sdhdf_obsParamsStruct,ut_date),stid);
	      H5Tinsert(val_tid,"LOCAL_TIME",HOFFSET(sdhdf_obsParamsStruct,local_time),stid);
				// changed in v4
				//H5Tinsert(val_tid,"RA_STR",HOFFSET(sdhdf_obsParamsStruct,raStr),stid);
				H5Tinsert(val_tid,"RIGHT_ASCENSION",HOFFSET(sdhdf_obsParamsStruct,raStr),stid);
				// changed in v4
	      //H5Tinsert(val_tid,"DEC_STR",HOFFSET(sdhdf_obsParamsStruct,decStr),stid);
				H5Tinsert(val_tid,"DECLINATION",HOFFSET(sdhdf_obsParamsStruct,decStr),stid);
				// deprecated in v4
				//H5Tinsert(val_tid,"RA_DEG",HOFFSET(sdhdf_obsParamsStruct,raDeg),H5T_NATIVE_DOUBLE);
	      //H5Tinsert(val_tid,"DEC_DEG",HOFFSET(sdhdf_obsParamsStruct,decDeg),H5T_NATIVE_DOUBLE);
	      //H5Tinsert(val_tid,"RA_OFFSET",HOFFSET(sdhdf_obsParamsStruct,raOffset),H5T_NATIVE_DOUBLE);
	      //H5Tinsert(val_tid,"DEC_OFFSET",HOFFSET(sdhdf_obsParamsStruct,decOffset),H5T_NATIVE_DOUBLE);
				// added in v4
				H5Tinsert(val_tid,"DRIVE_STATUS",HOFFSET(sdhdf_obsParamsStruct,fstat),H5T_NATIVE_INT);
				// changed in v4
	      //H5Tinsert(val_tid,"GL",HOFFSET(sdhdf_obsParamsStruct,gl),H5T_NATIVE_DOUBLE);
				H5Tinsert(val_tid,"GALACTIC_LONGITUDE",HOFFSET(sdhdf_obsParamsStruct,gl),H5T_NATIVE_DOUBLE);
				// changed in v4
	      //H5Tinsert(val_tid,"GB",HOFFSET(sdhdf_obsParamsStruct,gb),H5T_NATIVE_DOUBLE);
				H5Tinsert(val_tid,"GALACTIC_LATITUDE",HOFFSET(sdhdf_obsParamsStruct,gb),H5T_NATIVE_DOUBLE);
				// changed in v4
	      //H5Tinsert(val_tid,"AZ",HOFFSET(sdhdf_obsParamsStruct,az),H5T_NATIVE_DOUBLE);
				H5Tinsert(val_tid,"AZIMUTH_ANGLE",HOFFSET(sdhdf_obsParamsStruct,az),H5T_NATIVE_DOUBLE);
				// changed in v4
	      //H5Tinsert(val_tid,"ZE",HOFFSET(sdhdf_obsParamsStruct,ze),H5T_NATIVE_DOUBLE);
				H5Tinsert(val_tid,"ZENITH_ANGLE",HOFFSET(sdhdf_obsParamsStruct,ze),H5T_NATIVE_DOUBLE);
				// changed in v4
	      //H5Tinsert(val_tid,"EL",HOFFSET(sdhdf_obsParamsStruct,el),H5T_NATIVE_DOUBLE);
				H5Tinsert(val_tid,"ELEVATION_ANGLE",HOFFSET(sdhdf_obsParamsStruct,el),H5T_NATIVE_DOUBLE);
				// deprecated in v4
				//H5Tinsert(val_tid,"AZ_DRIVE_RATE",HOFFSET(sdhdf_obsParamsStruct,az_drive_rate),H5T_NATIVE_DOUBLE);
	      //H5Tinsert(val_tid,"ZE_DRIVE_RATE",HOFFSET(sdhdf_obsParamsStruct,ze_drive_rate),H5T_NATIVE_DOUBLE);
	      H5Tinsert(val_tid,"HOUR_ANGLE",HOFFSET(sdhdf_obsParamsStruct,hourAngle),H5T_NATIVE_DOUBLE);
				// changed in v4
	      //H5Tinsert(val_tid,"PARA_ANGLE",HOFFSET(sdhdf_obsParamsStruct,paraAngle),H5T_NATIVE_DOUBLE);
				H5Tinsert(val_tid,"PARALLACTIC_ANGLE",HOFFSET(sdhdf_obsParamsStruct,paraAngle),H5T_NATIVE_DOUBLE);
				// changed in v4
	      //H5Tinsert(val_tid,"WIND_DIR",HOFFSET(sdhdf_obsParamsStruct,windDir),H5T_NATIVE_DOUBLE);
	      //H5Tinsert(val_tid,"WIND_SPD",HOFFSET(sdhdf_obsParamsStruct,windSpd),H5T_NATIVE_DOUBLE);
				H5Tinsert(val_tid,"WIND_DIRECTION",HOFFSET(sdhdf_obsParamsStruct,windDir),H5T_NATIVE_DOUBLE);
	      H5Tinsert(val_tid,"WIND_SPEED",HOFFSET(sdhdf_obsParamsStruct,windSpd),H5T_NATIVE_DOUBLE);
				// added in v4
				H5Tinsert(val_tid,"PRESSURE",HOFFSET(sdhdf_obsParamsStruct,pressure),H5T_NATIVE_DOUBLE);
				H5Tinsert(val_tid,"PRESSURE_MSL",HOFFSET(sdhdf_obsParamsStruct,pressureMSL),H5T_NATIVE_DOUBLE);
				H5Tinsert(val_tid,"RELATIVE_HUMIDITY",HOFFSET(sdhdf_obsParamsStruct,relHumidity),H5T_NATIVE_DOUBLE);
				H5Tinsert(val_tid,"TEMPERATURE",HOFFSET(sdhdf_obsParamsStruct,temp),H5T_NATIVE_DOUBLE);

	      if (type==1)
				{
					status  = H5Dread(header_id,val_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,inFile->beam[i].bandData[j].astro_obsHeader);
					// Load attributes
				  sdhdf_loadAttributes(inFile, label, &inFile->beam[i].bandData[j].astro_obsHeaderAttr);
				}
				else if (type==2)
				{
					status  = H5Dread(header_id,val_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,inFile->beam[i].bandData[j].cal_obsHeader);
					// Load attributes
				  sdhdf_loadAttributes(inFile, label, &inFile->beam[i].bandData[j].cal_obsHeaderAttr);
				}

				status  = H5Tclose(val_tid);
	      status  = H5Tclose(stid);
	    }
	  status  = H5Dclose(header_id);
	}
    }
}

void sdhdf_initialise_spectralDumps(sdhdf_spectralDumpsStruct *in)
{
  in->freqAllocatedMemory=0;
  in->flagAllocatedMemory=0;
  in->dataWeightsAllocatedMemory=0;
  in->pol1AllocatedMemory=0;
  in->pol2AllocatedMemory=0;
  in->pol3AllocatedMemory=0;
  in->pol4AllocatedMemory=0;
  in->nFreqAttributes=0;
  in->nDataAttributes=0;
}

void sdhdf_initialise_bandHeader(sdhdf_bandHeaderStruct *header)
{
  strcpy(header->label,"unset");
  header->fc = -1;
  header->f0 = -1;
  header->f1 = -1;
  header->nchan = -1;
  header->npol = -1;
  header->ndump = -1;
}

void sdhdf_initialise_obsHeader(sdhdf_obsParamsStruct *obs)
{
  obs->timeElapsed=-1;
  obs->dtime = -1;
  // deprecated in v4
	//strcpy(obs->timedb,"unset");
  obs->mjd = -1;
  strcpy(obs->utc,"unset");
	// deprecated in v4
  //strcpy(obs->ut_date,"unset");
  strcpy(obs->local_time,"unset");
  strcpy(obs->raStr,"unset");
  strcpy(obs->decStr,"unset");
	// new for v4
	obs->fstat = -1;
	//
	obs->raDeg = -1;
  obs->decDeg = -1;
	// deprecated v4
  //obs->raOffset = -1;
  //obs->decOffset = -1;
  obs->gl = -1;
  obs->gb = -1;
  obs->az = -1;
  obs->ze = -1;
  obs->el = -1;
	// deprecated v4
  //obs->az_drive_rate = -1;
  //obs->ze_drive_rate = -1;
  obs->hourAngle = -1;
  obs->paraAngle = -1;
  obs->windDir = -1;
  obs->windSpd = -1;
	// added for v4
	obs->pressure = -1;
	obs->pressureMSL = -1;
	obs->relHumidity = -1;
	obs->temp = -1;
}

void sdhdf_writeHistory(sdhdf_fileStruct *outFile,sdhdf_historyStruct *outParams,int n)
{
  hid_t dset_id,datatype_id,group_id;
  herr_t status;
  hid_t dataspace_id,stid;
  hsize_t dims[1];
  char groupName[1024];
	char history[1024];

  dims[0] = n;
  datatype_id = H5Tcreate (H5T_COMPOUND, sizeof(sdhdf_historyStruct));
  stid = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(stid,64);
  dataspace_id = H5Screate_simple(1,dims,NULL);
  status = H5Tset_size(stid,20); H5Tinsert(datatype_id,DATE,HOFFSET(sdhdf_historyStruct,date),stid);
  status = H5Tset_size(stid,64); H5Tinsert(datatype_id,PROC,HOFFSET(sdhdf_historyStruct,proc_name),stid);
  status = H5Tset_size(stid,64); H5Tinsert(datatype_id,PROC_DESCR,HOFFSET(sdhdf_historyStruct,proc_descr),stid);
  status = H5Tset_size(stid,64); H5Tinsert(datatype_id,PROC_ARGS,HOFFSET(sdhdf_historyStruct,proc_args),stid);
  status = H5Tset_size(stid,64); H5Tinsert(datatype_id,PROC_HOST,HOFFSET(sdhdf_historyStruct,proc_host),stid);

  sprintf(groupName,"/%s",METADATA_GRP);
  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
    {
      group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Gclose(group_id);
    }

	sprintf(history,"%s/%s",METADATA_GRP,HISTORY);
  if (sdhdf_checkGroupExists(outFile,history) == 0) // was 0
    {
      printf("THE HISTORY TABLE ALREADY EXISTS\n");
			//printf("THE HISTORY TABLE DOES NOT EXIST\n");
      dset_id  = H5Dopen2(outFile->fileID,history,H5P_DEFAULT);
    }
  else
    {
			dset_id = H5Dcreate2(outFile->fileID,history,datatype_id,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
		}

  status  = H5Dwrite(dset_id,datatype_id,H5S_ALL,H5S_ALL,H5P_DEFAULT,outParams);
  status  = H5Dclose(dset_id);
  status  = H5Sclose(dataspace_id);
  status  = H5Tclose(datatype_id);
  status  = H5Tclose(stid);
}

void sdhdf_loadHistory(sdhdf_fileStruct *inFile)
{
  hid_t header_id,headerT,space,stid,val_tid;
  herr_t status;
  int ndims;
  hsize_t dims[2];
	char history[MAX_STRLEN];

	sprintf(history,"%s/%s",METADATA_GRP,HISTORY);
  if (sdhdf_checkGroupExists(inFile,history) == 1)
    {
      printf("WARNING: no history information\n");
      inFile->nHistory = 0;
      inFile->historyAllocatedMemory = 0;
    }
  else
    {
      header_id  = H5Dopen2(inFile->fileID,history,H5P_DEFAULT);
      headerT    = H5Dget_type(header_id);
      space      = H5Dget_space(header_id);
      ndims      = H5Sget_simple_extent_dims(space,dims,NULL);

      if (inFile->historyAllocatedMemory == 0)
			{
	  		inFile->history = (sdhdf_historyStruct *)malloc(sizeof(sdhdf_historyStruct)*MAX_HISTORY);
      	if (dims[0] > MAX_HISTORY)
				{
	  			printf("ERROR: need to increase MAX_HISTORY\n");
	  			exit(1);
				}
      	inFile->nHistory = dims[0];
      	inFile->historyAllocatedMemory = 1;
			}

      val_tid = H5Tcreate(H5T_COMPOUND,sizeof(sdhdf_historyStruct));
      stid = H5Tcopy(H5T_C_S1);

      status = H5Tset_size(stid,20); // Should set to value defined in sdhdf_v1.9.h
      H5Tinsert(val_tid,DATE,HOFFSET(sdhdf_historyStruct,date),stid);
      status = H5Tset_size(stid,64); // Should set to value defined in sdhdf_v1.9.h
      H5Tinsert(val_tid,PROC,HOFFSET(sdhdf_historyStruct,proc_name),stid);
      H5Tinsert(val_tid,PROC_DESCR,HOFFSET(sdhdf_historyStruct,proc_descr),stid);
      H5Tinsert(val_tid,PROC_ARGS,HOFFSET(sdhdf_historyStruct,proc_args),stid);
      H5Tinsert(val_tid,PROC_HOST,HOFFSET(sdhdf_historyStruct,proc_host),stid);

      status  = H5Dread(header_id,val_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,inFile->history);

			// Load attributes
			sdhdf_loadAttributes(inFile, history, &inFile->historyAttr);

      status = H5Tclose(val_tid);
      status = H5Tclose(stid);
      status = H5Dclose(header_id);
    }
}

void sdhdf_loadSchedule(sdhdf_fileStruct *inFile)
{
  hid_t header_id,headerT,space,stid,val_tid;
  herr_t status;
  int ndims;
  hsize_t dims[2];
	char sched[MAX_STRLEN];

	sprintf(sched,"%s/%s",METADATA_GRP,SCHEDULE);
  if (sdhdf_checkGroupExists(inFile,sched) == 1)
    {
      printf("WARNING: no schedule information\n");
      inFile->schedAllocatedMemory = 0;
    }
  else
    {
      header_id  = H5Dopen2(inFile->fileID,sched,H5P_DEFAULT);
      headerT    = H5Dget_type(header_id);
      space      = H5Dget_space(header_id);
      ndims      = H5Sget_simple_extent_dims(space,dims,NULL);

      if (inFile->schedAllocatedMemory == 0)
	    {
	  	 	inFile->sched = (sdhdf_scheduleStruct *)malloc(sizeof(sdhdf_scheduleStruct));
      	inFile->schedAllocatedMemory = 1;
			}

      val_tid = H5Tcreate(H5T_COMPOUND,sizeof(sdhdf_scheduleStruct));
      stid = H5Tcopy(H5T_C_S1);

      H5Tinsert(val_tid,DATE,HOFFSET(sdhdf_scheduleStruct,date),stid);
      H5Tinsert(val_tid,SCHED_HOST,HOFFSET(sdhdf_scheduleStruct,host),stid);
      H5Tinsert(val_tid,SCHED_HOSTNAME,HOFFSET(sdhdf_scheduleStruct,host_name),stid);
      H5Tinsert(val_tid,SCHED_API,HOFFSET(sdhdf_scheduleStruct,api),stid);
      H5Tinsert(val_tid,SCHED_VERSION,HOFFSET(sdhdf_scheduleStruct,version),stid);
			H5Tinsert(val_tid,SCHED_METADATA,HOFFSET(sdhdf_scheduleStruct,meta),stid);

      status  = H5Dread(header_id,val_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,inFile->sched);

			// Load attributes
			sdhdf_loadAttributes(inFile, sched, &inFile->schedAttr);

      status = H5Tclose(val_tid);
      status = H5Tclose(stid);
      status = H5Dclose(header_id);
    }
}

void sdhdf_addHistory(sdhdf_historyStruct *history,int n,char *procName,char *descr,char *args)
{
  time_t t = time(NULL);
  struct tm tm = *localtime(&t);
  char timeDate[1024];
  char host[1024];

  struct utsname unameData;
  uname(&unameData);
  strcpy(host,unameData.nodename);

  sprintf(timeDate,"%d-%02d-%02d-%02d:%02d:%02d",tm.tm_year+1900,tm.tm_mon+1,tm.tm_mday,tm.tm_hour,tm.tm_min,tm.tm_sec);
  strcpy(history[n].date,timeDate);
  strcpy(history[n].proc_name,procName);
  strcpy(history[n].proc_descr,descr);
  strcpy(history[n].proc_args,args);
  strcpy(history[n].proc_host,host);
}

void sdhdf_loadSoftware(sdhdf_fileStruct *inFile)
{
  hid_t header_id,headerT,space,stid,val_tid;
  herr_t status;
  int ndims;
  hsize_t dims[2];
	char dsetPthName[1024];

  sprintf(dsetPthName,"%s/%s",METADATA_GRP,SOFTWARE_VERSIONS);
  header_id  = H5Dopen2(inFile->fileID,dsetPthName,H5P_DEFAULT);
  headerT    = H5Dget_type(header_id);
  space      = H5Dget_space(header_id);
  ndims      = H5Sget_simple_extent_dims(space,dims,NULL);

  if (inFile->softwareAllocatedMemory == 0)
    {
      inFile->software = (sdhdf_softwareVersionsStruct *)malloc(sizeof(sdhdf_softwareVersionsStruct)*dims[0]);
      inFile->nSoftware = dims[0];
      inFile->softwareAllocatedMemory = 1;
    }

  val_tid = H5Tcreate(H5T_COMPOUND,sizeof(sdhdf_softwareVersionsStruct));
  stid = H5Tcopy(H5T_C_S1);

  status = H5Tset_size(stid,64); // Should set to value defined in sdhdf_v1.9.h
  H5Tinsert(val_tid,PROC,HOFFSET(sdhdf_softwareVersionsStruct,proc_name),stid);
  H5Tinsert(val_tid,SOFTWARE,HOFFSET(sdhdf_softwareVersionsStruct,software),stid);
  H5Tinsert(val_tid,SOFTWARE_DESCR,HOFFSET(sdhdf_softwareVersionsStruct,software_descr),stid);
  H5Tinsert(val_tid,SOFTWARE_VERSION,HOFFSET(sdhdf_softwareVersionsStruct,software_version),stid);

  status  = H5Dread(header_id,val_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,inFile->software);

	// Load attributes
	sdhdf_loadAttributes(inFile, dsetPthName, &inFile->softwareAttr);

  status = H5Tclose(val_tid);
  status = H5Tclose(stid);
  status = H5Dclose(header_id);
}

/*int sdhdf_getNattributes(sdhdf_fileStruct *inFile,char *dataName)
{
  hid_t dataset_id;
  herr_t status;
  H5O_info_t object_info;

	//printf("\nRetrieving number of attributes for: %s\n", dataName);
  dataset_id   = H5Dopen2(inFile->fileID,dataName,H5P_DEFAULT);
  #if H5_VERS_MINOR == 10
   status = H5Oget_info(dataset_id,&object_info);
  #else
   status = H5Oget_info(dataset_id,&object_info,H5O_INFO_NUM_ATTRS);
  #endif
  status = H5Dclose(dataset_id);
	//printf("Number of attributes found: %d\n",object_info.num_attrs);

  return object_info.num_attrs;
}*/

/*void sdhdf_readAttributeFromNum(sdhdf_fileStruct *inFile,char *dataName,int num,sdhdf_attributes_struct *attribute)
{
  hid_t dataset_id,attr_id,atype,atype_mem,aspace,space;
  H5T_class_t type_class;
  hid_t type,stid;
  herr_t status;
  int rank;
  hsize_t dims[5];
  int ndims;

  // Open the attribute based on its number
  attr_id      = H5Aopen_by_idx(inFile->fileID,dataName, H5_INDEX_NAME, H5_ITER_NATIVE,num,H5P_DEFAULT,H5P_DEFAULT);

  if (attr_id < 0)  // Unable to open attribute
    {
      strcpy(attribute->key,"uncertain");
      strcpy(attribute->value,"not defined");
      attribute->attributeType=0;
      return;
    }

  // Find out if this is a multi-dimensional array
  aspace       = H5Aget_space(attr_id);
  ndims        = H5Sget_simple_extent_dims(aspace,dims,NULL);
  if (ndims != 0) // 1 && ndims != 0)
    {
      hid_t dataspace_id;
      hsize_t storeSize;
      char *buf;

      H5Aget_name(attr_id,MAX_STRLEN,attribute->key);
      strcpy(attribute->value,"NOT SET");
      attribute->attributeType=0;

    }
  else // Here the number of dimensions is zero as storing a scalar
    {
      atype        = H5Aget_type(attr_id); // Retrieves a copy of the datatype for an attribute
      type_class   = H5Tget_class(atype);

      if (type_class == H5T_STRING)        // Is the attribute a string
	{
	  char string_out[MAX_STRLEN];
	  char *buffer;
	  buffer = (char *)malloc(sizeof(char)*MAX_STRLEN);

	  H5Aget_name(attr_id,MAX_STRLEN,attribute->key);

	  attribute->attributeType=0;
	  if (strcmp(attribute->key,"CLASS")==0)
	    {
	      strcpy(attribute->value,"FIX ME");
	    }
	  else
	    {
	      atype_mem = H5Tcopy(H5T_C_S1);
	      if (H5Tget_cset(atype) == H5T_CSET_UTF8) // UTF-8 STRING
		{
		  status    = H5Tset_cset(atype_mem,H5T_CSET_UTF8);
		  status    = H5Tset_size(atype_mem,H5T_VARIABLE);

		  status = H5Aread(attr_id,atype_mem,&buffer);
		  strcpy(attribute->value,buffer);
		}
	      else // ASCII STRING
		{
		  H5A_info_t attribute_info;
		  strcpy(buffer,"");

		  status = H5Aget_info(attr_id,&attribute_info);
		  status = H5Tset_cset(atype_mem,H5T_CSET_ASCII);
		  status = H5Tset_size(atype_mem,attribute_info.data_size);
		  status = H5Aread(attr_id,atype_mem,buffer);
		  buffer[attribute_info.data_size]='\0';

		  strcpy(attribute->value,buffer);

		}
	      status = H5Tclose(atype_mem);
	    }
	  free(buffer);
	}
      else if (type_class == H5T_FLOAT)
	{
	  double fval;
	  attribute->attributeType=1;
	  H5Aget_name(attr_id,MAX_STRLEN,attribute->key);
	  status  = H5Aread(attr_id, atype, &fval);
	  attribute->fvalue = fval;
	}
      else if (type_class == H5T_INTEGER)
	{
	  printf("HAVE AN INTEGER TYPE ATTRIBUTE. ERROR: DON'T KNOW WHAT TO DO\n");
	  attribute->attributeType=2;
	  exit(1);
	}
      else
	{
	  printf("CANNOT READ ATTRIBUTE\n");
	  strcpy(attribute->key,"unset");
	  strcpy(attribute->value,"unset");
	  attribute->attributeType=1;
	}
      status = H5Tclose(atype);
    }
  status = H5Sclose(aspace);
  status = H5Aclose(attr_id);

}*/

// checked OK
void sdhdf_copyBandHeaderStruct(sdhdf_bandHeaderStruct *in,sdhdf_bandHeaderStruct *out,int n)
{
  int i;
  for (i=0;i<n;i++)
    {
      strcpy(out[i].label,in[i].label);
      out[i].fc = in[i].fc;
      out[i].f0 = in[i].f0;
      out[i].f1 = in[i].f1;
      out[i].nchan = in[i].nchan;
      out[i].npol = in[i].npol;
      strcpy(out[i].pol_type,in[i].pol_type);
      out[i].dtime = in[i].dtime;
      out[i].ndump = in[i].ndump;
    }

}

void sdhdf_writeObsParams(sdhdf_fileStruct *outFile,char *bandLabel,char *beamLabel,int iband,sdhdf_obsParamsStruct *obsParams,int ndump,int type)
{
  hid_t dset_id,datatype_id,group_id;
  herr_t status;
  hid_t dataspace_id,stid;
  hsize_t dims[1];
  char name[1024];
  char groupName[1024];

	if (DEBUG == 1) printf("In sdhdf_writeObsParams...\n");

  // Need to allocate memory if I need the next line
  //  outFile->beam[ibeam].bandData[iband].nAstro_obsHeader = ndump;
  dims[0] = ndump;
  datatype_id = H5Tcreate (H5T_COMPOUND, sizeof(sdhdf_obsParamsStruct));
  stid = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(stid,64); // Should set to value defined in sdhdf_v1.9.h

  H5Tinsert(datatype_id,"ELAPSED_TIME",HOFFSET(sdhdf_obsParamsStruct,timeElapsed),H5T_NATIVE_DOUBLE);
  H5Tinsert(datatype_id,"INTEGRATION_TIME",HOFFSET(sdhdf_obsParamsStruct,dtime),H5T_NATIVE_DOUBLE);
  //H5Tinsert(datatype_id,"TIME_DB",HOFFSET(sdhdf_obsParamsStruct,timedb),stid);
  H5Tinsert(datatype_id,"MJD",HOFFSET(sdhdf_obsParamsStruct,mjd),H5T_NATIVE_DOUBLE);
  H5Tinsert(datatype_id,"UTC",HOFFSET(sdhdf_obsParamsStruct,utc),stid);
  //H5Tinsert(datatype_id,"UT_DATE",HOFFSET(sdhdf_obsParamsStruct,ut_date),stid);
  H5Tinsert(datatype_id,"LOCAL_TIME",HOFFSET(sdhdf_obsParamsStruct,local_time),stid);
  H5Tinsert(datatype_id,"RIGHT_ASCENSION",HOFFSET(sdhdf_obsParamsStruct,raStr),stid);
  H5Tinsert(datatype_id,"DECLINATION",HOFFSET(sdhdf_obsParamsStruct,decStr),stid);
	// added in v4
	H5Tinsert(datatype_id,"DRIVE_STATUS",HOFFSET(sdhdf_obsParamsStruct,fstat),H5T_NATIVE_INT);
  //H5Tinsert(datatype_id,"RA_DEG",HOFFSET(sdhdf_obsParamsStruct,raDeg),H5T_NATIVE_DOUBLE);
  //H5Tinsert(datatype_id,"DEC_DEG",HOFFSET(sdhdf_obsParamsStruct,decDeg),H5T_NATIVE_DOUBLE);
  //H5Tinsert(datatype_id,"RA_OFFSET",HOFFSET(sdhdf_obsParamsStruct,raOffset),H5T_NATIVE_DOUBLE);
  //H5Tinsert(datatype_id,"DEC_OFFSET",HOFFSET(sdhdf_obsParamsStruct,decOffset),H5T_NATIVE_DOUBLE);
  H5Tinsert(datatype_id,"GALACTIC_LONGITUDE",HOFFSET(sdhdf_obsParamsStruct,gl),H5T_NATIVE_DOUBLE);
  H5Tinsert(datatype_id,"GALACTIC_LATITUDE",HOFFSET(sdhdf_obsParamsStruct,gb),H5T_NATIVE_DOUBLE);
  H5Tinsert(datatype_id,"AZIMUTH_ANGLE",HOFFSET(sdhdf_obsParamsStruct,az),H5T_NATIVE_DOUBLE);
  H5Tinsert(datatype_id,"ZENITH_ANGLE",HOFFSET(sdhdf_obsParamsStruct,ze),H5T_NATIVE_DOUBLE);
  H5Tinsert(datatype_id,"ELEVATION_ANGLE",HOFFSET(sdhdf_obsParamsStruct,el),H5T_NATIVE_DOUBLE);
  //H5Tinsert(datatype_id,"AZ_DRIVE_RATE",HOFFSET(sdhdf_obsParamsStruct,az_drive_rate),H5T_NATIVE_DOUBLE);
  //H5Tinsert(datatype_id,"ZE_DRIVE_RATE",HOFFSET(sdhdf_obsParamsStruct,ze_drive_rate),H5T_NATIVE_DOUBLE);
	H5Tinsert(datatype_id,"HOUR_ANGLE",HOFFSET(sdhdf_obsParamsStruct,hourAngle),H5T_NATIVE_DOUBLE);
  H5Tinsert(datatype_id,"PARALLACTIC_ANGLE",HOFFSET(sdhdf_obsParamsStruct,paraAngle),H5T_NATIVE_DOUBLE);
  H5Tinsert(datatype_id,"WIND_DIRECTION",HOFFSET(sdhdf_obsParamsStruct,windDir),H5T_NATIVE_DOUBLE);
  H5Tinsert(datatype_id,"WIND_SPEED",HOFFSET(sdhdf_obsParamsStruct,windSpd),H5T_NATIVE_DOUBLE);
	// added in v4
	H5Tinsert(datatype_id,"PRESSURE",HOFFSET(sdhdf_obsParamsStruct,pressure),H5T_NATIVE_DOUBLE);
	H5Tinsert(datatype_id,"PRESSURE_MSL",HOFFSET(sdhdf_obsParamsStruct,pressureMSL),H5T_NATIVE_DOUBLE);
	H5Tinsert(datatype_id,"RELATIVE_HUMIDITY",HOFFSET(sdhdf_obsParamsStruct,relHumidity),H5T_NATIVE_DOUBLE);
	H5Tinsert(datatype_id,"TEMPERATURE",HOFFSET(sdhdf_obsParamsStruct,temp),H5T_NATIVE_DOUBLE);


  dataspace_id = H5Screate_simple(1,dims,NULL);
  // Do we need to create the groups

  sprintf(groupName,"%s",beamLabel);
  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
    {
      group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Gclose(group_id);
    }
  sprintf(groupName,"%s/%s",beamLabel,bandLabel);
  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
    {
      group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Gclose(group_id);
    }
  sprintf(groupName,"%s/%s/%s",beamLabel,bandLabel,METADATA_GRP);
  if (sdhdf_checkGroupExists(outFile,groupName) == 1)
    {
      group_id = H5Gcreate2(outFile->fileID,groupName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Gclose(group_id);
    }
  if (type==1)
    sprintf(name,"%s/%s/%s/%s",beamLabel,bandLabel,METADATA_GRP,OBS_PARAMS);
  else if (type==2)
    sprintf(name,"%s/%s/%s/%s",beamLabel,bandLabel,METADATA_GRP,CAL_OBS_PARAMS);
  dset_id = H5Dcreate2(outFile->fileID,name,datatype_id,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  status  = H5Dwrite(dset_id,datatype_id,H5S_ALL,H5S_ALL,H5P_DEFAULT,obsParams);
  status  = H5Dclose(dset_id);

  status  = H5Sclose(dataspace_id);
  status  = H5Tclose(datatype_id);

  status  = H5Tclose(stid);
}

void sdhdf_copySingleObsParams(sdhdf_fileStruct *inFile,int ibeam,int iband,int idump,sdhdf_obsParamsStruct *obsParam)
{
	if (DEBUG == 1) printf("In sdhdf_copySingleObsParams...beam %d band %d dump %d\n", ibeam, iband, idump);
  obsParam->timeElapsed = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].timeElapsed;
  obsParam->dtime = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].dtime;
  //strcpy(obsParam->timedb,inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].timedb);
  obsParam->mjd = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].mjd;
  strcpy(obsParam->utc,inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].utc);
  //strcpy(obsParam->ut_date,inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].ut_date);
  strcpy(obsParam->local_time,inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].local_time);
  strcpy(obsParam->raStr,inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].raStr);
  strcpy(obsParam->decStr,inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].decStr);
	obsParam->fstat = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].fstat;
	//strcpy(obsParam->fstat,inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].fstat);
	//obsParam->raDeg = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].raDeg;
  //obsParam->decDeg = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].decDeg;
  //obsParam->raOffset = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].raOffset;
  //obsParam->decOffset = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].decOffset;
  obsParam->gl = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].gl;
  obsParam->gb = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].gb;
  obsParam->az = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].az;
  obsParam->ze = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].ze;
  obsParam->el = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].el;
  //obsParam->az_drive_rate = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].az_drive_rate;
  //obsParam->ze_drive_rate = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].ze_drive_rate;
  obsParam->hourAngle = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].hourAngle;
  obsParam->paraAngle = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].paraAngle;
  obsParam->windDir = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].windDir;
  obsParam->windSpd = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].windSpd;
	obsParam->pressure = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].pressure;
	obsParam->pressureMSL = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].pressureMSL;
	obsParam->relHumidity = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].relHumidity;
	obsParam->temp = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].temp;
}

void sdhdf_copySingleObsParamsCal(sdhdf_fileStruct *inFile,int ibeam,int iband,int idump,sdhdf_obsParamsStruct *obsParam)
{
	if (DEBUG == 1) printf("In sdhdf_copySingleObsParamsCal...beam %d band %d dump %d\n", ibeam, iband, idump);
  obsParam->timeElapsed = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].timeElapsed;
  obsParam->dtime = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].dtime;
  //strcpy(obsParam->timedb,inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].timedb);
  obsParam->mjd = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].mjd;
  strcpy(obsParam->utc,inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].utc);
  //strcpy(obsParam->ut_date,inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].ut_date);
  strcpy(obsParam->local_time,inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].local_time);
  strcpy(obsParam->raStr,inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].raStr);
  strcpy(obsParam->decStr,inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].decStr);
	obsParam->fstat = inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].fstat;
	//strcpy(obsParam->fstat,inFile->beam[ibeam].bandData[iband].astro_obsHeader[idump].fstat);
	//obsParam->raDeg = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].raDeg;
  //obsParam->decDeg = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].decDeg;
  //obsParam->raOffset = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].raOffset;
  //obsParam->decOffset = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].decOffset;
  obsParam->gl = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].gl;
  obsParam->gb = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].gb;
  obsParam->az = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].az;
  obsParam->ze = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].ze;
  obsParam->el = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].el;
  //obsParam->az_drive_rate = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].az_drive_rate;
  //obsParam->ze_drive_rate = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].ze_drive_rate;
  obsParam->hourAngle = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].hourAngle;
  obsParam->paraAngle = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].paraAngle;
  obsParam->windDir = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].windDir;
  obsParam->windSpd = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].windSpd;
	obsParam->pressure = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].pressure;
	obsParam->pressureMSL = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].pressureMSL;
	obsParam->relHumidity = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].relHumidity;
	obsParam->temp = inFile->beam[ibeam].bandData[iband].cal_obsHeader[idump].temp;
}

void sdhdf_loadPersistentRFI(sdhdf_rfi *rfi,int *nRFI,int maxRFI,char *tel)
{
  FILE *fin;
  char fname[MAX_STRLEN];
  char runtimeDir[MAX_STRLEN];
  char str[4096];
  char dirName[1024];

  *nRFI = 0;
  sdhdf_getTelescopeDirName(tel,dirName);

  if (getenv("SDHDF_RUNTIME")==0)
    {
      printf("=======================================================================\n");
      printf("Error: require that the SDHDF_RUNTIME directory is set\n");
      printf("=======================================================================\n");
      exit(1);
    }
  strcpy(runtimeDir,getenv("SDHDF_RUNTIME"));

  sprintf(fname,"%s/observatory/%s/rfi/persistentRFI.dat",runtimeDir,dirName);
  if (!(fin = fopen(fname,"r")))
    {
      printf("ERROR: unable to open filename >%s<\n",fname);
    }
  else
    {
      while (!feof(fin))
	{
	  if (fgets(str,4096,fin)!=NULL)
	    {
	      if (str[0]!='#' && strlen(str) > 2) // Not a comment line
		{
		  sscanf(str,"%d %s %s %lf %lf %lf %lf",&(rfi[*nRFI].type),
			 rfi[*nRFI].observatory,rfi[*nRFI].receiver,&(rfi[*nRFI].f0),&(rfi[*nRFI].f1),
			 &(rfi[*nRFI].mjd0),&(rfi[*nRFI].mjd1));
		  (*nRFI)++;
		  // SHOULD CHCEK IF EXCEEDS MAXRFI -- FIX ME
		}
	    }
	}
    }
  fclose(fin);
}

void sdhdf_loadTransientRFI(sdhdf_transient_rfi *rfi,int *nRFI,int maxRFI,char *tel)
{
  FILE *fin;
  char fname[MAX_STRLEN];
  char runtimeDir[MAX_STRLEN];
  char str[4096];
  char dirName[1024];

  *nRFI = 0;


  sdhdf_getTelescopeDirName(tel,dirName);
  if (getenv("SDHDF_RUNTIME")==0)
    {
      printf("=======================================================================\n");
      printf("Error: require that the SDHDF_RUNTIME directory is set\n");
      printf("=======================================================================\n");
      exit(1);
    }
  strcpy(runtimeDir,getenv("SDHDF_RUNTIME"));

  sprintf(fname,"%s/observatory/%s/rfi/transientRFI.dat",runtimeDir,dirName);
  if (!(fin = fopen(fname,"r")))
    {
      printf("ERROR: unable to open filename >%s<\n",fname);
    }
  else
    {
      while (!feof(fin))
	{
	  if (fgets(str,4096,fin)!=NULL)
	    {
	      if (str[0]!='#' && strlen(str) > 2) // Not a comment line
		{
		  sscanf(str,"%d %s %s %lf %lf %lf %lf %lf %lf %lf",&(rfi[*nRFI].type),
			 rfi[*nRFI].observatory,rfi[*nRFI].receiver,&(rfi[*nRFI].f0),&(rfi[*nRFI].f1),
			 &(rfi[*nRFI].mjd0),&(rfi[*nRFI].mjd1),&rfi[*nRFI].threshold,&rfi[*nRFI].f2,&rfi[*nRFI].f3);
		  (*nRFI)++;
		  // SHOULD CHCEK IF EXCEEDS MAXRFI -- FIX ME
		}
	    }
	}
    }
  fclose(fin);
}

// Obtain the directory name corresponding to a telescope name
int sdhdf_getTelescopeDirName(char *tel,char *dir)
{
  char fname1[1024],fname2[1024];
  char loadLine[1024], observatoryDir[1024]="NULL";
  FILE *fin;
  char runtimeDir[1024];

  if (getenv("SDHDF_RUNTIME")==0)
    {
      printf("[sdhdf_getTelescopeDirName] ERROR: SDHDF_RUNTIME is not set\n");
      exit(1);
    }
  strcpy(runtimeDir,getenv("SDHDF_RUNTIME"));
  sprintf(fname1,"%s/observatory/observatories.list",runtimeDir);
  fin = fopen(fname1,"r");
  while (!feof(fin))
    {
      if (fgets(loadLine,1024,fin)!=NULL)
	{
	  if (loadLine[0]!='#')
	    {
	      if (strstr(loadLine,tel)!=NULL)
		{
		  sscanf(loadLine,"%s",observatoryDir);
		  break;
		}
	    }
	}
    }
  fclose(fin);
  if (strcmp(observatoryDir,"NULL")==0)
    {
      printf("ERROR: in sdhdfProc_metadata.c - unable to find observatory %s in %s\n",tel,fname1);
      exit(1);
    }
  strcpy(dir,observatoryDir);

  return 0;

}
