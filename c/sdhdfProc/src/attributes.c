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
#include <ctype.h>

char* sdhdf_getAttribute(sdhdf_attributes_struct *attrStruct, char *attr_name)
{
	int i;
	char* def = malloc(MAX_STRLEN * sizeof(char));

	for (i=0;i<MAX_ATTRIBUTES;i++)
		{
			if (strcmp(attrStruct[i].name,attr_name)==0)
				{
					printf("\nFound matching attribute: %s\n", attrStruct[i].name);
					printf("DESCRIPTION: %s\n", attrStruct[i].key);
					printf("UNIT:        %s\n", attrStruct[i].value);
					printf("DEFAULT:     %s\n", attrStruct[i].def);
					strcpy(def, attrStruct[i].def); break;
				}
		}
	return def;

}

void sdhdf_readAttributes(int num_attrs, sdhdf_attributes_struct *attrStruct)
{
	int i;
	for (i=0;i<MAX_ATTRIBUTES;i++)
		{
			if (strcmp(&attrStruct[i].name,"")==0) break;
			printf("NAME:        %s\n", &attrStruct[i].name);
			printf("DESCRIPTION: %s\n", &attrStruct[i].key);
			printf("UNIT:        %s\n", &attrStruct[i].value);
			printf("DEFAULT:     %s\n\n", &attrStruct[i].def);
		}

}

void sdhdf_copyAttributes(sdhdf_attributes_struct *in,int n_in,sdhdf_attributes_struct *out,int *n_out)
{
  int i;
  *n_out = n_in;
  for (i=0;i<n_in;i++)
    {
      printf("In sdhdf_copyAttributes: %s %s\n",in[i].key,in[i].value);
      strcpy(out[i].key,in[i].key);
      out[i].attributeType = in[i].attributeType;
      if (in[i].attributeType==0)
	strcpy(out[i].value,in[i].value);
      else if (in[i].attributeType==1)
	out[i].fvalue = in[i].fvalue;
      else if (in[i].attributeType==2)
	out[i].ivalue = in[i].ivalue;
    }
}

void sdhdf_copyAttributes2(sdhdf_attributes_struct *in, sdhdf_attributes_struct *out)
{
	int i;
	for (i=0;i<MAX_ATTRIBUTES;i++)
    {
			printf("In sdhdf_copyAttributes...\n");
			// hmm inFile->metaAttr appears to be empty!!
			if (strcmp(in[i].name,"")==0) break;
			strcpy(out[i].name,in[i].name);
			strcpy(out[i].key,in[i].key);
      strcpy(out[i].value,in[i].value);
			strcpy(out[i].def,in[i].def);
    }

}

void sdhdf_loadGroupAttributes(sdhdf_fileStruct *inFile, char *dataName, sdhdf_attributes_struct *attrStruct)
{
    hid_t  attr, attr_id, s1_tid, dataset_id, group_id, dtype;
    herr_t status;
		//char name[MAX_STRLEN];
		char attr_name[MAX_STRLEN];
		int i;
		int num_attrs;

		printf("dataName: %s\n", dataName);
		printf("to here 1\n");

		if (DEBUG == 1) printf("Loading attributes for %s...\n", dataName);
		group_id = H5Gopen(inFile->fileID, dataName, H5P_DEFAULT);
		num_attrs = H5Aget_num_attrs(group_id);

		for (i=0;i<num_attrs;i++)
			{
				printf("start loop %d\n", i);

				attr_id = H5Aopen_idx(group_id, i);
				H5Aget_name(attr_id,MAX_STRLEN,attr_name);
				printf("Name: %s\n", attr_name);

				attr = H5Aopen(group_id, attr_name, H5P_DEFAULT);
				dtype = H5Tcopy(H5T_C_S1);
				status = H5Tset_size(dtype, MAX_STRLEN);

				s1_tid = H5Tcreate(H5T_COMPOUND, sizeof(sdhdf_attributes_struct));

				H5Tinsert(s1_tid, "name", HOFFSET(sdhdf_attributes_struct, name), dtype);
				H5Tinsert(s1_tid, "description", HOFFSET(sdhdf_attributes_struct, key), dtype);
				H5Tinsert(s1_tid, "unit", HOFFSET(sdhdf_attributes_struct, value), dtype);
				H5Tinsert(s1_tid, "default", HOFFSET(sdhdf_attributes_struct, def), dtype);

				//printf("size of s1_tid in bytes : %d\n", sizeof(s1_tid));

				// BREAKS ON THIS LINE on the second loop with abort trap 6
				//status = H5Aread(attr, s1_tid, &attrStruct[i]);
				//

				strcpy(&attrStruct[i].name, attr_name);

				if (DEBUG == 1) printf("NAME:        %s\n", &attrStruct[i].name);
				if (DEBUG == 1) printf("DESCRIPTION: %s\n", &attrStruct[i].key);
				if (DEBUG == 1) printf("UNIT:   		 %s\n", &attrStruct[i].value);
				if (DEBUG == 1) printf("DEFAULT: 	   %s\n", &attrStruct[i].def);

				printf("to here 2\n");

				H5Tclose(s1_tid);
				H5Aclose(attr);

			}
			H5Gclose(group_id);
			printf("to here 3\n");
			printf("size of struct in bytes : %d\n", sizeof(attrStruct));
}

void sdhdf_loadAttributes(sdhdf_fileStruct *inFile, char *dataName, sdhdf_attributes_struct *attrStruct)
{
    hid_t  attr, attr_id, s1_tid, dataset_id, group_id, space, dtype1, space_type, type, type_class, dset, rootGroup_id;
    herr_t status, hErrVal;
		char name[MAX_STRLEN];
		char attr_name[MAX_STRLEN];
		char space_name[MAX_STRLEN];
		H5G_stat_t objStatInfo;
		int i;
		int num_attrs;

		//printf("dataName: %s\n", dataName);

		if (DEBUG == 1) printf("Loading attributes for %s...\n", dataName);
		rootGroup_id = H5Gopen(inFile->fileID, "/", H5P_DEFAULT);
		num_attrs = H5Aget_num_attrs(rootGroup_id);
		hErrVal = H5Gget_objinfo(rootGroup_id, dataName, 0, &objStatInfo);

		if (strcmp(dataName, "/") != 0)
		  {
		  if (objStatInfo.type == H5G_DATASET)
			  {
			  // dataset
			  if (DEBUG == 1) printf("Object %s is an HDF dataset\n", dataName);
			  dataset_id = H5Dopen2(inFile->fileID, dataName, H5P_DEFAULT);
				num_attrs = H5Aget_num_attrs(dataset_id);
		    }
		  else
		    {
			  // group
			  if (DEBUG == 1) printf("Object %s is an HDF group\n", dataName);
			  dataset_id = H5Dopen2(rootGroup_id, dataName, H5P_DEFAULT);
				num_attrs = H5Aget_num_attrs(dataset_id);
		    }
		  }
		/*else
		  {
			  num_attrs = H5Aget_num_attrs(dataset_id);
		  }*/

		printf("Num attrs: %d\n", num_attrs);
		//exit(1);

		for (i=0;i<MAX_ATTRIBUTES;i++)
			{
				printf("start loop %d\n", i);
				printf("dataName: %s\n", dataName); //dataName "/" gets reset to " " after first loop

				if (strcmp(dataName, "/") == 0)
				{
					//if (i == 2) break;
					attr_id = H5Aopen_idx(rootGroup_id, i);
					//H5Aget_name(attr_id,MAX_STRLEN,attr_name);
					//attr = H5Aopen(rootGroup_id, attr_name, H5P_DEFAULT);
					//attr_id = H5Aopen_by_idx(inFile->fileID, dataName, H5_INDEX_NAME, H5_ITER_NATIVE, i, H5P_DEFAULT, H5P_DEFAULT);
				}
				else
				{
				  attr_id = H5Aopen_by_idx(inFile->fileID, dataName, H5_INDEX_NAME, H5_ITER_NATIVE, i, H5P_DEFAULT, H5P_DEFAULT);
					//H5Aget_name(attr_id,MAX_STRLEN,attr_name);
					//attr = H5Aopen(dataset_id, attr_name, H5P_DEFAULT);
				}
				//else
			//	{
				  //attr_id = H5Aopen_by_idx(inFile->fileID, dataName, H5_INDEX_NAME, H5_ITER_NATIVE, i, H5P_DEFAULT, H5P_DEFAULT);
				//}
				if (attr_id == -1)
				{
					break;
				}
				printf("to here 1\n");
				H5Aget_name(attr_id,MAX_STRLEN,attr_name);

				// get attribute type
				space = H5Aget_space(attr_id);
				space_type = H5Sget_simple_extent_type(space);
				type = H5Aget_type(attr_id);
				type_class   = H5Tget_class(type);
			  if (space_type == H5S_SCALAR)
				  strcpy(space_name, "H5S_SCALAR");
			  else if (space_type == H5S_SIMPLE)
				  strcpy(space_name, "H5S_SIMPLE");
			  else if (space_type == H5S_NULL)
				  strcpy(space_name, "H5S_NULL");
			  else
				  strcpy(space_name, "UNKNOWN");

				// type_class = 6 COMPOUND
				// type_class = 3, space = 119, space_type = 0, space_name = H5S_SCALAR, attr_name = CLASS
				// type_class = 6, space = 120, space_type = 1, space_name = H5S_SIMPLE, attr_name = REFERENCE_LIST
				// type_class = 3, space = 121, space_type = 1, space_name = H5S_SIMPLE, attr_name = DIMENSION_LABELS
				// type_class = 9, space = 141, space_type = 1, space_name = H5S_SIMPLE, attr_name = DIMENSION_LIST

				if (DEBUG == 1) printf("attr_id = %d, type = %d, type_class = %d, space = %d, space_type = %d, space_name = %s, attr_name = %s\n\n",
								(int)attr_id,
								(int)type, (int)type_class,
								(int)space, (int)space_type,
								space_name, attr_name);

				if (type_class == 3) {
					// H5T_STRING datatype
					if (space_type == 0) {
						// H5S_SCALAR
					}
					else if (space_type == 1) {
						// H5S_SIMPLE
					}
					else {
						if (DEBUG == 1) printf("UNKNOWN dataspace: %d", space_type);
					}
				}
				else if (type_class == 6) {
					// H5T_COMPOUND datatype
					if (space_type == 0) {
						// H5S_SCALAR
					}
					else if (space_type == 1) {
						// H5S_SIMPLE
						// dataset
						/*if (objStatInfo.type == H5G_DATASET) {
							dataset = H5Dopen2(inFile->fileID, dataName, H5P_DEFAULT);
						}
						else
						{
						  // group
							// TODO this is not getting group attributes!
							// think we need to open this once not in a loop!
							dataset = H5Gopen(inFile->fileID, dataName, H5P_DEFAULT);
						}*/
							attr = H5Aopen(dataset_id, attr_name, H5P_DEFAULT);
	    		  	dtype1 = H5Tcopy(H5T_C_S1);
	    		  	status = H5Tset_size(dtype1, MAX_STRLEN);

					  	s1_tid = H5Tcreate(H5T_COMPOUND, sizeof(sdhdf_attributes_struct));

					  	H5Tinsert(s1_tid, "name", HOFFSET(sdhdf_attributes_struct, name), dtype1);
	    		  	H5Tinsert(s1_tid, "description", HOFFSET(sdhdf_attributes_struct, key), dtype1);
	    		  	H5Tinsert(s1_tid, "unit", HOFFSET(sdhdf_attributes_struct, value), dtype1);
	    		  	H5Tinsert(s1_tid, "default", HOFFSET(sdhdf_attributes_struct, def), dtype1);

					  	status = H5Aread(attr, s1_tid, &attrStruct[i]);

					  	strcpy(&attrStruct[i].name, attr_name);

					  	if (DEBUG == 1) printf("NAME:        %s\n", &attrStruct[i].name);
					  	if (DEBUG == 1) printf("DESCRIPTION: %s\n", &attrStruct[i].key);
					  	if (DEBUG == 1) printf("UNIT:   		 %s\n", &attrStruct[i].value);
					  	if (DEBUG == 1) printf("DEFAULT: 	   %s\n", &attrStruct[i].def);

							printf("to here 2\n");
	    		  	H5Tclose(s1_tid);
							/*if (objStatInfo.type == H5G_DATASET) {
	    		  		H5Dclose(dataset);
							}
							else
							{
								H5Gclose(dataset);
							}*/
							printf("to here 3\n");
							H5Aclose(attr_id);
							printf("to here 4\n");

					}
					else {
						if (DEBUG == 1) printf("UNKNOWN dataspace: %d", space_type);
					}
				}
				else if (type_class == 9) {
					// ??
					if (space_type == 0) {
						// H5S_SCALAR
					}
					else if (space_type == 1) {
						// H5S_SIMPLE
					}
					else {
						if (DEBUG == 1) printf("UNKNOWN dataspace: %d", space_type);
					}
				}
				else {
					if (DEBUG == 1) printf("UNKNOWN datatype: %d", type_class);
				}
			printf("end loop %d\n", i);
		}
		// need to FIX close!!
		/*if (objStatInfo.type == H5G_DATASET) {
			H5Dclose(dataset);
		}
		else
		{
			H5Gclose(dataset);
		}*/
		H5Gclose(rootGroup_id);
		printf("size of struct in bytes : %d\n", sizeof(attrStruct));
}

void sdhdf_loadDataFreqAttributes(sdhdf_fileStruct *inFile,int beam,int band,sdhdf_attributes_struct *dataAttributes,int *nDataAttributes,
				  sdhdf_attributes_struct *freqAttributes,int *nFreqAttributes)
{
  int j;
  char dataName[MAX_STRLEN];
  char beamLabel[MAX_STRLEN];

  strcpy(beamLabel,inFile->beamHeader[beam].label);
  sprintf(dataName,"%s/%s/%s/%s",beamLabel,inFile->beam[beam].bandHeader[band].label,ASTRONOMY_DATA_GRP,FREQUENCY);
  //*nFreqAttributes = sdhdf_getNattributes(inFile,dataName);
  //      printf("Number of attributes = %d FREQUENCY SECTION\n",inFile->beam[beam].bandData[band].nAstro_obsHeaderAttributes_freq);
  sdhdf_loadAttributes(inFile, dataName, &freqAttributes);
	//for (j=0;j<*nFreqAttributes;j++)
  //    sdhdf_readAttributeFromNum(inFile,dataName,j,&freqAttributes[j]);

  sprintf(dataName,"%s/%s/%s/%s",beamLabel,inFile->beam[beam].bandHeader[band].label,ASTRONOMY_DATA_GRP,DATA);

  // Read attributes
  //*nDataAttributes = sdhdf_getNattributes(inFile,dataName);
	sdhdf_loadAttributes(inFile, dataName, &dataAttributes);
  //for (j=0;j< *nDataAttributes; j++)
    //sdhdf_readAttributeFromNum(inFile,dataName,j,&dataAttributes[j]);

}

void sdhdf_writeAttribute(sdhdf_fileStruct *outFile,char *dataName,sdhdf_attributes_struct *attribute) //,char *attrName,char *result)
{
  hid_t memtype,dataset_id,attr,space,stid;
  herr_t status;
  hsize_t dims[1] = {1}; // Has only 1 string
  char *result;

  result = attribute->value;
  //  printf("In writeAttribute with dataName = >%s<, attrName = >%s< and result = >%s<\n",dataName,attribute->key,attribute->value);
  //  printf("AttributeType = %d\n",attribute->attributeType);
  // Python strings are written as variable length strings


  dataset_id   = H5Dopen2(outFile->fileID,dataName,H5P_DEFAULT);

  // This produces a multi dimensional attribute - which is not the same as the python code
  /*
  stid    = H5Tcopy(H5T_C_S1);
  status  = H5Tset_size(stid,H5T_VARIABLE);
  space   = H5Screate_simple(1,dims,NULL);
  attr    = H5Acreate(dataset_id,attrName,stid,space,H5P_DEFAULT,H5P_DEFAULT);
  printf("Trying to write >%s<\n",result);
  status = H5Awrite(attr,stid,&result);
  printf("The status after writing is: %d\n",status);
  H5Dclose(dataset_id);
  H5Aclose(attr);
  H5Tclose(memtype);
  */

  // This produces an attribute which is a string
  if (attribute->attributeType==0)
    {
      stid    = H5Tcopy(H5T_C_S1);
      //  status  = H5Tset_size(stid,H5T_VARIABLE);
      status  = H5Tset_size(stid,strlen(result));
      space   = H5Screate(H5S_SCALAR);
      attr    = H5Acreate(dataset_id,attribute->key,stid,space,H5P_DEFAULT,H5P_DEFAULT);
      //      printf("writeAttribute: Trying to write >%s<\n",result);
      status  = H5Awrite(attr,stid,result);
      //      printf("writeAttribute: The status after writing is: %d\n",status);
      H5Dclose(dataset_id);
      H5Tclose(memtype);

      H5Aclose(attr);
    }
  else if (attribute->attributeType==1)
    {
      float val = 1.0;
      int ival= 1;

      //      stid    = H5Tcopy(H5T_C_S1);
      //  status  = H5Tset_size(stid,H5T_VARIABLE);
      //      status  = H5Tset_size(stid,strlen(result));
      //      space   = H5Screate(H5S_SCALAR);

      space = H5Screate(H5S_SCALAR);

      //      attr    = H5Acreate(dataset_id,attribute->key,stid,space,H5P_DEFAULT,H5P_DEFAULT);
      attr    = H5Acreate(dataset_id,attribute->key,H5T_NATIVE_FLOAT,space,H5P_DEFAULT,H5P_DEFAULT);

      //      printf("writeAttribute: Trying to write >%s<\n",result);
      status  = H5Awrite(attr,H5T_NATIVE_FLOAT,&(attribute->fvalue));
      //      status  = H5Awrite(attr,stid,result);
      //      printf("writeAttribute: The status after writing is: %d\n",status);

      H5Dclose(dataset_id);
      H5Tclose(memtype);

      H5Aclose(attr);

      /*
      printf("GEORGE: WRITING A FLOATING VALUE INTO %s\n",attribute->key);
      space   = H5Screate(H5S_SCALAR);
      printf("GEORGE A %d\n",space);
      //      attr    = H5Acreate(dataset_id,attribute->key,H5T_NATIVE_FLOAT,space,H5P_DEFAULT,H5P_DEFAULT);
      //      attr    = H5Acreate(dataset_id,attribute->key,H5T_NATIVE_FLOAT,space,H5P_DEFAULT,H5P_DEFAULT);
      //      attr    = H5Acreate(dataset_id,attribute->key,H5T_FLOAT,space,H5P_DEFAULT,H5P_DEFAULT);
      attr    = H5Acreate(dataset_id,attribute->key,H5T_NATIVE_INT,space,H5P_DEFAULT,H5P_DEFAULT);
      printf("GEORGE B %f\n",attribute->fvalue);
      //      status  = H5Awrite(attr,H5T_NATIVE_FLOAT,&val);
      //      status  = H5Awrite(attr,H5T_FLOAT,&val);
      status  = H5Awrite(attr,H5T_NATIVE_INT,&ival);
      printf("GEORGE C >%d<\n",status);
      H5Aclose(attr);
      printf("COMPLETE WRITE WITH STATUS: %s\n",status);
      */
    }

  //  memtype = H5Tcopy(H5T_C_S1);
  //  status  = H5Tset_size(memtype,MAX_STRLEN);
  //  space   = H5Screate(H5S_SCALAR);

  //  status = H5Tset_size(memtype,H5T_VARIABLE);


  //
  //  space = H5Screate_simple(1,dims,NULL);
  //  space = H5Screate(H5S_SCALAR);


  //  attr = H5Acreate(dataset_id,attrName,memtype,space,H5P_DEFAULT,H5P_DEFAULT);
  // attr = H5Acreate(dataset_id,attrName,memtype,space,H5P_DEFAULT,H5P_DEFAULT);

  //  status = H5Awrite(attr,space,result);


}
