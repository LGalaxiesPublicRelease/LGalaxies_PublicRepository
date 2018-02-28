#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "allvars.h"
#include "proto.h"
#include "io_hdf5.h"

#ifdef HDF5_OUTPUT

void open_hdf5_file( int filenr){

rowsize=0;

char output_file[80];
#ifdef GALAXYTREE
sprintf(output_file,"%s/SA_galtree_%d.h5",OutputDir,filenr);
#else
sprintf(output_file,"%s/SA_output_%d.h5",OutputDir,filenr);
#endif

file_id = H5Fcreate(output_file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

output_offsets=(size_t*)calloc(nfields,sizeof(output_offsets));
output_sizes=(size_t*)calloc(nfields,sizeof(output_sizes));
field_types=(hid_t*)calloc(nfields,sizeof(field_types));

 printf("\n nfields=%d\n",nfields);
 
 //This sets the types and offsets
 for(ifield=0;ifield<nfields;ifield++){
 
   output_offsets[ifield]=rowsize;
   
   if(types[ifield]=='i'){
     field_types[ifield]=H5T_NATIVE_INT;
     output_sizes[ifield]=sizeof(int);
     rowsize+=sizeof(int);
   }
   else if(types[ifield]=='f'){
     field_types[ifield]=H5T_NATIVE_FLOAT;
    output_sizes[ifield]=sizeof(float);
     rowsize+=sizeof(float);
   }
   else if(types[ifield]=='l'){
     field_types[ifield]=H5T_NATIVE_LLONG;
     output_sizes[ifield]=sizeof(long long);
     rowsize+=sizeof(long long);
   }
   else if(types[ifield]=='3'){
     field_types[ifield]=H5Tarray_create(H5T_NATIVE_FLOAT,1,float_dims);
     output_sizes[ifield]=3*sizeof(float);
     rowsize+=3*sizeof(float);
   }
#ifdef COMPUTE_SPECPHOT_PROPERTIES
   else if(types[ifield]=='o'){
     field_types[ifield]=H5Tarray_create(H5T_NATIVE_FLOAT,1,nmag_dims);
     output_sizes[ifield]=NMAG*sizeof(float);
     rowsize+=NMAG*sizeof(float);
   }
#endif
#ifdef STAR_FORMATION_HISTORY
   else if(types[ifield]=='s'){
     field_types[ifield]=H5Tarray_create(H5T_NATIVE_FLOAT,1,sfh_dims);
     output_sizes[ifield]=SFH_NBIN*sizeof(float);
     rowsize+=SFH_NBIN*sizeof(float);
   }
#ifdef INDIVIDUAL_ELEMENTS
    else if(types[ifield]=='m'){
     field_types[ifield]=H5Tarray_create(H5T_NATIVE_FLOAT,1,float_dims);
     output_sizes[ifield]=sizeof(struct metals); 
     rowsize+=sizeof(struct metals);
   }
   else if(types[ifield]=='M'){
     field_types[ifield]=H5Tarray_create(H5T_NATIVE_FLOAT,2,mag_sfh_dims);
     output_sizes[ifield]=sizeof(struct metals)*SFH_NBIN; 
     rowsize+=sizeof(struct metals)*SFH_NBIN;
   }
   else if(types[ifield]=='e'){
     field_types[ifield]=H5Tarray_create(H5T_NATIVE_FLOAT,1,numele_dims);
     output_sizes[ifield]=NUM_ELEMENTS*sizeof(float);
     rowsize+=NUM_ELEMENTS*sizeof(float);
   }
   else if(types[ifield]=='E'){
     field_types[ifield]=H5Tarray_create(H5T_NATIVE_FLOAT,2,numele_sfh_dims);
     output_sizes[ifield]=SFH_NBIN*NUM_ELEMENTS*sizeof(float);
     rowsize+=SFH_NBIN*NUM_ELEMENTS*sizeof(float);
     }
#endif //INDIVIDUAL_ELEMENTS
#ifdef DETAILED_DUST 			//Scott 15/11/2015
   else if(types[ifield]=='d'){
     field_types[ifield]=H5Tarray_create(H5T_NATIVE_FLOAT,1,dustrate_dims);
     output_sizes[ifield]=sizeof(struct DustRates);
     rowsize+=sizeof(struct DustRates);
   }
   else if(types[ifield]=='D'){
     field_types[ifield]=H5Tarray_create(H5T_NATIVE_FLOAT,2,dustrate_sfh_dims);
     output_sizes[ifield]=SFH_NBIN*sizeof(struct DustRates);
     rowsize+=SFH_NBIN*sizeof(struct DustRates);
   }
#endif
#endif // STAR_FORMATION_HISTORY

    else {
     char err[100];
     sprintf(err,"***open_hdf5_file: undefined datatype %c \n", types[ifield]);
     terminate(err);  
    }
 }
}

void create_hdf5_table(int n){

    printf("\n Entering create_hdf5_table");

  char table_name[20];

  output_size=sizeof(struct GALAXY_OUTPUT);

  struct GALAXY_OUTPUT galaxy_output;

#ifdef GALAXYTREE
  sprintf(table_name,"GalTree"); 
#else
  sprintf(table_name,"%d",ListOutputSnaps[n]);
  printf("\n Making table %s\n",table_name);
#endif
#ifdef DEBUG_HDF5
  printf("table_name=%s\n",table_name);
  printf("file_id=%d\n",(int)file_id);
  printf("nfields=%d\n",nfields);
  printf("NRECORDS=%d\n",(int)NRECORDS);
  printf("output_size=%d\n",(int)output_size);
  int ifield;
  for (ifield=0;ifield<nfields;ifield++) {
      printf("field_names[%d]=%s\n",ifield,field_names[ifield]);
      printf("output_offsets[%d]=%d\n",ifield,(int)output_offsets[ifield]);
      printf("field_types[%d]=%d\n",ifield,(int)field_types[ifield]);
  }
  //exit(1);
#endif
  H5TBmake_table(table_name, file_id, table_name,nfields,NRECORDS,
           output_size,field_names, output_offsets, field_types,
           chunk_size, fill_data, COMPRESS, &galaxy_output  );

    printf("Leaving create_hdf5_table\n");

}

void hdf5_append_data(int n, struct GALAXY_OUTPUT * galaxy_output,int nrecords_app){

  char table_name[20];
  // Create the tablename from the snapshots
#ifdef GALAXYTREE
  sprintf(table_name,"GalTree");
#else
  sprintf(table_name,"%d",ListOutputSnaps[n]);
#endif

  H5TBappend_records(file_id,table_name,nrecords_app, output_size, output_offsets, output_sizes,galaxy_output);

}

void hdf5_close(){
  //Write tables from the input files, use: write_input_table(directory,filename)
#ifdef DEBUG_HDF5
  printf("Closing HDF5_file...\n");
#endif
  write_input_table("My_Makefile_options");
  write_input_table(inputFile);
#ifdef DEBUG_HDF5
  printf("wrote_input_table inputFile\n");
#endif
  write_prop_table();
#ifdef DEBUG_HDF5
  printf("wrote_prop_table\n");
#endif

  printf("\nClosing hdf5 file \n");
  H5Fclose( file_id ); 
#ifdef DEBUG_HDF5
  printf("HDF5 file closed\n");
#endif

  free(output_offsets);
  free(field_types);
  free(output_sizes);
}


void write_prop_table(void ){

  //Setup a Prop table to give a brief description about the data

  struct PropTable{
    char Name[HDF5_STRING_SIZE];
    char Units[HDF5_STRING_SIZE];
    char Description[HDF5_STRING_SIZE];
  }proptable;

  size_t proptable_output_size = sizeof(struct PropTable);

  hsize_t PropTable_nfields=3;
  const char* proptable_field_names[3]={"Field","Units","Description"};

  // Sets a string typ econtaining CSIZE characters
  hid_t string_type;
  size_t proptable_output_offsets[3]={HOFFSET(struct PropTable,Name),
				      HOFFSET(struct PropTable,Units),
				      HOFFSET(struct PropTable,Description)};
  hid_t   proptable_field_types[3];
  size_t proptable_output_sizes[3]={sizeof(proptable.Name),
				    sizeof(proptable.Units),
				    sizeof(proptable.Description)};

  struct PropTable prop_table[127];

  // struct PropTable *proptable;
  string_type=H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type,HDF5_STRING_SIZE);
  proptable_field_types[0]=string_type;
  proptable_field_types[1]=string_type;
  proptable_field_types[2]=string_type;
 
  H5TBmake_table("Prop Table", file_id, "Prop Table",PropTable_nfields,NRECORDS,
		 proptable_output_size ,proptable_field_names, proptable_output_offsets, proptable_field_types,
		 chunk_size, fill_data, COMPRESS, &prop_table  );

  FILE *fp; // The file pointer 
  char line[HDF5_STRING_SIZE];
  char *token;
  const char d[2]=",";
  char filename[30]="input/hdf5_field_props.txt";

    fp=fopen(filename,"r");

  if(fp==NULL){
    printf("could not open the file %s \n", filename);
    exit(1);
  }
  // Put data into the struct and append to the table

#ifdef DEBUG_HDF5
  int icount=0;
#endif
  while(fgets(line,sizeof(line),fp)){
    token=strtok(line,d);
    if (token!=NULL) 
	strcpy(prop_table->Name,token);
    else
	strcpy(prop_table->Name," ");
    token=strtok(NULL,d);
    if (token!=NULL) 
	strcpy(prop_table->Units,token);
    else
	strcpy(prop_table->Units," ");
    token=strtok(NULL,d);
    if (token!=NULL) 
	strcpy(prop_table->Description,token);
    else
	strcpy(prop_table->Description," ");
#ifdef DEBUG_HDF5
    printf("CP3.%d\n finished prop_table assignment",icount);
    printf("prop_table->Name = %s\n",prop_table->Name);
    printf("prop_table->Description = %s\n",prop_table->Description);
    printf("prop_table->Units = %s\n",prop_table->Units);
#endif
    H5TBappend_records(file_id,"Prop Table",1, proptable_output_size, proptable_output_offsets, proptable_output_sizes,&prop_table);
#ifdef DEBUG_HDF5
    icount++;
#endif
  } 

  H5Tclose( string_type ); 

}

void write_input_table(char *filename){

  //Setup a Prop table to give a brief description about the data

  struct InputTable{
    char Input[2048];
  }inputtable;

  char table_name[HDF5_STRING_SIZE];

  // Table filename after last slash (otherwise HDF5 crashes)
  char* delim;
  strcpy(table_name,filename);
  delim = strchr(filename,'/');
  while (delim!=NULL) {
      sprintf(table_name,"%s",++delim);
      delim = strchr(table_name,'/');
  }

  size_t inputtable_output_size = sizeof(struct InputTable);

  hsize_t InputTable_nfields=1;
  const char* inputtable_field_names[1]={filename};

  // Sets a string typ econtaining CSIZE characters
  hid_t string_type;
  size_t inputtable_output_offsets[1]={HOFFSET(struct InputTable,Input)};
  hid_t   inputtable_field_types[1];
  size_t inputtable_output_sizes[1]={sizeof(inputtable.Input)};

  struct InputTable input_table[1];

  // struct InputTable *inputtable;
  string_type=H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type,2048);
  inputtable_field_types[0]=string_type;

 

   H5TBmake_table(table_name, file_id, table_name,InputTable_nfields,NRECORDS,inputtable_output_size ,inputtable_field_names, inputtable_output_offsets, inputtable_field_types,chunk_size, fill_data, COMPRESS, &input_table  );

  FILE *fp; // The file pointer 
  char line[2048];

  fp=fopen(filename,"r");
  if(fp==NULL){
    printf("could not open the file %s \n", filename);
    exit(1);
  }
  // Put data into the struct and append to the table

 
  while(fgets(line,sizeof(line),fp)){
    strcpy(input_table->Input,line);
    H5TBappend_records(file_id,table_name,1, inputtable_output_size, inputtable_output_offsets, inputtable_output_sizes,&input_table);
  } 

  H5Tclose( string_type ); 
   

}


#endif //HDF5


