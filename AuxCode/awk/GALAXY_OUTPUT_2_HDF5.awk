# from a cleaned list of variables defining the GALAXY_OUPUT C struct, define hdf5 output datasiz struct definition
BEGIN{
    print "#include \"allvars.h\""
    print "#include \"proto.h\""
    print "#include \"hdf5.h\""
    print "#include \"hdf5_hl.h\""
    print " "
    print "#define NRECORDS (hsize_t) 0"
    print "#define HDF5_STRING_SIZE 2048"
    print " "
    print "/* File to set all the HDF5 table properties */ "
    print "int ifield;"
    print "int rowsize=0;"
    print "hsize_t chunk_size=CHUNK_SIZE;"
    print "int * fill_data=NULL;"
    print "hid_t file_id;"
    print "size_t * output_offsets;"
    print "hid_t * field_types;"
    print "size_t * output_sizes;"
    print "size_t output_size;"
    print " "
    print "// Write the datatype for HDF5 to write out the data"
    print " char types[]={"
    n=0
    size=0
}
{
    line=$0
    split(line,fields)

    type=fields[1]
    name[n]=fields[2]
    arraysize=1
    if(name[n]=="GALAXY_OUTPUT")
        next

    if(type == "long" || type == "double") 
    {
	dsize = 8
    } else if(type=="short")
    {
	dsize = 2
    } else 
    {
	dsize = 4
    }
    
    if(type == "float") {type= "'f"} 
    else if(type == "int"){type= "'i"}
    else if(type == "struct"){type="'s"}
    else if(type == "long" && name[n] == "long") {
    	type="'l"
    	name[n]=fields[3]
    }
    if(name[n]=="metals"){
        name[n]=fields[3]
        if(substr(name[n],0,3)=="sfh"){type="'M"}
        if(type=="'s"){type="'m"}
    }
    if(name[n]=="elements"){
        name[n]=fields[3]
        type="'e"
        if(substr(name[n],0,3)=="sfh"){type="'E"}
    } 
    ia=match(line,/\[.*\]/)
    if(ia>0) 
    {
    	arraysize=substr(line,ia+1)
    	ic=index(arraysize,"]")
    	if(ic > 0){
    	    arraysize=substr(arraysize,0,ic-1)
            if(arraysize==3){type="'3"}
            }
    }

    ia=match(name[n],/\[.*\]/)
    if(ia>0) 
	   name[n]=substr(name[n],0,ia-1)
    ia=match(name[n],";")

    if(ia>0) 
	   name[n]=substr(name[n],0,ia-1)
	if(substr(name[n],0,3)== "Mag")
	{
		type="'o"
	}
	if(substr(name[n],0,3)== "sfh" && type=="'f")
	{
		type="'s"
	}


    print type "',"
    n++
}
END{
    print "}; "
    print " "
    print "const char * field_names[]={"
    for(i=0;i<n;i++){
        print "\"" name[i] "\","
    }
    print "};"
    print " "
    print "//The number of fields in the data"
    print "int nfields=" n ";"
    print ""
    print "// Define the dimensions for the HDF5 table"
    print "hsize_t    float_dims[1]={3};"
    print "#ifdef COMPUTE_SPECPHOT_PROPERTIES"
    print "hsize_t    nmag_dims[1]={NMAG};"
    print "#endif"
    print "#ifdef STAR_FORMATION_HISTORY"
    print "hsize_t    sfh_dims[1]={SFH_NBIN};"
    print "hsize_t    sfh_3_dims[1]={3*SFH_NBIN};"
    print "#ifdef INDIVIDUAL_ELEMENTS"
    print "hsize_t    numele_dims[1]={NUM_ELEMENTS};"
    print "hsize_t    numele_sfh_dims[2]={NUM_ELEMENTS,SFH_NBIN};"
    print "hsize_t    mag_sfh_dims[2]={3,SFH_NBIN};"
    print "#endif //INDIVIDUAL_ELEMENTS"
    print "#endif //STAR_FORMATION_HISTORY"

}


