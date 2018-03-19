# From a cleaned list of variables defining the GALAXY_OUPUT C struct, define hdf5 output datasiz struct definition
# Here is a lookup table of types:
# f - float
# i - integer
# l - long long
# d - struct DustRates
# e - struct elements (no longer used)
# m - struct metals

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

# Set length of variable in bytes
    if(type == "long" || type == "double") dsize = 8
    else if(type=="short") dsize = 2
    else dsize = 4

# Short names for variable types.  These will get over-written below for more conplex types
    if(type == "float") {type= "f"} 
    else if(type == "int"){type= "i"}
    else if(type == "long" && name[n] == "long") {
    	type="'l"
    	name[n]=fields[3]
    }

# Flags that we will need later to assign more complex data types
# Is this a metal struct?
    flagMetals=0
    if(name[n]=="metals"){
        name[n]=fields[3]
	flagMetals=1
	type="m"
    }
# Is this a DustRate struct?
    flagDustRates=0
    if(name[n]=="DustRates"){
	name[n]=fields[3]
	flagDustRates=1
	type="d"
    }
# Flag to see if this has an array of size 3
    flag3=match(line,/\[3\]/)
# Flag to see if this is a magnitude
# Note: at this point NMAG has been substituted, but not any of the other parameterised arrays sizes - why not?
    flagMag=match(name[n],/Mag/)
# Flag to decide whether or not this entry contains Rings
    flagRings=match(name[n],/\[RNUM\]/)
# Flag to decide whether or not this entry contains star formation history
    flagSFH=match(name[n],/\[SFH_NBIN\]/)

    # if(name[n]=="metals"){
    #     name[n]=fields[3]
    #     if(substr(name[n],0,3)=="sfh"){type="M"}
    #     if(type=="'s"){type="m"}
    # }
    # if(name[n]=="elements"){
    #     name[n]=fields[3]
    #     type="e"
    #     if(substr(name[n],0,3)=="sfh"){type="E"}
    # } 
    # ia=match(line,/\[.*\]/)
    # if(ia>0) 
    # {
    # 	arraysize=substr(line,ia+1)
    # 	ic=index(arraysize,"]")
    # 	if(ic > 0){
    # 	    arraysize=substr(arraysize,0,ic-1)
    #         if(arraysize==3){type="3"}
    # 	}
    # }

# Now set type of variable.  These type names are used in io_hdf5.c to set data table sizes
# If this is an array of size 3 then set type equal to 3.  Assumes that can never be part of larger structure/array
    if(flag3>0) {
	if( (flagMag>0) || (flagRings>0) || (flagSFH>0) ) {
	    print "variable type not yet implemented"
	    type="X"
	}
	else type="3"
    }
# If this is a magnitude then set type to 'o'.  Assumes that can never be part of a larger structure/array (i.e. no rings)
    if(flagMag>0) {
	if( (flag3>0) || (flagRings>0) || (flagSFH>0) ) {
	    print "variable type not yet implemented"
	    type="X"
	}
	else type="o"
    }
# If this is a SFH array then 
    if(flagSFH>0) {
        if( (flag3>0) || (flagMag>0) || (flagRings>0) ) {
	    print "variable type not yet implemented"
	    type="X"
	}
	else if (flagMetals>0) type="M"
	else if (flagElements>0) type="E"
	else if (flagDustRates>0) type="D"
	else type="s"
    }
# If this is a rings array then
    if(flagRings>0) {
        if( (flag3>0) || (flagMag>0) || (flagSFH>0) ) {
	    print "variable type not yet implemented"
	    type="X"
	}
	else if (flagMetals>0) type="M"
	else if (flagElements>0) type="E"
	else if (flagDustRates>0) type="D"
	else type="r"
    }

# Print type to h_iodefs.h
    print "'" type "',"

# Extract the name.
# Strip an array expression off the end of name
    ia=match(name[n],/\[.*\]/)
    if(ia>0) 
	   name[n]=substr(name[n],0,ia-1)
# This strips off any trailing semi-colon
    ia=match(name[n],";")
    if(ia>0) 
	name[n]=substr(name[n],0,ia-1)

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
