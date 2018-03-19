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
    print "hid_t field_type;"
    print "hid_t * field_types;"
    print "size_t output_size;"
    print "size_t * output_sizes;"
    print " "
    print "// Define the dimensions for the HDF5 table"
    print "// Should be as large as the number of possible array dimensions+1"
    print "hsize_t dims[6];"
    print " "
    n=0
    size=0
}
{
    line=$0
    split(line,fields)

    type[n]=fields[1]
    name[n]=fields[2]
    arraysize=1
    if(name[n]=="GALAXY_OUTPUT")
        next

# Short names for variable types.
    if(type[n] == "float") {type[n]= "f"} 
    else if(type[n] == "int"){type[n]= "i"}
    else if(type[n] == "long" && name[n] == "long") {
    	type[n]="'l"
    	name[n]=fields[3]
    }
# Is this a metal struct?
    if(name[n]=="metals"){
        name[n]=fields[3]
	type[n]="m"
    }
# Is this a DustRate struct?
    if(name[n]=="DustRates"){
	name[n]=fields[3]
	type[n]="d"
    }

# IMPORTANT: this ordering specifies the order in which arrays should appear
# in variable definitions.
# Flag to decide whether or not this entry contains Rings
    flagRings[n]=match(name[n],/\[RNUM\]/)
    if (flagRings[n]>0) flagRings[n]=1
# Flag to decide whether or not this entry contains star formation history
    flagSFH[n]=match(name[n],/\[SFH_NBIN\]/)
    if (flagSFH[n]>0) flagSFH[n]=1
# Flag to decide whether or not this entry contains an element array
    flagElements[n]=match(name[n],/\[NUM_ELEMENTS\]/)
    if (flagElements[n]>0) flagElements[n]=1
# Flag to see if this has an array of size 3
    flag3[n]=match(line,/\[3\]/)
    if (flag3[n]>0) flag3[n]=1
# Flag to see if this is a magnitude
# Note: at this point NMAG has been substituted, but not any of the other parameterised arrays sizes - why not?
    flagMag[n]=match(name[n],/Mag/)
    if (flagMag[n]>0) flagMag[n]=1

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
    print "// The number of fields in the data"
    print "int nfields=" n ";"
    print " "

    print "// The field names"
    print "const char * field_names[]={"
    for(i=0;i<n;i++) print "\"" name[i] "\","
    print "};"
    print " "

    print "// Information describing the datatypes and array flags"
    print "// that we will use to construct HDF5 field_types."
    print " "
	
    print "char types[]={"
    for(i=0;i<n;i++) print "'" type[i] "',"
    print "}; "
    print " "

    print "int flag3[]={"
    for(i=0;i<n;i++) print flag3[i] ","
    print "}; "
    print " "

    print "int flagMag[]={"
    for(i=0;i<n;i++) print flagMag[i] ","
    print "}; "
    print " "

    print "int flagRings[]={"
    for(i=0;i<n;i++) print flagRings[i] ","
    print "}; "
    print " "

    print "int flagSFH[]={"
    for(i=0;i<n;i++) print flagSFH[i] ","
    print "}; "
    print " "

    print "int flagElements[]={"
    for(i=0;i<n;i++) print flagElements[i] ","
    print "}; "

}
