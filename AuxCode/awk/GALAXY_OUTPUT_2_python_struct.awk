# from a cleaned list of variables defining the GALAXY_OUPUT C struct, define a default IDL struct definition
BEGIN{
    print "# numpy dtype for LGAL_GAL_STRUCT"
    print "import numpy"
    print "struct_dtype = numpy.dtype(["
    n=0
    size=0
}
{
    line=$0
    split(line,fields)

    type=fields[1]
    name=fields[2]
    arraysize=1

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
    
    if(type == "float") {type= "numpy.float32"} 
    else if(type == "int"){type= "numpy.int32"}
    else if(type == "long" && name == "long") {
    	type="numpy.int64"
    	name=fields[3]
    }
    
    ia=match(line,/\[.*\]/)
    if(ia>0) 
    {
	arraysize=substr(line,ia+1)
	ic=index(arraysize,"]")
	if(ic > 0)
	    arraysize=substr(arraysize,0,ic-1)
    }
    ia=match(name,/\[.*\]/)
    if(ia>0) 
	name=substr(name,0,ia-1)
    ia=match(name,";")
    if(ia>0) 
	name=substr(name,0,ia-1)

    size += arraysize * dsize
    print "('" name "'," type "," arraysize "),"
    n+=1
}
END{
    print "('ending','i4',0)"
    print "])"
    print "properties_used = {}"
    print "for el in struct_dtype.names:"
    print "\tproperties_used[el] = False"
}
