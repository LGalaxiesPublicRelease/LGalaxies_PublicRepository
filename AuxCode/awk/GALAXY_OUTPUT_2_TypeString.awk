# from a cleaned list of variables defining the GALAXY_OUPUT C struct, construct a string whose characters
# indicate the size
BEGIN{
d=""
}
{

    line=$0
    split(line,fields)

	type=tolower(fields[1])
	name=tolower(fields[2])
	arraysize=1

    
    if(type == "float") {type= "F"} 
    else if(type == "int"){type= "I"}
    else if(type == "long" && name == "long") {type="L"}
    else if(type == "double") {type="D"}
    else if(type == "short") {type="L"}
    
    ia=match(line,/\[.*\]/)
	if(ia>0) 
	{
	  arraysize=substr(line,ia+1)
	  ic=index(arraysize,"]")
	  if(ic > 0)
	    arraysize=substr(arraysize,0,ic-1)
	}
	arraysize=arraysize+0 # make it explicitly numeric, 
	                      # or break if arraysize is not numeric, say due to missing NMAG declaration
	for(ia=1;ia<=arraysize;ia++){d=d""type}
}
END{
print d
}