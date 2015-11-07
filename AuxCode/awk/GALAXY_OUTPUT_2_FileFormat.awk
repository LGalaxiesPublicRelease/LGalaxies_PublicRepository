# from a cleaned list of variables defining the GALAXY_OUPUT 
# C struct, define a CSV file describing all columns
BEGIN{
; define an LGalaxy struct
n=0
}
{
    line=$0
    split(line,fields)

	type=fields[1]
	name=fields[2]
	arraysize=1

    if(type == "long" && name == "long") {
    	type="long"
    	name=fields[3]
    }

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
		
    print name "," type "," dsize "," arraysize		
}
END{
}