# from a cleaned list of variables defining the GALAXY_OUPUT C struct, define a default IDL struct definition
BEGIN{
print ";;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;"
; define an LGalaxy struct
print "PRO LGalaxy__define"
print "tmp = {LGalaxy $"
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
    
    if(type == "float") {type= "0.0"} 
    else if(type == "int"){type= "0L"}
    else if(type == "long" && name == "long") {
    	type="0LL"
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
	if(arraysize > 1)
	{
	  if(type=="0.0"){type="fltarr(" arraysize ")"}
	}
	print ", " name " : " type " $ "
	n+=1
}
END{
print "}"
print "end"
}
