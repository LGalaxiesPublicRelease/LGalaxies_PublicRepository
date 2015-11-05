# from a cleaned list of variables defining the GALAXY_OUPUT C struct, define a default table definition
BEGIN{
print "CREATE TABLE GALAXIES ("
n=0
}
{
    if(n == 0)
      prefix = ""
    else
      prefix= ", "
    line=$0
    split(line,fields)

	type=fields[1]
	name=fields[2]
	arraysize=1


    if(type == "float") 
    {
    	type= "REAL"
    } 
    else if(type == "int")
    {
    	type= "INTEGER"
    }
    else if(type == "short")
    {
    	type= "SMALLINT"
    }
    else if(type == "long" && name == "long") {
    	type="BIGINT"
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
	name=tolower(substr(name,0,1))""substr(name,2)
	if(arraysize == 1)
	{
		print prefix " " name " " type " NOT NULL "
		n+=1
	}
	else {
	  arraysize= arraysize+0 # ensure next comparison between ia and arraysize is done numerical iso alphabetical 
      for(ia=1;ia <= arraysize;ia++)
		{
			print prefix " " name "_" ia " " type " NOT NULL "
			n+=1
		}
	}
}
END{
print ")"
}