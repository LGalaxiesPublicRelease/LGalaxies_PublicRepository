# from a cleaned list of variables defining the GALAXY_OUPUT C struct, define a default table definition
BEGIN{

print ";;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;"
print "PRO LGalaxy_hist, LGs "
print "; plot for each field in the input LGalaxy struct a histogram"
}
{
    line=$0
    split(line,fields)

	type=fields[1]
	t=type
	name=fields[2]
	if(type == "float") {t= "F"} 
    else if(type == "int"){t= "I"}
    else if(type == "long" && name == "long") {t="L"}
    else if(type == "double") {t="D"}
    else if(type == "short") {t="L"}
    
	arraysize=1
    
	if(type == "long" && name == "long") {
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

	arraysize=arraysize+0 # make it explicitly numeric, 
	                      # or break if arraysize is not numeric, say due to missing NMAG declaration
		
	if(arraysize == 1)
	{
		print " print, 'plotting  " name "'"
		print "   plotHist, LGs." name ",'" name "','" t "'"
	}
	else {
		for(ia=0;ia<arraysize;ia++)
		{
			print " print, 'plotting  " name "[" ia "]'"
			print "   plotHist, LGs." name "(" ia "),'" name "[" ia "]','" t "'"
		}
	}
}
END{
print "end"
}