# from a cleaned list of variables defining the GALAXY_OUPUT C struct, define a default table definition
BEGIN{

print ";;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;"
print "PRO LGalaxy_zerofloats, LGs "
print "; test whether floats are NaN or too small for SQLServer"
print "; if so, set offending values to 0"
print "; assumes the existence of a function testFloat accepting an array of floats"
}
{
    line=$0
    split(line,fields)

	type=fields[1]
	t=type
	name=fields[2]
	if(type == "float") {
		arraysize=1
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
				print " sel = testFloat(LGs." name ")"
				print " if(sel(0) gt -1) then begin"
				print "     LGs[sel]." name " = 0"
				print " endif" 
		}
		else {
			for(ia=0;ia<arraysize;ia++)
			{
				print " sel = testFloat(LGs." name "(" ia "))"
				print " if(sel(0) gt -1) then begin"
				print "     LGs[sel]." name "(" ia ") = 0"
				print " endif" 
			}
		}
	}
}
END{
print "end"
}