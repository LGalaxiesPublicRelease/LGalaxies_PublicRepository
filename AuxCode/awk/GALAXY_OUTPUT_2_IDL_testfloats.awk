# from a cleaned list of variables defining the GALAXY_OUPUT C struct, define a default table definition
BEGIN{

print ";;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;"
print "FUNCTION LGalaxy_testfloats, LGs, nstart "
print "; test whether floats are NaN or too small for SQLServer"
print "; assumes the existence of a function testFloat accepting an array of floats"
print " badranges = []"
print " bad = 0"
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
				print "     bad=1"
#				print "     print, '" name " --- ', size(sel,/dimensions),':',nstart+sel"
				print "     print, '" name " --- ', nstart+sel"
				print "     print, '" name " --- ', LGs[sel]." name 
				print "     badranges=[badranges,sel]"
				print " endif" 
		}
		else {
			for(ia=0;ia<arraysize;ia++)
			{
				print " sel = testFloat(LGs." name "(" ia "))"
				print " if(sel(0) gt -1) then begin"
				print "     bad=1"
#				print "     print, '" name "[" ia "] --- ', size(sel,/dimensions),':',nstart+sel"
				print "     print, '" name "[" ia "] --- ', nstart+sel"
				print "     print, '" name " --- ', LGs[sel]." name "(" ia ")" 
				print "     badranges=[badranges,sel]"
				print " endif" 
			}
		}
	}
}
END{
print "if(bad) then begin "
print "     print, 'badranges found: ',badranges" 
#print "     badranges=badranges[uniq(badranges,sort(badranges))]"
print "endif"
print "return, badranges"
print "end"
}